
"""pas_seq.py
    Module for running poly-A site sequencing for alternative polyadenylation
    experiments.
"""


import itertools
import subprocess

from ruffus import (transform, split, collate, regex, suffix, add_inputs,
                        follows, files, jobs_limit)
from ruffus.task import active_if
import fisher
import pyper
import scipy.stats
from matplotlib import pyplot

from hts_waterworks.utils.ruffus_utils import (sys_call, main_logger as log,
                                               main_mutex as log_mtx)
from hts_waterworks.utils.makeGeneStructure import parse_gene_line
import hts_waterworks.mapping as mapping
from hts_waterworks.annotation import get_refseq_genes
from hts_waterworks.bootstrap import cfg
from hts_waterworks.utils.pas_seq_expression import (group_reads_by_gene,
                                                     group_adjacent_reads)
from hts_waterworks.utils.common import breakpoint


@active_if(False)
@files(None, '%s.polyA_DB' % cfg.get('DEFAULT', 'genome'), cfg.get('DEFAULT', 'genome'))
def get_polyA_DB(_, out_db, genome_build):
    cmd = r"curl 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/polyaDb.txt.gz' | gunzip - | cut -d $'\t' -f 2- > %s"
    cmd = cmd % (genome_build, out_db)
    sys_call(cmd, file_log=False)


@active_if(cfg.getint('PAS-Seq', 'min_read_count') > 0)
@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@transform(mapping.all_mappers_output, suffix('.mapped_reads'),
           '.overlap.mapped_reads', cfg.getint('PAS-Seq', 'min_read_count'))
def remove_nonoverlapping_reads(in_bed, out_bed, min_read_count):
    """
    Remove mapped reads that don't overlap with at least *min_read_count* reads
    """
    cmd = "intersectBed -wa -c -a %s -b %s | awk '$(NF) >= %s' |" \
          r"cut -f 1,2,3,4,5,6 > %s" % (in_bed, in_bed, min_read_count + 1,
                                        out_bed)
    sys_call(cmd, file_log=False)


@active_if(cfg.getboolean('PAS-Seq', 'merge_adjacent_reads'))
#@split(mapping.all_mappers_output, regex('(.*).mapped_reads$'),
@split(remove_nonoverlapping_reads, regex('(.*).mapped_reads$'),
           [r'\1.merged.mapped_reads', r'\1.merged.pileup_reads'],
           cfg.getint('PAS-Seq', 'merge_window_width'),
           cfg.getint('PAS-Seq', 'merge_num_iterations'),
           r'\1.merged.mapped_reads', r'\1.merged.pileup_reads',
           cfg.getint('PAS-Seq', 'min_read_count'))
def merge_adjacent_reads(in_bed, out_pattern, window_width, iterations,
                         out_merged, out_pileup, min_read_count):
    """Reassign read ends to a weighted average of adjacent reads"""
    # helper functions for parsing bed files
    filter_lines = lambda l: l.strip() and (not l.startswith('#') or \
                                            l.startswith('"'))
    read_bed_lines = lambda infile: itertools.ifilter(filter_lines, infile)
    
    # sort the input by chrom, stop
    tmpfile = in_bed + '.merged_adjacent_sorted'
    cmd = r"sort -t $'\t' -k 1,1 -k 3g,3 %s > %s" % (in_bed, tmpfile)
    print cmd
    sys_call(cmd, file_log=False)
    p_file = tmpfile
    outfile_pileup = None  # used on last iteration to generate the final pileup
    
    for i in range(iterations):
        print 'merge iteration %s' % i
        # read in from output of previous iteration
        infile = read_bed_lines(open(p_file))
        
        # output to a temp file except on the last iteration
        if i != iterations - 1:
            p_file = in_bed + '.merge_adjacent_%s' % i
        else:
            p_file = out_merged
            outfile_pileup = open(out_pileup, 'w')
        outfile = open(p_file, 'w')

        # parse first line
        (chrom, start, stop, name,
                        score, strand) = infile.next().strip().split('\t')[:6]
        if strand == '+':
            p_chrom, p_stops, p_names, p_strands = (chrom, [int(stop)],
                                                    [name], [strand])
        else:
            p_chrom, p_stops, p_names, p_strands = (chrom, [int(start)],
                                                    [name], [strand])
        print 'first line:', chrom, start, stop, name, score, strand
        
        for index, line in enumerate(infile):
            try:
                (chrom, start, stop,
                    name, score, strand) = line.strip().split('\t')[:6]
            except:
                print index, 'this line:', line
                raise
            if strand == '+':
                stop = int(stop)
            else:
                stop = int(start) + 1
            # is next read too far from first recorded?
            if p_chrom != chrom or (len(p_stops) > 0 and
                                    abs(p_stops[0] - stop) > window_width):
                if len(p_stops) == 0 or len(p_names) == 0:
                    print 'error!'
                    print line
                    print p_stops, p_names, p_strands
                    raise
                if len(p_stops) > min_read_count:
                    avg = int(round(sum(p_stops) / float(len(p_stops))))
                    # write out reads in this cluster, using avg as coordinate
                    outfile.writelines('\t'.join([p_chrom, str(max(0, avg-1)), str(avg),
                                             n_name, '0', n_strand]) + '\n'
                                  for n_name, n_strand in zip(p_names, p_strands))
                    if outfile_pileup is not None:
                        outfile_pileup.write('\t'.join([p_chrom, str(max(0, avg-1)), str(avg),
                                               p_names[0], str(len(p_stops)),
                                               p_strands[0]]) + '\n')
                # reset our record
                p_chrom = chrom
                p_stops = [stop]
                p_names =  [name]
                p_strands = [strand]
            # otherwise, the next read is within the window, on same chrom
            else:
                p_stops.append(stop)
                p_names.append(name)
                p_strands.append(strand)

        # output anything left in queue after EOF
        if len(p_stops) > 0:
            avg = int(round(sum(p_stops) / float(len(p_stops))))
            # write out reads in this cluster, using avg as coordinate
            outfile.writelines('\t'.join([chrom, str(max(0, avg-1)), str(avg),
                                     n_name, '0', n_strand]) + '\n'
                          for n_name, n_strand in zip(p_names, p_strands))
            if outfile_pileup is not None:
                outfile_pileup.write('\t'.join([chrom, str(max(0, avg-1)), str(avg),
                                           p_names[0], str(len(p_stops)),
                                           p_strands[0]]) + '\n')
        if outfile_pileup is not None:
            outfile_pileup.close()
        outfile.close()


@transform(merge_adjacent_reads , regex(r'(.*)\.(.*$)'), r'\1.no_prime.\2')
def remove_internal_priming_again(in_bed, out_bed):
    mapping.remove_internal_priming(in_bed, out_bed)

#@transform('hg19.refseq_genes', suffix('hg19.refseq_genes'),
#           'hg19.refseq_genes.middle_exons')
@transform(get_refseq_genes, suffix('.refseq_genes'),
           '.refseq_genes.middle_exons')
def make_middle_exons(in_refseq, out_exons):
    with open(out_exons, 'w') as outfile:
        for line in open(in_refseq):
            (name, chrom, strand, txStart, txEnd, cdsStart,
                cdsEnd, exons, name2, noncoding) = parse_gene_line(line)
            # remove the last exon
            if strand == '+':
                exons = exons[:-1]
            else:
                exons = exons[1:]
            for ex_start, ex_end in exons:
                outfile.write('\t'.join(map(str, [chrom, ex_start, ex_end]))
                                + '\n')

@follows(make_middle_exons)
@transform(remove_internal_priming_again, regex(r'(.*)\.(.*$)'),
           add_inputs(make_middle_exons), r'\1.no_exons.\2')
def remove_terminal_exon(in_files, out_bed):
    """Remove all exons but the last one using intersectBed"""
    in_bed, exon_file = in_files
    cmd = 'intersectBed -v -a %s -b %s > %s' % (in_bed, exon_file, out_bed)
    sys_call(cmd, file_log=False)


@active_if(cfg.getboolean('visualization', 'normalize_per_million'))
@transform(remove_terminal_exon, regex(r'(.*)\.(pileup_reads$)'),
           r'\1.norm_mil.\2')
def pileup_normalize_per_million(in_pileup, out_pileup):
    """Normalize pileup reads to tags per million mapping"""
    total_reads = sum(float(l.strip().split('\t')[4])
                                    for l in open(in_pileup))
    with open(in_pileup) as infile:
        with open(out_pileup, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                norm_score = float(fields[4]) / (total_reads / 1e6)
                outfile.write('\t'.join(fields[:4] + [str(norm_score)] +
                                        fields[5:]) + '\n')


short_name = lambda x: x.replace('hg19.refseq_genes.', '').split('.trim_regex')[0] + ('.plus' if 'plus' in x else '.minus')

@active_if(cfg.getboolean('PAS-Seq', 'test_differential_polya'))
@follows(remove_terminal_exon, pileup_normalize_per_million)
@split(get_refseq_genes, regex(r'(.*)'),
       add_inputs(pileup_normalize_per_million if
                  cfg.getboolean('visualization', 'normalize_per_million') else
                  remove_terminal_exon),
#@split(get_refseq_genes, regex(r'(.*)'), add_inputs('*.no_prime.norm_mil.pileup_reads'),
#@split('hg19.refseq_genes.extend3utr', regex(r'(.*)'), add_inputs('*.pileup_reads'),
           r'\1.*.polya.*.*_test',
           cfg.getint('PAS-Seq', 'compare_window_width'),
           r'\1.%s.vs.%s.polya.%s.fisher_test',
           r'\1.%s.vs.%s.polya.%s.t_test',
           r'\1.%s.vs.%s.polya.%s.avg_fisher_test',
           cfg.getfloat('PAS-Seq', 'min_score_for_site'))
def test_differential_polya(in_files, out_pattern, max_dist, out_template,
                            ttest_template, avg_fisher_template, min_score):
    """Test for differential poly-adenylation from PAS-seq pileups.
    
    Performs all pairwise tests of merged poly-A sites across all experiments.
    """
    print in_files
    in_genes, all_reads = in_files[0], in_files[1:]
    all_reads = filter(lambda f: f.endswith('pileup_reads'), all_reads)
    all_reads = sorted(all_reads)
    #read_counts = map(lambda f: sum(1 for i in open(f)), all_reads)
    if len(all_reads) == 0:
        raise RuntimeError('differential polyadenylation requires multiple '
                           'input datasets! I saw %s ', cur_reads)
    out_files = {}
    out_ttest = {}
    out_avg_fisher = {}
    total_sense = 0
    total_antisense = 0
    skipped_sense = 0
    skipped_antisense = 0
    loc_sense = 0
    loc_antisense = 0
    for compare_type in ['non_coding', 'independent', 'dependent', 'non_utr', 'utr']:
        out_files[compare_type] = {}
        out_ttest[compare_type] = {}
        out_avg_fisher[compare_type] = {}
    for strand in ['plus', 'minus']:
        strand_sign = '+' if strand == 'plus' else '-'
        cur_reads = filter(lambda f: strand in f, all_reads)
        if len(cur_reads) == 0:
            continue
        out_summary = open('fisher.test.%s.summary' % strand, 'w')
        out_summary.write('\t'.join(['%s_%s' % (c.split('.')[0], strand)
                                     for c in cur_reads]) + '\n')
        
        # collect the samples that belong to the same group (TM_CA and TM_GA are same)
        t_similar_samples = {}
        for index, exp_name in enumerate(cur_reads):
            similar_exp_name = exp_name.replace('CA','').replace('GA','')
            if similar_exp_name not in t_similar_samples:
                t_similar_samples[similar_exp_name] = []
            t_similar_samples[similar_exp_name].append(index)
        print '\n\n\n\here goes'
        #for g, inds in t_similar_samples.items():
        #    print g, [cur_reads[i] for i in inds]
        # group each set of reads by the gene it overlaps
        reads_by_gene = [group_reads_by_gene(in_genes, r) for r in cur_reads]
        #print reads_by_gene
        try:
            while True:
                # get the reads from each group associated with the next gene
                read_groups = []
                p_gene_fields = None
                for gene_iter in reads_by_gene:
                    gene_fields, read_fields = gene_iter.next()
                    same_strand_as_gene = 'sense' if gene_fields[2] == strand_sign else 'antisense'
                    #print same_strand_as_gene, gene_fields[2], strand_sign
                    read_groups.append(read_fields)
                    if p_gene_fields is None:
                        p_gene_fields = gene_fields
                    elif p_gene_fields != gene_fields:
                        raise RuntimeError('Parse error! gene fields %s and %s '
                                           'should match!', p_gene_fields,
                                           gene_fields)
                # gene_fields: (name, chrom, strand, txStart, txEnd, cdsStart,
                #                               cdsEnd, exons, name2, noncoding)
                #print read_groups
                expr_vals = dict(group_adjacent_reads(read_groups, max_dist,
                                                      min_score=min_score))
                
                if any(map(len, read_groups)):
                    if same_strand_as_gene == 'sense':
                        total_sense += 1
                    else:
                        total_antisense += 1
                        if (total_antisense / 1000) == 1:
                            breakpoint()
                #print expr_vals
                # sort locs s.t. those closest to the gene TSS come first
                locs = sorted(expr_vals.keys(), reverse=gene_fields[2] == '-')
                
                # save expr vals to disk
                out_summary.writelines('%s\t%s\n' % (l,
                                            '\t'.join(map(str, expr_vals[l])))
                                       for l in locs)
                
                
                # do a t-test on each location, pooling estimates from the same experiments but different libraries
                for cur_loc in locs:
                    if same_strand_as_gene == 'sense':
                        loc_sense += 1
                    else:
                        loc_antisense += 1
                    if gene_fields[9]:  # non-coding gene
                        compare_type = 'non_coding'
                    else:
                        # check if cur_loc is in the 3'utr
                        if gene_fields[2] == '+':
                            utr3_left = gene_fields[6]  # cdsEnd
                            utr3_right = gene_fields[4]  # txEnd
                        else:
                            utr3_left = gene_fields[3]  # txStart
                            utr3_right = gene_fields[5]  # cdsStart
                        in_utr = (utr3_left <= cur_loc[1] < utr3_right or
                                            (cur_loc[1] <= utr3_left and
                                                cur_loc[2] >= utr3_right))
                        if in_utr:
                            compare_type = 'utr'
                        else:
                            compare_type = 'non_utr'
                    #if (expr_vals[cur_loc][exp1] < min_score and 
                    #        expr_vals[cur_loc][exp2] < min_score):
                    #    continue  # both experiments have very little expression
                    #if (('plus' in cur_reads[exp1] and
                    #        'minus' in cur_reads[exp2]) or
                    #    ('minus' in cur_reads[exp1] and
                    #        'plus' in cur_reads[exp2])):
                    #    continue  # skip comparisons between different strands
                    #if cur_reads[exp1] == cur_reads[exp2]:
                    #    continue  # skip comparisons with same file



                    # t-test on pooled libraries from each experiment
                    for group1, group2 in itertools.combinations(
                                        sorted(t_similar_samples.keys()), 2):
                    
                        g1_vals = [expr_vals[cur_loc][cur_exp]
                                    for cur_exp in t_similar_samples[group1]]
                        g2_vals = [expr_vals[cur_loc][cur_exp]
                                    for cur_exp in t_similar_samples[group2]]
                        # skip if any of the expression values are < min_score
                        if any([v < min_score for v in g1_vals + g2_vals]):
                            #print 'skipping', same_strand_as_gene
                            if same_strand_as_gene == 'sense':
                                skipped_sense += 1
                            else:
                                skipped_antisense += 1
                            continue
                        t, p = scipy.stats.ttest_ind(g1_vals, g2_vals)
                    
                        outline = '\t'.join(map(str, [gene_fields[0],
                                                cur_loc, compare_type,
                                                same_strand_as_gene,
                                                t_similar_samples[group1],
                                                g1_vals,
                                                t_similar_samples[group2],
                                                g2_vals,
                                                t, p])) + '\n'
                        # create the output file if it's not ready
                        out_name = frozenset([group1, group2])
                        #print out_name
                        if out_name not in out_ttest[compare_type]:
                            file_name = ttest_template % ('%s.%s' % (group1.split('.')[0], strand),
                                                        '%s.%s' % (group2.split('.')[0], strand),
                                                        compare_type)
                            #file_name = ttest_template % ('%s.%s' % (short_name(group1), strand),
                            #                            '%s.%s' % (short_name(group2), strand),
                            #                            compare_type)
                            out_ttest[compare_type][out_name] = open(file_name, 'w')
                            out_ttest[compare_type][out_name].write('\t'.join([
                                'gene_name', 'loc', 'compare_type',
                                'sense_with_gene', 'exp1_name',
                                'exp1_counts',  'exp2_name', 'exp2_counts',
                                'ttest_stat', 'ttest_pvalue']) + '\n')
                        out_ttest[compare_type][out_name].write(outline)




                # get average expression across similar experiments, then do a fisher's test
                for group1, group2 in itertools.combinations(
                                                sorted(t_similar_samples.keys()), 2):
                    for ups_loc, downs_loc in itertools.combinations(locs, 2):
                        g1_up_val = scipy.mean([expr_vals[ups_loc][cur_exp]
                                    for cur_exp in t_similar_samples[group1]])
                        g1_down_val = scipy.mean([expr_vals[downs_loc][cur_exp]
                                    for cur_exp in t_similar_samples[group1]])
                        g2_up_val = scipy.mean([expr_vals[ups_loc][cur_exp]
                                    for cur_exp in t_similar_samples[group2]])
                        g2_down_val = scipy.mean([expr_vals[downs_loc][cur_exp]
                                    for cur_exp in t_similar_samples[group2]])

                        if gene_fields[9]:  # non-coding gene
                            compare_type = 'non_coding'
                        else:
                            # check if both ups_loc and downs_loc are in the 3'utr
                            if gene_fields[2] == '+':
                                utr3_left = gene_fields[6]  # cdsEnd
                                utr3_right = gene_fields[4]  # txEnd
                            else:
                                utr3_left = gene_fields[3]  # txStart
                                utr3_right = gene_fields[5]  # cdsStart
                            ups_in_utr = (utr3_left <= ups_loc[1] < utr3_right or
                                                (ups_loc[1] <= utr3_left and
                                                    ups_loc[2] >= utr3_right))
                            downs_in_utr = (utr3_left <= downs_loc[1] < utr3_right or
                                                (downs_loc[1] <= utr3_left and
                                                    downs_loc[2] >= utr3_right))
                            if ups_in_utr and downs_in_utr:
                                compare_type = 'independent'
                            elif ups_in_utr or downs_in_utr:
                                compare_type = 'dependent'
                            else:
                                compare_type = 'non_utr'

                        if (g1_up_val < min_score or
                            g1_down_val < min_score or
                            g2_up_val < min_score or
                            g2_down_val < min_score):
                            # skip where any site for any experiment has too
                            # low expression at both locations
                            continue
                        if ((g1_up_val < min_score and
                                g1_down_val < min_score) or
                            (g2_up_val < min_score and
                                g2_down_val < min_score)):
                            # skip where experiments have too low expression at
                            # both locations
                            continue
                        if ((g1_up_val < min_score and
                                g2_up_val < min_score) or
                            (g1_down_val < min_score and
                                g2_down_val < min_score)):
                            # skip where locations have too low expression in
                            # both experiments
                            continue
                        #if (('plus' in cur_reads[exp1] and
                        #        'minus' in cur_reads[exp2]) or
                        #    ('minus' in cur_reads[exp1] and
                        #        'plus' in cur_reads[exp2])):
                        #    continue  # skip comparisons between different strands
                        #if cur_reads[exp1] == cur_reads[exp2]:
                        #    continue  # skip comparisons with same file
                        p = fisher.pvalue(g1_up_val,
                                          g1_down_val,
                                          g2_up_val,
                                          g2_down_val)
                        #if gene_fields[0] == 'NM_001020':
                        #    breakpoint()
                        #if downs_loc == (('chr17', 34053421, 34053422)):
                        #    breakpoint()
                        outline = '\t'.join(map(str, [gene_fields[0],
                                                    ups_loc, downs_loc,
                                                    compare_type, same_strand_as_gene,
                                                    t_similar_samples[group1],
                                                    g1_up_val,
                                                    g1_down_val,
                                                    t_similar_samples[group2],
                                                    g2_up_val,
                                                    g2_down_val,
                                                    p.left_tail, p.right_tail,
                                                    p.two_tail])) + '\n'
                        # create the output file if it's not ready
                        out_name = frozenset([cur_reads[i] for i in t_similar_samples[group1]] + [cur_reads[i] for i in t_similar_samples[group2]])
                        if out_name not in out_avg_fisher[compare_type]:
                            #print "new created", out_name, compare_type
                            file_name = avg_fisher_template % ('%s' % short_name(group1),
                                                        '%s' % short_name(group2),
                                                        compare_type)
                            out_avg_fisher[compare_type][out_name] = open(file_name, 'w')
                            out_avg_fisher[compare_type][out_name].write('\t'.join([
                                'gene_name', 'ups_loc', 'downs_loc',
                                'compare_type', 'sense_with_gene', 'exp1_name',
                                'exp1_upstream_count', 'exp1_downstream_count',
                                'exp2_name', 'exp2_upstream_count',
                                'exp2_downstream_count', 'fisher_p_left',
                                'fisher_p_right', 'fisher_p_two_sided']) + '\n')
                        out_avg_fisher[compare_type][out_name].write(outline)




                
                # for each combination of experiment, do a fisher's test
                for exp1, exp2 in itertools.combinations(
                                                    range(len(cur_reads)), 2):
                    # do fisher-tests on all combinations of read locations    
                    for ups_loc, downs_loc in itertools.combinations(locs, 2):
                        # is the poly-adenylation splice-dependent?
                        #       i.e., in same annotated 3'utr?
                        if gene_fields[9]:  # non-coding gene
                            compare_type = 'non_coding'
                        else:
                            # check if both ups_loc and downs_loc are in the 3'utr
                            if gene_fields[2] == '+':
                                utr3_left = gene_fields[6]  # cdsEnd
                                utr3_right = gene_fields[4]  # txEnd
                            else:
                                utr3_left = gene_fields[3]  # txStart
                                utr3_right = gene_fields[5]  # cdsStart
                            ups_in_utr = (utr3_left <= ups_loc[1] < utr3_right or
                                                (ups_loc[1] <= utr3_left and
                                                    ups_loc[2] >= utr3_right))
                            downs_in_utr = (utr3_left <= downs_loc[1] < utr3_right or
                                                (downs_loc[1] <= utr3_left and
                                                    downs_loc[2] >= utr3_right))
                            if ups_in_utr and downs_in_utr:
                                compare_type = 'independent'
                            elif ups_in_utr or downs_in_utr:
                                compare_type = 'dependent'
                            else:
                                compare_type = 'non_utr'

                        if (expr_vals[ups_loc][exp1] < min_score or
                            expr_vals[downs_loc][exp1] < min_score or
                            expr_vals[ups_loc][exp2] < min_score or
                            expr_vals[downs_loc][exp2] < min_score):
                            # skip where any site for any experiment has too
                            # low expression at both locations
                            continue
                        if ((expr_vals[ups_loc][exp1] < min_score and
                                expr_vals[downs_loc][exp1] < min_score) or
                            (expr_vals[ups_loc][exp2] < min_score and
                                expr_vals[downs_loc][exp2] < min_score)):
                            # skip where experiments have too low expression at
                            # both locations
                            continue
                        if ((expr_vals[ups_loc][exp1] < min_score and
                                expr_vals[ups_loc][exp2] < min_score) or
                            (expr_vals[downs_loc][exp1] < min_score and
                                expr_vals[downs_loc][exp2] < min_score)):
                            # skip where locations have too low expression in
                            # both experiments
                            continue
                        if (('plus' in cur_reads[exp1] and
                                'minus' in cur_reads[exp2]) or
                            ('minus' in cur_reads[exp1] and
                                'plus' in cur_reads[exp2])):
                            continue  # skip comparisons between different strands
                        if cur_reads[exp1] == cur_reads[exp2]:
                            continue  # skip comparisons with same file
                        p = fisher.pvalue(expr_vals[ups_loc][exp1],
                                          expr_vals[downs_loc][exp1],
                                          expr_vals[ups_loc][exp2],
                                          expr_vals[downs_loc][exp2])
                        #if gene_fields[0] == 'NM_001020':
                        #    breakpoint()
                        #if downs_loc == (('chr17', 34053421, 34053422)):
                        #    breakpoint()
                        outline = '\t'.join(map(str, [gene_fields[0],
                                                    ups_loc, downs_loc,
                                                    compare_type, same_strand_as_gene,
                                                    cur_reads[exp1],
                                                    expr_vals[ups_loc][exp1],
                                                    expr_vals[downs_loc][exp1],
                                                    cur_reads[exp2],
                                                    expr_vals[ups_loc][exp2],
                                                    expr_vals[downs_loc][exp2],
                                                    p.left_tail, p.right_tail,
                                                    p.two_tail])) + '\n'
                        # create the output file if it's not ready
                        out_name = frozenset([cur_reads[exp1], cur_reads[exp2]])
                        #print out_name
                        if out_name not in out_files[compare_type]:
                            file_name = out_template % ('%s.%s' % (cur_reads[exp1].split('.')[0], strand),
                                                        '%s.%s' % (cur_reads[exp2].split('.')[0], strand),
                                                        compare_type)
                            out_files[compare_type][out_name] = open(file_name, 'w')
                            out_files[compare_type][out_name].write('\t'.join([
                                'gene_name', 'ups_loc', 'downs_loc',
                                'compare_type', 'sense_with_gene', 'exp1_name',
                                'exp1_upstream_count', 'exp1_downstream_count',
                                'exp2_name', 'exp2_upstream_count',
                                'exp2_downstream_count', 'fisher_p_left',
                                'fisher_p_right', 'fisher_p_two_sided']) + '\n')
                        out_files[compare_type][out_name].write(outline)
        except StopIteration:  # no more genes
            pass
        
        out_summary.close()
        
    for filedict in out_files.values():
        for fileout in filedict.values():
            fileout.close()
    
 

@collate(test_differential_polya,
         #regex(r'(.*)\.(?P<strand>plus|minus)\.(.*)\.(?P=strand)\.(.*)'), r'\1.\3.\5')
         regex(r'(.*)\.(?P<strand>plus|minus)\.(.*)\.(?P=strand)\.(.*)'), r'\1.\3.\4')
def merge_strands(in_files, out_merged):
    """concatenate the strand-specific analyses for plotting"""
    # output the first file in its entirety
    cmd = 'cat %s > %s ' % (in_files[0], out_merged)
    sys_call(cmd, file_log=False)
    # skip the header for remaining files
    for f in in_files[1:]:
        cmd = 'sed 1d %s >> %s' % (f, out_merged)
        sys_call(cmd, file_log=False)

@collate(merge_strands,
         regex(r'(.*\.polya)\.(.*)\.(t_test|fisher_test|avg_fisher_test)'), r'\1.\3')
def merge_comparison_types(in_files, out_merged):
    """concatenate the comparison types together for plotting"""
    cmd = 'cat %s > %s ' % (in_files[0], out_merged)
    sys_call(cmd, file_log=False)
    # skip the header for remaining files
    for f in in_files[1:]:
        cmd = 'sed 1d %s >> %s' % (f, out_merged)
        sys_call(cmd, file_log=False)


@collate(merge_strands,
         regex(r'(.*refseq_genes)\.(\w+)_(?P<exp>\w+)\.vs\.(\w+)_(?P=exp)\.(.*fisher_test)'),
         r'\1.\2.vs.\4.agreement.\5')
def intersect_comparison_types(in_files, out_merged):
    """
    Report those sites that are significant in ALL libraries for a given
    experiment
    
    The collate is formed s.t. significance must be attained across library
    types.  i.e.,
    Job = [[hg19.refseq_genes.CHX_CA.vs.DMSO_CA.polya.independent.fisher_test,
            hg19.refseq_genes.CHX_GA.vs.DMSO_GA.polya.independent.fisher_test]
    -> hg19.refseq_genes.CHX.vs.DMSO.agreement.polya.independent.fisher_test]
    
    which does not include CHX_CA.vs.DMSO_GA etc
    """
    site_pvals = {}
    CA_index = 0
    GA_index = 1
    for exp in in_files:
        for line in open(exp):
            fields = line.strip().split('\t')
            if fields[-1] == 'fisher_p_two_sided':
                continue
            loc = tuple(fields[1:3])
            pval = float(fields[-1])
            if loc not in site_pvals:
                site_pvals[loc] = []
            site_pvals[loc].append(float(pval))
    with open(out_merged, 'w') as outfile:
        for loc, pvals in site_pvals.iteritems():
            if all(p < .05 for p in pvals) and len(pvals) == len(in_files):
                outfile.write('\t'.join(['%s:%s-%s' %tuple(l.replace("'",'').replace(' ','').split(',')) for l in loc] + [str(pvals)]) + '\n')


@split([merge_strands, merge_comparison_types], regex(r'(.*)'),
    r'\1.ratio_sites.*.png', r'\1.ratio_sites.*.%s.png')
def plot_differential_polya(in_fisher, out_pattern, out_template):
    """plot a scatter plot of the log-expression difference of polya sites"""
    lhs, rhs = in_fisher.split('vs.')
    lhs, rhs = short_name(lhs), short_name(rhs)
    try:
        compare_type = in_fisher.replace('.fisher_test', '').split('polya.')[1]
    except:
        compare_type = 'ALL'
    for max_pval in [.05, .01, .001]:
        out_png = out_template % ('pval_%s' % max_pval)
        R_script = r"""
library(lattice)
d<-read.table(file="%(in_fisher)s", header=TRUE, sep="\t")
png('%(out_png)s')
sig_sites <- d$fisher_p_two_sided < %(max_pval)s

exp1_proximal = d$exp1_upstream_count[sig_sites]
exp1_distal = d$exp1_downstream_count[sig_sites]
exp2_proximal = d$exp2_upstream_count[sig_sites]
exp2_distal = d$exp2_downstream_count[sig_sites]

plot(log2(d$exp1_upstream_count/d$exp2_upstream_count), log2(d$exp1_downstream_count/d$exp2_downstream_count), cex=.8, col='lightgray', pch=20, xlab="%(xlab)s", ylab="%(ylab)s", main="Significant sites for %(lhs)s vs. %(rhs)s in %(compare_type)s", sub=paste("Significant sites:", sum(sig_sites), "/", dim(d)[1]))
points(log2(exp1_proximal/exp2_proximal), log2(exp1_distal/exp2_distal), col='red', cex=.8, pch=20)

dev.off()
""" % dict(plot_label=r'Differential Poly-A for\n%s' % in_fisher,
               xlab="log2(%s/%s)-proximal" % (lhs, rhs),
               ylab="log2(%s/%s)-distal" % (lhs, rhs),
               in_fisher=in_fisher, out_png=out_png, lhs=lhs, rhs=rhs,
               compare_type=compare_type, max_pval=max_pval)
        # print R_script
        r = pyper.R()
        r(R_script)


@transform([merge_strands, merge_comparison_types], suffix(''), '.scatter.png')
def plot_scatter_polya(in_fisher, out_png):
    """plot a scatter plot of the expression of polya sites"""
    lhs, rhs = in_fisher.split('vs.')
    lhs, rhs = short_name(lhs), short_name(rhs)
    try:
        compare_type = in_fisher.replace('.fisher_test', '').split('polya.')[1]
    except:
        compare_type = 'ALL'
    R_script = r"""
library(lattice)
d<-read.table(file="%(in_fisher)s", header=TRUE, sep="\t")
png('%(out_png)s')

exp1 <- c(d$exp1_upstream_count, d$exp1_downstream_count)
exp2 <- c(d$exp2_upstream_count, d$exp2_downstream_count)

sig_sites <- d$fisher_p_two_sided < .05

exp1_sig = c(d$exp1_upstream_count[sig_sites], d$exp1_downstream_count[sig_sites])
exp2_sig = c(d$exp2_upstream_count[sig_sites], d$exp2_downstream_count[sig_sites])

plot(log2(exp1), log2(exp2), cex=.8, col='lightgray', pch=20, xlab="%(xlab)s", ylab="%(ylab)s", main="All sites for %(lhs)s vs. %(rhs)s in %(compare_type)s", sub=paste("R^2 is ", cor(exp1, exp2)))

dev.off()
""" % dict(plot_label=r'Poly-A for\n%s' % in_fisher,
           xlab="log2(%s)" % (lhs),
           ylab="log2(%s)" % (rhs),
           in_fisher=in_fisher, out_png=out_png, lhs=lhs, rhs=rhs,
           compare_type=compare_type)
    # print R_script
    r = pyper.R()
    r(R_script)

@transform([merge_strands, merge_comparison_types],
    regex(r'(.*)\.t_test'), r'\1.t_test.ttest_scatter.png')
def plot_ttest_polya(in_ttest, out_png):
    """plot the t-test averages used as a scatter plot of the expression of polya sites"""
    lhs, rhs = in_ttest.split('vs.')
    lhs, rhs = short_name(lhs), short_name(rhs)
    try:
        compare_type = in_ttest.replace('.t_test', '').split('polya.')[1]
    except:
        compare_type = 'ALL'
    R_script = r"""
library(lattice)
d<-read.table(file="%(in_ttest)s", header=TRUE, sep="\t")
png('%(out_png)s')

exp1 <- unlist(lapply(lapply(strsplit(gsub("]", "", gsub("[", "", d$exp1_count, fixed=TRUE), fixed=TRUE), ", ", fixed=TRUE), as.numeric), mean))
exp2 <- unlist(lapply(lapply(strsplit(gsub("]", "", gsub("[", "", d$exp2_count, fixed=TRUE), fixed=TRUE), ", ", fixed=TRUE), as.numeric), mean))


sig_sites <- d$ttest_pvalue < .05
# upregulated means t-statistic is positive => exp1 < exp2
exp1_bigger <- d$ttest_pvalue < .05 & d$ttest_stat > 0
exp2_bigger <- d$ttest_pvalue < .05 & d$ttest_stat < 0

exp1_sig = exp1[sig_sites]
exp2_sig = exp2[sig_sites]

plot(log2(exp1), log2(exp2), cex=.8, col='lightgray', pch=20, xlab="", ylab="%(ylab)s", main="All sites for %(lhs)s vs. %(rhs)s in %(compare_type)s", sub=paste("R^2 is ", cor(exp1, exp2),"\nSig sites x > y: ", sum(exp1_bigger), "\nSig sites x < y: ", sum(exp2_bigger)))
points(log2(exp1_sig), log2(exp2_sig), col='red', cex=.8, pch=20)

dev.off()
""" % dict(plot_label=r'Poly-A for\n%s' % in_ttest,
           xlab="log2(%s)" % (lhs),
           ylab="log2(%s)" % (rhs),
           in_ttest=in_ttest, out_png=out_png, lhs=lhs, rhs=rhs,
           compare_type=compare_type)
    # print R_script
    r = pyper.R()
    r(R_script)


@transform(remove_terminal_exon, suffix('.pileup_reads'),
           add_inputs(get_polyA_DB),
           [r'.closest_polyA.png',
            r'.count_near_polyA.png'], [20, 40, 60, 80, 100])
def plot_closest_polyA_db(in_files, out_files, max_dists):
    """Find the nearest known polyA site for each read"""
    print 'plot_closest_polyA_db', in_files, out_files
    in_reads, in_polyA_DB = in_files
    out_dist, out_count = out_files
    cmd = r"closestBed -a {in_reads} -b {in_polyA_DB} -d -t first | cut -d $'\t' -f 5,15"
    cmd = cmd.format(in_reads=in_reads, in_polyA_DB=in_polyA_DB)
    p = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)
    counts_at_threshold = dict([(max_d, []) for max_d in max_dists])
    all_dists = []
    furthest = max(max_dists)
    for line in p.stdout:
        try:
            count, dist = map(lambda x: int(round(float(x))),
                                            line.strip().split('\t'))
        except Exception as e:
            print 'error on line:'
            print line
            raise e
        for max_d in max_dists:
            if dist <= max_d:
                # note how many reads are nearby
                counts_at_threshold[max_d].append(count)
        #print count, dist, furthest, all_dists, ([dist] * count)
        if dist <= furthest:  # only plot reads within window
            all_dists.extend([dist] * count)
    # plot the distance from reads to polyA site
    pyplot.figure()
    pyplot.hist(all_dists, bins=min(furthest, 50))
    pyplot.xlabel('Distance to PolyA')
    pyplot.ylabel('Count')
    pyplot.title('Distance for %s' % in_reads)
    pyplot.savefig(out_dist)
    
    # plot the count of reads at each polyA site for several thresholds
    pyplot.figure(figsize=(8,13))
    pyplot.subplot(len(max_dists), 1, 1)
    pyplot.title('Number of reads at PolyA sites')
    pyplot.ylabel('Frequency (Count)')
    for index, max_d in enumerate(max_dists):
        pyplot.subplot(len(max_dists), 1, index + 1)
        pyplot.hist(counts_at_threshold[max_d],
                    bins=min(counts_at_threshold[max_d], 50),
                    range=(0,1000))
        pyplot.xlabel('# reads, max_dist=%s' % max_d)
    pyplot.savefig(out_count)


#    R_script = """
#d<-read.table(file="%(out_data)s", header=FALSE, sep="\t")
#png('%(out_png)s')
#plot(density(d), xlab="Distance (bp)", ylab="Fraction", main="Distance to known polyA sites from %(in_reads)s")
#dev.off()
#"""
#    r = pyper.R()
#    r(R_script)
    
    

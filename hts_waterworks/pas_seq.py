
"""pas_seq.py
    Module for running poly-A site sequencing for alternative polyadenylation
    experiments.
"""

#  Current Version: 0.1-2-gce84a33
#  Last Modified: 2011-07-31 00:49

import itertools 

from ruffus import (transform, split, regex, suffix, add_inputs)
from ruffus.task import active_if
import fisher

from hts_waterworks.utils.ruffus_utils import (
                                           sys_call, main_logger as log,
                                           main_mutex as log_mtx)
import hts_waterworks.mapping as mapping
from hts_waterworks.annotation import get_refseq_genes
from hts_waterworks.bootstrap import cfg
from hts_waterworks.utils.pas_seq_expression import (group_reads_by_gene,
                                                     group_adjacent_reads)




@active_if(cfg.getboolean('PAS-Seq', 'merge_adjacent_reads'))
@split(mapping.all_mappers_output, regex('(.*).mapped_reads$'),
           [r'\1.merged.mapped_reads', r'\1.merged.pileup_reads'],
           cfg.getint('PAS-Seq', 'merge_window_width'),
           cfg.getint('PAS-Seq', 'merge_num_iterations'),
           r'\1.merged.mapped_reads', r'\1.merged.pileup_reads')
def merge_adjacent_reads(in_bed, out_pattern, window_width, iterations,
                         out_merged, out_pileup):
    """Reassign read ends to a weighted average of adjacent reads"""
    # helper functions for parsing bed files
    filter_lines = lambda l: l.strip() and (not l.startswith('#') or \
                                            l.startswith('"'))
    read_bed_lines = lambda infile: itertools.ifilter(filter_lines, infile)
    
    # sort the input by chrom, stop
    tmpfile = in_bed + '.merged_adjacent_sorted'
    cmd = 'sort -k1 -k3g %s > %s' % (in_bed, tmpfile)
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
        p_chrom, p_stops, p_names, p_strands = (chrom, [int(stop)],
                                                [name], [strand])
        print 'first line:', chrom, start, stop, name, score, strand
        
        for index, line in enumerate(infile):
            try:
                (chrom, start, stop,
                    name, score, strand) = line.strip().split('\t')[:6]
            except:
                print index, 'this line:', line
                raise
            start, stop = int(start), int(stop)
            # is next read too far from first recorded?
            if p_chrom != chrom or (len(p_stops) > 0 and
                                    abs(p_stops[0] - stop) > window_width):
                if len(p_stops) == 0 or len(p_names) == 0:
                    print 'error!'
                    print line
                    print p_stops, p_names, p_strands
                    raise
                avg = int(round(sum(p_stops) / float(len(p_stops))))
                # write out reads in this cluster, using avg as coordinate
                outfile.writelines('\t'.join([chrom, str(avg), str(avg+1),
                                         n_name, '0', n_strand]) + '\n'
                              for n_name, n_strand in zip(p_names, p_strands))
                if outfile_pileup is not None:
                    outfile_pileup.write('\t'.join([chrom, str(avg), str(avg+1),
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
            outfile.writelines('\t'.join([chrom, str(avg), str(avg+1),
                                     n_name, '0', n_strand]) + '\n'
                          for n_name, n_strand in zip(p_names, p_strands))
            if outfile_pileup is not None:
                outfile_pileup.write('\t'.join([chrom, str(avg), str(avg+1),
                                           p_names[0], str(len(p_stops)),
                                           p_strands[0]]) + '\n')
        if outfile_pileup is not None:
            outfile_pileup.close()
        outfile.close()

@transform(get_refseq_genes, suffix(''), add_inputs(merge_adjacent_reads),
           '.polyadenylation.fisher_test',
           cfg.getint('PAS-Seq', 'merge_window_width'))
def test_differential_polya(in_files, out_table, max_dist):
    """Test for differential poly-adenylation from PAS-seq pileups.
    
    Performs all pairwise tests of merged poly-A sites across all experiments.
    """
    print in_files
    in_genes, all_reads = in_files[0], in_files[1:]
    all_reads = filter(lambda f: f.endswith('pileup_reads'), all_reads)
    print all_reads
    if len(all_reads) == 0:
        raise RuntimeError('differential polyadenylation requires multiple '
                           'input datasets! I saw %s ', all_reads)
    # group each set of reads by the gene it overlaps
    reads_by_gene = [group_reads_by_gene(in_genes, r) for r in all_reads]
    
    outfile = open(out_table, 'w')
    outfile.write('\t'.join(['gene_name', 'left_loc', 'right_loc',
                             'compare_type','exp1_name', 'exp1_upstream_count',
                             'exp1_downstream_count', 'exp2_name',
                             'exp2_upstream_count', 'exp2_downstream_count',
                             'fisher_p_left', 'fisher_p_right',
                             'fisher_p_two_sided']) + '\n')
    try:
        while True:
            # get the reads from each group associated with the next gene
            read_groups = []
            p_gene_fields = None
            for gene_iter in reads_by_gene:
                gene_fields, read_fields = gene_iter.next()
                read_groups.append(read_fields)
                if p_gene_fields is None:
                    p_gene_fields = gene_fields
                elif p_gene_fields != gene_fields:
                    raise RuntimeError('Parse error! gene fields %s and %s '
                                       'should match!', p_gene_fields,
                                       gene_fields)
            # gene_fields: (name, chrom, strand, txStart, txEnd, cdsStart,
            #                               cdsEnd, exons, name2, noncoding)
            expr_vals = dict(group_adjacent_reads(read_groups, max_dist))
            # sort locs s.t. those closest to the gene TSS come first
            locs = sorted(expr_vals.keys(), reverse=gene_fields[2] == '-')
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
                # for each combination of experiment, do a fisher's test
                for exp1, exp2 in itertools.combinations(
                                                    range(len(all_reads)), 2):
                    p = fisher.pvalue(expr_vals[ups_loc][exp1],
                                      expr_vals[downs_loc][exp1],
                                      expr_vals[ups_loc][exp2],
                                      expr_vals[downs_loc][exp2])
                    outline = '\t'.join(map(str, [gene_fields[0],
                                                ups_loc, downs_loc,
                                                compare_type,
                                                all_reads[exp1],
                                                expr_vals[ups_loc][exp1],
                                                expr_vals[downs_loc][exp1],
                                                all_reads[exp2],
                                                expr_vals[ups_loc][exp2],
                                                expr_vals[downs_loc][exp2],
                                                p.left_tail, p.right_tail,
                                                p.two_tail])) + '\n'
                    outfile.write(outline)
    except StopIteration:  # no more genes
        pass
    outfile.close()


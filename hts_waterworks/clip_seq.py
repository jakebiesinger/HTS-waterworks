
import itertools
import tempfile
import re
import random
import cPickle as pickle
from collections import defaultdict, Counter
import os
import urllib
import StringIO
import gzip
import subprocess
import shlex


from ipdb import set_trace as breakpoint
from matplotlib import pyplot
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import pybedtools
from ruffus import (transform, follows, merge, split, collate,
                    add_inputs, regex, suffix, jobs_limit, files)
from ruffus.task import active_if

from hts_waterworks.bootstrap import cfg, genome_path, get_genome
from hts_waterworks.utils.ruffus_utils import sys_call
from hts_waterworks.utils.common import reverseComplement, consensus_to_regex
from hts_waterworks.utils.motif_significance import zscore_normal
import hts_waterworks.mapping as mapping

import scipy as sp
import scipy.stats
import bx.bbi.bigwig_file

@active_if(cfg.getboolean('CLIP-seq', 'truncate_to_starts'))
#@transform(mapping.all_mappers_output, suffix('.mapped_reads'),
#           '.read_starts.mapped_reads', cfg.getint('CLIP-seq', 'truncate_to_starts_offset'))
@split(mapping.all_mappers_output, regex(r'(.*)\.mapped_reads'),
           [r'\1.plus.read_starts.mapped_reads', r'\1.minus.read_starts.mapped_reads'],
            cfg.getint('CLIP-seq', 'truncate_to_starts_offset'))
def truncate_to_starts(in_bed, out_bed, offset):
    """Truncate reads as only the start + offset base.  For use in CLIP-seq
    experiments
    """
    out_plus, out_minus = open(out_bed[0], 'w'), open(out_bed[1], 'w')

    for line in open(in_bed):
        chrom,start,stop,name,score,strand = line.strip().split('\t')
        start, stop = int(start), int(stop)
        assert strand in ['+', '-']
        if strand == '+':
            out_plus.write('\t'.join(map(str, [chrom, max(0, start + offset),
                        max(0, start + offset + 1), name, score, strand])) + '\n')
        else:
            out_minus.write('\t'.join(map(str, [chrom, max(0, stop - 1 - offset),
                        max(0, stop - 1 - offset + 1), name, score, strand])) + '\n')

@transform(truncate_to_starts, suffix('.mapped_reads'),
           '.pileup_reads')
def pileup_starts(in_bed, out_pileup):
    """Truncate reads as only the start + offset base.  For use in CLIP-seq
    experiments
    """
    cmd = 'bedItemOverlapCount %s -chromSize=%s.chrom.sizes %s > %s' % (
                            genome_path(), genome_path(),
                            in_bed, out_pileup)
    sys_call(cmd, file_log=False)

@active_if(cfg.getboolean('CLIP-seq', 'truncate_to_starts'))
@merge(pileup_starts, 'reproducible.pkl')
def reads_to_pickle(in_pileups, out_pkl):
    count_by_position = {}  # position -> file index -> count
    print in_pileups
    in_pileups = sorted(in_pileups)
    for findex, infile in enumerate(map(open, in_pileups)):
        strand = '+' if 'plus' in in_pileups[findex] else '-'
        for line in infile:
            chrom, start, stop, count = line.strip().split('\t')
            count = int(count)
            if (chrom, start, strand) not in count_by_position:
                count_by_position[(chrom, start, strand)] = sp.array([0] * (len(in_pileups)/2))
            assert count_by_position[(chrom, start, strand)][findex//2] == 0
            count_by_position[(chrom, start, strand)][findex//2] = count
    outfile = open(out_pkl, 'w')
    no_strands = list(set(map(lambda x: x.replace('.plus', '').replace('.minus',''), in_pileups)))
    pickle.dump(no_strands, outfile, protocol=-1)
    pickle.dump(count_by_position, outfile, protocol=-1)


@transform(reads_to_pickle, suffix('.pkl'), '.positions')
def reproducible_positions(in_pkl, out_data):
    """Assess the number of reads at each posiiton in each experiment.
    What percent of positions exactly X reads in both (all three?) experiments
    """
    infile = open(in_pkl)
    in_pileups = pickle.load(infile)
    count_by_position = pickle.load(infile)
    # make summary table
    with open(out_data, 'w') as outfile:
        outfile.write('{|\n|-\n| ')
        outfile.write(' || '.join(['min_count'] + ['at_least_%s_exp' % i for i in range(1,len(in_pileups)+1)] + ['total_exp1']) + '\n')
        for mincount in range(1,11):
            outfile.write('|-\n| ' + str(mincount) + ' || ')
            total_passing = [0] * len(in_pileups)
            for pos in count_by_position:
                count_passing = sum(count_by_position[pos] >= mincount)
                for i in range(count_passing):
                    total_passing[i] += 1
            outfile.write(' || '.join(map(str, total_passing)) + '\n|-\n|')
            #
            #total_passing = 0
            #total_passing_other_exp = [0] * (len(in_pileups)+1)
            #for pos in count_by_position:
            #    if count_by_position[pos][0] == mincount:
            #        total_passing += 1
            #        other_exp_passing = sum(count_by_position[pos] >= mincount)
            #        total_passing_other_exp[other_exp_passing] += 1
            #for exp_count in range(1, len(in_pileups)+1):
            #    sum_exp_count = sum(total_passing_other_exp[exp_count:])
            #    outfile.write(str(sum_exp_count) + ' || ')
            #outfile.write(str(total_passing) + '\n')
        #
        #for mincount in [10]:
        #    outfile.write('|-\n| ' + str(mincount) + '+ || ')
        #    total_passing = [0] * len(in_pileups)
        #    for pos in count_by_position:
        #        count_passing = sum(count_by_position[pos] >= mincount)
        #        for i in range(count_passing):
        #            total_passing[i] += 1
        #    outfile.write(' || '.join(map(str, total_passing)) + '\n|-\n|')
            #total_passing = 0
            #total_passing_other_exp = [0] * (len(in_pileups)+1)
            #for pos in count_by_position:
            #    if count_by_position[pos][0] >= mincount:
            #        total_passing += 1
            #        other_exp_passing = sum(count_by_position[pos] >= mincount)
            #        total_passing_other_exp[other_exp_passing] += 1
            #for exp_count in range(1, len(in_pileups)+1):
            #    sum_exp_count = sum(total_passing_other_exp[exp_count:])
            #    outfile.write(str(sum_exp_count) + ' || ')
            #outfile.write(str(total_passing) + '\n')
        outfile.write('|}\n')


@split(reads_to_pickle, regex(r'(.*)'), r'*.repro.peaks', r'%s.repro.peaks')
def reproducible_positions_to_peaks(in_pkl, _out, out_peaks_template):
    """Assess the number of reads at each position in each experiment.
    What percent of positions exactly X reads in both (all three?) experiments
    """
    infile = open(in_pkl)
    in_pileups = pickle.load(infile)
    count_by_position = pickle.load(infile)
    #count_by_position = {}  # position -> file index -> count

    # write out reads that pass min read criteria in all datasets
    combined_name = re.sub('64_CLIP\d+', '64_CLIP_all', in_pileups[0])
    for mincount in [1]:
        out_passing = open(out_peaks_template % (combined_name + '_%sminreads' % mincount), 'w')
        for pos in count_by_position:
            total_reads = sum(count_by_position[pos])
            if sum(count_by_position[pos] >= mincount) >= 3:
                out_passing.write('\t'.join(map(str, pos[:2] + (int(pos[1]) + 1, '.', total_reads, pos[2]))) + '\n')

@split(reads_to_pickle, regex(r'(.*)'), r'*.pileup_together.peaks', r'%s.pileup_together.peaks')
def all_positions_pileup_together(in_pkl, _out, out_peaks_template):
    """Write out all peak positions, taking the summation across experiments
    """
    infile = open(in_pkl)
    in_pileups = pickle.load(infile)
    count_by_position = pickle.load(infile)
    #count_by_position = {}  # position -> file index -> count

    # write out reads that pass min read criteria in all datasets
    combined_name = re.sub('64_CLIP\d+', '64_CLIP_all', in_pileups[0])
    out_passing = open(out_peaks_template % (combined_name + '.pileup_together.peaks'), 'w')
    for pos in count_by_position:
        total_reads = sum(count_by_position[pos])
        out_passing.write('\t'.join(map(str, pos[:2] + (int(pos[1]) + 1, '.', total_reads, pos[2]))) + '\n')


@split(reads_to_pickle, regex(r'(.*)'), r'*.all.peaks', r'%s.all.peaks')
def all_positions_to_peaks(in_pkl, _out, out_peaks_template):
    """Output the pickle file as a bed pileup
    """
    for f in _out:
        os.unlink(f)
    infile = open(in_pkl)
    in_pileups = pickle.load(infile)
    count_by_position = pickle.load(infile)
    #count_by_position = {}  # position -> file index -> count


    for findex in range(len(in_pileups)):
        combined_name = in_pileups[findex].replace('.plus', '').replace('.minus','')
        out_file = open(out_peaks_template % (combined_name + '_%sminreads' % 0), 'a')
        for pos in count_by_position:
            total_reads = count_by_position[pos][findex]
            if total_reads > 0:
                out_file.write('\t'.join(map(str, pos[:2] + (int(pos[1]) + 1, '.', total_reads, pos[2]))) + '\n')



@transform(reproducible_positions_to_peaks, suffix('.peaks'), '.near_polyA_sites.peaks', 'hela.hg19.no_internal_priming.pas_seq_sites')
def find_peaks_near_polyA_sites(in_pileup, out_pileup, in_polya_sites):
    """Filter out peaks that are not near (80nt) polyA sites"""
    # DONE: liftover polyA sites from hg18
    # DONE: input polyA sites, get distance relative to polyA sites, filter out those not <80bp downstream

    polya_sites = pybedtools.BedTool(in_polya_sites)
    pileup = pybedtools.BedTool(in_pileup)

    # get sites that are downstream of polyA sites, on the same strand
    close_pileup_with_polya = pileup.closest(polya_sites, D='b', t='first', s=True, iu=True)
    close_pileup_with_polya = close_pileup_with_polya.filter(lambda a: 0 <= int(a[-1]) <= 80)
    close_pileup = close_pileup_with_polya.cut([0,1,2,3,4,5])
    close_pileup.saveas(out_pileup)



@transform('hela.hg19.no_internal_priming.pas_seq_sites', suffix('.pas_seq_sites'),
           '.random.pas_seq_sites')
def random_polya_sites(in_polya, out_random_sites):
    posns = pybedtools.BedTool(in_polya)
    genes = pybedtools.BedTool('%s.refseq_genes.all' % genome_path())
    posns_with_genes = intersect_same_strand_keep_shortest(posns, genes)
    
    def shuffle_within_genes(features):
        gene_start, gene_end = int(features[7]), int(features[8])
        features.start = random.randrange(gene_start, gene_end)
        features.stop = features.start + 1
        return features
    rand_posns = [posns_with_genes.each(shuffle_within_genes).cut(range(6)) for i in range(100)]
    rand_bed = '\n'.join(map(str, rand_posns))
    with open(out_random_sites, 'w') as outfile:
        outfile.write(rand_bed)



@split(['hela.hg19.no_internal_priming.pas_seq_sites', random_polya_sites],
    regex(r'(.*)\.pas_seq_sites'),
       add_inputs(reproducible_positions_to_peaks),
       [r'\1.clip_dependent.pas_seq.peaks',
        r'\1.clip_independent.pas_seq.peaks'])
def find_clip_dep_and_indep_polya_sites(in_files, out_files):
    """Filter out peaks that are near (80nt) polyA sites"""
    def make_index(features, _index=[0]):
        features.name = str(_index[0])
        _index[0] += 1
        return features
    polya_sites = pybedtools.BedTool(in_files[0]).filter(lambda f: int(f.score) > 20).each(make_index).saveas()
    print 'there are %s valid sites in %s' % (len(polya_sites), len(open(in_files[0]).readlines()))
    clip_pileup = pybedtools.BedTool(in_files[1])
    out_dependent = out_files[0]
    out_independent = out_files[1]

    # slop polya downstream by 80nt, then intersect with clip sites, only considering
    # sites that are on the same strand. collapse the read counts at the clip sites
    polya_with_clip_reads = polya_sites.slop(s=True, r=80, l=0,
                                       g='%s.chrom.sizes' % genome_path())\
                .intersect(clip_pileup, wao=True, s=True)\
                .groupby(g=[1,2,3,4,5,6], c=11, ops=['collapse']).saveas()
    # count the number of clip sites we intersected with.  If 5+, mark site as
    # dependent. Otherwise, mark as independent.
    #import ipdb; ipdb.set_trace()
    def has_enough_tall_reads(features):
        reads = features[-1].split(',')
        return len([r for r in reads if int(r) > 2]) >= 1
    clip_dep_sites = set()
    for p in polya_with_clip_reads.filter(has_enough_tall_reads):
        clip_dep_sites.add(p.name)

    # finally filter the original list based on dependence or independence
    dep_sites = polya_sites.filter(lambda f: f.name in clip_dep_sites).saveas(out_dependent)
    indep_sites = polya_sites.filter(lambda f: f.name not in clip_dep_sites).saveas(out_independent)




@split(['hela.hg19.no_internal_priming.pas_seq_sites'],
    regex(r'(.*)\.pas_seq_sites'),
       add_inputs(reproducible_positions_to_peaks),
       #[r'\1.clip_ratios.pas_seq.peaks'])
       [r'\1.clip_ratios.top1000.pas_seq.peaks',
        r'\1.clip_ratios.bottom1000.pas_seq.peaks'])
def find_clip_dep_and_indep_polya_sites_by_ratio(in_files, out_files):
    """Filter out peaks that are near (80nt) polyA sites"""
    def make_index(features, _index=[0]):
        features.name = str(_index[0])
        _index[0] += 1
        return features
    polya_sites = pybedtools.BedTool(in_files[0]).filter(lambda f: int(f.score) > 20).each(make_index).saveas()
    print 'there are %s valid sites in %s' % (len(polya_sites), len(open(in_files[0]).readlines()))
    clip_pileup = pybedtools.BedTool(in_files[1])

    # slop polya downstream by 80nt, then intersect with clip sites, only considering
    # sites that are on the same strand. collapse the read counts at the clip sites
    polya_with_clip_reads = polya_sites.slop(s=True, r=30, l=0,
                                       g='%s.chrom.sizes' % genome_path())\
                .intersect(clip_pileup, wao=True, s=True)\
                .groupby(g=[1,2,3,4,5,6], c=11, ops=['collapse']).saveas()
    utr3_region = pybedtools.BedTool('hg19.refseq_genes.utr3')
    print '20+ polya %s, 20+ polya in 3utr %s, 20+ polya in 3utr with clip %s' % (len(polya_sites),
                                                    len(polya_sites.intersect(utr3_region, s=True, u=True)),
                                                    len(polya_with_clip_reads.intersect(utr3_region, s=True, u=True)))
    # count the number of clip sites we intersected with.  If 5+, mark site as
    # dependent. Otherwise, mark as independent.
    #import ipdb; ipdb.set_trace()
    def get_ratio_pas_to_clip_reads(features):
        clip_reads = sum(int(r) if r != '-1' else 0 for r in features[-1].split(','))
        ratio = float(features.score) / clip_reads if clip_reads > 0 else float('inf')
        features.append(str(ratio))
        return features
    print 'there are %s PAS sites without a CLIP site' % len(polya_with_clip_reads.filter(lambda f: f[-1] == '-1'))
    with_ratio = polya_with_clip_reads.each(get_ratio_pas_to_clip_reads).saveas()
    sorted_ratio = pybedtools.BedTool('\n'.join(sorted(map(str, with_ratio), key=lambda f: (float(f.strip().split('\t')[-1]), float(f.strip().split('\t')[4])))), from_string=True)
    sorted_final = sorted_ratio.slop(s=True, r=-30,l=0,g='%s.chrom.sizes' % genome_path()).cut(range(6)).saveas()
    print 'final sites: %s' % len(sorted_final)
    final_sites = list(map(str, sorted_final))
    #final_sites = list(map(str, sorted_ratio))
     
    with open(out_files[0], 'w') as outfile:
        out_str = '\n'.join(map(str, final_sites[-1001:]))
        print '%s in 3utr: %s' % (out_files[0], len(pybedtools.BedTool(out_str, from_string=True).intersect(utr3_region, s=True, u=True)))
        outfile.write(out_str)
    with open(out_files[1], 'w') as outfile:
        out_str = '\n'.join(map(str, final_sites[:1001]))
        print '%s in 3utr: %s' % (out_files[1], len(pybedtools.BedTool(out_str, from_string=True).intersect(utr3_region, s=True, u=True)))
        outfile.write(out_str)
    


@merge(find_clip_dep_and_indep_polya_sites_by_ratio, 'pas_seq_to_clip_ratio_hist.png')
def plot_ratio_hists(in_reads, out_png):
    ratios = []
    for findex, filein in enumerate(in_reads):
        ratios.append([])
        for l in open(filein):
            value = float(l.strip().split('\t')[-1])
            if value == 0:
                value = 1
            ratios[findex].append(sp.log2(value))

    hist_range = min(min(r) for r in ratios), max(max(r) for r in ratios)
    hists = [sp.histogram(ratios[i], range=hist_range, bins=50)[0] / float(len(ratios[i])) for i in range(len(in_reads))]
    lines = []
    xs = sp.linspace(hist_range[0], hist_range[1], 50, endpoint=False)
    pyplot.figure()
    for i, h in enumerate(hists):
        lines.append(pyplot.plot(xs, h, label=in_reads[i]))
    pyplot.legend(lines, ['random' if 'random' in f else 'PAS sites' for f in in_reads])
    pyplot.xlabel('$\log_2$(PAS / CLIP) in 80nt downstream')
    pyplot.ylabel('fraction of all sites')
    pyplot.savefig(out_png)


@transform(find_clip_dep_and_indep_polya_sites, suffix('.peaks'), '.peaks.fasta')
def get_80nt_downstream_pas_seq(in_polya, out_fasta):
    """snag the 80bp downstream of the polya sites"""
    genome = get_genome(None, None, False)
    with open(out_fasta, 'w') as out_fasta:
        for i, line in enumerate(sorted(open(in_polya), key=lambda l: int(l.strip().split('\t')[4]), reverse=True)):
            #chrom, start, stop, count = line.strip().split('\t')
            chrom, start, stop, name, count, strand = line.strip().split('\t')[:6]
            if strand == '+':
                seq = genome[chrom][int(start):int(start)+80+1]
            else:
                seq = -genome[chrom][int(start)-80:int(start)+1]
            seq = str(seq).upper().replace('T', 'U')
            out_fasta.write('>seq_%s_%s_reads\n%s\n' % (i, count, seq))
            if i >= 10000:
                break


@transform(get_80nt_downstream_pas_seq, suffix('.fasta'), '.fasta.png')
def make_weblogo_downstream_pas_seq(in_fasta, out_png):
    """Generate a weblogo of the alignment of sequences around a clip binding site"""
    cmd = "weblogo -f {in_fasta} -F png_print -o {out_png} --fineprint '' --errorbars NO -Y NO --composition 'H. sapiens' -s large --scale-width NO -U probability"
    cmd = cmd.format(**locals())
    sys_call(cmd)


@transform(get_80nt_downstream_pas_seq, suffix('.fasta'), '.meme.discovered.motifs')
def meme_motif_search_downstream_pas_seq(in_fasta, out_motifs):
    import hts_waterworks.motif_discovery as motif_discovery
    motif_discovery.discover_meme_motifs(in_fasta, out_motifs)


# 2. Are AWTAAA overrepresented within 100nt upstream of intronic
#    CstF64 binding sites (read count threshold 1-10)?
#    Z score of AWTAAA comparing 100nt upstream of intronic CstF64 binding sites
#    and same-sized fragment randomly selected from the introns of the same gene.

def check_intronic_consensus_enrichment(in_pileup, in_consensus, in_introns, out_summary):
    # ID all positions
    def make_index(features, _index=[0]):
        features.name = str(_index[0])
        _index[0] += 1
        return features
    clip_sites = pybedtools.BedTool(in_pileup).each(make_index)

    # Find clip binding sites in introns of genes.
    introns = pybedtools.BedTool(in_introns)
    clip_in_introns = clip_sites.intersect(introns, wo=True, s=True)

    #   Count each site once, select the shortest transcript that it overlaps with
    # TODO: select only the shortest transcript I overlap with
    # TODO: map clip ID to what gene it overlaps with

    # check 100nt upstream of these regions for AWTAAA, get rate of occurrence
    cons_sites = pybedtools.BedTool(in_consensus)
    clip_near_cons = clip_in_introns.cut([0,1,2,3,4,5]).closest(cons_sites, D='a', id=True, s=True)
    clip_near_cons = clip_near_cons.filter(lambda p: -100 <= int(p[-1]) <= 0)

    fg_rate = len(clip_in_introns) / float(clip_near_cons)

    for p in clip_near_cons:
        pass
        # for each intronic binding site,
        #    shuffle the site around the introns of the current gene 100 times
        # TODO: get the complete collection of introns and sample from them for rand location
        #for i in range(100):
            # create a single bed file with random positions in it
        #    see how many times the AWTAAA is present
        # TODO: apply clip_near_cons filter to these random locations, add to total occurrence

    # calculate Zscore



# TODO: plot of nearest polyA sites (<80nt downstream of polyA sites)


# TODO: take polyA sites and find those that don't have CstF64 binding site nearby
# Calculate the % of polyA sites without CstF64 binding sites: Take all PAS-seq peaks with at least 20 reads, look at the 80nt downstream of the designated polyA sites, if less than 5 peaks (each with 2 reads or less) of CstF64 iCLIP tags are found, this is counted as a CstF64-independent polyA sites. All other sites are CstF64-dependent sites.



# Zscore/hexamer of poly-A downstream sites.
##  2. Sequence comparison between CstF64-dependent and -independent sites. Take 0-80nt downstream of polyA sites and perform MEME or Z score of hexamers.


@transform('hg19.refseq_genes.noncoding', suffix('.noncoding'),
           '.noncoding_with_wgRna_tRNA_lincRNA')
def get_noncoding_rna_tracks(in_refseq_noncoding, out_all_noncoding):
    """combine refseq noncoding with downloaded and combined wgRNA, lincRNA, and tRNA
    """

    with open(out_all_noncoding, 'w') as outfile:
        outfile.writelines(open(in_refseq_noncoding))
        data = urllib.urlopen('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgRna.txt.gz').read()
        data = gzip.GzipFile(fileobj=StringIO.StringIO(data))
        for line in data:
            f = line.strip().split('\t')
            outfile.write('\t'.join([f[1],f[2],f[3],f[4],f[9],f[6]]) + '\n')
        data = urllib.urlopen('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/lincRNAsTranscripts.txt.gz').read()
        data = gzip.GzipFile(fileobj=StringIO.StringIO(data))
        for line in data:
            f = line.strip().split('\t')
            outfile.write('\t'.join([f[2],f[4],f[5],f[1],f[0],f[3]]) + '\n')
        data = urllib.urlopen('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/tRNAs.txt.gz').read()
        data = gzip.GzipFile(fileobj=StringIO.StringIO(data))
        for line in data:
            f = line.strip().split('\t')
            outfile.write('\t'.join([f[1],f[2],f[3],f[4],f[5],f[6]]) + '\n')


@follows(get_noncoding_rna_tracks)
@transform(reproducible_positions_to_peaks, suffix('.peaks'), '.nogenes.peaks_in_noncoding.peaks')
def get_sites_nongenic_in_noncoding(in_pileup, out_noncoding_pileup):
    """Find pileup sites that *aren't* in refseq genes but are in other noncoding genes.
    """
    pileup = pybedtools.BedTool(in_pileup)
    genes = pybedtools.BedTool('hg19.refseq_genes.all')
    coding_genes = genes.filter(lambda f: not f[3].startswith('NR_'))
    # remove coding genes
    pileup_no_coding = pileup.intersect(coding_genes, s=True, v=True, wa=True)
    # get those that overlap with noncoding genes
    noncoding = pybedtools.BedTool('hg19.refseq_genes.noncoding_with_wgRna_tRNA_lincRNA')
    pileup_no_coding.intersect(noncoding, s=True, u=True).saveas(out_noncoding_pileup)


@follows(get_noncoding_rna_tracks)
@transform(reproducible_positions_to_peaks, suffix('.peaks'), '.nogenes.peaks_in_noncoding.table')
def nongenic_in_noncoding_table(in_pileup, out_noncoding_pileup):
    """Find pileup sites that *aren't* in refseq genes but are in other noncoding genes.
    Make a table of the number of counts
    """
    # id each pileup read
    def make_index(features, _index=[0]):
        features.name = str(_index[0])
        _index[0] += 1
        return features
    pileup = pybedtools.BedTool(in_pileup).each(make_index).saveas()
    genes = pybedtools.BedTool('hg19.refseq_genes.all')
    coding_genes = genes.filter(lambda f: not f[3].startswith('NR_'))
    # remove coding genes
    pileup_no_coding = pileup.intersect(coding_genes, s=True, v=True, wa=True)
    # get those that overlap with noncoding genes
    noncoding = pybedtools.BedTool('hg19.refseq_genes.noncoding_with_wgRna_tRNA_lincRNA')
    noncoding_genes_with_pileup = noncoding.intersect(pileup_no_coding, s=True, wa=True, wb=True)
    # aggregate the number of reads and count the number of sites that overlap
    reads_per_nc_gene = noncoding_genes_with_pileup.groupby(g=[1,2,3,4,5,6], c='10,11', ops='distinct,sum')
    def distinct_id_to_num_sites(features):
        features[6] = str(len(features[6].split(',')))
        return features
    reads_per_nc_gene.each(distinct_id_to_num_sites).saveas(out_noncoding_pileup)
    
    



@merge(all_positions_to_peaks, 'reproducible.motifs', 5)
#@merge(reproducible_positions_to_peaks, ['reproducible.1minreads.motifs', 'reproducible.10minreads.motifs'], 5)
def reproducible_motifs(in_pileups, out_mer_data, mer_length):
    """Examine the occurrence of pentamers overlapping pileup sites"""

    for index, fileset in enumerate([sorted(in_pileups)]):
        # track the total number of sites in each experiment
        total_sites = [0] * len(fileset)
        # track how many sites each kmer is present in
        # num_sites_with_kmer[file][kmer] -> num_sites
        num_sites_with_kmer = [Counter() for i in range(len(fileset))]
        genome = get_genome(None, None, False)
        for findex, infile in enumerate(map(open, fileset)):
            for lindex, line in enumerate(infile):
                total_sites[findex] += 1
                chrom, start, stop, name, count, strand = line.strip().split('\t')
                cur_seqs = set()
                # get all kmers overlapping the binding site
                for slice_start in range(-mer_length+1, 1):
                    seq = genome[chrom][int(start) + slice_start : int(start) + slice_start + mer_length]
                    assert strand in ['+', '-']
                    if strand == '+':
                        seq = str(seq).upper()
                    else:
                        seq = str(-seq).upper()
                    cur_seqs.add(seq)
                for seq in cur_seqs:
                    num_sites_with_kmer[findex][seq] += 1

        # sort the kmers by their occurence rate in expermient 0, largest first
        all_kmers = set()
        for kmers in num_sites_with_kmer:
            all_kmers.update(kmers.iterkeys())
        all_kmers = sorted(all_kmers, key=lambda k: num_sites_with_kmer[0][k], reverse=True)

        # reformat kmer frequencies to array and plot a 3d array R2 values
        kmer_freqs = sp.array([[num_sites_with_kmer[f][k] / float(total_sites[f]) for f in range(len(fileset))] for k in all_kmers])
        sp.save(out_mer_data + '.npy', kmer_freqs)

        # make a 3d plot of first 3 correlations and experiment kmer freqs
        oldval, matplotlib.rcParams['font.size'] = matplotlib.rcParams['font.size'], 10
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d')
        zeros = sp.zeros_like(kmer_freqs[:,0])
        ax.scatter(kmer_freqs[:,0], kmer_freqs[:,1], zeros, color='r')
        ax.scatter(kmer_freqs[:,0], zeros, kmer_freqs[:,2], color='g')
        ax.scatter(zeros, kmer_freqs[:,1], kmer_freqs[:,2], color='b')
        ax.scatter(kmer_freqs[:,0], kmer_freqs[:,1], kmer_freqs[:,2], color='k')
        ax.set_xlabel('\n'+fileset[0].split('.')[0], linespacing=2)
        ax.set_ylabel('\n'+fileset[1].split('.')[0], linespacing=2)
        ax.set_zlabel('\n'+fileset[2].split('.')[0], linespacing=2)
        ax.view_init(elev=10, azim=35)
        pyplot.savefig(out_mer_data + '.3d_scatter.png', dpi=240)
        ax.view_init(elev=40, azim=72)
        pyplot.savefig(out_mer_data + '.3d_scatter2.png', dpi=240)
        matplotlib.rcParams['font.size'] = oldval

        # make a wiki table of kmers
        with open(out_mer_data, 'w') as outfile:
            # table for R^2 and pvalues
            outfile.write('{|\n| comparison || r^2 || pvalue ')
            for lhs, rhs in itertools.combinations(range(len(fileset)), 2):
                lhs_name, rhs_name = fileset[lhs], fileset[rhs]
                r2, pval = scipy.stats.pearsonr(kmer_freqs[:,lhs], kmer_freqs[:,rhs])
                outfile.write('{lhs_name} vs {rhs_name} || {r2} || {pval} \n|-\n'.format(**locals()))
            outfile.write('|}\n')

            outfile.write('{|\n| ')
            outfile.write(' || '.join('kmer count in %s' % f for f in fileset) + '\n')
            for kmer in all_kmers:
                outfile.write('|-\n| ' + kmer + '|| ' + ' || '.join(map(str, [num_sites_with_kmer[i][kmer] / float(total_sites[i]) for i in range(len(fileset))])) + '\n')
            outfile.write('|-\n| total sites || ' + ' || '.join(map(str, total_sites)))
            outfile.write('|}\n')





#@transform([reproducible_positions_to_peaks, find_peaks_near_polyA_sites,
#            all_positions_pileup_together, find_clip_dep_and_indep_polya_sites,
#            find_clip_dep_and_indep_polya_sites_by_ratio],
@transform([find_clip_dep_and_indep_polya_sites_by_ratio],
    suffix('.peaks'), '.peaks_in_genes')
def reads_within_same_gene(in_pileup, out_sites_with_gene):
    """For each read in pileup, find the transcript that it overlaps with.
    When a read overlaps with multiple transcripts, keep only the shortest one.
    """
    posns = pybedtools.BedTool(in_pileup)
    genes = pybedtools.BedTool('%s.refseq_genes.all' % genome_path())
    posns_with_genes = intersect_same_strand_keep_shortest(posns, genes)
    posns_with_genes.saveas(out_sites_with_gene)



@transform(reads_within_same_gene, suffix(''), '.fdr_table')
def calculate_read_gene_fdr(in_pileup_with_gene, out_fdr_table):
    """Given reads and the genes they overlap with, calculate the FDR by
    reshuffling the reads around the gene uniformly, then taking the ratio
          (bg heights >= x) / (fg heights >= x)
    """
    # calculate how many reads are in each transcript
    print '# sorting by gene position'
    # sort by gene position
    def sort_by_gene(line):
        f = line.strip().split('\t')
        return [f[6]] + map(int, f[7:9])
    sortedlines = sorted(open(in_pileup_with_gene), key=sort_by_gene)
    posns_in_genes = pybedtools.BedTool('\n'.join(sortedlines), from_string=True)

    # group by gene position, count number of reads in that position
    reads_per_gene = posns_in_genes.groupby(g=[7,8,9], c=5, ops='sum')
    rand_heights = []
    with open(out_fdr_table + '.diag', 'w') as outfile:
        for p in reads_per_gene:
            outfile.write('{readcount}\t{gene_length}\n'.format(readcount=p[-1], gene_length=int(p[2])-int(p[1])))
    all_heights, all_lengths = zip(* [map(int, line.strip().split('\t')) for line in open(out_fdr_table + '.diag')])
    biggest_h = max(all_heights)
    print biggest_h
    biggest_l = max(all_lengths)
    print biggest_l

    pyplot.figure()
    pyplot.plot(all_heights, all_lengths, 'o')
    expected_height = sp.logspace(1,sp.log10(biggest_h))
    expected_length = sp.logspace(1,sp.log10(biggest_l))
    for i in [1, 10, 100]:
        pyplot.plot(expected_height * i, expected_length, '.')
    pyplot.xlabel("Reads within gene")
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.ylabel("Gene Lengths")
    pyplot.title("Reads per gene for FDR calculation")
    pyplot.savefig(out_fdr_table + '.diag.png', dpi=240)

    rand_zeros = 0
    for gene in reads_per_gene:
        read_count = int(gene[-1])
        gene_length = int(gene[2]) - int(gene[1])
        for i in range(100):
            # count how many times each position is sampled, keep the heights
            positions = sp.floor(sp.round_(sp.rand(read_count) * (gene_length - 1)))
            heights = Counter(positions).values()
            rand_heights.extend(heights)
            # don't forget to count all those heights of 0
            if len(heights) < gene_length:
                rand_zeros += gene_length - len(heights)

    print '# getting fg heights'
    # get fg height distribution
    fg_heights = [int(p.score) for p in posns_in_genes]

    print '# calculatin fdr'
    # calculate FDR
    fg_cumsum, bg_cumsum = make_cumsum_dicts(fg_heights, rand_heights, rand_zeros)
    with open(out_fdr_table, 'w') as outfile:
        outfile.write('x reads\tbg>=x\tfg>=x\tbg_rate\tfg_rate\tfdr\n')
        for readcount in sorted(fg_cumsum):
            print readcount
            bg_passing = bg_cumsum[readcount]
            fg_passing = fg_cumsum[readcount]
            bg_rate = float(bg_passing) / (len(rand_heights) + rand_zeros)
            fg_rate = float(fg_passing) / len(fg_heights)
            if bg_rate == 0 and fg_rate == 0:
                fdr = float('nan')
            elif fg_rate == 0:
                fdr = float('inf')
            elif bg_rate == 0:
                fdr = 0.
            else:
                fdr = bg_rate / fg_rate
            outfile.write('{readcount}\t{bg_passing}\t{fg_passing}\t{bg_rate:.3f}\t{fg_rate:.3f}\t{fdr:.8f}\n'.format(**locals()))



@merge([reproducible_positions_to_peaks, pileup_starts], 'reproducilbility.fdr_table')
def calculate_read_gene_reproducibility_fdr(in_all_pileups, out_fdr_table):
    """Given reads and the genes they overlap with, calculate the FDR by
    reshuffling the reads around the gene uniformly, then taking the ratio
          (bg heights >= x) / (fg heights >= x)
    with a reproducibility minimum threshold
    """
    read_count_in_genes = {}
    def make_index(features, _index=[0]):
        features.append(str(_index[0]))
        _index[0] += 1
        return features
    genes = pybedtools.BedTool('%s.refseq_genes.all' % genome_path()).each(make_index)    
    id_col = 6 + len(genes.__iter__().next())
    for index, in_pileup in enumerate(in_all_pileups):
        print in_pileup
        posns = pybedtools.BedTool(in_pileup)
        posns_in_genes = intersect_same_strand_keep_shortest(posns, genes)
        # sort by gene id
        sortedlines = sorted(open(posns_in_genes.fn), key=lambda f: f.strip().split('\t')[-1])
        posns_in_genes = pybedtools.BedTool('\n'.join(sortedlines), from_string=True)

        # calculate how many reads are in current transcript from current dataset
        reads_per_gene = posns_in_genes.groupby(g=[id_col], c=5, ops='sum')
        #read_count_in_genes[]
        

    # create a new random pileup for each experiment in the current transcript
    # get reproducible heights, see how many pass threshold in all three experiments
    




    
    
    print '# sorting by gene position'
    

    # group by gene position, count number of reads in that position
    reads_per_gene = posns_in_genes.groupby(g=[7,8,9], c=5, ops='sum')
    rand_heights = []
    with open(out_fdr_table + '.diag', 'w') as outfile:
        for p in reads_per_gene:
            outfile.write('{readcount}\t{gene_length}\n'.format(readcount=p[-1], gene_length=int(p[2])-int(p[1])))
    all_heights, all_lengths = zip(* [map(int, line.strip().split('\t')) for line in open(out_fdr_table + '.diag')])
    biggest_h = max(all_heights)
    print biggest_h
    biggest_l = max(all_lengths)
    print biggest_l

    pyplot.figure()
    pyplot.plot(all_heights, all_lengths, 'o')
    expected_height = sp.logspace(1,sp.log10(biggest_h))
    expected_length = sp.logspace(1,sp.log10(biggest_l))
    for i in [1, 10, 100]:
        pyplot.plot(expected_height * i, expected_length, '.')
    pyplot.xlabel("Reads within gene")
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.ylabel("Gene Lengths")
    pyplot.title("Reads per gene for FDR calculation")
    pyplot.savefig(out_fdr_table + '.diag.png', dpi=240)

    rand_zeros = 0
    for gene in reads_per_gene:
        read_count = int(gene[-1])
        gene_length = int(gene[2]) - int(gene[1])
        for i in range(100):
            # count how many times each position is sampled, keep the heights
            positions = sp.floor(sp.round_(sp.rand(read_count) * (gene_length - 1)))
            heights = Counter(positions).values()
            rand_heights.extend(heights)
            # don't forget to count all those heights of 0
            if len(heights) < gene_length:
                rand_zeros += gene_length - len(heights)

    print '# getting fg heights'
    # get fg height distribution
    fg_heights = [int(p.score) for p in posns_in_genes]

    print '# calculatin fdr'
    # calculate FDR
    fg_cumsum, bg_cumsum = make_cumsum_dicts(fg_heights, rand_heights, rand_zeros)
    with open(out_fdr_table, 'w') as outfile:
        outfile.write('x reads\tbg>=x\tfg>=x\tbg_rate\tfg_rate\tfdr\n')
        for readcount in sorted(fg_cumsum):
            print readcount
            bg_passing = bg_cumsum[readcount]
            fg_passing = fg_cumsum[readcount]
            bg_rate = float(bg_passing) / (len(rand_heights) + rand_zeros)
            fg_rate = float(fg_passing) / len(fg_heights)
            if bg_rate == 0 and fg_rate == 0:
                fdr = float('nan')
            elif fg_rate == 0:
                fdr = float('inf')
            elif bg_rate == 0:
                fdr = 0.
            else:
                fdr = bg_rate / fg_rate
            outfile.write('{readcount}\t{bg_passing}\t{fg_passing}\t{bg_rate:.3f}\t{fg_rate:.3f}\t{fdr:.8f}\n'.format(**locals()))






def make_cumsum_dicts(fg_counts, bg_counts, bg_zeros):
    """Creates two dictionaries with the cumulative sum of counts >= key
    """
    fg_cumsum, bg_cumsum = {}, {}
    items = set(fg_counts).union(bg_counts)
    fg_above, bg_above = 0, 0
    fg_index, bg_index = len(fg_counts) - 1, len(bg_counts) - 1
    fg_counts = sorted(fg_counts)
    bg_counts = sorted(bg_counts)
    # for each level, keep track of how many in fg and bg pass
    for level in sorted(items, reverse=True):
        while fg_counts[fg_index] >= level:
            if fg_index == 0:
                break
            fg_index -= 1
            fg_above += 1
        while bg_counts[bg_index] >= level:
            if bg_index == 0:
                break
            bg_index -= 1
            bg_above += 1
        fg_cumsum[level] = fg_above
        bg_cumsum[level] = bg_above
    fg_cumsum[0] = fg_above
    bg_cumsum[0] = bg_zeros
    return fg_cumsum, bg_cumsum





@active_if(False)
@transform(reads_within_same_gene,
    regex(r'(.*(?<!\.pileup_together))\.peaks_in_genes'), r'\1.peaks_in_genes.6mer.vs_same_gene', 6, 10)
def kmer_zscore_compare_within_same_gene(in_pileup_with_gene, out_zscores, kmer_length, windowsize):
    """count kmers in the [-windowsize : windowsize+1] region around read positions.
    Keep only the reads falling within genes. For each kmer, sample the gene region intersected
    by the current read to generate a background occurrence rate, then calculate the z-score for the kmer
    vs this local current-gene background.
    """
    samples_per_overlap = 100
    genome = get_genome(None, None, None)
    fg_mers = defaultdict(set)  # kmer -> which overlaps contain it
    num_fg_regions = 0
    bg_mers = defaultdict(set)  # kmer -> which overlaps contain it
    num_posns_fields = 6    
    
    posns_in_genes = pybedtools.BedTool(in_pileup_with_gene)

    for real_index, overlap in enumerate(posns_in_genes):
        assert overlap.strand in ['+', '-']
        num_fg_regions += 1
        # get the kmers in this region
        seq = genome[overlap.chrom][overlap.start - windowsize : overlap.start + windowsize + 1]
        seq = str(seq) if overlap.strand == '+' else str(-seq)
        seq = seq.lower()
        kmers = [seq[start:start+kmer_length] for start in range(len(seq)-kmer_length+1)]
        # count the real kmers
        for kmer in kmers:
            fg_mers[kmer].add(real_index)

        # get 100 random regions in the same gene, get their kmer counts
        gene_start = int(overlap[num_posns_fields + 1])
        gene_end = int(overlap[num_posns_fields + 2])
        rand_starts = [random.randrange(gene_start, gene_end) for i in range(samples_per_overlap)]
        for random_index, rand_start in enumerate(rand_starts):
            rand_seq = genome[overlap.chrom][rand_start - windowsize : rand_start + windowsize + 1]
            rand_seq = str(rand_seq) if overlap.strand == '+' else str(-rand_seq)
            rand_seq = rand_seq.lower()
            rand_kmers = [rand_seq[start:start+kmer_length] for start in range(len(rand_seq)-kmer_length+1)]
            # count the bg kmers
            for kmer in rand_kmers:
                bg_mers[kmer].add('%s_%s' % (real_index,random_index))
    print 'done with overlap'

    with open(out_zscores, 'w') as outfile:
        # for each kmer we saw in fg, calculate zscore-- (x - mu) / sigma
        for kmer in fg_mers:
            z = zscore_normal(len(fg_mers[kmer]), num_fg_regions, len(bg_mers[kmer]), num_fg_regions * samples_per_overlap)
            fg_rate = len(fg_mers[kmer]) / float(num_fg_regions)
            bg_rate = len(bg_mers[kmer]) / float(num_fg_regions * samples_per_overlap)
            print >>outfile, '%s\t%s\t%s\t%.5f\t%s\t%s\t%.5f\t%.5f' % (kmer, len(fg_mers[kmer]), num_fg_regions, fg_rate, len(bg_mers[kmer]), num_fg_regions * samples_per_overlap, bg_rate, z)


@transform(reads_within_same_gene,
    #regex(r'(.*(?<=dependent)\.pas_seq)\.peaks_in_genes'), r'\1.peaks_in_genes.80nt_down.6mer.vs_same_gene', 6, 80)
    regex(r'(.*(?<=1000)\.pas_seq)\.peaks_in_genes'), r'\1.peaks_in_genes.30nt_down.6mer.vs_same_gene', 6, 30)
def kmer_zscore_compare_within_same_gene_pas_downstream(in_pileup_with_gene, out_zscores, kmer_length, windowsize):
    """count kmers in the 80nt downstream of the pas-
    Keep only the reads falling within genes. For each kmer, sample the gene region intersected
    by the current read to generate a background occurrence rate, then calculate the z-score for the kmer
    vs this local current-gene background.
    """
    samples_per_overlap = 100
    genome = get_genome(None, None, None)
    fg_mers = defaultdict(set)  # kmer -> which overlaps contain it
    num_fg_regions = 0
    bg_mers = defaultdict(set)  # kmer -> which overlaps contain it
    num_posns_fields = 6
    
    filenames = [in_pileup_with_gene.replace('clip_independent', 'clip_dependent'), in_pileup_with_gene.replace('clip_dependent', 'clip_independent')]
    most_sites_allowed = min(len(pybedtools.BedTool(filenames[0])), len(pybedtools.BedTool(filenames[1])))
    
    #import ipdb; ipdb.set_trace()
    best_sites = sorted(open(in_pileup_with_gene), key=lambda l: int(l.strip().split('\t')[4]), reverse=True)
    best_sites = best_sites[:most_sites_allowed]
    posns_in_genes = pybedtools.BedTool('\n'.join(map(str, best_sites)), from_string=True)
    

    for real_index, overlap in enumerate(posns_in_genes):
        assert overlap.strand in ['+', '-']
        num_fg_regions += 1
        # get the kmers in this region
        if overlap.strand == '+':
            seq = str(genome[overlap.chrom][overlap.start : overlap.start + windowsize + 1])
        else:
            seq = str(- genome[overlap.chrom][max(0, overlap.start - windowsize) : overlap.start + 1])
        seq = seq.lower()
        kmers = [seq[start:start+kmer_length] for start in range(len(seq)-kmer_length+1)]
        # count the real kmers
        for kmer in kmers:
            fg_mers[kmer].add(real_index)

        # get 100 random regions in the same gene, get their kmer counts
        gene_start = int(overlap[num_posns_fields + 1])
        gene_end = int(overlap[num_posns_fields + 2])
        rand_starts = [random.randrange(gene_start, gene_end) for i in range(samples_per_overlap)]
        for random_index, rand_start in enumerate(rand_starts):
            if overlap.strand == '+':
                rand_seq = str(genome[overlap.chrom][rand_start : rand_start + windowsize + 1])
            else:
                rand_seq = str(- genome[overlap.chrom][max(0, rand_start - windowsize) : rand_start + 1])
            rand_seq = rand_seq.lower()
            rand_kmers = [rand_seq[start:start+kmer_length] for start in range(len(rand_seq)-kmer_length+1)]
            # count the bg kmers
            for kmer in rand_kmers:
                bg_mers[kmer].add('%s_%s' % (real_index,random_index))
    print 'done with overlap'

    with open(out_zscores, 'w') as outfile:
        # for each kmer we saw in fg, calculate zscore-- (x - mu) / sigma
        for kmer in fg_mers:
            z = zscore_normal(len(fg_mers[kmer]), num_fg_regions, len(bg_mers[kmer]), num_fg_regions * samples_per_overlap)
            fg_rate = len(fg_mers[kmer]) / float(num_fg_regions)
            bg_rate = len(bg_mers[kmer]) / float(num_fg_regions * samples_per_overlap)
            print >>outfile, '%s\t%s\t%s\t%.5f\t%s\t%s\t%.5f\t%.5f' % (kmer, len(fg_mers[kmer]), num_fg_regions, fg_rate, len(bg_mers[kmer]), num_fg_regions * samples_per_overlap, bg_rate, z)





@transform([kmer_zscore_compare_within_same_gene, kmer_zscore_compare_within_same_gene_pas_downstream],
    suffix('mer.vs_same_gene'), 'mer.vs_same_gene.png')
def plot_kmer_zscore_distribution(in_zscores, out_png):
    """Plot a histogram of the kmer Zscores"""
    zs = [float(l.strip().split('\t')[-1]) for l in open(in_zscores)]
    pyplot.figure()
    pyplot.hist(zs, bins=50)
    pyplot.savefig(out_png, dpi=240)


@transform([reproducible_positions_to_peaks, find_peaks_near_polyA_sites], suffix('.peaks'), '.hist_nearby.png')
#@transform(find_peaks_near_polyA_sites, suffix('.peaks'), '.hist_nearby.png')
def cluster_tall_reads_to_histogram(in_pileup, out_histogram):
    """Cluster tall reads together, then use the tallest ones as the center of
    a histogram of read signal
    """
    raise DeprecationWarning("cluster_tall_reads_to_histogram is deprecated. Use histogram_dist_to_downstream_read instead!")
    # ID all positions
    def make_index(features, _index=[0]):
        features.name = str(_index[0])
        _index[0] += 1
        return features
    id_pileup = pybedtools.BedTool(in_pileup).each(make_index).saveas()
    # sort reads by height
    #id_pileup = id_pileup.sort(chrThenScoreD=True)  # this sorts scores alphabetically!
    def sort_by_score(line):
        f = line.strip().split('\t')
        return float(f[4])
    sortedlines = sorted(open(id_pileup.fn).readlines(), key=sort_by_score, reverse=True)
    id_pileup = pybedtools.BedTool('\n'.join(sortedlines), from_string=True)

    # select >= 20 reads
    tall_pileup = id_pileup.filter(lambda f: int(f.score) >= 20)
    # slop them +/-30nt and intersect them with original, preserving all >=20.
    tall_with_nearby = tall_pileup.slop(g='%s.chrom.sizes' % genome_path(),
                                    b=30).intersect(id_pileup, wao=True, s=True)
    # from tallest to shortest, Keep track of what reads are already included and histogram those that haven't
    seen_reads = set()
    read_counts = defaultdict(int)
    for ovlp in tall_with_nearby:
        tall_site = ovlp[:6]
        other_site = ovlp[6:12]
        tall_id, tall_count = tall_site[3], int(tall_site[4])
        # count the tall read at position 0
        if tall_id not in seen_reads:
            read_counts[0] += tall_count
            seen_reads.add(tall_id)
        # count other reads at their distance from tall_site
        if other_site[0] != '.':  # some read overlaps with tall_site
            other_id, other_count = other_site[3], int(other_site[4])
            if other_id not in seen_reads:
                dist = (int(other_site[1]) + int(other_site[2])) / 2 - (int(tall_site[1]) + int(tall_site[2])) / 2
                if tall_site[5] == '-':
                    dist = -dist
                read_counts[dist] += other_count
                seen_reads.add(other_id)

    pyplot.figure()
    vals = [read_counts[i] for i in range(-30, 30 + 1)]
    pyplot.bar(range(-30,30+1), vals, align='center')
    pyplot.savefig(out_histogram, dpi=240)



    # TODO: cluster analysis: Take peaks with at least 20 reads that are higher all other peaks located within -30nt to +30nt. Make histograms of read counts at each position within this region.
    ##  that means we need the raw signal at these locations
    ## TODO: convert wig file to raw signal here-- bx-python?





@transform([reproducible_positions_to_peaks, find_peaks_near_polyA_sites], suffix('.peaks'), '.hist_distance.png')
#@transform(find_peaks_near_polyA_sites, suffix('.peaks'), '.hist_nearby.png')
def histogram_dist_to_downstream_read(in_pileup, out_histogram):
    """Histogram the distance between reads
    """
    # ID all positions
    def make_index(features, _index=[0]):
        features.name = str(_index[0])
        _index[0] += 1
        return features
    id_pileup = pybedtools.BedTool(in_pileup).each(make_index).saveas()

    # select >= 20 reads
    tall_pileup = id_pileup.filter(lambda f: int(f.score) >= 20)
    # slop them 100nt downstream and intersect them with original reads.
    # Keep only those that intersect (are 100nt downstream of a >=20 read)
    tall_with_nearby = tall_pileup.slop(g='%s.chrom.sizes' % genome_path(),
                                l=0, r=100, s=True).intersect(id_pileup, wo=True, s=True)
    # Keep track of comparisons are already included and
    seen_reads = set()
    seen_compares = set()
    read_counts = defaultdict(int)
    for ovlp in tall_with_nearby:
        tall_site = ovlp[:6]
        other_site = ovlp[6:12]
        tall_id, tall_count = tall_site[3], int(tall_site[4])
        other_id, other_count = other_site[3], int(other_site[4])
        # ignore self-overlaps
        if tall_id == other_id:
            continue
        # count the tall read at position 0 only once
        if tall_id not in seen_reads:
            read_counts[0] += tall_count
            seen_reads.add(tall_id)
        # count other reads at their distance from tall_site (before tall_site was slopped)
        compare = frozenset([tall_id, other_id])
        if compare not in seen_compares:
            if tall_site[5] == '+':
                dist = int(other_site[1]) - int(tall_site[1])
            elif tall_site[5] == '-':
                dist = int(tall_site[2]) - 1 - int(other_site[1])
            else:
                raise ValueError("Unrecognized strand in intersection: %s" % tall_site[5])
            read_counts[dist] += other_count
            seen_compares.add(compare)

    pyplot.figure()
    vals = [read_counts[i] for i in range(0, 100 + 1)]
    pyplot.bar(range(0,100+1), vals, align='center')
    #pyplot.gca().set_yscale('log')
    pyplot.savefig(out_histogram, dpi=240)







####@transform(pileup_starts, suffix('.pileup_reads'), '.meme.discovered.motifs')
#@transform([reproducible_positions_to_peaks, find_peaks_near_polyA_sites], suffix('.peaks'), '.meme.discovered.motifs')
@transform([reproducible_positions_to_peaks, find_peaks_near_polyA_sites], regex(r'(.*)\.peaks'), r'\1.meme.discovered.motifs')
def find_meme_motifs_around_sites(in_pileup, out_motifs):
    """Find motifs within the 21bp window around each binding site"""
    import hts_waterworks.motif_discovery as motif_discovery
    in_fasta = tempfile.NamedTemporaryFile(delete=False)
    genome = get_genome(None, None, False)
    for i, line in enumerate(sorted(open(in_pileup), key=lambda l: int(l.strip().split('\t')[4]), reverse=True)):
        #chrom, start, stop, count = line.strip().split('\t')
        chrom, start, stop, name, count, strand = line.strip().split('\t')[:6]
        seq = genome[chrom][int(start)-10 : int(start)+10+1]
        if strand == '+':
            seq = str(seq).upper()
        else:
            seq = str(-seq).upper()
        in_fasta.write('>seq_%s_%s_reads\n%s\n' % (i, count, seq))
        if i >= 5000:
            break
    in_fasta.close()
    motif_discovery.discover_meme_motifs(in_fasta.name, out_motifs)


@transform([reproducible_positions_to_peaks, find_peaks_near_polyA_sites], regex(r'(.*)\.peaks'), r'\1.peaks.fasta')
def get_sequence_around_sites(in_pileup, out_fasta):
    genome = get_genome(None, None, False)
    with open(out_fasta, 'w') as out_fasta:
        for i, line in enumerate(sorted(open(in_pileup), key=lambda l: int(l.strip().split('\t')[4]), reverse=True)):
            #chrom, start, stop, count = line.strip().split('\t')
            chrom, start, stop, name, count, strand = line.strip().split('\t')[:6]
            seq = genome[chrom][int(start)-10:int(start)+10+1]
            if strand == '+':
                seq = str(seq).upper()
            else:
                seq = str(-seq).upper()
            seq = seq.replace('T', 'U')
            out_fasta.write('>seq_%s_%s_reads\n%s\n' % (i, count, seq))
            #if i >= 100:
            #    break


@transform(get_sequence_around_sites, suffix('.fasta'), '.fasta.png')
def make_weblogo_around_sites(in_fasta, out_png):
    """Generate a weblogo of the alignment of sequences around a clip binding site"""
    windowsize = 10
    cmd = "weblogo -f {in_fasta} -F png_print -o {out_png} --fineprint '' --errorbars NO -Y NO --composition 'H. sapiens' -i -{windowsize} -s large --scale-width NO -U probability"
    cmd = cmd.format(**locals())
    sys_call(cmd)


@active_if(len(cfg.get('motifs', 'consensus_sequence')) > 0)
@split(None, ['%s.consensus.sites' % c for c in cfg.get('motifs', 'consensus_sequence').split(',')],
       cfg.get('motifs', 'consensus_sequence').split(','))
def search_genome_consensus(_, out_sites_list, consensus_list):
    """Search the genome (both strands) for matches to a consensus sequence"""
    genome = get_genome(None, None, False)
    out_sites = map(lambda f: open(f, 'w'), out_sites_list)
    for chrom in genome:
        seq = str(genome[chrom])
        revseq = reverseComplement(seq)
        for cons_seq, outfile in zip(consensus_list, out_sites):
            cons_re = re.compile(consensus_to_regex(cons_seq))
            try:
                i = -1
                while i < len(seq):
                    i = cons_re.search(seq, i+1).start()
                    outfile.write('\t'.join([chrom, str(i), str(i+len(cons_seq)), cons_seq, '0', '+']) + '\n')
            except AttributeError:
                pass
            try:
                i = -1
                while i < len(revseq):
                    i = cons_re.search(revseq, i+1).start()
                    outfile.write('\t'.join([chrom, str(len(seq) - (i+len(cons_seq))), str(len(seq) - i), cons_seq, '0', '-']) + '\n')
            except AttributeError:
                pass

@active_if(len(cfg.get('motifs', 'consensus_sequence')) > 0)
@follows(search_genome_consensus)
@split(pileup_starts, regex(r'(.*)\.pileup_reads'),
           add_inputs(search_genome_consensus),
           r'\1.pileup_reads.dist_to.consensus.*.dist', r'\1.pileup_reads.dist_to.consensus.%s.dist')
def closest_consensus_sites(in_files, _, out_dist_pattern):
    import hts_waterworks.annotation as annotation
    in_pileup, in_consensus_sites = in_files[0], in_files[1:]
    annotation.get_nearest_features((in_pileup, '%s.chrom.sizes' % genome_path()) + in_consensus_sites, _, out_dist_pattern)



@transform(reproducible_positions_to_peaks, suffix('.peaks'), '.intronic.peaks', 'hg19.refseq_genes.intron')
def get_intronic_peaks(in_pileup, out_intron_pileup, in_introns):
    """filter out non-intronic peaks from reproducible positions"""
    clip_sites = pybedtools.BedTool(in_pileup)
    introns = pybedtools.BedTool(in_introns)
    # keep the clip sites that overlap any introns
    intronic_clip = clip_sites.intersect(introns, s=True, u=True)
    intronic_clip.saveas(out_intron_pileup)


@transform(reproducible_positions_to_peaks, suffix('.peaks'), '.peaks_in_3utr.peaks', 'hg19.refseq_genes.utr3')
def get_3utr_peaks(in_pileup, out_3utr_pileup, in_3utr):
    """filter out non-intronic peaks from reproducible positions"""
    clip_sites = pybedtools.BedTool(in_pileup)
    utr3 = pybedtools.BedTool(in_3utr)
    # keep the clip sites that overlap any 3utr regions
    utr3_clip = clip_sites.intersect(utr3, s=True, u=True)
    utr3_clip.saveas(out_3utr_pileup)


@files(None, 'hg19.phastcons.conserved.elements')
def import_conservation_tracks(_in, out_elements):
    """grab the phastcons conserved elements from UCSC"""
    data = urllib.urlopen('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/phastConsElements46wayPlacental.txt.gz')
    data = gzip.GzipFile(fileobj=StringIO.StringIO(data))
    pybedtools.BedTool(data, from_string=True).cut([1,2,3,4,5]).saveas(out_elements)

@files(None, 'all_chroms.phyloP46way.vertebrates.bigwig')
def import_phylop_wig(_in, out_file):
    """download the entire phylop directory for the current genome"""
    print 'cd phylop; wget -m ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate/'
    print 'mv hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate/* .'
    print 'zcat * | gzip -c > all_chroms.wigFix.gz'
    print 'wigToBigWig all_chroms.wigFix.gz ../hg19.chrom.sizes ../%s' % out_file


@transform(reproducible_positions_to_peaks, suffix('.peaks'), '.peaks_in_introns')
def intronic_peaks_with_intron(in_pileup, out_sites_with_intron):
    """filter out non-intronic peaks, but also retain the intron we overlap with.
    When multiple introns are overlapped with, keep the shortest.
    """
    posns = pybedtools.BedTool(in_pileup)
    introns = pybedtools.BedTool('%s.refseq_genes.intron' % genome_path())
    posns_with_introns = intersect_same_strand_keep_shortest(posns, introns)
    posns_with_introns.saveas(out_sites_with_intron)



@transform(reproducible_positions_to_peaks, suffix('.peaks'), '.peaks_in_utr3')
def utr3_peaks_with_utr3(in_pileup, out_sites_with_utr3):
    """filter out non-3' utr peaks, but also retain the 3'utr we overlap with.
    When multiple 3' utrs are overlapped with, keep the shortest.
    """
    posns = pybedtools.BedTool(in_pileup)
    utr3 = pybedtools.BedTool('%s.refseq_genes.utr3' % genome_path())
    posns_with_utr3 = intersect_same_strand_keep_shortest(posns, utr3)
    posns_with_utr3.saveas(out_sites_with_utr3)


def intersect_same_strand_keep_shortest(posns, others):
    """Intersect the bedtools in_sites and in_other, requiring same strand and keeping only
    intersecting reads. When a read intersects multiple times, keep only the
    intersection with the shortest in_other
    """
    # report posns and genes they overlap with; only report same-strand overlap
    posns_in_other = posns.intersect(others, wo=True, s=True)

    # for each duplicate pileup, keep only intersection with the shortest transcript
    pileups_shortest_others = {}  # pileup -> line
    for p in posns_in_other:
        pileup = tuple(p[:6])
        other = p[6:]
        other_length = int(other[2]) - int(other[1])
        try:
            old_other = pileups_shortest_others[pileup][6:]
            old_length = int(old_other[2]) - int(old_other[1])
            if other_length < old_length:
                pileups_shortest_others[pileup] = p
        except KeyError:  # first time pileup seen
            pileups_shortest_others[pileup] = p

    # only shortest transcripts are saved; output results
    outstr = '\n'.join(str(p) for p in pileups_shortest_others.itervalues())
    final_sites = pybedtools.BedTool(outstr, from_string=True)
    return final_sites



@split([intronic_peaks_with_intron, utr3_peaks_with_utr3], regex(r'(.*\.peaks_in_\w+$)'),
    r'\1.pileup_reads.random_dist_to.consensus.*.dist',
    r'\1.pileup_reads.random_dist_to.consensus.%s.dist')
def random_dist_to_consensus_within_region(in_pileup, _out, out_template):
    pass
    

@transform('all_chroms.phyloP46way.vertebrates.bigwig', suffix('.bigwig'), '.array.pkl')
def bigwig_to_array(in_bigwig, out_array_dict):
    chrom_sizes = pybedtools.chromsizes('hg19')
    scores_bigwig = bx.bbi.bigwig_file.BigWigFile(open(in_bigwig))
    all_scores = dict((chrom, scores_bigwig.get_as_array(chrom, start, stop))
                        for chrom, (start, stop) in chrom_sizes.items())
    pickle.dump(all_scores, open(out_array_dict, 'w'), -1)


@follows(bigwig_to_array)
@transform([intronic_peaks_with_intron, utr3_peaks_with_utr3], regex(r'(.*\.peaks_in_\w+$)'), r'\1.conservation.png', 100, 'all_chroms.phyloP46way.vertebrates.array.pkl')
def phylop_scores_around_sites(in_sites, out_png, windowsize, phylop_array):
    """Draw the mean phylop score around in_sites, -/+ windowsize"""
    pileup = pybedtools.BedTool(in_sites)
    phylo_scores = pickle.load(open(phylop_array))

    scores_around = []
    rand_scores_around = []
    chrom_sizes = pybedtools.chromsizes('hg19')

    for index, p in enumerate(pileup):
        chrom = p.chrom
        score = phylo_scores[chrom][p.start - windowsize : p.start + windowsize + 1]
        if len(score) != windowsize * 2 + 1:
            print 'found a score near end of chrom?'
            print midpoint, len(score)
            continue
        scores_around.append(score)

        # for each real score, do 10 random windows around the currently
        # intersected feature-- (either intron on utr3)
        other_start, other_end = int(p[7]), int(p[8])
        other_length = other_end - other_start
        rand_posns = sp.floor(sp.round_(sp.rand(10) * (other_length - 1))) + other_start

        for rindex, r_pos in enumerate(rand_posns):
            r_pos = int(r_pos)
            score = phylo_scores[chrom][r_pos - windowsize : r_pos + windowsize + 1]
            if len(score) != windowsize * 2 + 1:
                print 'found a random score near end of chrom?'
                print midpoint, len(score)
                continue
            rand_scores_around.append(score)

    scores_around = sp.ma.masked_invalid(sp.array(scores_around))
    mean_per_base = scores_around.mean(axis=0)
    std_per_base = scores_around.std(axis=0)
    # standard error of the mean is std / sqrt(n)
    sem_per_base = std_per_base / sp.sqrt(index - 1)

    rand_scores_around = sp.ma.masked_invalid(sp.array(rand_scores_around))
    rand_mean_per_base = rand_scores_around.mean(axis=0)
    rand_std_per_base = rand_scores_around.std(axis=0)
    # standard error of the mean is std / sqrt(n)
    rand_sem_per_base = rand_std_per_base / sp.sqrt(index - 1 * 10)

    xs = sp.arange(-windowsize, windowsize+1)
    fig, ax = pyplot.subplots(1)
    ax.plot(xs, mean_per_base, 'k')
    ax.fill_between(xs, mean_per_base + sem_per_base, mean_per_base - sem_per_base,
                    facecolor='gray', alpha=.7, lw=0)

    ax.plot(xs, rand_mean_per_base, 'gray')
    ax.fill_between(xs, rand_mean_per_base + rand_sem_per_base, rand_mean_per_base - rand_sem_per_base,
                    facecolor='gray', alpha=.7, lw=0)

    ax.set_xlabel('distance from crosslink nt')
    ax.set_ylabel('conservation')
    pyplot.xlim((-windowsize - 10, windowsize + 10))
    pyplot.savefig(out_png, dpi=240)



@transform([reads_within_same_gene, intronic_peaks_with_intron, utr3_peaks_with_utr3],
    regex(r'(64_CLIP_all.*\.peaks_in_(\w+)$)'), r'\1.distribution_within_\2.png')
def read_distribution_within_gene(in_pileup_with_gene, out_png):
    """Given reads and their position within the gene, create a gene-normalized
    (0-100% of length) distribution...
    """
    posns_with_genes = pybedtools.BedTool(in_pileup_with_gene)
    if in_pileup_with_gene.endswith('.peaks_in_genes'):
        total_genes = len(pybedtools.BedTool('%s.refseq_genes.all' % genome_path()))
    elif in_pileup_with_gene.endswith('.peaks_in_introns'):
        total_genes = len(pybedtools.BedTool('%s.refseq_genes.intron' % genome_path()))
    elif in_pileup_with_gene.endswith('.peaks_in_utr3'):
        total_genes = len(pybedtools.BedTool('%s.refseq_genes.utr3' % genome_path()))
    else:
        import ipdb; ipdb.set_trace()
        raise RuntimeError('Name Error %s' % in_pileup_with_gene)
    all_starts = []
    for p in posns_with_genes:
        gene_start, gene_stop = int(p[7]), int(p[8])
        gene_length = gene_stop - gene_start
        relative_start = (p.start - gene_start) / float(gene_length) * 100
        all_starts.append(relative_start)
    total_per_base, _ = sp.histogram(all_starts, range=[min(0, *all_starts), max(100 *all_starts)], bins=1000)
    mean_per_base = total_per_base / float(total_genes)
    # get the variance for each site:  avg( (X - mu)^2 )
    def calc_variance(index):
        var_present = (1 - mean_per_base[index])**2 * total_per_base[index]
        var_missing = (0 - mean_per_base[index])**2 * (total_genes - total_per_base[index])
        # average variance
        return (var_present + var_missing) / total_genes
    var_per_base = sp.array([calc_variance(i) for i in range(len(total_per_base))])
    std_per_base = sp.sqrt(var_per_base)
    sem_per_base = std_per_base / sp.sqrt(total_genes)

    pyplot.figure()
    xs = sp.linspace(0,100,len(total_per_base), endpoint=False)
    fig, ax = pyplot.subplots(1)
    ax.plot(xs, mean_per_base, 'k')
    ax.fill_between(xs, mean_per_base + sem_per_base, mean_per_base - sem_per_base,
                    facecolor='gray', alpha=.7, lw=0)
    #ax.fill_between(xs, mean_per_base + std_per_base, mean_per_base - std_per_base,
    #                facecolor='gray', alpha=.7, lw=0)

    pyplot.xlim((-10, 110))
    pyplot.savefig(out_png, dpi=240)


#@transform(reads_within_same_gene, regex(r'(64_CLIP_all.*\.repro\.peaks_in_(\w+)$)'), r'\1.rna_fold_probs.vs_\2.png', r'seq_%s_%s_\1.\2')
#@transform([reads_within_same_gene, utr3_peaks_with_utr3, intronic_peaks_with_intron],
@transform([utr3_peaks_with_utr3],
    regex(r'(64_CLIP_all.*\.repro\.peaks_in_(\w+)$)'), r'\1.rna_fold_probs_w30.vs_\2.png')
def get_rna_pairing_probabilities(in_pileup_with_gene, out_png):
    """Use RNAplfold to generate the per-base probability of RNA folding in a given region
    """
    posns_with_genes = pybedtools.BedTool(in_pileup_with_gene)
    if in_pileup_with_gene.endswith('.peaks_in_genes'):
        total_genes = len(pybedtools.BedTool('%s.refseq_genes.all' % genome_path()))
    elif in_pileup_with_gene.endswith('.peaks_in_introns'):
        total_genes = len(pybedtools.BedTool('%s.refseq_genes.intron' % genome_path()))
    elif in_pileup_with_gene.endswith('.peaks_in_utr3'):
        total_genes = len(pybedtools.BedTool('%s.refseq_genes.utr3' % genome_path()))
    else:
        import ipdb; ipdb.set_trace()
        raise RuntimeError('Name Error %s' % in_pileup_with_gene)

    print len(posns_with_genes)
    windowsize = 100
    genome = get_genome(None, None, False)
    window_seqs = []
    rand_seqs = []
    for index, p in enumerate(posns_with_genes):
        gene_start, gene_stop = int(p[7]), int(p[8])
        gene_length = gene_stop - gene_start
        window_seq = genome[p.chrom][p.start - windowsize : p.start + windowsize + 1]
        if p.strand == '-':
            window_seq = -window_seq
        window_seq = str(window_seq).upper().replace('T', 'U')
        window_seqs.append(window_seq)

        # randomly sample this gene region
        rand_posns = sp.floor(sp.round_(sp.rand(10) * (gene_length - 1))) + gene_start
        for rindex, r_pos in enumerate(rand_posns):
            r_pos = int(r_pos)
            window_seq = genome[p.chrom][r_pos - windowsize : r_pos + windowsize + 1]
            if p.strand == '-':
                window_seq = -window_seq
            window_seq = str(window_seq).upper().replace('T', 'U')
            rand_seqs.append(window_seq)

    # seqs are generated; now run RNAplfold on them
    print '# folding foreground'
    pair_probs = get_RNAplfold_probabilities(window_seqs)
    print '# folding background'
    rand_pair_probs = get_RNAplfold_probabilities(rand_seqs)

    sp.save(out_png.replace('.png','.true_probs.npy'), pair_probs)
    sp.save(out_png.replace('.png','.random_probs.npy'), rand_pair_probs)

    # now sampling dance is done so let's plot it
    pair_probs = 1 - sp.array(pair_probs)
    mean_per_base = pair_probs.mean(axis=0)
    std_per_base = pair_probs.std(axis=0)
    # standard error of the mean is std / sqrt(n)
    sem_per_base = std_per_base / sp.sqrt(index - 1)

    rand_pair_probs = 1 - sp.array(rand_pair_probs)
    rand_mean_per_base = rand_pair_probs.mean(axis=0)
    rand_std_per_base = rand_pair_probs.std(axis=0)
    # standard error of the mean is std / sqrt(n)
    rand_sem_per_base = rand_std_per_base / sp.sqrt(index - 1 * 10)

    xs = sp.arange(-windowsize, windowsize+1)
    fig, ax = pyplot.subplots(1)
    ax.plot(xs, mean_per_base, 'k')
    ax.fill_between(xs, mean_per_base + sem_per_base, mean_per_base - sem_per_base,
                    facecolor='gray', alpha=.7, lw=0)
    ax.plot(xs, rand_mean_per_base, 'gray')
    ax.fill_between(xs, rand_mean_per_base + rand_sem_per_base, rand_mean_per_base - rand_sem_per_base,
                    facecolor='gray', alpha=.7, lw=0)

    ax.set_xlabel('distance from crosslink nt')
    ax.set_ylabel('average base pairing probability')
    pyplot.xlim((-windowsize - 10, windowsize + 10))
    pyplot.savefig(out_png, dpi=240)




def get_RNAplfold_probabilities(dna_seqs):
    """call RNAplfold for the given DNA sequence and return pairing probabilities
    for each dna position
    """
    all_probs = []
    with tempfile.NamedTemporaryFile() as tmp:
        proc = subprocess.Popen(['RNAplfold','-u','1','-W','30'], stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        tmp_name = os.path.split(tmp.name)[1]
        out_str = ''
        for i, seq in enumerate(dna_seqs):
            seq_name = tmp_name + '.file_%s' % i
            out_str += '>%s' % seq_name + '\n%s\n' % seq
        stdout, stderr = proc.communicate(out_str)
        for i, seq in enumerate(dna_seqs):
            seq_name = tmp_name + '.file_%s' % i
            probs_filename = seq_name + '_lunp'
            with open(probs_filename) as infile:
                probs = [float(l.strip().split('\t')[1]) for l in infile if '#' not in l]
            os.unlink(probs_filename)
            os.unlink(seq_name + '_dp.ps')
            all_probs.append(probs)
    return all_probs



@active_if(len(cfg.get('motifs', 'consensus_sequence')) > 0)
@follows(search_genome_consensus)
#@split(pileup_starts, regex(r'(.*)\.pileup_reads'),
@split([reproducible_positions_to_peaks, get_intronic_peaks,
        get_sites_nongenic_in_noncoding, get_3utr_peaks], regex(r'(.*)\.peaks'),
           add_inputs(search_genome_consensus),
           [r'\1.pileup_reads.dist_to.consensus.*.dist'],
           r'\1.pileup_reads.dist_to.consensus.%s.dist')
def closest_upstream_consensus_sites(in_files, _, out_dist_pattern):
    in_pileup, in_consensus_sites = in_files[0], in_files[1:]
    pileup = pybedtools.BedTool(in_pileup)
    def fiveprime_to_start_dist(feature):
        assert feature.strand in ['+', '-']
        if feature.strand == '+':
            feature.append(str(int(feature[7]) - int(feature[1])))
        else:
            feature.append(str(int(feature[1]) - int(feature[8]) - 1))
        return feature

    for c in in_consensus_sites:
        # get the pileups that have consensus sites in the 100bp upstream region, and are on the same strand
        consensus_sites = pybedtools.BedTool(c)
        pileup_to_cons = pileup.closest(consensus_sites, D='a', t='first', id=True, s=True)
        pileup_to_cons_real_dist = pileup_to_cons.each(fiveprime_to_start_dist)
        # filter on -100 to 0 window
        pileup_with_near_cons = pileup_to_cons_real_dist.filter(lambda a: -100 <= int(a[-1]) <= 0)
        dists = [int(p[-1]) for p in pileup_with_near_cons]
        with open(out_dist_pattern % c, 'w') as outfile:
            outfile.writelines('\n'.join(map(str, dists)))



@active_if(len(cfg.get('motifs', 'consensus_sequence')) > 0)
@follows(search_genome_consensus)
#@split(pileup_starts, regex(r'(.*)\.pileup_reads'),
#@split([reproducible_positions_to_peaks, get_intronic_peaks,
@split([get_sites_nongenic_in_noncoding, get_3utr_peaks],
    regex(r'(.*)\.peaks'),
           add_inputs(search_genome_consensus),
           [r'\1.pileup_reads.random_dist_to.consensus.*.dist'],
           r'\1.pileup_reads.random_dist_to.consensus.%s.dist')
def random_upstream_consensus_sites_same_genome_distribution2(in_files, _, out_dist_pattern):
    """ we want to shuffle the pileups, but need to keep the distribution in the
    genome the same (same count in introns, exons, etc). So we do weighted
    sampling of the regions, shuffling the original reads within the appropriate
    regions an equal number of times as there are percents...
    """

    in_pileup, in_consensus_sites = in_files[0], in_files[1:]
    pileup = pybedtools.BedTool(in_pileup)
    def fiveprime_to_start_dist(feature):
        assert feature.strand in ['+', '-']
        if feature.strand == '+':
            feature.append(str(int(feature[7]) - int(feature[1])))
        else:
            feature.append(str(int(feature[1]) - int(feature[8]) - 1))
        return feature
    
    if '.peaks_in_' in in_pileup:
        restrict_region = re.search('\.peaks_in_(\w+)', in_pileup).groups()[0]
        regions = {'introns' : 'hg19.refseq_genes.intron',
                   '3utr' : 'hg19.refseq_genes.utr3',
                   'noncoding' : 'hg19.refseq_genes.noncoding_with_wgRna_tRNA_lincRNA'
                   }
        shuffle_region = regions[restrict_region]
        all_shuffles = []
        print 'shuffling 100 times in %s' % shuffle_region
        for i in range(100):
            all_shuffles.append(pileup.shuffle(genome='hg19', incl=shuffle_region))
        all_shuffles = pybedtools.BedTool('\n'.join(map(str, all_shuffles)), from_string=True)
    else:
        # get all the genic regions...
        utr3 = pybedtools.BedTool('hg19.refseq_genes.utr3')
        utr5 = pybedtools.BedTool('hg19.refseq_genes.utr5')
        exon = pybedtools.BedTool('hg19.refseq_genes.exon')
        introns = pybedtools.BedTool('hg19.refseq_genes.intron')
        noncoding = pybedtools.BedTool('hg19.refseq_genes.noncoding_with_wgRna_tRNA_lincRNA')
    
        categories = [utr3, utr5, exon, introns, noncoding]
        genic_regions = '\n'.join(map(str, categories))
        genic_regions = pybedtools.BedTool(genic_regions, from_string=True).sort().merge()
        intergenic_regions = genic_regions.complement(genome='hg19')
        categories.append(intergenic_regions)
        
        # get the occurence of sites in each category
        leftovers = [pileup]
        count_per_category = [len(pileup)]
        for c in categories:
            print 'removing', c.fn
            # remove current category, then count how many were removed
            leftovers.append(leftovers[-1].intersect(c, v=True))
            count_per_category.append(count_per_category[0] - sum(count_per_category[1:]) - len(leftovers[-1]))
            print count_per_category
        
        percent_per_category = [count_per_category[i] / float(count_per_category[0]) * 100
                                for i in range(1,len(count_per_category))]
        # now shuffle the categories according to the genomic distribution
        # this will do 100 total shuffles
        all_shuffles = []
        for count, category in zip(percent_per_category, categories):
            count = int(round(count))
            print 'shuffling %s times for %s' % (count, category.fn)
            for i in range(count):
                all_shuffles.append(pileup.shuffle(genome='hg19', incl=category.fn))
        all_shuffles = pybedtools.BedTool('\n'.join(map(str, all_shuffles)), from_string=True)
    
    for c in in_consensus_sites:
        print 'closest for consensus sites %s' % c
        # get the pileups that have consensus sites in the 100bp upstream region, and are on the same strand
        consensus_sites = pybedtools.BedTool(c)
        pileup_to_cons = all_shuffles.closest(consensus_sites, D='a', t='first', id=True, s=True)
        pileup_to_cons_real_dist = pileup_to_cons.each(fiveprime_to_start_dist)
        # filter on -100 to 0 window
        pileup_with_near_cons = pileup_to_cons_real_dist.filter(lambda a: -100 <= int(a[-1]) <= 0)
        dists = [int(p[-1]) for p in pileup_with_near_cons]
        with open(out_dist_pattern % c, 'w') as outfile:
            outfile.writelines('\n'.join(map(str, dists)))





@active_if(len(cfg.get('motifs', 'consensus_sequence')) > 0)
@follows(search_genome_consensus)
#@split(pileup_starts, regex(r'(.*)\.pileup_reads'),
#@split([reproducible_positions_to_peaks, get_intronic_peaks,
@split([get_sites_nongenic_in_noncoding, get_3utr_peaks],
    regex(r'(.*)\.peaks'),
           add_inputs(search_genome_consensus),
           [r'\1.pileup_reads.random_dist_to.consensus.*.dist'],
           r'\1.pileup_reads.random_dist_to.consensus.%s.dist')
def random_upstream_consensus_sites_same_genome_distribution(in_files, _, out_dist_pattern):
    """ we want to shuffle the pileups, but need to keep the distribution in the
    genome the same (same count in introns, exons, etc). So we do weighted
    sampling of the regions, shuffling the original reads within the appropriate
    regions an equal number of times as there are percents...
    """

    in_pileup, in_consensus_sites = in_files[0], in_files[1:]
    pileup = pybedtools.BedTool(in_pileup)
    def fiveprime_to_start_dist(feature):
        assert feature.strand in ['+', '-']
        if feature.strand == '+':
            feature.append(str(int(feature[7]) - int(feature[1])))
        else:
            feature.append(str(int(feature[1]) - int(feature[8]) - 1))
        return feature
    
    if '.peaks_in_' in in_pileup:
        restrict_region = re.search('\.peaks_in_(\w+)', in_pileup).groups()[0]
        regions = {'introns' : 'hg19.refseq_genes.intron',
                   '3utr' : 'hg19.refseq_genes.utr3',
                   'noncoding' : 'hg19.refseq_genes.noncoding_with_wgRna_tRNA_lincRNA'
                   }
        shuffle_region = regions[restrict_region]
        all_shuffles = []
        print 'shuffling 100 times in %s' % shuffle_region
        for i in range(100):
            all_shuffles.append(pileup.shuffle(genome='hg19', incl=shuffle_region))
        all_shuffles = pybedtools.BedTool('\n'.join(map(str, all_shuffles)), from_string=True)
    else:
        # get all the genic regions...
        utr3 = pybedtools.BedTool('hg19.refseq_genes.utr3')
        utr5 = pybedtools.BedTool('hg19.refseq_genes.utr5')
        exon = pybedtools.BedTool('hg19.refseq_genes.exon')
        introns = pybedtools.BedTool('hg19.refseq_genes.intron')
        noncoding = pybedtools.BedTool('hg19.refseq_genes.noncoding_with_wgRna_tRNA_lincRNA')
    
        categories = [utr3, utr5, exon, introns, noncoding]
        genic_regions = '\n'.join(map(str, categories))
        genic_regions = pybedtools.BedTool(genic_regions, from_string=True).sort().merge()
        intergenic_regions = genic_regions.complement(genome='hg19')
        categories.append(intergenic_regions)
        
        # get the occurence of sites in each category
        leftovers = [pileup]
        count_per_category = [len(pileup)]
        for c in categories:
            print 'removing', c.fn
            # remove current category, then count how many were removed
            leftovers.append(leftovers[-1].intersect(c, v=True))
            count_per_category.append(count_per_category[0] - sum(count_per_category[1:]) - len(leftovers[-1]))
            print count_per_category
        
        percent_per_category = [count_per_category[i] / float(count_per_category[0]) * 100
                                for i in range(1,len(count_per_category))]
        # now shuffle the categories according to the genomic distribution
        # this will do 100 total shuffles
        all_shuffles = []
        for count, category in zip(percent_per_category, categories):
            count = int(round(count))
            print 'shuffling %s times for %s' % (count, category.fn)
            for i in range(count):
                all_shuffles.append(pileup.shuffle(genome='hg19', incl=category.fn))
        all_shuffles = pybedtools.BedTool('\n'.join(map(str, all_shuffles)), from_string=True)
    
    for c in in_consensus_sites:
        print 'closest for consensus sites %s' % c
        # get the pileups that have consensus sites in the 100bp upstream region, and are on the same strand
        consensus_sites = pybedtools.BedTool(c)
        pileup_to_cons = all_shuffles.closest(consensus_sites, D='a', t='first', id=True, s=True)
        pileup_to_cons_real_dist = pileup_to_cons.each(fiveprime_to_start_dist)
        # filter on -100 to 0 window
        pileup_with_near_cons = pileup_to_cons_real_dist.filter(lambda a: -100 <= int(a[-1]) <= 0)
        dists = [int(p[-1]) for p in pileup_with_near_cons]
        with open(out_dist_pattern % c, 'w') as outfile:
            outfile.writelines('\n'.join(map(str, dists)))


@active_if(len(cfg.get('motifs', 'consensus_sequence')) > 0)
@follows(search_genome_consensus)
@split(reproducible_positions_to_peaks, regex(r'(.*)\.peaks'),
           add_inputs(search_genome_consensus),
           [r'\1.pileup_reads.dist_to.cons_and_polya.*.dist.png'],
           r'\1.pileup_reads.dist_to.cons_and_polya.%s.dist.png',
           'hela.hg19.no_internal_priming.pas_seq_sites',
           'hg19.refseq_genes.utr3')
def closest_upstream_consensus_sites_with_polyA_in_utr3(in_files, _out, out_dist_png, pas_seq_sites, utr3_region):
    """Plot the distance from the consensus site to the nearest polyA site
    as well as from the consensus site to clip binding sites
    """

    in_pileup, in_consensus_sites = in_files[0], in_files[1:]

    pileup = pybedtools.BedTool(in_pileup)
    polya_sites = pybedtools.BedTool(pas_seq_sites)
    utr3 = pybedtools.BedTool(utr3_region)

    pileup_in_utr3 = pileup.intersect(utr3, u=True, s=True)
    polya_in_utr3 = polya_sites.intersect(utr3, u=True, s=True)
    
    # distance from closest is minimum between all feature edges. fix to be end - start
    def nearest_edge_distance(feature):
        # distance from consensus 3' end to binding site start
        assert feature.strand in ['+', '-']
        if feature.strand == '+':
            feature.append(str(int(feature[7]) - int(feature[1])))
        else:
            feature.append(str(int(feature[2]) - 1 - int(feature[7])))
        return feature

    # consensus is slopped 100nt. Get distance from orignal consensus' furthest point
    for c in in_consensus_sites:

        # get the pileups that have consensus sites in the 100bp upstream region, and are on the same strand
        consensus_sites = pybedtools.BedTool(c)
        len_consensus = consensus_sites.__iter__().next().length
        consensus_sites_downstream = consensus_sites.slop(s=True, r=100 - len_consensus,
                            l=0, g='%s.chrom.sizes' % genome_path())
        # get sites downstream of consensus, keep both consensus and clip site
        valid_consensus_with_pileups = consensus_sites_downstream.intersect(pileup_in_utr3, s=True, wb=True, wa=True)
        pileup_dists = [int(f[-1]) for f in valid_consensus_with_pileups.each(nearest_edge_distance)]

        # take the slopped consensus sites that have pileups and intersect them
        # to get the nearest polyA site
        valid_consensus_with_polya = consensus_sites_downstream.intersect(polya_in_utr3, s=True, wb=True, wa=True)
        polya_dists = [int(f[-1]) for f in valid_consensus_with_polya.each(nearest_edge_distance)]

        # TODO: rerun this
        pileup_dists = sp.histogram(pileup_dists, range=[0,100], bins=100)[0] / float(len(valid_consensus_with_pileups))
        polya_dists = sp.histogram(polya_dists, range=[0,100], bins=100)[0] / float(len(valid_consensus_with_polya))
        xs = sp.linspace(0,100,100,endpoint=False)

        #fig = pyplot.figure()
        #clip_ax = fig.add_subplot(111)
        #clip_lines = clip_ax.plot(xs, pileup_dists, 'b', label='CLIP sites')
        #clip_ax.set_xlabel('distance from AWTAAA consensus')
        #clip_ax.set_ylabel("fraction of 3'utr CLIP sites")
        #
        #pas_ax = clip_ax.twinx()
        #pas_lines = pas_ax.plot(xs, polya_dists, 'g', label='PAS sites')
        #pas_ax.set_ylabel("fraction of 3'utr PAS sites")
        #pyplot.legend([clip_lines, pas_lines], [clip_lines.get_label(), pas_lines.get_label()])
        #pyplot.savefig(out_dist_png % c)

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.plot(xs, pileup_dists, 'b', label='CLIP sites', lw=2)
        ax.plot(xs, polya_dists, 'g', label='PAS sites', lw=2)
        ax.set_xlabel('distance from AWTAAA consensus')
        ax.set_ylabel("fraction of 3'utr sites with AWTAAA")
        pyplot.legend()
        pyplot.savefig(out_dist_png % c)



@transform([closest_upstream_consensus_sites,
           random_upstream_consensus_sites_same_genome_distribution],
    suffix('.dist'), '.dist.png')
def plot_individual_closest_upstream_consensus_sites(in_dists, out_png):
    """Plot the distance distribution for each clip dataset to the nearest
    upstream consensus site.
    """
    all_dists = map(int, open(in_dists))
    pyplot.figure()
    pyplot.hist(all_dists, bins=50)
    pyplot.savefig(out_png, dpi=240)

#@collate([closest_upstream_consensus_sites, random_upstream_consensus_sites_same_genome_distribution],
#    )

@merge(closest_upstream_consensus_sites, r'all.pileup_reads.dist_to.all.consensus.dist.png')
def plot_all_closest_upstream_consensus_sites(in_all_dists, out_png):
    """Plot the distance distribution of all upstream consensus sites"""
    all_dists = []
    for f in map(open, in_all_dists):
        all_dists.extend(map(int, f.readlines()))
    pyplot.figure()
    pyplot.hist(all_dists, bins=50)
    pyplot.savefig(out_png, dpi=240)


#@merge(closest_upstream_consensus_sites_with_polyA_in_utr3, r'all.pileup_reads.dist_to.all.cons_and_polya.dist.png')
#def plot_from_upstream_consensus_with_polya_sites(in_all_dists, out_png):
#    clip_site_dists = []
#    polya_site_dists = []
#    for f in map(open, in_all_dists):
#        clip_site_dists.extend(map(int, f.readline().strip().split('\t')))
#        polya_site_dists.extend(map(int, f.readline().strip().split('\t')))
#    pyplot.figure()
#    pyplot.hist(clip_site_dists, bins=100)
#    pyplot.hist(polya_site_dists, bins=100)
#    pyplot.savefig(out_png, dpi=240)
#


@transform(closest_consensus_sites, suffix('.dist'),
           '.dist.png')
def plot_closest_consensus_sites(in_dist, out_png):
    """plot distance to nearest consensus sites"""
    import hts_waterworks.annotation as annotation
    annotation.plot_nearest_features(in_dist, out_png, window_size=5)

#@transform(annotation.refseq_genes_to_regions, suffix('.intron'), '.intron.5p_base')
@transform('hg19.refseq_genes.intron', suffix('.intron'), '.intron.5p_base')
def get_5p_introns(in_introns, out_5p_introns):
    """create a feature for the 5' base of each intron region"""
    in_introns = pybedtools.BedTool(in_introns)
    with open(out_5p_introns, 'w') as outfile:
        for r in in_introns:
            if r.strand == '+':
                outfile.write('\t'.join([r.chrom, str(r.start), str(r.start+1)] + r.fields[3:]) + '\n')
            else:
                outfile.write('\t'.join([r.chrom, str(r.end-1), str(r.end)] + r.fields[3:]) + '\n')

#@transform(annotation.refseq_genes_to_regions, suffix('.intron'), '.intron.3p_base')
@transform('hg19.refseq_genes.intron', suffix('.intron'), '.intron.3p_base')
def get_3p_introns(in_introns, out_3p_introns):
    """create a feature for the 3' base of each intron region"""
    in_introns = pybedtools.BedTool(in_introns)
    with open(out_3p_introns, 'w') as outfile:
        for r in in_introns:
            if r.strand == '+':
                outfile.write('\t'.join([r.chrom, str(r.end-1), str(r.end)] + r.fields[3:]) + '\n')
            else:
                outfile.write('\t'.join([r.chrom, str(r.start), str(r.start+1)] + r.fields[3:]) + '\n')


#@split([pileup_starts, '*.custom.peaks', reproducible_positions_to_peaks], regex('(.*\.(peaks|pileup_reads))'),
@split(reproducible_positions_to_peaks, regex(r'(.*)\.peaks'),
#@split([reproducible_positions],
    add_inputs('hg19.refseq_genes.intron', get_5p_introns, get_3p_introns),
    [r'\1.dist_to.intron_5p.distances', r'\1.dist_to.intron_3p.distances'])
def dist_within_introns(in_files, out_files):
    """For reads within introns, check preference for 5' or 3' introns,
    and compare vs. reads randomly distributed within the *same* introns"""
    in_reads = in_files[0]
    intron5, intron3 = in_files[-2:]
    intron_regions = filter(lambda x: x.endswith('.intron'), in_files)[0]
    out_5, out_3 = out_files
    print in_reads, intron5, intron3, intron_regions, out_5, out_3

    in_reads = pybedtools.BedTool(in_reads)
    intron_regions = pybedtools.BedTool(intron_regions)
    intron5 = pybedtools.BedTool(intron5)
    intron3 = pybedtools.BedTool(intron3)
    dists_5, dists_3 = [], []

    # get reads within introns only
    print '# getting intronic read distances'
    reads_in_introns = in_reads.intersect(intron_regions, wo=True)
    dists_5 = [e.fields[-1] for e in reads_in_introns.closest(intron5, D='a')]
    dists_3 = [e.fields[-1] for e in reads_in_introns.closest(intron3, D='a')]
    # for each read, generate random sites uniformly within the overlapping intron
    print '# getting random sites'
    random_sites = []
    for r in reads_in_introns:
        o_intron = pybedtools.Interval(r.fields[6], int(r.fields[7]), int(r.fields[8]), *r.fields[9:-1])
        for i in xrange(10):
            rand_start = random.randrange(o_intron.start, o_intron.stop - r.length + 1)
            rand_end = rand_start + r.length
            random_sites.append('\t'.join([o_intron.chrom, str(rand_start), str(rand_end), '.', '0', r.strand]))
    random_sites = '\n'.join(random_sites)
    random_sites = pybedtools.BedTool(random_sites, from_string=True)
    rand_dists_5 = [e.fields[-1] for e in random_sites.closest(intron5, D='a')]
    rand_dists_3 = [e.fields[-1] for e in random_sites.closest(intron3, D='a')]

    # save for R plotting
    with open(out_5, 'w') as outfile:
        outfile.write("intron_5p\tRandom\n")
        outfile.writelines('\t'.join(map(str, l)) + '\n' \
                    for l in itertools.izip_longest(dists_5, rand_dists_5, fillvalue='NA'))
    with open(out_3, 'w') as outfile:
        outfile.write("intron_5p\tRandom\n")
        outfile.writelines('\t'.join(map(str, l)) + '\n' \
                    for l in itertools.izip_longest(dists_3, rand_dists_3, fillvalue='NA'))






#Fig 4. Intronic CstF64-binding sites.
#1. Total number of intronic sites (using different read count threshold 1-10)
#2. Are AWTAAA overrepresented within 100nt upstream of intronic CstF64 binding sites (read count threshold 1-10)? Z score of AWTAAA comparing 100nt upstream of intronic CstF64 binding sites and same-sized fragment randomly selected from the introns of the same gene.
#3. Do these sites cause Pol II pause: compare with HeLa Pol II ChIP-seq data (visualize data: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12781)

@transform('hela_gene_names-071410.txt', suffix('.txt'), '.pas_seq_sites')
def make_pas_seq_sites(in_datafile, out_pas_sites):
    with open(in_datafile) as infile:
        with open(out_pas_sites, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue
                gene = fields[3:6]
                strand = fields[6]
                pas_sites = fields[7:]
                for p in pas_sites:
                    if not p:
                        continue
                    print p, p.split('-')
                    chrom, start, count = p.split('-')
                    print >>outfile, '\t'.join([chrom, start, str(int(start)+1), '.', count, strand])


@transform('hela.hg19.no_internal_priming.pas_seq_sites.extended', suffix(''), '.bedgraph')
def bed_score_to_bedgraph(in_pileup, out_bedgraph):
    """Convert a pileup (bed file with score as read count) to a wiggle file"""
    pileups = pybedtools.BedTool(in_pileup)
    with open(out_bedgraph, 'w') as outfile:
        #outfile.write('track type=bedGraph name="%s"\n' % out_bedgraph)
        for p in pileups:
            outfile.write('\t'.join(p[:3] + [p[4] if p.strand == '+' else str(-int(p[4]))]) + '\n')


@transform('hela.hg19.no_internal_priming.pas_seq_sites.extended', suffix(''), '.bedgraph')
def bed_score_to_bed_reads(in_pileup, out_bedgraph):
    """Convert a pileup (bed file with score as read count) to one line per read"""
    pileups = pybedtools.BedTool(in_pileup)
    with open(out_bedgraph, 'w') as outfile:
        #outfile.write('track type=bedGraph name="%s"\n' % out_bedgraph)
        for p in pileups:
            for i in range(int(p[-2])):
                outfile.write('\t'.join(p[:3] + [p[4] if p.strand == '+' else str(-int(p[4]))]) + '\n')



@follows(plot_kmer_zscore_distribution, find_meme_motifs_around_sites,
         plot_all_closest_upstream_consensus_sites,
         plot_individual_closest_upstream_consensus_sites,
         dist_within_introns,
         make_weblogo_around_sites, histogram_dist_to_downstream_read,
         reproducible_positions, reproducible_motifs, calculate_read_gene_fdr)
def clip_done():
    pass

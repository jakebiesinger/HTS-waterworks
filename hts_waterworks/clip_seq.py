
import itertools
import tempfile
import re

from ipdb import set_trace as breakpoint

from ruffus import (transform, follows, merge, split, collate,
                    add_inputs, regex, suffix, jobs_limit)
from ruffus.task import active_if

from hts_waterworks.bootstrap import cfg, genome_path, get_genome
from hts_waterworks.utils.ruffus_utils import sys_call
from hts_waterworks.utils.common import reverseComplement, consensus_to_regex
import hts_waterworks.mapping as mapping

import scipy as sp

@active_if(cfg.getboolean('CLIP-seq', 'truncate_to_starts'))
@transform(mapping.all_mappers_output, suffix('.mapped_reads'),
           '.read_starts.mapped_reads', cfg.getint('CLIP-seq', 'truncate_to_starts_offset'))
def truncate_to_starts(in_bed, out_bed, offset):
    """Truncate reads as only the start + offset base.  For use in CLIP-seq
    experiments
    """
    with open(out_bed, 'w') as outfile:
        for line in open(in_bed):
            chrom,start,stop,name,score,strand = line.strip().split('\t')
            start = int(start)
            outfile.write('\t'.join(map(str, [chrom, max(0, start + offset),
                            max(0, start + offset + 1), name, score, strand])) + '\n')

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

#@merge(pileup_starts, regex(r'(.*)\.(pileup_reads)'), ['reproducible.positions', '*.pileup.repro.peaks'], r'\1.%s.pileup.repro.peaks')
@collate(pileup_starts, regex(r'(.*)\.pileup_reads'), ['reproducible.positions', '*.pileup.repro.peaks'], r'%s.repro.peaks')
def reproducible_positions(in_pileups, out_data, out_pattern):
    """Assess the number of reads at each posiiton in each experiment.
    What percent of positions exactly X reads in both (all three?) experiments
    """
    out_data = out_data[0]
    count_by_position = {}  # position -> file index -> count
    for findex, infile in enumerate(map(open, in_pileups)):
        for line in infile:
            chrom, start, stop, count = line.strip().split('\t')
            count = int(count)
            if (chrom, start) not in count_by_position:
                count_by_position[(chrom, start)] = sp.array([0] * len(in_pileups))
            count_by_position[(chrom, start)][findex] = count
    # make summary table
    with open(out_data, 'w') as outfile:
        outfile.write('{|\n|-\n| ')
        outfile.write(' || '.join(['min_count'] + ['at_least_%s_exp' % i for i in range(1,len(in_pileups)+1)] + ['total_exp1']) + '\n')
        for mincount in range(1,10):
            outfile.write('|-\n| ' + str(mincount) + ' || ')
            total_passing = 0
            total_passing_other_exp = [0] * (len(in_pileups)+1)
            for pos in count_by_position:
                if count_by_position[pos][0] == mincount:
                    total_passing += 1
                    other_exp_passing = sum(count_by_position[pos] >= mincount)
                    total_passing_other_exp[other_exp_passing] += 1
            for exp_count in range(1, len(in_pileups)+1):
                sum_exp_count = sum(total_passing_other_exp[exp_count:])
                outfile.write(str(sum_exp_count) + ' || ')
            outfile.write(str(total_passing) + '\n')

        for mincount in [10]:
            outfile.write('|-\n| ' + str(mincount) + '+ || ')
            total_passing = 0
            total_passing_other_exp = [0] * (len(in_pileups)+1)
            for pos in count_by_position:
                if count_by_position[pos][0] >= mincount:
                    total_passing += 1
                    other_exp_passing = sum(count_by_position[pos] >= mincount)
                    total_passing_other_exp[other_exp_passing] += 1
            for exp_count in range(1, len(in_pileups)+1):
                sum_exp_count = sum(total_passing_other_exp[exp_count:])
                outfile.write(str(sum_exp_count) + ' || ')
            outfile.write(str(total_passing) + '\n')
        outfile.write('|}\n')
    
    # write out reads that pass min read criteria in all datasets
    min_count = 10
    out_passing = {}
    for exp in in_pileups:
        out_passing[exp] = open(out_pattern % exp, 'w')
    for pos in count_by_position:
        if sum(count_by_position[pos] >= mincount) >= 3:
            for expindex, exp in enumerate(in_pileups):
                if count_by_position[pos][expindex] > 0:
                    out_passing[exp].write('\t'.join(map(str, pos + ('.', count_by_position[pos], '+'))) + '\n')


@merge(pileup_starts, 'reproducible.motifs', 5)
def reproducible_motifs(in_pileups, out_mer_data, mer_length):
    """Examine the occurrence of pentamers overlapping pileup sites"""
    mer_counts = {}  # kmer -> exp counts
    total_sites = [0] * len(in_pileups)
    genome = get_genome(None, None, False)
    for findex, infile in enumerate(map(open, in_pileups)):
        for line in infile:
            total_sites[findex] += 1
            chrom, start, stop, count = line.strip().split('\t')
            seq = str(genome[chrom][int(start) - int(round(mer_length/2.)) : int(start) + mer_length/2]).upper()
            for seq_start in range(len(seq) - mer_length + 1):
                seq_slice = seq[seq_start : seq_start + mer_length]
                if seq_slice not in mer_counts:
                    mer_counts[seq_slice] = [0] * len(in_pileups)
                mer_counts[seq_slice][findex] += 1
    with open(out_mer_data, 'w') as outfile:
        outfile.write('{|\n| ')
        outfile.write(' || '.join('kmer count in %s' % f for f in in_pileups) + '\n')
        for kmer in mer_counts:
            outfile.write('|-\n| ' + kmer + '|| ' + ' || '.join(map(str, mer_counts[kmer])) + '\n')
        outfile.write('|-\n| total sites || ' + ' || '.join(map(str, total_sites)))
        outfile.write('|}\n')

@transform(pileup_starts, suffix('.pileup_reads'), '.meme.discovered.motifs')
def find_meme_motifs_around_sites(in_pileup, out_motifs):
    import hts_waterworks.motif_discovery as motif_discovery
    in_fasta = tempfile.NamedTemporaryFile(delete=False)
    genome = get_genome(None, None, False)
    for i, line in enumerate(sorted(open(in_pileup), key=lambda l: int(l.strip().split('\t')[3]), reverse=True)):
        chrom, start, stop, count = line.strip().split('\t')
        in_fasta.write('>seq_%s_%s_reads\n%s\n' % (i, count, str(genome[chrom][int(start)-10:int(start)+10])))
        if i >= 1000:
            break
    in_fasta.close()
    motif_discovery.discover_meme_motifs(in_fasta.name, out_motifs)

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

@transform(closest_consensus_sites, suffix('.dist'),
           '.dist.png')
def plot_closest_consensus_sites(in_dist, out_png):
    """plot distance to nearest consensus sites"""
    import hts_waterworks.annotation as annotation
    annotation.plot_nearest_features(in_dist, out_png, window_size=5)



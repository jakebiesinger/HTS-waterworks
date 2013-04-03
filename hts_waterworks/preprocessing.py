"""preprocessing.py
    module for filtering out low quality reads, trimming poly-A tails, etc
"""

import os
import re
from os.path import join
from subprocess import Popen, PIPE, check_call, CalledProcessError
import gzip
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from Bio import SeqIO
from ruffus import (transform, follows, collate, files, split, merge,
                    suffix, mkdir, jobs_limit, output_from)
from ruffus.task import active_if

from hts_waterworks.utils.ruffus_utils import (sys_call, main_logger as log,
                                               main_mutex as log_mtx)
from hts_waterworks.bootstrap import cfg
from hts_waterworks.utils.common import parseFastq

# filtering
original_reads = '*.fastq'
prev_output = original_reads
prev_suffix = '.fastq'

@active_if(cfg.getboolean('filtering', 'convert_sanger_to_illumina'))
@transform(prev_output, suffix(prev_suffix), '.fastq_illumina')
def convert_fastq(in_fastq, out_fastq):
    'convert sanger fastq format (phred-33) to illumina format (phred-64)'
    base_out = os.path.splitext(out_fastq)[0]
    records = SeqIO.parse(in_fastq, "fastq")
    with open(base_out, 'w') as outfile:
        SeqIO.write(records, outfile, "fastq-illumina")
    check_call('gzip %s' % base_out, shell=True)
if cfg.getboolean('filtering', 'convert_sanger_to_illumina'):
    prev_output = convert_fastq
    prev_suffix = ''



@active_if(cfg.getboolean('filtering', 'clip_adapter'))
@transform(prev_output, suffix(prev_suffix), '.noAdapter')
def clip_adapter(in_fastq, out_fastq):
    'remove adapter sequence from raw reads'
    cmd1 = 'cat %s' % in_fastq
    cmd2 = 'fastx_clipper -o %s -a %s' % (out_fastq,
                                    cfg.get('filtering', 'adapter_sequence'))
    p1 = Popen([cmd1], stdout=PIPE, shell=True)
    p2 = Popen([cmd2], stdin=p1.stdout, shell=True)
    p2.communicate()
    if p1.returncode:
        raise CalledProcessError(p1.returncode, cmd1)
    if p2.returncode:
        raise CalledProcessError(p2.returncode, cmd2)
if cfg.getboolean('filtering', 'clip_adapter'):
    prev_output = clip_adapter
    prev_suffix = ''


@active_if(cfg.getboolean('filtering', 'trim_reads'))
@transform(prev_output, suffix(prev_suffix), '.trimmed')
def trim_reads(in_fastq, out_fastq):
    'trim leading and/or trailing bases from all reads'
    cmd1 = 'cat %s' % in_fastq
    cmd2 = 'fastx_trimmer -o %s -f %s -l %s' % (out_fastq,
                                    cfg.get('filtering', 'trim_start'),
                                    cfg.get('filtering', 'trim_end'))
    p1 = Popen([cmd1], stdout=PIPE, shell=True)
    p2 = Popen([cmd2], stdin=p1.stdout, shell=True)
    p2.communicate()
    if p1.returncode:
        raise CalledProcessError(p1.returncode, cmd1)
    if p2.returncode:
        raise CalledProcessError(p2.returncode, cmd2)
if cfg.getboolean('filtering', 'trim_reads'):
    prev_output = trim_reads
    prev_suffix = ''


@active_if(cfg.get('filtering', 'trim_regex') != '')
@transform(prev_output, suffix(prev_suffix), '.trim_regex',
           cfg.get('filtering', 'trim_regex'))
def trim_regex(in_fastq, out_fastq, trim_pattern):
    """Search the reads for a regex, and trim everything matching the pattern
        and all succeeding sequence.
    
    """
    pattern = re.compile(trim_pattern)
    with gzip.open(in_fastq) as infile:
        with gzip.open(out_fastq, 'w') as outfile:
            for header, seq, qual in parseFastq(infile):
                matches = [m.span() for m in pattern.finditer(seq)]
                if len(matches) > 0:
                    # match to re found--
                    #   trim the right-most hit and add the trimmed sequence to the read ID
                    m = matches[-1]
                    header = seq[m[0]:] + '_' + header
                    seq = seq[:m[0]]
                    qual = qual[:m[0]]
                if len(matches) > 0 or not cfg.getboolean('filtering', 'require_regex'):
                    if len(seq) >= 10:  # TODO: add adjustable min length
                        outfile.write('@%s\n%s\n+%s\n%s\n' % (header, seq,
                                                              header, qual))
if cfg.get('filtering', 'trim_regex') != '':
    prev_output = trim_regex
    prev_suffix = ''


@active_if(cfg.getboolean('filtering', 'filter_artifacts'))
@transform(prev_output, suffix(prev_suffix), '.noArtifacts')
def filter_artifacts(in_fastq, out_fastq):
    """Remove sequences with only 3/4 nucleotides (no As or no Ts or ...)"""
    cmd1 = 'zcat %s' % in_fastq
    cmd2 = 'fastx_artifacts_filter -o %s -z' % (out_fastq)
    p1 = Popen([cmd1], stdout=PIPE, shell=True)
    p2 = Popen([cmd2], stdin=p1.stdout, shell=True)
    p2.communicate()
    if p1.returncode:
        raise CalledProcessError(p1.returncode, cmd1)
    if p2.returncode:
        raise CalledProcessError(p2.returncode, cmd2)
if cfg.getboolean('filtering', 'filter_artifacts'):
    prev_output = filter_artifacts
    prev_suffix = ''


@active_if(cfg.getboolean('filtering', 'filter_quality'))
@transform(prev_output, suffix(prev_suffix), '.min_qual',
        cfg.getint('filtering', 'filter_min_quality'),
        cfg.getint('filtering', 'filter_percent_bases_at_min'))
def filter_min_quality(in_fastq, out_fastq, min_qual, min_percent):
    """Remove sequences that have < min_precent bases with quality < min_qual"""
    cmd1 = 'cat %s' % in_fastq
    cmd2 = 'fastq_quality_filter -o %s -q %s -p %s' % (out_fastq,
                                                        min_qual, min_percent)
    p1 = Popen([cmd1], stdout=PIPE, shell=True)
    p2 = Popen([cmd2], stdin=p1.stdout, shell=True)
    p2.communicate()
    if p1.returncode:
        raise CalledProcessError(p1.returncode, cmd1)
    if p2.returncode:
        raise CalledProcessError(p2.returncode, cmd2)
if cfg.getboolean('filtering', 'filter_quality'):
    prev_output = filter_min_quality
    prev_suffix = ''


@transform(prev_output, suffix(prev_suffix), '.read_lengths.png')
def read_length_histogram(in_fastq, out_hist):
    """draw histogram of read lengths"""
    length_count = defaultdict(int)
    for header, seq, qual in parseFastq(gzip.open(in_fastq)):
        length_count[len(seq)] += 1
    pyplot.figure()
    pyplot.bar(range(1,max(length_count)+1),
               [length_count[l] for l in range(1, max(length_count)+1)],
               align='center')
    pyplot.savefig(out_hist)


@transform([original_reads, prev_output],
        suffix(prev_suffix), '.qual_stats')
def quality_stats(in_fastq, out_stats):
    """aggregate quality score statistics for original and filtered reads"""
    cmd1 = 'cat %s' % in_fastq
    cmd2 = 'fastx_quality_stats -o %s' % (out_stats)
    p1 = Popen([cmd1], stdout=PIPE, shell=True)
    p2 = Popen([cmd2], stdin=p1.stdout, shell=True)
    p2.communicate()
    if p1.returncode:
        raise CalledProcessError(p1.returncode, cmd1)
    if p2.returncode:
        raise CalledProcessError(p2.returncode, cmd2)


@transform(quality_stats, suffix('.qual_stats'), '.qual_boxplot.png')
def quality_boxplot(in_stats, out_boxplot):
    """draw a boxplot for the quality scores"""
    cmd = 'fastq_quality_boxplot_graph.sh -t %s -i %s -o %s' % (in_stats,
                                                in_stats, out_boxplot)
    check_call(cmd)


@transform(quality_stats, suffix('.qual_stats'), '.qual_nuc_dist.png')
def quality_nuc_dist(in_stats, out_dist):
    'show the nucleotide distribution across the reads'
    cmd = 'fastx_nucleotide_distribution_graph.sh -t %s -i %s -o %s' % (
        in_stats, in_stats, out_dist)
    check_call(cmd)


@follows(mkdir('summaries'))
@merge([original_reads, clip_adapter, trim_reads, trim_regex,
        filter_artifacts, filter_min_quality],
    #join('..', 'summaries', 'fastq.wikisummary'))
    'fastq.wikisummary')
def summarize_fastq_reads(in_fastq, out_summary):
    """Summarize fastq line counts"""
    with open(out_summary, 'w') as outfile:
        outfile.write("""
{| class="wikitable"
|+ Summary of raw read counts
!scope="col" | Dataset
!scope="col" | Number of raw reads
|-
""")
        for infile in in_fastq:
            for count, line in enumerate(gzip.open(infile)):
                pass
            outfile.write("| %s || %s\n|-\n" % (infile, count//4))
        outfile.write('|}\n')


final_output = prev_output

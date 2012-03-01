"""preprocessing.py
    module for filtering out low quality reads, trimming poly-A tails, etc
"""


import re

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
final_output = '*.fastq'
prev_suffix = '.fastq'

@active_if(cfg.getboolean('filtering', 'convert_sanger_to_illumina'))
@transform(final_output, suffix(prev_suffix), '.fastq_illumina')
def convert_fastq(in_fastq, out_fastq):
    'convert sanger fastq format to illumina format'
    records = SeqIO.parse(in_fastq, "fastq")
    with open(out_fastq, 'w') as outfile:
        SeqIO.write(records, outfile, "fastq-illumina")
if cfg.getboolean('filtering', 'convert_sanger_to_illumina'):
    final_output = convert_fastq
    prev_suffix = ''

@active_if(cfg.getboolean('filtering', 'clip_adapter'))
@transform(final_output, suffix(prev_suffix), '.noAdapter')
def clip_adapter(in_fastq, out_fastq):
    'remove adapter sequence from raw reads'
    cmd = 'fastx_clipper -i %s -o %s -a %s' % (in_fastq, out_fastq,
                                    cfg.get('filtering', 'adapter_sequence'))
    sys_call(cmd)
if cfg.getboolean('filtering', 'clip_adapter'):
    final_output = clip_adapter
    prev_suffix = ''

@active_if(cfg.getboolean('filtering', 'trim_reads'))
@transform(final_output, suffix(prev_suffix), '.trimmed')
def trim_reads(in_fastq, out_fastq):
    'trim leading and/or trailing bases from all reads'
    cmd = 'fastx_trimmer -i %s -o %s -f %s -l %s' % (in_fastq, out_fastq,
                                    cfg.get('filtering', 'trim_start'),
                                    cfg.get('filtering', 'trim_end'))
    sys_call(cmd)
if cfg.getboolean('filtering', 'trim_reads'):
    final_output = trim_reads
    prev_suffix = ''


@active_if(cfg.get('filtering', 'trim_regex') != '')
@transform(final_output, suffix(prev_suffix), '.trim_regex',
           cfg.get('filtering', 'trim_regex'))
def trim_regex(in_fastq, out_fastq, trim_pattern):
    """Search the reads for a regex, and trim everything matching the pattern
        and all succeeding sequence.
    
    """
    pattern = re.compile(trim_pattern)
    with open(in_fastq) as infile:
        with open(out_fastq, 'w') as outfile:
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
    final_output = trim_regex
    prev_suffix = ''

@active_if(cfg.getboolean('filtering', 'filter_artifacts'))
@transform(final_output, suffix(prev_suffix), '.noArtifacts')
def filter_artifacts(in_fastq, out_fastq):
    """Remove sequences with only 3/4 nucleotides (no As or no Ts or ...)"""
    cmd = 'fastx_artifacts_filter -i %s -o %s' % (in_fastq, out_fastq)
    sys_call(cmd)
if cfg.getboolean('filtering', 'filter_artifacts'):
    final_output = filter_artifacts
    prev_suffix = ''

@active_if(cfg.getboolean('filtering', 'filter_quality'))
@transform(final_output, suffix(prev_suffix), '.min_qual',
        cfg.getint('filtering', 'filter_min_quality'),
        cfg.getint('filtering', 'filter_percent_bases_at_min'))
def filter_min_quality(in_fastq, out_fastq, min_qual, min_percent):
    """Remove sequences that have < min_precent bases with quality < min_qual"""
    cmd = 'fastq_quality_filter -i %s -o %s -q %s -p %s' % (in_fastq, out_fastq,
                                                          min_qual, min_percent)
    sys_call(cmd)
if cfg.getboolean('filtering', 'filter_quality'):
    final_output = filter_min_quality
    prev_suffix = ''

@transform(final_output, suffix(''), '.read_lengths.png')
def read_length_histogram(in_fastq, out_hist):
    'draw histogram of read lengths'
    cmd = 'fasta_clipping_histogram.pl %s %s' % (in_fastq, out_hist)
    # lengths = Counter(imap(lambda l: len(l[1]),
    #                       parseFastq(open('./CHX_CA.fastq'))))
    # pyplot.bar(lengths.keys(), lengths.values())
    sys_call(cmd)

@transform(['*.fastq_illumina' if cfg.getboolean('filtering', 'convert_sanger_to_illumina') else '*.fastq', final_output],
        suffix(''), '.qual_stats')
def quality_stats(in_fastq, out_stats):
    'aggregate quality score statistics for reads'
    cmd = 'fastx_quality_stats -i %s -o %s' % (in_fastq, out_stats)
    sys_call(cmd)

@transform(quality_stats, suffix('.qual_stats'), '.qual_boxplot.png')
def quality_boxplot(in_stats, out_boxplot):
    'draw a boxplot for the quality scores'
    cmd = 'fastq_quality_boxplot_graph.sh -t %s -i %s -o %s' % (in_stats,
                                                in_stats, out_boxplot)
    sys_call(cmd)

@transform(quality_stats, suffix('.qual_stats'), '.qual_nuc_dist.png')
def quality_nuc_dist(in_stats, out_dist):
    'show the nucleotide distribution across the reads'
    cmd = 'fastx_nucleotide_distribution_graph.sh -t %s -i %s -o %s' % (
        in_stats, in_stats, out_dist)
    sys_call(cmd)

@merge(['*.fastq_illumina' if cfg.getboolean('filtering', 'convert_sanger_to_illumina') else '*.fastq',
        clip_adapter, trim_reads, trim_regex,
        filter_artifacts, filter_min_quality], 'fastq.wikisummary')
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
            for count, line in enumerate(open(infile)):
                pass
            outfile.write("| %s || %s\n|-\n" % (infile, count//4))
        outfile.write('|}\n')


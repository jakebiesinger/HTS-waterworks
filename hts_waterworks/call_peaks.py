"""call_peaks.py
    Module for calling peaks using a variety of peak callers
"""

#  Current Version: 0.0
#  Last Modified: 2011-07-22 16:52

import re
import itertools

from ruffus import (transform, follows, collate, files, split, merge,
                    add_inputs, regex, suffix, mkdir, jobs_limit, output_from)
from ruffus.task import active_if

from hts_waterworks.utils.ruffus_utils import (sys_call, main_logger as log,
                                           main_mutex as log_mtx)
from hts_waterworks.bootstrap import cfg
import hts_waterworks.mapping as mapping
from hts_waterworks.utils.common import (bedCommentFilter, readBedLines,
                                         parse_ucsc_range)


@active_if(cfg.getboolean('peaks', 'run_macs'))
@collate(mapping.all_mappers_output, regex(r'(.+)\.treat(.*)\.mapped_reads'), 
         add_inputs(r'\1.control\2.mapped_reads'), r'\1.treat\2.macs.peaks',
         cfg.getfloat('peaks', 'max_FDR'))
def run_macs(in_files, out_peaks, max_fdr):
    """Call peak with MACS (v1.3).
    Apply a maximum FDR threshold and treat centers as peak summits
    
    """
    in_treat, in_control = in_files[0]
    matches = re.search(r'(.*\.treat)(.*)\.mapped_reads', in_treat).groups()
    name = matches[0] + matches[1] + '.macs.peaks'
    max_fdr = cfg.getfloat('peaks', 'max_FDR')
    cmd = 'macs -t %s -c %s --name=%s %s' % (in_treat, in_control, name,
                                               cfg.get('peaks', 'macs_params'))
    sys_call(cmd)
    
    # convert to proper bedfile- ints for score and + for strand
    with open(out_peaks, 'w') as outfile:
        with open(name + '_peaks.xls') as infile:
            for index, line in enumerate(itertools.ifilter(
                                        bedCommentFilter, infile)):
                fields = line.split('\t')
                if fields[0] == 'chr':
                    continue # skip header
                start = str(max(0, int(fields[1])))
                score = str(max(0, min(1000, int(float(fields[6])))))
                fdr = float(fields[8])
                if fdr <= max_fdr:
                    outfile.write('\t'.join([fields[0], start, fields[2],
                                        'MACS_peak_%s' % (index + 1), score]) +
                                    '\t+\n')
    # take region surrounding the peak center as the summit
    summit_size = cfg.getint('peaks', 'peak_summit_size')
    with open(out_peaks + '_summits.%s_around' % \
                        cfg.get('peaks', 'peak_summit_size'), 'w') as outfile:
        with open(name + '_peaks.xls') as infile:
            for index, line in enumerate(itertools.ifilter(bedCommentFilter,
                                                                    infile)):
                fields = line.strip().split('\t')
                if fields[0] == 'chr':
                    continue # skip header
                score = str(max(0, min(1000, int(float(fields[6])))))
                p_start, p_stop = max(0, int(fields[1])), int(fields[2])
                p_center = p_start + (p_stop - p_start) / 2
                s_start = p_center - summit_size / 2
                s_stop = p_center + summit_size / 2
                fdr = float(fields[8])
                if fdr <= max_fdr:
                    outfile.write('\t'.join([fields[0], str(s_start),
                                    str(s_stop),
                                    'MACS_peak_%s' % (index + 1), score])
                                        + '\t+\n')

@active_if(cfg.getboolean('peaks', 'run_macs14'))
@collate(mapping.all_mappers_output, regex(r'(.+)\.treat(.*)\.mapped_reads'), 
         add_inputs(r'\1.control\2.mapped_reads'), r'\1.treat\2.macs14.peaks',
         cfg.getfloat('peaks', 'max_FDR'))
def run_macs14(in_files, out_peaks, max_fdr):
    """Call peaks using MACS (v1.4). Apply a maximum FDR threshold."""
    in_treat, in_control = in_files[0]
    matches = re.search(r'(.*\.treat)(.*)\.mapped_reads', in_treat).groups()
    name = matches[0] + matches[1] + '.macs14.peaks'
    cmd = 'macs14 -t %s -c %s --name=%s %s' % (in_treat, in_control, name,
                                             cfg.get('peaks', 'macs14_params'))
    sys_call(cmd)
    peaks_to_keep = set()
    # convert to proper bedfile- ints for score and + for strand
    with open(out_peaks, 'w') as outfile:
        with open(name + '_peaks.xls') as infile:
            for index, line in enumerate(itertools.ifilter(bedCommentFilter,
                                                                    infile)):
                fields = line.split('\t')
                if fields[0] == 'chr':
                    continue # skip header
                start = str(max(0, int(fields[1])))
                score = str(max(0, min(1000, int(float(fields[6])))))
                fdr = float(fields[8])
                if fdr <= max_fdr:
                    outfile.write('\t'.join([fields[0], start, fields[2],
                                        'MACS14_peak_%s' % (index + 1), score])
                                                + '\t+\n')
                    peaks_to_keep.add(index)
    # take region surrounding the peak summit
    summit_size = cfg.getint('peaks', 'peak_summit_size')
    with open(out_peaks + '_summits.%s_around' % \
                        cfg.get('peaks', 'peak_summit_size'), 'w') as outfile:
        with open(name + '_summits.bed') as infile:
            for index, line in enumerate(itertools.ifilter(bedCommentFilter,
                                                                    infile)):
                fields = line.strip().split('\t')
                if fields[0] == 'chr':
                    continue # skip header
                # score is number of reads at summit
                score = str(max(0, min(1000, int(float(fields[-1])))))
                start = str(max(0, int(fields[1]) - summit_size / 2))
                stop = str(int(fields[2]) + summit_size / 2)
                if index in peaks_to_keep:
                    outfile.write('\t'.join([fields[0], start, stop,
                                        'MACS_peak_%s' % (index + 1), score])
                                            + '\t+\n')

@active_if(cfg.getboolean('peaks', 'run_arem'))
@collate(mapping.all_mappers_output, regex(r'(.+)\.treat(.*)\.mapped_reads'), 
         add_inputs(r'\1.control\2.mapped_reads'), r'\1.treat\2.arem.peaks',
         cfg.getfloat('peaks', 'max_FDR'))
def run_arem(in_files, out_peaks, max_fdr):
    """Call peaks using AREM, applying a maximum FDR threshold"""
    in_treat, in_control = in_files[0]
    matches = re.search(r'(.*\.treat)(.*)\.mapped_reads', in_treat).groups()
    name = matches[0] + matches[1] + '.arem.peaks'
    cmd = 'arem -t %s -c %s --name=%s %s' % (in_treat, in_control, name,
                                               cfg.get('peaks', 'arem_params'))
    sys_call(cmd)
    # convert to proper bedfile- ints for score and + for strand
    peaks_to_keep = set()
    with open(out_peaks, 'w') as outfile:
        with open(name + '_peaks.xls') as infile:
            for index, line in enumerate(itertools.ifilter(bedCommentFilter,
                                                                    infile)):
                fields = line.split('\t')
                if fields[0] == 'chr':
                    continue # skip header
                start = str(max(0, int(fields[1])))
                score = str(max(0, min(1000, int(float(fields[6])))))
                fdr = float(fields[8])
                if fdr <= max_fdr:
                    outfile.write('\t'.join([fields[0], start, fields[2],
                                        'AREM_peak_%s' % (index + 1), score])
                                                + '\t+\n')
                    peaks_to_keep.add(index)
    # take region surrounding the peak summit
    summit_size = cfg.getint('peaks', 'peak_summit_size')
    with open(out_peaks + '_summits.%s_around' % \
                        cfg.get('peaks', 'peak_summit_size'), 'w') as outfile:
        with open(name + '_summits.bed') as infile:
            for index, line in enumerate(itertools.ifilter(bedCommentFilter,
                                                                    infile)):
                fields = line.strip().split('\t')
                if fields[0] == 'chr':
                    continue # skip header
                score = str(max(0, min(1000, int(float(fields[-1])))))
                start = str(max(0, int(fields[1]) - summit_size / 2))
                stop = str(int(fields[2]) + summit_size / 2)
                if index in peaks_to_keep:
                    outfile.write('\t'.join([fields[0], start, stop,
                                        'AREM_peak_%s' % (index + 1), score])
                                                + '\t+\n')


@active_if(cfg.getboolean('peaks', 'run_glitr'))
@transform(mapping.all_mappers_output,
    suffix('.mapped_reads'), '.mapped_reads_starts')
def bed_to_glitr(in_bed, out_starts):
    """Convert reads to (chrom, start, strand) for GLITR."""
    with open(in_bed) as infile:
        with open(out_starts, 'w') as outfile:
            for chrom, start, stop, strand in readBedLines(infile):
                outfile.write('\t'.join([chrom, str(start), strand]) + '\n')

@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@collate(bed_to_glitr, regex(r'(.+)\.treat\.mapped_reads_starts'), 
        add_inputs(r'\1.control.mapped_reads_starts'), r'\1.treat.glitr.ranges')
def run_glitr(in_files, out_peaks):
    """Call peaks with GLITR"""
    in_treat, in_control = in_files[0]
    glitr_dir = in_treat + '.GLITR_out'
    cmd = ('rm -r %s; mkdir %s; cd %s; GLITR.pl --CHIP=../%s ' + \
            '--CONTROL=../%s --GENOME=%s %s ') % (
                glitr_dir, glitr_dir, glitr_dir, in_treat, in_control,
                cfg.get('DEFAULT', 'genome').upper(),
                cfg.get('peaks', 'glitr_params'))
    sys_call(cmd)
    sys_call('cp %s/allChIP.FDR_1.5PercentFDR %s' % (glitr_dir, out_peaks))

@transform(run_glitr, suffix('.glitr.ranges'), '.glitr.peaks')
def glitr_range_to_bed(in_range, out_bed):
    """Convert GLITR ranges to BED format, use peak centers as summits"""
    summit_size = cfg.get('peaks', 'peak_summit_size')
    with open(in_range) as infile:
        with open(out_bed, 'w') as outfile:
            with open(out_bed + '_summits.%s_around' % summit_size, 'w') \
                                                            as outfile_summits:
                for i, line in enumerate(infile):
                    fields = line.strip('\n').split('\t')
                    chrom, start, stop = parse_ucsc_range(fields[0])
                    start = max(0, start)
                    foldchange = fields[3]
                    outfile.write('\t'.join([chrom, str(start), str(stop),
                                             'GLITR_peak_%s'%(i+1),
                                             str(int(float(foldchange))),'+'])
                                                + '\n')
                    # take bases around center as summit
                    center = start + (stop - start) / 2
                    center_start = center - summit_size / 2
                    center_stop = center + summit_size / 2
                    outfile_summits.write('\t'.join([chrom, str(center_start),
                                    str(center_stop), 'GLITR_peak_%s'%(i+1),
                                    str(int(float(foldchange))),'+']) + '\n')


all_peak_caller_functions = [run_macs, run_macs14, run_arem, glitr_range_to_bed]

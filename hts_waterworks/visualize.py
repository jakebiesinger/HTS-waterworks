"""visualize.py
    a module for visualizing read placement via bed, bedgraph, etc

"""

#  Current Version: 0.0
#  Last Modified: 2011-07-22 16:52

import tempfile
import re

from ruffus import (transform, follows, merge,
                    add_inputs, regex, suffix, jobs_limit)
from ruffus.task import active_if

from hts_waterworks.utils.ruffus_utils import sys_call, touch
from hts_waterworks.bootstrap import genome_path, cfg
import hts_waterworks.bootstrap as bootstrap
import hts_waterworks.mapping as mapping
import hts_waterworks.call_peaks as call_peaks


@jobs_limit(cfg.get('DEFAULT', 'max_throttled_jobs'), 'throttled')
@follows(bootstrap.get_chrom_sizes)
@transform(call_peaks.all_peak_caller_functions + mapping.all_mappers_output,
           suffix(''), '.clipped.sorted')
def bed_clip_and_sort(in_bed, out_sorted):
    """Sort the bed file and constrain bed regions to chromosome sizes"""
    with tempfile.NamedTemporaryFile() as tmp_clipped:
        cmd = 'bedClip %s %s.chrom.sizes %s' % (in_bed, genome_path(),
                                                tmp_clipped.name)
        sys_call(cmd)
        #cmd = 'bedSort %s %s' % (out_clipped, out_sorted)
        cmd = 'sort -k1,1 -k2,2n -S 2G %s > %s' % (tmp_clipped.name, out_sorted)
        sys_call(cmd)

@active_if(cfg.getboolean('visualization', 'uniquefy_track'))
@transform(bed_clip_and_sort, suffix(''), '.unique',
           cfg.getint('visualization', 'uniquefy_track_max_reads'))
def bed_uniquefy(in_bed, out_bed, max_reads):
    'Given a sorted bed file, remove tags that are on the same start, strand'
    with open(in_bed) as infile:
        with open(out_bed, 'w') as outfile:
            prev_start, prev_chrom = None, None
            plus_seen, minus_seen = 0, 0
            for line in infile:
                fields = line.split('\t')
                chrom, start, stop = fields[:3]
                if prev_start is None or prev_start != start or \
                                                       prev_chrom != chrom:
                    prev_start, prev_chrom = start, chrom
                    plus_seen, minus_seen = 0, 0
                if len(fields) < 6 or fields[5] == '+':
                    if plus_seen <= max_reads:
                        outfile.write(line)
                    plus_seen += 1
                else:
                    if minus_seen <= max_reads:
                        outfile.write(line)
                    minus_seen += 1

@transform([bed_uniquefy, bed_clip_and_sort], suffix(''), '.colored')
def bed_color_strand(in_bed, out_bed):
    """Color the regions of the BED file by thier strand"""
    with open(in_bed) as infile:
        with open(out_bed, 'w') as outfile:
            for line in infile:
                fields = line.strip('\n').split('\t')
                chrom, start, stop = fields[:3]
                strand = fields[5] if len(fields) >= 6 else '+'
                # +:RED, -:GREEN
                color = '255,0,0' if strand == '+' else '0,255,0'
                outfile.write('\t'.join(fields + [start, stop, color]) + '\n')


@transform(bed_color_strand, suffix(''), '.bigbed')
def bed_to_bigbed(in_bed, out_bigbed):
    """Convert a BED file to .bigbed for viewing on UCSC browser"""
    cmd = 'bedToBigBed %s %s.chrom.sizes %s' % (in_bed,
                                                genome_path(), out_bigbed)
    sys_call(cmd)

@transform([bed_uniquefy, bed_clip_and_sort],
    regex('(.*mapped_reads).clipped.sorted(.unique|)'),
    add_inputs(bootstrap.get_chrom_sizes),
    r'\1\2.bedgraph')
def bed_to_bedgraph(in_files, out_bedgraph):
    'extend reads to the full fragment length and create a bedgraph from them'
    in_bed, in_chrom_sizes = in_files
    cmd = ('slopBed -i %s -s -r %s -l 0 -g %s | ' + \
            'bedItemOverlapCount %s -chromSize=%s.chrom.sizes stdin > %s') % (
                        in_bed,
                        cfg.getint('DEFAULT','fragment_size') - \
                                            cfg.getint('DEFAULT','tag_size'),
                        in_chrom_sizes, cfg.get('DEFAULT', 'genome'),
                        genome_path(), out_bedgraph)
    sys_call(cmd)

@active_if(cfg.getboolean('visualization', 'normalize_per_million'))
@follows(bed_to_bedgraph)
@transform([bed_uniquefy, bed_clip_and_sort],
        regex('(.*mapped_reads).clipped.sorted(.unique|)'),
        add_inputs(r'\1\2.bedgraph'), r'\1\2.normalized.bedgraph')
def bedgraph_normalize_per_million(in_files, out_bedgraph):
    """Normalize bedgraph heights to tags per million mapping"""
    in_bed, in_bedgraph = in_files
    total_bed = sum(1 for l in open(in_bed))
    with open(in_bedgraph) as infile:
        with open(out_bedgraph, 'w') as outfile:
            for line in infile:
                fields = line.rstrip('\n').split('\t')
                norm_score = float(fields[-1]) / (total_bed / 1e6)
                outfile.write('\t'.join(fields[:-1] + [str(norm_score)]) + '\n')

#def bed_to_bedgraph_histogram(in_bed, out_bedgraph):
#    'create a histogram of read depth at a given window width'
#    with open(in_bed) as infile:
#        start, chrom = None, None
#        line = infile.next()
#        for line in infile:
            

@transform([bedgraph_normalize_per_million, bed_to_bedgraph],
    suffix('.bedgraph'), '.bigwig')
def bedgraph_to_bigwig(in_bedgraph, out_bigwig):
    """Convert the bedgraph file to .bigwig for viewing on UCSC"""
    cmd = 'bedGraphToBigWig %s %s.chrom.sizes %s' % (in_bedgraph, genome_path(),
                                                     out_bigwig)
    sys_call(cmd)

@merge([bedgraph_to_bigwig, bed_to_bigbed], 'ucsc.track.headers')
def make_track_headers(in_files, out_header):
    """For all the visualization files, create UCSC track headers"""
    with open(out_header,'w') as outfile:
        for in_track in in_files:
            if in_track.endswith('.bigwig'):
                track_extras = 'type=bigWig'
            elif in_track.endswith('.bigbed'):
                track_extras = 'type=bigBed itemRgb="On"'
            else:
                raise RuntimeError("Unrecognized file type: %s" % in_track)
            url = cfg.get('visualization', 'public_url_base') + '/' + in_track
            # remove cruft from the names
            short_name = re.sub(r'(mapped_reads|clipped|sorted|colored|\.)+',
                                                            ' ', in_track)
            track_str = 'track %s name="%s" description="%s" ' \
                        'bigDataUrl=%s\n' % (track_extras, short_name,
                                             short_name, url)
            outfile.write(track_str)

@follows(make_track_headers)
@merge([bedgraph_to_bigwig, bed_to_bigbed], 'ucsc.track.ready')
def deploy_track_files(in_files, out_header):
    """Copy UCSC tracks to public url"""
    remote = cfg.get('visualization', 'remote_ssh_dir')
    remote_host = remote.split(':')[0]
    remote_dir = remote.split(':')[1]
    for in_track in in_files:
        sys_call('ssh %s mkdir -p %s' % (remote_host, remote_dir))
        sys_call('scp %s %s' % (in_track, remote))
    touch(out_header)

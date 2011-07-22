
"""pas_seq.py
    Module for running poly-A site sequencing for alternative polyadenylation
    experiments.
"""

#  Current Version: 0.0
#  Last Modified: 2011-07-22 16:52

import itertools 

from ruffus import (transform, split, regex, suffix)
from ruffus.task import active_if


from hts_waterworks.utils.ruffus_utils import (
                                           sys_call, main_logger as log,
                                           main_mutex as log_mtx)
import hts_waterworks.mapping as mapping
from hts_waterworks.bootstrap import cfg




@active_if(cfg.getboolean('PAS-Seq', 'merge_adjacent_reads'))
@transform(mapping.all_mappers_output, suffix('.mapped_reads'),
           '.merged.mapped_reads',
           cfg.getint('PAS-Seq', 'merge_window_width'),
           cfg.getint('PAS-Seq', 'merge_num_iterations'))
def merge_adjacent_reads(in_bed, out_merged, window_width, iterations):
    """Reassign read ends to a weighted average of adjacent reads"""
    # helper functions for parsing bed files
    filter_lines = lambda l: l.strip() and (not l.startswith('#') or \
                                            l.startswith('"'))
    read_bed_lines = lambda infile: itertools.ifilter(filter_lines, infile)
    
    # sort the input by strand, chrom, stop
    tmpfile = in_bed + '.merged_adjacent_sorted'
    cmd = 'sort -k6 -k1 -k3g %s > %s' % (in_bed, tmpfile)
    print cmd
    sys_call(cmd, file_log=False)
    p_file = tmpfile
    
    for i in range(iterations):
        print 'merge iteration %s' % i
        # read in from output of previous iteration
        infile = read_bed_lines(open(p_file))
        
        # output to a temp file except on the last iteration
        if i != iterations - 1:
            p_file = in_bed + '.merge_adjacent_%s' % i
        else:
            p_file = out_merged
        outfile = open(p_file, 'w')

        # parse first line
        chrom, start, stop, name, score, strand = infile.next().split('\t')[:6]
        p_chrom, p_stops, p_names, p_strands = (chrom, [int(stop)],
                                                [name], [strand])
        print 'first line:', chrom, start, stop, name, score, strand
        
        for index, line in enumerate(infile):
            try:
                (chrom, start, stop,
                    name, score, strand) = line.rstrip('\n').split('\t')[:6]
            except:
                print index, 'this line:', line
                raise
            start, stop = int(start), int(stop)
            # is next read too far from first recorded?
            if p_chrom != chrom or (len(p_stops) > 0 and
                                    abs(p_stops[0] - stop) > window_width):
                if len(p_stops) == 0:
                    print 'error!'
                    print line
                    print p_stops, p_names, p_strands
                    raise
                avg = int(round(sum(p_stops) / float(len(p_stops))))
                # write out reads in this cluster, using avg as coordinate
                outfile.writelines('\t'.join([chrom, str(avg), str(avg+1),
                                         n_name, '0', n_strand]) + '\n'
                              for n_name, n_strand in zip(p_names, p_strands))
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
        outfile.close()

@split(merge_adjacent_reads, regex(r'(.*\.merged)\.mapped_reads$'),
       '\1.*.pileup', r'\1')
def pileup_strands(in_bed, out_pattern, out_base):
    "split the strands from the merged reads and give a pileup of the reads"
    # TODO finish up here
    pass

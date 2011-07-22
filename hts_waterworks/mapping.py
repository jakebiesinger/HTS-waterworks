
"""mapping.py
    module for mapping reads to the genome.
"""

#  Current Version: 0.0
#  Last Modified: 2011-07-22 16:52

import re
import random
import os
import tempfile

from ruffus import (transform, follows, files, split, merge, add_inputs,
                    regex, suffix, jobs_limit)
from ruffus.task import active_if
from pygr import worldbase

from hts_waterworks.utils.ruffus_utils import (sys_call, main_logger as log,
                                           main_mutex as log_mtx)
from hts_waterworks.bootstrap import genome_path, get_genome, cfg
import hts_waterworks.preprocessing as preprocessing


@active_if(cfg.getboolean('mapping', 'run_mosaik'))
@files(genome_path(), genome_path() + '.mosaik_dat')
def run_mosaik_build_reference(in_genome, out_bin):
    'convert reference to mosaik binary'
    cmd = 'MosaikBuild -fr %s -oa %s' % (in_genome, out_bin)
    sys_call(cmd)

mosaik_suffix_base = r'\1.mosaik_jump_%s' % cfg.getint('mapping',
                                                       'mosaik_hash_size')
@split(run_mosaik_build_reference, regex('(.*)\.mosaik_dat'),
       [mosaik_suffix_base + '_keys.jmp',
        mosaik_suffix_base + '_meta.jmp',
        mosaik_suffix_base + '_positions.jmp'], mosaik_suffix_base)
def run_mosiak_jump_reference(in_dat, _, out_jump_base):
    'create mosaik jump db on reference'
    cmd = 'MosaikJump -ia %s -out %s -hs %s' % (in_dat, out_jump_base,
                                    cfg.getint('mapping', 'mosaik_hash_size'))
    sys_call(cmd)


@active_if(cfg.getboolean('mapping', 'run_mosaik'))
@transform(preprocessing.final_output, suffix(''), '.mosaik_reads_dat')
def run_mosaik_build_reads(in_fastq, out_dat):
    'convert reads to mosaik binary'
    cmd = 'MosaikBuild -q %s -out %s -st illumina' % (in_fastq, out_dat)
    sys_call(cmd)

@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@transform(run_mosaik_build_reads, suffix('.mosaik_reads_dat'),
           add_inputs(run_mosaik_build_reference, run_mosiak_jump_reference),
           '.mosaik_align_dat')
def run_mosaik_align(in_files, out_align):
    'align reads to reference using MosaikAligner'
    # MosaikAligner -in sequence_archives/c_elegans_chr2_test.dat -out sequence_archives/c_elegans_chr2_test_aligned.dat -ia reference/c.elegans_chr2.dat -hs 14 -act 17 -mm 2 -m unique
    in_reads, in_genome_dat, in_genome_jump, _, _ = in_files
    in_genome_jump = in_genome_jump.replace('_keys.jmp', '')
    cmd = 'MosaikAligner -in %s -ia %s -j %s -out %s -hs %s  %s'
    cmd = cmd % (in_reads, in_genome_dat, in_genome_jump, out_align,
                   cfg.getint('mapping', 'mosaik_hash_size'),
                   cfg.get('mapping', 'mosaik_params'))
    sys_call(cmd)

@transform(run_mosaik_align, suffix('.mosaik_align_dat'), '.mosaik_align_sam')
def mosaik_to_sam(in_dat, out_sam):
    'convert mosaik alignments to SAM format'
    cmd = 'MosaikText -in %s -sam %s' % (in_dat, out_sam)
    sys_call(cmd)
    cmd = 'gunzip %s.gz' % out_sam
    sys_call(cmd)

@transform(mosaik_to_sam, suffix('.mosaik_align_sam'), '.mosaik.mapped_reads')
def mosaik_to_bed(in_sam, out_bed):
    'convert mosaik alignments to BED format'
    cmd = 'sam2bed %s %s' % (in_sam, out_bed)
    sys_call(cmd)


@active_if(cfg.getboolean('mapping', 'run_pash'))
@transform(preprocessing.final_output, suffix(''), '.pash_reads')
def run_pash(in_fastq, out_pash):
    'align reads using PASH'
    cmd = 'pash-3.0lx.exe -h %s -v %s -o %s  %s '
    cmd = cmd % (genome_path(), in_fastq, out_pash,
                        cfg.get('mapping', 'pash_params'))
    sys_call(cmd)

@transform(run_pash, suffix('.pash_reads'), '.pash.mapped_reads')
def pash_to_bed(in_pash, out_bed):
    'convert PASH output to BED format'
    with open(in_pash) as infile:
        with open(out_bed, 'w') as outfile:
            for line in infile:
                fields = line.rstrip('\n').split('\t')
                outfile.write('\t'.join(fields[:3] +
                                        ['0', '0', fields[6]]) + '\n')


@active_if(cfg.getboolean('mapping', 'run_bowtie'))
@split(None, '%s.*.ebwt' % genome_path())
def get_bowtie_index(_, out_pattern):
    'retrieve or build bowtie index'
    cmd = 'wget -N ftp://ftp.cbcb.umd.edu/pub/data/bowtie_indexes/%s.ebwt.zip'
    cmd = cmd % cfg.get('DEFAULT', 'genome')
    try:
        sys_call(cmd)
        cmd = 'unzip -o %s.ebwt.zip -d %s' % (cfg.get('DEFAULT', 'genome'),
                                              os.path.split(genome_path())[0])
        sys_call(cmd)
    except:
        # downloading the file failed, so we have to build the index ourselves
        with log_mtx:
            log.info('Bowtie index download failed. '
                     'Trying to build a new index')
        cmd = 'bowtie-build %s %s' % (genome_path(), genome_path())
        sys_call(cmd)

@active_if(cfg.getboolean('mapping', 'run_bowtie'))
@follows(get_bowtie_index)
@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@transform(preprocessing.final_output, suffix(''), '.bowtie_reads')
def run_bowtie(in_fastq, out_bowtie):
    'align reads to reference using Bowtie'
    cmd = 'bowtie %s %s %s %s' % (genome_path(),
                                  cfg.get('mapping', 'bowtie_params'),
                                  in_fastq, out_bowtie)
    sys_call(cmd)

@transform(run_bowtie, suffix('.bowtie_reads'), '.bowtie.mapped_reads')
def bowtie_to_bed(in_bowtie, out_bed):
    'Convert Bowtie-formatted output into bed format'
    with open(in_bowtie) as infile:
        # use first ten reads to determine read length
        read_lengths = [len(infile.readline().split('\t')[4])
                                            for i in range(10)]
        read_lengths = sum(read_lengths) / len(read_lengths)
        infile.seek(0)
        with open(out_bed, 'w') as outfile:
            for line in infile:
                name, strand, chrom, start = line.split('\t')[:4]
                stop = int(start) + read_lengths + 1  # stop is fencepost after
                outfile.write('\t'.join([chrom, start, str(stop), name, '0',
                                         strand]) + '\n')


# Mapping with ssaha2
@active_if(cfg.getboolean('mapping', 'run_ssaha2'))
@split(genome_path(), suffix(''),
       [r'\1.body', r'\1.base', r'\1.head', r'\1.name', r'\1.size'])
def get_ssaha2_hashtable(in_genome, out_ssaha2):
    """Use ssaha2Build to generate a hash table for the genetic
       sequences stored in an input .fasta file
        
       ssaha2Build writes five files to disk, each preceded by
       the hash name. Their file extensions are:
       .base, .body, .head, .name, .size
    """
    
    #TODO: add useful parameters to cmd and config file
    cmd = 'ssaha2Build -save %s %s' % (out_ssaha2, in_genome)
    sys_call(cmd)

''' LOG
* Added hash_name to cfg
* removed sentinel file in get_ssaha2_hashtable
* Cleaned up output file names for run_ssaha2
*
'''

@active_if(cfg.getboolean('mapping', 'run_ssaha2'))
@follows(get_ssaha2_hashtable)
@transform(preprocessing.final_output, suffix(''), '.ssaha2_reads')
def run_ssaha2(in_fastq, out_ssaha2):
    """ Runs ssaha2 command using the prebuilt hash table from
        get_ssaha2_hashtable. 
    
        The ssaha2 command maps DNA sequence reads onto a genomic 
        reference sequence using a combination of word hashing and 
        dynamic programming. (From ssaha2 manual)
    """
    #TODO: add useful parameters to cmd and config file
    #cmd = 'ssaha2 -outfile %s -save %s %s' % (out_ssaha2, hash_name, in_fastq)
    cmd = 'ssaha2 -outfile %s -disk 1 -save %s %s'
    cmd = cmd % (out_ssaha2, (cfg.get('mapping', 'ssaha2_hash_name')), in_fastq)
    sys_call(cmd)

# Mapping with Maq
@active_if(cfg.getboolean('mapping', 'run_maq'))
@files(genome_path(), genome_path() + '.maq.bfa')
def maq_index_reference(in_fasta, out_bfa):
    """ Use maq fasta2bfa to convert reference sequences in .fasta
    format to BFA format, which is a binary representation.
    """
    cmd = 'maq fasta2bfa %s %s' % (in_fasta, out_bfa)
    sys_call(cmd)

@active_if(cfg.getboolean('mapping', 'run_maq'))
@follows(maq_index_reference) 
@transform(preprocessing.final_output, suffix(''), '.maq.bfq') 
def maq_index_reads(in_fastq, out_bfq):
    """ Use maq fastq2bfq to convert read sequences in .fastq
    format to BFQ format, which is a binary representation.
    """
    cmd = 'maq fastq2bfq %s %s' % (in_fastq, out_bfq)
    sys_call(cmd)
    
# TODO: Can merge preserve file names?
@active_if(cfg.getboolean('mapping', 'run_maq'))
@merge([maq_index_reference, maq_index_reads], '.maq.binary.map')
def maq_map_reads(in_files, out_map):
    """ Use maq match to align the reads to the reference.
    Input files are in binary format and output is in .map format.
    """
    cmd = 'maq match %s %s %s' % (out_map, in_files[0], in_files[1])
    sys_call(cmd)
    

@transform(maq_map_reads, suffix('binary.map'), '.readable.map')
def maq_view_reads(in_map, out_map):
    """ Use maq mapview to generate a human readable .map
    format.
    """
    cmd = 'maq mapview %s > %s' % (in_map, out_map)
    sys_call(cmd)
    
@transform(maq_view_reads, suffix('readable.map'), '.bed')
def maq_map_to_bed(in_map, out_bed):
    """ Convert maq map file to BED format """
    with open(in_map) as infile:
        # use first ten reads to determine read length
        read_lengths = [len(infile.readline().split('\t')[14])
                                                    for i in range(10)]
        read_lengths = sum(read_lengths) / len(read_lengths)
        infile.seek(0)
        with open(out_bed, 'w') as outfile:
            for line in infile:
                fields = line.strip().split('\t')
                chrom, start, strand = fields[1], fields[2], fields[3]
                name = cfg.get('mapping', 'maq_bed_name')
                score = cfg.get('mapping', 'maq_bed_score')
                stop = int(start) + read_lengths + 1  # stop is fencepost after
                outfile.write('\t'.join([chrom, str(start), str(stop),
                                         str(name), str(score), str(strand)])
                              + '\n')


all_mappers_output = [bowtie_to_bed, pash_to_bed, mosaik_to_bed, run_ssaha2,
                      maq_map_to_bed]


@active_if(cfg.getboolean('peaks', 'downsample_reads'))
@split(all_mappers_output, regex(r'(.+)\.treat(.*)\.mapped_reads'),
       add_inputs(r'\1.control\2.mapped_reads'),
       [r'\1.treat\2.matched_size_*.mapped_reads',
        r'\1.control\2.matched_size_*.mapped_reads'])
def uniquefy_downsample_reads(in_files, out_files):
    """Uniquefy sequence reads then downsample so the total unique tag count in
    treatment and control is the same.  This may generate many downsampled datasets.
    """
    # WARNING: this is a circular dependency.  It has to be included at runtime
    #    Top-level import will cause this module to load only 1/2 way
    #    we import here because we need to call this function directly,
    #    and not just when using ruffus
    from hts_waterworks.visualize import bed_clip_and_sort, bed_uniquefy
    if not cfg.getboolean('peaks', 'downsample_reads'):
        with log_mtx:
            log.debug('NOT downsampling the sequence reads!')
    else:
        in_treat, in_control = in_files
        out_treat_template = re.sub(r'mapped_reads$',
                                    'matched_size_%s.mapped_reads', in_treat)
        out_control_template = re.sub(r'mapped_reads$',
                                    'matched_size_%s.mapped_reads', in_control)
        if out_treat_template == in_treat:
            raise RuntimeError('regex substitution failed from %s to %s' % (
                                                in_treat, out_treat_template))
        if out_control_template == in_control:
            raise RuntimeError('regex substitution failed from %s to %s' % (
                                            in_control, out_control_template))
        
        # sort the reads
        bed_clip_and_sort(in_treat, in_treat+'.tmp_sorted')
        bed_clip_and_sort(in_control, in_control+'.tmp_sorted')
        
        # uniquefy the reads
        bed_uniquefy(in_treat+'.tmp_sorted', in_treat+'.tmp_unique',
                     cfg.getint('visualization', 'uniquefy_track_max_reads'))
        bed_uniquefy(in_control+'.tmp_sorted', in_control+'.tmp_unique',
                     cfg.getint('visualization', 'uniquefy_track_max_reads'))
        
        total_treat = sum(1 for l in open(in_treat+'.tmp_unique'))
        total_control = sum(1 for l in open(in_control+'.tmp_unique'))
        if total_treat == total_control:
            with log_mtx:
                log.debug('No downsampling required-- tag counts identical')
        else:
            # downsample num_down_sample times
            for i in xrange(cfg.getint('peaks', 'num_down_samples')):
                out_treat = out_treat_template % i
                out_control = out_control_template % i
                if total_treat > total_control:
                    # reduce number of treatment reads
                    inds_to_keep = set(random.sample(xrange(total_treat),
                                                                total_control))
                    in_orig, out_orig = in_control+'.tmp_unique', out_control
                    in_subset, out_subset = in_treat+'.tmp_unique', out_treat
                else:
                    # reduce number of control reads
                    inds_to_keep = set(random.sample(xrange(total_control),
                                                     total_treat))
                    in_orig, out_orig = in_treat+'.tmp_unique', out_treat
                    in_subset, out_subset = in_control+'.tmp_unique', out_control
                sys_call('cp %s %s' % (in_orig, out_orig))
                # subset the tags
                with open(in_subset) as infile:
                    with open(out_subset, 'w') as outfile:
                        outfile.writelines(line for i, line in enumerate(infile) 
                                                        if i in inds_to_keep)
if cfg.getboolean('peaks', 'downsample_reads'):
    all_mappers_output.append(uniquefy_downsample_reads)


@active_if(cfg.getboolean('mapping', 'collapse_mapped_read_IDs'))
@transform(all_mappers_output, suffix('.mapped_reads'),
           '.collapsedID.mapped_reads')
def collapse_mapped_read_IDs(in_bed, out_bed):
    """Reads with the same ID and start,strand are collapsed into a single read
    The ID is defined as any characters following the final underscore in the
    4th column of the BED file ("name")
    """
    with tempfile.NamedTemporaryFile() as tmp_sorted:
        # sort the file by chrom, start positions
        cmd = 'sort -k1,1 -k2,2n -S 2G %s > %s' % (in_bed, tmp_sorted.name)
        sys_call(cmd)
        with open(out_bed, 'w') as outfile:
            p_chrom, p_start, p_name, p_strand = None, None, None, None
            for line in open(tmp_sorted.name):
                chrom,start,stop,name,score,strand = line.split('\t')
                name = name.split('_')[-1] if '_' in name else name
                if not (p_name == name and p_start == start and \
                        p_strand == strand and p_chrom == chrom):
                    # current read does NOT have the same position or ID
                    outfile.write(line)
                p_chrom, p_start, p_name, p_strand = chrom, start, name, strand
if cfg.getboolean('mapping', 'collapse_mapped_read_IDs'):
    all_mappers_output = [collapse_mapped_read_IDs]


@follows(get_genome)
@active_if(cfg.getboolean('mapping', 'remove_internal_priming'))
@transform(all_mappers_output, suffix('.mapped_reads'),
           '.no_internal_priming.mapped_reads')
def remove_internal_priming(in_bed, out_bed):
    """Reads that map to genomic locations with 7+ A's in the 10nt downstream of
    the end of the read should be filtered out
    """
    wb_genome = worldbase(cfg.get('DEFAULT', 'worldbase_genome'))
    with open(out_bed, 'w') as outfile:
        for line in open(in_bed):
            chrom,start,stop,name,score,strand = line.split('\t')
            start, stop = int(start), int(stop)
            if strand == '+':
                try:
                    seq = str(wb_genome[chrom][start:stop])
                except IndexError:
                    seq = ''
                seq_A = seq.count('AAAAAA')
                try:
                    downstream = str(wb_genome[chrom][stop:stop+10])
                except IndexError:
                    downstream = ''
                down_A = downstream.count('A')
            else:
                try:
                    seq = str(wb_genome[chrom][start:stop])
                except IndexError:
                    seq = ''
                seq_A = seq.count('TTTTTT')
                try:
                    downstream = str(wb_genome[chrom][max(0,start-10):start])
                except IndexError:
                    downstream = ''
                down_A  = downstream.count('T')
            #filter if 6+ consecutive A's in sequence or 7+ A's downstream
            if seq_A < 1 and down_A < 7:
                outfile.write(line)
if cfg.getboolean('mapping', 'remove_internal_priming'):
    all_mappers_output = [remove_internal_priming]

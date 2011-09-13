"""bootstrap.py
    module to set up a pipeline program
"""

#  Current Version: 0.1-4-g6334e5a
#  Last Modified: 2011-09-12 20:18

import StringIO
from ConfigParser import ConfigParser

from ruffus import files
from pygr import worldbase

from hts_waterworks.utils.ruffus_utils import (touch, sys_call,
                                           main_logger as log,
                                           main_mutex as log_mtx)

# --- config and options---
CONFIG_FILENAME = 'pipeline.cfg'
CONFIG_TEMPLATE = """
[DEFAULT]
# Lab name and project name are used for uploading files to a webserver
lab_name = test
project_name = pipeline

# Max number cpu or memory heavy jobs to run simultaneously
max_throttled_jobs = 3
# Maximum number of threads per job
max_threads = 3
WORLDBASEPATH=.

worldbase_genome = Bio.Seq.Genome.HUMAN.%(genome)s
genome = hg19
R_annotation = org.Hs.eg.db
R_mart = hsapiens_gene_ensembl

#worldbase_genome = Bio.Seq.Genome.MOUSE.%(genome)s
#genome = mm9
#R_annotation = org.Mm.eg.db
#R_mart = mmusculus_gene_ensembl

#worldbase_genome = Bio.Seq.Genome.XENTR.%(genome)s
#genome = xenTro2
#R_annotation = org.Xl.eg.db
#R_mart = xtropicalis_gene_ensembl

fragment_size = 200
tag_size = 50


[filtering]
convert_sanger_to_illumina = False

clip_adapter = False
adapter_sequence =

trim_reads = False
trim_start = 1
trim_end = 50


# trim everything to the right of this regular expression, including the re match
# leave blank if no trimming is to take place
trim_regex =

filter_artifacts = False

filter_quality = True
filter_min_quality = 20
filter_percent_bases_at_min = 90


[mapping]
run_bowtie = True
# Keep only uniquely mapping reads
bowtie_params = -n 2 -m 1 --phred64-quals -s 1 -p %(max_threads)s
# Keep up to 10 alignments for each sequence read
#bowtie_params = -n 2 -k 10 --best --phred64-quals -s 1 -p %(max_threads)s

run_tophat = False
tophat_params = -p %(max_threads)s

run_pash = False
pash_params = -G 6 -k 13 -n 21 -s 30 -d 800 -S .

run_mosaik = False
mosaik_hash_size = 15
mosaik_params = -pd -act 17 -mm 4 -mhp 100 -p %(max_threads)s

run_ssaha2 = False
ssaha2_hash_name = htab

run_maq = False
maq_bed_name = 0
maq_bed_score = 0


collapse_mapped_read_IDs = False
remove_internal_priming = False

map_to_transcriptome = False
separate_by_strand = False

[peaks]
downsample_reads = True
num_down_samples = 1

# All peaks reported must be <= this FDR
max_FDR = 5

run_macs = False
macs_params = --tsize=%(tag_size)s

run_macs14 = False
macs14_params = -g hs

run_arem = False
arem_params = -g hs

run_glitr = False
glitr_params = --NUMCPUS=%(max_threads)s --GLITR=/opt --FRAGLEN=%(fragment_length) --SEQLEN=%(tag_size)s --NN=70 --KNN=100 --FDR=%(max_FDR)s

run_quest = False
quest_params =
quest_directory = ''
quest_params_generate_step1 = y
quest_params_generate_step2 = 1
quest_params_generate_step3 = 2
quest_params_generate_step4 = y

peak_summit_size = 50


[visualization]
uniquefy_track = True
# maximum number of reads to show at each position
uniquefy_track_max_reads = 2

normalize_per_million = True

remote_ssh_dir = localhost:~/public_html/projects/%(lab_name)s/%(project_name)s
public_url_base = http://www.ics.uci.edu/~wbiesing/projects/%(lab_name)s/%(project_name)s/


[motifs]
run_meme = True
meme_params = -dna -minw 7 -maxw 20 -mod zoops -nmotifs 3 -p %(max_threads)s -maxsize 10000000
run_nmica = True
nmica_params = -minLength 7 -maxLength 20 -numMotifs 3 -maxCycles 5000 -revComp -threads %(max_threads)s

# number of sampled chunks to run motif discovery against
motif_num_chunks = 1
# max number of sequences input to motif discovery
motif_chunk_size = 500
# genomic samples for threshold score
motif_threshold_sample_size = 500000
# genomic samples for significance estimate
motif_significance_sample_size = 100000
sampling_exclude_repeats = False
sampling_exclude_N = True

motif_zscores = 4.29


[genes]
nearby_genes_max_dist = 20000

promoter_size = 2000
promoter_extend = 0
downstream_size = 500
downstream_extend = 0

gene_overlap_order = ['promoter', 'utr5', 'exon', 'intron', 'utr3', 'down', 'intergenic']
gene_profile_on_raw_reads = False

download_refseq = True

ks_test_default_value = 1.0

[PAS-Seq]
merge_adjacent_reads = False
merge_window_width = 40
merge_num_iterations = 2

compare_window_width = 24
"""

# global configuration
cfg = ConfigParser()
cfg.readfp(StringIO.StringIO(CONFIG_TEMPLATE))  # set up defaults
cfg.read(CONFIG_FILENAME)  # add user-specified settings

def genome_path():
    'returns the path to the genome fasta file (and downloads it if necessary)'
    genome = worldbase(cfg.get('DEFAULT', 'worldbase_genome'), download=True)
    return genome.filepath


@files(None, genome_path())
def get_genome(_, out_genome_path, touch_file=True):
    'download the worldbase genome'
    genome = worldbase(cfg.get('DEFAULT', 'worldbase_genome'), download=True)
    if touch_file:
        touch(out_genome_path)
    return genome


@files(None, '%s.chrom.sizes' % genome_path())
def get_chrom_sizes(_, out_sizes):
    'retrieve the chromosome sizes for the current genome from UCSC'
    cmd = 'fetchChromSizes %s > %s' % (cfg.get('DEFAULT', 'genome'), out_sizes)
    sys_call(cmd, file_log=False)

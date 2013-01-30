#!/usr/bin/env python
"""waterworks.py
   Run a chip-seq pipeline, with output in the current directory.
   Place the raw reads at *.treat and *.control

"""

# --- imports ---
from ruffus import follows, files

from hts_waterworks.utils.ruffus_utils import (ruffus_main, ruffus_opt_parser,
                                               touch)
import hts_waterworks.preprocessing as preprocessing
import hts_waterworks.mapping as mapping
import hts_waterworks.call_peaks as call_peaks
import hts_waterworks.visualize as visualize
import hts_waterworks.annotation as annotation
import hts_waterworks.pas_seq as pas_seq
import hts_waterworks.motif_discovery as motif_discovery
import hts_waterworks.clip_seq as clip_seq
from hts_waterworks.bootstrap import CONFIG_FILENAME, CONFIG_TEMPLATE


@files(None, CONFIG_FILENAME)
def bootstrap(_, out_config):
    'create a configuration template'
    with open(out_config, 'w') as outfile:
        outfile.write(CONFIG_TEMPLATE)

@follows(preprocessing.final_output, preprocessing.read_length_histogram,
         preprocessing.quality_boxplot, preprocessing.quality_nuc_dist,
         preprocessing.summarize_fastq_reads)
def preprocess():
    """preprocessing/filtering steps complete"""
    pass


@follows(*(mapping.all_mappers_output + [mapping.summarize_mapped_reads]))
def map_reads():
    """read mapping steps complete"""
    pass


@follows(*call_peaks.all_peak_caller_functions)
def peak_calling():
    """peak calling steps complete"""
    pass


# visualization
@follows(visualize.deploy_track_files)
def visualization():
    """visualization steps complete"""
    pass


@follows(annotation.gene_overlap, annotation.gene_ontology,
         annotation.find_nearby_genes, annotation.plot_nearest_features,
         annotation.draw_expression_correlation, annotation.draw_expression_ks)
def expression():
    """expression/annotation steps complete"""
    pass


@follows(motif_discovery.motif_enrichment_genomic,
         motif_discovery.motif_enrichment_control,
         motif_discovery.consensus_enrichment,
         motif_discovery.motif_presence_sorted_peaks,
         motif_discovery.make_seq_logo
         #motif_discovery.run_stamp,
         )
def motifs():
    """motif discovery complete"""
    pass

@follows(pas_seq.plot_closest_polyA_db, pas_seq.plot_ttest_polya,
         pas_seq.plot_scatter_polya, pas_seq.plot_differential_polya,
         pas_seq.intersect_comparison_types, pas_seq.test_differential_polya)
def do_pas_seq():
    """Run complete pas_seq pipeline"""

@follows(clip_seq.reproducible_motifs, clip_seq.reproducible_positions,
         clip_seq.find_meme_motifs_around_sites)
def do_clip_seq():
    """Run complete CLIP-seq pipeline"""

@follows(preprocess, map_reads, peak_calling, visualization, expression, motifs,
         do_pas_seq, do_clip_seq)
@files(None, 'all_complete.ready')
def all_complete(_, out_sentinel):
    """HTS workflow complete"""
    touch(out_sentinel)


# --- main ---
if __name__ == '__main__':
    parser = ruffus_opt_parser()
    opts, args = parser.parse_args()
    ruffus_main(opts, args)

"""motif_discovery.py
    A module for performing motif discovery on called peaks.
    
    Discovery performed with MEME or NestedMICA, with plots and significance.

"""


import shlex
import glob
import re
import pickle
import random
from array import array
from collections import defaultdict
import math

from pygr import worldbase
import numpy
from matplotlib import pyplot
import matplotlib
from ruffus import (transform, suffix, files, follows, regex, split, collate,
                    add_inputs, jobs_limit)
from ruffus.task import active_if

from hts_waterworks.utils.ruffus_utils import (sys_call,
                                               main_logger as log)
from hts_waterworks.utils.common import (readBedLines,
                                         parseFastaLines,
                                         fastaToSequenceList)
import hts_waterworks.utils.get_bed_sequence as get_bed_sequence
import hts_waterworks.utils.sequence_motif as sequence_motif
import hts_waterworks.utils.sampling as sampling
import hts_waterworks.utils.motif_significance as motif_significance
from hts_waterworks.bootstrap import cfg, get_genome, genome_path
import hts_waterworks.call_peaks as call_peaks
import hts_waterworks.annotation as annotation


#from ipdb import set_trace as breakpoint

# motif setup

@transform(call_peaks.all_peak_caller_functions + 
          ['*.peaks_summits.%s_around' % cfg.get('peaks', 'peak_summit_size')],
        regex(r'(.*\.peaks$|.*\..*_around$|_genes.promoter.*_ext[\d]+$)'),
        r'\1.top%s.peaks' % cfg.getint('motifs', 'motif_chunk_size'),
        cfg.getint('motifs', 'motif_chunk_size'))
def get_top_peaks(in_peaks, out_subset, num_peaks_to_keep):
    """keep only the top peaks as input to motif discovery"""
    with open(in_peaks) as infile:
        seqs = list(readBedLines(infile, dataOnly=False))
        # sort by score, highest first
        seqs.sort(key=lambda x: int(x[4]), reverse=True)
        with open(out_subset, 'w') as outfile:
            subset = seqs[:num_peaks_to_keep]
            outfile.writelines('\t'.join(map(str, s)) + '\n' for s in subset)

#@follows(get_genome)
@transform([get_top_peaks], suffix(''), '.fasta')
def get_peak_sequence(in_peaks, out_fasta):
    """Get fasta file for peak summits
    """
    in_summits = out_fasta.replace('.fasta', '')
    args = shlex.split('''--genome=%s %s %s''' % (
                                        cfg.get('DEFAULT', 'worldbase_genome'),
                                        in_summits, out_fasta))
    get_bed_sequence.main(args)

@split(get_peak_sequence, regex(r'(.+peaks_summits.+).fasta'), r'\1.small_sample.*.fasta')
def motif_select_random_seqs(in_fasta, out_pattern):
    """Split a fasta file into several chunks so motif discovery is easier"""
    name = name = re.search('(.*).fasta', in_fasta).groups()[0]
    with open(in_fasta) as infile:
        seqs = list(parseFastaLines(infile))
        if len(seqs) <= cfg.get('motifs', 'motif_chunk_size'):
            num_chunks = 1
        else:
            num_chunks = cfg.get('motifs', 'motif_num_chunks')
        # get a random sample of peaks
        for i in xrange(num_chunks):
            with open(name + '.small_sample.%s.fasta' % i, 'w') as outfile:
                subset = random.sample(seqs, min(len(seqs),
                                    cfg.getint('motifs', 'motif_chunk_size')))
                outfile.writelines('>%s\n%s\n' % (s[0].strip(), s[1].strip())
                                                                for s in subset)

# motif discovery
@active_if(cfg.getboolean('motifs', 'run_meme'))
@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@transform(motif_select_random_seqs,
           #suffix('.fasta'), '.meme.discovered.motifs')
           #regex(r'(.*(?=_around).*(?=top).*).fasta$'),
           regex(r'(.*(?=top).*).fasta$'),
           r'\1.meme.discovered.motifs')
def discover_meme_motifs(in_fasta, out_motifs):
    """Discover sequence motifs in peaks by running meme"""
    cmd = 'meme %s %s -oc %s_meme_out ' % (in_fasta,
                                           cfg.get('motifs', 'meme_params'),
                                           out_motifs)
    #if 'top' in in_fasta and 'around' in in_fasta:
    sys_call(cmd)
    motifs = sequence_motif.parseMemeMotifs('%s_meme_out/meme.txt' % out_motifs)
    pickle.dump(motifs, open(out_motifs, 'w'))

@active_if(cfg.getboolean('motifs', 'run_nmica'))
@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@transform(motif_select_random_seqs, suffix('.fasta'),
           '.nmica.discovered.motifs')
def discover_nmica_motifs(in_fasta, out_motifs):
    """Discover sequence motifs in peaks by running nestedMICA"""
    cmd = 'nminfer -seqs %s %s ' % (in_fasta, cfg.get('motifs', 'nmica_params'))
    sys_call(cmd)
    motifs_name = in_fasta.replace('.fasta', '.motifs.xms')
    sys_call('mv motifs.xms %s' % motifs_name)
    motifs = sequence_motif.parse_xms_motifs(motifs_name)
    pickle.dump(motifs, open(out_motifs, 'w'))
    #args = shlex.split('%s %s' % (motifs_name, out_motifs))
    #parse_nmica_motifs.main(args)

# motif enrichment
@follows(get_genome)
@files(None, '%s.genome_samples.size30.num%s.fasta' % (genome_path(),
                            cfg.get('motifs', 'motif_threshold_sample_size')))
def sample_genome_short(_, out_samples):
    """Genomic sampling for threshold score"""
    args = shlex.split('''%s --genome=%s --sample_length=30 --num_samples=%s
                       ''' % (out_samples,
                              cfg.get('DEFAULT', 'worldbase_genome'),
                              cfg.get('motifs', 'motif_threshold_sample_size')))
    sampling.main(args)

@transform('*.known.motifs.transfac', suffix('.transfac'), '')
def convert_transfac_motifs(in_transfac, out_pickle):
    """Convert text files with motifs into our pickled format"""
    transfac_str = open(in_transfac).read()
    m = sequence_motif.parseMotifsFromTransfac(transfac_str)
    pickle.dump(m, open(out_pickle, 'w'))

@transform('*.known.motifs.dreme', suffix('.dreme'), '')
def convert_dreme_motifs(in_dreme, out_pickle):
    """Convert text files with motifs into our pickled format"""
    dreme_str = open(in_dreme).read()
    m = sequence_motif.parseMotifsFromTransfac(dreme_str)
    pickle.dump(m, open(out_pickle, 'w'))
    
@follows(convert_transfac_motifs)
@transform([discover_meme_motifs, discover_nmica_motifs, '*.known.motifs'],
        suffix('.motifs'),
        add_inputs(sample_genome_short), '.with_mean_sd.motifs')
def motif_mean_sd(in_files, out_motifs):
    """calculate the motifs' background score distributions, with mean and sd.
    
    """
    in_motif, genome_samples = in_files
    # Convert the .fasta file to a list of sequences
    sequenceList = fastaToSequenceList(genome_samples)

    with open(in_motif) as infile:
        all_motifs = pickle.load(infile)

    # maybe need to create new empty mutable dictionary
    # and populate it with items in the for loop

    for motif in all_motifs.values():
        motif.calculate_background(None, samples=sequenceList)

    # Save the motifs that were in the motifs input file
    # genome_samples are in .fasta format
    with open(out_motifs, 'w') as outfile:
        pickle.dump(all_motifs, outfile)

@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@split(call_peaks.all_peak_caller_functions + [get_top_peaks, '*.custom.peaks'],
        regex(r'(.*)\.peaks$'),
        [r'\1.peaks.similar.genomic.sample',
         r'\1.peaks.similar.genomic.sample.locations'])
def sample_genome_like_peaks(in_peaks, out_files):
    """Sample from the genome, keeping the sample widths the same as peaks"""
    out_sample, out_locations = out_files[:2]
    peak_lengths = array('i', (stop - start for chrom, start, stop, strand in
                        readBedLines(open(in_peaks))))
    if len(peak_lengths) == 0:
        raise RuntimeError("Peaks file %s is empty!" % in_peaks)
    
    wb_genome = worldbase(cfg.get('DEFAULT', 'worldbase_genome'))
    s = sampling.sample_genome(wb_genome, peak_lengths,
                               sampleSize=cfg.getint('motifs',
                                        'motif_significance_sample_size'),
                               excludeRepeat=cfg.getboolean('motifs',
                                                'sampling_exclude_repeats'),
                               excludeN=cfg.getboolean('motifs',
                                                'sampling_exclude_N'),
                               ignoreCharacters='_', weighted=True)
    with open(out_sample, 'w') as outfile:
        with open(out_locations, 'w') as outlocations:
            for index, line in enumerate(s):
                outfile.write('>%s\n%s\n' % (index, line))
                outlocations.write('\t'.join([line.id, str(line.start),
                                             str(line.stop), str(index), '0',
                                '+' if line.orientation == 1 else '-']) + '\n')

@active_if(len(glob.glob('*control*fastq')) > 0)
@jobs_limit(cfg.getint('DEFAULT', 'max_throttled_jobs'), 'throttled')
@split(call_peaks.all_peak_caller_functions + [get_top_peaks],
        regex(r'(.*).peaks$'),
        [r'\1.peaks.similar.control.sample',
         r'\1.peaks.similar.control.sample.locations'])
def sample_control_like_peaks(in_peaks, out_files):
    """Sample from the control IgG, with similar widths as the peaks"""
    out_sample, out_locations = out_files[:2]
    peak_lengths = array('i', (stop - start for chrom, start, stop, strand in
                                                readBedLines(open(in_peaks))))
    if len(peak_lengths) == 0:
        raise RuntimeError("Peaks file %s is empty!" % in_peaks)
    wb_genome = worldbase(cfg.get('DEFAULT', 'worldbase_genome'))
    # do the dance to map peaks back to their control raw reads
    control_bed = re.sub(r'treat', 'control', in_peaks)
    control_bed = re.sub(r'\.top[\d]+\.peaks$', '', control_bed)
    control_bed = re.sub(r'_summits\.[\d]+_around', '', control_bed)
    control_bed = re.sub(r'peaks', 'mapped_reads', control_bed)
    control_bed = re.sub(r'\.(macs(14)*|arem|glitr)', '', control_bed)
    with open(control_bed) as control_file:
        with open(out_locations, 'w') as outlocations:
            s = sampling.sample_middles(wb_genome, peak_lengths, control_file,
                                sampleSize=cfg.getint('motifs',
                                            'motif_significance_sample_size'))
            with open(out_sample, 'w') as outfile:
                for index, seq in enumerate(s):
                    # repr() gives location, str() gives sequence
                    outfile.write('>%s_%s\n%s\n' % (index, repr(seq), str(seq)))
                    outlocations.write('\t'.join([seq.id, str(seq.start),
                                                 str(seq.stop), str(index), '0',
                                '+' if seq.orientation == 1 else '-']) + '\n')

@follows(sample_genome_like_peaks)
#@collate(motif_mean_sd, regex(r'(.*)\.(.*)\.peaks(.*)\.with_mean_sd.motifs$'), 
#         add_inputs(r'\1.\2.peaks', r'\1.\2.peaks.similar.genomic.sample'), 
#         r'\1.\2.peaks\3.with_mean_sd.motifs.genomic.enrichment')

# Work in progress...
###@collate([call_peaks.all_peak_caller_functions +
###                    [get_top_peaks, '*.custom.peaks'],
###          sample_genome_like_peaks,
###          motif_mean_sd],
###    regex(), )
@active_if(False)  #whats going on with sampling???
@split(motif_mean_sd, regex(r'(.*)\.with_mean_sd\.motifs'),
       add_inputs([call_peaks.all_peak_caller_functions +
                  [get_top_peaks, '*custom.peaks'], sample_genome_like_peaks]),
       r'\1.with_mean_sd.*.vs.*.enrichment', r'\1.with_mean_sd.vs.similar.genomic.sample.zscore_%s.enrichment')
def motif_enrichment_genomic(in_files, out_pattern, out_template):
    """Determine a motif's enrichment vs. genomic samples"""
    in_motifs = in_files[0]
    in_peaks = in_files[1][0]
    in_control_samples = filter(lambda x: x.endswith('sample'), in_files[1][1:])
    
    for peak_file in in_peaks:
        # get the similar control data
        cur_control = filter(lambda x: x == (peak_file + '.similar.genomic.sample'),
                             in_control_samples)
        for c in cur_control:
            short_control = c.split(peak_file)[1][1:]
            for zscore in cfg.get('motifs', 'motif_zscores').split(','):
                outfile = out_template % (zscore)
                args = shlex.split( '%s --motif_file=%s --bg_samples=%s '
                                   '--genome=%s --output_file=%s --zscore=%s' %
                                        (peak_file, in_motifs, c,
                                         cfg.get('DEFAULT', 'worldbase_genome'),
                                         outfile, zscore))
                print args
                motif_significance.main(args)


@active_if(len(glob.glob('*control*fastq')) > 0)
@follows(sample_control_like_peaks)
@collate(motif_mean_sd, regex(r'(.*)\.(.*)\.peaks(.*)\.with_mean_sd.motifs$'), 
         add_inputs(r'\1.\2.peaks', r'\1.\2.peaks.similar.control.sample'), 
         r'\1.\2.peaks\3.with_mean_sd.motifs.control.enrichment')
def motif_enrichment_control(in_files, out_enrichment):
    """Determine a motif's enrichment vs. control data"""
    in_motifs, in_peaks, in_control_sample = in_files[0]
    for zscore in cfg.get('motifs', 'motif_zscores').split(','):
        args = shlex.split('''%s --motif_file=%s --bg_samples=%s --genome=%s
                              --output_file=%s --zscore=%s''' % (
                                        in_peaks, in_motifs, in_control_sample,
                                        cfg.get('DEFAULT', 'worldbase_genome'),
                                        out_enrichment, zscore))
        motif_significance.main(args)

@active_if(len(glob.glob('*.consensus.motifs')) > 0)
@transform([sample_control_like_peaks, sample_genome_like_peaks],
            regex(r'(.*)\.(.*)\.peaks.similar.(.*).sample$'),
            add_inputs(r'\1.\2.peaks', '*.consensus.motifs'),
            r'\1.\2.peaks.similar.\3.enrichment')
def consensus_enrichment(in_files, out_enrichment):
    """Determine a consensus motif's enrichment vs. genomic samples"""
    in_samples, in_peaks = in_files[:2]
    in_consensuses = in_files[2:]
    for in_con in in_consensuses:
        args = shlex.split('''%s --consensus_file=%s --bg_samples=%s --genome=%s
                              --output_file=%s ''' % (
                                        in_peaks, in_con, in_samples,
                                        cfg.get('DEFAULT', 'worldbase_genome'),
                                        out_enrichment))
        motif_significance.main(args)

@active_if(False)  #whats going on with sampling???
@follows(motif_mean_sd)
@split(call_peaks.all_peak_caller_functions +
                [get_top_peaks, '*.custom.peaks'] +
                [sample_control_like_peaks, sample_genome_like_peaks],
           regex(r'(.*)((\.peaks$)|(.locations$))'),
           add_inputs('*.consensus.motifs', motif_mean_sd),
           [r'\1\2*.peak_motif_presence', r'\1\2*.peak_motif_presence.png',
            r'\1*.peak_motif_locations',
            r'\1*.peak_motif_locations.bed',], r'\1', r'\2')
           #r'\1.test')
def motif_presence_sorted_peaks(in_files, out_patterns, in_prefix, in_suffix):
    """Plot the running motif presence, starting at most significant peaks"""
    in_peaks, in_motifs = in_files[0], in_files[1:]
    out_summary = in_prefix + in_suffix + '.%s.peak_motif_presence'
    out_png = in_prefix + in_suffix + '.%s.peak_motif_presence.png'
    out_locations = in_prefix + in_suffix + '.%s.peak_motif_locations'
    out_locations_bed = in_prefix + in_suffix + '.%s.peak_motif_locations.bed'
    wb_genome = worldbase(cfg.get('DEFAULT', 'worldbase_genome'))
    old_size = matplotlib.rcParams['font.size']
    matplotlib.rcParams['font.size'] = 6
    # read in the peaks file, sorting it by *score*
    print in_peaks
    print open(in_peaks).readline()
    try:
        peaks = [float(l.strip().split('\t')[4]) for l in open(in_peaks)]
        print peaks
        peaks = sorted([l.strip().split('\t') for l in open(in_peaks)],
                        key=lambda line:float(line[4]), reverse=True)
    except ValueError:
        print 'here is the error!', l.strip(), float(l.strip().split('\t')[4])
        raise
    motifs_in_peaks = dict((tuple(p), defaultdict(list)) for p in peaks)
    for m_file in in_motifs:
        cur_motifs = {}
        m_file_short = re.sub(r'((treat|fastq|fastq_illumina|min_qual|bowtie|' +
                                    r'maq|peaks|with_mean_sd|discovered|' +
                                    r'motifs_meme_out|motifs|matched_size_[0-9]|sorted|[0-9]+_around|small_sample)\.)+(motifs\.*)*',
                              '', m_file)
        #print m_file_short
        with open(m_file) as infile:
            try:
                cur_motifs.update(pickle.load(infile))
            except:
                infile.seek(0)
                for line in infile:
                    #print line,
                    name, consensus = line.strip('\n').split('\t')
                    cur_motifs.update({name:
                                    sequence_motif.makePWMFromIUPAC(consensus)})
        #print m_file, cur_motifs
        all_motif_percent = {}
        for zscore in cfg.get('motifs','motif_zscores').strip().split(','):
            for name, pwm in cur_motifs.items():
                with_motif = 0
                percent_with = []  # percent with motif at each peak
                for total, p in enumerate(peaks):
                    chrom, start, stop = p[0], int(p[1]), int(p[2])
                    region = wb_genome[chrom][start:stop]
                    # extend peaks to at least pwm length
                    while len(region) < len(pwm):
                        region = wb_genome[chrom][region.start-5:region.stop+5]
                        # catch nasty infinite loops for very short scaffolds
                        if len(region) == len(wb_genome[chrom]):
                            break
                    # check if the motif occurs in the region
                    try:
                        hits = list(pwm.find_in_region(region,
                                                       zscore=float(zscore)))
                    except Exception as e:
                        log.debug('issue with sequence', repr(region),
                                        name, e.message)
                        hits = []
                    if len(hits) > 0:
                        with_motif += 1
                        # add all peak locations to the list
                        motifs_in_peaks[tuple(p)][name].extend((
                                            h[0] + start, h[1] + start,
                                            '+' if h[2] == 1 else '-')
                                                                for h in hits)
                    percent_with.append(float(with_motif) / (total+1))
                
                #print all_motif_percent, name, percent_with
                all_motif_percent[name] = percent_with
            # having calculated for all motifs in all files,
            # plot a figure and give a summary
            with open(out_summary % ('z' + zscore), 'w') as outfile:
                outfile.writelines('%s\t%s\n' % (name, percent)
                                for name, percent in all_motif_percent.items())

            # write the peak locations along with the motif instances
            # that occur in them
            with open(out_locations % ('z' + zscore), 'w') as outfile:
                with open(out_locations_bed % ('z' + zscore), 'w') as out_bed:
                    # header is 6 columns of peak info, then motif info
                    outfile.write('\t'.join(['p_chrom', 'p_start', 'p_stop',
                                             'p_name', 'p_score', 'p_strand']))
                    for motif_name in sorted(cur_motifs):
                        outfile.write('\t%s\t#instances_%s' % (motif_name,
                                                               motif_name))
                    outfile.write('\n')
                    
                    # write one line per peak, then the motif counts and
                    # instances in the peak
                    # instances for each motif are all in one column
                    for p in peaks:
                        outfile.write('\t'.join(map(str, p)))
                        for motif_name in sorted(cur_motifs):
                            hits = motifs_in_peaks[tuple(p)][motif_name]
                            outfile.write('\t%s\t%s' % (len(hits), hits))
                            for h in hits:
                                out_bed.write('\t'.join(map(str, [p[0], h[0],
                                                        h[1], motif_name, 1000,
                                                        h[2]])) + '\n')
                        outfile.write('\n')
                    
            all_motif_percent_dict = sorted(all_motif_percent.items())
            names = [k for k, v in all_motif_percent_dict]
            datapoints = numpy.array([v for k, v in all_motif_percent_dict]).T
            
            # plot original data
            pyplot.plot(datapoints)
            pyplot.legend(names)
            pyplot.title('Motifs from\n%s\nPresence in\n%s' % (m_file_short,
                                                               in_peaks))
            pyplot.savefig(out_png % ('z'+zscore))
            pyplot.close()
            
            # plot top 10% of data
            plot_top = len(datapoints) / 10
            #print datapoints
            #print datapoints[:plot_top, :]
            # check if the slice is the right dimension
            pyplot.plot(datapoints[:plot_top, :])
            pyplot.legend(names)
            pyplot.title('Top 10%% of Motifs from\n%s\nPresence in\n%s' % (
                                                        m_file_short, in_peaks))
            pyplot.savefig(out_png % ('z' + zscore + '.top10percent'))
            pyplot.close()
        
    matplotlib.rcParams['font.size'] = old_size 




#
#sirna3_pval.txt.gene.expression.with_peaks.pygo2.treat.min_qual.bowtie.sorted.macs14.peaks.nearby.genes.reversed.ks_data
#with
#pygo2.treat.min_qual.bowtie.sorted.macs14.peaks.z4.29.peak_motif_locations

@follows(motif_presence_sorted_peaks)
@split(annotation.make_expression_ks,
    regex(r'(.*)\.with_peaks\.(.*)\.nearby\.genes\.ks_data$'),
    add_inputs(r'\2.*.peak_motif_locations'),
    r'\1.with_peaks.\2.ks_data_with.%s.peak_motif_locations.with_expression',
    r'\1.with_peaks.\2.ks_data_with.%s.peak_motif_locations.with_expression', r'\2')
def ks_with_motifs(in_files, out_pattern, out_template, peak_file):
    """join ks expression data with peak information and motif information"""
    in_ks_data, in_motifs = in_files[0], in_files[1:]
    #print 'in_data', in_ks_data
    in_motifs = filter(lambda f: 'similar' not in f and f.startswith(peak_file) and 'top' not in f.split('.z')[0], in_motifs)
    
    #print '\n\n\n\n', in_motifs
    
    ks_data = {}  # index by peak name
    for i, line in enumerate(open(in_ks_data)):
        fields = line.strip().split('\t')
        peak_id = fields[-1]
        if peak_id == 'NA':
            peak_id = 'NA_%s' % i
        ks_data[peak_id] = fields

    for i, motif_file in enumerate(in_motifs):
        #out_file = open(out_template % f.split(peak_file)[1] for f in in_motifs])
        out_file = open(out_template % i, 'w')
        out_file.write('\t'.join(['motifs_from_file', motif_file]) + '\n')
        out_file.write('\t'.join(['gene_id', 'expr_val', 'has_peak', 'peak_loc', 'peak_score', 'peak_name'] + open(motif_file).readline().strip().split('\t')[6:]) + '\n\n')
        # index motifs by peak id
        motif_data = {}
        for line in open(motif_file):
            fields = line.strip().split('\t')
            peak_id = fields[3]
            if peak_id == 'NA':
                raise RuntimeError("peak_id must be not be NA!")
            elif peak_id in motif_data:
                raise RuntimeError("peak_id must be unique! %s %s" % (peak_id, motif_data))
            else:
                motif_data[peak_id] = fields
        
        # join ks data and motif data
        for peak_id in ks_data:
            out_file.write('\t'.join(ks_data[peak_id]))
            if peak_id in motif_data:
                out_file.write('\t' + '\t'.join(motif_data[peak_id][6:]))
            out_file.write('\n')
        out_file.close()


@split([discover_meme_motifs, discover_nmica_motifs, r'*.known.motifs'],
            regex(r'(.*)\.motifs$'), r'\1.motif_*.logo.png', r'\1')
def make_seq_logo(in_pwm, eps_pattern, out_base, to_generate=1000,
                                                            out_format="PNG"):
    "convert a pwm into a sequence logo"
    import urllib, urllib2
    url = 'http://weblogo.berkeley.edu/logo.cgi'
    values = {#'sequence' : al,
          'format' : out_format,
          'logowidth' : '18',
          'logoheight' : '5',
          'logounits' : 'cm',
          'kind' : 'AUTO',
          'firstnum' : "1",
          'command' : 'Create Logo',
          'smallsamplecorrection' : "on",
          'symbolsperline' : 32,
          'res' : '96',
          'res_units' : 'ppi',
          'antialias' : 'on',
          'title' : '',
          'barbits' : '',
          'xaxis': 'on',
          'xaxis_label'  : '',
          'yaxis': 'on',
          'yaxis_label' : '',
          'showends' : 'on',
          'shrink' : '0.5',
          'fineprint' : 'on',
          'ticbits' : '1',
          'colorscheme' : 'DEFAULT',
          'color1' : 'green',
          'color2' : 'blue',
          'color3' : 'red',
          'color4' : 'black',
          'color5' : 'purple',
          'color6' : 'orange',
          }
    pwms = pickle.load(open(in_pwm))
    for name, pwm in pwms.items():
        columns = []
        for row in pwm.matrix:
            #print row
            if min(row) < 0:
                row = [2**(el + math.log(.25, 2)) for el in row]
                row = [el / sum(row) for el in row]
            #print row
            l = list(sampling.weighted_sample(zip(row, 'ACGT'), to_generate))
            random.shuffle(l)
            columns.append(l)
        instances = '\n'.join('>%s\n%s' % (i, ''.join(seq))
                                for i,seq in enumerate(zip(*columns)))
        values.update(sequence=instances)
        data = urllib.urlencode(values)
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req)
        im = response.read()
        print 'saving logo to', (out_base + '.motif_%s.logo.png' % name)
        with open(out_base + '.motif_%s.logo.png' % name, 'w') as outfile:
            outfile.write(im)

@transform(['*.consensus.motifs', motif_mean_sd], suffix(''), '.stamp_out.pdf')
def run_stamp(in_pwm, out_stamp_pdf):
    """Check for similar motifs using the STAMP webservice"""
    import urllib, urllib2
    url = 'http://www.benoslab.pitt.edu/stamp/run.php'
    to_generate = 1000
    values = {#'input':fasta_seq,
              'MAX_FILE_SIZE':'1000000',
              'mtrim':'on',
              'ICbox':'0.4',
              'motiffile':'',
              'filename':'',
              'num_hits':'10',
              'match_db':'ALL',
              'metric':'PCC',
              'align':'SWU',
              'mulalign':'IR',
              'tree':'UPGMA'
            }
    pwms = pickle.load(open(in_pwm))
    instances = '\n'.join('>%s\n%s' % (
                            name, '\n'.join(pwm.generate_sites(to_generate)))
                          for name,pwm in pwms.items())
    values.update(input=instances)
    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    im = response.read()
    print im
    job_id = re.search(r'results/(\d+)_results.html', im).groups()[0]
    pdf_url = re.sub(r'.*(?=http)', '', job_id)
    pdf_url = re.sub(r'_results\.html.*', '.pdf', pdf_url)
    req = urllib2.Request(pdf_url)
    response = urllib2.urlopen(req)
    im = response.read()
    with open(out_stamp_pdf, 'w') as outfile:
        outfile.write(im)

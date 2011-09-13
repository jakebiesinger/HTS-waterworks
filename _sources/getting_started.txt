Getting Started
===============

Expose the reads
----------------
Once you have hts-waterworks and its dependencies installed, you will need to
rename name your sequencing reads so that the hts-waterworks can see your files.

All of your fastq-formatted short reads should end in `.fastq`. If you have
ChIP-Seq reads, the treatment (immunoprecipitated reads) and control (IgG or
whole cell extract) reads should start have the same name, but treatment should
include the name `treat` while control should include the name `control`, as
follows::

    SREBP-1.treat.fastq     # immunoprecipitated reads
    SREBP-1.control.fastq   # IgG control

You can perform several comparisons using the same IgG control (or using the
same treatment data with different controls) by creating symbolic links to the
files and changing the prefix. For example::

    SREBP-1_vs_IgG.treat.fastq       # immunoprecipitated reads
    SREBP-1_vs_IgG.control.fastq     # IgG control
    SREBP-1_vs_Input.treat.fastq     # symbolic link to immunoprecipitated reads
    SREBP-1_vs_Input.control.fastq   # Whole cell extract control

The prefix used in these files will follow throughout the analysis steps and
a separate analysis will be performed on for each prefix.


Configure the analysis
----------------------

To create the initial configuration file, you will need to run::

    hts_waterworks.py -t bootstrap

This will create a new configuration file as pipeline.cfg.

There are many options available in pipeline.cfg, including selecting which one
(or many) read mappers, peak callers, or motif discovery programs will be used.
You can also update experimental details in this file (such as library tag
size).


Run the analysis
----------------

To run all of the analysis that you've specified in the config file, you can
use::

    waterworks.py -t all_complete

`all_complete` will run every job in the analysis pipeline.  There are several
other targets you may specify, such as:

.. automodule:: waterworks
    :members:

* `preprocess` which performs minimum quality thresholding, read clipping,
        and generates summary graphs of the reads
* `map_reads` which maps the reads using the mapper(s) of choice
* `peak_calling` 

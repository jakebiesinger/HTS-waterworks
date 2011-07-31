Preliminaries
=============

Installation
---------------

The simplest way to install hts-waterworks is to use the python program pip::

    pip install hts_waterworks

This will download and install hts_waterworks as well as the required python
package dependencies.

Prerequisites
-------------

hts-waterworks depends on many third-party tools to run properly and completely.
Unfortunately, installing these can be difficult and error prone and the
complexity of the task goes well beyond the scope of this document.  Some
packages are included and built automatically during installation, but others
must be installed by you.

The following python packages are required and are automatically installed
if you use pip as described above.

Required Python Packages
^^^^^^^^^^^^^^^^^^^^^^^^

* Custom Ruffus http://github.com/jakebiesinger/ruffus
* Biopython http://biopython.org/wiki/Biopython
* matplotlib http://matplotlib.sourceforge.net/
* Custom Motility http://github.com/jakebiesinger/motility/
* numpy http://numpy.scipy.org/
* pygr http://github.com/cjlee112/pygr
* pyper http://www.webarray.org/softwares/PypeR/

Required external programs
^^^^^^^^^^^^^^^^^^^^^^^^^^^
* BEDTools http://code.google.com/p/bedtools/
* Kent Admin Binaries  http://hgdownload.cse.ucsc.edu/admin/exe/
    * fetchChromSizes
    * bedClip
    * bedToBigBed
    * bedGraphToBigWig
* Fastx Toolkit http://hannonlab.cshl.edu/fastx_toolkit/

Optional external programs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Short read mappers:*

* Bowtie http://bowtie-bio.sourceforge.net/
* BWA (soon) http://bio-bwa.sourceforge.net/
* Pash (soon) http://brl.bcm.tmc.edu/pash/
* Mosaik http://bioinformatics.bc.edu/marthlab/Mosaik
* SSAHA2 http://www.sanger.ac.uk/resources/software/ssaha2/
* MAQ http://maq.sourceforge.net/

*Peak Calling:*

* MACS http://liulab.dfci.harvard.edu/MACS/
* MACS 1.4 http://liulab.dfci.harvard.edu/MACS/
* AREM http://arem.sourceforge.net/
* GLITR http://web.me.com/kaestnerlab1/GLITR/
* QuEST (soon) http://mendel.stanford.edu/SidowLab/downloads/quest/

*Motif discovery:*

* MEME http://meme.sdsc.edu/
* NestedMICA http://www.sanger.ac.uk/resources/software/nestedmica/


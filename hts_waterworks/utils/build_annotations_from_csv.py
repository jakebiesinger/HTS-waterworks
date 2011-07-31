#!/usr/bin/env python
""" Create an NLMSA and an AnnotDB based on a sqlite version of the given csv file

"""

import re
import os
import optparse
import csv
import sys
import itertools
from pygr import worldbase, cnestedlist, annotation, mapping, sqlgraph

from hts_waterworks.utils.common import (getGenome, readBedLines,
                                     bedCommentFilter, peakIter)

def main():
    """ Load the given csv file into an sqlite table, saving an
        annotationDB and an NLMSA version of the original file """

    parser = optparse.OptionParser("%prog [options] infile.csv\n"+main.__doc__)
    parser.add_option("--datapath", '-p', dest="datapath", type="string",
                      default='/home/shared/pygrdata/annotations/HUMAN/hg18',
                      help="""Sets the datafile path.  Default=%default""")
    parser.add_option("--table_name", '-t', dest="table_name", type="string",
                      help="""The resource table's name and data stem, e.g.,
                      refGene => datapath/refGene.sqlite """)
    parser.add_option("--genome", '-g', dest="genome_resource", type="string", default='hg18',
                      help="""The pygr resource for the genome, default=%default""")
    parser.add_option("--save_resource", '-r', dest="save_resource", type="string",
                      help="""Where to save the created annotationDB and NLMSA. eg, 
                      Bio.Annotation.HUMAN.hg18.MotifMap.M0001""")
    parser.add_option("--bind_attribute", '-b', dest="bind_attribute", type="string", 
                      help="""The attribute to access annotationDB from genome region, eg, 
                      'officialGenes' would be accessible via triCas3['ChLG2'][100:200].officialGenes 
                      Default is not to bind an attribute to genome""")
    parser.add_option("--slice_attrs", '-s', dest="slice_attrs", type="string",
                      default='dict(id="chromosome", start="start", stop="stop", orientation="orientation")',
                      help="""dictionary providing aliases in csv file for id, start, stop, etc. 
                      default=%default'""")
    parser.add_option("--bed_format", dest="bed_format", action='store_true',
                      help="""csv file is in BED file format, without headers.""")
    opts, args = parser.parse_args()
    if len(args) < 1: 
        parser.print_help()
        print 'Please specify at least one csv file to read'
        sys.exit(-1)
    if None in [opts.save_resource, opts.table_name]:
        parser.print_help()
        print 'Required options: save_resource, table_name'
        sys.exit(-1)
    
    fileIn = open(args[0])
    if not opts.bed_format:
        reader = csv.DictReader(fileIn, delimiter='\t')
    else:
        fileIn = itertools.ifilter(bedCommentFilter, fileIn)
        reader = csv.DictReader(fileIn, delimiter='\t', fieldnames=['chromosome', 'start', 'stop'], restkey='junkData')
    fieldnames = reader.fieldnames
    print fieldnames
    
    print '# Loading genome %s' % opts.genome_resource
    genome = getGenome(opts.genome_resource)
    
    opts.table_name = opts.table_name.replace('.','_')      # SQL interprets . as membership
    tablePath = os.path.join(opts.datapath,opts.table_name + '.sqlite')
    print '# Creating sqlite table for %s at %s' % (opts.table_name, tablePath)
    dataTable = convertBedToSQLite(reader, opts.table_name, fieldNames=fieldnames)
 
 
    
    print '# Making AnnotationDB and NLMSA...'
    annotDB = annotation.AnnotationDB(dataTable, genome, annotationType=opts.table_name+':',
                                      sliceAttrDict=eval(opts.slice_attrs))
    annotDB.__doc__ = 'AnnotationDB for %s on %s' % (opts.table_name, opts.genome_resource)
    
    msaName = os.path.join(opts.datapath, opts.table_name + '_')
    annotMap = makeNLMSA([annotDB], dataPath=msaName)

    print '# Saving results to worldbase as %s and %s...' % (opts.save_resource,
                                                             opts.save_resource+'_db')
    worldbase.add_resource(opts.save_resource, annotMap)
    worldbase.add_resource(opts.save_resource+'_db', annotDB)
    worldbase.commit()

def makeResourceFromBed(fileLines, genome, docstring='Temp Resource From BED', dataPath='memory'):
    'Generate a sqlite table, annotDB, and NLMSA from the given bed lines'
    bedLines = readBedLines(fileLines)
    bedDict = makeDictFromBed(bedLines)
    tableName = os.path.split(dataPath)[1]
    sqlDataPath = dataPath if dataPath != 'memory' else ':memory:'  # SQLite has special name for in-memory tables
    dataTable = convertDictToSQLite(bedDict, tableName, sqlDataPath)
    annotDB = annotation.AnnotationDB(dataTable, genome,
                                      sliceAttrDict=eval(defaultSliceAttrs))
    annotMap = makeNLMSA([annotDB], dataPath)
    return dataTable, annotDB, annotMap


def makeNLMSA(annotDBList, dataPath='memory'):
    if dataPath == 'memory':
        annotMap = cnestedlist.NLMSA(dataPath, 'memory', pairwiseMode=True)
    else:
        annotMap = cnestedlist.NLMSA(dataPath, 'w', pairwiseMode=True)
    annotMap.__doc__ = 'NLMSA built against '
    for annotDB in annotDBList:
        annotMap.__doc__ += ' %s, ' % (annotDB.__doc__)
        print '# Adding annotations to NLMSA from %s...' % annotDB.__doc__
        for annot in annotDB.values():
            annotMap.addAnnotation(annot)
    print '# Building annotation map...'
    if dataPath == 'memory':
        annotMap.build()
    else:
        annotMap.build(saveSeqDict=True)
    return annotMap


def convertDictToSQLite(bedDict, tableName, dataPath=':memory:'):
    liteserver = sqlgraph.SQLiteServerInfo(dataPath)
    liteserver.cursor().execute('BEGIN TRANSACTION;')
    liteserver._connection.text_factory = str
    # generate a list of all attributes defined in the file-- orientation, start, stop are all INT's
    firstLine, bedDict = peakIter(bedDict)
    fieldNames = firstLine.keys()
    attributes = [('%s INT' if field in ['start', 'stop', 'orientation'] else '%s TEXT') % field
                    for field in fieldNames]
    attributeStr = 'rowid INT PRIMARY KEY, ' + ', '.join(attributes)
    dataTable = sqlgraph.SQLTable(tableName, serverInfo=liteserver, writeable=True,
                                    dropIfExists=True,
                                    createTable='CREATE TABLE %s (%s);' % 
                                                (tableName, attributeStr))
    for index, lineData in enumerate(bedDict):
        # add entry to table
        if 'hap' in lineData['chromosome']:
            continue
        lineData['rowid'] = index
        outputAnnot = dataTable.new(**lineData)
    liteserver._connection.commit()
    return dataTable

def makeDictFromBed(bedLines):
    return ( dict(zip(['chromosome', 'start', 'stop', 'strand'], line)) for line in bedLines )

defaultSliceAttrs = 'dict(id="chromosome", start="start", stop="stop", orientation="orientation")'

if __name__ == '__main__':
    main()

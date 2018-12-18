#slin 201707  added lib path and others.

import os
import sys

SCRIPTS_DIR = os.path.dirname(__file__)
PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

import numpy as np
from collections import defaultdict

import ArgLib
import DataFile
import Validator

import glob


def getOutputName(outputBase, prog, ext):
    return outputBase + '_' + prog + ext


def countTargetsAndDecoys(proc_tag_graph):
    decoys = defaultdict(int)

    for item in proc_tag_graph.values():
        decoys[item['Decoy Status']] += 1

    return decoys

def filterToFDR(proc_tag_graph, desired_fdr=0.01, write_roc=False):
    scores, crit = [], []

    for item in proc_tag_graph.values():
        scores += [float(item['Alignment Score'])]
        crit += [0] if item['Decoy Status'] == 'Yes' else [1]

    sens, prec, score_cutoffs = Validator.getSensitivityAndPrecisionArrs(np.array(scores), np.array(crit), numPoints=200)
    
    if write_roc:
        Validator.writeSensitivityAndPrecisionToFile(np.array(scores), sens, prec, getOutputName(outBase, 'ROC', '.tdv'))
    
    # This logic is not very pythonic
    # Maximize sensitivity if under the desired FDR, else just try to minimize the FDR
    sens_under, under_ind = 0, -1
    fdr_over, over_ind = 1, -1
    for i in range(len(prec)):
        if sens[i] > 0:
            fdr = 1 - prec[i]
            
            if fdr > desired_fdr and fdr - desired_fdr < fdr_over:
                fdr_over = fdr - desired_fdr
                over_ind = i
            elif fdr <= desired_fdr and sens_under < sens[i]:
                sens_under = sens[i]
                under_ind = i

    if under_ind < 0:
        fdr = 1 - prec[over_ind]
        score_cutoff = score_cutoffs[over_ind]
        print 'Could not filter to desired FDR %f %%. Filtering to FDR of %f %% instead...'%(100*desired_fdr, 100*fdr)
    else:
        fdr = 1 - prec[under_ind]
        score_cutoff = score_cutoffs[under_ind]
        print 'Filtering to FDR of %f %%...'%(100*fdr,)

    # Filtering...
    for scanF in proc_tag_graph.keys():
        if float(proc_tag_graph[scanF]['Alignment Score']) < score_cutoff:
            del proc_tag_graph[scanF]
            

def writeOutput(outfile_name, proc_tag_graph):
    outFile = open(outfile_name, 'w')
    cols = ['ScanF', 'Alignment Score', 'Context', 'Modifications', 'Proteins', 'Matching Tag Length', 'PEAKS ALC (%)', 'PEAKS Peptide', 'Decoy?']
    
    outFile.write('\t'.join(cols) + '\n')
    for item in proc_tag_graph.values():
        scanData = {'ScanF': item['ScanF'], 'Alignment Score': item['Alignment Score'], 'Matching Tag Length': item['Matching Tag Length'], 'PEAKS ALC (%)': item['De Novo Score'], 'PEAKS Peptide': item['De Novo Peptide'], 'Decoy?': item['Decoy Status']}
        for context in item['Context']:
            scanData['Context'] = context[0]
            scanData['Modifications'] = context[1]
            scanData['Proteins'] = context[2]
            
            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
            
    outFile.close()
                                                                        

if __name__ == '__main__':
    print 'dtadir argument points to parent directory of de_novo and mzML files'
    print 'output is a string of values to append to the de novo filename in the output'
    options = ArgLib.parse(['init', 'dtadir', 'ppmstd', 'modtolerance', 'unimoddict', 'modmaxcounts', 'maxcounts', 'fmindex', 'denovo', 'model', 'config', 'output'], [{'opts': ('-F', '--fraction'), 'attrs': {'type': 'int', 'dest': 'fraction', 'help': 'Fraction to run TAG_GRAPH on'}}, {'opts': ('-x', '--splittaxon'), 'attrs': {'dest': 'splittaxon', 'action': 'store_true', 'default': False, 'help': 'Flag. For searches of metaproteomic databases, split identical context entries by taxon for accurate consideration via EM.'}}])

    fileFound = False
    
    outDir = os.path.join(options.dtadir, 'taggraph')
    peaksDir = os.path.join(options.dtadir, 'de_novo')
    dataDir = os.path.join(options.dtadir, 'data')
    localDtaDir = ''
    try:
        # This should throw an exception if the directory for dta files does not yet exist (it will then be created)
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        dtaDir = glob.glob(dataDir + os.path.sep + '*f%02d'%options.fraction)[0] + os.path.sep
        '''
        localDtaDir = glob.glob(dataDir + '/' + '*f%02d'%options.fraction)[0] + '/'
        print localDtaDir+"\n"
    except IndexError:
        # Only unpack dta files if not already unpacked

        ''' Check whether we have mzML or mzXML files, and extract appropriately '''

        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        mzMLFile = os.path.abspath(glob.glob(dataDir + os.path.sep + '*f%02d.mzML'%options.fraction)[0])
        '''
        mzMLFiles  = glob.glob(dataDir + '/' + '*f%02d.mzML'%options.fraction)
        mzXMLFiles = glob.glob(dataDir + '/' + '*f%02d.mzXML'%options.fraction)

        if len(mzMLFiles) > 0:
            ## Found mzML file for this fraction
            fileFound = True
            mzMLFile = os.path.abspath(mzMLFiles[0])
            
            ## The directory for DTA files will be the same as the mzML file without the .mzML extension
            mzml_name_base = mzMLFile[:-5]
            print 'mzMLFile: "%s"' % (mzMLFile)
            print 'mzml_name_base: "%s"' % (mzml_name_base)
        
            # Unpack DTAs
            DataFile.executeProcess(SCRIPTS_DIR, 'mzml2dta.py', ['-o', mzml_name_base, mzMLFile])
            ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
            dtaDir = glob.glob(dataDir + os.path.sep + '*f%02d'%options.fraction)[0] + os.path.sep
            '''
            localDtaDir = glob.glob(dataDir + '/' + '*f%02d'%options.fraction)[0] + '/'
            print 'Found mzML, setting dtaDir to %s' % (localDtaDir)
        elif len(mzXMLFiles) > 0:
            ## Found mzXML file for this fraction
            fileFound = True
            mzXMLFile = os.path.abspath(mzXMLFiles[0])
            
            '''
            ## The directory for DTA files will be the same as the mzML file without the .mzML extension
            mzml_name_base = mzMLFile[:-5]
            ## print 'mzml_name_base: "%s"' % (mzml_name_base)
            '''
            
            # Unpack DTAs
            mzXML2SearchDir =  os.path.abspath( os.path.join(os.path.join( PAR_DIR, 'lib'), 'MzXML2Search') )
            DataFile.executeProcess(mzXML2SearchDir, 'MzXML2Search', ['-dta', '-P1', '-T10000', '-B200', mzXMLFile], interpreter = False)
            ##DataFile.executeProcess(SCRIPTS_DIR, 'mzml2dta.py', ['-o', mzml_name_base, mzMLFile])
            ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
            dtaDir = glob.glob(dataDir + os.path.sep + '*f%02d'%options.fraction)[0] + os.path.sep
            '''
            localDtaDir = glob.glob(dataDir + '/' + '*f%02d'%options.fraction)[0] + '/'
            print 'Found mzXML, setting dtaDir to %s' % (localDtaDir)
        else:
            print '*** FAILURE *** No mzML or mzXML data file found for fraction %02d' % options.fraction


    if fileFound == True:
        # Run TAG-GRAPH
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        deNovoFile = os.path.abspath(glob.glob(peaksDir + os.path.sep + '*F%i.tdv'%options.fraction)[0])
        '''
        deNovoFile = os.path.abspath(glob.glob(peaksDir + '/' + '*F%i.tdv'%options.fraction)[0])
        args = DataFile.getArgs(options, ['init', 'ppmstd', 'modtolerance', 'unimoddict', 'modmaxcounts', 'maxcounts', 'fmindex', 'model', 'config'])
        args += ['--splittaxon'] if options.splittaxon else []
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        args += ['--dtadir', dtaDir, '--denovo', deNovoFile, '--output', os.path.abspath(outDir + os.path.sep + os.path.splitext(os.path.basename(deNovoFile))[0] + '_' + options.output + '.tdv')]
        '''
        args += ['--dtadir', localDtaDir, '--denovo', deNovoFile, '--output', os.path.abspath(outDir + '/' + os.path.splitext(os.path.basename(deNovoFile))[0] + '_' + options.output + '.tdv')]
        DataFile.executeProcess(LIB_DIR, 'TAG_GRAPH_PROB_SCORE.py', args)

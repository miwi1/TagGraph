'''
Created on Sep 15, 2011

@author: Arun
'''
#slin 201707  added lib path and others.

import os
import sys

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

import ArgLib
import Constants
import DataFile

import numpy as np
import pickle


def getPeptideData(progName, progDict, scanF):
    try:
        seq = processedInfo[progName][scanF]['Peptide']
        if 'Ambiguous Edges' in processedInfo[progName][scanF]:
            return (seq, processedInfo[progName][scanF]['Ambiguous Edges'])
        else:
            return (seq, None)
    except KeyError:
        return False

def parseScans(fDict, prog, seqMap, dbDict, delimiter=',', srchID = None, seqDelimLen=2):
    processedInfo = {}
    for csvfile in fDict.keys():
        cols, data = DataFile.getScanInfo(csvfile, dbDict[prog]['fields'] + (['SrchID'] if srchID != None else []), delimiter=delimiter)
        processedInfo[fDict[csvfile]] = DataFile.preprocessDatabaseScanInfo(data, seqMap[fDict[csvfile]], dbDict[prog]['fieldmap'], seqDelimLen=seqDelimLen)
    
    return processedInfo

def updateProgDict(fDict, progDict, prog):
    for name in fDict.values():
        progDict[name] = prog

def getAllScanF(processedInfo):
    scanFs = np.array([], dtype=np.dtype('int'))
    for progName in processedInfo.keys():
        scanFs = np.append(scanFs, np.array(processedInfo[progName].keys(), dtype=np.dtype('int')))

    return np.unique(scanFs)


if __name__ == '__main__':
    options = ArgLib.parse(['init', 'lads', 'sequest', 'mascot', 'pepnovo', 'output', 'database', 'symbolmap', 'pnovo', 'peaks', 'combined'])
    
    paramsDict = ArgLib.parseInitFile(options.init, options)
    progDict = ArgLib.getProgDict(DataFile.searchprogs, options)
    
    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)

    seqMap = DataFile.generateSeqMap(progDict, symbolMap, paramsDict)

    dbDict = DataFile.getDBInfo(options.database)
    processedInfo = {}  
    if options.lads:
        LADSdict = eval(options.lads)
        for tdvfile in LADSdict.keys():
            LADSScanInfo = DataFile.getScanInfo(tdvfile, dbDict['LADS']['fields'], delimiter='\t')
            processedInfo[LADSdict[tdvfile]] = DataFile.preprocessLADSScanInfo(LADSScanInfo, seqMap[LADSdict[tdvfile]], paramsDict['Pair Configurations'], dbDict['LADS']['fieldmap'])

    if options.pepnovo:
        pepNovoDict = eval(options.pepnovo)
        processedInfo.update(parseScans(pepNovoDict, 'PepNovo', seqMap, dbDict, delimiter='\t', seqDelimLen=0))

    if options.mascot:
        MASCOTdict = eval(options.mascot)
        processedInfo.update(parseScans(MASCOTdict, 'MASCOT', seqMap, dbDict))
        
    if options.sequest:
        SEQUESTdict = eval(options.sequest)
        processedInfo.update(parseScans(SEQUESTdict, 'SEQUEST', seqMap, dbDict))
        
    if options.pnovo:
        pNovoDict = eval(options.pnovo)
        processedInfo.update(parseScans(pNovoDict, 'pNovo', seqMap, dbDict, delimiter='\t', seqDelimLen=0))

    if options.peaks:
        peaksDict = eval(options.peaks)
        processedInfo.update(parseScans(peaksDict, 'PEAKS', seqMap, dbDict, delimiter='\t', seqDelimLen=0))
        
    if options.combined:
        combinedDict = eval(options.combined)
        processedInfo.update(parseScans(combinedDict, 'Combined', seqMap, dbDict, delimiter='\t', seqDelimLen=0))
    
    cols = ['ScanF']
    progNames = processedInfo.keys()
    for progName in progNames:
        cols.extend([progName + ' ' + val for val in dbDict[progDict[progName]]['cols']])
    
    consensi = {}
    for i in range(len(progNames)):
        consensi[progNames[i]] = []
        for j in range(i + 1, len(progNames)):
            cols.extend([progNames[i] + ' ' + progNames[j] + ' ' + item for item in ('Precision', 'Accuracy')])
            consensi[progNames[i]] += [progNames[j]]

    outFile = open(options.output, 'w')
    outFile.write('\t'.join([col for col in cols]) + '\n')
    
    for i in getAllScanF(processedInfo):
        scanData = {}
        scanData['ScanF'] = i
        for progName in progNames:
            for col in dbDict[progDict[progName]]['cols']:
                try:
                    scanData[progName + ' ' + col] = processedInfo[progName][i][col]
                except KeyError:
                    scanData[progName + ' ' + col] = None
            progPept = getPeptideData(progName, progDict, i)

            for sprogName in consensi[progName]:
                sProgPept = getPeptideData(sprogName, progDict, i)
                if progPept and sProgPept:
                    comp = Constants.comparePeptideResults(progPept[0], sProgPept[0], ambigEdges1=progPept[1], ambigEdges2=sProgPept[1])
                else:
                    comp = (None, None, None)
                for j, item in enumerate(['Precision', 'Accuracy', 'Consensus']):
                    scanData[progName + ' ' + sprogName + ' ' + item] = comp[j]

        outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
        print '\nScan Number %i \n' % i 
        print '\t'.join(['Program', 'Peptide', 'Reference', 'Score', 'Obs M+H'])
        
        for progName in progNames:
            peptide = str(scanData[progName + ' ' + dbDict['infoMap'][progDict[progName]]['Peptide']])
            score = dbDict['infoMap'][progDict[progName]]['Score'] + ': ' + str(scanData[progName + ' ' + dbDict['infoMap'][progDict[progName]]['Score']])
            try:
                MH = str(scanData[progName + ' ' + dbDict['infoMap'][progDict[progName]]['Obs M+H']])
            except KeyError:
                MH = str(None)
            try:
                reference = str(scanData[progName + ' ' + dbDict['infoMap'][progDict[progName]]['Reference']])
            except KeyError:
                reference = str(None)
            print '\t'.join([progName, peptide, reference, score, MH])
        
        compHead = []
        compLineA = []
        compLineB = []
        for progName in consensi.keys():
            for sprogName in consensi[progName]:
                compHead += [progName + ' ' + sprogName +  ' P: %s A: %s' % (str(scanData[progName + ' ' + sprogName + ' Precision']), str(scanData[progName + ' ' + sprogName + ' Accuracy']))]
                try:
                    compLineA += [str(scanData[progName + ' ' + sprogName + ' Consensus'][0])]
                    compLineB += [str(scanData[progName + ' ' + sprogName + ' Consensus'][1])]
                except TypeError:
                    compLineA += ['N/A']
                    compLineB += ['N/A']
        print '\n' + '\t'.join(compHead) + '\n' + '\t'.join(compLineA) + '\n' + '\t'.join(compLineB) + '\n'
            
            
        
        

                   

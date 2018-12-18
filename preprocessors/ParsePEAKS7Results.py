import sys
import os

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

import Constants
import ArgLib
import DataFile

import csv
import pickle
import glob


def getMatchingModSymbol(aa, index, pept_length, seqMap, mod):
    pos_locs = []
    if index == 0:
        pos_locs += ['N-term']
    elif index + 1 == pept_length:
        pos_locs += ['C-term']

    pos_locs += [aa]
    for loc in pos_locs:
        try:
            return seqMap['PEAKS']['Mods'][(mod, loc)]
        except KeyError:
            pass
    
    raise KeyError('%s with possible localizations %s not found in seqMap'%(mod, pos_locs))

def getPEAKSProcSeq(seq, seqMap, pept_length, LCs, cutOff = 0, ambigAA='@'):
    #print 'cutoff', cutOff
    procSeq = []
    mod = ''
    i = -1
    curr_aa = None
    seqList = iter(seq)
    #print LCs
    try:
        while True:
            aa = seqList.next()
            if aa in seqMap['PEAKS']['AAs']:
                procSeq += [seqMap['PEAKS']['AAs'][aa]]
                cur_aa = aa
                i += 1
            else:
                modLetter = aa
                if modLetter == '(':
                    mod = ''
                    while modLetter != ')':
                        mod += modLetter
                        modLetter = seqList.next()
                    mod += modLetter
                    procSeq[-1] = procSeq[-1] + getMatchingModSymbol( cur_aa, i, pept_length, seqMap, ''.join(mod) )
                    mod = ''
    except StopIteration:
        if mod != '':
            procSeq[-1] = procSeq[-1] + getMatchingModSymbol( cur_aa, i, pept_length, seqMap, ''.join(mod) )

    #print procSeq
    if cutOff == 0:
        return ''.join(procSeq), []
    else:
        ambig_edges = []
        for i, aa in enumerate(procSeq):
            if LCs[i] < cutOff:
                procSeq[i] = ambigAA
                ambig_edges += [(0, Constants.aminoacids[aa][2])]
        return ''.join(procSeq), ambig_edges
                    
if __name__ == '__main__':
    print 'In this program, the PEAKS argument is just the location of the PEAKS output to parse. Number argument indicates ALC cutoff to form ambig edges (set to 0 to not form any amibiguous edges'
    options = ArgLib.parse(['init', 'output', 'symbolmap', 'peaks', 'cutoff'])

    AMBIG_AA = '@'
    
    paramsDict = ArgLib.parseInitFile(options.init, options)
    with open(options.symbolmap, 'r') as fin:
        symbolMap = pickle.load(fin)
    seqMap = DataFile.generateSeqMap({'PEAKS': 'PEAKS'}, symbolMap, paramsDict)
    #print seqMap
    
    scanInfo = DataFile.getScanInfo(options.peaks, delimiter=',')[1]
    if 'Peptide' in scanInfo[0]:
        seq_col = 'Peptide'
    else:
        seq_col = 'Sequence'
    
    outFile = open(options.output, 'w')
    
    #print 'cutoff_arge', options.cutoff
    cols = ['ScanF', 'Charge', 'RT', 'Obs M+H', 'Peptide', 'ALC (%)', 'LC']
    alc_cutoff = options.cutoff if options.cutoff else 0
    if alc_cutoff > 0:
        cols.insert(-2, 'Ambig Edges')
    #print 'alccutoff', alc_cutoff
    outFile.write('\t'.join([col for col in cols]) + '\n')

    scanData = {}
    for scan in scanInfo:
        
        if 'Fraction' in scan:
            scanData['ScanF'] = 'F' + scan['Fraction'] + ':' + scan['Scan']
        else:
            scanData['ScanF'] = scan['Scan']

            
        scanData['Charge'] = scan['z']
        scanData['RT'] = scan['RT']
        scanData['Obs M+H'] = float(scan['m/z']) * int(scan['z']) - ((int(scan['z']) - 1) * Constants.mods['H+'])

        pept_length = int(scan['Tag length']) if ('Tag length' in scan) else int(scan['length'])
        peptide, ambig_edges = getPEAKSProcSeq(scan[seq_col], seqMap, pept_length, [int(LC) for LC in scan['local confidence (%)'].split(' ')], alc_cutoff)

        # Don't write the peptide if it is just ambig amino acids (will cause seg fault in tag-graph)
        if all([aa == AMBIG_AA for aa in peptide]):
            continue
        
        scanData['Peptide'] = peptide
        if alc_cutoff > 0:
            scanData['Ambig Edges'] = ambig_edges
        scanData['ALC (%)'] = scan['ALC (%)']
        scanData['LC'] = scan['local confidence (%)']
        outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

    outFile.close()

import sys
import os
import networkx as nx


import numpy as np
import time
import pickle
import glob
from collections import defaultdict, deque

import anydbm

CUR_DIR = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(CUR_DIR,'fmindex'))

import fmindex as fm
import ArgLib
import DataFile
import SixFrameDNATranslator as SFDT
import Constants



###################
# KNOWN BUGS
# Adding arginine mass mod (+156) instead of matching an isobaric sub to arginine (due to scoring matching AA +1), penalize mods more?
# Allowing insertion and deletion of AAs instead of matching isobaric sub (due to scoring matching AA + 1), penalize mods more?
####################

EXACT_MATCH_SCORE = 100

def grabGreedyTags(diags):
    acceptedTagRegions = []
    for diag in diags:
        diagX, diagY = zip(*diag)
        
        acceptTag = True
        for tag in acceptedTagRegions:
            tagX, tagY = zip(*tag)
            print (not set(tagX) & set(diagX)), (not set(tagY) & set(diagY))
            if set(tagX) & set(diagX) or set(tagY) & set(diagY):
                acceptTag = False

        if acceptTag:
            acceptedTagRegions += [diag]

    return acceptedTagRegions
    
                                                 

def assembleSequenceTags(indicesDict):
    overlapRegions = []
    tagInds = sorted(indicesDict.keys())


    # Get overlapping sequence tag match regions
    diagHash = defaultdict(deque)
    for tagInd in tagInds:
        for match in indicesDict[tagInd]:
            diagHash[tagInd - match] += [(tagInd, match)]

    #print 'Diag Hash', diagHash
    #Break Apart non-contiguous diagonals
    diags = []
    for diag in diagHash.values():
        tempDiag = [diag.popleft()]
        while diag:
            pair = diag.popleft()
            if tempDiag[-1][1] + 1 == pair[1]:
                tempDiag += [pair]
            else:
                diags += [tempDiag]
                tempDiag = [pair]

        diags += [tempDiag]
            

    return diags
            
                
def generateSequenceTags(deNovoPept, dbPept, tagLength = 3, ambigAA = 'X'):
    indicesDict = {}

    for i in range(len(deNovoPept) - tagLength + 1):
        tag = deNovoPept[i:i+tagLength]
        if ambigAA not in tag:
            matches = list(findAll(tag, dbPept))
            if matches:
                indicesDict[i] = matches

    #print 'Indices Dict', indicesDict
    # Assign tags to peptide positions in a greedy fashion
    acceptedTagRegions = assembleSequenceTags(indicesDict)

    #print 'Accepted Tag Regions', acceptedTagRegions
    acceptedTags = []
    for tagRegion in sorted(acceptedTagRegions, key = lambda t: t[0][0]):
        acceptedTags += [((tagRegion[0][0], tagRegion[-1][0] + tagLength), (tagRegion[0][1], tagRegion[-1][1] + tagLength))]

    return acceptedTags

def greedyTagsToAlignment(deNovoPept, dbPept, acceptedTags):
    alignment = [deNovoPept[:acceptedTags[0][0][0]], dbPept[:acceptedTags[0][0][1]]]

    for i, tag in enumerate(acceptedTags):
        alignment[0] += '*%s*' % (deNovoPept[tag[0][0]:tag[1][0]],)
        alignment[1] += '*%s*' % (dbPept[tag[0][1]:tag[1][1]],)

        if i < len(acceptedTags) - 1:
            alignment[0] += deNovoPept[tag[1][0]:acceptedTags[i+1][0][0]]
            alignment[1] += dbPept[tag[1][1]:acceptedTags[i+1][0][1]]

    alignment[0] += deNovoPept[tag[1][0]:]
    alignment[1] += dbPept[tag[1][1]:]

    return alignment


def generateStartAndEndTags(deNovoPept, dbPept):
    startTags = []
    endTags = []

    for i in range(len(dbPept) + 1):
        startTags += (((0,0), (i,i)),)
        endTags += (((len(deNovoPept), len(deNovoPept)), (i,i)),)

    return startTags, endTags

def getSequenceTagGraph(startTags, endTags, sequenceTags):
    tagGraph = nx.DiGraph()

    #print 'Sequence Tags', sequenceTags
    for tag in startTags:
        tagGraph.add_node(tag, position="start")
    for tag in endTags:
        tagGraph.add_node(tag, position="end")
    for tag in sequenceTags:
        tagGraph.add_node(tag, position="internal")

    for startTag in startTags:
        for internalTag in sequenceTags:
            if startTag[0][1] <= internalTag[0][0] and startTag[1][1] <= internalTag[1][0]:
                tagGraph.add_edge(startTag, internalTag)

    for endTag in endTags:
        for internalTag in sequenceTags:
            if internalTag[0][1] <= endTag[0][0] and internalTag[1][1] <= endTag[1][0]:
                tagGraph.add_edge(internalTag, endTag)

    for i in range(len(sequenceTags)):
        for j in range(i):
            tag1 = sequenceTags[i]
            tag2 = sequenceTags[j]

            if tag1[0][1] <= tag2[0][0] and tag1[1][1] <= tag2[1][0]:
                tagGraph.add_edge(tag1, tag2)
            elif tag2[0][1] <= tag1[0][0] and tag2[1][1] <= tag1[1][0]:
                tagGraph.add_edge(tag2, tag1)

    return tagGraph

def alignDeNovoToDBSequence(deNovoPRMLadder, deNovoPept, dbPept, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges = None, tagLength=2, isobaricPenalty=-0.5, defModPenalty=-1, inDelPenalty=-2, undefModPenalty=-3, defaultScore=0, cutOffDiff=2):
    #print deNovoPRMLadder

    #print 'De Novo', deNovoPept
    #print 'DB', dbPept
    
    dbPRMLadder = Constants.getPRMLadder(dbPept, addEnds=True)

    startTags, endTags = generateStartAndEndTags(deNovoPept, dbPept)
    sequenceTags = generateSequenceTags(deNovoPept, dbPept, tagLength=tagLength)

    tagGraph = getSequenceTagGraph(startTags, endTags, sequenceTags)

    #print sorted(tagGraph.nodes(data=True))
    #print sorted(tagGraph.edges(data=True))
    topo_sort = nx.topological_sort(tagGraph)
    maxLScore, maxEndScoringTag = scoreAndInsertLScores(tagGraph, topo_sort, deNovoPRMLadder, deNovoPept, dbPRMLadder, dbPept, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges, tagLength, isobaricPenalty, defModPenalty, inDelPenalty, undefModPenalty, defaultScore)
    maxRScore, maxStartScoringTag = insertRScores(tagGraph, topo_sort, defaultScore)

    max_scoring_paths = getBestAlignments(tagGraph, maxLScore-cutOffDiff, startTags)
    pept_mod_data = [((path[0][1][0], path[-1][1][0]), [tagGraph.edge[ path[i] ][ path[i+1] ]['mods'] for i in range(len(path) - 1) if len(tagGraph.edge[ path[i] ][ path[i+1] ]['mods']) > 0], score) for score, path in max_scoring_paths]

    return pept_mod_data

def scoreAndInsertLScores(tagGraph, topo_sort, deNovoPRMLadder, deNovoPept, dbPRMLadder, dbPept, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges = None, tagLength=2, isobaricPenalty=-0.5, defModPenalty=-1, inDelPenalty=-2, undefModPenalty=-3, defaultScore=0):

    maxScore, maxScoringTag = None, None
    
    for tag in topo_sort:
        nodeScore = tag[0][1] - tag[0][0]
        #print 'Tag', tag
        for prevTag in tagGraph.predecessors(tag):
            nModSymbol = None
            # Define terminus of peptide for modification annotation
            if tagGraph.node[prevTag]['position'] == 'start':
                term = 'N-term'
            elif tagGraph.node[tag]['position'] == 'end':
                term = 'C-term'
            else:
                term = None

            
            refMass = dbPRMLadder[tag[1][0]] - dbPRMLadder[prevTag[1][1]]
            deNovoMass = deNovoPRMLadder[tag[0][0]] - deNovoPRMLadder[prevTag[0][1]]
            refSubSeq = dbPept[prevTag[1][1]:tag[1][0]]
            deNovoSubSeq = deNovoPept[prevTag[0][1]:tag[0][0]]

            mods = resolveInterval(refMass, deNovoMass, refSubSeq, deNovoSubSeq, hashedUnimodDict, unimodDict, paramsDict, term=term, nModSymbol=nModSymbol)

            if not mods:
                modPenalty = 0
            elif 'Isobaric Substitution' == mods[0][0][0]:
                modPenalty = isobaricPenalty
            elif 'Insertion' == mods[0][0][0] or 'Deletion' == mods[0][0][0]:
                modPenalty = inDelPenalty
            elif 'Undefined Mass Shift' == mods[0][0][0]:
                modPenalty = undefModPenalty
            else:
                modPenalty = defModPenalty

            edge_score = nodeScore + modPenalty
            tagGraph.edge[prevTag][tag]['edge_score'] = edge_score
            tagGraph.edge[prevTag][tag]['mods'] = mods

            #print prevTag, tag, deNovoSubSeq, refSubSeq, mods
            
            if 'lscore' not in tagGraph.node[prevTag]:
                tagGraph.node[prevTag]['lscore'] = defaultScore

            try:
                tagGraph.node[tag]['lscore'] = max(tagGraph.node[tag]['lscore'], tagGraph.node[prevTag]['lscore'] + edge_score)
            except KeyError:
                tagGraph.node[tag]['lscore'] = tagGraph.node[prevTag]['lscore'] + edge_score

        try:
            if tagGraph.node[tag]['position'] == 'end' and tagGraph.node[tag]['lscore'] > maxScore:
                maxScore = tagGraph.node[tag]['lscore']
                maxScoringTag = tag
        except KeyError:
            pass


    return maxScore, maxScoringTag

def insertRScores(tagGraph, topo_sort, defaultScore = 0):

    maxScore, maxScoringTag = None, None

    for tag in topo_sort[::-1]:
        if 'rscore' not in tagGraph.node[tag]:
            tagGraph.node[tag]['rscore'] = defaultScore

        for prev_tag in tagGraph.predecessors(tag):
            edge_score = tagGraph.edge[prev_tag][tag]['edge_score']
            try:
                tagGraph.node[prev_tag]['rscore'] = max(tagGraph.node[prev_tag]['rscore'], tagGraph.node[tag]['rscore'] + edge_score)
            except KeyError:
                tagGraph.node[prev_tag]['rscore'] = tagGraph.node[tag]['rscore'] + edge_score

        if tagGraph.node[tag]['position'] == 'start' and tagGraph.node[tag]['rscore'] > maxScore:
            maxScore = tagGraph.node[tag]['rscore']
            maxScoringTag = tag

    return maxScore, maxScoringTag

def getBestAlignments(tagGraph, cutOff, startTags):

    paths_passing_cutoff = []
    paths = deque([(0, [tag]) for tag in startTags])

    while paths:
        path_score, current_path = paths.pop()
        current_node = current_path[-1]
        
        for s in tagGraph.successors(current_node):
            s_node = tagGraph.node[s]
            s_score =  path_score + tagGraph.edge[current_node][s]['edge_score']
            
            max_attainable_score = s_score + s_node['rscore']
            if max_attainable_score >= cutOff:
                if s_node['position'] == 'end':
                    paths_passing_cutoff += [(s_score, current_path + [s])]
                else:
                    paths.extend([(s_score, current_path + [s])])

    return paths_passing_cutoff
                
    

def calculateMassDeltaFrequencies(precMassArr, epsilon=0.01):
    precMassMatrix = np.reshape(np.tile(precMassArr, precMassArr.size), (precMassArr.size, precMassArr.size))
    deltaMatrix = precMassMatrix - precMassMatrix.transpose()

    freqHash = defaultdict(int)
    for entry in deltaMatrix.flat:
        hMass = hashMass(entry, epsilon)
        for m in range(hMass-1, hMass+2):
            freqHash[m] += 1

    return freqHash

def getScoringModel(freqHash):
    return 0

    
def getBestAlignment(tagGraph, dbPept, maxScore, maxScoringTag):
    modListReversed = []

    endInd = maxScoringTag[1][0]
    currentNode = maxScoringTag
    currentScore = maxScore
    while True:
        for node in tagGraph.predecessors(currentNode):
            if currentScore - tagGraph.edge[node][currentNode]['edgeScore'] == tagGraph.node[node]['score']:
                modListReversed += [tagGraph.edge[node][currentNode]['mods']] if len(tagGraph.edge[node][currentNode]['mods']) > 0 else []
                if tagGraph.node[node]['position'] == 'start':
                    return (node[1][0], endInd), modListReversed[::-1], maxScore

                currentNode = node
                currentScore = tagGraph.node[node]['score']
                break
            


def findAll(sub, string):
    index = -1
    try:
        while True:
            index = string.index(sub, index + 1)
            yield index
    except ValueError:
        pass

def updateSeqCount(seqCountArr, seqModArr, proteinSeq, peptide, modList = None, countUp = 1, maxModLength=5):
    #print proteinSeq, peptide
    startPos = proteinSeq.index(peptide)
    #print startPos, len(proteinSeq)
    for i in range(startPos, startPos + len(peptide)):
        seqCountArr[i] += countUp

    
    if modList != None:
        for mod in modList:
            if 'Isobaric' not in mod[0][0] and len(mod[2]) < maxModLength:
                seqModArr[startPos + peptide.index(mod[2])] += 1
    

def getConnectedDisagreementRegions(disagreeArr, trueVal = 0):
    intervals = []

    startInd = -1
    endInd = -1
    for i, elem in enumerate(disagreeArr):
        if elem == trueVal:
            if startInd == -1:
                startInd = i
            endInd = i+1
        else:
            if startInd != -1:
                intervals += [(startInd, endInd)]
            startInd = endInd = -1

    if startInd != -1:
        intervals += [(startInd, endInd)]
        
    return intervals
                

def hashUnimodDict(unimodDict, epStep=0.0025, maxEp=0.1):
    hUnimodDict = {}

    maxIntEp = int(np.ceil(maxEp/epStep))
    minIntEp = int(np.floor(-maxEp/epStep))

    for mod in unimodDict:
        massKey = np.round(unimodDict[mod]['mass']/epStep)

        for intEp in range(minIntEp, maxIntEp+1):
            try:
                hUnimodDict[intEp + massKey][mod] = epStep*intEp
            except KeyError:
                hUnimodDict[intEp + massKey] = {mod: epStep*intEp}

    return hUnimodDict

def hashMass(modMass, epStep=0.0025):
    return int(np.round(modMass/epStep))

def getAlignedIndsMap(alignment):
    alignedIndsMap = {'De Novo': {}, 'Ref': {}}
    numAADeNovo = 0
    numAARef = 0
    
    for i in range(len(alignment[0])):
        alignedIndsMap['De Novo'][i] = numAADeNovo
        alignedIndsMap['Ref'][i] = numAARef
        if alignment[0][i] != '-':
            numAADeNovo += 1
        if alignment[1][i] != '-':
            numAARef += 1

    alignedIndsMap['De Novo'][len(alignment[0])] = numAADeNovo
    alignedIndsMap['Ref'][len(alignment[0])] = numAARef
    return alignedIndsMap

def isValidModLocation(unmodSeq, term, location):
    if location[0] in unmodSeq:
        if location[1] == 'Anywhere':
            return True
        elif location[1] == term:
            if term == 'N-term' and unmodSeq[0] == location[0]:
                return True
            elif term == 'C-term' and unmodSeq[-1] == location[0]:
                return True
            else:
                return False
        else:
            return False
    elif location[0] == term:
        return True
    else:
        return False

def resolveModification(hUnimodDict, unimodDict, deNovoMass, refMass, deNovoSeq, unmodSeq, paramsDict, term = None, nModSymbol=None, epsilon=0.02):

    possUnmodMasses = [refMass]
    for modEntry in paramsDict['Static Mods']:
        if modEntry[1] in unmodSeq:
            possUnmodMasses += [refMass - float(paramsDict['Static Mods'][modEntry])]

    for modEntry in paramsDict['Diff Mods'].values():
        if modEntry[1] in unmodSeq or (term != None and modEntry[1].lower() == term.lower()):
            possUnmodMasses += [refMass + float(modEntry[2])]
            
    mods = []
    for unmodMass in possUnmodMasses:

        # see if iso substitution was masked by diff mod
        if abs(unmodMass - deNovoMass) < epsilon:
            return ([('Isobaric Substitution', 0, deNovoMass - unmodMass)], deNovoMass, deNovoSeq, unmodSeq)
        
        try:
            possMods = hUnimodDict[hashMass(deNovoMass - unmodMass)]
        except KeyError:
            continue

        for mod in possMods:
            if abs(possMods[mod]) < epsilon:
                locations = unimodDict[mod]['locations']
                for location in locations:
                    if isValidModLocation(unmodSeq, term, location):
                        mods += [(mod, unimodDict[mod]['mass'], possMods[mod])]
                        break

    if mods:
        return (mods, deNovoMass, deNovoSeq, unmodSeq)
    else:
        return False

# Returns a tuple of modifications found in the given interval
def resolveInterval(refMass, deNovoMass, refSubSeq, deNovoSubSeq, hUnimodDict, unimodDict, paramsDict, term=None, nModSymbol=None, epsilon=0.02):
    if not refSubSeq and not deNovoSubSeq:
        return tuple()
    elif not refSubSeq:
        return ([('Insertion', deNovoMass, 0)], deNovoMass, deNovoSubSeq, refSubSeq)
    elif not deNovoSubSeq:
        return ([('Deletion', -refMass, 0)], deNovoMass, deNovoSubSeq, refSubSeq)
    elif abs(deNovoMass - refMass) < epsilon:
        return ([('Isobaric Substitution', 0, deNovoMass - refMass)], deNovoMass, deNovoSubSeq, refSubSeq)
    else:
        matchedMods = resolveModification(hUnimodDict, unimodDict, deNovoMass, refMass, deNovoSubSeq, refSubSeq, paramsDict, term = term, nModSymbol=nModSymbol, epsilon=epsilon)
        if matchedMods:
            return matchedMods
        else:
            return ([('Undefined Mass Shift', None, deNovoMass - refMass)], deNovoMass, deNovoSubSeq, refSubSeq)

def loadFASTAs(fastaFile, fromBLAST=True):
    seqDict = {}
    seqGen = SFDT.sequenceGenerator(fastaFile)

    for seqName, seq in seqGen:
        seqNameTruncated = seqName.split(' ')[0][1:]
        if fromBLAST and len(seqNameTruncated) >= 63:
            seqNameTruncated = seqNameTruncated[:63] + '...'
        seqDict[seqNameTruncated] = seq

    return seqDict

def getSequence(fastaFile, reference):
    cmd = '/lab/seq_db/user/adevabhaktuni/ABBLAST/xdget -p %s %s' % (fastaFile, reference)

    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)

    result = proc.stdout.read()

    if len(result) == 0:
        return ''
    
    seq = ''
    for line in result.split('\n')[1:]:
        if len(line) > 0 and line[0] != '>':
            seq += line
        else:
            return seq

# Designed to throw a KeyError if uncharacterized amino acid found in extended_sequence
def trimExtendedSequence(extended_sequence, tag_interval, add_n=0, add_c=0, seq_delimiter='-'):

    base_seq = extended_sequence[tag_interval[0]:tag_interval[1]]
    #print 'Starting trim... Base', extended_sequence, base_seq, add_n, add_c
    n_index, c_index = 0, 0

    if add_n > 0:
        n_mass = 0
        for aa in extended_sequence[:tag_interval[0]][::-1]:
            if aa == seq_delimiter:
                break
            
            n_index += 1            
            n_mass += Constants.aminoacids[aa][2]
            if n_mass >= add_n:
                break

    if add_c > 0:
        c_mass = 0
        for aa in extended_sequence[tag_interval[1]:]:
            if aa == seq_delimiter:
                break
            
            c_index += 1
            c_mass += Constants.aminoacids[aa][2]
            if c_mass >= add_c:
                break
            
    #print n_index, c_index, tag_interval

    # Add one more amino acid to the N- and C- terminus so that sequence context can be read out later
    return extended_sequence[tag_interval[0]-n_index-1:tag_interval[1]+c_index+1]

# max counts is maximum number of matching sequences
# NOTE: This method only returns first found match in the case of a tie in match_length (i.e., one with shortest match_end_offset)
def getMaxMatchingSubstring(seq_index, de_novo_peptide, pept_length):
    match_counts, max_length_match, match_end_offset = 0, 0, []
    
    for end_offset in range(pept_length):
        # Don't need to keep going if it's impossible to find a longer match
        if (pept_length - end_offset) < max_length_match:
            break

        counts, match_length = seq_index.count_max_suffix(de_novo_peptide[:(pept_length-end_offset)])
        # TODO: Does not consider case of tie, want to return both end_indexes in this case?
        if match_length > max_length_match:
            max_length_match = match_length
            match_end_offset = [end_offset]
            match_counts = counts
        elif match_length == max_length_match:
            match_counts += counts
            match_end_offset += [end_offset]
            
    return match_counts, max_length_match, match_end_offset

def getCounts(counts_per_index):
    max_length = max([item[1] for item in counts_per_index])
    counts = 0
    indices = []

    for i, index_counts in enumerate(counts_per_index):
        if index_counts[1] == max_length:
            counts += index_counts[0]
            indices += [i]

    return counts, max_length, indices
                                                                                                                                                                                        
def getProteinName(protein_names, protein_offsets, position):
    return protein_names[position[0]][str(protein_offsets[position[0]][np.searchsorted(protein_offsets[position[0]], np.uint32(position[1])) - 1])]

def getPeptideContext(extended_sequence, start=None, end=None):
    return extended_sequence[start-1] + '.' + extended_sequence[start:end] + '.' + extended_sequence[end]

def getStartAndEndInds(proteinSeq, peptide):
    startInd = proteinSeq.index(peptide)
    return startInd, startInd+len(peptide)


if __name__ == '__main__':
    start_time = time.time()

    options = ArgLib.parse(['init', 'output', 'ppmstd', 'denovo', 'unimoddict', 'modmaxcounts', 'maxcounts', 'fmindex', 'scans'])

    # Amount of mass to add to ends of peptides when fetching sequence context for incomplete masses, may move to command line arguments later
    end_buffer = 250
    # Minimum AA mass, used to decide how many amino acids to request from sequence index when attempting to resolve an inexact match
    min_aa_mass = 50
    # Length of matched tags to tile across de novo and database sequences
    tag_length = 2
    # Maximum number of matches to return per scan
    max_candidates = 100

    # Set to Filtered mode
    filter_for_scans = False
    if options.scans != None:
        filter_for_scans = True
        scans = eval(options.scans)
        print 'Only sequencing scans', scans

    # Load mod annotation data(for de novo sequence results)
    params_dict = ArgLib.parseInitFile(options.init, options)

    # Load mod info from unimod
    with open(options.unimoddict) as fin:
        unimod_dict = pickle.load(fin)
    hashed_unimod_dict = hashUnimodDict(unimod_dict)

    index_basename = os.path.splitext(options.fmindex)[0]
    
    # Load positions of sequence seperators
    protein_names = []
    for seqnames_file in sorted(glob.glob(index_basename + '.seqnames*')):
        protein_names += [anydbm.open(seqnames_file)]

    with open(index_basename + '.offsets') as fin:
        protein_offsets = pickle.load(fin)

    # Load FM Index
    seq_indices = []
    for seqnames_file in sorted(glob.glob(index_basename + '.fm*')):
        seq_indices += [fm.FM_load(seqnames_file)]

    # Load de novo sequence info
    de_novo_cols, de_novo_results = DataFile.getScanInfo(options.denovo, delimiter='\t')
    # Initialize ambiguous edges to empty list if data is not present in de novo results (i.e., not running TAG-GRAPH on LADS)
    ambig_edges_present = 'Ambig Edges' in de_novo_cols
    ambig_edges = []

    # Prep output file
    outFile = open(options.output, 'w')
    cols = ['ScanF', 'Alignment Score', 'Rank', 'Context', 'Modifications', 'Proteins', 'Matching Tag Length', 'De Novo Peptide', 'Unmod De Novo Peptide', 'De Novo Score', 'Time Taken']
    if ambig_edges_present:
        cols.insert(-3, 'Ambig Edges')
    outFile.write('\t'.join(cols) + '\n')

    # Here we go!
    num_skipped = 0
    for entry in de_novo_results:

        if filter_for_scans and entry['ScanF'] not in scans:
            continue
        
        # For Timing Purposes
        t1 = time.time()
        
        scanData = defaultdict(lambda: None)
        scanData['ScanF'] = entry['ScanF']
        scanData['De Novo Peptide'] = entry['Peptide']
        scanData['Unmod De Novo Peptide'] = Constants.stripModifications(scanData['De Novo Peptide'], noRemove=[])
        scanData['De Novo Score'] = entry['ALC (%)']

        if ambig_edges_present:
            ambig_edges = eval(entry['Ambig Edges'])
            scanData['Ambig Edges'] = ambig_edges

        # Find matches in sequence index
        de_novo_peptide = scanData['Unmod De Novo Peptide']
        pept_length = len(de_novo_peptide)

        counts_per_index = []
        for i, seq_index in enumerate(seq_indices):
            # Each call to getMaxMatchingSubstring returns a tuple of the form counts, match_length, match_end_offset
            counts_per_index += [getMaxMatchingSubstring(seq_index, de_novo_peptide, pept_length)]

        # Find counts for longest match
        counts, match_length, indices_to_search = getCounts(counts_per_index)
        scanData['Matching Tag Length'] = match_length

        match_data = []
        match_list = {}

        if counts > options.maxcounts:
            # de novo peptide does not narrow down database candidates enough
            print 'Counts exceed max counts for peptide %s at scan number %s - counts: %i match length: %i '%(de_novo_peptide, scanData['ScanF'], counts, match_length)
            num_skipped += 1
        elif match_length == pept_length:
            # Found exact match in sequence index
            locations = []
            for i in indices_to_search:
                locations += [{'index_locations': seq_indices[i].locate(de_novo_peptide)[0], 'seq_index': i}]

            matches = []
            for location in locations:
                seq_index = location['seq_index']
                
                for index_location in location['index_locations']:
                    extended_sequence = seq_indices[seq_index].extract(index_location-1, index_location+pept_length+1)
                    if extended_sequence in match_list:
                        match_list[extended_sequence] += [(seq_index, index_location)]
                    else:
                        match_data += [([(0, pept_length), [], EXACT_MATCH_SCORE], extended_sequence)]
                        match_list[extended_sequence] = [(seq_index, index_location)]

        elif counts > options.modmaxcounts:
            # de novo peptide does not narrow down database candidates enough for mod search
            print 'Counts exceed mod max counts for inexact matching peptide %s at scan number %s - counts: %i match length: %i '%(de_novo_peptide, scanData['ScanF'], counts, match_length)
            num_skipped += 1
        else:

            # Did not find exact match, resolve mismatches with mods or substitutions
            # Need to store data on which start_ind, end_ind, and sequence index to use for each location
            locations = []
            for i in indices_to_search:
                for match_end_offset in counts_per_index[i][-1]:
                    match_start_ind, match_end_ind = pept_length - ( match_end_offset + match_length ), pept_length - match_end_offset
                    index_locations = seq_indices[i].locate(de_novo_peptide[ match_start_ind:match_end_ind ])[0]

                    locations += [{'index_locations': index_locations, 'seq_index': i, 'match_start_ind': match_start_ind, 'match_end_ind': match_end_ind}]

            # Get and score modifications at all matching locations
            de_novo_prmladder = Constants.getPRMLadder(scanData['De Novo Peptide'], ambigEdges=ambig_edges)
            for location in locations:
                add_n, add_c = 0, 0
                seq_index = location['seq_index']
                if location['match_start_ind'] > 0:
                    add_n = de_novo_prmladder[ location['match_start_ind'] ] + end_buffer
                if location['match_end_ind'] < pept_length:
                    add_c = de_novo_prmladder[-1] - de_novo_prmladder[ location['match_end_ind'] ] + end_buffer

                # Pull out sequence surrounding match and trim to match mass buffers
                # Retreive one more aa than necessary from termini to ensure that sequence context can be inferred
                for index_location in location['index_locations']:
                    add_aa_n, add_aa_c = int(add_n/min_aa_mass) + 1, int(add_c/min_aa_mass) + 1
                    if index_location - add_aa_n < 0:
                        add_aa_n = index_location
                    
                    extended_sequence = seq_indices[seq_index].extract( index_location-add_aa_n, index_location+match_length+add_aa_c)
                    # TODO: Is this even necessary? See the extra processing time to make sure mass addition constraints are satisfied this makes things faster or slower as opposed to just simply truncating at protein boundaries
                    extended_sequence = trimExtendedSequence(extended_sequence, (add_aa_n, len(extended_sequence)-add_aa_c ), add_n, add_c)

                    if extended_sequence not in match_list:
                        matches = alignDeNovoToDBSequence(de_novo_prmladder, de_novo_peptide, extended_sequence[1:-1], hashed_unimod_dict, unimod_dict, params_dict, deNovoAmbigEdges = ambig_edges, tagLength=tag_length, isobaricPenalty=0, defModPenalty=-1, inDelPenalty=-3, undefModPenalty=-3, defaultScore=0)
                        match_data += [(match, extended_sequence) for match in matches]
                        match_list[extended_sequence] = [(seq_index, index_location)]
                    else:
                        match_list[extended_sequence] += [(seq_index, index_location)]

        match_data.sort(key=lambda d: -d[0][2])
        scanData['Time Taken'] = time.time() - t1
        
        # TODO: Check if this is way to write data is SLOW
        for i in range(min(max_candidates, len(match_data))):
            match, extended_sequence = match_data[i]

            scanData['Rank'] = i+1
            scanData['Context'] = getPeptideContext(extended_sequence, match[0][0]+1, match[0][1]+1)
            scanData['Modifications'] = match[1]
            scanData['Proteins'] = [getProteinName(protein_names, protein_offsets, location) for location in match_list[extended_sequence]]
            scanData['Alignment Score'] = match[2]
            
            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

    outFile.close()
            
    print 'Finished. Total Time Taken: %f. Num skipped due to matches exceeding max_counts: %i' % (time.time() - start_time, num_skipped)

import sys
import os
import networkx as nx


import numpy as np
import time
import pickle
import glob
from collections import defaultdict, deque
import itertools
import anydbm
import copy
import math

CUR_DIR = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(CUR_DIR,'fmindex'))

import fmindex as fm
import ArgLib
import DataFile
import Constants
import ProbNetwork as PN
import Database
import Validator

###################
# KNOWN BUGS
# Adding arginine mass mod (+156) instead of matching an isobaric sub to arginine (due to scoring matching AA +1), penalize mods more?
# Allowing insertion and deletion of AAs instead of matching isobaric sub (due to scoring matching AA + 1), penalize mods more?
####################

EXACT_MATCH_SCORE = 100
inf = float("inf")

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
    
def uniqueModsKey(mods):
    mods = [mod for mod in mods if mod[0][0] != 'Isobaric Substitution']
    return Validator.getUniqueModTuple(mods)

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

def getMatchLength(path):
    return sum([tag[0][1] - tag[0][0] for tag in path])

# Returns list of mods in math (and all info required for scoring)
def getModData(tagGraph, path):
    return [tagGraph.edge[ path[i] ][ path[i+1] ]['mods'] + (path[i][1][1] - path[0][1][0],) for i in range(len(path) - 1) if len(tagGraph.edge[path[i] ][ path[i+1] ]['mods']) > 0]

# Returns a tuple of the form (refseq_start_ind, refseq_end_ind), mods
# Mods is a list of mods in the form ([List of matching mod annotations], mod_mass, de_novo_subseq, ref_subseq, ref_mod_ind)
# ref_mod_ind is index of mod location starting from start of returned peptide (i.e., 'true' start ind in dbPept (extended_sequence) is refseq_start_ind + ref_mod_ind
def alignDeNovoToDBSequence(deNovoPRMLadder, deNovoPept, dbPept, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges = None, tagLength=2, isobaricPenalty=-0.5, defModPenalty=-1, inDelPenalty=-2, undefModPenalty=-3, defaultScore=0, cutOffDiff=2, modTolerance=0.02, align_to_length_ratio = 0.2):
    #print deNovoPRMLadder

    #print 'De Novo', deNovoPept
    #print 'DB', dbPept
    
    dbPRMLadder = Constants.getPRMLadder(dbPept, considerTerminalMods = False, addEnds=True)
    #print dbPept, dbPRMLadder
    
    startTags, endTags = generateStartAndEndTags(deNovoPept, dbPept)
    sequenceTags = generateSequenceTags(deNovoPept, dbPept, tagLength=tagLength)

    #print 'Tags', sequenceTags
    tagGraph = getSequenceTagGraph(startTags, endTags, sequenceTags)

    topo_sort = nx.topological_sort(tagGraph)
    maxLScore, maxEndScoringTag = scoreAndInsertLScores(tagGraph, topo_sort, deNovoPRMLadder, deNovoPept, dbPRMLadder, dbPept, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges, tagLength, isobaricPenalty, defModPenalty, inDelPenalty, undefModPenalty, defaultScore, modTolerance)
    maxRScore, maxStartScoringTag = insertRScores(tagGraph, topo_sort, defaultScore)

    #print 'Scores', maxLScore, maxEndScoringTag, maxRScore, maxStartScoringTag
    #print sorted(tagGraph.nodes(data=True))
    #for edge in  sorted(tagGraph.edges(data=True), key=lambda k: -k[2]['edge_score']):
    #    print edge

    #if maxLScore/( len(deNovoPRMLadder) - 1 ) < align_to_length_ratio:
    #    cutOff = maxLScore
    #else:
    #    cutOff = maxLScore - cutOffDiff
    try:
        cutOff = maxLScore - cutOffDiff
    except TypeError:
        # TypeError occurs when no tags are found between database and de novo sequence (i.e., if a given db pept at a location has no two adjacent amino acids in common with de novo peptide
        # This is clearly an error case, unknown how often locate() returns an erroneous result and/if locate is not returning correct results (i.e., creating false negatives). Further testing is needed
        # For now, just throw error and print error message so this can be monitored
        raise TypeError
    
    max_scoring_paths = getBestAlignments(tagGraph, cutOff, startTags)
    pept_mod_data = [(path, getModData(tagGraph, path), score) for score, path in max_scoring_paths]

    return pept_mod_data

def scoreAndInsertLScores(tagGraph, topo_sort, deNovoPRMLadder, deNovoPept, dbPRMLadder, dbPept, hashedUnimodDict, unimodDict, paramsDict, deNovoAmbigEdges = None, tagLength=2, isobaricPenalty=-0.5, defModPenalty=-1, inDelPenalty=-2, undefModPenalty=-3, defaultScore=0, modTolerance=0.02):

    maxScore, maxScoringTag = None, None
    num_denovo_prms = len(deNovoPRMLadder)
    
    for tag in topo_sort:
        nodeScore = tag[0][1] - tag[0][0]
        #print 'Tag', tag
        for prevTag in tagGraph.predecessors(tag):
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

            # If there are static mods on the peptide termini, and we are trying to match a disagreement region which spans one of the de novo peptide termini
            # add the static mass to the reference sequence mass here
            if prevTag[0][1] == 0:
                refMass += Constants.NTermMods['static']
            if tag[0][0] == num_denovo_prms:
                refMass += Constants.CTermMods['static']
                
            #print refMass, deNovoMass, refSubSeq, deNovoSubSeq
            mods = resolveInterval(refMass, deNovoMass, refSubSeq, deNovoSubSeq, hashedUnimodDict, unimodDict, paramsDict, term=term, modTolerance=modTolerance)

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

def isValidModLocation(unmodSeq, term, location, static_mod_aas = [], static_aa = False):
    if not static_aa and location[0] in static_mod_aas:
        return False

    if static_aa and location[0] != static_aa:
        return False
    
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


def getModLocations(unmodSeq, term, locations):
    indexes = []

    for location in locations:
        
        if location[1] == 'Anywhere':
            indexes += [(index, location) for index in list(findAll(location[0], unmodSeq))]
        elif location[1] == term:
            if term == 'N-term' and (unmodSeq[0] == location[0] or term == location[0]):
                indexes += [(0, location)]
            elif term == 'C-term' and (unmodSeq[-1] == location[0] or term == location[0]):
                indexes += [(len(unmodSeq) - 1, location)]

    return set(indexes)
            

def resolveModification(hUnimodDict, unimodDict, deNovoMass, refMass, deNovoSeq, unmodSeq, paramsDict, term = None, modTolerance=0.02):

    # TODO: Add reporting of what static mods were assumed to be substracted (and diff mods assumed to be added) to mod annotation
    # Substract static mods from unmod mass (covers case where AA mod blocked static mod, happens with cysteines sometimes)
    possUnmodMasses = [(refMass, 0, False)]

    # Sometimes, a static modification in the search parameters can 'mask' a modification on that residue due to other means.
    # This logic searches these possibilities explicitly
    for modEntry in paramsDict['Static Mods']:
        if modEntry[1] in unmodSeq:
            possUnmodMasses += [( refMass - float(paramsDict['Static Mods'][modEntry]), float(paramsDict['Static Mods'][modEntry]), modEntry[1] )]

    if Constants.NTermMods['static'] > 0 and term == 'N-term':
        possUnmodMasses += [( refMass - Constants.NTermMods['static'], Constants.NTermMods['static'], 'N-term')]
    if Constants.CTermMods['static'] > 0 and term == 'C-term':
        possUnmodMasses += [( refMass - Constants.CTermMods['static'], Constants.CTermMods['static'], 'C-term')]

    # TODO: Add this back in if many mods are being masked by close proximity to user specified differential mods (should now show up as an undefined mod)
    # Reference sequence may have differential mods on amino acids which may mask true mod identity on another residue
    # This case mainly arises with LADS samples that have N-terminal dimethylation, may be a source of error in many other cases
    #for modEntry in paramsDict['Diff Mods'].values():
    #    if modEntry[1] in unmodSeq or (term != None and modEntry[1].lower() == term.lower()):
    #        possUnmodMasses += [refMass + float(modEntry[2])]

    #print [deNovoMass - unmodMass[0] for unmodMass in possUnmodMasses], unmodSeq
            
    mods = []
    for unmodMass in possUnmodMasses:

        # AD - commented this out for now for easier scoring, seems to only be useful for LADS samples
        # see if iso substitution was masked by diff mod
        #if abs(unmodMass - deNovoMass) < epsilon:
        #    return ([('Isobaric Substitution', 0, deNovoMass - unmodMass)], deNovoMass, deNovoSeq, unmodSeq)
        
        try:
            possMods = hUnimodDict[hashMass(deNovoMass - unmodMass[0])]
        except KeyError:
            continue

        for mod in possMods:
            valid_locs = []
            if abs(possMods[mod]) < modTolerance:
                locations = unimodDict[mod]['locations']
                #print mod, locations, unmodMass
                for location in locations:
                    if isValidModLocation(unmodSeq, term, location, static_mod_aas = paramsDict['Static AAs'], static_aa = unmodMass[2]):
                        valid_locs += [location]

                # Subtracts mass of assumed static mod from mass of modification (mod found in this case assumed to be on unmodified residue where static mod not present)
                # TODO: Return all Mods or just mod with minimum mass error?
                if len(valid_locs) > 0:
                    mods += [(mod, unimodDict[mod]['mass'] - unmodMass[1], possMods[mod], valid_locs)]

    if mods:
        return (mods, deNovoMass, deNovoSeq, unmodSeq)
    else:
        return False

# Returns a tuple of modifications found in the given interval
def resolveInterval(refMass, deNovoMass, refSubSeq, deNovoSubSeq, hUnimodDict, unimodDict, paramsDict, term=None, modTolerance=0.02):
    if not refSubSeq and not deNovoSubSeq:
        return tuple()
    elif not refSubSeq:
        return ([('Insertion', deNovoMass, 0)], deNovoMass, deNovoSubSeq, refSubSeq)
    elif not deNovoSubSeq:
        return ([('Deletion', -refMass, 0)], deNovoMass, deNovoSubSeq, refSubSeq)
    elif abs(deNovoMass - refMass) < modTolerance:
        return ([('Isobaric Substitution', 0, deNovoMass - refMass)], deNovoMass, deNovoSubSeq, refSubSeq)
    else:
        matchedMods = resolveModification(hUnimodDict, unimodDict, deNovoMass, refMass, deNovoSubSeq, refSubSeq, paramsDict, term = term, modTolerance=modTolerance)
        if matchedMods:
            return matchedMods
        else:
            # TODO: Also find a way to report when undefined mod lies on AA with static mod (may be an undefined mod on unmodified form)?
            return ([('Undefined Mass Shift', deNovoMass - refMass, None)], deNovoMass, deNovoSubSeq, refSubSeq)


def getDiffModMap(diff_mods):
    mod_map = {'AA': {}, 'N-term': {}, 'C-term': {} }
    for mod_symbol in diff_mods:
        mod_data = diff_mods[mod_symbol]
        if mod_data[1] == 'N-term' or mod_data[1] == 'C-term':
            mod_map['N-term'][mod_symbol] = ((mod_data[0], float(mod_data[2]), 0), (mod_data[1], mod_data[1]))
        else:
            for aa in mod_data[1]:
                mod_map['AA'][aa + mod_symbol] = ((mod_data[0], float(mod_data[2]), 0), (aa, 'Anywhere'))

    return mod_map

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
    #print tag_interval[0]-n_index-1, tag_interval[1]+c_index+1, tag_interval
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

def getModifiedSubsequenceVariants(mod, term):
    modified_sequences = {}
    
    for mod_type in mod[0]:
        if mod_type[0] == 'Undefined Mass Shift':
            modified_sequences[mod_type] = []
            for index in range(len(mod[-2])):
                modSeq = getModSeq(mod[-2], mod_type[1], index)
                if modSeq:
                    modified_sequences[mod_type] += [(modSeq, (index, (mod[-2][index], 'Anywhere')) )]
        else:
            mod_hash = mod_type[:-1]
            modified_sequences[mod_hash] = []
            index_locs = getModLocations(mod[-2], term, mod_type[-1])
            for index_loc in index_locs:
                modSeq = getModSeq(mod[-2], mod_type[1], index_loc[0], min_mod_AA_mass = 10)
                if modSeq:
                    modified_sequences[mod_hash] += [(modSeq, index_loc)]

    return modified_sequences

def getModSeq(unmodSeq, modMass, index, modAASymbol = '-', min_mod_AA_mass = 10):
    modAA = unmodSeq[index]
    s = list(unmodSeq)
    s[index] = modAASymbol

    if Constants.aminoacids[modAA][2] + modMass > min_mod_AA_mass:
        return (''.join(s), [(0, Constants.aminoacids[modAA][2] + modMass)])
    else:
        return False

def getSeqFromDeNovoTag(de_novo_aas, start_ind, end_ind, context_start_ind, diff_mod_map, modAASymbol = '-'):
    de_novo_seq = []
    mod_ranges = []
    mod_ambig_edges = []
    mods = []

    #print diff_mod_map, de_novo_aas, start_ind, end_ind
    for i, aa in enumerate(de_novo_aas['AA'][start_ind:end_ind]):
        if len(aa) > 1 or (i + context_start_ind == 0 and de_novo_aas['N-term'] != None) or (i + context_start_ind == len(de_novo_aas['AA']) and de_novo_aas['C-term'] != None):
            de_novo_seq += [modAASymbol]
            mod_ranges += [(context_start_ind+i, context_start_ind+i+1)]
            
            mod_aa_mass = Constants.aminoacids[aa][2]
            if len(aa) > 1:
                mods += [diff_mod_map['AA'][aa] + (context_start_ind+i,)]
            
            if i + context_start_ind == 0 and de_novo_aas['N-term'] != None:
                mod_aa_mass += Constants.NTermMods[de_novo_aas['N-term']]
                mods += [diff_mod_map['N-term'][de_novo_aas['N-term']] + (context_start_ind+i,)]
            if i + context_start_ind == len(de_novo_aas['AA']) and de_novo_aas['C-term'] != None:
                mod_aa_mass += Constants.CTermMods[de_novo_aas['C-term']]
                mods += [diff_mod_map['C-term'][de_novo_aas['C-term']] + (context_start_ind+i,)]

            mod_ambig_edges += [(0.0, mod_aa_mass)]
        else:
            de_novo_seq += [aa]

    return ''.join(de_novo_seq), mods, mod_ambig_edges, mod_ranges

def expandMatch(peptide, mods, path, de_novo_aas, diff_mod_map):
    #print '----Expand match-----', peptide, mods, de_novo_aas

    if len(mods) == 0:
        # Return de novo sequence in case of exact match, so that all diff mods are properly placed in peptide
        de_novo_seq, de_novo_mods, de_novo_ambig_edges, de_novo_mod_ranges = getSeqFromDeNovoTag(de_novo_aas, 0, len(de_novo_aas['AA']), 0, diff_mod_map)
        return [(de_novo_seq, de_novo_ambig_edges, de_novo_mods)], de_novo_mod_ranges

    # Enumerate sequence variants for annotated mods for scoring

    # Each item in the replace_subseqs list has 2 elements:
    # (start_ind, end_id) - start and end index of region in peptide to replace
    # list of peptide subsequences to replace region with. Each peptide item has 3 elements:
    # replacement_sequence - this is a two-member tuple of the form (sequence, ambig_edges). The sequence will have an 'X' where the modification is found (if not a deletion, insertion, or isobaric substitution)
    # mod_name - name of modification
    # index_of_mod - location of modification in parent peptide, used for mod reporting annotation later on
    mod_ranges = []
    all_replacement_pepts = []
    for mod in mods:
        replacement_pepts = []
        if mod[0][0][0] == 'Insertion':
            # TODO: Currently, inserted sequences do not preserve diff mods on de novo peptide subsequence, WILL affect scoring
            mod_range = (mod[-1], mod[-1])
            replacement_pepts = [( (mod[-3], []), (mod[0][0], mod[-3], mod[-1]) )]
            
        elif mod[0][0][0] == 'Deletion':
            mod_range = (mod[-1], mod[-1] + len(mod[-2]))
            replacement_pepts = [( ('', []), (mod[0][0], mod[-2], mod[-1]) )]
            
        elif mod[0][0][0] == 'Isobaric Substitution':
            mod_range = (mod[-1], mod[-1] + len(mod[-2]))
            replacement_pepts = [( (mod[-2], []), (mod[0][0], tuple(), mod[-1]) )]
            
        else:
            # Figure out if mod occurs on terminus
            if mod[-1] == 0:
                term = 'N-term'
            elif mod[-1] + len(mod[-2]) == len(peptide):
                term = 'C-term'
            else:
                term = None
            
            modified_sequences = getModifiedSubsequenceVariants(mod, term)
            mod_range = (mod[-1], mod[-1] + len(mod[-2]))
            replacement_pepts = []
            #print modified_sequences
            for mod_type in modified_sequences:
                for sequence in modified_sequences[mod_type]:
                    replacement_pepts += [( sequence[0], (mod_type, sequence[1][1], sequence[1][0] + mod[-1]) )]


        if len(replacement_pepts) > 0:
            mod_ranges += [mod_range]
            all_replacement_pepts += [replacement_pepts]
        else:
            # One of the mods is invalid, return
            # Happens when mass of mod + mass of modified AA is less than the minimum mass for a modified amino acid (currently set at 10 Daltons)
            # TODO: exclude these matches at the graph candidate enumeration step instead of here? (don't draw edge if mod is implausible)
            #print 'ERROR: mod invalid for peptide %s: %s. Will not expand match.'%(peptide,str(mod)) 
            return [], []

    
    # Mapping of shared tags to de novo sequences (to properly account for differential mod annotations in de novo peptide)
    de_novo_seq_map = {}
    return_mod_ranges = []
    i = 0
    for node in path:
        de_novo_seq, de_novo_mods, de_novo_ambig_edges, de_novo_mod_ranges = getSeqFromDeNovoTag(de_novo_aas, node[0][0], node[0][1], node[1][0] - path[0][1][0], diff_mod_map)
        
        # Normalize indices of tag to start and end of peptide context (previously index referred to index relative to extended_sequence)
        de_novo_seq_map[(node[1][0] - path[0][1][0], node[1][1] - path[0][1][0])] = (de_novo_seq, de_novo_mods, de_novo_ambig_edges)

        # Some nodes might not yield a modification (i.e., connection between start node and first shared tag)
        if i < len(mod_ranges) and node[1][0]-path[0][1][0] == mod_ranges[i][1]:
            return_mod_ranges += [mod_ranges[i]]
            i += 1
        return_mod_ranges += de_novo_mod_ranges
    
    # Now create list of peptides from all expanded mod regions
    expanded_peptide_list = []
    for mod_chain in itertools.product(*all_replacement_pepts):
        new_pept, ambig_edges, mods = '', [], []
        prev_mod_range = (0,0)
        for i, mod in enumerate(mod_chain):
            mod_range = mod_ranges[i]

            # Add data for shared tag region (between de novo and db pept)
            de_novo_seq, de_novo_mods, de_novo_ambig_edges = de_novo_seq_map[(prev_mod_range[1], mod_range[0])]
            new_pept = new_pept + de_novo_seq + mod[0][0]
            mods += de_novo_mods
            ambig_edges += de_novo_ambig_edges

            # Add mod data
            ambig_edges += mod[0][1]
            mods += [mod[1]]
            prev_mod_range = mod_range

        # Cap peptide with remaining shared tag region
        de_novo_seq, de_novo_mods, de_novo_ambig_edges = de_novo_seq_map[(prev_mod_range[1], len(peptide))]
        new_pept = new_pept + de_novo_seq
        mods += de_novo_mods
        ambig_edges += de_novo_ambig_edges

        # Add peptide to candidate list
        expanded_peptide_list += [(new_pept, ambig_edges, mods)]

    return expanded_peptide_list, return_mod_ranges

def scoreMatch(spec, charge, peptide, ambig_edges=[], mods=[], scanF = -1, no_spec_default_score = -5):
    score = 0

    try:
        pm = Constants.getPM(peptide, ambigEdges=ambig_edges)
    except KeyError:
        # Happens when insertion contains ambiguous amino acid
        #print 'ERROR: mod invalid for peptide %s, %s, %s. Will not score'%(peptide, str(ambig_edges), str(mods))
        return -10, 0

    if not spec:
        return no_spec_default_score, pm
        
    for node in Constants.nodeInfoGen(peptide, considerTerminalMods = True, ambigEdges=ambig_edges):
        score += spec.getNodeScore( node, pm, str( int(max(min(float(charge), 4), 2) )) )

    """
    for mod in mods:
        if mod[0][0] == 'Insertion' or mod[0][0] == 'Deletion':
           score += indel_penalty
        elif mod[0][0] == 'Undefined Mass Shift':
            score += undefined_penalty
        elif mod[0][0] == 'Isobaric Substitution':
            score += isobaric_penalty
        else:
            score += def_mod_penalty
    """

    if score == inf:
        # Happens sometimes when noise model returns 0 prob and spectrum model returns positive prob. Unknown why this occurs, must monitor
        print "CRITICAL ERROR: Score of infinity returned for peptide %s with ambiguous edges %s in scan %i"%(peptide, str(ambig_edges), scanF)
        return -10, pm
    else:
        return score/len(peptide), pm

def byAlignAndProteaseSpecificity(match_data):
    match, extended_sequence = match_data
    
    context = getPeptideContext(extended_sequence, match[0][0][1][0]+1, match[0][-1][1][0]+1)
    specificity = Database.getProteaseContext(context, params_dict['Enzyme']['specificity'])
    if specificity == 'full':
        spec_num = 0
    elif specificity == 'nonspecific':
        spec_num = 2
    else:
        spec_num = 1

    return (-match[2], spec_num)

# retreives the results whose scores match the max_composite_score for a given (context, mod_key) pair.
# Sometimes, duplicates can be returned by tag-graph for a variety of reasons:
# (1) two locations have different extended_sequences, but the actual contexts are identical after alignment
# (2) tag_graph can return the same result multiple times. This occurs if tag-graph decides to 'skip' using a matching tag of minimal length, thereby extending a mod range and adding in an isobaric substitution. These results should have different max_composite_scores than the top scoring result and not be considered.
# TODO: Is there ANY situation in which two results would have the same context, mod context, and mod key, but in fact represent distinct results and should be written separately?
def getMaxUniqueResults(candidates, write_data, prob_scores, pms, composite_scores, expanded_match_map, max_composite_score):
    mod_context_map = {}

    for i, expanded_match in enumerate(candidates):
        if round(composite_scores[i], 2) == max_composite_score:
            j = expanded_match_map[i]
            if expanded_match[0] not in mod_context_map:
                mod_context_map[expanded_match[0]] = [ expanded_match, prob_scores[i], pms[i], composite_scores[i], set(write_data[j]['Proteins']), j ]
            else:
                # fixes issue where protein localizations with same context but different extended sequences get split into two separate lines
                mod_context_map[expanded_match[0]][4] |= set(write_data[j]['Proteins'])

    return mod_context_map.values()


if __name__ == '__main__':
    start_time = time.time()

    options = ArgLib.parse(['init', 'output', 'dtadir', 'model', 'config', 'ppmstd', 'denovo', 'unimoddict', 'modmaxcounts', 'modtolerance', 'maxcounts', 'fmindex', 'scans', 'excludescans'], optArgs = [{'opts': ('-x', '--splittaxon'), 'attrs': {'dest': 'splittaxon', 'action': 'store_true', 'default': False, 'help': 'Flag. For searches of metaproteomic databases, split identical context entries by taxon for accurate consideration via EM.'}}])

    # Amount of mass to add to ends of peptides when fetching sequence context for incomplete masses, may move to command line arguments later
    end_buffer = 250
    # Minimum AA mass, used to decide how many amino acids to request from sequence index when attempting to resolve an inexact match
    min_aa_mass = 50
    # Length of matched tags to tile across de novo and database sequences
    tag_length = 2
    # Maximum number of matches to return per scan in alignment step
    max_candidates_alignment = 5
    # Maximum  numbe of matches to return per scan in scoring step
    # max_candidates = 100
    # Default score if spectrum is not present (for MzXMLToSearch, occurs when spectrum has less than 5 peaks)
    NO_SPEC_DEFAULT_SCORE = -5
    # Minimum alignment score to peptide length ratio under which results are not enumerated exhaustively
    MIN_ALIGN_TO_LENGTH_RATIO = 0.2
    # Minimum spectrum score required to continue expanding candidates (helps to not spend unnecessary time enumerating likely incorrect candidates)
    MIN_SPEC_SCORE_TO_CONTINUE = -1

    # max difference between max score and score under consideration to expand match (past max_candidates)
    # RATIONALE FOR SETTING AT 2: All other things equal (i.e., # of defined, undefined mods, etc., this lets us return candidates which choose to bypass the use of a single matching tag of length 2)
    cut_off_diff = 2

    # Set to Filtered mode
    filter_for_scans = False
    if options.scans != None:
        filter_for_scans = True
        scans = set(eval(options.scans))
        print 'Only sequencing scans', scans

    exclude_scans = set()
    if options.excludescans != None:
        exclude_scans = set(eval(options.excludescans))

    # Load mod annotation data(for de novo sequence results)
    params_dict = ArgLib.parseInitFile(options.init, options)
    diff_mod_map = getDiffModMap(params_dict['Diff Mods'])
    #print params_dict['Diff Mods']
    #print diff_mod_map

    # Load Probabilistic Scoring Model
    PNet = PN.ProbNetwork(options.config, options.model)

    # Load mod info from unimod
    with open(options.unimoddict) as fin:
        unimod_dict = pickle.load(fin)
    hashed_unimod_dict = hashUnimodDict(unimod_dict, maxEp = options.modtolerance)

    index_basename = os.path.splitext(options.fmindex)[0]

    # Load scan data
    dtaList = glob.glob(options.dtadir + os.path.sep + '*.dta')
    scanFDict = DataFile.getScanFDict(dtaList)
    
    # Load positions of sequence seperators
    protein_names = []
    for seqnames_file in sorted(glob.glob(index_basename + '.seqnames*'), key = lambda k: int(os.path.splitext(k)[1][1:])):
        protein_names += [anydbm.open(seqnames_file)]

    with open(index_basename + '.offsets') as fin:
        protein_offsets = pickle.load(fin)

    # Load FM Index
    seq_indices = []
    for seqnames_file in sorted(glob.glob(index_basename + '.fm*'), key = lambda k: int(os.path.splitext(k)[1][1:])):
        seq_indices += [fm.FM_load(seqnames_file)]

    # Load de novo sequence info
    de_novo_cols, de_novo_results = DataFile.getScanInfo(options.denovo, delimiter='\t')
    # Initialize ambiguous edges to empty list if data is not present in de novo results (i.e., not running TAG-GRAPH on LADS)
    ambig_edges_present = 'Ambig Edges' in de_novo_cols
    ambig_edges = []

    # Prep output file
    outFile = open(options.output, 'w')
    cols = ['ScanF', 'RT', 'Charge', 'Theo M+H', 'Obs M+H', 'PPM', 'Alignment Score', 'Spectrum Probability Score', 'Composite Score', 'Rank', 'Context', 'Mod Context', 'Mod Ambig Edges', 'Match Modifications', 'Mod Ranges', 'Proteins', 'Matching Tag Length', 'De Novo Peptide', 'Unmod De Novo Peptide', 'De Novo Score', 'Num Matches', 'Time Taken']
    if ambig_edges_present:
        cols.insert(-3, 'De Novo Ambig Edges')
    outFile.write('\t'.join(cols) + '\n')

    # Here we go!
    num_skipped = 0
    for entry in de_novo_results:
        # print entry['ScanF']
        if filter_for_scans and entry['ScanF'] not in scans:
            continue

        if entry['ScanF'] in exclude_scans:
            continue
        
        # For Timing Purposes
        t1 = time.time()

        scanF = int(entry['ScanF'])
        scanData = defaultdict(lambda: None)
        scanData['ScanF'] = entry['ScanF']
        scanData['De Novo Peptide'] = entry['Peptide']
        scanData['Unmod De Novo Peptide'] = Constants.stripModifications(scanData['De Novo Peptide'], noRemove=[], ambigAA='@')
        scanData['De Novo Score'] = entry['ALC (%)']
        scanData['Obs M+H'] = float(entry['Obs M+H'])
        scanData['Charge'] = entry['Charge']
        scanData['RT'] = entry['RT']

        if ambig_edges_present:
            ambig_edges = eval(entry['Ambig Edges'])
            scanData['De Novo Ambig Edges'] = ambig_edges

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
        scanData['Num Matches'] = counts

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

            for location in locations:
                seq_index = location['seq_index']
                
                for index_location in location['index_locations']:
                    extended_sequence = seq_indices[seq_index].extract(index_location-1, index_location+pept_length+1)
                    if extended_sequence in match_list:
                        match_list[extended_sequence] += [(seq_index, index_location)]
                    else:
                        match_data += [([ [((0,0),(0,0)), ((pept_length, pept_length),(pept_length, pept_length))], [], match_length ], extended_sequence)]
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
            #print scanData['De Novo Peptide'], ambig_edges
            de_novo_prmladder = Constants.getPRMLadder(scanData['De Novo Peptide'], ambigEdges=ambig_edges, ambigAA='@')
            #print scanData['ScanF'], scanData['De Novo Peptide'], de_novo_prmladder
            for location in locations:
                # print 'LOCATION', location
                add_n, add_c = 0, 0
                seq_index = location['seq_index']
                if location['match_start_ind'] > 0:
                    add_n = de_novo_prmladder[ location['match_start_ind'] ] + end_buffer
                if location['match_end_ind'] < pept_length:
                    add_c = de_novo_prmladder[-1] - de_novo_prmladder[ location['match_end_ind'] ] + end_buffer

                # Pull out sequence surrounding match and trim to match mass buffers
                # Retreive one more aa than necessary from termini to ensure that sequence context can be inferred
                for index_location in location['index_locations']:
                    add_aa_n, add_aa_c = int(add_n/min_aa_mass) + 2, int(add_c/min_aa_mass) + 2
                    #print add_n, add_aa_n, add_c, add_aa_c
                    
                    # TODO: How much faster is a try: catch: block than these checks that the sequence lies within the DB bounds (these are checks for INCREDIBLY rare errors)
                    # Only occurs if match is at start of DB
                    if index_location - add_aa_n < 0:
                        add_aa_n = index_location

                    # print 'Extracting extended sequence', index_location-add_aa_n, index_location+match_length+add_aa_c
                    extended_sequence = seq_indices[seq_index].extract( index_location-add_aa_n, index_location+match_length+add_aa_c)

                    # Only occurs if sequence matches to the end of the last DB entry
                    if add_aa_n + match_length + add_aa_c > len(extended_sequence):
                        add_aa_c = len(extended_sequence) - (match_length + add_aa_n)

                    # TODO: Is this even necessary? See the extra processing time to make sure mass addition constraints are satisfied this makes things faster or slower as opposed to just simply truncating at protein boundaries
                    # print 'Extended', index_location, extended_sequence, (add_aa_n, len(extended_sequence)-add_aa_c )
                    extended_sequence = trimExtendedSequence(extended_sequence, (add_aa_n, len(extended_sequence)-add_aa_c ), add_n, add_c)
                    # print 'Trimmed', extended_sequence
                    
                    if extended_sequence not in match_list:
                        try:
                            matches = alignDeNovoToDBSequence(de_novo_prmladder, de_novo_peptide, extended_sequence[1:-1], hashed_unimod_dict, unimod_dict, params_dict, deNovoAmbigEdges = ambig_edges, tagLength=tag_length, isobaricPenalty=0, defModPenalty=-1, inDelPenalty=-3, undefModPenalty=-3, defaultScore=0, modTolerance=options.modtolerance, align_to_length_ratio = MIN_ALIGN_TO_LENGTH_RATIO, cutOffDiff = cut_off_diff)
                        except TypeError:
                            print 'CRITICAL ERROR: No matching tags found between de novo peptide %s and db peptide %s at location %i for scan number %s'%(de_novo_peptide, extended_sequence[1:-1], index_location, scanData['ScanF'])
                            continue
                        
                        match_data += [(match, extended_sequence) for match in matches]
                        match_list[extended_sequence] = [(seq_index, index_location)]
                    else:
                        match_list[extended_sequence] += [(seq_index, index_location)]

        if len(match_data) > 0:

            if scanF in scanFDict:
                #precMass = float(entry['Obs M+H'])
                epSTD = options.ppmstd * scanData['Obs M+H'] * 10**-6
                spec = PN.Spectrum(PNet, scanData['Obs M+H'], Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=DataFile.getMassIntPairs(scanFDict[scanF]['dta']), useMemo=True)
                spec.initializeNoiseModel()
            else:
                spec = False
                print 'No corresponding dta file found for scan %i'%scanF

            # Expand and score matches        
            # TODO: Check if this is way to write data is SLOW
            expanded_match_map = defaultdict(list)
            expanded_matches = defaultdict(list)
            write_data = []
            de_novo_aas = Constants.getAllAAs(scanData['De Novo Peptide'], ambigEdges=ambig_edges, ambigAA='@')
            
            match_data.sort(key = byAlignAndProteaseSpecificity)

            i = 0
            max_score = match_data[0][0][2]
            align_ratio_over_cut = float(max_score) / len(scanData['Unmod De Novo Peptide']) >= MIN_ALIGN_TO_LENGTH_RATIO
            #print 'Ratio', align_ratio_over_cut, max_score, len(scanData['Unmod De Novo Peptide']), max_score / len(scanData['Unmod De Novo Peptide'])

            # This step is the slowest so we use a couple of heuristics to abort expansion of candidates if it doesn't look like it will yield a right answer
            # align_ratio_over_cut: is the ratio of the alignment score to the length of the de novo peptide over a defined cutoff (MIN_ALIGN_TO_LENGTH_RATIO)
            # is the max returned spectrum score for a given match over a defined cutoff (MIN_SPEC_SCORE_TO_CONTINUE)
            # We still write a few candidates in case of a likely wrong answer, as these will be used to define a null distribution during the EM phase
            # NOTE: These cutoffs need to be extremely liberal! We do NOT want to abort expansion prematurely and miss out on a right answer because of it!
            mod_ranges = []
            prob_scores = defaultdict(list)
            composite_scores = defaultdict(list)
            pms = defaultdict(list)
            while i < len(match_data) and ( (align_ratio_over_cut and match_data[i][0][2] == max_score) or i < max_candidates_alignment):
                match, extended_sequence = match_data[i]
                context = getPeptideContext(extended_sequence, match[0][0][1][0]+1, match[0][-1][1][0]+1)
                match_candidates, mod_range = expandMatch(context[2:-2], match[1], match[0], de_novo_aas, diff_mod_map)

                try:
                    new_prob_scores, new_pms = zip(*[scoreMatch(spec, scanData['Charge'], candidate[0], candidate[1], candidate[2], scanF, no_spec_default_score = NO_SPEC_DEFAULT_SCORE) for candidate in match_candidates])
                except ValueError:
                    # Occurs when zero match candidates are returned (i.e. putative mod annotation was invalid)
                    new_prob_scores, new_pms = [], []

                scanData['Alignment Score'] = match[2]                
                scanData['Proteins'] = [getProteinName(protein_names, protein_offsets, location) for location in match_list[extended_sequence]]            
                scanData['Context'] = context
                write_data += [copy.copy(scanData)]

                # Separate matches by (context, mod_key) pair. Will only write top scoring candidates for each set of candidates with the same (context, mod_key)
                for j, expanded_match in enumerate(match_candidates):
                    mod_key = uniqueModsKey(expanded_match[2])
                    expanded_matches[(context, mod_key)] += [expanded_match]
                    expanded_match_map[(context, mod_key)] += [i]
                    prob_scores[(context, mod_key)] += [new_prob_scores[j]]
                    composite_scores[(context, mod_key)] += [match[2] + new_prob_scores[j]]
                    pms[(context, mod_key)] += [new_pms[j]]
                    
                mod_ranges += [mod_range]
                i += 1

                # stop enumerating candidates if spectrum score of returned hits for top match is too low (top matches are most likely to be correct based on alignment score)
                if len(new_prob_scores) > 0 and max(new_prob_scores) < MIN_SPEC_SCORE_TO_CONTINUE:
                    break

            # Sort results by composite score and also compute maximum composite score for a particular (context, mod_key) pair
            max_composite_scores = sorted([( item_key, round(max(composite_scores[item_key]), 2) ) for item_key in composite_scores], key = lambda k: -k[1])
            #print max_composite_scores
            #print expanded_matches
            #print prob_scores
            #print composite_scores
            #print expanded_match_map
            #print write_data
            # write results
            time_taken = time.time() - t1
            i = 0

            #for item_key in max_composite_scores:
            #    print item_key, expanded_matches[item_key[0]], expanded_match_map[item_key[0]], write_data[expanded_match_map[item_key[0]][0]]
            
            for item_key, max_composite_score in max_composite_scores:
                results = getMaxUniqueResults(expanded_matches[item_key], write_data, prob_scores[item_key], pms[item_key], composite_scores[item_key], expanded_match_map[item_key], max_composite_score)
                #print 'ITEM', item_key, results
                for match, spectrum_score, pm, composite_score, proteins, write_data_ind in results:
                    scanData = copy.deepcopy(write_data[write_data_ind])
                    scanData['Rank'] = i + 1
                    scanData['Time Taken'] = time_taken
                    scanData['Match Modifications'] = match[2]
                    scanData['Mod Ambig Edges'] = match[1]
                    scanData['Spectrum Probability Score'] = spectrum_score
                    scanData['Composite Score'] = composite_score
                    scanData['Mod Ranges'] = mod_ranges[write_data_ind]
                    scanData['Proteins'] = list(proteins)
                    scanData['Theo M+H'] = pm + Constants.mods['H2O'] + Constants.mods['H+']
                    scanData['PPM'] = (scanData['Obs M+H'] - scanData['Theo M+H'])/scanData['Theo M+H'] * 1000000
    
                    scanData['Mod Context'] = scanData['Context'][:2] + match[0] + scanData['Context'][-2:]

                    if not options.splittaxon:
                        outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')
                    else:
                        #print scanData['Context'], match, proteins
                        taxon_prot_map = defaultdict(list)
                        for protein in scanData['Proteins']:
                            taxon_prot_map[protein.split('TAXON=')[1]] += [protein]

                        #print taxon_prot_map.keys()
                        for taxon in taxon_prot_map:
                            scanData['Proteins'] = taxon_prot_map[taxon]
                            outFile.write('\t'.join([str(scanData[col]) for col in cols]) + '\n')

                    i += 1
                
                
    outFile.close()
            
    print 'Finished. Total Time Taken: %f. Num skipped due to matches exceeding max_counts: %i' % (time.time() - start_time, num_skipped)

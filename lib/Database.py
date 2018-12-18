import DataFile
import ArgLib
import Constants

import dbm
import pickle
import numpy as np
from collections import defaultdict
import decimal
decimal.getcontext().prec = 100

import random
import re

# set both leucine and isoleucine codons to translate to isoleucine
codonDict = {

'TTT': 'F',
'TTC': 'F',
'TTA': 'L',
'TTG': 'L',
'CTT': 'L',
'CTC': 'L',
'CTA': 'L',
'CTG': 'L',
'ATT': 'I',
'ATC': 'I',
'ATA': 'I',
'ATG': 'M',
'GTT': 'V',
'GTC': 'V',
'GTA': 'V',
'GTG': 'V',
'TCT': 'S',
'TCC': 'S',
'TCA': 'S',
'TCG': 'S',
'CCT': 'P',
'CCC': 'P',
'CCA': 'P',
'CCG': 'P',
'ACT': 'T',
'ACC': 'T',
'ACA': 'T',
'ACG': 'T',
'GCT': 'A',
'GCC': 'A',
'GCA': 'A',
'GCG': 'A',
'TAT': 'Y',
'TAC': 'Y',
'TAA': '*',
'TAG': '*',
'CAT': 'H',
'CAC': 'H',
'CAA': 'Q',
'CAG': 'Q',
'AAT': 'N',
'AAC': 'N',
'AAA': 'K',
'AAG': 'K',
'GAT': 'D',
'GAC': 'D',
'GAA': 'E',
'GAG': 'E',
'TGT': 'C',
'TGC': 'C',
'TGA': '*',
'TGG': 'W',
'CGT': 'R',
'CGC': 'R',
'CGA': 'R',
'CGG': 'R',
'AGT': 'S',
'AGC': 'S',
'AGA': 'R',
'AGG': 'R',
'GGT': 'G',
'GGC': 'G',
'GGA': 'G',
'GGG': 'G',

}



# mapping onto antisense strand, includes N:N to map undetermined nucleotides onto undetermined nucleotides
antisenseMap = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A', 'N': 'N'}

def calculateRandomMatchPeptideProbabilities(fastaFileName, pepLengths=range(5,20), ILEqual=False):

    probDict = {}
    numAA = 20
    if ILEqual:
        numAA -= 1

    for pepLength in pepLengths:
        seqGen = sequenceGenerator(fastaFileName)
        uniquePeptSet = set()
        for seqName, sequence in seqGen:
            startInd = 0
            while startInd <= len(sequence) - pepLength:
                subSeq = sequence[startInd:startInd+pepLength]
                if all(aa in Constants.aminoacids for aa in subSeq):
                    uniquePeptSet.add(subSeq)
                startInd += 1

        probDict[pepLength] = '%.5g' % ((float(len(uniquePeptSet))/(numAA**pepLength)),)

    return probDict

# frames are 1,2,3,-1,-2,-3
def getTranslation(sequence, frame):
    if frame < 0:
        # get antisense strand and reverse orientation
        sequence = ''.join([antisenseMap[nucleotide] for nucleotide in sequence])[::-1]
        frame = frame * -1

    translatedSequence = []
    for i in range(frame-1, len(sequence), 3):
        try:
            translatedSequence += [codonDict[sequence[i:i+3]]]
        except KeyError:
            translatedSequence += ['X']

    return ''.join(translatedSequence)
    
    
def sequenceGenerator(fastaFileName):
    fastaFile = open(fastaFileName)
    sequence = [] # Added by smp to fix error: UnboundLocalError: local variable 'sequence' referenced before assignment at line 143
    seqName = None
    for line in fastaFile:
        if line[0] == '>':
            if seqName != None:
                yield seqName, ''.join(sequence)

            seqName = line.strip()
            sequence = []
        else:
            sequence += [line.strip()]

    yield seqName, ''.join(sequence)

# pattern is of the form (match1_length, gap_length, match2_length) 
# i.e., (4,2,4) matches a 4-mer, with a 2-aa gap (can be anything) followed by another 4-mer
# calculates the number and proportion (out of all possible combinations) of this pattern
def calculate_pattern_frequency(pattern, fasta_file, AAs = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']) ):
    seq_gen = sequenceGenerator(fasta_file)
    found_seqs = {}

    pattern_length = sum(pattern)
    second_match_offset = pattern[0] + pattern[1]
    for seq_name, sequence in seq_gen:
        for i in range(len(sequence) - pattern_length):
            pept1 = sequence[i:i+pattern[0]]
            pept2 = sequence[i+second_match_offset:i+pattern_length]
            if not ( any([aa not in AAs for aa in pept1]) or any([aa not in AAs for aa in pept2]) ):
                found_seqs[(pept1, pept2)] = 0

    return len(found_seqs)

# Calculates number of unique sequences present in database of length N
def calculate_peptlength_frequency(pept_length, fasta_file, AAs = set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']) ):
    seq_gen = sequenceGenerator(fasta_file)
    found_seqs = {}

    for seq_name, sequence in seq_gen:
        for i in range(len(sequence) - pept_length):
            peptide = sequence[i:i+pept_length]
            if not any([aa not in AAs for aa in peptide]):
                found_seqs[peptide] = 0

    return len(found_seqs)    


def calculate_db_stats(fasta_file):
    seq_gen = sequenceGenerator(fasta_file)

    num_prot = 0
    num_aa = 0
    aa_freq = defaultdict(int)

    for seq_name, sequence in seq_gen:
        num_prot += 1
        num_aa += len(sequence)
        for aa in sequence:
            aa_freq[aa] += 1

    return {'Proteins': num_prot, 'AAs': num_aa, 'AA Frequency': aa_freq}

def get_taxons_at_score_percent_cutoff(get_taxons_file, score_percent_cutoff = 0.001):
    taxons = []
    all_pepts = set()
    cols, data = DataFile.getScanInfo(get_taxons_file, delimiter='\t')

    for item in data:
        all_pepts |= eval(item['Peptide Cover'])

    for item in data:
        if float(item['Score']) / len(all_pepts) >= score_percent_cutoff:
            taxons += [item['Taxon']]

    return taxons

# calculates first order markov model based on sequence database
# aa_conditional_probs of form P(A_n | A_n-1) = aa_conditional_probs[A_n-1][A_n]
def create_first_order_markov_model(fasta_file, initial_count = 0.1):
    seq_gen = sequenceGenerator(fasta_file)

    num_prot = 0
    num_aa = 0
    aa_freq = {}
    aa_conditional_probs = {}

    for seq_name, sequence in seq_gen:
        num_prot += 1

        for i in range(len(sequence)):
            num_aa += 1

            if sequence[i] not in aa_freq:
                aa_freq[sequence[i]] = 1

                # Initialize counts in markov model to initial_count (to avoid probability of 0)
                for aa in aa_conditional_probs:
                    aa_conditional_probs[aa][sequence[i]] = initial_count

                aa_conditional_probs[sequence[i]] = {}
                for aa in aa_conditional_probs:
                    aa_conditional_probs[sequence[i]][aa] = initial_count
            else:
                aa_freq[sequence[i]] += 1
                if i > 0:
                    aa_conditional_probs[sequence[i-1]][sequence[i]] += 1

    # calculate probs
    total_aa = float(sum(aa_freq.values()))
    for aa in aa_freq:
        aa_freq[aa] = aa_freq[aa]/total_aa

    for aa in aa_conditional_probs:
        total_counts = float(sum(aa_conditional_probs[aa].values()))
        for aa2 in aa_conditional_probs[aa]:
            aa_conditional_probs[aa][aa2] = aa_conditional_probs[aa][aa2]/total_counts

    return {'Proteins': num_prot, 'AAs': num_aa, 'AA Probabilities': aa_freq, 'Markov Model': aa_conditional_probs}


# NOTE: Due to the way this is coded, insertions and deletions are still penalized more harshly than mods
# This is because the specificity is calculated off of the context, not the mod_context, and the context does not include the inserted amino acids (or the removal of the deleted ones)
def getProteaseContext(context, motifs, seq_delim='-'):

    N_satisfied, C_satisfied = False, False
    if context[0] == seq_delim:
        N_satisfied = True
    else:
        for motif in motifs:
            if bool(re.match(motif[0], context[0])) and bool(re.match(motif[1], context[2])):
                N_satisfied = True
                break

    if context[-1] == seq_delim:
        C_satisfied = True
    else:
        for motif in motifs:
            if bool(re.match(motif[0], context[-3])) and bool(re.match(motif[1], context[-1])):
                C_satisfied = True
                break

    if N_satisfied and C_satisfied:
        return 'full'
    elif not N_satisfied and C_satisfied:
        return 'n-ragged'
    elif N_satisfied and not C_satisfied:
        return 'c-ragged'
    else:
        return 'nonspecific'
                                                                                                                                                                            
def getNumMissedCleavages(context, mod_context, mods, motifs, modAASymbol = '-'):

    num_missed = 0
    for i in range(2, len(mod_context) - 3):
        missed_cleavage = False
        for motif in motifs:
            if re.match(motif[0], mod_context[i]) and re.match(motif[1], mod_context[i+1]):
                if not mod_context[i] == modAASymbol:
                    missed_cleavage = True
                else:
                    # If mod in N- or C-terminal, then it is still a missed cleavage
                    if i == 2 and 'N-term' in [mod[-1][0] for mod in mods]:
                        missed_cleavage = True
                    if i == len(mod_context) - 4 and 'C-term' in [mod[-1][0] for mod in mods]:
                        missed_cleavage = True

        if missed_cleavage:
            num_missed += 1

    # TODO: Subtract missed cleavages introduced by insertions (can show up for transpeptidation events, not really missed cleavages)?
    #print(context, mod_context, mods, num_missed)
    return num_missed



# calculates probability of sequence (relative to all other sequences of same length) using markov model
def calculate_sequence_probability_markov(markov_model, sequence):
    if len(sequence) == 0:
        return 1

    prob =  decimal.Decimal(str(markov_model['AA Probabilities'][sequence[0]]))
    for i in range(1, len(sequence)):
        prob = prob * decimal.Decimal(str(markov_model['Markov Model'][sequence[i-1]][sequence[i]]))

    return prob

def calculate_match_probability_markov(markov_model, sequence):
    if len(sequence) == 0:
        return 1

    seq_prob = calculate_sequence_probability_markov(markov_model, sequence)
    num_tries = markov_model['AAs'] - markov_model['Proteins'] * (len(sequence) - 1)
    return decimal.Decimal(1) - (decimal.Decimal(1) - seq_prob)**num_tries

def calculate_sequence_probability_constant(sequence, num_aas = 19):
    if len(sequence) == 0:
        return 1

    return decimal.Decimal(1)/num_aas**len(sequence)

def calculate_match_probability_constant(db_stats, sequence):
    if len(sequence) == 0:
        return 1

    seq_prob = calculate_sequence_probability_constant(sequence)
    num_tries = db_stats['AAs'] - db_stats['Proteins'] * len(sequence) - 1
    return decimal.Decimal(1) - (decimal.Decimal(1) - seq_prob)**num_tries

# Uses constant probability model
def calculate_sequence_length_probabilities(db_stats, lengths, default_prob = 1.0):
    probabilities = defaultdict(lambda: default_prob)

    for length in lengths:
        probabilities[length] = calculate_match_probability_constant(db_stats, ['A'] * length)

    return probabilities

# seqStart and seqEnd are defined as distance away from start of chromosome with first position indexed as one, seqStart is inclusive and seqEnd is exclusive
def chunkAndWriteSequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, lineLength=60, chunkSize=None, dontReadThrough=['*', 'X'], minPeptLength=10, overhang=100):
    i = 0
    while i < len(fullSequence):
        if fullSequence[i] in dontReadThrough:
            i += 1
        else:
            lowestInd = min(fullSequence.find(forbidAA,i) for forbidAA in dontReadThrough)
            if lowestInd == -1:
                lowestInd = len(fullSequence)
                for forbidAA in dontReadThrough:
                    if fullSequence.find(forbidAA,i) != -1:
                        lowestInd = min(lowestInd, fullSequence.find(forbidAA,i))

            if lowestInd - i > minPeptLength and (not chunkSize or lowestInd - i <= chunkSize):
                prepAndWriteSubsequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, i, lowestInd, lineLength=lineLength)
            elif chunkSize and lowestInd - i > chunkSize:
                for startInd in range(i, lowestInd, chunkSize-overhang):
                    prepAndWriteSubsequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, startInd, min(startInd+chunkSize, lowestInd), lineLength=60)

            i = lowestInd


def prepAndWriteSubsequence(outFile, baseSeqName, fullSequence, frame, ntSeqLength, startInd, endInd, lineLength=60):
    seqStart = 1 + abs(frame) + 3*startInd
    seqEnd = seqStart + 3*(endInd - startInd)
    if frame < 0:
        seqEnd = ntSeqLength - seqStart + 2
        seqStart = seqEnd - 3*(endInd - startInd)

    seqName = getTransSeqNameForGFY(baseSeqName, frame, '_start_Pos%i_end_Pos%i' % (seqStart, seqEnd))
    writeSequence(outFile, seqName, fullSequence[startInd:endInd], lineLength=lineLength)

def writeSequence(outFile, seqName, sequence, lineLength=60):
    outFile.write(seqName + '\n')
    for i in range(0, len(sequence), lineLength):
        outFile.write(sequence[i:i+lineLength] + '\n')

def getTransSeqNameForGFY(seqName, frame, addStringToBase=''):
    seqNameList = seqName.split(' ')
    seqNameList[0] = seqNameList[0] + ('_+' if frame > 0 else '_') + str(frame) + addStringToBase
    return ' '.join(seqNameList) + ' frame=%i' %(frame,)

# Generates decoy database from target (reversed + shuffle)
# 1:1 mapping of same subsequences (to avoid decoy factor exceeding 1:1)
# Keeps certain AAs fixed so that similar numbers of peptides are returned from target and decoy at a given precmass
def generateDecoyDB(seq_gen, keep_fixed =['K']):
    decoy_seqs = {}
    decoy_pept_mapping = {}

    for seq_name, sequence in seq_gen:
        fixed_locs = [-1]
        new_seq = []
        for i in range(len(sequence)):
            if sequence[i] in keep_fixed:
                fixed_locs += [i]

        fixed_locs += [len(sequence)]

        for i in range(len(fixed_locs) - 1):
            replace_pept = sequence[fixed_locs[i]+1:fixed_locs[i+1]]
            if replace_pept not in decoy_pept_mapping:
                r_list = list(replace_pept)
                random.shuffle(r_list)
                decoy_pept_mapping[replace_pept] = r_list

            new_seq += decoy_pept_mapping[replace_pept]

            if fixed_locs[i+1] < len(sequence):
                new_seq += [ sequence[fixed_locs[i+1]] ]

        decoy_seqs[seq_name] = ''.join(new_seq[::-1])

    return decoy_seqs

def writeTrueAndDecoyDB(source_fasta, outFileName, keep_fixed = ['K'], decoy_delim='##'):
    seq_gen = sequenceGenerator(source_fasta)
    decoy_db = generateDecoyDB(seq_gen, keep_fixed=keep_fixed)

    outFile = open(outFileName, 'w')

    seq_gen = sequenceGenerator(source_fasta)
    for seq_name, sequence in seq_gen:
        writeSequence(outFile, seq_name, sequence, lineLength=60)

    for seq_name, sequence in decoy_db.items():
        seq_name = seq_name[0] + decoy_delim + seq_name[1:]
        writeSequence(outFile, seq_name, sequence, lineLength=60)

    outFile.close()



# Define maximum size of index to be slightly less than UINT32_MAX to allow suffix array indexing
# Also saves offsets file defining sequence names offsets in db
def makeDBForFMIndexFromFASTA(fastaFile, outFileBase, transformLtoI = True, seqSeperator='-', maxIndexSize = 2**31-2):
    print(outFileBase)
    print("\n")
    print(fastaFile)
    print("\n")
    seqGen = sequenceGenerator(fastaFile)

    i = 1
    indSize = 0

    # Initialize first index partion
    nameDB = dbm.open(outFileBase + '.seqnames.%i'%i, 'n')
    offSet = 0
    outFile = open(outFileBase + '_fmFormatted.txt.%i'%i, 'w')
    offsets = []

    offsets_arr = []
    for seqName, sequence in seqGen:
        # Sometimes FASTAs have weird whitespace in them, can screw up index
        sequence = clean_sequence(sequence)

        indSize += len(sequence) + 1
        if indSize > maxIndexSize:
            print('Index size %i greater than max %i, creating new index partition'%(indSize, maxIndexSize))
            nameDB.close()
            outFile.close()
            offsets_arr += [np.array(offsets, dtype=np.uint32)]

            i += 1
            indSize = len(sequence) + 1
            nameDB = dbm.open(outFileBase + '.seqnames.%i'%i, 'n')
            offSet = 0
            outFile = open(outFileBase + '_fmFormatted.txt.%i'%i, 'w')
            offsets = []


        parsed_seq_name = seqName[1:]
        nameDB[str(offSet)] = parsed_seq_name
        offsets += [offSet]
        offSet += len(sequence) + 1

        outFile.write(seqSeperator)
        if transformLtoI:
            outFile.write(sequence.replace('L', 'I'))
        else:
            outFile.write(sequence)

    offsets_arr += [np.array(offsets, dtype=np.uint32)]
    with open(outFileBase + '.offsets', 'w') as fout:
        pickle.dump(offsets_arr, fout)

    outFile.close()
    nameDB.close()

def clean_sequence(sequence):
    return sequence.replace(" ", "").replace("*", "")

def format_seqs_for_motifX(seqs_file, out_file_name, motif_len = 15, seq_len = 41, seq_end_aa = 'X', non_canon_aas = ['X', 'B', 'Z', 'U', 'O', 'J'], unique = False, max_size = 10000, max_num_crossvals = 1):
    seqs = open(seqs_file)
    new_seqs = []
    for seq in seqs:
        clip = (seq_len - motif_len)/2
        seq = seq.strip().upper()[clip:-clip]

        if any([aa in seq for aa in non_canon_aas]):
            continue

        # replace sequence ends with wildcard
        new_seqs += [seq.replace('*', 'X')]

    print('get unique')
    if unique:
        new_seqs = set(new_seqs)

    if len(new_seqs) > max_size:
        print('newseqs g')
        #for i in range(min(max_num_crossvals, len(new_seqs)/float(max_size) * 2)):
        rand_smpl = [ new_seqs[i] for i in sorted(random.sample(range(len(new_seqs)), max_size)) ]
        #out_file = open(out_file_name + str(i+1), 'w')
        out_file = open(out_file_name, 'w')
        for seq in rand_smpl:
            out_file.write(seq + '\n')
        out_file.close()
    else:
        print('maxsize g')
        out_file = open(out_file_name, 'w')
        for seq in new_seqs:
            out_file.write(seq + '\n')
        out_file.close()

if __name__ == '__main__':
    print('This program will take a FASTA file from the --lads argument and output a six-frame translation of the file to output. Number refers to maximum size of sequence in resulting FASTA file. If a chromosomal region exceeds this length with no stop codons, the sequence will be chunked with a 100 aa overhang at each edge. Minimum Length of peptide in FASTA file is 5.')
    options = ArgLib.parse(['lads', 'output', 'number'])

    outFile = open(options.output, 'w')
    #chunkSize = int(options.number)

    for seqName, sequence in sequenceGenerator(options.lads):
        #for frame in [1, 2, 3, -1, -2, -3]:
        for frame in [1,2,3]:
            #print(seqName, frame)
            transSeq = getTranslation(sequence.upper(), frame)
            #chunkAndWriteSequence(outFile, seqName, transSeq, frame, len(sequence), lineLength=60, chunkSize=2000, dontReadThrough=['X'], minPeptLength=10, overhang=50)
            #transSeqName = seqName + ('_+' if frame > 0 else '_') + str(frame)
            transSeqName = getTransSeqNameForGFY(seqName, frame)
            writeSequence(outFile, transSeqName, transSeq)
            
    outFile.close()
            

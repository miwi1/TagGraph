from collections import defaultdict
import itertools

import numpy as np
import copy
import math
import re

import StringFoldingWrapper as SFW

import Constants
import Database

aaBrevs = set([aaInfo[0].lower() for aaInfo in Constants.aminoacids.values()])
aaBrevs.add('xle')

# Checks to see if mod is an amino acid substitution, we won't use these to expand plausible mod candidates
def isNotAASub(mod_name):
    possAAs = mod_name.split('->')
    return not ( len(possAAs) == 2 and all([possAA in aaBrevs for possAA in possAAs]) )
        

# Retrieves peptides which match the following criteria:
#
# Score bonus is added to modified peptides matched to LysC-consisten site
#
#
def getConfidentPeptides(indexed_info, mod_freqs, score_key="Composite Score",  score_bonus=20):
    confident_items = {}

    for scanF in indexed_info:
        unmod = False
        best_mod_count = 0
        confident_items[scanF] = []
        for item in indexed_info[scanF]:
            mods = [mod for mod in eval(item['Match Modifications']) if mod[0][0] != 'Isobaric Substitution']
            context = getLysCContext(item['Mod Context'])
            item[score_key] = float(item[score_key])

            if len(mods) == 0:
                if context == 'full':
                    item[score_key] += score_bonus

                # Reset matched scans if mod peptide hasn't been found yet
                if not unmod:
                    unmod = True
                    confident_items[scanF] = []

                confident_items[scanF] += [item]                    
            elif not unmod:
                mod_tuple = getUniqueModTuple(mods)
                if context == 'full' and mod_tuple in mod_freqs:
                    item[score_key] += score_bonus

                current_mod_count = 0
                if mod_tuple in mod_freqs:
                    current_mod_count = mod_freqs[mod_tuple]
                    
                if current_mod_count > best_mod_count:
                    best_mod_count = current_mod_count
                    confident_items[scanF] = []

                if current_mod_count == best_mod_count:
                    confident_items[scanF] += [item]

        confident_items[scanF] = sorted(confident_items[scanF], key = lambda k: -k[score_key])

    return confident_items

def redundify(mods, redundant_mod_map):

    new_mod_tuple = []
    trimmed_mods = [mod for mod in mods if mod[0][0] != 'Isobaric Substitution']
    for mod in getUniqueModTuple(trimmed_mods):
        if (mod,) in redundant_mod_map:
            new_mod_tuple += [ redundant_mod_map[(mod,)][0] ]
        else:
            new_mod_tuple += [mod]

    return tuple(sorted(new_mod_tuple))

# TODO: mod freqs is set to only find LysC peptides, switch to general protease later
# Calculates mod frequencies from top scoring peptide in 
def getModFreqs(indexed_tag_graph_info, score_key="Composite Score", protease_specific = True):
    mod_freqs = defaultdict(int)
    
    for scanF in indexed_tag_graph_info:
        top_score = float(indexed_tag_graph_info[scanF][0][score_key])
        for item in indexed_tag_graph_info[scanF]:
            if float(item[score_key]) != top_score:
                break
            
            mods = [mod for mod in eval(item['Match Modifications']) if mod[0][0] != 'Isobaric Substitution']

            # Check to make sure mods are present and are not all isobaric mods
            if len(mods) > 0:
                if not protease_specific or getLysCContext(item['Mod Context']) == 'full':
                    mod_tuple = getUniqueModTuple(mods)
                    mod_freqs[mod_tuple] += 1

    return mod_freqs

# Mod classes are Defined Mod: 1, AA Sub: 2, Undefined Mod: 3, Insertion or Deletion: 4 (0 for no mods, but the distribution does not consider this cases
# For multiply modified peptides, mod_class is max of all individual mod classes
def getModClass(mod_tuple):
    mod_classes = []
    for mod in mod_tuple:
        if mod[0] == 'insertion' or mod[0] == 'deletion':
            mod_classes += [4]
        elif mod[0] == 'undefined mass shift':
            mod_classes += [3]
        elif isNotAASub(mod[0]):
            mod_classes += [1]
        else:
            mod_classes += [2]

    try:
        return max(mod_classes)
    except ValueError:
        return 0
    
# Round undefined mass shifts to first place and don't include AA specificity for easier matching to one another
# TODO: ADD N-term and C-term specificity in to undefined mass shift (should rescue sensitivity and allow aa-calling)
def getUniqueModTuple(mods, undef_mod_round_precision = 1):
    mod_tuple = []
    for mod in mods:
        if mod[0][0] == 'Undefined Mass Shift':
            # TODO: Just lowercase everything in Unimod so that time isn't wasted here
            mod_tuple += [(mod[0][0].lower(), round(mod[0][1], undef_mod_round_precision), ())]
        else:
            mod_tuple += [(mod[0][0].lower(), round(mod[0][1], 2), mod[1])]
    return tuple(sorted(mod_tuple))


# Precalculates attributes on which EM will optimize
# Currently calculates: mod tuples, protease context (only LysC)
#                                         0    1              2        3            4          5          6                   7            8               9         10                11          12             13            14                  15                                          16
# key is scanF, value is list of format [(id, spectrum_score, context, mod_context, mod_tuple, proteins, matching_tag_length, specificity, ppm mass error, num_mods, missed_cleavages, mod_class, context_length, max_mod_size, modified (0 or 1)), set(mod_masses rounded to nearest integer)), charge]
#    0                     1                    2                    3
#  (unique_sibling_var, unique_mods_on_context, num_mod_occurrences, num_single_mods_found_in_context (for multi-modded peptides)),
#   0                      1
# (em_probability, 1-lg10 em probability)]
def precalculateEMAttributes(indexed_taggraph_results, motifs, ep_step = 0.0025, mod_size_bin = 100, max_mod_size_bin = 10):
    for scanF in indexed_taggraph_results:
        for i, item in enumerate(indexed_taggraph_results[scanF]):

            indexed_taggraph_results[scanF][i] = precalculateEMAttributesForItem(item, motifs, ep_step = ep_step, mod_size_bin = mod_size_bin, max_mod_size_bin = max_mod_size_bin)

def precalculateEMAttributesForItem(item, motifs, ep_step = 0.0025, mod_size_bin = 250, max_mod_size_bin = 8, ppm_bin=0.25, max_ppm = 15, min_ppm = -15):
    all_mods = eval(item[5])
    mods = [mod for mod in all_mods if mod[0][0] != 'Isobaric Substitution']
    
    mod_tuple = getUniqueModTuple(mods)
    num_mods = len(mod_tuple)
    
    """
    # TODO: Currently DOES NOT include mass error of isobaric substitutions in distribution, try with and without
    # /DISREGARD Perhaps dangerous because all existing mods are already binned by ep_step, could underestimate prob for some mods with low but nonzero mass errors
    # Current layout:
    # Isobaric substitutions and all defined mods get put into a mass error bin (as defined by ep_step)
    # Items with no mod get put into their own special bin ('No Mod')
    # Items with Insertion, Deletion and Undefined Mass Shift get put into their own bin ('Indeterminate')
    # TODO: Use absolute value of mass error instead of mass error?
    if len(all_mods) == 0:
        mass_error = 'Exact'
    else:
        mass_error = 'Indeterminate'
        for mod in all_mods:
            if mod[0][0] not in ['Insertion', 'Deletion', 'Undefined Mass Shift']:
                if mass_error == 'Indeterminate':
                    mass_error = round(mod[0][-1]/ep_step)
                elif abs(round(mod[0][-1]/ep_step)) > abs(mass_error):
                    mass_error = round(mod[0][-1]/ep_step)
    """

    ppm_error = round(item[2]/ppm_bin) * ppm_bin
    ppm_error = max( min(ppm_error, max_ppm), min_ppm )
    
    proteins = SFW.fold_list(eval(item[6]))
            
    specificity = Database.getProteaseContext(item[3], motifs)
    missed_cleavages = Database.getNumMissedCleavages(item[3], item[4], item[5], motifs)
    mod_class = getModClass(mod_tuple)
    modified = min( len(mod_tuple), 1 )

    if len(mod_tuple) == 0:
        max_mod_size = 'Indeterminate'
    else:
        max_mod_size = min(max_mod_size_bin, int(max([abs(mod[1]) for mod in mod_tuple]))/mod_size_bin)

    # Charge State-specific scoring models only trained for charges 2, 3, and 4
    charge = max(min(item[8], 4), 2)

    return [( item[0], item[1], item[3], item[4], mod_tuple, proteins, item[7], specificity, ppm_error, num_mods, missed_cleavages, mod_class, len(item[3]) - 4, max_mod_size, modified, tuple(sorted(round(mod[1]) for mod in mod_tuple)), charge ), tuple(), tuple()]


def recalculateNetworkAttributes(indexed_taggraph_results, top_only = True, first_iteration = False):

    # For each protein, only adding top ranked context for a particular item to avoid overly inflating protein counts
    protein_map = defaultdict(set)
    
    # Add mod tuples here (for seeing whether or not single mods found for multi-modded pepts)
    context_map = defaultdict(set)
    
    # Add mod masses here (to avoid inflating context count)
    context_mass_map = defaultdict(set)
    mod_map = defaultdict(set)

    for scanF in indexed_taggraph_results:
        protein_context_map = defaultdict(list)
        top_item = indexed_taggraph_results[scanF][0]

        for protein in top_item[0][5]:
            # protein,            context
            protein_context_map[protein] += [(0, top_item[0][2])]

        # context,              mod_tuple
        context_map[top_item[0][2]].add(top_item[0][4])
        context_mass_map[top_item[0][2]].add(top_item[0][15])
        # mod_tuple            context
        mod_map[top_item[0][4]].add(top_item[0][2])

        # Do not calculate network attributes based on rest of results if top_only is true
        if not top_only:
            top_em_prob = top_item[2][1]
            for i, item in enumerate(indexed_taggraph_results[scanF][1:]):
                for protein in item[0][5]:
                    # protein,            context
                    protein_context_map[protein] += [(i+1, item[0][2])]
                
                # context,              mod_tuple
                context_mass_map[item[0][2]].add(item[0][15])

                # Don't add mod counts (or exact mod identity to context_map) for lower probability items even if top_only is set to false (or if em probability is 0)
                item_em_prob = item[2][1]
                #if item[0][2] == top_item[0][2]:
                #    print scanF, item[0][3], item[0][4], top_item[0][3], top_item[0][4], item_em_prob, top_em_prob, item_em_prob == top_em_prob
                if item[0][2] == top_item[0][2] and item_em_prob == top_em_prob:
                    mod_map[item[0][4]].add(item[0][2])
                    context_map[item[0][2]].add(item[0][4])



                
        # Add top ranked context for each protein to protein_map (avoids overcounting bay adding multiple similar contexts to same protein for a given scan)
        for protein in protein_context_map:
            protein_map[protein].add(sorted(protein_context_map[protein])[0][1])
                                
    protein_count = defaultdict(lambda: 1)
    context_count = defaultdict(lambda: 1)
    mod_count = defaultdict(lambda: 1)

    for context in context_mass_map:
        context_count[context] = len(context_mass_map[context])

    for protein in protein_map:
        protein_count[protein] = len(protein_map[protein])

    for mod in mod_map:
        mod_count[mod] = len(mod_map[mod])

    # Write network level attributes
    for scanF in indexed_taggraph_results:
        for item in indexed_taggraph_results[scanF]:
            # TODO: Split context mod variants feature into 'number of mods which contain this mod' and 'number of mods which do not contain this mod' (make sure these are independent), used to give additional evidence for single mods which have multiple modified peptides with the same single mod at the same context
            #          Unique sibling variants                                  context mod variants       num mod occurrences    # num single mods found in context (for mult-modded peptides)
            if len(item[0][4]) == 1:
                num_single = 1
            else:
                num_single = len([1 for mod in item[0][4] if (mod,) in context_map[item[0][2]]])
            item[1] = (max([protein_count[protein] for protein in item[0][5]]), context_count[item[0][2]], mod_count[item[0][4]], num_single)
            
                      
        

# Sets initial EM Probabilities of taggraph results based on initial models for positive and negative gaussians
def initializeEMProbabilities(indexed_taggraph_results, spectrum_match_models, db_match_models, prior_probs = (0.5, 0.5)):
    for scanF in indexed_taggraph_results:
        for item in indexed_taggraph_results[scanF]:
            # print scanF, item
            pos_model_prob = spectrum_match_models[0](item[0][1], item[0][12]) * db_match_models[0][item[0][6]]
            neg_model_prob = spectrum_match_models[1](item[0][1], item[0][12]) * db_match_models[1][item[0][6]]
            em_prob = pos_model_prob*prior_probs[0] / (pos_model_prob * prior_probs[0] + neg_model_prob * prior_probs[1])
            lg_em_prob = calculateLogEMProb(em_prob, neg_model_prob, pos_model_prob, prior_probs)
            item[2] = (em_prob, lg_em_prob)

# Calculates P(Context | +) and P(Context | -) for indexed taggraph_results
# Calculated using MLE for parzen windows (using soft-EM assumption)
def calculateProteaseContextProbabilities(indexed_taggraph_results, num_scans, prior_probs, min_prob=1e-40):
    pos_probs = defaultdict(lambda: min_prob)
    neg_probs = defaultdict(lambda: min_prob)

    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]

        pos_probs[item[0][7]] += item[2][0]
        neg_probs[item[0][7]] += 1 - item[2][0]

    for specificity in pos_probs:
        pos_probs[specificity] = max( min_prob, pos_probs[specificity] / (num_scans * prior_probs[0]) )

    for specificity in neg_probs:
        neg_probs[specificity] = max( min_prob, neg_probs[specificity] / (num_scans * prior_probs[1]) )

    return pos_probs, neg_probs

def calculateMissedCleavageProbabilities(indexed_taggraph_results, num_scans, prior_probs, min_prob=1e-40):
    pos_probs = defaultdict(lambda: min_prob)
    neg_probs = defaultdict(lambda: min_prob)

    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        
        pos_probs[item[0][10]] += item[2][0]
        neg_probs[item[0][10]] += 1 - item[2][0]
        
    for num_missed in pos_probs:
        pos_probs[num_missed] = max( min_prob, pos_probs[num_missed] / (num_scans * prior_probs[0]) )
        
    for num_missed in neg_probs:
        neg_probs[num_missed] = max( min_prob, neg_probs[num_missed] / (num_scans * prior_probs[1]) )
        
    return pos_probs, neg_probs
                                                            
def calculatePPMErrorProbabilities(indexed_taggraph_results, num_scans, prior_probs, min_prob=1e-40):
    pos_probs = defaultdict(lambda: min_prob)
    neg_probs = defaultdict(lambda: min_prob)

    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        
        pos_probs[item[0][8]] += item[2][0]
        neg_probs[item[0][8]] += 1 - item[2][0]
        
    for ppm in pos_probs:
        pos_probs[ppm] = max( min_prob, pos_probs[ppm] / (num_scans * prior_probs[0]) )
        
    for ppm in neg_probs:
        neg_probs[ppm] = max( min_prob, neg_probs[ppm] / (num_scans * prior_probs[1]) )
        
    return pos_probs, neg_probs

# Calculates P(Mod | +) and P(Mod | -) for indexed taggraph results
# Zero probabilities shouldn't theoretically arise unless P(+|D) or P(-|D) == 1 for a particular mod
# Calculated using MLE for parzen windows (using soft-EM assumption)
# pos_default and neg_default set default bias against mods not ranked in top 1 (and thus not present in dict). Currently: 0.01
# pos_default and neg_default are also used as the minimum values any entry in the two model databases can take (to avoid entries with 0 probability showing up)
def calculateModProbabilities(indexed_taggraph_results, num_scans, prior_probs):

    pos_modified, neg_modified, pos_num_mods, neg_num_mods, pos_context_probs, neg_context_probs, pos_mod_size_probs, neg_mod_size_probs, pos_occurrence_probs, neg_occurrence_probs, pos_single_mod_of_multi, neg_single_mod_of_multi, pos_class_probs, neg_class_probs = getModProbDistributions(indexed_taggraph_results, num_scans, prior_probs)
    
    def pos_mod_prob(item):
        return pos_modified[item[0][14]] * pos_num_mods[item[0][14]][item[0][9]] * pos_context_probs[item[1][1]] * pos_mod_size_probs[item[0][14]][item[0][13]] * pos_occurrence_probs[item[0][9]][item[1][2]] * pos_single_mod_of_multi[item[0][9]][item[1][3]] * pos_class_probs[item[0][14]][item[0][11]]
    

    def neg_mod_prob(item):
        return neg_modified[item[0][14]] * neg_num_mods[item[0][14]][item[0][9]] * neg_context_probs[item[1][1]] * neg_mod_size_probs[item[0][14]][item[0][13]] * neg_occurrence_probs[item[0][9]][item[1][2]] * neg_single_mod_of_multi[item[0][9]][item[1][3]] * neg_class_probs[item[0][14]][item[0][11]]

    return pos_mod_prob, neg_mod_prob

def getModProbDistributions(indexed_taggraph_results, num_scans, prior_probs, min_prob = 1e-10):

    pos_context_probs = defaultdict(lambda: min_prob)
    neg_context_probs = defaultdict(lambda: min_prob)

    pos_num_mods = defaultdict(lambda: defaultdict(lambda: min_prob))
    neg_num_mods = defaultdict(lambda: defaultdict(lambda: min_prob))

    pos_single_mod_of_multi = defaultdict(lambda: defaultdict(lambda: min_prob))
    neg_single_mod_of_multi = defaultdict(lambda: defaultdict(lambda: min_prob))

    pos_occurrence_probs = defaultdict(lambda: defaultdict(lambda: min_prob))
    neg_occurrence_probs = defaultdict(lambda: defaultdict(lambda: min_prob))

    pos_class_probs = defaultdict(lambda: defaultdict(lambda: min_prob))
    neg_class_probs = defaultdict(lambda: defaultdict(lambda: min_prob))

    pos_mod_size_probs = defaultdict(lambda: defaultdict(lambda: min_prob))
    neg_mod_size_probs = defaultdict(lambda: defaultdict(lambda: min_prob))

    # This is just to save this sum so it can be used later in calculating model probabilities
    pos_num_mods_sum = defaultdict(lambda: min_prob)
    neg_num_mods_sum = defaultdict(lambda: min_prob)

    pos_modified = defaultdict(float)
    neg_modified = defaultdict(float)

    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]

        pos_context_probs[item[1][1]] += item[2][0]
        neg_context_probs[item[1][1]] += 1 - item[2][0]

        pos_num_mods[item[0][14]][item[0][9]] += item[2][0]
        neg_num_mods[item[0][14]][item[0][9]] += 1 - item[2][0]

        pos_mod_size_probs[item[0][14]][item[0][13]] += item[2][0]
        neg_mod_size_probs[item[0][14]][item[0][13]] += 1 - item[2][0]

        pos_occurrence_probs[item[0][9]][item[1][2]] += item[2][0]
        neg_occurrence_probs[item[0][9]][item[1][2]] += 1 - item[2][0]

        pos_class_probs[item[0][14]][item[0][11]] += item[2][0]
        neg_class_probs[item[0][14]][item[0][11]] += 1 - item[2][0]

        pos_single_mod_of_multi[item[0][9]][item[1][3]] += item[2][0]
        neg_single_mod_of_multi[item[0][9]][item[1][3]] += 1 - item[2][0]

        pos_num_mods_sum[item[0][9]] += item[2][0]
        neg_num_mods_sum[item[0][9]] += 1 - item[2][0]

        pos_modified[item[0][14]] += item[2][0]
        neg_modified[item[0][14]] += 1 - item[2][0]

    for context_count in pos_context_probs:
        pos_context_probs[context_count] = max(min_prob, pos_context_probs[context_count]/(num_scans * prior_probs[0]))
        neg_context_probs[context_count] = max(min_prob, neg_context_probs[context_count]/(num_scans * prior_probs[1]))

    for modified in pos_num_mods:
        for num_mods in pos_num_mods[modified]:
            pos_num_mods[modified][num_mods] = max(min_prob, pos_num_mods[modified][num_mods]/pos_modified[modified])
            neg_num_mods[modified][num_mods] = max(min_prob, neg_num_mods[modified][num_mods]/neg_modified[modified])

    for modified in pos_mod_size_probs:
        for mod_size in pos_mod_size_probs[modified]:
            pos_mod_size_probs[modified][mod_size] = max(min_prob, pos_mod_size_probs[modified][mod_size]/pos_modified[modified])
            neg_mod_size_probs[modified][mod_size] = max(min_prob, neg_mod_size_probs[modified][mod_size]/neg_modified[modified])

    for num_mods in pos_occurrence_probs:
        for num_occ in pos_occurrence_probs[num_mods]:
            pos_occurrence_probs[num_mods][num_occ] = max(min_prob, pos_occurrence_probs[num_mods][num_occ]/pos_num_mods_sum[num_mods])
            neg_occurrence_probs[num_mods][num_occ] = max(min_prob, neg_occurrence_probs[num_mods][num_occ]/neg_num_mods_sum[num_mods])

    for num_mods in pos_single_mod_of_multi:
        for num_single in pos_single_mod_of_multi[num_mods]:
            pos_single_mod_of_multi[num_mods][num_single] = max(min_prob, pos_single_mod_of_multi[num_mods][num_single]/pos_num_mods_sum[num_mods])
            neg_single_mod_of_multi[num_mods][num_single] = max(min_prob, neg_single_mod_of_multi[num_mods][num_single]/neg_num_mods_sum[num_mods])
            
    for modified in pos_class_probs:
        for mod_class in pos_class_probs[modified]:
            pos_class_probs[modified][mod_class] = max(min_prob, pos_class_probs[modified][mod_class]/pos_modified[modified])
            neg_class_probs[modified][mod_class] = max(min_prob, neg_class_probs[modified][mod_class]/neg_modified[modified])

    for modified in pos_modified:
        pos_modified[modified] = max(min_prob, pos_modified[modified]/ (num_scans * prior_probs[0]))
        neg_modified[modified] = max(min_prob, neg_modified[modified]/ (num_scans * prior_probs[1]))

    return pos_modified, neg_modified, pos_num_mods, neg_num_mods, pos_context_probs, neg_context_probs, pos_mod_size_probs, neg_mod_size_probs, pos_occurrence_probs, neg_occurrence_probs, pos_single_mod_of_multi, neg_single_mod_of_multi, pos_class_probs, neg_class_probs

# min_pos_count and min_neg_count are used to make sure that zero values don't appear in match probs
def calculateDBMatchProbabilities(indexed_taggraph_results, num_scans, prior_probs, min_prob=1e-40):
    pos_probs = defaultdict(float)
    neg_probs = defaultdict(float)

    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]

        pos_probs[item[0][6]] += item[2][0]
        neg_probs[item[0][6]] += 1 - item[2][0]

    for specificity in pos_probs:
        pos_probs[specificity] = max(min_prob, pos_probs[specificity] / (num_scans * prior_probs[0]))
        neg_probs[specificity] = max(min_prob, neg_probs[specificity] / (num_scans * prior_probs[1]))

    return pos_probs, neg_probs

def calculateUniqueSiblingCountProbabilities(indexed_taggraph_results, num_scans, prior_probs, min_prob=1e-5):
    pos_probs = defaultdict(lambda: min_prob)
    neg_probs = defaultdict(lambda: min_prob)

    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]

        pos_probs[item[1][0]] += item[2][0]
        neg_probs[item[1][0]] += 1 - item[2][0]

    for count in pos_probs:
        pos_probs[count] = max(min_prob, pos_probs[count] / (num_scans * prior_probs[0]))
        neg_probs[count] = max(min_prob, neg_probs[count] / (num_scans * prior_probs[1]))

    return pos_probs, neg_probs


def getPriorProbabilities(indexed_taggraph_results, num_scans):
    pos_prob_total = 0
    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        pos_prob_total += item[2][0]

    print 'Prior', (pos_prob_total/num_scans, 1 - pos_prob_total/num_scans)
    return (pos_prob_total/num_scans, 1 - pos_prob_total/num_scans)

# item[0][1], item[0][12]
# returns P(z|+), P(z|-) as well as multivariate gaussian parameters for P(S,L|z,+), P(S,L|z,-) [parameters returned for a gaussian for each charge state]
def getSpectrumMatchParametersMulti(indexed_taggraph_results, num_scans, prior_probs):

    pos_spec_score_mean, pos_pept_length_mean, neg_spec_score_mean, neg_pept_length_mean = defaultdict(float), defaultdict(float), defaultdict(float), defaultdict(float)
    charge_pos_sums, charge_neg_sums = defaultdict(float), defaultdict(float)
    
    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        #print item[0][14], item[2][0]

        pos_spec_score_mean[item[0][16]] += (item[0][1] * item[2][0])
        pos_pept_length_mean[item[0][16]] += (item[0][12] * item[2][0])
        neg_spec_score_mean[item[0][16]] += (item[0][1] * (1 - item[2][0]))
        neg_pept_length_mean[item[0][16]] += (item[0][12] * (1 - item[2][0]))
        #print pos_mean, neg_mean, item[0][14] * item[2][0], item[0][14] * (1 - item[2][0])

        # Sum up probabilities for charge state distributons
        charge_pos_sums[item[0][16]] += item[2][0]
        charge_neg_sums[item[0][16]] += 1 - item[2][0]


    #print pos_mean, neg_mean
    for charge in charge_pos_sums:
        pos_spec_score_mean[charge] = pos_spec_score_mean[charge] / charge_pos_sums[charge]
        pos_pept_length_mean[charge] = pos_pept_length_mean[charge] / charge_pos_sums[charge]
        neg_spec_score_mean[charge] = neg_spec_score_mean[charge] / charge_neg_sums[charge]
        neg_pept_length_mean[charge] = neg_pept_length_mean[charge] / charge_neg_sums[charge]

    # COORDINATES (ROW, COL)(1,1)(1,2)(2,1)(2,2)
    pos_sigma, neg_sigma = {}, {}
    for charge in charge_pos_sums:
        pos_sigma[charge] = [0.0, 0.0, 0.0, 0.0]
        neg_sigma[charge] = [0.0, 0.0, 0.0, 0.0]
        
    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        spec_score_mu = item[0][1] - pos_spec_score_mean[item[0][16]]
        pept_length_mu = item[0][12] - pos_pept_length_mean[item[0][16]]
        covar = item[2][0] * spec_score_mu * pept_length_mu
        pos_sigma[item[0][16]][0] += item[2][0] * spec_score_mu**2
        pos_sigma[item[0][16]][1] += covar
        pos_sigma[item[0][16]][2] += covar
        pos_sigma[item[0][16]][3] += item[2][0] * pept_length_mu**2
        
        
        spec_score_mu = item[0][1] - neg_spec_score_mean[item[0][16]]
        pept_length_mu = item[0][12] - neg_pept_length_mean[item[0][16]]
        covar = (1 - item[2][0]) * spec_score_mu * pept_length_mu
        neg_sigma[item[0][16]][0] += (1 - item[2][0]) * spec_score_mu**2
        neg_sigma[item[0][16]][1] += covar
        neg_sigma[item[0][16]][2] += covar
        neg_sigma[item[0][16]][3] += (1 - item[2][0]) * pept_length_mu**2
        
    #print pos_sigma, neg_sigma
    for charge in pos_sigma:
        for i in range(4):
            pos_sigma[charge][i] = pos_sigma[charge][i] / charge_pos_sums[charge]
            neg_sigma[charge][i] = neg_sigma[charge][i] / charge_neg_sums[charge]

    for charge in charge_pos_sums:
        charge_pos_sums[charge] = charge_pos_sums[charge] / (num_scans * prior_probs[0])
        charge_neg_sums[charge] = charge_neg_sums[charge] / (num_scans * prior_probs[1])

    return charge_pos_sums, [pos_spec_score_mean, pos_pept_length_mean], pos_sigma, charge_neg_sums, [neg_spec_score_mean, neg_pept_length_mean], neg_sigma
     
def getSpectrumMatchParameters(indexed_taggraph_results, num_scans, prior_probs):
    pos_mean, neg_mean = 0, 0
    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        pos_mean += item[0][1] * item[2][0]
        neg_mean += item[0][1] * (1 - item[2][0])


    pos_mean = pos_mean/ (prior_probs[0] * num_scans)
    neg_mean = neg_mean/ (prior_probs[1] * num_scans)
    
    pos_sigma, neg_sigma = 0, 0
    for scanF in indexed_taggraph_results:
        item =  indexed_taggraph_results[scanF][0]
        x_mu = item[0][1] - pos_mean
        pos_sigma += item[2][0] * x_mu**2

        x_mu = item[0][1] - neg_mean
        neg_sigma += (1 - item[2][0]) * x_mu**2

    pos_sigma = pos_sigma / (prior_probs[0] * num_scans)
    neg_sigma = neg_sigma / (prior_probs[1] * num_scans)

    return pos_mean, pos_sigma, neg_mean, neg_sigma

def calculateSpectrumMatchProbabilities(indexed_taggraph_results, num_scans, prior_probs):

    pos_charge, pos_mean, pos_sigma, neg_charge, neg_mean, neg_sigma = getSpectrumMatchParametersMulti(indexed_taggraph_results, num_scans, prior_probs)
    print 'charge pos', pos_charge, 'neg', neg_charge, 'gaussian', 'pos', pos_mean, pos_sigma, 'neg', neg_mean, neg_sigma

    pos_multi_pdf, neg_multi_pdf = {}, {}
    for charge in pos_charge:
        pos_multi_pdf[charge] = getMultivariateNormalPDF( (pos_mean[0][charge], pos_mean[1][charge]), pos_sigma[charge])
        neg_multi_pdf[charge] = getMultivariateNormalPDF( (neg_mean[0][charge], neg_mean[1][charge]), neg_sigma[charge])

    def pos_spec_prob(item):
        return pos_charge[item[0][16]] * pos_multi_pdf[item[0][16]](item[0][1], item[0][12])

    def neg_spec_prob(item):
        return neg_charge[item[0][16]] * neg_multi_pdf[item[0][16]](item[0][1], item[0][12])
    
    return pos_spec_prob, neg_spec_prob

def calculateLogEMProb(em_prob, neg_model_prob, pos_model_prob, prior_probs):
    if em_prob == 1.0:
        lg_em_prob = -math.log10( (neg_model_prob*prior_probs[1]) / (pos_model_prob*prior_probs[0]) )
    elif em_prob == 0.0:
        lg_em_prob = (pos_model_prob*prior_probs[0]) / (neg_model_prob*prior_probs[1])
    else:
        lg_em_prob = -math.log10( 1 - em_prob )

    return lg_em_prob
    
    
# Calculate EM Probabilities P(+|S,D,M) [i.e., hidden variables] based on inferred model parameters from previous step
# inferred model parameters: P(context|+), P(context|-), P(mod|+), P(mod|-), mixture gaussian parameters for spectrum probability score and database probability score
# Also re-sorts the returned values for each scanF by their EM Probabilities
def calculateEMProbabilities(indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs, rerank = True, score_all = True):

    for scanF in indexed_taggraph_results:
        # If we're not reranking, only update the EM Probability for the top ranked result
        if score_all:
            items = indexed_taggraph_results[scanF]
        else:
            items = [indexed_taggraph_results[scanF][0]]

        for item in items:
            # Note: For items which have mods that weren't present in the top mod slot, these values default to the chosen pos_bias and neg_bias (0 and 1 currently)
            #print item
            pos_model_prob, neg_model_prob = calculateClassProbabilities(item, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models)

            em_prob = pos_model_prob*prior_probs[0] / (neg_model_prob*prior_probs[1] + pos_model_prob*prior_probs[0])
            lg_em_prob = calculateLogEMProb(em_prob, neg_model_prob, pos_model_prob, prior_probs)

            item[2] = (em_prob, lg_em_prob)
        if rerank:
            indexed_taggraph_results[scanF] = sorted(indexed_taggraph_results[scanF], key = lambda item: -item[2][1])

def calculateClassProbabilities(item, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models):
    pos_probs = calculateModelProbabilities(item, spectrum_match_models[0], mod_models[0], context_models[0], db_match_models[0], protein_count_models[0], missed_cleavage_models[0], ppm_error_models[0])
    neg_probs = calculateModelProbabilities(item, spectrum_match_models[1], mod_models[1], context_models[1], db_match_models[1], protein_count_models[1], missed_cleavage_models[1], ppm_error_models[1])

    pos_model_prob = 1.0
    neg_model_prob = 1.0
    for i in range(len(pos_probs)):
        pos_model_prob *= pos_probs[i]
        neg_model_prob *= neg_probs[i]

    #print pos_model_prob, pos_probs
    #print neg_model_prob, neg_probs
    return pos_model_prob, neg_model_prob

def calculateModelProbabilities(item, spectrum_match_model, mod_model, context_model, db_match_model, protein_count_model, missed_cleavage_model, ppm_error_model):
    return [ spectrum_match_model(item), mod_model(item), context_model[item[0][7]], db_match_model[item[0][6]], protein_count_model[item[1][0]], missed_cleavage_model[item[0][10]], ppm_error_model[item[0][8]] ]



def writeEMProbabilities(outFileName, indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs, topOnly=False, score_cut = 0.99):
    cols = ['ScanF', 'Charge', 'Matching Tag Length', 'Spectrum Probability Score', 'Specificity', 'Modified?', 'PPM Error', 'Mod Size', 'Num Mods', 'Unique Siblings', 'Context Mod Variants', 'Num Mod Occurrences', 'Mod Class', 'Num Single Mods Found', 'Context', 'Mod Context', 'Mod Tuple', 'EM Probability', '1-lg10 EM', 'Pos Spectrum Prob', 'Neg Spectrum Prob', 'Pos Mod Prob', 'Neg Mod Prob', 'Pos Context Prob', 'Neg Context Prob', 'Pos DB Match Prob', 'Neg DB Match Prob', 'Pos Protein Count Prob', 'Neg Protein Count Prob', 'Pos Missed Cleavage Prob', 'Neg Missed Cleavage Prob', 'Pos PPM Error Prob', 'Neg PPM Error Prob', 'Pos Prob', 'Neg Prob']

    outFile = open(outFileName, 'w')
    outFile.write('\t'.join(cols) + '\n')

    for scanF in indexed_taggraph_results:
        top_item = indexed_taggraph_results[scanF][0]
        items = [top_item]

        for item in indexed_taggraph_results[scanF][1:]:
            if not topOnly:
                items += [item]
            elif top_item[2][0] > score_cut and round(item[0][1], 1) == round(top_item[0][1], 1) and item[0][4] == top_item[0][4] and item[0][2] == top_item[0][2]:
                items += [item]
            
            
        for item in items:
            
            pos_probs = calculateModelProbabilities(item, spectrum_match_models[0], mod_models[0], context_models[0], db_match_models[0], protein_count_models[0], missed_cleavage_models[0], ppm_error_models[0])
            neg_probs = calculateModelProbabilities(item, spectrum_match_models[1], mod_models[1], context_models[1], db_match_models[1], protein_count_models[1], missed_cleavage_models[1], ppm_error_models[1])
                
            pos_model_prob, neg_model_prob = calculateClassProbabilities(item, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models)

            em_prob = pos_model_prob*prior_probs[0] / (neg_model_prob*prior_probs[1] + pos_model_prob*prior_probs[0])

            if em_prob == 1.0:
                lg_em_prob = -math.log10( (neg_model_prob*prior_probs[1]) / (pos_model_prob*prior_probs[0]) )
            elif em_prob == 0.0:
                lg_em_prob = (pos_model_prob*prior_probs[0]) / (neg_model_prob*prior_probs[1])
            else:
                lg_em_prob = -math.log10( 1 - em_prob )

            write_info = {}
            write_info['ScanF'] = scanF
            write_info['Charge'] = item[0][16]
            write_info['Matching Tag Length'] = item[0][6]
            write_info['Spectrum Probability Score'] = item[0][1]
            write_info['Specificity'] = item[0][7]
            write_info['PPM Error'] = item[0][8]
            write_info['Modified?'] = item[0][14]
            write_info['Num Mods'] = item[0][9]
            write_info['Mod Class'] = item[0][11]
            write_info['Mod Size'] = item[0][13]
            write_info['Unique Siblings'] = item[1][0]
            write_info['Num Single Mods Found'] = item[1][3]
            write_info['Context Mod Variants'] = item[1][1]
            write_info['Num Mod Occurrences'] = item[1][2]
            write_info['Context'] = item[0][2]
            write_info['Mod Context'] = item[0][3]
            write_info['Mod Tuple'] = item[0][4]
            write_info['EM Probability'] = em_prob
            write_info['1-lg10 EM'] = lg_em_prob

            write_info['Pos Spectrum Prob'] = pos_probs[0]
            write_info['Pos Mod Prob'] = pos_probs[1]
            write_info['Pos Context Prob'] = pos_probs[2]
            write_info['Pos DB Match Prob'] = pos_probs[3]
            write_info['Pos Protein Count Prob'] = pos_probs[4]
            write_info['Pos Missed Cleavage Prob'] = pos_probs[5]
            write_info['Pos PPM Error Prob'] = pos_probs[6]
            
            write_info['Neg Spectrum Prob'] = neg_probs[0]
            write_info['Neg Mod Prob'] = neg_probs[1]
            write_info['Neg Context Prob'] = neg_probs[2]
            write_info['Neg DB Match Prob'] = neg_probs[3]
            write_info['Neg Protein Count Prob'] = neg_probs[4]
            write_info['Neg Missed Cleavage Prob'] = neg_probs[5]
            write_info['Neg PPM Error Prob'] = neg_probs[6]
                
            write_info['Pos Prob'] = pos_model_prob
            write_info['Neg Prob'] = neg_model_prob
                
            outFile.write('\t'.join([str(write_info[col]) for col in cols]) + '\n')

    outFile.close()

def writeEMProbabilitiesOnlyProbs(outFileName, indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs, topOnly=False, score_cut = 0.99):
    cols = ['ScanF', 'EM Probability', '1-lg10 EM']

    outFile = open(outFileName, 'w')
    outFile.write('\t'.join(cols) + '\n')

    for scanF in indexed_taggraph_results:
        top_item = indexed_taggraph_results[scanF][0]
        items = [top_item]

        for item in indexed_taggraph_results[scanF][1:]:
            if not topOnly:
                items += [item]
            elif top_item[2][0] > score_cut and round(item[0][1], 1) == round(top_item[0][1], 1) and item[0][4] == top_item[0][4] and item[0][2] == top_item[0][2]:
                items += [item]
            
            
        for item in items:
            
            pos_probs = calculateModelProbabilities(item, spectrum_match_models[0], mod_models[0], context_models[0], db_match_models[0], protein_count_models[0], missed_cleavage_models[0], ppm_error_models[0])
            neg_probs = calculateModelProbabilities(item, spectrum_match_models[1], mod_models[1], context_models[1], db_match_models[1], protein_count_models[1], missed_cleavage_models[1], ppm_error_models[1])
                
            pos_model_prob, neg_model_prob = calculateClassProbabilities(item, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models)

            em_prob = pos_model_prob*prior_probs[0] / (neg_model_prob*prior_probs[1] + pos_model_prob*prior_probs[0])

            if em_prob == 1.0:
                lg_em_prob = -math.log10( (neg_model_prob*prior_probs[1]) / (pos_model_prob*prior_probs[0]) )
            elif em_prob == 0.0:
                lg_em_prob = (pos_model_prob*prior_probs[0]) / (neg_model_prob*prior_probs[1])
            else:
                lg_em_prob = -math.log10( 1 - em_prob )

            write_info = {}
            write_info['ScanF'] = scanF
            write_info['EM Probability'] = em_prob
            write_info['1-lg10 EM'] = lg_em_prob

            outFile.write('\t'.join([str(write_info[col]) for col in cols]) + '\n')
    outFile.close()

def writeModels(outFileName, indexed_taggraph_results, cutOff=0.99, ep_step = 0.0025):
    num_scans = len(indexed_taggraph_results)

    outFile = open(outFileName, 'w')

    prior_probs = getPriorProbabilities(indexed_taggraph_results, num_scans)
    outFile.write('PRIOR PROBABILITIES\n')
    outFile.write(str(prior_probs) + '\n\n')

    pos_charge, pos_mean, pos_sigma, neg_charge, neg_mean, neg_sigma = getSpectrumMatchParametersMulti(indexed_taggraph_results, num_scans, prior_probs)
    outFile.write('CHARGE STATE MODELS\n')
    outFile.write('\t'.join(['Charge', 'Pos', 'Neg', 'Ratio']) + '\n')
    for i in sorted(pos_charge):
        outFile.write('\t'.join([str(i), str(pos_charge[i]), str(neg_charge[i]), str(pos_charge[i]/neg_charge[i])]) + '\n')
    
    outFile.write('\nSPECTRUM MATCH MODELS (MULTIVARIATE GAUSSIAN OVER SPECTRUM MATCH SCORE, PEPTIDE LENGTH)\n')
    for charge in sorted(pos_charge):
        outFile.write('Charge %i Pos Mean %s, Pos Sigma %s\n'%(charge, str((pos_mean[0][charge], pos_mean[1][charge])),str(pos_sigma[charge])))
        outFile.write('Charge %i Neg Mean %s, Neg Sigma %s\n'%(charge, str((neg_mean[0][charge], neg_mean[1][charge])),str(neg_sigma[charge])))

    context_models = calculateProteaseContextProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    outFile.write('\nCONTEXT MODELS\n')
    outFile.write('Pos ' + str(context_models[0]) + '\n')
    outFile.write('Neg ' + str(context_models[1]) + '\n\n')

    missed_cleavage_models = calculateMissedCleavageProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    outFile.write('\nMISSED CLEAVAGE MODELS\n')
    outFile.write('\t'.join(['Num Missed', 'Pos', 'Neg', 'Ratio']) + '\n')
    for i in sorted(missed_cleavage_models[0].keys()):
        outFile.write('\t'.join([str(i), str(missed_cleavage_models[0][i]), str(missed_cleavage_models[1][i]), str(missed_cleavage_models[0][i]/missed_cleavage_models[1][i])]) + '\n')

    db_match_models = calculateDBMatchProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    outFile.write('\nDB MATCH MODELS\n')
    outFile.write('\t'.join(['Acc', 'Pos', 'Neg', 'Ratio']) + '\n')
    for i in sorted(db_match_models[0].keys()):
        outFile.write('\t'.join([str(i), str(db_match_models[0][i]), str(db_match_models[1][i]), str(db_match_models[0][i]/db_match_models[1][i])]) + '\n')

    protein_count_models = calculateUniqueSiblingCountProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    outFile.write('\nUNIQUE SIBLING COUNT MODELS\n')
    outFile.write('\t'.join(['Num Siblings', 'Pos', 'Neg', 'Ratio']) + '\n')
    for i in sorted(protein_count_models[0].keys()):
        outFile.write('\t'.join([str(i), str(protein_count_models[0][i]), str(protein_count_models[1][i]), str(protein_count_models[0][i]/protein_count_models[1][i])]) + '\n')

    ppm_error_models = calculatePPMErrorProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    outFile.write('\nPPM Error DISTRIBUTIONS\n')
    outFile.write('\t'.join(['PPM', 'Pos', 'Neg', 'Ratio']) + '\n')
    for i in sorted(ppm_error_models[0].keys()):
        outFile.write('\t'.join([str(i), str(ppm_error_models[0][i]), str(ppm_error_models[1][i]), str(ppm_error_models[0][i]/ppm_error_models[1][i])]) + '\n')
    
    # Each function is should never produce a probability that is zero due to floating point inaccuracy
    pos_modified, neg_modified, pos_num_mods, neg_num_mods, pos_context_probs, neg_context_probs, pos_mod_size_probs, neg_mod_size_probs, pos_occurrence_probs, neg_occurrence_probs, pos_single_mod_of_multi, neg_single_mod_of_multi, pos_class_probs, neg_class_probs = getModProbDistributions(indexed_taggraph_results, num_scans, prior_probs)

    outFile.write('MOD MODELS\n')

    modified_map = {0: 'No', 1: 'Yes'}
    outFile.write('\nMODIFIED DISTRIBUTIONS\n')
    outFile.write('\t'.join(['Modified?', 'Pos', 'Neg', 'Ratio']) + '\n')
    for modified in pos_modified:
        outFile.write('\t'.join([modified_map[modified], str(pos_modified[modified]), str(neg_modified[modified]), str(pos_modified[modified]/neg_modified[modified])]) + '\n')
                    
                    
    outFile.write('\nNUM MODS PER PEPTIDE DISTRIBUTIONS\n')
    outFile.write('\t'.join(['Modified?', 'Num Mods', 'Pos', 'Neg', 'Ratio']) + '\n')
    for m in sorted(pos_num_mods.keys()):
        for i in sorted(pos_num_mods[m].keys()):
            outFile.write('\t'.join([modified_map[m], str(i), str(pos_num_mods[m][i]), str(neg_num_mods[m][i]), str(pos_num_mods[m][i]/neg_num_mods[m][i])]) + '\n')

    outFile.write('\nMOD CONTEXT VARIANT DISTRIBUTIONS\n')
    outFile.write('\t'.join(['Num Variants', 'Pos', 'Neg', 'Ratio']) + '\n')
    for i in sorted(pos_context_probs.keys()):
        outFile.write('\t'.join([str(i), str(pos_context_probs[i]), str(neg_context_probs[i]), str(pos_context_probs[i]/neg_context_probs[i])]) + '\n')

    outFile.write('\nMOD SIZE DISTRIBUTIONS\n')
    outFile.write('\t'.join(['Modified?', 'Size Bin', 'Pos', 'Neg', 'Ratio']) + '\n')
    for m in sorted(pos_mod_size_probs.keys()):
        for i in sorted(pos_mod_size_probs[m].keys()):
            outFile.write('\t'.join([modified_map[m], str(i), str(pos_mod_size_probs[m][i]), str(neg_mod_size_probs[m][i]), str(pos_mod_size_probs[m][i]/neg_mod_size_probs[m][i])]) + '\n')
    
    outFile.write('\nNUM MODS OCCURRENCES DISTRIBUTIONS\n')
    for num in pos_occurrence_probs:
        outFile.write('NUM MODS %i\n'%num)
        outFile.write('\t'.join(['Num Occurrences', 'Pos', 'Neg', 'Ratio']) + '\n')
        for i in sorted(pos_occurrence_probs[num].keys()):
            outFile.write('\t'.join([ str(i), str(pos_occurrence_probs[num][i]), str(neg_occurrence_probs[num][i]), str(pos_occurrence_probs[num][i]/neg_occurrence_probs[num][i]) ]) + '\n')

    outFile.write('\nNUMBER SINGLE MODS ON SAME CONTEXT (FOR MULTI-MODDED PEPTIDES)\n')
    outFile.write('\t'.join(['Num Mods', 'Num Single Mods Found', 'Pos', 'Neg', 'Ratio']) + '\n')
    for num_mods in sorted(pos_single_mod_of_multi.keys()):
        for i in sorted(pos_single_mod_of_multi[num_mods].keys()):
            outFile.write('\t'.join([str(num_mods), str(i), str(pos_single_mod_of_multi[num_mods][i]), str(neg_single_mod_of_multi[num_mods][i]), str(pos_single_mod_of_multi[num_mods][i]/neg_single_mod_of_multi[num_mods][i])]) + '\n')

    class_map = {0: 'Unmodified', 1: 'Defined Mod', 2: 'AA Sub', 3: 'Undefined Mod', 4: 'Insertion/Deletion'}
    outFile.write('\nMOD CLASS DISTRIBUTIONS\n')
    outFile.write('\t'.join(['Modified?', 'Mod Class', 'Pos', 'Neg', 'Ratio']) + '\n')
    for m in sorted(pos_class_probs.keys()):
        for i in sorted(pos_class_probs[m].keys()):
            outFile.write('\t'.join([modified_map[m], class_map[i], str(pos_class_probs[m][i]), str(neg_class_probs[m][i]), str(pos_class_probs[m][i]/neg_class_probs[m][i])]) + '\n')

    mod_counts = {}
    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        if item[0][4] not in mod_counts:
            mod_counts[item[0][4]] = {'Unique': item[1][2], 'Total': 1, 'Num Mods': len(item[0][4]), 'Over Cutoff': 0, 'Over Cutoff Unique': set(), 'Error Sum': 0}
        else:
            mod_counts[item[0][4]]['Total'] += 1

        if item[2][0] > cutOff:
            mod_counts[item[0][4]]['Over Cutoff'] += 1
            mod_counts[item[0][4]]['Over Cutoff Unique'].add(item[0][2])

            try:
                mod_counts[item[0][4]]['Error Sum'] += item[0][8]
            except TypeError:
                # If mass error is indeterminate
                pass

    outFile.write('\nNUM MODS CONSIDERED %i\n\n'%len(mod_counts) )
    cols = ['Mod', 'Length', 'Over Cutoff (%f)'%cutOff, 'Over Cutoff Unique', 'Total', 'Unique', 'Average Error']
    outFile.write('\t'.join(cols) + '\n')
    for mod in sorted(mod_counts.keys(), key= lambda k: -mod_counts[k]['Over Cutoff']):
        if mod_counts[mod]['Over Cutoff'] == 0:
            break
        
        outFile.write('\t'.join([ str(mod), str(mod_counts[mod]['Num Mods']), str(mod_counts[mod]['Over Cutoff']), str(len(mod_counts[mod]['Over Cutoff Unique'])), str(mod_counts[mod]['Total']),  str(mod_counts[mod]['Unique']), str(mod_counts[mod]['Error Sum']/mod_counts[mod]['Over Cutoff'] * ep_step) ]) + '\n')

    outFile.close()
    

# Calculate model parameters based on maximizing expectated value of results relative to hidden variable probabilities given previous iteration model parameters
def updateModelParameters(indexed_taggraph_results, num_scans, recalculate_network = True, first_iteration = False):

    # Recalculates # of other unique pepts/protein and # of unique contexts/
    if recalculate_network:
        recalculateNetworkAttributes(indexed_taggraph_results, top_only = False, first_iteration = first_iteration)
    
    # Each function is should never produce a probability that is zero due to floating point inaccuracy
    prior_probs = getPriorProbabilities(indexed_taggraph_results, num_scans)
    spectrum_match_models = calculateSpectrumMatchProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    mod_models = calculateModProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    context_models = calculateProteaseContextProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    db_match_models = calculateDBMatchProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    protein_count_models = calculateUniqueSiblingCountProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    missed_cleavage_models = calculateMissedCleavageProbabilities(indexed_taggraph_results, num_scans, prior_probs)
    ppm_error_models = calculatePPMErrorProbabilities(indexed_taggraph_results, num_scans, prior_probs)

    return spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs
    

def getUnivariateNormalPDF(mu, sigma):
    mu = float(mu)
    sigma = float(sigma)
    
    norm_const = 1.0/ math.pow( 2*np.pi * sigma, 1.0/2 )

    def prob_norm_calc(x):
        x_mu = x - mu
        return norm_const * math.pow(math.e, -x_mu**2 / (2 * sigma))

    return prob_norm_calc

# Returns multivariate normal PDF as a function of input datapoint x (given parameters mu and sigma of distribution)
def getMultivariateNormalPDF(mu, sigma):

    size = len(mu)
    sigma_mat = np.matrix(sigma).reshape((2,2))
    det = np.linalg.det(sigma_mat)
    norm_const = 1.0/ ( math.pow((2*np.pi),float(size)/2) * math.pow(det,1.0/2) )

    inv = sigma_mat.I
    inv_list = [inv[0,0], inv[0,1], inv[1,0], inv[1,1]]
    spec_score_mean, pept_length_mean = mu

    def prob_norm_calc(spec_score, pept_length):
        spec_score_mu = spec_score - spec_score_mean
        pept_length_mu = pept_length - pept_length_mean

        # Matrix multiply x.T * E * x
        conjugate = (spec_score_mu*inv_list[0] + pept_length_mu*inv_list[2])*spec_score_mu + (spec_score_mu*inv_list[1] + pept_length_mu*inv_list[3])*pept_length_mu 
        
        return norm_const * math.pow(math.e, -0.5 * conjugate)

    return prob_norm_calc

def createInitialGuess(pos_mean=(5, 10), pos_sigma=(10, -5, -5, 100), neg_mean=(-2, 10), neg_sigma=(10, -5, -5, 100), match_lengths=range(1,40), match_models = False):
    spec_match_models = (getMultivariateNormalPDF(pos_mean, pos_sigma), getMultivariateNormalPDF(neg_mean, neg_sigma))

    if not match_models:
        pos_match = defaultdict(lambda: 0.5)
        neg_match = defaultdict(lambda: 0.5)
        num_lengths = len(match_lengths)
        step_up = 0.98/num_lengths
        for i, ml in enumerate(match_lengths):
            pos_match[ml] = 0.01 + (i-1)*step_up
            neg_match[ml] = 0.99 - (i-1)*step_up

        match_models = (pos_match, neg_match)
        
    return spec_match_models, match_models

# Performs EM until difference between previous and current results is less than cutoff
# Difference is a concatenated vector of EM Probabilities and ranks of top-scoring results (from initial ranking)
# This convergence condition ensures that re-ranking has stabilized and Probability scores are converged
def performEM(indexed_taggraph_results, initial_spectrum_match_models, initial_db_match_models, params, cutOff=0.00001, max_iter_all=100, max_iter_after_rerank=100, outBase = '', write_data = True, write_model_every_time = False, rerank=True):
    num_scans = len(indexed_taggraph_results)

    #precalculateEMAttributes(indexed_taggraph_results, params['Enzyme']['specificity'])
    initializeEMProbabilities(indexed_taggraph_results, initial_spectrum_match_models, initial_db_match_models)

    num_iterations = 1
    prev_vec = np.array([0.0] * 2 * num_scans)
    # This counts as first iteration
    spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs = updateModelParameters(indexed_taggraph_results, num_scans, recalculate_network = True, first_iteration = True)
    calculateEMProbabilities(indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs)
    curr_vec = []
    for scanF in indexed_taggraph_results:
        curr_vec += [indexed_taggraph_results[scanF][0][2][0], float(indexed_taggraph_results[scanF][0][0][0])]
    curr_vec = np.array(curr_vec)
    dist = np.linalg.norm(prev_vec - curr_vec)

    if write_model_every_time:
        writeModels(outBase + '_MODELS_iter1.log', indexed_taggraph_results)
        writeEMProbabilitiesOnlyProbs(outBase + '_EMProbs_iter1.tdv', indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs, topOnly=True)

    while num_iterations < max_iter_all and dist > cutOff:
        spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs = updateModelParameters(indexed_taggraph_results, num_scans, recalculate_network = True)
        calculateEMProbabilities(indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs)

        num_iterations += 1
        prev_vec = curr_vec
        curr_vec = []
        for scanF in indexed_taggraph_results:
            curr_vec += [indexed_taggraph_results[scanF][0][2][0], float(indexed_taggraph_results[scanF][0][0][0])]
        curr_vec = np.array(curr_vec)
        dist = np.linalg.norm(prev_vec - curr_vec)

        print 'Iteration %i, Distance %f'%(num_iterations, dist)

    if write_data:
        writeModels(outBase + '_MODELS_BEFORERERANK.log', indexed_taggraph_results)
        writeEMProbabilities(outBase + '_EMProbs_BEFORERERANK.tdv', indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs)
    
    # Correction - corrects spectra where lower ranked annotation has much higher mod prevalence than highest ranked annotation if the scores are within a certain ratio (i.e., 0.001)
    # This step is done twice in order to correctly rerank multi-modded peptides (counts for single modded peptides are corrected in first step, then used to correct multi-modded peptides n second step)
    if rerank:
        rankMostPrevalentModOnTop(indexed_taggraph_results, report_file = outBase + '_RERANKSTATS_ONE.tdv', write_data = write_data)
        recalculateNetworkAttributes(indexed_taggraph_results)
        rankMostPrevalentModOnTop(indexed_taggraph_results, report_file = outBase + '_RERANKSTATS_TWO.tdv', write_data = write_data)

    recalculateNetworkAttributes(indexed_taggraph_results, top_only = True)

    num_iterations = 1
    while num_iterations < max_iter_after_rerank and dist > cutOff:
        spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs = updateModelParameters(indexed_taggraph_results, num_scans, recalculate_network = False)
        calculateEMProbabilities(indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs, rerank = False, score_all = False)

        num_iterations += 1
        prev_vec = curr_vec
        curr_vec = []
        for scanF in indexed_taggraph_results:
            curr_vec += [indexed_taggraph_results[scanF][0][2][0], float(indexed_taggraph_results[scanF][0][0][0])]
        curr_vec = np.array(curr_vec)
        dist = np.linalg.norm(prev_vec - curr_vec)

        if write_model_every_time:
            writeModels(outBase + '_MODELS_iter%i.log'%(num_iterations), indexed_taggraph_results)
            writeEMProbabilitiesOnlyProbs(outBase + '_EMProbs_iter%i.tdv'%(num_iterations), indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs, topOnly=True)

        print 'Iteration %i, Distance %f'%(num_iterations, dist)

    # Kludge, will probably remove this conditional after finishing cross-validation analysis
    if write_data or write_model_every_time:
        # calculate probs for all items (so that items with the same prob as top 'ranked' item can also be marked as top
        # calculateEMProbabilities(indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, prior_probs, rerank = False)
        writeModels(outBase + '_MODELS_END.log', indexed_taggraph_results)
        writeEMProbabilities(outBase + '_EMProbs_END_TOPONLY.tdv', indexed_taggraph_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, missed_cleavage_models, ppm_error_models, prior_probs, topOnly=True)

    return spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, ppm_error_models, prior_probs

def getTopRankAndProbVectors(indexed_taggraph_results):

    items, ranks, probs = [], [], []
    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]
        items += [copy.copy(item)]
        ranks += [int(item[0][0])]
        probs += [float(item[2][0])]

    return items, ranks, probs

def updateModContextWinCounts(top_item, items, mod_context_win_counts):
    # Get best scoring single mods for context, will be used later to identify most likely mod localization site on a given context
    top_item_context = top_item[0][2]
    mod_score_dict = defaultdict(lambda: -1000)
    for item in items:
        if item[0][2] == top_item_context and len(item[0][4]) == 1:
            mod_score_dict[item[0][4]] = max(mod_score_dict[item[0][4]], round(item[0][1], 2))


    sorted_mods = sorted(mod_score_dict.items(), key = lambda k: -k[1])
    if len(sorted_mods) > 0:
        max_score = sorted_mods[0][1]
        for mod, score in sorted_mods:
            if score < max_score:
                break
            mod_context_win_counts[top_item_context][mod] += 1

def rankMostPrevalentModOnTop(indexed_taggraph_results, min_log_score_diff = 1, report_file = False, write_data = False):

    mod_stats = defaultdict(lambda: {'Prob Sum': 0, 'Total': 0, 'Unique': set(), 'Over 99': 0, 'Unique Over 99': set(), 'Total Top': 0, 'Prob Sum Top': 0, 'Over 99 Top': 0, 'Unique Over 99 Top': set()})
    mod_without_loc_counts = defaultdict(int)
    mod_context_win_counts = defaultdict(lambda: defaultdict(int))
    
    for scanF in indexed_taggraph_results:
        item = indexed_taggraph_results[scanF][0]

        if len(item[0][4]) == 1:
            updateModContextWinCounts(item, indexed_taggraph_results[scanF], mod_context_win_counts)
        
        top_item_em_prob = item[2][1]
        top_item_context = item[0][2]
        
        for item in indexed_taggraph_results[scanF]:

            if item[0][2] == top_item_context and item[2][1] == top_item_em_prob:
                mod_stats[item[0][4]]['Total Top'] += 1
                mod_stats[item[0][4]]['Prob Sum Top'] += item[2][0]
                
                if item[2][0] >= 0.99:
                    mod_stats[item[0][4]]['Over 99 Top'] += 1
                    mod_stats[item[0][4]]['Unique Over 99 Top'].add(item[0][2])
                
                    # Count occurences of mods in peptide without localization (i.e., to know to prioritize carbamidomethylation over gly, etc.)
                    for mod in item[0][4]:
                        mod_without_loc_counts[mod[:2]] += 1
        

            mod_stats[item[0][4]]['Prob Sum'] += item[2][0]
            
            mod_stats[item[0][4]]['Total'] += 1
            mod_stats[item[0][4]]['Unique'].add(item[0][2])
            if item[2][0] >= 0.99:
                mod_stats[item[0][4]]['Over 99'] += 1
                mod_stats[item[0][4]]['Unique Over 99'].add(item[0][2])

      
    for scanF in indexed_taggraph_results:
        # First do switchWithTop with all items which have the same context as top item (to get the 'best' mod annotation for that context), then do switchWithTop with all other mods
        items = indexed_taggraph_results[scanF]
        top_item_context = items[0][0][2]
        top_item_mods = items[0][0][4]
        for i in range(1, len(items)):
            if items[i][0][2] != top_item_context or len(items[i][0][4]) != len(top_item_mods):
                continue
            print 'CANDIDATE', scanF, items[0][0][3], items[0][0][4], items[i][0][3], items[i][0][4]
            if switchWithTop(items[0], items[i], mod_stats, mod_without_loc_counts, mod_context_win_counts, min_log_score_diff):
                print 'SWITCHING'
                items[i], items[0] = items[0], items[i]

        for i in range(1, len(items)):
            if items[i][0][2] != top_item_context or len(items[i][0][4]) == len(top_item_mods):
                continue
            print 'CANDIDATE', scanF, items[0][0][3], items[0][0][4], items[i][0][3], items[i][0][4]
            if switchWithTop(items[0], items[i], mod_stats, mod_without_loc_counts, mod_context_win_counts, min_log_score_diff):
                print 'SWITCHING'
                items[i], items[0] = items[0], items[i]

        for i in range(1, len(items)):
            if items[i][0][2] == top_item_context:
                continue
            print 'CANDIDATE', scanF, items[0][0][3], items[0][0][4], items[i][0][3], items[i][0][4]
            if switchWithTop(items[0], items[i], mod_stats, mod_without_loc_counts, mod_context_win_counts, min_log_score_diff):
                print 'SWITCHING'
                items[i], items[0] = items[0], items[i]
        
    if report_file and write_data:
        outFile = open(report_file, 'w')
        cols = ['Mod Tuple', 'Prob Sum', 'Total', 'Unique', 'Total Top', 'Prob Sum Top', 'Over 99', 'Unique Over 99', 'Over 99 Top', 'Unique Over 99 Top', 'Without Loc Counts']
        outFile.write('\t'.join(cols) + '\n')
        for mod_tuple in sorted(mod_stats.keys(), key = lambda k: -len(mod_stats[k]['Unique Over 99']))[:10000]:
            outFile.write('\t'.join([
                str(mod_tuple),
                str(mod_stats[mod_tuple]['Prob Sum']),
                str(mod_stats[mod_tuple]['Total']),
                str(len(mod_stats[mod_tuple]['Unique'])),
                str(mod_stats[mod_tuple]['Total Top']),
                str(mod_stats[mod_tuple]['Prob Sum Top']),
                str(mod_stats[mod_tuple]['Over 99']),
                str(len(mod_stats[mod_tuple]['Unique Over 99'])),
                str(mod_stats[mod_tuple]['Over 99 Top']),
                str(len(mod_stats[mod_tuple]['Unique Over 99 Top'])),
                str(sum([mod_without_loc_counts[mod[:2]] for mod in mod_tuple]))
                ]) + '\n')
            
            
        outFile.close()

def hasAASub(mod_tuple):
    return not all([isNotAASub(item[0]) for item in mod_tuple])

# Switch with top if:
# top is a modified peptide and item is unmodified (always)
# TODO: Add constraint to switch with top mod if same mod is found on other residues (i.e., carbamidomethylation and oxidation)
# TODO: For multiply modded peptides, rerank based on some heuristic if mods on peptide are also found as highly prevalent single mods (i.e., mod prevalence is max of mod prevalences of single mods or combo mod) except if score difference is greater than min_log_score_diff?
def switchWithTop(top_item, item, mod_stats, mod_without_loc_counts, mod_context_win_counts, min_log_score_diff = 1, min_diff_for_context_win = 5, min_spec_score_diff_same_context = 0.1, mod_error_diff = 10, min_unique_over_cut_generic = 10, min_unique_over_cut_same_context = 20, min_unique_count_diff_to_switch = 1, ratio_for_same_context_switch = 50.0):
    # Don't switch if current item has an AA substitution (NOTE: This biases against AA subs)
    #if hasAASub(item[0][4]):
    #    return False

    top_spec_score = round(top_item[0][1], 2)
    item_spec_score = round(item[0][1], 2)

    # Handles cases where one of the items being compared is an exact match - want to always choose the best exact match in this case
    if len(item[0][4]) == 0:
        if len(top_item[0][4]) > 0:
            return True
        elif top_spec_score < item_spec_score:
            return True
        else:
            return False

    if len(item[0][4]) > 0 and len(top_item[0][4]) == 0:
        return False

    # Sometimes a multiply modded pept outscores a singly modified peptide for the same scanF, even though all of the mods in the less modified peptide are also present in the more modified result. The less modified result is almost universally better. This occurs because multiply modified peptides get a boost if single modded counterparts are found elsewhere
    #print item[0][15], top_item[0][15], set(item[0][15]) in set(top_item[0][15])
    if len(top_item[0][4]) > len(item[0][4]) and set(item[0][15]).issubset(set(top_item[0][15])):
        print 'Switching based on item being subset of top'
        return True

    item_mods_only = tuple([mod for mod in item[0][4] if mod not in top_item[0][4]])
    top_item_mods_only = tuple([mod for mod in top_item[0][4] if mod not in item[0][4]])
    indiv_item_mod_counts = [ len(mod_stats[(mod,)]['Unique Over 99 Top']) for mod in item[0][4] ]
    indiv_top_item_mod_counts = [ len(mod_stats[(mod,)]['Unique Over 99 Top']) for mod in top_item[0][4] ]
    
    try:
        min_mods_not_in_top_item_counts = min([len(mod_stats[(mod,)]['Unique Over 99 Top']) for mod in item_mods_only])
    except ValueError:
        # Both item and top item have same mods, but probably in different places. Will likely never need to switch in this case unless mod localization is somehow important
        min_mods_not_in_top_item_counts = min_unique_over_cut_same_context
        
    try:
        min_mods_only_in_top_item_counts = min([len(mod_stats[(mod,)]['Unique Over 99 Top']) for mod in top_item_mods_only])
    except ValueError:
        min_mods_only_in_top_item_counts = min_unique_over_cut_same_context

    # Deals with cases where top mod is on top just due to prevalence erroneously (i.e., AA sub such as Gly->Ala). Hopefully does not pick suboptimal mod when both mods are on same residue (i.e., Phospho vs Sulfo, Guanidinyl vs Acetyl)
    # Biases against AA substitutions, relaxed switching constraint if top item has AA subs, does not trigger if current item has AA sub
    # Only execute this code if unmod contexts are the same and the two items to compare have the same number of mods
    # Spectrum score is often not enough to overwhelm other scoring attributes in this case, so we give it a little 'help' (could also be wrong localization due to spurious peak matches, so we'll need to watch this to make sure it doesn't return junk 
    if item[0][2] == top_item[0][2]:
        print 'entering same context switch'
        # Get minimum mod counts over all mods in the current item which are not in top item (normalizes for case when both top_item and item share a 'rare' mod such as GlycerylPE)

        # For mass error constraint, make sure that item mod error is less than mod error diff or that top item mod error is less than mod error diff
        # Also, mod class of item is 0 or 1 or mod class of top item is greater than 1
        # print item[0][11] < 2 or top_item[0][11] >= 2, len(item_mods_only) == 1 and len(top_item_mods_only) == 1, (top_item[0][8] == 'Indeterminate' or item[0][8] == 'Indeterminate' or abs(item[0][8]) < mod_error_diff or abs(top_item[0][8]) > mod_error_diff)
        if len(item_mods_only) == 1 and len(top_item_mods_only) == 1:
            print 'Both mods counts one', top_item_mods_only, mod_context_win_counts[item[0][2]][top_item_mods_only], item_mods_only, mod_context_win_counts[item[0][2]][item_mods_only], sum(indiv_top_item_mod_counts), sum(indiv_item_mod_counts), min_mods_only_in_top_item_counts, min_mods_not_in_top_item_counts, top_item[0][10], item[0][10]
            if (top_item[0][8] == 'Indeterminate' or item[0][8] == 'Indeterminate' or abs(item[0][8]) <= mod_error_diff or abs(top_item[0][8]) > mod_error_diff):
                if mod_context_win_counts[item[0][2]][item_mods_only] - abs(round((min_mods_only_in_top_item_counts - min_mods_not_in_top_item_counts)/ratio_for_same_context_switch)) - 2*(item[0][10] - top_item[0][10]) > mod_context_win_counts[item[0][2]][top_item_mods_only] and (top_item[0][8] == 'Indeterminate' or (item[0][8] != 'Indeterminate' and abs(item[0][8]) <= abs(top_item[0][8])) ):
                    # Makes sure that the current mod isn't switched for a REALLY abundant one (i.e., Xle->Arg for Carbamyl, Insertion of G for Carbamidomethyl, etc.)
                    # Also add two to # of context wins required to switch for each missed cleavage more which is present in the item relative to the top item
                    print 'Switching based on context win'
                    return True

                if mod_context_win_counts[item[0][2]][item_mods_only] >= mod_context_win_counts[item[0][2]][top_item_mods_only]:
                    # Make sure that mod context win counts for item mod at least equals that for top mod even for prevalence based switching
                    if item_spec_score >= top_spec_score and sum(indiv_item_mod_counts) - min_unique_count_diff_to_switch >= sum(indiv_top_item_mod_counts):
                        print 'Switching based on mod counts'
                        return True

                    # To distinguish Gly from carbamidomethyl, etc.
                    # Only when both items have same mod context
                    if item[0][3] == top_item[0][3] and item_spec_score >= top_spec_score and item_mods_only[0][1] == top_item_mods_only[0][1] and sum([ mod_without_loc_counts[mod[:2]] for mod in item[0][4] ]) > sum([ mod_without_loc_counts[mod[:2]] for mod in top_item[0][4] ]):
                        print 'Switching based on mod without loc counts'
                        return True

            return False
        elif item_spec_score >= top_spec_score and min_mods_not_in_top_item_counts >= min_unique_over_cut_same_context and sum(indiv_item_mod_counts)/len(indiv_item_mod_counts) > sum(indiv_top_item_mod_counts)/len(indiv_top_item_mod_counts):
            print 'Switching based on average mod count greater', min_mods_not_in_top_item_counts >= min_unique_over_cut_same_context,  sum(indiv_item_mod_counts)/len(indiv_item_mod_counts) > sum(indiv_top_item_mod_counts)/len(indiv_top_item_mod_counts)
            # TODO: For comparing a more modified protein with a less modified protein, add constraint that spectrum score must be strictly greater or that number of single mods found is greater than 0?
            # For comparing multiply modded peptides, or peptides with different numbers of mods on same context
            return True
        elif len(item_mods_only) != 0 and len(top_item_mods_only) > len(item_mods_only) and item_spec_score >= top_spec_score and (min_mods_only_in_top_item_counts < min_unique_over_cut_same_context or min_mods_not_in_top_item_counts > min_mods_only_in_top_item_counts) and ( (len(item[0][4]) == 1 and top_item[1][3] == 0) or (len(item[0][4]) > 1 and item[1][3] >= top_item[1][3]) ):
            # For case when combo mod is on top and single mod is candidate for rerank. Want to be really careful here, but allow rerank if top combo mod is highly implausible (i.e., one of the mods in the combo is 'rare', or the item to rerank has more prevalent minumum mods than the top mod and no single mod from combo mod has been found at same context). Must also pass a spectrum score cutoff
            print 'Switching based on top mod min less prevalent than cutoff or item mod', min_mods_only_in_top_item_counts < min_unique_over_cut_same_context or min_mods_not_in_top_item_counts > min_mods_only_in_top_item_counts, len(item[0][4]) == 1 and top_item[1][3] == 0, (len(item[0][4]) == 1 and top_item[1][3] == 0) or (len(item[0][4]) > 1 and item[1][3] >= top_item[1][3])
            return True
            #elif len(top_item_mods_only) != 0 and len(item_mods_only) > len(top_item_mods_only) and min_mods_not_in_top_item_counts > min_unique_over_cut_same_context and ( (len(top_item[0][4]) == 1 and item[1][3] > 0) or item[1][3] > top_item[1][3] ):
            # For case when expanded combo mod is candidate for rerank and top_item has less mods, only rerank if all mods in combo mod are prevalent and more single mods have been found in context for item than top_item (different from above scenario as this does not require the spectrum score of the combo mod to be greater)
            #return True
        else:
            print 'same context switch false'
            return False

    
    print 'entering diff context switch'
    # Generic switching comparator based on number of unique hits with EM prob > 0.99 for mod, Only triggers in cases where above comparators fail to return (i.e., different contexts, no exact match, mods for one item not entirely contained in mods for the other item)
    # Only switches if context for one of the peptide candidates is within context for another peptide candidate (i.e., Succinyl vs Carbamyl at same 'site'), does not switch if contexts are on different proteins entirely
    more_prevalent = False
    if (len(mod_stats[item[0][4]]) >= min_unique_over_cut_generic and len(mod_stats[top_item[0][4]]['Unique Over 99 Top']) < len(mod_stats[item[0][4]]['Unique Over 99 Top'])):
        more_prevalent = True
    if len(top_item_mods_only) > 0 and len(item[0][4]) <= len(top_item[0][4]) and min_mods_not_in_top_item_counts >= min_unique_over_cut_generic and min_mods_not_in_top_item_counts > min_mods_only_in_top_item_counts:
        more_prevalent = True

    same_context = (top_item[0][2].replace('.', '') in item[0][2].replace('.', '')) or (item[0][2].replace('.', '') in top_item[0][2].replace('.', ''))

    # Comparing singly to doubly modded peptides sometimes has large differences in EM Prob. EM Probs have to be close here UNLESS top item has similar mass error to item or top item is multiply modified
    if same_context and more_prevalent and (top_item[2][1]-item[2][1] < min_log_score_diff or len(top_item[0][4]) > 1 or (item[0][8] != 'Indeterminate' and (top_item[0][8] == 'Indeterminate' or abs(top_item[0][8]) >= abs(item[0][8]))) ):
        print 'Switching based on diff context more counts', len(mod_stats[item[0][4]]) >= min_unique_over_cut_generic and len(mod_stats[top_item[0][4]]['Unique Over 99 Top']) < len(mod_stats[item[0][4]]['Unique Over 99 Top']), len(top_item_mods_only) > 0 and min_mods_not_in_top_item_counts >= min_unique_over_cut_generic and min_mods_not_in_top_item_counts > min_mods_only_in_top_item_counts
        return True
    else:
        return False


# Performs one step of EM algorithm
def EMStep(indexed_taggraph_results, num_scans):
    gaussian_models, mod_models, context_models, prior_probs = updateModelParameters(indexed_taggraph_results, num_scans)
    calculateEMProbabilities(indexed_taggraph_results, gaussian_models, mod_models, context_models, prior_probs)

    return gaussian_models, mod_models, context_models, prior_probs

def getSensitivityAndPrecisionArrs(scoreArr, critArr, critCutOff=0.9, numPoints=100):
    truePosInds = np.where(critArr >= critCutOff)
    trueScores = scoreArr[truePosInds]
    falseScores = np.delete(scoreArr, truePosInds)

    print scoreArr
    maxScore = np.amax(scoreArr)
    minScore = np.amin(scoreArr)
    scoreCutoffs = np.linspace(minScore, maxScore, num=numPoints)

    sens, prec = [], []
    for cutoff in scoreCutoffs:
        TP = len(np.where(trueScores >= cutoff)[0])
        FP = len(np.where(falseScores >= cutoff)[0])
        sens += [float(TP)/len(trueScores)]
        try:
            prec += [float(TP)/(TP + FP)]
        except ZeroDivisionError:
            prec += [1]

    return sens, prec, scoreCutoffs

# outFile can be file handler object or string indicating file to write to
def writeSensitivityAndPrecisionToFile(scoreArr, sensArr, precArr, outFile, allScoresArr = None):
    
    cols = ['Score Cutoff', 'Sensitivity', 'Precision', 'Num Results']
    if allScoresArr != None:
        cols += ['Num Total Results']

    try:
        outFile.write('\t'.join([col for col in cols]) + '\n')
    except:
        outFile = open(outFile, 'w')
        outFile.write('\t'.join([col for col in cols]) + '\n')

    cutOffs = np.linspace(np.amin(scoreArr), np.amax(scoreArr), num=len(sensArr))
    for i in range(cutOffs.size):
        outFile.write('\t'.join([str(cutOffs[i]), str(sensArr[i]), str(precArr[i]), str( np.where(scoreArr >= cutOffs[i])[0].size ), str( np.where(allScoresArr >= cutOffs[i])[0].size ) if allScoresArr != None else 'None']) + '\n')
    
    outFile.close()

def getTrueAndFalseScores(scoreArr, critArr, critCutOff=0.9):
    truePosInds = np.where(critArr >= critCutOff)
    trueScores = scoreArr[truePosInds]
    falseScores = np.delete(scoreArr, truePosInds)

    return trueScores, falseScores

def calculateFDRArray(scoreArr, critArr, critCutOff=0.9, numPoints=100):
    trueScores, falseScores = getTrueAndFalseScores(scoreArr, critArr, critCutOff)

    FDRs = []
    scoreCutoffs = np.linspace(np.amin(scoreArr), np.amax(scoreArr), num=numPoints)
    for cutoff in scoreCutoffs:
        TP = len(np.where(trueScores >= cutoff)[0])
        FP = len(np.where(falseScores >= cutoff)[0])
        try:
            FDR = float(FP)/(TP + FP)
        except ZeroDivisionError:
            if FP > 0:
                FDR = 1
            else:
                FDR = 0

        FDRs += [FDR]
    
    return scoreCutoffs, FDRs

def writeFDRArr(scoreCutoffs, FDRArr, outFileName):
    outFile = open(outFileName, 'w')
    
    cols = ['Score Cutoff', 'FDR']
    outFile.write('\t'.join([col for col in cols]) + '\n')

    for i in range(scoreCutoffs.size):
        outFile.write('\t'.join([str(scoreCutoffs[i]), str(FDRArr[i])]) + '\n')

    outFile.close()



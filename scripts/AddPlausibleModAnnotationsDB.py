'''
Created on Sep 29, 2011

@author: Arun
'''

#slin 201707  added lib path and others.

import os
import sys
import numpy as np
import glob
import copy
PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1,PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1,LIB_DIR)
import ArgLib
import DataFile
import Validator
import StringFoldingWrapper as SFW
import Constants
import TAG_GRAPH_PROB_SCORE as TG
import ProbNetwork as PN
sys.path.insert(2, os.path.abspath(os.path.join(PAR_DIR, 'database')))
from Models import Experiment, Result, Fraction
from sqlalchemy import create_engine, desc
from sqlalchemy.sql import bindparam, select, and_
from collections import defaultdict

experiments = Experiment.__table__
results = Result.__table__
fractions = Fraction.__table__

def getModInfo(connection, fracs, spectrum_score_cutoff = 0):
    indexed_results = defaultdict(lambda: -1000)
    inf = float("inf")
    top_mods_info = []
    all_mods_info = defaultdict(list)
    for fraction_id, fraction_name in fracs:
        stmt = select([results.c.scan, results.c.composite_score, results.c.mods, results.c.spectrum_score]).where(results.c.fraction_id == fraction_id).where(results.c.spectrum_score > spectrum_score_cutoff).order_by(results.c.scan).order_by(desc(results.c.composite_score))
        response = connection.execution_options(stream_results=True).execute(stmt)
        for row in SFW.string_folding_wrapper(response):
            scan_name = fraction_name + ':' + str(row[0])
            mods = [mod for mod in eval(row[2]) if mod[0][0] != 'Isobaric Substitution']
            mod_tuple = Validator.getUniqueModTuple(mods, undef_mod_round_precision = 2)
            if len(mod_tuple) == 0:
                # Make sure no values can be added for a scan if an unmodified peptide has been found
                indexed_results[scan_name] = inf
                all_mods_info[scan_name] = tuple()
                continue

            if includeInHash(mod_tuple):
                try:
                    all_mods_info[scan_name] += [(mods, mod_tuple, row[3], round(mods[0][0][1]))]
                except TypeError:
                    continue
                
                score = round(row[1], 2)
                if score >= indexed_results[scan_name]:
                    top_mods_info += [(mods, mod_tuple)]
                    indexed_results[scan_name] = score

    for scan_name in all_mods_info.keys():
        if len(all_mods_info[scan_name]) == 0:
            del all_mods_info[scan_name]

    return top_mods_info, all_mods_info, len(indexed_results)

def getUndefModMassesFromTopResults(connection, fracs, prob_cutoff = 0.99):
    frac_ids = [fracs[0] for frac in fracs]

    undef_mods = defaultdict(list)

    stmt = select([results.c.mods]).where(results.c.fraction_id.in_(frac_ids)).where(results.c.em_probability >= prob_cutoff)
    response = connection.execute(stmt)

    for row in response:
        mods = [mod for mod in eval(row[0]) if mod[0][0] != 'Isobaric Substitution']

        if len(mods) == 1 and mods[0][0][0] == 'Undefined Mass Shift':
            undef_mods[Validator.getUniqueModTuple(mods)[0]].add(mods[0])

    # Now calculate masses of mods as average of undef mods binned together
    for undef_mod in undef_mods:
        mod_mass = sum([mod[0][1] for mod in undef_mods[undef_mod]])/len(undef_mods[undef_mod])
        undef_mods[undef_mod] = ('Undefined Mass Shift', mod_mass, ())

    return undef_mods

#TODO: Do a small number (~5) of EM iterations with rerank before this to get more accurate estimations of most prevalent mods and abundances?
def getModExpansionCandidatesList(top_mods_info, all_mods_info, num_scans, mod_tolerance = 0.1, ep_step = 0.01, mod_combo_e_value = 0.1, undef_mod_combo_e_value = 0.1, min_counts_for_add = 10, out_base = None):
    
    expand_mods = defaultdict(int)
    initial_mod_counts = defaultdict(int)
    hashed_mod_combos = defaultdict(list)
    all_single_mods = set()

    max_int_ep = int( np.round(mod_tolerance/ep_step) )

    # Get initial mod counts
    for mods, mod_tuple in top_mods_info:
        initial_mod_counts[mod_tuple[0]] += 1

    # Pick most parsimonious candidates for each scan (if count of items is 1, just pick top scoring result)
    # This is to avoid mods like succinyl, etc. not being expanded due to also being shadow mods
    all_undef_mods = defaultdict(set)
    for scan_name in all_mods_info:
        max_score = max([initial_mod_counts[mod_item[1][0]] for mod_item in all_mods_info[scan_name]])
        items = all_mods_info[scan_name]
        # Add mod to all_single_mods set for later expansion using combo mods
        for mods, mod_tuple, spectrum_score, mod_mass in items:
            all_single_mods.add(mod_tuple[0])

        # Record highest scoring mod mass and score so that same mass mods with same score can also be counted (i.e., count Methyl, N-term as well as Methyl, H, etc.)
        if max_score > 1:
            max_items = [mod_item for mod_item in items if initial_mod_counts[mod_item[1][0]] == max_score]
            top_mod_mass = max_items[0][3]
            max_spec_score = max([mod_item[2] for mod_item in items if mod_item[3] == top_mod_mass])
        else:
            items = [items[0]]
            max_spec_score = items[0][2]
            top_mod_mass = items[0][3]

        for mods, mod_tuple, spectrum_score, mod_mass in items:
            if mod_mass == top_mod_mass and spectrum_score == max_spec_score:
                expand_mods[mod_tuple[0]] += 1

    
    # Get 'true' mass of undefined modifications (original mod_tuple function simply rounds to nearest indicated decimal precision (default: tenth)
    # TODO: Is binning necessay here, or can we just record mod tuples with a stricter precision (i.e, hundredths. Only thing is that this would thorw off counts in expand_mods because undef mods wouldn't be binned as often, decreasing their ability to act as candidates in expand_combo_mods and expand_single_mods, or making them be expanded unnecessarily in expand_combo_mods because counts are artificially low)
    # NOTE: Currently, undefined mods are not allowed to be candidates in expand_single_mods
    """
    all_undef_mod_tuple_map = {}
    for mod_tuple in all_undef_mods:
        true_mass = sum([ mod[0][1] for mod in all_undef_mods[mod_tuple] ]) / len(all_undef_mods[mod_tuple])
        new_mod_tuple = ('Undefined Mass Shift', true_mass, ())
        all_single_mods.add(new_mod_tuple)
        all_undef_mod_tuple_map[new_mod_tuple] = mod_tuple

    undefined_mod_tuple_map = {}
    for mod_tuple in undef_mods:
        true_mass = sum([ mod[0][1] for mod in undef_mods[mod_tuple] ]) / len(undef_mods[mod_tuple])
        new_mod_tuple = ('Undefined Mass Shift', true_mass, ())
        expand_mods[new_mod_tuple] = len(undef_mods[mod_tuple])
        undefined_mod_tuple_map[new_mod_tuple] = mod_tuple
    """

    # hash mods candidates for expansion in dictionary for easy lookup later
    combo_mod_cutoff = mod_combo_e_value * num_scans
    undef_combo_mod_cutoff = undef_mod_combo_e_value * num_scans
    expand_mod_list = sorted(expand_mods.items(), key=lambda k: -k[1])
    print 'Hashing mods by mass for %i plausible mods. Combo mod cutoff is %f'%( len(expand_mod_list), combo_mod_cutoff)

    hashed_mod_combos = getHashedModCombos(all_single_mods, expand_mods, combo_mod_cutoff, undef_combo_mod_cutoff, min_counts_for_add, max_int_ep, ep_step)
    
    # Remake hashed single mods, this time only using mods with a mod count higher than min_counts_for_add
    hashed_single_mods = defaultdict(list)
    for mod_item in expand_mod_list:
        mod_1 = mod_item[0]

        # Don't add mod as a candidate for single mod substitution if it is an AA Sub
        #if not Validator.isNotAASub(mod_1[0]):
        #    continue
        
        if expand_mods[mod_1] < min_counts_for_add:
            break
        
        h_mass = hashMass(mod_1[1], ep_step)
        for int_ep in range(-max_int_ep, max_int_ep+1):
            hashed_single_mods[h_mass+int_ep] += [(mod_1, int_ep*ep_step)]


    # Add mod_tuples for undefined mods to expand_mods (currently replaced by accurate mass mod_tuples for the purposes of proper mod enumeration)
    # This will allow proper prioritization and filtering of add mod candidates
    #for new_mod_tuple in undefined_mod_tuple_map:
    #    expand_mods[undefined_mod_tuple_map[new_mod_tuple]] = expand_mods[new_mod_tuple]


    out_file = open(out_base + '_poss_combo_mods.tdv', 'w')
    out_file.write('\t'.join(['Mod', 'Combo Candidates']) + '\n')
    for mod in hashed_mod_combos:
        # Note: won't output proper value for mod counts for mod variable if it is an undefined mass shift (but reflects accurately how the program treats undefined mods)
        out_file.write('\t'.join([str( (mod, expand_mods[mod]) ), str([(item, expand_mods[item[0]], expand_mods[item[1]]) for item in hashed_mod_combos[mod]])]) + '\n')
    out_file.close()

    out_file = open(out_base + '_poss_single_mods.tdv', 'w')
    out_file.write('\t'.join(['Mass', 'Poss Mods']) + '\n')
    for hash_mass in sorted(hashed_single_mods.keys()):
        out_file.write('\t'.join([str(hash_mass * ep_step), str([(item, expand_mods[item[0]]) for item in hashed_single_mods[hash_mass]])]) + '\n')
    out_file.close()
    
    return expand_mods, hashed_single_mods, hashed_mod_combos

def getHashedModCombos(all_single_mods, expand_mods, combo_mod_cutoff, undef_combo_mod_cutoff, min_counts_for_add, max_int_ep, ep_step):
    # Initially make hashed_single_mods without using min_counts_for_add cutoff, because we need to hash all mod masses to calculate mod combos to expand
    expand_mod_list = sorted(expand_mods.items(), key=lambda k: -k[1])
    all_single_mods_hashed = defaultdict(list)
    hashed_mod_combos = defaultdict(list)
    
    for mod_tuple in all_single_mods:
        h_mass = hashMass(mod_tuple[1], ep_step)
        for int_ep in range(-max_int_ep, max_int_ep+1):
            all_single_mods_hashed[h_mass+int_ep] += [(mod_tuple, int_ep*ep_step)]
            
    for i, mod_item in enumerate(expand_mod_list):
        for j in range(i, len(expand_mod_list)):
            mod_1, counts_1 = expand_mod_list[i]
            mod_2, counts_2 = expand_mod_list[j]
            if counts_1 * counts_2 < combo_mod_cutoff:
                break

            # Don't add mod combo if count for mod is less than min_counts_for_add
            if counts_1 < min_counts_for_add or counts_2 < min_counts_for_add:
                continue

            # Don't add combo mod as candidate if mod_1 or mod_2 are AA subs
            #if not (Validator.isNotAASub(mod_1[0]) and Validator.isNotAASub(mod_2[0])):
            #    continue

            h_mass = hashMass(mod_1[1] + mod_2[1], ep_step)
            for single_mod, error in all_single_mods_hashed[h_mass]:
                # Only expand mod into combo if both mods in combo are more prevalent than single mod in question
                if expand_mods[single_mod] < counts_1 and expand_mods[single_mod] < counts_2:
                    hashed_mod_combos[single_mod] += [(mod_1, mod_2, error)]

    return hashed_mod_combos

# Replace mod rules
# Try to replace all single mods with double mods (OF DEF AND UNDEF MODS ONLY, NOT AA SUBS AND INDELS)
# Try and replace single AA subs, insertions, deletions, undef, and def mods with DEF AND UNDEF MODS of similar mass ones IF not already present (DO NOT ADD AA SUBS AND INDELS HERE)
def addPlausibleCandidatesFromModList(connection, fracs, expand_mods, data_dir, hashed_single_mods, hashed_mod_combos, prob_network, ep_step = 0.01, mod_tolerance = 0.1, ppmSTD = 10, isobaric_mod_penalty = -0.5, def_mod_penalty = -1, indel_penalty = -3, undef_mod_penalty = -3, spectrum_score_cutoff = 0, max_per_scan = 10):
    
    for fraction_id, fraction_name in fracs:
        # Load in dta info
        frac_num = int(fraction_name[1:])
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        dta_dir = glob.glob(data_dir + os.path.sep + '*f%02d'%frac_num)[0] + os.path.sep
        '''
        dta_dir = glob.glob(data_dir + '/' + '*f%02d'%frac_num)[0] + '/'
        dtaList = glob.glob(dta_dir + '*.dta')
        scanFDict = DataFile.getScanFDict(dtaList)

        # Get TAG-GRAPH results
        stmt = select([results.c.scan, results.c.alignment_score, results.c.context, results.c.mod_context, results.c.mods, results.c.mod_ranges, results.c.mod_ambig_edges, results.c.proteins, results.c.matching_tag_length, results.c.de_novo_peptide, results.c.unmod_de_novo_peptide, results.c.de_novo_score, results.c.num_matches, results.c.obs_mh, results.c.retention_time]).where(results.c.fraction_id == fraction_id).where(results.c.spectrum_score > spectrum_score_cutoff).order_by(results.c.scan).order_by(desc(results.c.composite_score))
        response = connection.execution_options(stream_results=True).execute(stmt)
        
        indexed_results = defaultdict(list)
        for row in SFW.string_folding_wrapper(response):
            indexed_results[row[0]] += [row[1:]]

        new_scan_items = {}
        for scanF in indexed_results:
            # TODO: Don't add mod candidates for crappy results to save time (use spectrum prob score for this or protease specificity of localization?)
            # Can also do initial rounds of EM and use a probability cutoff (see above idea for gating mod candidates)

            # Eval mod lists and mod ranges once (will be using them over and over)
            mod_lists = []
            mod_ranges_list = []
            mod_tuples_list = []
            enumerated_mods = defaultdict(set)
            exact_match = False
            for item in indexed_results[scanF]:
                mods = eval(item[3])
                mod_lists += [mods]
                mod_ranges = tuple(eval(item[4]))
                mod_ranges_list += [ mod_ranges ]
                if len(mods) == 0 or all([mod[0][0] == 'Isobaric Substitution' for mod in mods]):
                    exact_match = True
                    break

                mod_tuples = []
                for j, mod in enumerate(mods):
                    mod_tuple = Validator.getUniqueModTuple([mod], undef_mod_round_precision = 2)[0]
                    enumerated_mods[(item[1], mod_ranges[j])].add( mod_tuple )
                    mod_tuples += [mod_tuple]
                    
                mod_tuples_list += [mod_tuples]

            #print fraction_name, scanF, exact_match
            # Don't add mod candidates if exact match is found
            if exact_match:
                continue

            #if scanF != 5841 or fraction_name != 'F12':
            #    continue

            #print '-----------------------'
            #print indexed_results[scanF]
            # Add mod candidates which can plausibly be the sum of two separate mods
            new_combo_mods = addComboModCandidates(scanF, indexed_results[scanF], mod_lists, mod_ranges_list, mod_tuples_list, enumerated_mods, scanFDict, expand_mods, hashed_mod_combos, prob_network, ep_step, mod_tolerance, ppmSTD)
            #print 'Num Combo Mods after getUniqueCandidates', sum([len(val[1]) for val in new_combo_mods])
            #print '---Combo Mods---'
            #print new_combo_mods

            #print enumerated_mods
            new_single_mods = addSingleModCandidates(scanF, indexed_results[scanF], mod_lists, mod_ranges_list, mod_tuples_list, enumerated_mods, scanFDict, expand_mods, hashed_single_mods, prob_network, ep_step, mod_tolerance, ppmSTD)
            #print 'Num Single Mods after getUniqueCandidates', sum([len(val[1]) for val in new_single_mods])
            #print '---Single Mods---'
            #print new_single_mods

            
            
            new_scan_items[scanF] = new_single_mods + new_combo_mods
            #print scanF, new_scan_items[scanF]

        #print 'Indexed results scans', sorted(indexed_results.keys())
        #print 'new_scan_items scans', sorted(new_scan_items.keys())
        # Import new candidates into DB
        values = []
        for scanF in indexed_results:
            #print scanF, scanF in new_scan_items, int(scanF) in new_scan_items, str(scanF) in new_scan_items
            if scanF in new_scan_items:
                #print 'NUM', len(new_scan_items[scanF])
                indexed_item = indexed_results[scanF][0]
                i = 0
                # Sort candidates by sum of prevalences of mods
                # No need to sort by composite_score, only top scoring mods for each (context, mod_tuple) pair are included in the new_scan_items (filtered using getUniqueCandidates)
                for item in sorted( new_scan_items[scanF], key=lambda k: -sum([expand_mods[mod] for mod in k[0][1]])/len(k[0][1]) ):
                    for candidate in item[1]:
                        candidate.update({
                            "scan": scanF,
                            "charge": scanFDict[scanF]['charge'],
                            "matching_tag_length": indexed_item[7],
                            "time_taken": 0,
                            "de_novo_peptide": indexed_item[8],
                            "unmod_de_novo_peptide": indexed_item[9],
                            "de_novo_score": indexed_item[10],
                            "num_matches": indexed_item[11],
                            "obs_mh": indexed_item[12],
                            "retention_time": indexed_item[13],
                            "ppm": (candidate["theo_mh"] - indexed_item[12])/candidate["theo_mh"] * 1000000,
                            "fraction_id": fraction_id                
                            })
                        
                        values += [candidate]   
                        i += 1

                    if i > max_per_scan:
                        break
                    
        print 'Adding %i candidates for fraction %s'%(len(values), fraction_name)        
        res = connection.execute(results.insert(), values)

    return new_scan_items

# aa_spread argument indicates how far from mod_range to search for viable candidates
def addSingleModCandidates(scanF, scan_items, mod_lists, mod_ranges_list, mod_tuples_list, enumerated_mods, scanFDict, expand_mods, hashed_single_mods, prob_network, ep_step=0.01, mod_tolerance=0.1, ppmSTD=10, aa_spread=3):
    #print enumerated_mods
    #for i in range(len(scan_items)):
    #    print scan_items[i]
    #    print mod_lists[i]
    #    print mod_tuples_list[i]
    #    print '-----------------'

    candidates_map = {}
    mod_range_map = {}
    # print enumerated_mods
    for i, item in enumerate(scan_items):
        # print 'Scan', scan_items[i]
        for j, mod in enumerate(mod_lists[i]):

            if mod[0][0] == 'Isobaric Substitution':
                continue

            context = item[1]
            seq_without_ends = context[2:-2]
            mod_range = mod_ranges_list[i][j]

            if (context, mod_range, mod[0][0]) in candidates_map:
                continue

            if j > 0:
                start_cut = max([mod_ranges_list[i][j-1][1], mod_range[0] - aa_spread])
            else:
                start_cut = max([0, mod_range[0] - aa_spread])                

                
            try:
                end_cut = min([mod_range[1] + aa_spread, mod_ranges_list[i][j+1][0]])
            except IndexError:
                end_cut =  min([mod_range[1] + aa_spread, len(seq_without_ends)])


            repl_mods = []
            replace_seqs = []
            mod_class = Validator.getModClass([mod_tuples_list[i][j]])
            #has_defined_mod = any([ Validator.getModClass([mod_tuple]) == 1 for mod_tuple in enumerated_mods[(context, mod_range)] ])
            # Don't expand if mod is an AASub and there already exists a defined mod for the same mod_range/context
            #if mod_class > 1 and has_defined_mod:
            #    continue
                                  
            # print 'Mod', j, mod, start_cut, end_cut
            for repl_mod, error in hashed_single_mods[ hashMass(mod[0][1], ep_step) ]:
                # Don't enumerate mod if already enumerated or if mod error exceeds mod_tolerance
                repl_mod_error = (0 if not mod[0][2] else mod[0][2]) + error
                # print repl_mod, repl_mod_error, repl_mod[0],  mod_tuples_list[i][j], mod_tuples_list[i][j][0], expand_mods[repl_mod], expand_mods[mod_tuples_list[i][j]]
                # print enumerated_mods[(context, mod_range)], context, mod_range
                if repl_mod[0] == 'undefined mass shift' or repl_mod in enumerated_mods[(context, mod_range)] or abs(repl_mod_error) > mod_tolerance:
                    # Don't add candidate if candidate already exists for scan item or mod_error of replacement exceeds tolerance or if it is an undefined mod
                    continue
                elif mod_class == 1 and not (repl_mod[0] == mod_tuples_list[i][j][0] and expand_mods[repl_mod] > expand_mods[mod_tuples_list[i][j]]):
                    # If mod is a defined mod, only expand if the replacement mod has the same name but a different and more prevalent localization
                    continue

                # print 'replace candidate', mod, repl_mod, error, repl_mod_error, context, seq_without_ends[start_cut:end_cut]
                repl_mods += [(repl_mod, repl_mod_error)]

            
            term = getTerminus(start_cut, end_cut, seq_without_ends)
            subseq = seq_without_ends[start_cut:end_cut]
            for repl_mod, error in repl_mods:
                # add candidates
                try:
                    locs = getModLocs(subseq, term, repl_mod)
                    for loc in locs[1]:
                        replace_seqs += getSingleModSeq(subseq, repl_mod, loc, error)
                except IndexError:
                    # Index error occurs when subseq is blank. This happens if an insertion is sandwhiched between two mods exactly or is sandwhiched between a mod and a terminus. Doesn't happen normally, but can happen in special cases, such as when the first amino acid in the context is an oxidized methionine and the Insertion is N-terminal to this
                    # TODO: See if this case happens frequently enough with carbamidomethyls to mess up the results
                    pass

            if len(replace_seqs) > 0:
                candidates_map[(context, mod_range, mod[0][0])] = replace_seqs
                mod_range_map[(context, mod_range)] = (start_cut, end_cut)

    # print candidates_map
    new_scan_items = []
    if len(candidates_map) > 0 and scanF in scanFDict:
        precMass = scanFDict[scanF]['precMass']
        epSTD = ppmSTD * precMass * 10**-6
        spec = PN.Spectrum(prob_network, precMass, Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=DataFile.getMassIntPairs(scanFDict[scanF]['dta']), useMemo=True)
        spec.initializeNoiseModel()
                                        
        # Generate new peptide candidates
        for i, item in enumerate(scan_items):
            for j, mod_range in enumerate(mod_ranges_list[i]):
                # add candidate
                    
                if (item[1], mod_range, mod_lists[i][j][0][0]) in candidates_map:
                    # print '----------------'
                    # print item[1]
                    # print mod_range
                    # print mod_lists[i][j]
                    # print mod_ranges_list
                    # for candidate in candidates_map[(item[1], mod_range, mod_lists[i][j][0][0])]:
                    #     print candidate

                    new_scan_items += [getSingleModCandidate(spec, scanFDict[scanF]['charge'], item, mod_lists[i], mod_ranges_list[i], j, mod_range_map[(item[1], mod_range)], candidate[0], candidate[1], candidate[2], candidate[3]) for candidate in candidates_map[(item[1], mod_range, mod_lists[i][j][0][0])]]
                
    #print 'Num Single Mods Before getUniqueCandidates', len(new_scan_items)
    return getUniqueCandidates(new_scan_items)
                                                             

def getSingleModCandidate(spec, charge, item, mods, mod_ranges, replace_index, replace_interval, replace_subseq, replace_ambig_edges, replace_mod, mass_error, modAASymbol = '-'):
    unmod_pept = list(item[1][2:-2])
    #print item, mods, replace_index, replace_interval, replace_subseq, replace_mod
    
    # Need to keep separate index for ambiguous edges since insertion and deletion mods don't generate them
    ambig_edges = eval(item[5])
    ambig_edge_index = 0
    
    new_mod_pept = []
    new_mods = []
    new_ambig_edges = []
    
    mod_ranges = list(copy.copy(mod_ranges))
    mod_ranges[replace_index] = replace_interval
    
    prev_mod_range = (0, 0)
    for i, mod_range in enumerate(mod_ranges):
        new_mod_pept += unmod_pept[prev_mod_range[1]:mod_range[0]]
        prev_mod_range = mod_range
        
        if i == replace_index:
            new_mod_pept += [replace_subseq]
            
            mod, mod_index = replace_mod

            # Construct new mod annotation string
            new_mods += [((mod[0], mod[1], mass_error), mod[2], mod_range[0]+mod_index)]
            new_ambig_edges += replace_ambig_edges
            
            if mods[i][0][0] != 'Insertion' and mods[i][0][0] != 'Deletion':
                ambig_edge_index += 1
        else:

            new_mods += [mods[i]]
            if mods[i][0][0] == 'Insertion':
                new_mod_pept += [mods[i][1]]
            elif mods[i][0][0] == 'Deletion':
                pass
            else:
                if mods[i][0][0] != 'Isobaric Substitution':
                    unmod_pept[mods[i][2]] = modAASymbol
                    new_ambig_edges += [ambig_edges[ambig_edge_index]]
                    ambig_edge_index += 1
                new_mod_pept += unmod_pept[mod_range[0]:mod_range[1]]

                    
    new_mod_pept += unmod_pept[mod_range[1]:]
    
    #print ''.join(new_mod_pept), new_ambig_edges, new_mods
    spectrum_score, pm = TG.scoreMatch(spec, charge, ''.join(new_mod_pept), new_ambig_edges, new_mods)
    
    return {
        'alignment_score': item[0],
        'spectrum_score': spectrum_score,
        'composite_score': item[0] + spectrum_score,
        'context': item[1],
        'theo_mh': pm + Constants.mods['H2O'] + Constants.mods['H+'],
        'mod_context': item[1][:2] + ''.join(new_mod_pept) + item[1][-2:],
        'mods': new_mods,
        'mod_ranges': mod_ranges,
        'proteins': item[6],
        'mod_ambig_edges': new_ambig_edges
    }            
    



def addComboModCandidates(scanF, scan_items, mod_lists, mod_ranges_list, mod_tuples_list, enumerated_mods, scanFDict, expand_mods, hashed_mod_combos, prob_network, ep_step = 0.01, mod_tolerance= 0.1, ppmSTD = 10):
    add_mods_map = defaultdict(set)

    # Go through entire candidate list, identify alternate combo mod interpretations for given mod ranges
    for i, item in enumerate(scan_items):
        mod_ranges = mod_ranges_list[i]
        
        for j, mod in enumerate(mod_lists[i]):
            # Continue if mod has already been expanded
            # Format of key in add_mods_map is (context, mod_range)
            if mod[0][0] == 'Insertion' or mod[0][0] == 'Deletion' or mod[0][0] == 'Isobaric Substitution':
                continue

            # print j, mod, mod_ranges, mod_ranges[j], item[1]
            # Initialize set so that this can be skipped if it comes up in the future (and no candidates are found)
            add_mods_map[(item[1], mod_ranges[j], mod_tuples_list[i][j])] = set()
            # now hash mass of mods in peptides to see if alternate combo candidates can be found
            for mod_combo_candidate in hashed_mod_combos[ mod_tuples_list[i][j] ]:

                # ModError is defined as mass of de_novo_seq - mass of modified reference_seq
                mod_error = (0 if not mod[0][2] else mod[0][2]) - mod_combo_candidate[-1]
                # make sure that mod_error is within tolerance and no expanded mod for given mod interval has greater prevalence than either mod in mod combo
                if not (abs(mod_error) > mod_tolerance or any([ expand_mods[enum_mod] > expand_mods[mod_combo_candidate[0]] or expand_mods[enum_mod] > expand_mods[mod_combo_candidate[1]] for enum_mod in enumerated_mods[(item[1], mod_ranges[j])] ])):
                     
                    add_mods_map[ (item[1], mod_ranges[j], mod_tuples_list[i][j]) ].add( (mod_combo_candidate, mod_error) )

    #print 'Add mods', add_mods_map
    # Get Sequence candidates for mod ranges which have valid combo mods
    candidates_map = {}
    # print add_mods_map.keys()
    for context, mod_range, mod_tuple in add_mods_map:
        candidates = []
        term = getTerminus(mod_range[0], mod_range[1], context)
            
        for mod_combo in add_mods_map[(context, mod_range, mod_tuple)]:
            mod_1, mod_2 = mod_combo[0][:2]
            subseq = context[2:-2][mod_range[0]:mod_range[1]]
            locs_1 = getModLocs(subseq, term, mod_1)
            locs_2 = getModLocs(subseq, term, mod_2)
            # TODO: Only add mod combo with most prevalent single mods for a given (context, mod_range) combination (as opposed to all valid mod combos as is done now)?
            for loc_1 in locs_1[1]:
                for loc_2 in locs_2[1]:
                    if loc_1 < loc_2:
                        candidates += getComboModSeq(subseq, mod_1, loc_1, mod_2, loc_2, mod_combo[1])
                    elif loc_1 > loc_2:
                        candidates += getComboModSeq(subseq, mod_2, loc_2, mod_1, loc_1, mod_combo[1])
                    elif locs_1[0] != locs_2[0] and (mod_1[0] == 'undefined mass shift' or mod_2[0] == 'undefined mass shift' or mod_1[2][1] != mod_2[2][1]):
                        # Second part of if statement guards against things like putting a trimethyl (A, N-term) with a Carbamyl (N-term, N-term) on the same residue
                        candidates += getComboModSeq(subseq, mod_1, loc_1, mod_2, loc_2, mod_combo[1])

        if len(candidates) > 0:
            candidates_map[(context, mod_range, mod_tuple)] = candidates

    #print candidates_map
    # Note: the way this code is written now, this method might produce LOTS of duplicate entries, particularly if the peptide in question is multiply modified
    # This is because the mod_range which has a valid combo mod is expanded for each scan item (in each scan item, the position of the mod within the mod_range may be different, but it expands out to the same set of new candidates)
    # In the future (if this way is too slow), we can minimize this time by caching the position, mod for each enumerated candidate, and only expand if the set of all mods (not including the mod_range to expand) is unique
    # For now, redundant candidates are filtered out at return step
    new_scan_items = []
    if len(candidates_map) > 0 and scanF in scanFDict:
        precMass = scanFDict[scanF]['precMass']
        epSTD = ppmSTD * precMass * 10**-6
        spec = PN.Spectrum(prob_network, precMass, Nmod=0.0, Cmod=0.0, epsilon=2*epSTD, spectrum=DataFile.getMassIntPairs(scanFDict[scanF]['dta']), useMemo=True)
        spec.initializeNoiseModel()
        
        # Generate new entries for scan from peptides
        for i, item in enumerate(scan_items):
            for j, mod_range in enumerate(mod_ranges_list[i]):
                if (item[1], mod_range, mod_tuples_list[i][j]) in candidates_map:
                    new_scan_items += [getComboModCandidate(spec, scanFDict[scanF]['charge'], item, mod_lists[i], mod_ranges_list[i], j, candidate[0], candidate[1], candidate[2], candidate[3], candidate[4]) for candidate in candidates_map[(item[1], mod_range, mod_tuples_list[i][j])]]

    #print 'Num Combo Mods before getUniqueCandidates', len(new_scan_items)
    return getUniqueCandidates(new_scan_items)
                
def getComboModCandidate(spec, charge, item, mods, mod_ranges, replace_index, replace_subseq, replace_ambig_edges, replace_mod1, replace_mod2, mass_error, modAASymbol = '-'):

    unmod_pept = list(item[1][2:-2])
    #print item, mods, mod_ranges, replace_index
    #print replace_subseq, replace_ambig_edges, replace_mod1, replace_mod2
    
    # Need to keep separate index for ambiguous edges since insertion and deletion mods don't generate them
    ambig_edges = eval(item[5])
    ambig_edge_index = 0
    
    new_mod_pept = []
    new_mods = []
    new_ambig_edges = []
    
    prev_mod_range = (0, 0)
    for i, mod_range in enumerate(mod_ranges):
        new_mod_pept += unmod_pept[prev_mod_range[1]:mod_range[0]]
        prev_mod_range = mod_range
        
        if i == replace_index:
            new_mod_pept += [replace_subseq]

            mod1, index1 = replace_mod1
            mod2, index2 = replace_mod2
            # Construct new mod annotation string
            new_mods += [((mod1[0], mod1[1], mass_error), mod1[2], mod_range[0]+index1), ((mod2[0], mod2[1], mass_error), mod2[2], mod_range[0]+index2)]
            new_ambig_edges += replace_ambig_edges

            if mods[i][0][0] != 'Insertion' and mods[i][0][0] != 'Deletion':
                ambig_edge_index += 1
        else:
            new_mods += [mods[i]]
            if mods[i][0][0] == 'Insertion':
                new_mod_pept += [mods[i][1]]
            elif mods[i][0][0] == 'Deletion':
                pass
            else:
                if mods[i][0][0] != 'Isobaric Substitution':
                    unmod_pept[mods[i][2]] = modAASymbol
                    new_ambig_edges += [ambig_edges[ambig_edge_index]]
                    ambig_edge_index += 1
                new_mod_pept += unmod_pept[mod_range[0]:mod_range[1]]


    new_mod_pept += unmod_pept[mod_range[1]:]

    # Use old mods here so spectrum score isn't penalized for the inclusion of the new expanded modifications
    # TODO: Increase spectrum/composite score if new mod combination is prevalent (i.e., add E[new_mod] - 1 to score or some similar metric to make sure a good combo is ranked first)
    # print ''.join(new_mod_pept), new_ambig_edges, new_mods
    spectrum_score, pm = TG.scoreMatch(spec, charge, ''.join(new_mod_pept), new_ambig_edges, new_mods)

    return {
        'alignment_score': item[0],
        'spectrum_score': spectrum_score, 
        'composite_score': item[0] + spectrum_score,
        'context': item[1],
        'theo_mh': pm + Constants.mods['H2O'] + Constants.mods['H+'],
        'mod_context': item[1][:2] + ''.join(new_mod_pept) + item[1][-2:],
        'mods': new_mods, 
        'mod_ranges': mod_ranges,
        'proteins': item[6],
        'mod_ambig_edges': new_ambig_edges
        }

def getUniqueCandidates(candidate_list, spec_score_cut_for_list = 0):

    unique_candidates = {}
    for candidate in candidate_list:
        item_key = ( candidate['context'], Validator.getUniqueModTuple(candidate['mods']) )
        if item_key not in unique_candidates:
            unique_candidates[item_key] = [candidate]
        elif unique_candidates[item_key][0]['spectrum_score'] < candidate['spectrum_score']:
            unique_candidates[item_key] = [candidate]
        elif unique_candidates[item_key][0]['spectrum_score'] == candidate['spectrum_score'] and candidate['spectrum_score'] > spec_score_cut_for_list and candidate['mod_context'] not in [item['mod_context'] for item in unique_candidates[item_key]]:
            unique_candidates[item_key] += [candidate]

    return stringify(unique_candidates.items())

# Needed for import into DB to work properly
def stringify(candidate_list):

    for item in candidate_list:
        for candidate in item[1]:
            candidate['mods'] = str(candidate['mods'])
            candidate['mod_ranges'] = str(candidate['mod_ranges'])
            candidate['mod_ambig_edges'] = str(candidate['mod_ambig_edges'])

    return candidate_list

def getTerminus(start_ind, end_ind, context):
    if start_ind == 0:
        term = 'N-term'
    elif end_ind == len(context):
        term = 'C-term'
    else:
        term = None

    return term

# index1 must be less than or equal to index2
def getComboModSeq(unmodSeq, mod1, index1, mod2, index2, mass_error, modAASymbol = '-', min_mod_AA_mass = 10):
    seq = list(unmodSeq)
    mass1 = mod1[1]
    mass2 = mod2[1]
    
    # If both masses go on same residue
    if index1 == index2:
        mod_mass = mass1 + mass2
        modAA = seq[index1]
        seq[index1] = modAASymbol
        if Constants.aminoacids[modAA][2] + mod_mass > min_mod_AA_mass:
            return [(''.join(seq), [(0, Constants.aminoacids[modAA][2] + mod_mass)], (mod1, index1), (mod2, index2), mass_error)]
        else:
            return []
    else:
        modAA1 = seq[index1]
        modAA2 = seq[index2]
        seq[index1] = modAASymbol
        seq[index2] = modAASymbol
        
        if Constants.aminoacids[modAA1][2] + mass1 > min_mod_AA_mass and Constants.aminoacids[modAA2][2] + mass2 > min_mod_AA_mass:
            return [(''.join(seq), [(0, Constants.aminoacids[modAA1][2] + mass1), (0, Constants.aminoacids[modAA2][2] + mass2)], (mod1, index1), (mod2, index2), mass_error)]
        else:
            return []

def getSingleModSeq(unmodSeq, mod, index, mass_error, modAASymbol = '-', min_mod_AA_mass = 10):
    seq = list(unmodSeq)
    mod_mass = mod[1]

    modAA = seq[index]
    seq[index] = modAASymbol
    if Constants.aminoacids[modAA][2] + mod_mass > min_mod_AA_mass:
        return [(''.join(seq), [(0, Constants.aminoacids[modAA][2] + mod_mass)], (mod, index), mass_error)]
    else:
        return []

def getModLocs(peptide, term, mod):
    #print peptide, term, mod
    if mod[0] == 'undefined mass shift':
        return ( 'AA', range(len(peptide)) )
    else:
        if mod[2][0] == term:
            if term == 'N-term':
                locs = ('N-term', [0])
            else:
                locs = ('C-term', [len(peptide) - 1])
        else:
            locs = ( 'AA', [loc[0] for loc in TG.getModLocations(peptide, term, [mod[2]])] )
        
        return locs

def includeInHash(mod_tuple):
    return len(mod_tuple) == 1 and mod_tuple[0][0] not in ['insertion', 'deletion']

def hashMass(modMass, epStep=0.01):
    return int(np.round(modMass/epStep))
    


def writeModOccurrences(output_name, indexed_results):

    outFile = open(output_name, 'w')
    cols = ['Mod Tuple', 'Num', 'Num Unique']
    outFile.write('\t'.join(cols) + '\n')

    mod_counts = defaultdict(lambda: {'Total': 0, 'Unique': set()})

    for scanF in indexed_results:
        context, mod = indexed_results[scanF]
        mod_counts[mod]['Total'] += 1
        mod_counts[mod]['Unique'].add(context)

    for mod in sorted(mod_counts, key= lambda k: -len(mod_counts[k]['Unique'])):
        outFile.write('\t'.join([str(mod), str(mod_counts[mod]['Total']), str(len(mod_counts[mod]['Unique']))]) + '\n')

    outFile.close()
                    
    
if __name__ == '__main__':
    print 'modmaxcounts argument is number of iterations in initial EM over all results, maxcounts argument indicates maximum number of iterations before expectation-maximization terminates after reranking is completed (on the top ranked results only). Set fracs to \"all\" to run over all fractions for experiment or supply desired fracs to run EM over separated by commas'
    options = ArgLib.parse(['init', 'dtadir', 'sqlitedb', 'experimentname', 'fracs', 'model', 'config', 'modtolerance', 'ppmstd', 'output'])

    params_dict = ArgLib.parseInitFile(options.init, options)

    engine = create_engine('sqlite:///' + options.sqlitedb, echo=True)
    conn = engine.connect()
    conn.execute("PRAGMA max_page_count = max_page;");
    conn.execute("PRAGMA temp_store = 2;")
    conn.execute("PRAGMA page_size")

    out_base = os.path.splitext(options.output)[0]
    
    p_net = PN.ProbNetwork(options.config, options.model)

    try:
        experiment_id = conn.execute(select([experiments.c.id]).where(experiments.c.name == options.experimentname)).fetchone()[0]
    except TypeError:
        raise ValueError("ERROR: No experiment by name %s"%options.experimentname)

    fracs = []
    if options.fracs.lower() == 'all':
        fracs = conn.execute(select([fractions.c.id, fractions.c.name]).where(fractions.c.experiment_id == experiment_id)).fetchall()
    else:
        for fraction_name in options.fracs.split(','):
            try:
                fraction_id = conn.execute(select([fractions.c.id, fractions.c.name]).where(and_(fractions.c.name == fraction_name, fractions.c.experiment_id == experiment_id))).fetchone()[0]
                fracs += [(fraction_id, fraction_name)]
            except TypeError:
                raise ValueError("ERROR: Fraction %s not found in experiment"%fraction_name)

    
    top_mods_info, all_mods_info, num_scans = getModInfo(conn, fracs)
    expand_mod_list, hashed_single_mods, hashed_mod_combos = getModExpansionCandidatesList(top_mods_info, all_mods_info, num_scans, mod_tolerance = options.modtolerance, out_base = out_base)
    addPlausibleCandidatesFromModList(conn, fracs, expand_mod_list, options.dtadir, hashed_single_mods, hashed_mod_combos, p_net, mod_tolerance = options.modtolerance, ppmSTD = options.ppmstd)
    
    conn.close()

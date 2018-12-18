import numpy as np
import os
import sys
import csv
import copy
import pickle
import configparser
import glob
import re
from collections import defaultdict
import ftplib
import xml.etree.ElementTree as etree
import Constants
import shlex
import subprocess

csv.field_size_limit(sys.maxsize)
searchprogs = ['LADS', 'SEQUEST', 'MASCOT', 'SORCERER', 'PepNovo', 'NovoHMM', 'MaxQuant', 'LuteFisk', 'X!Tandem', 'SPIDER', 'PEAKS', 'OMSSA', 'InSpect', 'pNovo', 'Combined']

def getMotifZResults(results_file):
    return [item.split(' ')[0] for item in open(results_file).readlines()[5:]]

def getArgs(options, argList):
    args = []
    for arg in argList:
        try: 
            opt = getattr(options, arg)
            if opt != None:
                args.extend(['--' + arg, '\"' + str(opt) + '\"'])
        except AttributeError:
            print('AttributeError for arg %s'%arg)
            pass

    return args

def parse_models_file(models_file):
    in_file = open(models_file)

    in_file.readline()

    prior = eval(in_file.readline().strip())

    for i in range(2):
        in_file.readline()

    pos_gauss = in_file.readline()
    pos_mean, pos_sigma = pos_gauss.split(', Pos Sigma ')
    pos_sigma = eval(pos_sigma.strip())
    pos_mean = eval(pos_mean.split('Pos Mean ')[1])

    neg_gauss = in_file.readline()
    neg_mean, neg_sigma = neg_gauss.split(', Neg Sigma ')
    neg_sigma = eval(neg_sigma.strip())
    neg_mean = eval(neg_mean.split('Neg Mean ')[1])

    for i in range(2):
        in_file.readline()

    pos_context = eval('{' + in_file.readline().split(', {')[1].strip()[:-1])
    neg_context = eval('{' + in_file.readline().split(', {')[1].strip()[:-1])

    for i in range(5):
        in_file.readline()

    # missed cleavage models
    pos_missed_cleavage, neg_missed_cleavage = load_model_from_table(in_file)

    for i in range(2):
        in_file.readline()

    pos_match, neg_match = load_model_from_table(in_file)

    for i in range(2):
        in_file.readline()

    pos_sibling, neg_sibling = load_model_from_table(in_file, end_cond = 'MOD MODELS')

    for i in range(3):
        in_file.readline()

    pos_mod, neg_mod = load_model_from_table(in_file)

    """
    for i in range(2):
        in_file.readline()

    pos_num_mod, neg_num_mod = load_model_from_table(in_file)
    """

    return prior, pos_mean, pos_sigma, neg_mean, neg_sigma, pos_context, neg_context, pos_missed_cleavage, neg_missed_cleavage, pos_match, neg_match, pos_sibling, neg_sibling, pos_mod, neg_mod

def load_model_from_table(in_file, end_cond = ''):
    pos, neg = {}, {}
    line = in_file.readline().strip()
    while line != end_cond:
        bin, pos_v, neg_v = line.split('\t')[:3]
        pos[bin] = float(pos_v)
        neg[bin] = float(neg_v)
        line = in_file.readline().strip()

    return pos, neg


def calculateShannonEntropy(abundance_vec):
    probs = np.array(abundance_vec)/np.sum(abundance_vec)

    return -np.sum(probs * np.log(probs))

def load_david_functional_annotation_clustering(file):
    in_file = open(file).readlines()
    cols = ['Category', 'Term', 'Count', '%', 'PValue', 'Genes', 'List Total', 'Pop Hits', 'Pop Total', 'Fold Enrichment', 'Bonferroni', 'Benjamini', 'FDR']
    annotation_clusters = []

    i = 0
    while i < len(in_file):
        line = in_file[i]
        # New annotation cluster
        if 'Annotation Cluster' in line:
            cluster_num, score = line.strip().split('\t')
            cluster = {'Cluster': cluster_num.split(' ')[-1], 'Score': score.split(' ')[-1], 'Categories': []}
            i += 2
            while i < len(in_file) and in_file[i].strip() != '':
                items = in_file[i].strip().split('\t')
                cluster['Categories'] += [dict( (cols[j], items[j]) for j in range(len(items)) )]
                i += 1

            annotation_clusters += [cluster]

        i += 1

    return annotation_clusters

# Just a cut-paste from ipython, not yet argumentized
def combine_prot_lists_for_david():
    ptm_list = ['g-h1_39.99_R', 'crotonaldehyde_70.04_H', 'crotonaldehyde_70.04_K', 'cyano_-32.03_C', 'nitro_44.99_Y', 'hydroxymethyl_30.01_N', 'amino_15.01_Y', 'deamidated_0.98_R']

    ptm_map = {}

    for file in ptm_list:
        ptm_map[file] = open('./%s_prots_sponly.tdv'%file).readlines()
        out_file = open('./combined_go_list_gh1_R_crotonaldehyde_KH_cyano_C_nitro_Y_hydroxymethyl_N_amino_Y_deamidated_R.txt', 'w')

        out_file.write('\t'.join(ptm_list) + '\n')

    for i in range(max([len(k) for k in ptm_map.values()])):
        ptm_vec = []
        for file in ptm_list:
            if i >= len(ptm_map[file]):
                ptm_vec += ['']
            else:
                ptm_vec += [ptm_map[file][i].strip()]
        out_file.write('\t'.join(ptm_vec) + '\n')
    out_file.close()


def load_models_end(models_end_file, discard_unmod = True):

    # Get offset of mods table in models end file
    in_file = open(models_end_file)
    skip_lines = 0
    for line in in_file:
        if 'NUM MODS CONSIDERED' in line:
            break

        skip_lines += 1
    skip_lines += 2
    in_file.close()

    cols, mods_info = getScanInfo(models_end_file, skipLines = skip_lines, delimiter='\t')
    for mod_item in mods_info:
        mod_item['Mod'] = eval(mod_item['Mod'])
        mod_item['Average Error'] = float(mod_item['Average Error'])
        for col in ['Total', 'Unique', 'Length', 'Over Cutoff Unique', 'Over Cutoff (0.990000)']:
            mod_item[col] = int(mod_item[col])

    if discard_unmod:
        mods_info = [mod_item for mod_item in mods_info if len(mod_item['Mod']) > 0]

    return mods_info

def bin_by_tissue(accum_file, out_file_name, tissue_list_col_offset, col_map, delim='\t', combine_func = lambda vec: sum(vec)/len(vec)):
    cols = getCols(accum_file, delimiter = delim)

    new_cols = cols[:tissue_list_col_offset] + sorted(col_map.keys())
    out_file = open(out_file_name, 'w')
    out_file.write(delim.join(new_cols) + '\n')

    try:
        for item in getScanInfoIterator(accum_file, delimiter=delim):
            vec_map = defaultdict(list)
            for col in col_map:
                vecs = []
                for m_col in col_map[col]:
                    vecs += [float(item[m_col])]
                vec_map[col] = combine_func(vecs)

            scan_data = {}
            for col in cols[:tissue_list_col_offset]:
                scan_data[col] = item[col]

            for col in col_map.keys():
                scan_data[col] = vec_map[col]

            out_file.write(delim.join([str(scan_data[col]) for col in new_cols]) + '\n')
    except IndexError:
        # Some of the accumulated output files errored before they could finish, but don't really care about the super low abundant mods anyway so I just use this to gracefully exit
        pass

    out_file.close()

def transform_and_average_by_tissue(accum_file, out_file_name, tissue_list_col_offset, col_map, delim='\t', transform_func = np.log, zero_replace_func = min):
    # NOTE: BE CAREFUL WHEN USING THIS FOR STOICHIOMETRY - WAS A PEPTIDE OBSERVED FOR THE SITE OR WAS THE ZERO BECAUSE THERE WAS NO OBSERVATION
    cols = getCols(accum_file, delimiter = delim)

    new_cols = cols[:tissue_list_col_offset]
    for col in sorted(col_map.keys()):
        new_cols += ['%s Mean'%col,'%s Variance'%col]

    out_file = open(out_file_name, 'w')
    out_file.write(delim.join(new_cols) + '\n')

    exp_cols = []
    for item in col_map.values():
        exp_cols += item

    vecs = defaultdict(list)
    base_data = []
    try:
        for item in getScanInfoIterator(accum_file, delimiter=delim):
            for col in exp_cols:
                vecs[col] += [float(item[col])]

            base_data += [ [item[col] for col in cols[:tissue_list_col_offset]] ]
    except IndexError:
        # Some of the accumulated output files errored before they could finish, but don't really care about the super low abundant mods anyway so I just use this to gracefully exit
        pass

    # Replace zero values with min nonzero value for each col and transform using transform_func
    for col in vecs:
        vecs[col] = np.array(vecs[col])
        replace_val = zero_replace_func(vecs[col][np.where(vecs[col] > 0)])
        vecs[col][np.where(vecs[col] == 0)] = replace_val
        vecs[col] = transform_func(vecs[col])

    write_data = {}
    for col in col_map:
        consolidated_data = np.array([ vecs[exp_col] for exp_col in col_map[col] ])
        write_data['%s Mean'%col] = np.mean(consolidated_data, 0)
        write_data['%s Variance'%col] = np.var(consolidated_data, 0)

    for i in range(len(vecs.values()[0])):
        out_file.write('\t'.join(base_data[i] + [str(write_data[col][i]) for col in new_cols[tissue_list_col_offset:]]) + '\n')

    out_file.close()

# Note: this isn't the original method to do this, couldn't find it so I had to rewrite
def parse_unimod_xml(unimod_xml):

    unimod_dict = {}
    for event, elem in etree.iterparse(unimod_xml, events = ('start', 'end')):
        tag = strip_namespace(elem.tag)
        if event == 'start':
            if tag == 'mod':
                mod_name = elem.attrib['title']
                sites = []
                mass = 0
            elif tag == 'specificity':
                # TODO: Change to encode protein vs any term specificity in future. Currently does not consider distinction and overwrites some unimod entries as a result
                if 'N-term' in elem.attrib['position']:
                    pos = 'N-term'
                elif 'C-term' in elem.attrib['position']:
                    pos = 'C-term'
                else:
                    pos = elem.attrib['position']

                sites += [((elem.attrib['site'], pos), elem.attrib['classification'])]
            elif tag == 'delta':
                mass = float(elem.attrib['mono_mass'])
        elif event == 'end':
            if tag == 'mod':
                unimod_dict[mod_name] = {'mass': mass, 'sites': dict(sites)}

    return unimod_dict

# Note: in HUMAN PROTEOME pep.xml files, naming conventions of files are inconsistent
# some are f01,f02, etc. and others have different fraction identifiers like d01, d02, etc.
# some files are listed in sets of 12 with a preceding letter (i.e., A1, A2, A3... A12, B1, B2... B12). In order of letter then number, they map onto fractions
# sometimes a filename might appears twice with a slightly different name like _REDO, _1112334, appended if fraction was rerun. These become the next "fraction"
def convert_pepxml_to_tdv(xml_file_path, out_file_base, score_cols = {'MASCOT': ['IonScore', 'Exp Value'], 'SEQUEST': ['XCorr', 'SpScore', 'Probability']}):

    base_cols = ['Dataset', 'ScanF', 'Precursor Neutral Mass', 'Charge', 'Retention Time', 'Peptide', 'Proteins', 'Modifications']
    out_file_names = {}
    scan_data = {'Proteins': [], 'Modifications': []}
    in_results = False
    for event, elem in etree.iterparse(xml_file_path, events = ('start', 'end')):
        tag = strip_namespace(elem.tag)
        if event == 'start':
            if tag == 'msms_run_summary':
                in_results = True
            elif tag == 'search_summary' and in_results:
                search_engine = elem.attrib['search_engine']
                out_file = open(out_file_base + '_%s.tdv'%search_engine, 'w')
                out_file_names[search_engine] = out_file_base + '_%s.tdv'%search_engine

                cols = base_cols + score_cols[search_engine]
                out_file.write('\t'.join(cols) + '\n')
            elif tag == 'spectrum_query':
                scan_data['ScanF'] = elem.attrib['start_scan']
                scan_data['Precursor Neutral Mass'] = elem.attrib['precursor_neutral_mass']
                scan_data['Retention Time'] = elem.attrib['retention_time_sec']
                scan_data['Dataset'] = elem.attrib['spectrum'][:elem.attrib['spectrum'].find('.raw')]
                scan_data['Charge'] = elem.attrib['assumed_charge']
            elif tag == 'search_hit':
                scan_data['Peptide'] = elem.attrib['peptide']
                scan_data['Proteins'] += [elem.attrib['protein_descr']]
            elif tag == 'search_score':
                scan_data[elem.attrib['name']] = elem.attrib['value']
            elif tag == 'alternative_protein':
                scan_data['Proteins'] += [elem.attrib['protein_descr']]
            elif tag == 'mod_aminoacid_mass':
                scan_data['Modifications'] += [(elem.attrib['position'], elem.attrib['mass'])]
            elif tag == 'modification_info':
                try:
                    scan_data['Modifications'] += [('N-term', elem.attrib['mod_nterm_mass'])]
                except KeyError:
                    pass

        elif event == 'end':
            if tag == 'spectrum_query':
                elem.clear()
                out_file.write('\t'.join([str(scan_data[col]) for col in cols]) + '\n')
                scan_data = {'Proteins': [], 'Modifications': []}
            elif tag == 'msms_run_summary':
                in_results = False
                out_file.close()

    combine_tdv_results(out_file_names, out_file_base + '_COMBINED.tdv')

# NOTE: Currently, the summary file only reports accurate numbers in the case of TWO programs being combined
def combine_tdv_results(program_file_map, out_file_name, base_cols = ['Dataset', 'ScanF', 'Precursor Neutral Mass', 'Charge', 'Retention Time', 'Peptide', 'Proteins', 'Modifications'], prog_specific_cols = {'MASCOT': ['IonScore', 'Exp Value'], 'SEQUEST': ['XCorr', 'SpScore', 'Probability']}):

    cols = copy.copy(base_cols)
    for prog in program_file_map:
        cols += prog_specific_cols[prog]

    out_file = open(out_file_name, 'w')
    out_file.write('\t'.join(cols) + '\n')

    indexed_db_info = {}
    for prog in program_file_map:
        indexed_db_info[prog] = indexDataByCompositeKey(getScanInfo(program_file_map[prog], delimiter='\t')[1], key = lambda k: '%s:%s'%(k['Dataset'], k['ScanF']))

    all_scans = set()
    for prog in indexed_db_info:
        for scanF in indexed_db_info[prog].keys():
            all_scans.add(scanF)

    discarded = 0
    total_results = 0
    all_agree = 0
    for scan in sorted(list(all_scans)):
        prog_pepts = []
        scan_data = defaultdict(lambda: None)
        for prog in indexed_db_info:
            if scan in indexed_db_info[prog]:
                prog_pepts += [indexed_db_info[prog][scan]['Peptide']]
                for col in cols:
                    try:
                        scan_data[col] = indexed_db_info[prog][scan][col]
                    except KeyError:
                        pass

        if len(set(prog_pepts)) == 1:
            total_results += 1
            out_file.write('\t'.join([str(scan_data[col]) for col in cols]) + '\n')
            if len(prog_pepts) == len(indexed_db_info):
                all_agree += 1
        else:
            print('Discarded due to disagreement')
            for prog in indexed_db_info:
                if scan in indexed_db_info[prog]:
                    print(indexed_db_info[prog][scan])

            discarded += 1
    num_overlap = all_agree + discarded

    out_file.close()
    summary_out = open(os.path.splitext(out_file_name)[0] + '_summary.txt', 'w')
    summary_out.write('Total results: %i\n'%total_results)
    for prog in indexed_db_info:
        summary_out.write('%s Total: %i\n'%(prog, len(indexed_db_info[prog])))
    summary_out.write('All Overlap not discarded: %i\n'%all_agree)
    summary_out.write('Discarded due to disagreement: %i\n'%discarded)
    for prog in indexed_db_info:
        summary_out.write('%s Only: %i\n'%(prog, len(indexed_db_info[prog]) - num_overlap))
    summary_out.close()

def strip_namespace(tag):
    return tag[tag.find('}')+1:]

def map_human_proteome_dbresults(em_top_results_file, human_proteome_db_results, summary_out_path, num_per_frac = 100, em_cut = 0.99, min_ratio_to_accept_mapping = 0.5, score_key = 'IonScore', peptide_key = 'Peptide', include_db_mods = True, clip_flanking_aas = False, db_dataset_id = 'Dataset', db_file_delim='\t', discard_decoys = False, protein_key = 'Reference'):
    em_results = getScanInfo(em_top_results_file, delimiter='\t')[1]
    db_results = getScanInfo(human_proteome_db_results, delimiter=db_file_delim)[1]

    indexed_em = defaultdict(dict)
    for item in em_results:
        frac, scan = item['ScanF'].split(':')
        indexed_em[frac][scan] = item

    indexed_db = defaultdict(dict)
    for item in db_results:
        if discard_decoys and item[protein_key][0] == '#':
            continue

        indexed_db[item[db_dataset_id]][item['ScanF']] = item
        if clip_flanking_aas:
            item[peptide_key] = Constants.stripModifications(item[peptide_key][2:-2].replace('L', 'I'))
        else:
            item[peptide_key] = Constants.stripModifications(item[peptide_key].replace('L', 'I'))

    top_results = {}
    for dataset in indexed_db:
        top_results[dataset] = sorted( indexed_db[dataset].values(), key = lambda k: -float(k[score_key] if k[score_key] != 'None' else 0) )[:num_per_frac]

    frac_map = {}
    for dataset in top_results:
        items = top_results[dataset]
        #print(items)
        frac_agree_counts = defaultdict(int)
        for frac in indexed_em:
            for item in items:
                if item['ScanF'] in indexed_em[frac]:
                    #print(item[peptide_key], indexed_em[frac][item['ScanF']]['Context'][2:-2])
                    if item[peptide_key] == indexed_em[frac][item['ScanF']]['Context'][2:-2]:
                        frac_agree_counts[frac] += 1
        try:
            best_frac =  sorted(frac_agree_counts.items(), key = lambda k: -k[1])[0]
            if best_frac[1] >= min_ratio_to_accept_mapping * len(items):
                frac_map[dataset] = best_frac[0]
        except IndexError:
            # Happens if db fraction not present in TG data (possibly corrupted raw file or some other reason)
            print('INDEX ERROR', frac_agree_counts)
            pass

    if len(set(frac_map.values())) != len(frac_map.values()):
        print(sorted(frac_map.items(), key = lambda k: k[1]))
        #print(em_top_results_file, human_proteome_db_results)
        print('CRITICAL ERROR: REDUNDANT FRACTION MAPPING %s'%str( sorted(frac_map.items(), key = lambda k: k[1]) ))

    if len(set(frac_map.values())) != len(indexed_em.keys()):
        print(len(set(frac_map.values())), len(indexed_em.keys()))
        print('CRITICAL ERROR: Number of fractions in fraction map does not match number of fractions in TAG-GRAPH results for %s %s'%(em_top_results_file, human_proteome_db_results))

    # Get Summary Statistics
    total_spectra, agree, disagree, db_total, tg_total = 0, 0, 0, 0, 0
    disagree_items = []
    db_peptides, db_mod_peptides, tg_contexts, tg_mod_contexts = set(), set(), set(), set()
    for frac in indexed_em:
        for item in indexed_em[frac].values():
            total_spectra += 1
            if float(item['EM Probability']) >= em_cut:
                tg_total += 1
                tg_contexts.add(item['Context'][2:-2])
                tg_mod_contexts.add( (item['Context'][2:-2], eval(item['Mod Tuple'])) )

    for dataset in indexed_db:
        for item in indexed_db[dataset].values():
            db_total += 1
            db_peptides.add(item[peptide_key])
            if include_db_mods:
                db_mod_peptides.add( (item[peptide_key], tuple(eval(item['Modifications']))) )
            else:
                db_mod_peptides.add(item[peptide_key])

    for dataset in frac_map:
        frac = frac_map[dataset]
        for scanF in indexed_db[dataset]:
            if scanF in indexed_em[frac] and float(indexed_em[frac][scanF]['EM Probability']) >= em_cut:
                if indexed_db[dataset][scanF][peptide_key] == indexed_em[frac][scanF]['Context'][2:-2]:
                    agree += 1
                elif 'Met-loss+Acetyl' in indexed_em[frac][scanF]['Mod Tuple'] and indexed_db[dataset][scanF][peptide_key] in indexed_em[frac][scanF]['Context'][2:-2]:
                    # Covers case where TG reports N-terminal acetylation with initiator methionine and DB reports it as an Acetyl mod
                    agree += 1
                else:
                    disagree += 1
                    disagree_items += [(indexed_db[dataset][scanF], indexed_em[frac][scanF])]

    try:
        percent_disagree = float(disagree)/(agree + disagree) * 100
    except ZeroDivisionError:
        # Happens in Adult_Monocytes bRP_Elite because pep.xml was mislabeled (duplicate of bRP_Velos)
        percent_disagree = 'ERROR'
    db_only = db_total - (agree + disagree)
    tg_only = tg_total - (agree + disagree)
    all_total = db_only + tg_only + agree + disagree

    with open(os.path.splitext(summary_out_path)[0] + '_disagree_data.pck', 'w') as fout:
        pickle.dump(disagree_items, fout)

    return [total_spectra, db_total, tg_total, all_total, len(db_peptides), len(db_mod_peptides), len(tg_contexts), len(tg_mod_contexts), agree, disagree, percent_disagree, db_only, tg_only]



# Set interpreter to False for native programs
def executeProcess(progLoc, prog, args, interpreter='python'):

    if not interpreter:
        cmd = ' '.join([os.path.join(progLoc, prog)] + args)
    else:
        cmd = ' '.join([interpreter] + [os.path.join(progLoc, prog)] + args)
    print(cmd)

    '''
    By Default, Posix is True. This properly handles the double-quotes that ArgLib adds around arguments/paths,
    but it removes all backslashes as escape characters. This is a problem on Windows where backslash is the
    path separator character. Solution for now is to leave Posix as True, and make sure all paths use slashes
    instead of backslashes. Note that os.path.normpath will change slashes to backslashes on Windows, so that
    function cannot be used. --smp
    '''
    proc = subprocess.Popen(shlex.split(cmd, posix=True))

    proc.communicate()

def getFilesFTP(ftp_url, ftp_dir, wildcard):
    file_list = []
    ftp = ftplib.FTP(ftp_url)
    ftp.login()
    ftp.cwd(ftp_dir)

    resp = ftp.retrlines('NLST', file_list.append)
    print('Number of files in directory', len(file_list))
    ftp.quit()

    for file in file_list:
        if re.match(wildcard, file):
            # reopen ftp for each file since connection is unstable
            ftp = ftplib.FTP(ftp_url)
            ftp.login()
            ftp.cwd(ftp_dir)

            print('Retrieving file %s'%file)
            outfile = open(file, 'w')
            ftp.retrbinary("RETR " + file, outfile.write)
            outfile.close()

            ftp.quit()

'''
This method will return a tuple containing the precursor mass and charge read
from the first line of a .dta file
@param absPath: absolute path to .dta file
@return: tuple, first element is precursor mass, last element is charge
@author: Arun
'''
def getPrecMassAndCharge(absPath):
    fin = open(absPath, 'r')

    firstline = fin.readline()
    precMass, charge = firstline.split()

    fin.close()
    return float(precMass), float(charge)

def getDTAFNamesInDir(dirPath):
    dtas = []
    dirList = os.listdir(dirPath)
    for file in dirList:
        ext = os.path.splitext(file)[1]
        if ext == '.dta':
            dtas.extend([dirPath + file])

    return dtas

def getDTAByScanNum(dtaDir, scanNum):
    return glob.glob(dtaDir + '/*.%(scanNum)04i.*dta' % {'scanNum': scanNum})

def indexDataByKey(data, key='ScanF', overrideKey=None, dtyper=int):
    dataDict = {}
    for datum in data:
        if dtyper(datum[key]) in dataDict and overrideKey != None:
            if float(datum[overrideKey]) > float(dataDict[dtyper(datum[key])][overrideKey]):
                dataDict[dtyper(datum[key])] = datum
        else:
            dataDict[dtyper(datum[key])] = datum

    return dataDict

def indexDataByCompositeKey(data, key=lambda k: k['ScanF'], overrideKey=None):
    dataDict = {}
    for datum in data:
        if key(datum) in dataDict and overrideKey != None:
            if float(datum[overrideKey]) > float(dataDict[key(datum)][overrideKey]):
                dataDict[key(datum)] = datum
        else:
            dataDict[key(datum)] = datum

    return dataDict

def PEAKS7_split_by_fraction(peaks_file, frac_map = None, first_frac = 1, delim = '\t', scan_col="ScanF", renumber_fracs = True):
    base, ext = os.path.splitext(peaks_file)
    cols, data = getScanInfo(peaks_file, delimiter=delim)

    # no fraction number present in peaks file, write all data to F1
    # This allows us to run RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py on our data if we just want to run 1 raw file (or export data from peaks in a different way)
    if len(data[0][scan_col].split(':')) == 1:
        out_file = open(base + '_F1' + ext, 'w')
        out_file.write(delim.join(cols) + '\n')
        for item in data:
            out_file.write(delim.join([item[col] for col in cols]) + '\n')

        out_file.close()
        return

    # Get fractions
    if frac_map == None:
        frac_map = {}
        for item in data:
            frac_map[item[scan_col].split(':')[0]] = 0

        if renumber_fracs:
            frac_nums = sorted([int(frac[1:]) for frac in frac_map.keys()])
            frac_diff = first_frac - frac_nums[0]
            for frac_num in frac_nums:
                frac_map['F' + str(frac_num)] = open(base + '_F' + str(frac_num + frac_diff) + ext, 'w')

        else:
            for fraction in frac_map:
                frac_map[fraction] = open(base + '_' + fraction + ext, 'w')

    for fraction in frac_map:
        frac_map[fraction].write(delim.join(cols) + '\n')

    # Split Data
    for item in data:
        frac, scan = item[scan_col].split(':')
        item[scan_col] = scan
        frac_map[frac].write(delim.join([item[col] for col in cols]) + '\n')

    for fraction in frac_map:
        frac_map[fraction].close()


def combineDatafiles(data1, data2, key='ScanF', overrideKey1=None, overrideKey2=None, dtyper=int):
    indexedData1 = indexDataByKey(data1, key=key, overrideKey=overrideKey1, dtyper=dtyper)
    indexedData2 = indexDataByKey(data2, key=key, overrideKey=overrideKey2, dtyper=dtyper)

    cols1 = indexedData1.values()[0].keys()
    cols2 = indexedData2.values()[0].keys()
    try:
        cols1.remove(key)
        cols2.remove(key)
    except ValueError:
        pass

    combinedDataArr = []
    allKeyVals = set(indexedData1.keys()) | set(indexedData2.keys())
    for keyVal in allKeyVals:
        scanData = {}
        scanData[key] = keyVal

        if keyVal in indexedData1:
            for col in cols1:
                scanData[col] = indexedData1[keyVal][col]
        else:
            for col in cols1:
                scanData[col] = None

        if keyVal in indexedData2:
            for col in cols2:
                scanData[col] = indexedData2[keyVal][col]
        else:
            for col in cols2:
                scanData[col] = None

        combinedDataArr += [scanData]

    return combinedDataArr

def getCols(fname, skipLines = 0, delimiter=','):
    rec = csv.reader(open(fname, 'r'), delimiter=delimiter)
    for i in range(skipLines):
        rec.next()

    return rec.next()

def getGeneOntologyInfo(gene_ontology_file):
    return getScanInfo(gene_ontology_file, delimiter='\t', skipLines=12, cols = ['DB', 'Object ID', 'Gene Name', 'Qualifier', 'GO', 'GO_REF', 'Evidence Code', 'With (or) From', 'Aspect', 'DB Object Name', 'DB Object Synonym', 'DB Object Type', 'Taxon', 'Date', 'Assigned By', 'Annotation Extension', 'Gene Product Form ID'])[1]



def getScanInfo(fname, fields=None, delimiter=',', skipLines=0, cols = False):
    runInfo = []
    rec = csv.reader(open(fname, 'r'), delimiter=delimiter)
    for i in range(skipLines):
        rec.next()

    # Can supply cols if sheet does not have cols listed
    if not cols:
        cols = rec.next()

    if not fields:
        fields = [col for col in cols]

    fieldInds = []
    for i in range(len(cols)):
        if cols[i] in fields:
            fieldInds += [i]

    for row in rec:
        rowInfo = {}
        for ind in fieldInds:
            rowInfo[cols[ind]] = row[ind]

        runInfo += [rowInfo]

    return cols, runInfo

def getScanInfoIterator(fname, delimiter=','):
    rec = csv.reader(open(fname, 'r'), delimiter=delimiter)
    cols = rec.next()

    for row in rec:
        rowInfo = {}
        for i in range(len(cols)):
            rowInfo[cols[i]] = row[i]
        yield rowInfo

def getScanNum(dtaFName):
    fname = os.path.basename(dtaFName)
    base = os.path.splitext(fname)[0]
    base = os.path.splitext(base)[0]
    base = os.path.splitext(base)[0]
    scanNum = os.path.splitext(base)[1]
    return int(scanNum[1:])


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def getMassIntPairs(absPath):
    fin = open(absPath, 'r')

    massIntPairs = []
    fin.readline()

    for i, line in enumerate(fin):
        mass, intensity = [float(l) for l in line.split()]
        if intensity > 0:
            massIntPairs += [(mass, intensity)]

    fin.close()
    return np.array(massIntPairs)

def getMassIntData(absPath):
    fin = open(absPath, 'r')

    prec_charge = fin.readline().split('\t')
    prec_mass, charge = float(prec_charge[0]), int(prec_charge[1])

    massIntPairs = []
    for i, line in enumerate(fin):
        mass, intensity = [float(l) for l in line.split()]
        if intensity > 0:
            massIntPairs += [(mass, intensity)]

    fin.close()
    return prec_mass, charge, np.array(massIntPairs)

def getScoringMatrix(absPath):
    fin = open(absPath, 'r')
    line = fin.readline()
    while line[0] == '#':
        line = fin.readline()

    i = 0
    AAMap = {}
    for char in line:
        if char != ' ' and char != '\n':
            AAMap[char] = i
            i += 1

    dim = max(AAMap.values()) + 1
    matrix = np.zeros((dim, dim), dtype=np.dtype('int'))
    for i in range(dim):
        line = fin.readline()
        char = line[0]
        if AAMap[char] != i:
            print('ERROR: scoring matrix not properly formatted for AA %s' % char)
        else:
            data = line[2:]
            k = j = 0
            while j < dim:
                if data[k] == '-':
                    matrix[i][j] = int(data[k:k + 2])
                    k += 1
                    j += 1
                elif data[k] != ' ':
                    matrix[i][j] = int(data[k])
                    j += 1

                k += 1

    return AAMap, matrix

def getDBInfo(db, key=None):
    import anydbm as dbm
    db = dbm.open(db, 'r')
    if key:
        return pickle.loads(db[key])
    else:
        dbDict = {}
        for key in db.keys():
            dbDict[key] = pickle.loads(db[key])
        return dbDict

#    modDict = {'SEQUEST': {'silac(k)': '#', 'oxidation': '*'},
#                'MASCOT': {'silac(k)': '*', 'oxidation': '#'},
#               'LADS': {'silac(k)': '*', 'oxidation': '#'}}
#    AADict = {'SEQUEST': {'X': 'I'}, 'MASCOT': {}, 'LADS': {'-': 'X'}}
def generateSeqMap(progDict, symbolMap, paramsDict):
    seqMap = {}
    for prog in progDict:
        seqMap[prog] = {'AAs':{}, 'Mods': {}}
        for aa in Constants.origAAs.keys():
            seqMap[prog]['AAs'][aa] = aa

        if progDict[prog] != 'LADS':
            seqMap[prog]['AAs']['L'] = 'I'
        elif progDict[prog] == 'LADS':
            seqMap[prog]['AAs']['-'] = 'X'

        try:
            seqMap[prog]['AAs'].update(symbolMap[prog]['AAs'])
        except KeyError:
            print('AAs dictionary for program %s not present in symbol map' % (prog,))

    for mod in paramsDict['Diff Mods'].keys():
        modType = paramsDict['Diff Mods'][mod][0]
        for prog in progDict:
            try:
                progMod = symbolMap[prog]['Mods'][modType]
                if symbolMap[prog]['Mods'][modType] != 'static':
                    # (mod, aa_loc [n-term or c-term])
                    seqMap[prog]['Mods'][ (progMod, paramsDict['Diff Mods'][mod][1]) ] = mod
                else:
                    # AD 9/6/2016 Not sure what this clause does, leaving it in just in case it's important
                    # modify all amino acids with static mod
                    seqMap[prog]['AAs'][paramsDict['Diff Mods'][mod][1]] = seqMap[prog]['AAs'][paramsDict['Diff Mods'][mod][1]] + mod
            except KeyError:
                print('Modification of type %s unaccounted for in modDict for program %s' % (modType, prog))

    for modData in paramsDict['Static Mods'].keys():
        modType, aa_loc = modData
        for prog in progDict:
            try:
                progMod = symbolMap[prog]['Mods'][modType]
                seqMap[prog]['Mods'][(progMod, aa_loc)] = ''
            except KeyError:
                print('Modification of type %s unaccounted for in modDict for program %s' % (modType, prog))

    return seqMap

paramHandler = {

'Static Mods': lambda datum, args: Constants.addStaticMod(datum, args),
'Diff Mods': lambda datum, args: Constants.addDiffMod(datum, args),
'Models': lambda datum, args: parseModelInfo(datum, args),
'Enzyme': lambda datum, args: parseEnzymeInfo(datum, args),
'Amino Acids': lambda datum, args: Constants.addAA(datum, args),
'LADS Parameters': lambda datum, args: parseLADSParametersInfo(datum, args),
'Pair Configurations': lambda datum, args: parsePairConfiguration(datum, args),
'Cluster Configuration': lambda datum, args: parseClusterConfiguration(datum, args)
}

def parseParams(fname):
    #print(fname)
    paramsDict = {}
    params = configparser.configparser()
    params.read(fname)
    for section in params.sections():
        paramsDict[section] = {}
        for datum in params.options(section):
            try:
                paramsDict[section].update(paramHandler[section](datum, params.get(section, datum).split(' ')))
            except KeyError:
                raise KeyError('%s not a valid category. Valid categories are: %s' % (section, str(paramHandler.keys())))

    #print(paramsDict)
    paramsDict['Static AAs'] = [entry[1] for entry in paramsDict['Static Mods']]
    return paramsDict

paramHandler_v1 = {
'Static Mods': lambda datum, args: Constants.addStaticMod(datum, args),
'Diff Mods': lambda datum, args: Constants.addDiffMod_v1(datum, args),
'Models': lambda datum, args: parseModelInfo(datum, args),
'Enzyme': lambda datum, args: parseEnzymeInfo(datum, args),
'Amino Acids': lambda datum, args: Constants.addAA(datum, args),
'LADS Parameters': lambda datum, args: parseLADSParametersInfo(datum, args),
'Pair Configurations': lambda datum, args: parsePairConfiguration(datum, args),
'Cluster Configuration': lambda datum, args: parseClusterConfiguration(datum, args)
}
def parseParams_v1(fname):
    #print(fname)
    paramsDict = {}
    params = configparser.configparser()
    params.read(fname)
    for section in params.sections():
        paramsDict[section] = {}
        for datum in params.options(section):
            try:
                paramsDict[section].update(paramHandler_v1[section](datum, params.get(section, datum).split(' ')))
            except KeyError:
                raise KeyError('%s not a valid category. Valid categories are: %s' % (section, str(paramHandler_v1.keys())))

    paramsDict['Static AAs'] = [entry[1] for entry in paramsDict['Static Mods']]
    return paramsDict

def parseModelInfo(datum, args):
    return {datum: {'config': args[0], 'model': args[1]}}

def parseEnzymeInfo(datum, args):
    if datum == 'specificity':
        return {datum: parseEnzymeSpecificity(args)}
    else:
        return {datum: args[0]}

def parseLADSParametersInfo(datum, args):
    return {datum: args[0]}

def parseEnzymeSpecificity(args):
    motifs = []
    for arg in args:
        motifs += [arg.split(';')]

    return motifs



def parsePairConfiguration(datum, pairConfig):
    pairConfiguration = {'NStatic': np.array([float(pairConfig[0])]), 'CStatic': np.array([float(pairConfig[1])]), 'NMod': float(pairConfig[2]), 'CMod': float(pairConfig[3]), 'NModSymbol': pairConfig[4], 'CModSymbol': pairConfig[5], 'Model': pairConfig[6]}
    if pairConfiguration['NModSymbol'] in ['None', '0', 'False']:
        pairConfiguration['NModSymbol'] = ''
    if pairConfiguration['CModSymbol'] in ['None', '0', 'False']:
        pairConfiguration['CModSymbol'] = ''

    return {datum: pairConfiguration}

def parseClusterConfiguration(datum, args):
    return {datum: args[0]}

# Gets top-ranked result(s) plus some stats for each scan
def processTAGGRAPHInfo(scanInfo, dtyper=lambda x: x, score_key="Composite Score"):
    processedDict = {}

    for item in scanInfo:
        if dtyper(item['ScanF']) not in processedDict:
            processedDict[dtyper(item['ScanF'])] = item
            item['Context'] = [( (item['Context']), eval(item['Modifications']), eval(item['Proteins']) )]
            item['Proteins'] = set([prot.split(' ')[0] for prot in eval(item['Proteins'])])
            item['Num Contexts'] = 1
            item['Num Proteins'] = len(item['Proteins'])
            item[score_key] = float(item[score_key])

            item['Decoy Status'] = getDecoyStatus(item['Proteins'])

            del item['Modifications']
        elif float(item[score_key]) == processedDict[dtyper(item['ScanF'])][score_key]:
            p_item = processedDict[dtyper(item['ScanF'])]
            p_item['Context'] += [( (item['Context']), eval(item['Modifications']), eval(item['Proteins']) )]
            p_item['Proteins'] |= set([prot.split(' ')[0] for prot in eval(item['Proteins'])])
            p_item['Num Proteins'] = len(p_item['Proteins'])
            p_item['Num Contexts'] += 1

            p_item['Decoy Status'] = getDecoyStatus(p_item['Proteins'])
        else:
            p_item = processedDict[dtyper(item['ScanF'])]
            if 'Delta %s'%score_key not in p_item:
                p_item['Delta %s'%score_key] = p_item[score_key] - float(item[score_key])

    return processedDict

def indexedTAGGRAPHInfo(scanInfo):
    processedDict = defaultdict(list)

    for item in scanInfo:
        processedDict[item['ScanF']] += [item]

    return processedDict

# Only process items with defined mods
def processTAGGRAPHInfoDefinedOnly(scanInfo, score_key="Composite Score", min_tag_length = 0, min_context_length = 0, de_novo_score_cut = -10000):
    processedDict = {}

    for item in scanInfo:
        if 'Undefined Mass Shift' in item['Match Modifications'] or int(item['Matching Tag Length']) < min_tag_length or len(item['Context']) - 4 < min_context_length or float(item['De Novo Score']) < de_novo_score_cut:
            continue

        if item['ScanF'] not in processedDict:
            processedDict[item['ScanF']] = item
            item['Context'] = [( (item['Context']), eval(item['Match Modifications']), eval(item['Proteins']) )]
            item['Proteins'] = set([prot.split(' ')[0] for prot in eval(item['Proteins'])])
            item['Num Contexts'] = 1
            item['Num Proteins'] = len(item['Proteins'])
            item[score_key] = float(item[score_key])

            item['Decoy Status'] = getDecoyStatus(item['Proteins'])

            del item['Match Modifications']
        elif float(item[score_key]) == processedDict[item['ScanF']][score_key]:
            p_item = processedDict[item['ScanF']]
            p_item['Context'] += [( (item['Context']), eval(item['Match Modifications']), eval(item['Proteins']) )]
            p_item['Proteins'] |= set([prot.split(' ')[0] for prot in eval(item['Proteins'])])
            p_item['Num Proteins'] = len(p_item['Proteins'])
            p_item['Num Contexts'] += 1

            p_item['Decoy Status'] = getDecoyStatus(p_item['Proteins'])
        else:
            p_item = processedDict[item['ScanF']]
            if 'Delta %s'%score_key not in p_item:
                p_item['Delta %s'%score_key] = p_item[score_key] - float(item[score_key])

    return processedDict


# Given processed tag_graph results, this function returns the unique peptides and a representative element (the top scoring PSM for that peptide)
# Note that a single 'PSM' might yield multiple unique peptides if there are different peptides in the returned contexts with the same score
# unique_peptides: mapping of unique peptides and a representative element (the top scoring PSM for that peptide) and list of scans corresponding to this peptide
#
# TODO: Add scan_to_peptides_mapping? As mentioned above, scans might map to multiple, distinct, equally scoring protiens. Do we need to make sure that a spectrum isn't
# used 'twice' in that it generates two different protein mappings? Note that this case is rare and will become much rarer when we implement accurate PSM scoring
def getUniquePeptides(proc_tag_graph, score_key='De Novo Score'):
    unique_peptides = {}

    for scanF in proc_tag_graph:
        matching_pepts = defaultdict(list)
        for context in proc_tag_graph[scanF]['Context']:
            matching_pepts[context[0][2:-2]] += [context]

        if len(matching_pepts) > 1:
            print('Multiple matching peptides found for de novo peptide %s: %s'%(proc_tag_graph[scanF]['De Novo Peptide'], matching_pepts.keys()))

        for peptide in matching_pepts:
            element = copy.copy(proc_tag_graph[scanF])
            element['Context'] = matching_pepts[peptide]
            if peptide not in unique_peptides:
                unique_peptides[peptide] = {'representative': element, 'scans': [scanF]}
            else:
                unique_peptides[peptide]['scans'] += [scanF]
                if float(unique_peptides[peptide]['representative'][score_key]) < proc_tag_graph[scanF][score_key]:
                    unique_peptides[peptide]['representative'] = element

    return unique_peptides

def getUniquePeptidesFromEMResults(em_results_file, prob_cut = 0.5, default_prob = 0.0000000000001, min_pept_length = 0):
    cols, em_data = getScanInfo(em_results_file, delimiter='\t')

    unique_peptides = {}
    for item in em_data:
        context = item['Context'][2:-2]
        if len(context) < min_pept_length:
            continue

        if context not in unique_peptides:
            contexts_dict = defaultdict(lambda: {'mods': [], 'proteins': []})
            contexts_dict[item['Context']]  = {'mods': [eval(item['Mods'])], 'proteins': eval(item['Proteins'])}
            unique_peptides[context] = {'representative': item, 'scans': set([item['ScanF']]), 'contexts': contexts_dict}
            if item['1-lg10 EM'] == 'None':
                item['EM Probability'] = default_prob
                item['1-lg10 EM'] = default_prob
        else:
            unique_peptides[context]['scans'].add(item['ScanF'])
            unique_peptides[context]['contexts'][item['Context']]['mods'] += [eval(item['Mods'])]
            # Note, this overwrites the existing value for the proteins column if this context was already found for this unique peptide, but this should be the same if the peptide is the same anyway
            unique_peptides[context]['contexts'][item['Context']]['proteins'] = eval(item['Proteins'])
            if item['1-lg10 EM'] != 'None' and float(item['1-lg10 EM']) > float(unique_peptides[context]['representative']['1-lg10 EM']):
                unique_peptides[context]['representative'] = item

    # Do this after because not all elements for a given scan will have an associated score (if outputing all exact matches instead of top matches)
    if prob_cut != None:
        for peptide in unique_peptides.keys():
            unique_peptides[peptide]['scans'] = list(unique_peptides[peptide]['scans'])
            if float(unique_peptides[peptide]['representative']['EM Probability']) < prob_cut:
                del unique_peptides[peptide]

    return unique_peptides

def getDecoyStatus(prots):
    if any([prot[0] == '#' for prot in prots]):
        if all([prot[0] == '#' for prot in prots]):
            return 'Yes'
        else:
            return 'Both'
    else:
        return 'No'

def filterForExactsAndIsos(proc_tag_graph):
    for scanF in proc_tag_graph.keys():
        contexts = [context for context in proc_tag_graph[scanF]['Context'] if isIsobaricContext(context)]
        if len(contexts) > 0:
            proc_tag_graph[scanF]['Context'] = contexts
        else:
            del proc_tag_graph[scanF]

def isIsobaricContext(context):
    return len(context[1]) == 0 or all([item[0][0] == 'Isobaric Substitution' for item in context[1]])

def preprocessDatabaseScanInfo(scanInfo, seqMap, fieldMap, scanNumKey='ScanF', seqDelimLen=2):
    DBScanInfo = {}
    for scan in scanInfo:
        #print(scan)
        scanNum = int(scan[scanNumKey])
        info = {}

        if 'Ambiguous Edges' in scan:
            scan['Ambiguous Edges'] = eval(scan['Ambiguous Edges'])
            ambigEdges = scan['Ambiguous Edges']
        else:
            ambigEdges = None

        try:
            if seqDelimLen > 0:
                procSeq = Constants.preprocessSequence(scan['Peptide'][seqDelimLen:-seqDelimLen], seqMap, ambigEdges = ambigEdges)
            else:
                procSeq = Constants.preprocessSequence(scan['Peptide'], seqMap, ambigEdges = ambigEdges)
        #except KeyError:
        #    print('ERROR: Could not process sequence %s with current sequence map' % (scan['Peptide'][seqDelimLen:-seqDelimLen],))
        #    continue
        except IndexError:
            if scan['Peptide'] == 'NULL':
                continue
            else:
                raise IndexError

        info['Peptide'] = procSeq
        for field in fieldMap.keys():
            info[fieldMap[field]] = scan[field]

        addToInfo = True
        try:
            AAs = Constants.getAllAAs(procSeq, ambigEdges=ambigEdges)
        except KeyError:
            raise KeyError

        DBScanInfo[scanNum] = info

    return DBScanInfo

def getScanFDict(dtaList):
    scanFDict = {}
    for dta in dtaList:
        scanF = getScanNum(dta)
        precMass, charge = getPrecMassAndCharge(dta)
        scanFDict[scanF] = {'dta': dta, 'precMass': precMass, 'charge': charge}

    return scanFDict



"""
def plotSpectra(massIntPairs, precMass = 'None', charge = 'Undefined'):
    masses = massIntPairs[:,0]
    intensities = massIntPairs[:,1]
    minMass = masses.min()
    maxMass = masses.max()
    maxInt = intensities.max()

    for i in range(len(masses)):
        pl.axvline(x = masses[i], ymin = 0, ymax = intensities[i]/(maxInt * 1.1))

    pl.ylim([0, maxInt * 1.1])
    pl.xlabel( 'M/Z, NumSpectra: ' + str(np.size(masses)) )
    pl.ylabel('Intensity')
    pl.title('Precursor Mass: ' + str(precMass) + ', Charge: ' + str(charge))
    pl.xlim([minMass * 0.8,maxMass * 1.2])
    pl.show()
"""
if __name__ == '__main__':
    #print(getScanInfo('adevabhaktuni_1310166306.csv'))
    """
    dirPath = 'C:\\Users\\Arun\\Pythonprojects\\DeNovoSequencing\\src\\SpectraCorrelation\\LF2_short_HCD+CID_ath001862_244\\'
    names = getDTAFNamesInDir(dirPath)
    fname = 'C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3383.3383.1.dta'
    massIntPairs = getMassIntPairs(fname)
    precMass, charge = getPrecMassAndCharge(fname)
    plotSpectra(massIntPairs, precMass, charge)
    """
    #fields = ['LIGHT SCAN', 'HEAVY SCAN', 'SEQ', 'SCORE', 'AMBIGUOUS EDGES', 'M+H', 'EPSILON', 'SHARED PEAKS RATIO', 'SEQUENCING TIME']
    #print(getScanInfo('C:\\Users\\Arun\\Proteomics Files\\ath001862UPen10KPen15LRRestrictTest.tdv', fields, delimiter='\t'))
    paramsDict = parseParams('./Misc/LADS_SILAC_Trypsin.ini')
    print(paramsDict)
    seqMap = generateSeqMap(['SEQUEST', 'MASCOT'], paramsDict)
    print(seqMap)


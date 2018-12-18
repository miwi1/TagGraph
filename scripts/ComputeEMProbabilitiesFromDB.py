
'''
Created on Sep 29, 2011

@author: Arun
'''
#slin 201707  added lib path and others.

import os
import sys
import time

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)
import ArgLib
import DataFile
import Validator
import StringFoldingWrapper as SFW

sys.path.insert(2, os.path.abspath(os.path.join(PAR_DIR, 'database')))

from Models import Experiment, Result, Fraction, write_top_results

from sqlalchemy import create_engine, desc
from sqlalchemy.sql import bindparam, select, and_

from collections import defaultdict

experiments = Experiment.__table__
results = Result.__table__
fractions = Fraction.__table__

#                         0   1              2    3        4            5     6         7                    8
# key is scanF, value is id, spectrum_score, ppm, context, mod_context, mods, proteins, matching_tag_length, charge
def getIndexedResultsFromDB(connection, fracs, motifs):
    indexed_results = defaultdict(list)
    for fraction_id, fraction_name in fracs:
        stmt = select([results.c.scan, results.c.id, results.c.spectrum_score, results.c.ppm, results.c.context, results.c.mod_context, results.c.mods, results.c.proteins, results.c.matching_tag_length, results.c.charge]).where(results.c.fraction_id == fraction_id).order_by(results.c.scan).order_by(desc(results.c.composite_score))
        response = connection.execution_options(stream_results=True).execute(stmt)

        for row in SFW.string_folding_wrapper(response):
            indexed_results[fraction_name + ':' + str(row[0])] += [Validator.precalculateEMAttributesForItem(row[1:], motifs)]

    return indexed_results

def writeEMProbabilitiesAndData(connection, indexed_results, fracs, top_only=True, score_cut = 0.99):

    connection.execute(results.update(values={results.c.top_result: False}).where(results.c.fraction_id.in_([frac[0] for frac in fracs])))
                       

    print 'Building update statement to write EM probs'
    stmt = results.update().where(results.c.id == bindparam('b_id')).values({results.c.em_probability: bindparam('em'), results.c.log_em_probability: bindparam('log_em'), results.c.unique_sibling_peptides: bindparam('sib'), results.c.context_mod_variants: bindparam('con'), results.c.num_mod_occurrences: bindparam('occ'), results.c.top_result: bindparam('top')})
    update_data = []
    for scanF in indexed_results:
        top_item = indexed_results[scanF][0]
        update_data += [{'b_id': top_item[0][0], 'em': top_item[2][0], 'log_em': top_item[2][1], 'sib': top_item[1][0], 'con': top_item[1][1], 'occ': top_item[1][2], 'top': True}]

        for item in indexed_results[scanF][1:]:
            if not top_only:
                update_data += [{'b_id': item[0][0], 'em': item[2][0], 'log_em': item[2][1], 'sib': item[1][0], 'con': item[1][1], 'occ': item[1][2], 'top': False}]
            elif top_item[2][0] > score_cut and round(item[0][1], 1) == round(top_item[0][1], 1) and item[0][2] == top_item[0][2] and item[0][4] == top_item[0][4]:
                # Add new item if score is equal to that of top scoring item (i.e., if same localization/mods, but two or more valid sites on protein, etc.)
                update_data += [{'b_id': item[0][0], 'em': top_item[2][0], 'log_em': top_item[2][1], 'sib': item[1][0], 'con': item[1][1], 'occ': item[1][2], 'top': True}]
                

    print 'Executing update'
    connection.execute(stmt, update_data)
    
if __name__ == '__main__':
    print 'modmaxcounts argument is number of iterations in initial EM over all results, maxcounts argument indicates maximum number of iterations before expectation-maximization terminates after reranking is completed (on the top ranked results only). Set fracs to \"all\" to run over all fractions for experiment or supply desired fracs to run EM over separated by commas'
    options = ArgLib.parse(['init', 'sqlitedb', 'experimentname', 'fracs', 'maxcounts', 'modmaxcounts', 'output'])

    t1 = time.time()
    paramsDict = ArgLib.parseInitFile(options.init, options)
    outBase = os.path.splitext(options.output)[0]

    engine = create_engine('sqlite:///' + options.sqlitedb, echo=True)
    conn = engine.connect()

    conn.execute("PRAGMA max_page_count = max_page;");
    conn.execute("PRAGMA temp_store = 2;")
    conn.execute("PRAGMA page_size")

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

    print 'Fetching results from DB for experiment %s, fractions %s'%( str((experiment_id, options.experimentname)) , str(fracs) )
    indexed_results = getIndexedResultsFromDB(conn, fracs, paramsDict['Enzyme']['specificity'])
    
    print 'Results fetched, performing EM'
    initial_spec_match, initial_db_match = Validator.createInitialGuess(pos_mean=(5, 12), pos_sigma=(10, -5, -5, 100), neg_mean=(-2, 12), neg_sigma=(10, -5, -5, 100), match_lengths=range(1,40))
    #indexed_taggraph_info = DataFile.indexedTAGGRAPHInfo(DataFile.getScanInfo(options.taggraph, delimiter='\t')[1])
    spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, ppm_error_models, prior_probs = Validator.performEM(indexed_results, initial_spec_match, initial_db_match, paramsDict, max_iter_all=options.modmaxcounts, max_iter_after_rerank=options.maxcounts, outBase = outBase)

    t2 = time.time()
    writeEMProbabilitiesAndData(conn, indexed_results, fracs)
    print 'Database update completed. Time taken for DB update: %f'%(time.time() - t2,)

    print 'Writing Top Results to file.'
    write_top_results(conn, options.experimentname, os.path.join(os.path.dirname(options.output), '%s_TopResults.tdv'%options.experimentname))

    print 'Total Time taken for EM Step: %f'%(time.time() - t1,)
    #cols = ['ScanF', 'Alignment Score', 'Spectrum Probability Score', 'Composite Score', 'EM Probability', '1-lg10 EM', 'Rank', 'Context', 'Mod Context', 'Mod Ambig Edges', 'Match Modifications', 'Unique Sibling Peptides', 'Context Mod Variants', 'Num Mod Occurrences', 'Modifications', 'Proteins', 'Matching Tag Length', 'De Novo Peptide', 'Unmod De Novo Peptide', 'De Novo Score', 'Time Taken']
    #outFile = open(options.output, 'w')


    # CODE TO WRITE TOP RESULTS HERE
    #outFile = open(outBase + '_TopOnly.tdv', 'w')
    #outFile.write('\t'.join(cols) + '\n')
    #for scanF in indexed_taggraph_info:
    #    item = indexed_taggraph_info[scanF][0]
    #    outFile.write('\t'.join([str(item[col]) for col in cols]) + '\n')
    #outFile.close()

    #Validator.writeModels(outBase + '_MODELS.log', indexed_results)
    #Validator.writeEMProbabilities(outBase + '_EMProbs.tdv', indexed_results, spectrum_match_models, mod_models, context_models, db_match_models, protein_count_models, prior_probs)

    #session.close()
    conn.close()

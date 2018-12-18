import os
import sys
import glob

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

import DataFile
import ArgLib

from Models import Experiment, Result, Fraction, delete_experiment

#from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from sqlalchemy.sql import select, and_


experiment = Experiment.__table__
result = Result.__table__
fraction = Fraction.__table__

def importTAGGRAPHResults(connection, experiment_name, fraction_name, taggraph_files, max_batch_size = 500000):

    try:
        experiment_id = connection.execute(select([experiment.c.id]).where(experiment.c.name == experiment_name)).fetchone()[0]
    except TypeError:
        raise ValueError("ERROR: No experiment by name %s"%experiment_name)

    try:
        fraction_id = connection.execute(select([fraction.c.id]).where(and_(fraction.c.name == fraction_name, fraction.c.experiment_id == experiment_id))).fetchone()[0]
        print "Using existing fraction %s in database for experiment"%str((fraction_id, fraction_name))
    except TypeError:
        print 'FRAC NOT FOUND, CREATING NEW FRAC', fraction_name, experiment_id
        res = connection.execute(fraction.insert().values(name=fraction_name, experiment_id=experiment_id))
        fraction_id = res.inserted_primary_key[0]

    connection.execute(fraction.update().where(fraction.c.id=='fraction_id').values(taggraph_file=str(taggraph_files)))

    values = []
    for taggraph_file in taggraph_files:
        taggraph_info = DataFile.getScanInfoIterator(taggraph_file, delimiter='\t')
        for item in taggraph_info:
            values += [{
                "scan": item['ScanF'],
                "charge": item['Charge'],
                "obs_mh": item['Obs M+H'],
                "theo_mh": item['Theo M+H'],
                "ppm": item['PPM'],
                "retention_time": item['RT'],
                "alignment_score": item['Alignment Score'],
                "spectrum_score": item['Spectrum Probability Score'],
                "composite_score": item['Composite Score'],
                "context": item['Context'],
                "mod_context": item['Mod Context'],
                "mods": item['Match Modifications'],
                "mod_ranges": item['Mod Ranges'],
                "mod_ambig_edges": item['Mod Ambig Edges'],
                "proteins": item['Proteins'],
                "matching_tag_length": item['Matching Tag Length'],
                "time_taken": item['Time Taken'],
                "de_novo_peptide": item['De Novo Peptide'],
                "unmod_de_novo_peptide": item['Unmod De Novo Peptide'],
                "de_novo_score": item['De Novo Score'],
                "num_matches": item['Num Matches'],
                "fraction_id": fraction_id
                
                }]

            if len(values) > max_batch_size:
                res = connection.execute(result.insert(), values)
                values = []
        
        #fraction.results.extend([Result(scan=item['ScanF'], alignment_score=item['Alignment Score'], spectrum_score=item['Spectrum Probability Score'], composite_score=item['Composite Score'], context=item['Context'], mod_context=item['Mod Context'], mods=item['Match Modifications'], mod_ambig_edges=item['Mod Ambig Edges'], proteins=item['Proteins'], matching_tag_length=item['Matching Tag Length'], time_taken=item['Time Taken'], de_novo_peptide=item['De Novo Peptide'], unmod_de_novo_peptide=item['Unmod De Novo Peptide'], de_novo_score=item['De Novo Score'], num_matches=item['Num Matches'])])
    if len(values) > 0:
        res = connection.execute(result.insert(), values)
    
    return True
    
def importExperiment(connection, experiment_name, fmindex, params_file, prog_params, results_dir, fractions):

    experiment_id = connection.execute(select([experiment.c.id]).where(experiment.c.name == experiment_name)).fetchone()

    # Delete existing experiment of same name if present in database
    if experiment_id:
        delete_experiment(connection, experiment_name)

    res = connection.execute(experiment.insert().values(name=experiment_name, fmindex=fmindex, params=params_file, taggraph_parameters=str(prog_params)))
    experiment_id = res.inserted_primary_key[0]


    values = []
    for fraction_name in fractions:
        values += [{
            "name": fraction_name,
            "experiment_id": experiment_id
            }]
    res = connection.execute(fraction.insert(), values)

    #session.commit()

    for fraction_name in fractions:
        results_files = glob.glob(results_dir +  os.path.sep + '*_' + fraction_name + '_*tdv')
        print fraction_name, results_files
        if importTAGGRAPHResults(connection, experiment_name, fraction_name, results_files):
            print 'Fraction %s imported successfully'%fraction_name
        else:
            print 'ERROR: Unable to import results for fraction %s'%fraction_name

    return True



if __name__ == '__main__':
    options = ArgLib.parse(['sqlitedb', 'taggraph', 'init', 'fmindex', 'modtolerance', 'ppmstd', 'maxcounts', 'modmaxcounts', 'experimentname', 'fracs'], optArgs=[{ 'opts': ('-Y', '--type'), 'attrs': {'type': 'string', 'dest': 'type', 'default': None, 'help': 'Value is either experiment or single, defines where import is an experiment import or single sample import'} }])

    print 'If importing just a single fraction, ignore all arguments except sqlitedb, taggraph, type (set as single), experimentname, fracs (set as fraction name)'
    print 'If type of import is experiment, set taggraph argument to directory where results are located, and fracs to a tuple of fractions to import. TAGGRAPH parameters (init, fmindex, modtolerance, ppmstd, maxcounts, modmaxcounts) will be saved in db for record keeping'

    engine = create_engine('sqlite:///' + options.sqlitedb, echo=True)
    #Session = sessionmaker(bind=engine)

    #session = Session()

    conn = engine.connect()

    if options.type == "single":
        importTAGGRAPHResults(conn, options.experimentname, options.fracs, eval(options.taggraph))
    elif options.type == "experiment":

        prog_params = {'ppmstd': options.ppmstd, 'maxcounts': options.maxcounts, 'modmaxcounts': options.modmaxcounts, 'modtolerance': options.modtolerance}
        importExperiment(conn, options.experimentname, options.fmindex, options.init, str(prog_params), options.taggraph, options.fracs.split(','))
    else:
        print 'ERROR: Please specify type as single or experiment'

    conn.close()
    #session.close()

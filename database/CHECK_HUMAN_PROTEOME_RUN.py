#verify, create .ERROR file when encounter erros.

import os
import sys

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

import ArgLib
import DataFile

import glob

from Models import Experiment, Result, Fraction

from sqlalchemy import create_engine, func
from sqlalchemy.sql import select, and_

experiment = Experiment.__table__
fraction = Fraction.__table__
result = Result.__table__

from collections import defaultdict

if __name__ == '__main__':
    print 'output is the name of file to write the report to'
    print 'experimentname is the experiment to check for integrity'
    print 'dtadir is directory of parent project'
    options = ArgLib.parse(['dtadir', 'experimentname', 'output'])
    sqlitedb_loc = os.path.join(options.dtadir, 'results.db')
    engine = create_engine('sqlite:///' + sqlitedb_loc, echo=True)
    conn = engine.connect()
    experiment_name = options.experimentname
    append_string = ''
    try:
        experiment_id = conn.execute(select([experiment.c.id]).where(experiment.c.name == experiment_name)).fetchone()[0]
    except TypeError:
        append_string = '.ERROR'
        experiment_id = 0
    experiment_dir = os.path.join(options.dtadir, experiment_name)
    #Get all the fraction numbers
    peaks_dir = os.path.join(experiment_dir, 'de_novo')
    peaks_files = glob.glob(peaks_dir + os.path.sep + '*_parsed_F*.tdv')
    fracs = []
    for file in peaks_files:
        fracs += [int(file.split('_')[-1][1:-4])]
    fracs = sorted(fracs)
    frac_counts = {}
    for frac in fracs:
        fraction_name = 'F%i'%frac
        try:
            frac_id = conn.execute(select([fraction.c.id]).where(and_(fraction.c.name == fraction_name, fraction.c.experiment_id == experiment_id))).fetchone()[0]
            counts = conn.execute(select([func.count(result.c.id)]).where(result.c.fraction_id == frac_id)).fetchone()[0]
            frac_counts[frac] = counts
            if counts == 0:
                append_string = '.ERROR'
        except TypeError:
            append_string = '.ERROR'
            frac_counts[frac] = 'NOT PRESENT'
    conn.close()
    # Write output summary
    outFile = open(options.output + '.%i'%len(fracs) + append_string, 'w')
    outFile.write(options.dtadir + '\n')
    outFile.write('Experiment Name %s ID %i\n\n'%(experiment_name, experiment_id))
    outFile.write('Result Counts for %i fractions\n'%(len(fracs)))
    for frac in fracs:
        outFile.write('F%i: %s\n'%(frac, str(frac_counts[frac])))
    outFile.close()
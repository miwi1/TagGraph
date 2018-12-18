import os
import sys

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float, Text, Boolean
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy import create_engine, desc

import ArgLib
import Validator
import StringFoldingWrapper as SFW

Base = declarative_base()

from sqlalchemy.sql import select
import sqlalchemy

class Experiment(Base):
    __tablename__ = 'experiment'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    fmindex = Column(Text)
    params = Column(Text)
    taggraph_parameters = Column(Text)
    em_parameters = Column(Text)
    location = Column(Text)
    taggraph_finished = Column(Boolean, default=False)
    em_finished = Column(Boolean, default=False)
    

class Fraction(Base):
    __tablename__ = 'fraction'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    experiment_id = Column(Integer, ForeignKey('experiment.id'), index=True)
    em_finished = Column(Boolean, default=False)
    taggraph_file = Column(Text)
    
    experiment = relationship("Experiment", backref=backref('fractions', order_by=id))


class Result(Base):
    __tablename__ = 'result'
    id = Column(Integer, primary_key=True)
    scan = Column(Integer)
    charge = Column(Integer)
    obs_mh = Column(Float)
    theo_mh = Column(Float)
    retention_time = Column(Float)
    ppm = Column(Float)
    alignment_score = Column(Float)
    spectrum_score = Column(Float)
    composite_score = Column(Float)
    em_probability = Column(Float)
    log_em_probability = Column(Float)
    context = Column(Text)
    mod_context = Column(Text)
    modified = Column(Boolean)
    mods = Column(Text)
    mod_ambig_edges = Column(Text)
    mod_ranges = Column(Text)
    proteins = Column(Text)
    matching_tag_length = Column(Integer)
    time_taken = Column(Float)
    de_novo_peptide = Column(Text)
    unmod_de_novo_peptide = Column(Text)
    de_novo_score = Column(Float)
    num_matches = Column(Integer)
    top_result = Column(Boolean, default=False)
    unique_sibling_peptides = Column(Integer)
    context_mod_variants = Column(Integer)
    num_mod_occurrences = Column(Integer)
    spectral_counts = Column(Integer)
    abundance = Column(Float)
    fraction_id = Column(Integer, ForeignKey('fraction.id'), index=True)

    fraction = relationship("Fraction", backref=backref('results', order_by=id, lazy='dynamic'))

class Protein(Base):
    __tablename__ = 'protein'
    id = Column(Integer, primary_key = True)
    name = Column(Text)
    experiment_id = Column(Integer, ForeignKey('experiment.id'), index=True)
    counts = Column(Float)
    abundance = Column(Float)
    
class ProteinMapping(Base):
    __tablename__ = 'protein_mapping'
    id = Column(Integer, primary_key = True)
    offset = Column(Integer)
    # Note: for cases where there are multiple 'top result' peptides for a given scanF, each one will have a corresponding ProteinMapping Entry (will need to remember to only count UNIQUE result_ids when doing protein count calculations)
    # Indicates whether assigned protein is present in the linked result entry or not. 
    # Also, sometimes there are multiple 'top result' entries which may have different linked proteins. 
    # Sometimes most parsimonious protein can be from lower ranking entry. This is due to a bug in the EM algorithm as of 1/23/2015, may be fixed at a later date so that this would be 'tied' with top result or score more highly
    # If a peptide matches in two places on the same protein, it will have two corresponding (result_id, protein_id) data points with different offsets
    from_result = Column(Boolean)
    protein_id = Column(Integer, ForeignKey('protein.id'))
    result_id = Column(Integer, ForeignKey('result.id'))

class Modification(Base):
    __tablename__ = 'modification'
    id = Column(Integer, primary_key = True)
    # Mod format: ( (mod_name, mod_mass, localization), .... )
    # Will likely grab using LIKE query
    mod = Column(Text)
    mod_length = Column(Integer)
    accurate_mass_mod = Column(Text)
    average_ppm_error = Column(Float)
    average_ppm_error_mappings = Column(Float)
    experiment_id = Column(Integer, ForeignKey('experiment.id'), index=True)
    over_cutoff = Column(Integer)
    over_cutoff_single_counts = Column(Text)
    over_cutoff_mappings = Column(Integer)

class ModMapping(Base):
    __tablename__ = 'mod_mapping'
    id = Column(Integer, primary_key = True)
    spec_score = Column(Float)
    spec_score_diff = Column(Float)
    ppm = Column(Float)
    enumerated = Column(Boolean)
    alternate_explanation = Column(Boolean)
    
    result_id = Column(Integer, ForeignKey('result.id'), index=True)
    modification_id = Column(Integer, ForeignKey('modification.id'), index=True)

result = Result.__table__
experiment = Experiment.__table__
fraction = Fraction.__table__
protein = Protein.__table__
protein_mapping = ProteinMapping.__table__
modification = Modification.__table__
mod_mapping = ModMapping.__table__

def write_top_results(connection, experiment_name, out_file_name, with_same_pept_other_context = False):
    try:
        experiment_id = connection.execute(select([experiment.c.id]).where(experiment.c.name == experiment_name)).fetchone()[0]
    except TypeError:
        print 'No experiment by name %s present in database'%experiment_name

    cols = ['ScanF', 'Charge', 'Retention Time', 'Obs M+H', 'Theo M+H', 'PPM', 'EM Probability', '1-lg10 EM', 'Spectrum Score', 'Alignment Score', 'Composite Score', 'Unique Siblings', 'Context Mod Variants', 'Num Mod Occurrences', 'Context', 'Mod Context', 'Mods', 'Mod Ambig Edges', 'Mod Ranges', 'Proteins', 'De Novo Peptide', 'De Novo Score', 'Matching Tag Length', 'Num Matches']
    out_file = open(out_file_name, 'w')
    out_file.write('\t'.join(cols) + '\n')
    
    fracs = connection.execute(select([fraction.c.name, fraction.c.id]).where(fraction.c.experiment_id == experiment_id)).fetchall()
    for frac_name, frac_id in fracs:
        stmt = select([result.c.scan, result.c.charge, result.c.retention_time, result.c.obs_mh, result.c.theo_mh, result.c.ppm, result.c.em_probability, result.c.log_em_probability, result.c.spectrum_score, result.c.alignment_score, result.c.composite_score, result.c.unique_sibling_peptides, result.c.context_mod_variants, result.c.num_mod_occurrences, result.c.context, result.c.mod_context, result.c.mods, result.c.mod_ambig_edges, result.c.mod_ranges, result.c.proteins, result.c.de_novo_peptide, result.c.de_novo_score, result.c.matching_tag_length, result.c.num_matches, result.c.top_result]).where(result.c.fraction_id == frac_id).order_by(result.c.scan).order_by(desc(result.c.top_result))
        response = connection.execution_options(stream_results=True).execute(stmt)

        top_results = {}
        for row in response:
            if row[24]:
                write_info = (frac_name + ':' + str(row[0]),) + row[1:]
                out_file.write('\t'.join([str(write_info[i]) for i in range(len(cols))]) + '\n')
                top_results[row[0]] = row
            elif with_same_pept_other_context and row[14][2:-2] == top_results[row[0]][14][2:-2] and row[16] == top_results[row[0]][16]:
                write_info = (frac_name + ':' + str(row[0]),) + row[1:]
                out_file.write('\t'.join([str(write_info[i]) for i in range(len(cols))]) + '\n')

    out_file.close()

# set attributes to list of desired columns if you want to select specific columns from DB (default is all)
def fetch_top_results(connection, experiment_id, attributes = None):
    frac_ids = [id[0] for id in connection.execute(select([fraction.c.id]).where(fraction.c.experiment_id == experiment_id)).fetchall()]

    if attributes == None:
        return SFW.string_folding_wrapper(connection.execution_options(stream_results=True).execute(select([result]).where(result.c.fraction_id.in_(frac_ids)).where(result.c.top_result == True)))
    else:
        return SFW.string_folding_wrapper(connection.execution_options(stream_results=True).execute(select(attributes).where(result.c.fraction_id.in_(frac_ids)).where(result.c.top_result == True)))

def fetch_proteins(connection, experiment_id):
    return connection.execute(select([protein.c.id, protein.c.name, protein.c.counts, protein.c.abundance]).where(protein.c.experiment_id == experiment_id))

def fetch_modifications(connection, experiment_id, attributes = None):
    if attributes == None:
        cols = [modification]
    else:
        cols = attributes

    return connection.execute(select(cols).where(modification.c.experiment_id == experiment_id))

def fetch_mod_by_name(connection, experiment_id, name):
    return connection.execute(select(modification).where(modification.c.mod.like('%' + name + '%')))

def fetch_mod_mappings(connection, modification_id, em_prob_cutoff = 0.99):
    return connection.execute(select([result.c.id, mod_mapping.c.id, result.c.fraction_id, result.c.scan, result.c.mods, result.c.context, result.c.mod_context, result.c.em_probability, mod_mapping.c.spec_score, mod_mapping.c.spec_score_diff, mod_mapping.c.ppm, mod_mapping.c.alternate_explanation, mod_mapping.c.enumerated]).where(mod_mapping.c.modification_id == modification_id).where(result.c.id == mod_mapping.c.result_id).where(result.c.em_probability >= em_prob_cutoff)).fetchall()

def fetch_results_for_fraction(connection, fraction_id, top_only = True, attributes = None):
    if attributes == None:
        cols = [result]
    else:
        cols = attributes

    if top_only:
        return SFW.string_folding_wrapper(connection.execution_options(stream_results=True).execute(select([cols]).where(result.c.fraction_id == fraction_id).where(result.c.top_result == True)))
    else:
        return SFW.string_folding_wrapper(connection.execution_options(stream_results=True).execute(select(cols).where(result.c.fraction_id == fraction_id)))

def fetch_all_proteins_and_results_with_scores(connection, experiment_id, em_cut = 0.99):
    return connection.execute(select([result.c.fraction_id, result.c.scan, result.c.context, result.c.mod_context, result.c.mods, result.c.em_probability, protein_mapping.c.offset, protein.c.id, protein.c.name, result.c.ppm, result.c.spectrum_score]).where(result.c.id == protein_mapping.c.result_id).where(result.c.em_probability >= em_cut).where(protein_mapping.c.protein_id == protein.c.id).where(protein.c.experiment_id == experiment_id)).fetchall()

def fetch_all_proteins_and_results(connection, experiment_id, em_cut = 0.99):
    return connection.execute(select([result.c.fraction_id, result.c.scan, result.c.context, result.c.mod_context, result.c.mods, result.c.em_probability, protein_mapping.c.offset, protein.c.id, protein.c.name]).where(result.c.id == protein_mapping.c.result_id).where(result.c.em_probability >= em_cut).where(protein_mapping.c.protein_id == protein.c.id).where(protein.c.experiment_id == experiment_id)).fetchall()

def write_exact_matches(connection, experiment_name, out_file_name):
    try:
        experiment_id = connection.execute(select([experiment.c.id]).where(experiment.c.name == experiment_name)).fetchone()[0]
    except TypeError:
        print 'No experiment by name %s present in database'%experiment_name

    cols = ['ScanF', 'EM Probability', '1-lg10 EM', 'Spectrum Score', 'Alignment Score', 'Composite Score', 'Unique Siblings', 'Context Mod Variants', 'Num Mod Occurrences', 'Context', 'Mod Context', 'Mods', 'Mod Ambig Edges', 'Mod Ranges', 'Proteins', 'De Novo Peptide', 'De Novo Score', 'Matching Tag Length', 'Num Matches']
    out_file = open(out_file_name, 'w')
    out_file.write('\t'.join(cols) + '\n')
    
    fracs = connection.execute(select([fraction.c.name, fraction.c.id]).where(fraction.c.experiment_id == experiment_id)).fetchall()
    for frac_name, frac_id in fracs:
        stmt = select([result.c.scan, result.c.em_probability, result.c.log_em_probability, result.c.spectrum_score, result.c.alignment_score, result.c.composite_score, result.c.unique_sibling_peptides, result.c.context_mod_variants, result.c.num_mod_occurrences, result.c.context, result.c.mod_context, result.c.mods, result.c.mod_ambig_edges, result.c.mod_ranges, result.c.proteins, result.c.de_novo_peptide, result.c.de_novo_score, result.c.matching_tag_length, result.c.num_matches]).where(result.c.fraction_id == frac_id).order_by(result.c.scan)
        response = connection.execution_options(stream_results=True).execute(stmt)

        for row in response:
            if len([mod for mod in eval(row[11]) if mod[0][0] != 'Isobaric Substitution']) == 0:
                write_info = (frac_name + ':' + str(row[0]),) + row[1:]
                out_file.write('\t'.join([str(write_info[i]) for i in range(len(cols))]) + '\n')

    out_file.close()

def fetch_all_experiments(connection):
    return connection.execute(select([experiment.c.id, experiment.c.name])).fetchall()

def fetch_all_fractions(connection, experiment_id):
    return connection.execute(select([fraction.c.id, fraction.c.name]).where(fraction.c.experiment_id == experiment_id)).fetchall()

def get_experiment_id_from_name(connection, experiment_name):
    try:
        return connection.execute(select([experiment.c.id]).where(experiment.c.name == experiment_name)).fetchone()[0]
    except TypeError:
        raise ValueError("ERROR: No experiment by name %s"%experiment_name)

def get_peptides_for_protein(connection, protein_id, em_prob_cutoff = 0.99):
    return connection.execute(select([result.c.id, protein_mapping.c.id, result.c.fraction_id, result.c.scan, result.c.proteins, result.c.context, result.c.mod_context, result.c.mods, result.c.em_probability, protein_mapping.c.offset]).where(protein_mapping.c.protein_id == protein_id).where(result.c.id == protein_mapping.c.result_id).where(result.c.em_probability >= em_prob_cutoff)).fetchall()

# Note: this method does not recalculate the composite score correctly, just assumes that these scans are crap and makes sure that they are ranked at the bottom
# THIS HAPPENS VERY RARELY, MAYBE 1 out every million spectra analyzed
# New versions of TAG_GRAPH_PROB_SCORE.py will detect this case, set the score to -10, and output a CRITICAL ERROR for log purposes, so this method should no longer be necessary
def reset_scores_where_inf(connection, new_spectrum_score = -10, new_composite_score = -10):
    connection.execute(result.update(values={result.c.spectrum_score: new_spectrum_score, result.c.composite_score: new_composite_score}).where(result.c.spectrum_score == float("inf")))

def delete_proteins(connection):
    res = connection.execute(protein.delete())
    print 'Deleted %i protein entries'%res.rowcount

def delete_protein_mappings(connection):
    res = connection.execute(protein_mapping.delete())
    print 'Deleted %i protein_mapping entries'%res.rowcount

def delete_modifications(connection):
    res = connection.execute(modification.delete())
    print 'Deleted %i modification entries'%res.rowcount

def delete_mod_mappings(connection):
    res = connection.execute(mod_mapping.delete())
    print 'Deleted %i mod_mapping entries'%res.rowcount
    
def delete_experiment(connection, experiment_name):
    try:
        experiment_id = connection.execute(select([experiment.c.id]).where(experiment.c.name == experiment_name)).fetchone()[0]
    except TypeError:
        print 'No experiment by name %s present in database'%experiment_name

    frac_ids = [id[0] for id in connection.execute(select([fraction.c.id]).where(fraction.c.experiment_id == experiment_id)).fetchall()]

    # delete results
    res = connection.execute(result.delete().where(result.c.fraction_id.in_(frac_ids)))
    print 'Deleted %i results'%res.rowcount

    # delete fractions
    res = connection.execute(fraction.delete().where(fraction.c.experiment_id == experiment_id))
    print 'Deleted %i fractions'%res.rowcount

    # delete experiment
    res = connection.execute(experiment.delete().where(experiment.c.id == experiment_id))
    print 'Deleted %i experiments'%res.rowcount

def add_column(engine, table, column):
    table_name = table.description
    column_name = column.compile(dialect=engine.dialect)
    column_type = column.type.compile(engine.dialect)
    try:
        engine.execute('ALTER TABLE %s ADD COLUMN %s %s' % (table_name, column_name, column_type))
    except sqlalchemy.exc.OperationalError:
        print 'CRITICAL ERROR: column %s already exists in table'%column_name

def get_tissue_experiment_map(tissue_paths):
    tissue_map = defaultdict(list)
    for experiment_dir in tissue_paths:
        db_loc = os.path.join(experiment_dir, 'results.db')
        engine = create_engine('sqlite:///' + db_loc, echo=True)
        conn = engine.connect()

        experiments = Models.fetch_all_experiments(conn)
        for experiment_id, experiment_name in experiments:
            tissue_map[os.path.basename(os.path.normpath(path))] += [experiment_name]

        conn.close()
        
    return tissue_map
       
if __name__ == '__main__':
    print 'Creates DB from defined schema'
    options = ArgLib.parse(['sqlitedb'])

    engine = create_engine('sqlite:///' + options.sqlitedb, echo=True)
    Base.metadata.create_all(engine)
    

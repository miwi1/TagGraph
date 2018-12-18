#slin 201707  added lib path and others.

import os
import sys


PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

DATABASE_SCRIPT_DIR = os.path.join(PAR_DIR, 'database')
PREPROCESSOR_SCRIPT_DIR = os.path.join(PAR_DIR, 'preprocessors')
RESOURCES_DIR = os.path.join(PAR_DIR, 'resources')

import pepInput
import shutil
import ArgLib
import DataFile
import pepInput
import glob
import pickle

if __name__ == '__main__':
    print 'dtadir is the directory containing the mzXML files to analyze'
    print 'peaks is a dictionary mapping {experiment_name: peaks csv}'
    print 'output is the directory to move all files to and set up the project in'
    options = ArgLib.parse(['init', 'dtadir', 'peaks', 'output'])

    print 'options.output: %s' % (options.output)
    print 'normpath(options.output): %s' % (os.path.normpath(options.output))
    # Fails with an OSError if directory already exists
    os.makedirs(options.output)

    # Create database
    args = ['--sqlite', os.path.join(options.output, 'results.db')]
    print 'Models.py dir: %s' % (DATABASE_SCRIPT_DIR)
    DataFile.executeProcess(DATABASE_SCRIPT_DIR, 'Models.py', args)

    # Make experiment directories
    # Structure
    # /options.output
    # .../ExperimentName
    # ...../data
    # ...../de_novo
    # ...../taggraph
    for experiment, peaks_file in eval(options.peaks).items():
        experiment_dir = os.path.join(options.output, experiment)
        os.makedirs(experiment_dir)

        # Make the de_novo subdirectory
        peaks_dir = os.path.join(experiment_dir, 'de_novo')
        os.makedirs(peaks_dir)

        # Copy the peaks file to the de_novo subdirectory
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        shutil.copy(peaks_file, peaks_dir + os.path.sep + '%s_PEAKS.csv'%experiment)
        '''
        new_peaks_file = peaks_dir + '/' + experiment + '_PEAKS.csv'

        # Move the pickled fileFractionMapping file from the symLinkDir to the output directory
        shutil.move(options.dtadir + '/fileFractionMapping.pck', options.output + '/')

        if peaks_file.upper().endswith('.XML') or peaks_file.upper().endswith('.PEPXML'):
            # Copy the original pepXML to the de_novo output directory
            shutil.copy(peaks_file, peaks_dir + "/")
            mappingFile = open(options.output + '/fileFractionMapping.pck', 'rb')
            fileFractionMapping = pickle.load(mappingFile)
            pepInput.convertPepxmlToCSV(peaks_file, new_peaks_file, fileFractionMapping)
        else:
            # Note: This assumes it will have .CSV extension, if it doesn't have .[pep]XML
            # Might be good to add an additional check
            shutil.copy(peaks_file, new_peaks_file)

        # preprocess peaks result file (convert to tab-delimited .TDV file)
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        args = ['--init', options.init, '--symbolmap', os.path.join(RESOURCES_DIR, 'PEAKSSymbMap.txt'), '--peaks', peaks_dir+os.path.sep+'%s_PEAKS.csv'%experiment, '--output', peaks_dir+os.path.sep+'%s_PEAKS_parsed.tdv'%experiment]
        '''
        #args = ['--init', options.init, '--symbolmap', os.path.join(RESOURCES_DIR, 'PEAKSSymbMap.txt'), '--peaks', peaks_dir+'/'+'%s_PEAKS.csv'%experiment, '--output', peaks_dir+'/'+'%s_PEAKS_parsed.tdv'%experiment]
        args = ['--init', options.init, '--symbolmap', os.path.join(RESOURCES_DIR, 'PEAKSSymbMap.txt'), '--peaks', new_peaks_file, '--output', peaks_dir+'/'+'%s_PEAKS_parsed.tdv'%experiment]
        DataFile.executeProcess(PREPROCESSOR_SCRIPT_DIR, 'ParsePEAKS7Results.py', args)

        # Unpack results into fractions
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        DataFile.PEAKS7_split_by_fraction(peaks_dir+os.path.sep+'%s_PEAKS_parsed.tdv'%experiment)
        '''
        DataFile.PEAKS7_split_by_fraction(peaks_dir+'/'+'%s_PEAKS_parsed.tdv'%experiment)

        data_dir = os.path.join(experiment_dir, 'data')
        os.makedirs(data_dir)
        ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
        for mzXML in glob.glob(options.dtadir + os.path.sep + '*%s*mzML'%experiment):
        '''
        for mzXML in glob.glob(options.dtadir + '/' + '*%s*mz*ML'%experiment):
            ''' Replace os.path.sep with '/' to fix Windows backslash issues. --smp
            shutil.copy(mzXML, data_dir+os.path.sep+os.path.basename(mzXML))
            '''
            shutil.copy(mzXML, data_dir+'/'+os.path.basename(mzXML))

        os.makedirs(os.path.join(experiment_dir, 'taggraph'))
            
import os
import sys
import time
import ConfigParser
import getopt
import shutil
import glob
import csv
import pickle
import time, datetime
import tempfile

CUR_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1,CUR_DIR)
LIB_DIR = CUR_DIR+"/lib"
sys.path.insert(1,LIB_DIR)
SCRIPTS_DIR=CUR_DIR+"/scripts"

import ArgLib
import DataFile
import pepInput
import verifyEM

MSGBORDER="---------------------------------------------------------------------------------------"
TAGGRAPH_CONFIG_HEADER =  '\n\n***************************************\n'
TAGGRAPH_CONFIG_HEADER += '*** START TAGGRAPH FROM CONFIG FILE ***\n'
TAGGRAPH_CONFIG_HEADER += '***************************************\n'
TAGGRAPH_CONFIG_FOOTER =  '\n\n*************************************\n'
TAGGRAPH_CONFIG_FOOTER += '*** END TAGGRAPH FROM CONFIG FILE ***\n'
TAGGRAPH_CONFIG_FOOTER += '*************************************\n'

symLinkBaseDir='/tmp/'
ZTestPctToGrab=0.2

def write2FileStdout(fh,message):
    fh.write("%s\n" % message)
    print message

def printUsage():
    print '\nUsage: '
    print '   This HELP:'
    print '      $ python '+sys.argv[0]+' -h'
    print '      $ python '+sys.argv[0]+' --help'
    print '   Run TagGraph:'
    print '      $ python '+sys.argv[0]+' <TagGraph Config File>'

def ConfigSectionMap(Config, section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def isNumeric(value):
    try:
        float(value)
        return True
    except:
        return False

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h", ["help"])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        print 'opt: %s' % (opt)
        print 'arg: %s' % (arg)
        if opt in ("-h", "--help"):
            printUsage()
            sys.exit(1)
    ''' Now process the arguments (INI file path)'''
    if len(args) != 1:
        printUsage()
        sys.exit(1)
    configFileName = args[0]
    ### create a output file/handle:
    tmpFolder=tempfile.gettempdir()
    (tm_year,tm_mon,tm_mday,tm_hour,tm_min,tm_sec,tm_wday,tm_yday,tm_isdst)=time.localtime(time.time())
    runCapture=tmpFolder+'/RunTG'+str(tm_mon)+str(tm_mday)+str(tm_hour)+str(tm_min)+str(tm_sec)+'.txt'
    fh=open(runCapture,'w')
    write2FileStdout(fh,'**** start TagGraph process: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,TAGGRAPH_CONFIG_HEADER)
    write2FileStdout(fh,configFileName+"\n")
    if os.path.isfile(configFileName) and os.access(configFileName, os.R_OK):
        write2FileStdout(fh,MSGBORDER)
        write2FileStdout(fh,"Using Configuration File: %s"%configFileName)
        write2FileStdout(fh,MSGBORDER)
    else:
        #print ' ** FAILURE ** Could not read configuration file: \'%s\'' % (configFileName)
        write2FileStdout(fh,' ** FAILURE ** Could not read configuration file: \'%s\'' % (configFileName))
        sys.exit(1)
    theConfig = ConfigParser.ConfigParser()
    theConfig.optionxform = str
    theConfig.read(configFileName)
    #sectionNames = theConfig.sections()
    generalSectionMap = ConfigSectionMap(theConfig, "General")
    tagGraphSectionMap = ConfigSectionMap(theConfig, "TagGraph")
    ######    init', 'dtadir', 'peaks', 'output', 'ppmstd', 'modtolerance', 'unimoddict', 'maxcounts', 'modmaxcounts', 'fmindex', 'model', 'config'
    ## Define our Required Arguments ##
    fatalError = False
    ## Arguments that must exist, and be numbers ##
    requiredTagGraphNumericArgs = ['ppmstd','modtolerance','maxcounts','modmaxcounts']
    ## Arguments that must exist, and be paths that point to files that exist and are Readable ##
    requiredTagGraphExistingFiles = ['unimoddict','model','config','init','de_novo']
    ## Arguments that must exist, and be directories that can be created on the filesystem ##
    requiredTagGraphToCreate = ['output']
    ## Special Arguments:
    # ExperimentName must be a string
    # d must be a directory, with mzML/mzXML files in it that start with ExperimentName
    # f must be an fmindex name of the form <basepath>.fm, where <basepath> is the basename and the following files should exist: <basename>.fm.1, <basename>.seqnames.1, <basename>.offsets
    ## Arguments that must exist, and be numbers ##
    for currArg in requiredTagGraphNumericArgs:
        if currArg in tagGraphSectionMap:
            if isNumeric(tagGraphSectionMap[currArg]):
                write2FileStdout(fh,'* Found Required Numeric TagGraph Parameter \'%s\'  : \'%s\'' % (currArg, tagGraphSectionMap[currArg]))
            else:
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'%s\' must be a numeric value, found value \'%s\'' % (currArg,tagGraphSectionMap[currArg]))
        else:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'%s\' not found in config file' % (currArg))
    ## Arguments that must exist, and be paths that point to files that exist and are Readable ##
    for currArg in requiredTagGraphExistingFiles:
        if currArg in tagGraphSectionMap:
            if os.path.isfile(tagGraphSectionMap[currArg]) and os.access(tagGraphSectionMap[currArg], os.R_OK):
                write2FileStdout(fh,'* Found Required Readable File for TagGraph Parameter \'%s\' : \'%s\'' % (currArg, tagGraphSectionMap[currArg]))
            else:
                if not os.path.isfile(tagGraphSectionMap[currArg]):
                    fatalError = True
                    write2FileStdout(fh,'** FAILURE ** Could not find file for Required Parameter \'%s\' at \'%s\'' % (currArg, tagGraphSectionMap[currArg]))
                elif not os.access(tagGraphSectionMap[currArg], os.R_OK):
                    fatalError = True
                    write2FileStdout(fh,'** FAILURE ** Could not Read file for Required Parameter \'%s\' at \'%s\' (check permissions)' % (currArg, tagGraphSectionMap[currArg]))
        else:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'%s\' not found in config file' % (currArg))
    ## Arguments that must exist, and be directories that should not already exist but can be created on the filesystem ##
    for currArg in requiredTagGraphToCreate:
        if currArg in tagGraphSectionMap:
            dirToCreate = tagGraphSectionMap[currArg]
            if not os.path.exists(dirToCreate):
                try:
                    ## Should be able to make the directory, and then remove it ##
                    os.makedirs(dirToCreate)
                    os.rmdir(dirToCreate)
                    write2FileStdout(fh,'* Found Required Createable Directory for TagGraph Parameter \'%s\' : \'%s\'' % (currArg, dirToCreate))
                except OSError:
                    fatalError = True
                    write2FileStdout(fh,'** FAILURE ** Unable to Create Directory for Required Parameter \'%s\' at \'%s\'' % (currArg, dirToCreate))
            else:
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** File/Directory for Required Parameter \'%s\' at \'%s\' already exists! Should be created by TagGraph' % (currArg, dirToCreate))
        else:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'%s\' not found in config file' % (currArg))
    ## Now Lets Handle the Special Cases
    ## ExperimentName must be a string
    experimentName = ''
    if not 'ExperimentName' in tagGraphSectionMap:
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'ExperimentName\' not found in config file')
    else:
        experimentName = tagGraphSectionMap['ExperimentName']
        write2FileStdout(fh,'* Found Required TagGraph Parameter ExperimentName: \'%s\'' % (experimentName))
        
    ## New Method: numFractions = 2, fraction01 = <path to file 1>, fraction02 = <path to file 2>
    numFractions = 0
    foundNumFractions = False
    dataDirectory = ''
    symLinkDir = symLinkBaseDir
    if not 'numFractions' in tagGraphSectionMap:
        ## Check for dataDirectory and automatically finding data files from the de novo files
        if not 'dataDirectory' in tagGraphSectionMap:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required Directory TagGraph Parameter \'dataDirectory\' not found in config file')
        else:
            dataDirectory = tagGraphSectionMap['dataDirectory']
            if not dataDirectory.endswith('/'):
                dataDirectory += '/'            
            if ( not (dataDirectory.startswith('/'))):
                levelup=dataDirectory.count('../')
                if (levelup==0):
                    dataDirectory=CUR_DIR+'/'+ dataDirectory
                else:
                    splitDataDir=dataDirectory.split("/")
                    splitCurDir=CUR_DIR.split("/")
                    tmpD=''
                    for i in xrange(0,len(splitCurDir)-levelup):
                        tmpD = tmpD+splitCurDir[i]+"/"
                    for i in xrange(levelup,len(splitDataDir)-1):
                        tmpD = tmpD+splitDataDir[i]+"/"
                    dataDirectory=tmpD
            write2FileStdout(fh,"dataDirectory: %s" % dataDirectory)
            if not os.path.exists(dataDirectory):
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** Required Directory TagGraph Parameter \'dataDirectory\' does not exist at: \'%s\'' % (dataDirectory))
            elif not os.path.isdir(dataDirectory):
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** Required Directory TagGraph Parameter \'dataDirectory\' does not point to a directory at: \'%s\'' % (dataDirectory))
            else:
                ## We need to get the data file names from the de novo file, and check for them in the dataDirectory
                fileFractionMapping = []
                deNovoFile = tagGraphSectionMap['de_novo']
                if deNovoFile.upper().endswith('.XML') or deNovoFile.upper().endswith('.PEPXML') or deNovoFile.upper().endswith('.CSV'):
                    if deNovoFile.upper().endswith('.XML') or deNovoFile.upper().endswith('.PEPXML'):
                        fileFractionMapping = pepInput.getFileFractionMappingFromPepXML(deNovoFile)
                    else: ## deNovoFile.upper().endswith('.CSV'):
                        fileFractionMapping = pepInput.getFileFractionMappingFromCSV(deNovoFile)
                    ## We should now have fileMapping, a list of tuples: (2-Digit Fraction Num, FileName)
                    ## mz[X]ML Files should be located in the dataDirectory
                    write2FileStdout(fh,'fileFractionMapping: %s'%fileFractionMapping)
                    symLinkDir += experimentName + '_' + str(os.getpid()) + '/'
                    dataFileSuffix = "mzML"
                    try:
                        ## Should be able to make the directory, and then remove it ##
                        os.makedirs(symLinkDir)
                        write2FileStdout(fh,'* Created temporary sym-link Directory for TagGraph mz[X]ML files \'%s\'' % (symLinkDir))
                        ## Lets write out the fileFractionMapping, pickled for easy reading/writing
                        mappingFilename = 'fileFractionMapping.pck'
                        mappingFilePath = os.path.join(symLinkDir, mappingFilename)
                        mappingOutput = open(mappingFilePath, 'wb')
                        pickle.dump(fileFractionMapping, mappingOutput)
                        mappingOutput.close()
                        ##Create a symbolic link pointing to source named link_name.
                        for currFilledFractionNumber, currFilename in fileFractionMapping:
                            ## Check if source file exists
                            currFilePath = dataDirectory + currFilename
                            if not os.path.exists(currFilePath):
                                fatalError = True
                                write2FileStdout(fh,'** FAILURE ** Data File \'%s\' referenced in de novo file does not exist in dataDirectory \'%s\'' % (currFilename,dataDirectory))
                            elif not os.access(currFilePath, os.R_OK):
                                fatalError = True
                                write2FileStdout(fh,'** FAILURE ** Data file \'%s\' Not Readable' % (currFilePath))
                            else:
                                currFractionFile = currFilePath
                                if currFractionFile.endswith('mzML'):
                                    dataFileSuffix = 'mzML'
                                elif currFractionFile.endswith('mzXML'):
                                    dataFileSuffix = 'mzXML'
                                else:
                                    fatalError = True
                                    dataFileSuffix = ''
                                    write2FileStdout(fh,'** FAILURE ** Data file \'%s\' must end in .mzML or .mzXML!' % (currFractionFile))
                                if not dataFileSuffix == '':
                                    symLinkFile = symLinkDir + experimentName + '_f' + currFilledFractionNumber + '.' + dataFileSuffix
                                    os.symlink(currFractionFile, symLinkFile)
                                    write2FileStdout(fh,'   * Created symLink \'%s\' to data file \'%s\'' % (symLinkFile, currFractionFile))
                    except OSError:
                        fatalError = True
                        write2FileStdout(fh,'** FAILURE ** Unable to Create Directory for TagGraph mz[X]ML sym-links at \'%s\'' % (symLinkDir))
                else:
                    fatalError = True
                    write2FileStdout(fh,'** FAILURE ** Required de novo TagGraph Parameter \'de_novo\' must be named .CSV or .XML/.PEPXML, found \'%s\'' % (deNovoFile))
    else:
        numFractions = tagGraphSectionMap['numFractions']
        if isNumeric(numFractions):
            if float(numFractions).is_integer():
                foundNumFractions = True
                write2FileStdout(fh,'* Found Required integer TagGraph Parameter \'numFractions\'  : \'%s\'' % (numFractions))
                numFractions = int(numFractions)
            else:
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'numFractions\' must be an integer value, found value \'%s\'' % (numFractions))
        else:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'numFractions\' must be a numeric value, found value \'%s\'' % (numFractions))  
    ## If we found numFractions, lets get the paths to the data files and make sym-links to them in a new directory ##
    ## sym-links will be named <ExperimentName>_f01.mz[X]ML, etc.                                                   ##
    if True == foundNumFractions:
        symLinkDir += experimentName + '_' + str(os.getpid()) + '/'
        dataFileSuffix = "mzML"
        try:
            ## Should be able to make the directory, and then remove it ##
            os.makedirs(symLinkDir)
            write2FileStdout(fh,'* Created temporary sym-link Directory for TagGraph mz[X]ML files \'%s\'' % (symLinkDir))
        except OSError:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Unable to Create Directory for TagGraph mz[X]ML sym-links at \'%s\'' % (symLinkDir))
        ##Create a symbolic link pointing to source named link_name.
        for currFraction in xrange(1,numFractions+1):
            filledFractionNumber = str(currFraction).zfill(2)
            if not str('fraction'+filledFractionNumber) in tagGraphSectionMap:
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'fraction%s\' not found in config file' % (filledFractionNumber))
            currFractionFile = tagGraphSectionMap['fraction'+filledFractionNumber]
            if currFractionFile.endswith('mzML'):
                dataFileSuffix = 'mzML'
            elif currFractionFile.endswith('mzXML'):
                dataFileSuffix = 'mzXML'
            else:
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** Data file \'%s\' must end in mzML or mzXML!' % (currFractionFile))
            symLinkFile = symLinkDir + experimentName + '_f' + filledFractionNumber + '.' + dataFileSuffix
            os.symlink(currFractionFile, symLinkFile)
            write2FileStdout(fh,'   * Created symLink \'%s\' to data file \'%s\'' % (symLinkFile, currFractionFile))
    # f must be an fmindex name of the form <basepath>.fm, where <basepath> is the full file path without the .fm extension, and the following files should exist: <basename>.fm.1, <basename>.seqnames.1, <basename>.offsets
    fmindexBase = ''
    if not 'fmindex' in tagGraphSectionMap:
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'fmindex\' (should be the basename of the fmindex files, ending in \'.fm\') not found in config file')
    else:
        fmParam = tagGraphSectionMap['fmindex']
        write2FileStdout(fh,'* Found Required fmindex TagGraph Parameter \'%s\'' % (fmParam))
        if fmParam.endswith('.fm'):
            fmindexBase = fmParam[:-3]
        else:
            fmindexBase = fmParam
        # Now lets check for 3 fmIndex files ending in: .fm.1, .offsets, and .seqnames.1
        fmFile = fmindexBase + ".fm.1"
        fmOffsetFile = fmindexBase + ".offsets"
        fmSeqnamesFile = fmindexBase + ".seqnames.1"
        if not os.path.isfile(fmFile):
            fatalError = True
            write2FileStdout(fh,'    ** FAILURE ** Could not find required fmindex file at \'%s\'' % (fmFile))
        elif not os.access(fmFile, os.R_OK):
            fatalError = True
            write2FileStdout(fh,'    ** FAILURE ** Could not Read required fmindex file \'%s\' (check permissions)' % (fmFile))
        else:
            write2FileStdout(fh,'   * Found Required readable fmindex file at \'%s\'' % (fmFile))
        if not os.path.isfile(fmOffsetFile):
            fatalError = True
            write2FileStdout(fh,'    ** FAILURE ** Could not find required fmindex Offset file at \'%s\'' % (fmOffsetFile))
        elif not os.access(fmOffsetFile, os.R_OK):
            fatalError = True
            write2FileStdout(fh,'    ** FAILURE ** Could not Read required fmindex Offset file \'%s\' (check permissions)' % (fmOffsetFile))
        else:
            write2FileStdout(fh,'   * Found Required readable fmindex Offset file at \'%s\'' % (fmOffsetFile))
        if not os.path.isfile(fmSeqnamesFile):
            fatalError = True
            write2FileStdout(fh,'    ** FAILURE ** Could not find required fmindex Seqnames file at \'%s\'' % (fmSeqnamesFile))
        elif not os.access(fmSeqnamesFile, os.R_OK):
            fatalError = True
            write2FileStdout(fh,'    ** FAILURE ** Could not Read required fmindex Seqnames file \'%s\' (check permissions)' % (fmSeqnamesFile))
        else:
            write2FileStdout(fh,'   * Found Required readable fmindex Seqnames file at \'%s\'' % (fmSeqnamesFile))
    ### Now lets Check the EM step parameters that can be checked before TG runs ###
    expectationMaximizationSectionMap = ConfigSectionMap(theConfig, "EM")
    '''
    -i: same as TG -i parameter
    -F all
    -M 100
    -C 20
    -B = <-o parameter from TG>/results.db [checked after TG]
    -E: Same as TG ExperimentName parameter.
    -o: Output Prefix, will create files with the prefix <EM -o parameter> in the directory specified by the <TG -o parameter>
    '''
    ## Arguments that must exist, and be numbers
    # Special Case: EMFractions must be 'all' or a number. Note: EMFractions is now assumed to always be 'all'
    requiredEMNumericArgs = ['maxIterations','initIterations']#,'EMFractions']
    ## Special Arguments:
    ## -o must be a string, the file prefix for the EM Output files (often 'EM_Results')
    ## Arguments that must exist, and be numbers ('EMFractions' is special, as a number or 'all')
    for currArg in requiredEMNumericArgs:
        if currArg in expectationMaximizationSectionMap:
            if isNumeric(expectationMaximizationSectionMap[currArg]):
                write2FileStdout(fh,'* Found Required EM Numeric Parameter \'%s\'  : \'%s\'' % (currArg, expectationMaximizationSectionMap[currArg]))
            else:
                fatalError = True
                write2FileStdout(fh,'** FAILURE ** Required EM Parameter \'%s\' must be a numeric value, found value \'%s\'' % (currArg, expectationMaximizationSectionMap[currArg]))
        else:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required EM Parameter \'%s\' not found in config file' % (currArg))
    ## Now Lets Handle the Special Cases
    # resultsPrefix (Output Prefix) must be a string
    emResultsPrefix = ''
    if not 'resultsPrefix' in expectationMaximizationSectionMap:
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Required EM Parameter \'resultsPrefix\' not found in config file')
    else:
        emResultsPrefix = expectationMaximizationSectionMap['resultsPrefix']
        write2FileStdout(fh,'* Found Required EM Parameter \'resultsPrefix\': \'%s\'' % (emResultsPrefix))
    #options = ArgLib.parse(['init', 'dtadir', 'peaks', 'output', 'ppmstd', 'modtolerance', 'unimoddict', 'maxcounts', 'modmaxcounts', 'fmindex', 'model', 'config'], optArgs=[{'opts': ('-x', '--splittaxon'), 'attrs': {'dest': 'splittaxon', 'action': 'store_true', 'default': False, 'help': 'Flag. For searches of metaproteomic databases, split identical context entries by taxon for accurate consideration via EM.'}}])
    ### If a fatal error was thrown, do not proceed ###
    if fatalError == True:
        write2FileStdout(fh,'*****  HALTING DUE TO FATAL ERROR IN TAGGRAPH OR EM PARAMETERS, SEE OUTPUT ABOVE!!! ')
        sys.exit(1)
    ## Lets set up the args properly for RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py ##
    tg_ppmstd = str(tagGraphSectionMap['ppmstd'])
    tg_modtolerance = str(tagGraphSectionMap['modtolerance'])
    tg_maxcounts = str(tagGraphSectionMap['maxcounts'])
    tg_modmaxcounts = str(tagGraphSectionMap['modmaxcounts'])
    tg_config = tagGraphSectionMap['config']
    tg_init = tagGraphSectionMap['init']
    tg_dtadir = symLinkDir ## tagGraphSectionMap['d']
    tg_model = tagGraphSectionMap['model']
    tg_output = tagGraphSectionMap['output']
    tg_unimoddict = tagGraphSectionMap['unimoddict']
    tg_fmindex = tagGraphSectionMap['fmindex']
    tg_peaks = '{\'' + tagGraphSectionMap['ExperimentName'] + '\': \'' + tagGraphSectionMap['de_novo'] + '\'}' # K = "{'e009133': '/lab/samba/shared/Users/Sam/20160630_Pulldown_dcas9_in_gel_digest_test_DENOVO_5/de_novo_peptides.csv'}"
    ### tg_output directory will now end with a slash
    if not tg_output.endswith('/'):
        tg_output += '/'
    tgArgs = []
    tgArgs.extend(['-p', '\"' + tg_ppmstd + '\"'])
    tgArgs.extend(['-l', '\"' + tg_modtolerance + '\"'])
    tgArgs.extend(['-M', '\"' + tg_maxcounts + '\"'])
    tgArgs.extend(['-C', '\"' + tg_modmaxcounts + '\"'])
    tgArgs.extend(['-c', '\"' + tg_config + '\"'])
    tgArgs.extend(['-i', '\"' + tg_init + '\"'])
    tgArgs.extend(['-d', '\"' + tg_dtadir + '\"'])
    tgArgs.extend(['-m', '\"' + tg_model + '\"'])
    tgArgs.extend(['-o', '\"' + tg_output + '\"']) 
    tgArgs.extend(['-Q', '\"' + tg_unimoddict + '\"'])
    tgArgs.extend(['-f', '\"' + tg_fmindex + '\"'])
    tgArgs.extend(['-K', '\"' + tg_peaks + '\"'])
    write2FileStdout(fh,'\nTG ARGS: %s\n\n'%tgArgs)
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'*** CALLING RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py from runTG.py')
    write2FileStdout(fh,MSGBORDER+"\n")
    DataFile.executeProcess(SCRIPTS_DIR,'RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py', tgArgs)
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END CALLING RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py from runTG.py')
    write2FileStdout(fh,MSGBORDER)
    ### VERIFY TG RUN ###
    '''
    Now lets check the TG output to make sure it ran correctly. We'll check for:
    * <output_dir>/results.db should exist and have size > 0 (do actual db check?)
    * The files <output_dir>/<experiment_name>_addPlausibleMods_poss_[combo/single]_mods.tdv both exist and have reasonable sizes
    * Check that output_dir/<experiment_name>/data/ contains directories of DTA files named <experiment_name>_f01/ etc
    * Check that output_dir/<experiment_name>/de_novo/<experiment_name>_PEAKS.csv/PEAKS_parsed.tdv/PEAKS_parsed_F1.tdv etc exist
    * Check that output_dir/<experiment_name>/taggraph/<experiment_name>_PEAKS_parsed_F1_TAGGRAPH.tdv etc exist
    * output_dir/<experiment_name>_CHECK.txt.<numFractions> contains count numbers for each fraction:
    -------------------------------
    /lab/samba/shared/Users/Sam/newtest/diet60_output
    Experiment Name diet60 ID 1
    Result Counts for 4 fractions
    F1: 399878
    F2: 395964
    F3: 346932
    F4: 270693
    -------------------------------
    '''
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'*** VERIFYING TAGGRAPH OUTPUTS in runTG.py ')
    write2FileStdout(fh,MSGBORDER)
    minDBFileSize = 1000000 ## 1Megabyte minimum db size after TG runs?
    minAddPlausibleModsFileSize = 2000 ## 10kBytes min size for <experiment_name>_addPlausibleMods_[combo/single]_mods.tdv files
    ## <output_dir>/results.db should exist and have size > 0 (do actual db check?)
    dbFile = tg_output + 'results.db'
    if not os.path.exists(dbFile):
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Required SQLITE DB File \'%s\' does not exist!!' % (dbFile))
    else:
        dbFileSize = os.path.getsize(dbFile)
        if dbFileSize < minDBFileSize:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required SQLITE DB File \'%s\' is too small: %d Bytes!!' % (dbFile, dbFileSize))
        else:
            write2FileStdout(fh,'* Found Required SQLITE DB File \'%s\', size %d Bytes OK' % (dbFile, dbFileSize)) 
    ## The files <output_dir>/<experiment_name>_addPlausibleMods_poss_[combo/single]_mods.tdv both exist
    singleModsFile = tg_output + experimentName + '_addPlausibleMods_poss_single_mods.tdv'
    comboModsFile  = tg_output + experimentName + '_addPlausibleMods_poss_combo_mods.tdv'
    if not os.path.exists(singleModsFile):
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Required Single Mods File \'%s\' does not exist!!' % (singleModsFile))
    else:
        singleModsFileSize = os.path.getsize(singleModsFile)
        if singleModsFileSize < minAddPlausibleModsFileSize:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required Single Mods File \'%s\' is too small: %d Bytes!!' % (singleModsFile, singleModsFileSize))
        else:
            write2FileStdout(fh,'* Found Required Single Mods File \'%s\', size %d Bytes OK' % (singleModsFile, singleModsFileSize))
    if not os.path.exists(comboModsFile):
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Required Combo Mods File \'%s\' does not exist!!' % (comboModsFile))
    else:
        comboModsFileSize = os.path.getsize(comboModsFile)
        if comboModsFileSize < minAddPlausibleModsFileSize:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required Combo Mods File \'%s\' is too small: %d Bytes!!' % (comboModsFile, comboModsFileSize))
        else:
            write2FileStdout(fh,'* Found Required Combo Mods File \'%s\', size %d Bytes OK' % (comboModsFile, comboModsFileSize))
    ## Check that output_dir/<experiment_name>/data/ contains directories of DTA files named <experiment_name>_f01/ etc
    dataDir = tg_output + experimentName + '/data/'
    for currFraction in xrange(1,numFractions+1):
        filledFractionNumber = str(currFraction).zfill(2)
        currDtaDirName = dataDir + experimentName + '_f' + filledFractionNumber
        if not os.path.exists(currDtaDirName):
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Missing directory of DTA files at: \'%s\'' % (currDtaDirName))
        elif not os.path.isdir(currDtaDirName):
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** \'%s\' exists but is not a Directory!' % (currDtaDirName))
        else:
            write2FileStdout(fh,'* Found DTA directory: \'%s\'' % (currDtaDirName))
    ## Check that output_dir/<experiment_name>/de_novo/<experiment_name>_PEAKS.csv/PEAKS_parsed.tdv/PEAKS_parsed_F1.tdv etc exist
    deNovoDir = tg_output + experimentName + '/de_novo/'
    deNovoCSV = deNovoDir + experimentName + '_PEAKS.csv'
    peaksParsed = deNovoDir + experimentName + '_PEAKS_parsed.tdv'
    fractionsParsedBase = deNovoDir + experimentName + '_PEAKS_parsed_F'
    if not os.path.exists(deNovoCSV):
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Missing de novo CSV File \'%s\' !!' % (deNovoCSV))
    else:
        write2FileStdout(fh,'* Found Required de novo CSV File \'%s\'' % (deNovoCSV)) 
    if not os.path.exists(peaksParsed):
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Missing Parsed de novo File \'%s\' !!' % (peaksParsed))
    else:
        write2FileStdout(fh,'* Found Required Parsed de novo File \'%s\'' % (peaksParsed))
    for currFraction in xrange(1,numFractions+1):
        currParsedFractionFile = fractionsParsedBase + str(currFraction) + '.tdv'
        if not os.path.exists(currParsedFractionFile):
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Missing Parsed de novo Fraction File \'%s\' !!' % (currParsedFractionFile))
        else:
            write2FileStdout(fh,'* Found Required Parsed de novo Fraction File \'%s\'' % (currParsedFractionFile))
    ## Check that output_dir/<experiment_name>/taggraph/<experiment_name>_PEAKS_parsed_F1_TAGGRAPH.tdv etc exist
    taggraphDir = tg_output + experimentName + '/taggraph/'
    taggraphParsedBase = taggraphDir + experimentName + '_PEAKS_parsed_F'
    taggraphParsedSuffix = '_TAGGRAPH.tdv'
    for currFraction in xrange(1,numFractions+1):
        currTaggraphFractionFile = taggraphParsedBase + str(currFraction) + taggraphParsedSuffix
        if not os.path.exists(currTaggraphFractionFile):
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Missing Parsed TagGraph Fraction File \'%s\' !!'%(currTaggraphFractionFile))
        else:
            write2FileStdout(fh,'* Found Required Parsed TagGraph Fraction File \'%s\''%(currTaggraphFractionFile))
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END VERIFYING TAGGRAPH OUTPUTS in runTG.py')
    write2FileStdout(fh,MSGBORDER)
    ### END VERIFY TG RUN ###
    ### If a fatal error was thrown, do not proceed ###
    if fatalError == True:
        write2FileStdout(fh,'*****  HALTING DUE TO FATAL ERROR IN VERIFYING TAGGRAPH RUN, SEE OUTPUT ABOVE!!')
        sys.exit(1)
    ## Copy configuration file to output tree for safe keeping ##
    configFileBaseName = os.path.basename(configFileName)
    checkConfigDestination = tg_output
    if os.path.exists(checkConfigDestination + configFileBaseName):
        write2FileStdout(fh,'** WARNING ** config file \'%s\' already exists in output directory \'%s\'' % (configFileBaseName,checkConfigDestination))
    else:
        shutil.copy(configFileName, checkConfigDestination)
        write2FileStdout(fh,'* Successfully copied Configuration File \'%s\' to Output Directory \'%s\'' % (configFileName,checkConfigDestination))
    ## Lets set up the args properly for ComputeEMProbabilitiesFromDB.py ##
    '''
    -i: same as TG -i parameter
    -F all
    -M 100
    -C 20
    -B = <-o parameter from TG>/results.db [checked after TG runs]
    -E: Same as TG ExperimentName parameter.
    -o: Output Prefix, will create files with the prefix <EM -o parameter> in the directory specified by the <TG -o parameter>
    '''
    em_init = tg_init
    em_fractions = 'all' ## EMFractions is always 'all' now! ## = str(expectationMaximizationSectionMap['EMFractions'])
    em_maxIterations = str(expectationMaximizationSectionMap['maxIterations'])
    em_initIterations = str(expectationMaximizationSectionMap['initIterations'])
    em_dbLocation = tg_output + 'results.db'
    em_experimentName = tagGraphSectionMap['ExperimentName']
    em_output = tg_output
    if not em_output.endswith('/'):
        em_output += '/'
    em_output += emResultsPrefix
    emArgs = []
    emArgs.extend(['-i', '\"' + em_init + '\"'])
    emArgs.extend(['-F', '\"' + em_fractions + '\"']) 
    emArgs.extend(['-M', '\"' + em_maxIterations + '\"'])
    emArgs.extend(['-C', '\"' + em_initIterations + '\"'])
    emArgs.extend(['-B', '\"' + em_dbLocation + '\"']) 
    emArgs.extend(['-E', '\"' + em_experimentName + '\"'])
    emArgs.extend(['-o', '\"' + em_output + '\"']) 
    write2FileStdout(fh,'EM ARGS: %s\n'% emArgs)
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'*** CALLING ComputeEMProbabilitiesFromDB.py from runTG.py')
    write2FileStdout(fh,MSGBORDER+"\n")
    DataFile.executeProcess(SCRIPTS_DIR,'ComputeEMProbabilitiesFromDB.py', emArgs)
    write2FileStdout(fh,'*** command executed: python ComputeEMProbabilitiesFromDB.py %s'%emArgs)
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END CALLING ComputeEMProbabilitiesFromDB.py from runTG.py')
    write2FileStdout(fh,MSGBORDER)
    EMProbs_TOPONLY = tg_output + 'EM_Results_EMProbs_END_TOPONLY.tdv'
    if not os.path.exists(EMProbs_TOPONLY):
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Missing EMProbs END TOPONLY file \'%s\'.' % (EMProbs_TOPONLY))
        sys.exit(1)
    else:
        write2FileStdout(fh,'* Found EMProbs END TOPONLY file \'%s\'' % (EMProbs_TOPONLY))
    write2FileStdout(fh,"\n\n"+MSGBORDER)
    write2FileStdout(fh,'*** CALLING verify EM result tests from runTG.py')
    write2FileStdout(fh,"\ntime now: @ %s"%datetime.datetime.now())
    result=verifyEM.verifyEM(tg_output)
    write2FileStdout(fh,result)
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,"\ntime now: @ %s"%datetime.datetime.now())
    write2FileStdout(fh,'*** END CALLING verify EM result tests from runTG.py')
    write2FileStdout(fh,MSGBORDER)
    topResultsFile = tg_output + experimentName + '_TopResults.tdv'
    if not os.path.exists(topResultsFile):
        fatalError = True
        write2FileStdout(fh,'** FAILURE ** Missing TopResult file \'%s\'.' % (topResultsFile))
        sys.exit(1)
    else:
        write2FileStdout(fh,'* Found TopResult file \'%s\'' % (topResultsFile))
    outputPerFraction="No"
    write2FileStdout(fh,'**** start parseResultsDB process: %s'%(datetime.datetime.now()))
    FDRCutoff=0.01
    logEMCutoff=100
    DisplayProteinNum=5
    if "outputPerFraction" in generalSectionMap:
        if True==theConfig.getboolean('General','outputPerFraction'):
            outputPerFraction="Yes"
    if "FDRCutoff" in generalSectionMap:
        if isNumeric(generalSectionMap["FDRCutoff"]):
            write2FileStdout(fh,'* Found  Numeric TagGraph Parameter \'%s\'  : \'%s\'' % ("FDRCutoff", generalSectionMap["FDRCutoff"]))
        FDRCutoff=generalSectionMap['FDRCutoff']
    if "logEMCutoff" in generalSectionMap:
        if isNumeric(generalSectionMap["logEMCutoff"]):
            write2FileStdout(fh,'* Found  Numeric TagGraph Parameter \'%s\'  : \'%s\'' % ("logEMCutoff", generalSectionMap["logEMCutoff"]))
        logEMCutoff=generalSectionMap['logEMCutoff']
    if "DisplayProteinNum" in generalSectionMap:
        if isNumeric(generalSectionMap["DisplayProteinNum"]):
            write2FileStdout(fh,'* Found  Numeric TagGraph Parameter \'%s\'  : \'%s\'' % ("DisplayProteinNum", generalSectionMap["DisplayProteinNum"]))
        DisplayProteinNum=generalSectionMap['DisplayProteinNum']
    writeTopArgs=[]
    writeTopArgs.extend(['\"' + tg_output + '\"'])
    writeTopArgs.extend(['\"' + tg_init + '\"'])
    writeTopArgs.extend(['\"' + outputPerFraction + '\"'])
    writeTopArgs.extend(['\"' + str(FDRCutoff) + '\"'])
    writeTopArgs.extend(['\"' + str(logEMCutoff) + '\"'])
    writeTopArgs.extend(['\"' + str(DisplayProteinNum) + '\"'])
    ## Now lets parse the original TG tab-delimted format ##
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'*** CALLING parseResultsDB.py from runTG.py')
    write2FileStdout(fh,MSGBORDER+"\n")
    DataFile.executeProcess(SCRIPTS_DIR, 'parseResultsDB.py',writeTopArgs)
    write2FileStdout(fh,'*** command executed: python parseResultsDB.py %s'%writeTopArgs)
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END CALLING parseResultsDB.py from runTG.py')
    write2FileStdout(fh,'**** done parseResultsDB process: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER)
    topResultsFinalFile = tg_output + experimentName + '_TopResults*.txt'
    foundFile=0
    if len(glob.glob(topResultsFinalFile))>0:
        foundFile=1
    if foundFile==0:
        fatalError=True
        write2FileStdout(fh,'** FAILURE ** Missing result file \'%s\' from parseResultsDB.py process. Please check.' % (topResultsFinalFile))
        sys.exit(1)
    if 'generatePepXML' in generalSectionMap:
        if True == theConfig.getboolean('General','generatePepXML'):
            ## Now lets generate the output in PepXML format ##
            '''
            python -u /lab/scratch/taggraph_sarah/taggraphsourcecode/database/resultPepXML.py \
            tg_init-i /lab/scratch/taggraph_sarah/taggraphsourcecode/resources/TAG_GRAPH_Tryp_CysCarbam_MetOx.ini \
            tg_ppmstd-p 10 \
            tg_modtolerance-l 0.1 \
            tg_maxcounts-M 400 \
            tg_modmaxcounts-C 200 \
            tg_fmindex-f /var/www/html/TAG_GRAPH/lib/databases/20141209_UniHUMAN_cRAP_ILEq.fm \
            tg_model-m /lab/scratch/taggraph_sarah/taggraphsourcecode/resources/AllChargeDist_posOnlyDependence_20150808_HumanProt500000.pck \
            xxxx-c /lab/scratch/taggraph_sarah/taggraphsourcecode/resources/AllChargeDist_posOnlyDependence_20150808.txt \
            tg_unimoddict-Q /lab/scratch/taggraph_sarah/taggraphsourcecode/resources/unimodDict_noLabels.pck \
            tg_output-o /lab/samba/shared/Users/Sarah/taggraph/testmzml/output/ \
            tg_dtadir-d /lab/samba/shared/Users/Sarah/taggraph/testmzml \
            >& /lab/samba/shared/Users/Sarah/taggraph/testmzml/OutputpepXML.txt
            '''
            pepArgs = []
            pepArgs.extend(['\"' + tg_init + '\"'])
            pepArgs.extend(['\"' + tg_ppmstd + '\"'])
            pepArgs.extend(['\"' + tg_modtolerance + '\"'])
            pepArgs.extend(['\"' + tg_maxcounts + '\"'])
            pepArgs.extend(['\"' + tg_modmaxcounts + '\"'])
            pepArgs.extend(['\"' + tg_fmindex + '\"'])  # tagGraphSectionMap['fmindex']
            pepArgs.extend(['\"' + tg_model + '\"'])  # tagGraphSectionMap['model']
            #pepArgs.extend(['\"' + tg_config + '\"'])  # tagGraphSectionMap['config']
            pepArgs.extend(['\"' + tg_unimoddict + '\"'])  # tagGraphSectionMap['unimoddict']
            pepArgs.extend(['\"' + tg_output + '\"'])  # tagGraphSectionMap['output']
            pepArgs.extend(['\"' + tg_dtadir + '\"'])  # symLinkDir
            pepArgs.extend(['\"' + str(FDRCutoff) + '\"'])
            pepArgs.extend(['\"' + str(logEMCutoff) + '\"'])
            write2FileStdout(fh,MSGBORDER)
            write2FileStdout(fh,'*** CALLING generatePepXMLDBperFrac.py from runTG.py')
            write2FileStdout(fh,MSGBORDER+"\n")
            DataFile.executeProcess(SCRIPTS_DIR, 'generatePepXMLDBperFrac.py',pepArgs)
            write2FileStdout(fh,'*** command: python generatePepXMLDBperFrac.py %s'%pepArgs)
            write2FileStdout(fh,"\n"+MSGBORDER)
            write2FileStdout(fh,'*** END CALLING generatePepXMLDBperFrac.py from runTG.py')
            write2FileStdout(fh,MSGBORDER)
            '''
            Now lets clean up our temporary items and copied data files as configured! ###
            We need to:
            * Remove the sym-link directory in /tmp/ (symLinkDir)
            * If cleanMzDataFilesFromOutput is True, clean the dataDir (<output_dir>/<experiment_name>/data/)
            directory of mz[X]ML files and the DTA directories of the same name
            '''
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'***    CLEANING UP')
    write2FileStdout(fh,MSGBORDER)
    ### Remove the sym-link directory in /tmp/ (symLinkDir)
    shutil.rmtree(symLinkDir)
    if os.path.exists(symLinkDir):
        write2FileStdout(fh,'** FAILURE ** Failed to removed temporary symbolic link directory \'%s\'' % (symLinkDir))
    else:
        write2FileStdout(fh,'* Successfully removed temporary symbolic link directory \'%s\'' % (symLinkDir))
    if 'cleanInputDataFilesFromOutput' in generalSectionMap:
        if True == theConfig.getboolean('General','cleanInputDataFilesFromOutput'):
            shutil.rmtree(dataDir)
            #os.makedirs(dataDir)
            write2FileStdout(fh,'* Removed mz[X]ML and DTA files from data directory \'%s\' (cleanInputDataFilesFromOuput is True)' % (dataDir))
        else:
            write2FileStdout(fh,'* Leaving mz[X]ML and DTA files in data directory \'%s\' (cleanInputDataFilesFromOuput is False)' % (dataDir))
    if 'cleanIntermediateFiles' in generalSectionMap:
        denovoOutputDir = tg_output +'/'+ experimentName + '/de_novo/'
        taggraphOutputDir = tg_output +'/'+ experimentName + '/taggraph/'
        experimentOutputDir = tg_output +'/'+ experimentName
        if True == theConfig.getboolean('General','cleanIntermediateFiles'):
            shutil.rmtree(denovoOutputDir)
            shutil.rmtree(taggraphOutputDir)
            if os.path.exists(dataDir):
                shutil.rmtree(dataDir)
            shutil.rmtree(experimentOutputDir)
            files = os.listdir(tg_output)
            for file in files:
                if (file.endswith(".tdv") or (file.find("_CHECK.txt.") > 0) or file.endswith(".db") or file.endswith(".log")):
                    if (os.path.exists(os.path.join(tg_output,file))):
                        write2FileStdout(fh,"remove %s" % os.path.join(tg_output,file))
                        os.remove(os.path.join(tg_output,file))
                    else:
                        write2FileStdout(fh,"keeper %s" % os.path.join(tg_output,file))
            write2FileStdout(fh,'* Removed mz[X]ML and Intermediate files from output directory \'%s\' (cleanIntermediateFiles is True)' % (dataDir))
        else:
            write2FileStdout(fh,'* Leaving mz[X]ML and Intermediate files in output directory \'%s\' (cleanIntermediateFiles is False)' % (dataDir))
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'***  END CLEANING UP')
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'%s'%TAGGRAPH_CONFIG_FOOTER)
    write2FileStdout(fh,'**** end TagGraph process: %s'%(datetime.datetime.now()))
    fh.close()
    #move file back to output folder:
    toDest=tg_output+"runReport.log"
    shutil.move(runCapture,toDest)
    sys.exit(0)

if __name__ == '__main__':
    main(sys.argv[1:])

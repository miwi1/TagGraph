# Created on Sep 23, 2011
# @author: Arun
#
# 20160921    slin  modify to the take value from a init File:  init.config under project directory
# 20171113    slin  made init.config optinal for the project.
#
#


import sys
import os
from optparse import OptionParser

import Constants
import DataFile

def setupParams():
    fullPath=os.path.dirname(__file__)
    defaultVal = dict() 
    defaultVal["ppmstd"]=5
    defaultVal["modtolerance"]=0.02
    defaultVal["unimoddict"]=""
    defaultVal["maxcounts"]=300
    defaultVal["ambigpenalty"]=20
    defaultVal["ppmpenalty"]=20
    defaultVal["peaks"]=""
    defaultVal["subgraphcut"]=300
    defaultVal["unfiltered"]=False
    defaultVal["roc"]=False
    defaultVal["modmaxcounts"]=200
    defaultVal["alpha"]=0.90
    defaultVal["model"]=""
    defaultVal["config"]=""
    defaultVal["init"]=""
    defaultVal["dtadir"]=""
    defaultVal["fmindex"]=""
    defaultVal["output"]=""
    defaultVal["verbose"]=True
    defaultVal["unfiltered"]=False
    defaultVal["srchid"]=False
    defaultVal["sqlitedb"]=False
    defaultVal["experimentname"]=False
    defaultVal["fracs"]=False

    if os.path.exists(fullPath+"/init.config"):
        with open(fullPath+"/init.config") as f:
            for line in f:
                line = line.strip()
                if "=" in line:
                    splitLine = line.split("=")
                    if splitLine[0] in defaultVal.keys():
                        defaultVal.update({splitLine[0]:splitLine[1]})
                    else:
                        defaultVal[splitLine[0]]=splitLine[1]

    A = {
    'dtadir': {'opts': ('-d', '--dtadir'), 'attrs': {'type': 'string', 'dest': 'dtadir','default': defaultVal["dtadir"],'help': 'Directory of .dta files or dta.tgz file'}},  
    'config': {'opts': ('-c', '--config'), 'attrs': {'type': 'string', 'dest': 'config','default': defaultVal["config"], 'help': 'Path to model configuration file. If none is provided, will check the params file.'}},
    'model': {'opts': ('-m', '--model'), 'attrs': {'type': 'string', 'dest': 'model','default': defaultVal["model"],'help': 'Path to probabilistic model file. If none is provided, will check the params file.'}},
    'output': {'opts': ('-o', '--output'), 'attrs': {'type': 'string', 'dest': 'output','default': defaultVal["output"],'help': 'Name of output file'}},
    'database': {'opts': ('-D', '--database'), 'attrs': {'type': 'string', 'dest': 'database', 'help': 'Database containing formatting parameters of search programs'}},
    'taggraph': {'opts': ('-T', '--taggraph'), 'attrs': {'type': 'string', 'dest': 'taggraph', 'help': 'Location of TAG_GRAPH.py (.tdv) output (NOT Filtered or Parsed output!)'}},
    'columns': {'opts': ('-k', '--columns'), 'attrs': {'type': 'string', 'dest': 'columns', 'help': 'Path to pickled python tuple listing columns to include in LADS .tdv output, default is light scan number, heavy scan number (where applicable), [M+H], score, sequence'}},
    'verbose': {'opts': ('-v', '--verbose'), 'attrs': {'action': 'store_true', 'default': defaultVal["verbose"], 'dest': 'verbose', 'help': 'Print status updates to std out'}},
    'unfiltered': {'opts': ('-u', '--unfiltered'), 'attrs': {'action': 'store_true', 'default': defaultVal["unfiltered"], 'dest': 'unfiltered', 'help': 'Flag. Set to write unfiltered output of TAG-GRAPH'}},
    'roc': {'opts': ('-r', '--roc'), 'attrs': {'action': 'store_true', 'default': defaultVal["roc"], 'dest': 'roc', 'help': 'Flag. Set to write ROC curve used to derive score cutoffs'}},
    'ppmstd': {'opts': ('-p', '--ppmstd'), 'attrs': {'type': 'float', 'dest': 'ppmstd', 'default': defaultVal["ppmstd"], 'help': 'Expected standard deviation in ppm error distributions of fragment ions. Recommend 5 for HCD 30,000 resolution'}},
    'modtolerance': {'opts': ('-l', '--modtolerance'), 'attrs': {'type': 'float', 'dest': 'modtolerance', 'default': defaultVal["modtolerance"], 'help': 'Maximum deviation between expected and observed mod mass to consider mod as a candidate match.'}},
    'number': {'opts': ('-y', '--number'), 'attrs': {'type': 'int', 'dest': 'number', 'help': 'Just a number (integer). Use varies with program.'}},
    'ambigpenalty': {'opts': ('-a', '--ambigpenalty'), 'attrs': {'type': 'float', 'dest': 'ambigpenalty', 'default': defaultVal["ambigpenalty"], 'help': 'Score penalty for ambiguous edges in spectrum graph'}},
    'ppmpenalty': {'opts': ('-A', '--ppmpenalty'), 'attrs': {'type': 'float', 'dest': 'ppmpenalty', 'default': defaultVal["ppmpenalty"], 'help': 'Maximum score penalty for ppm deviation from true sequence mass in edges of spectrum graph'}},
    'scans': {'opts': ('-N', '--scans'), 'attrs': {'type': 'string', 'dest': 'scans','help': 'Optional. A list of scans to specifically target using tag-graph, separated by commas (e.g., scan1,scan2,etc.)'}},    
    'excludescans': {'opts': ('-s', '--excludescans'), 'attrs': {'type': 'string', 'dest': 'excludescans', 'help': 'Optional. A list of scans to exclude from tag-graph analysis, separated by commas (e.g., scan1,scan2,etc.)'}},
    'maxcounts': {'opts': ('-M', '--maxcounts'), 'attrs': {'type': 'int', 'dest': 'maxcounts', 'default': defaultVal["maxcounts"], 'help': 'Maximum number of times a maximum matching substring can be found in sequence database. Recommend 1,000 for single organism database (i.e., human uniprot) and 5,000 for nr or 6-frame translations'}},
    'modmaxcounts': {'opts': ('-C', '--modmaxcounts'), 'attrs': {'type': 'int', 'dest': 'modmaxcounts', 'default': defaultVal["modmaxcounts"], 'help': 'Maximum number of times a maximum matching substring can be found in sequence database for in-exact de novo peptide matches subsequently evaluated for presence of modifications. Recommend 200 for single organism database (i.e., human uniprot) and 1,000 for nr or 6-frame translations'}},
    'taxonomy': {'opts': ('-y', '--taxonomy'), 'attrs': {'type': 'string', 'dest': 'taxonomy', 'help': 'File containing pickled python taxonomy tree'}},
    'gimap': {'opts': ('-P', '--gimap'), 'attrs': {'type': 'string', 'dest': 'gimap', 'help': 'Location of gimap database (dbm)'}},
    'alpha': {'opts': ('-r', '--alpha'), 'attrs': {'type': 'float', 'dest': 'alpha', 'default': defaultVal["alpha"], 'help': 'Minimum score ratio between suboptimal and optimal paths'}},
    'lads': {'opts': ('-L', '--lads'), 'attrs': {'type': 'string', 'dest': 'lads', 'help': 'Dictionary mapping location of LADS search results (.tdv) to its name in the output'}},
    'mascot': {'opts': ('-t', '--mascot'), 'attrs': {'type': 'string', 'dest': 'mascot', 'help': 'Dictionary mapping location of MASCOT search results (.csv) to its name in the output'}},
    'pepnovo': {'opts': ('-O', '--pepnovo'), 'attrs': {'type': 'string', 'dest': 'pepnovo', 'help': 'Location of pepNOVO search results (.tdv)'}},
    'pnovo': {'opts': ('-N', '--pnovo'), 'attrs': {'type': 'string', 'dest': 'pnovo', 'help': 'Location of pNovo search results (.tdv)'}},
    'peaks': {'opts': ('-K', '--peaks'), 'attrs': {'type': 'string', 'dest': 'peaks', 'default': defaultVal["peaks"], 'help': 'Location of exported PEAKS search results'}},
    'fmindex': {'opts': ('-f', '--fmindex'), 'attrs': {'type': 'string', 'dest': 'fmindex','default': defaultVal["fmindex"], 'help': 'Location of fmindex to query against'}},
    'fasta': {'opts': ('-f', '--fasta'), 'attrs': {'type': 'string', 'dest': 'fasta', 'help': 'Location of FASTA file'}},
    'testvals':  {'opts': ('-V', '--testvals'), 'attrs': {'type': 'string', 'dest': 'testvals', 'help': 'Input of the form a,b,c, etc. which evals to a tuple. Used to test a program over multiple values of a parameter'}},
    'progdict': {'opts': ('-r', '--progdict'), 'attrs': {'type': 'string', 'dest': 'progdict', 'help': 'Dictionary which maps names of searches in CompareSearches.py output to the program which generated them'}},
    'denovo': {'opts': ('-n', '--denovo'), 'attrs': {'type': 'string', 'dest': 'denovo', 'help': 'location of de novo search results'}},
    'mainprogname': {'opts': ('-R', '--mainprogname'), 'attrs': {'type': 'string', 'dest': 'mainprogname', 'help': 'name of output to compare against others when analyzing compareSearches.py output'}},
    'subgraphcut': {'opts': ('-u', '--subgraphcut'), 'attrs': {'type': 'float', 'dest': 'subgraphcut', 'default': defaultVal["subgraphcut"], 'help': 'maximum edge length to be taken for complete sequencing (as opposed to creating a sub-spectrum graph) in resolving ambiguous edges'}},
    'cutoff': {'opts': ('-b', '--cutoff'), 'attrs': {'type': 'float', 'dest': 'cutoff', 'help': 'Number, use varies with program'}},
    'init': {'opts': ('-i', '--init'), 'attrs': {'type': 'string', 'dest': 'init','default': defaultVal["init"], 'help': 'Initiation file used to configure TAG-GRAPH'}},
    'symbolmap': {'opts': ('-S', '--symbolmap'), 'attrs': {'type': 'string', 'dest': 'symbolmap', 'help': 'Key mapping program symbols from search program output the their respective modification types or amino acids(path to pickled python dict).'}},
    'modalpha': {'opts': ('-q', '--modalpha'), 'attrs': {'type': 'float', 'dest': 'modalpha', 'help': 'suboptimal pathing cutoff alpha for modification amino acids'}},
    'unimoddict': {'opts': ('-Q', '--unimoddict'), 'attrs': {'type': 'string', 'dest': 'unimoddict','default': defaultVal["unimoddict"],'help': 'location of pickled unimod dictionary'}},
    'combined': {'opts': ('-C', '--combined'), 'attrs': {'type': 'string', 'dest': 'combined', 'help': 'Dictionary which maps location of combined sequest/mascot results to their name in the output'}},
    'srchid': {'opts': ('-H', '--srchid'), 'attrs': {'type': 'int', 'dest': 'srchid', 'default': defaultVal["srchid"], 'help': 'Dictionary which maps name of database search in output to the searchID of interest in aggregated exported database search file (.csv), optional'}},
    'sqlitedb':  {'opts': ('-B', '--sqlitedb'), 'attrs': {'type': 'string', 'dest': 'sqlitedb', 'default': defaultVal["sqlitedb"], 'help': 'Location of sqlite db to connect to, must be an absolute path.'}},
    'experimentname': {'opts': ('-E', '--experimentname'), 'attrs': {'type': 'string', 'dest': 'experimentname', 'default': defaultVal["experimentname"], 'help': 'Name of experiment to pull from database.'}},
    'fracs': {'opts': ('-F', '--fracs'), 'attrs': {'type': 'string', 'dest': 'fracs', 'default': defaultVal["fracs"], 'help': 'Fractions to use from database. Argument should be in form frac1,frac2,frac3,etc.'}}
    }
    return A

def parse(arglist, optArgs=[]):
    A=setupParams()
    parser = OptionParser(usage="%prog [options]")        
    for arg in arglist:
        parser.add_option(*A[arg]['opts'], **A[arg]['attrs'])
    for optArg in optArgs:
        parser.add_option(*optArg['opts'], **optArg['attrs'])
    (options, args) = parser.parse_args()
    return options

def parseInitFile(init, options):
    A=setupParams()
    paramsDict = DataFile.parseParams(init)
    for param in paramsDict['Parameters'].keys():
        try:
            paramType = A[param]['attrs']['type']
            val = init['Parameters'][param]
            if paramType != 'string':
                val = getattr('__builtin__', paramType)(val)
            setattr(options, param, val)
        except KeyError:
            pass
    return paramsDict

def getProgDict(progs, options, progDict={}):
    for prog in progs:
        if hasattr(options, prog.lower()) and getattr(options, prog.lower()) != None:
            fDict = eval(getattr(options, prog.lower()))
            for name in fDict.values():
                progDict[name] = prog
    return progDictogDict

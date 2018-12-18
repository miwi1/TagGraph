#slin 201707  added lib path and others.

import os
import sys
import glob

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

import Database
import DataFile
import ArgLib
import argparse

if __name__ == '__main__':
	#eg:  python Build_FMIndex_new.py /lab/samba/shared/Users/Sarah/HC_FMindex/HumanChimpRNAseq_V2.1.simplified.20170516.fasta
    parser = argparse.ArgumentParser(description='Generate FMIndex files from a FASTA file.')
    parser.add_argument('-o', '--output', default='', help='Prefix for FMIndex output files. Base of Fasta input file will be used if not supplied.')
    parser.add_argument('fasta', help='Input FASTA filename. Must end with .fasta')
    args = parser.parse_args()

    if not os.path.basename(args.fasta.lower()).endswith('.fasta'):
        print 'Error! FASTA input \'%s\' doesn\'t end with \'.fasta\'! Exiting...' % (args.fasta)
        sys.exit()

    print 'FASTA input file: %s' % (args.fasta)
    
    if args.output == '':
        print 'BLANK OUTPUT basename - using the FASTA input file  base: \'%s\'' % (os.path.basename(args.fasta[:-6]))
        args.output = os.path.basename(args.fasta[:-6])
    else:
        print 'Output basename: %s' % (args.output)

    Database.makeDBForFMIndexFromFASTA(args.fasta, args.output)

    fmbuild_loc  =  os.path.abspath( os.path.join(os.path.join( PAR_DIR, 'lib'), 'fmindex') )
    for fm_formatted in glob.glob(args.output + '*fmFormatted*'):
        DataFile.executeProcess(fmbuild_loc, 'fmbuild', ['-v', fm_formatted, args.output + '.fm%s'%os.path.splitext(fm_formatted)[1]], interpreter = False)
        os.remove(fm_formatted)

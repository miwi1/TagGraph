import os
import sys
import glob
import Database
import DataFile
import argparse

PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)



#eg:  python Build_FMIndex_new.py /lab/samba/shared/Users/Sarah/HC_FMindex/HumanChimpRNAseq_V2.1.simplified.20170516.fasta
parser = argparse.ArgumentParser(description='Generate FMIndex files from a FASTA file.')
parser.add_argument('-o', '--output', default='', help='Prefix for FMIndex output files. Base of Fasta input file will be used if not supplied.')
parser.add_argument('fasta', help='Input FASTA filename. Must end with .fasta')
args = parser.parse_args()

if not os.path.basename(args.fasta.lower()).endswith('.fasta'):
    raise FileNotFoundError("Error! FASTA input {} doesn't end with .fasta!".format(args.fasta))

if args.output == '':
    output_filename = os.path.basename(args.fasta[:-6])
    print('BLANK OUTPUT basename - using the FASTA input file  base: {}'.format(output_filename))
    args.output = output_filename

print(args)
sys.exit()

Database.makeDBForFMIndexFromFASTA(args.fasta, args.output)

fmbuild_loc  =  os.path.abspath( os.path.join(os.path.join( PAR_DIR, 'lib'), 'fmindex') )
for fm_formatted in glob.glob(args.output + '*fmFormatted*'):
    DataFile.executeProcess(fmbuild_loc, 'fmbuild', ['-v', fm_formatted, args.output + '.fm%s'%os.path.splitext(fm_formatted)[1]], interpreter = False)
    os.remove(fm_formatted)

# bitbucket antigeMS issue #22
# input file format:  
# {accession}
# [{peptide1}	{val1}]+
# --- cut here ---
# Q05086
# SKLPLAAQG	1
# FRRGFHMVTNES	0.2
# PFEEFINEP	3
# --- end here ---
#
#  eg: python proteinPeptidesHistogramMap.py input.txt output_smallfasta.txt test.fasta 
#      python proteinPeptidesHistogramMap.py input.txt output_fullfasta.txt /lab/samba/shared/2015-02-11_REV_Plasmo_acpGFP_bsaBiotin_HUMAN.fasta
#
#      python proteinPeptidesHistogramMap.py /lab/samba/shared/Users/Sarah/proteinPeptides/input.txt ~/output_1.txt test.fasta
# 20171006 slin script creation.
# 20171006 slin add error messages.

import sys
import subprocess
from itertools import islice
from decimal import Decimal
import os.path
MARKER=">"

def readReqInfo(inputFile):
	peptidesList=[]
	accession=""
	rowcnt=0
	with open(inputFile,'r') as lines:
		for line in lines:
			line=line.strip()
			if (rowcnt==0):
				accession=line
			else:
				peptideInfo=line.split("\t")
				peptidesList.append([peptideInfo[0],peptideInfo[1]])
			rowcnt=rowcnt+1
	return (accession,peptidesList)	
		
def findProteinViaAccession(fasta,accession):
	protein=""
	startLine=0
	endLine=0
	rowCnt=0
	with open(fasta) as lines:
		for line in lines:
			rowCnt=rowCnt+1
			line=line.strip()
			if (startLine>0) and (MARKER not in line):
				protein=protein+line
			if MARKER in line:  #only search when ">" exist
				if (startLine>0):
					endLine=rowCnt-1
					break
				if accession in line:
					startLine=rowCnt
	return protein

def createHistogram(protein,proteinArray,peptidesList):
	for peptideInfo in peptidesList:
		peptide=peptideInfo[0]
		assVal=(round(Decimal(peptideInfo[len(peptideInfo)-1]),6))  #not sure why peptideInfo[1] is not working.
		if (peptide in protein):
			startPos=protein.index(peptide)
			endPos=startPos+len(peptide)
			for i in range(startPos,endPos):
				proteinArray[i]=proteinArray[i]+assVal
	return proteinArray
	
if __name__ == '__main__':
	if (len(sys.argv)<4):
		print " Error! \n Usage: python proteinPeptidesHistogramMap.py inputFile outputFile fastaFile [checkfile]\n"
		print " eg: python proteinPeptidesHistogramMap.py input.txt output.txt test.fasta "
		print " eg: python proteinPeptidesHistogramMap.py input.txt output.txt test.fasta x "
		sys.exit(2)
	checkfile=0	
	inputFile=sys.argv[1]
	outputFile=sys.argv[2]
	fasta=sys.argv[3]
	if (len(sys.argv)==5):
		checkfile=sys.argv[4]
	inputNotExist=0
	inputNotExistMsg=""
	fastaNotExist=0
	fastaNotExistMsg=""
	if not (os.path.exists(inputFile)):
		inputNotExist=1
		inputNotExistMsg=("  input file:%s does not exist."%(inputFile))
	if not (os.path.exists(fasta)):
		fastaNotExist=1
		fastaNotExistMsg=("  fasta file:%s does not exist."%(fasta))
	if inputNotExist or fastaNotExist:
		print " Error:"
		if inputNotExist:
			print inputNotExistMsg
		if fastaNotExist:
			print fastaNotExistMsg
		print ""
		exit()
	(accession,peptidesList)=readReqInfo(inputFile)
	protein=findProteinViaAccession(fasta,accession)
	if (protein==""):
		print (" Error:\n  fasta file: %s got no matching accession: %s\n"%(fasta,accession))
		exit()
	proteinArray=[0]*len(protein)
	proteinArray=createHistogram(protein,proteinArray,peptidesList)
	fout=open(outputFile,'w')
	for i in range(len(proteinArray)): #position starts at 1
		fout.write("%s\t%s\n"%(i+1,proteinArray[i]))
	fout.close()
	print (" output file is located at %s\n"%(outputFile))
	if checkfile:
		fout1=open("checkfile.txt",'w')  #write out info for checking.
		fout1.write("%s\n"%(protein))
		fout1.write("%s\n"%(peptidesList))
		for j in range(len(proteinArray)): #position starts at 1
			fout1.write("%s\t%s\n"%(j+1,proteinArray[j]))
		fout1.close()
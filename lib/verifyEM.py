#20170915 slin #bitbucket issue#8

import sqlite3 as lite
import sys
import csv
import math
import itertools
import decimal
MAXLINE=50000
TESTCASES=[["SpectrumProbabilityScore",">=3"],["MatchingTagLength",">=8"],["ModSize",">0"]]

TABLEFIELDS="(ScanF varchat(50),Charge int,MatchingTagLength int,SpectrumProbabilityScore decimal(23,20),Specificity varchar(40),Modified int,PPMError decimal(14,10),ModSize varchar(50),NumMods int,UniqueSiblings int,ContextModVariants int,NumModOccurrences int,ModClass int,NumSingleModsFound int,Context varchar(300),ModContext varchar(300),ModTuple varcbar(300),EMProbability decimal(23,20),OneSubtractlg10EM decimal(23,20),PosSpectrumProb  decimal(23,20),NegSpectrumProb decimal(23,20),PosModProb decimal(23,20),NegModProb decimal(23,20),PosContextProb decimal(23,20),NegContextProb decimal(23,20),PosDBMatchProb decimal(23,20),NegDBMatchProb decimal(23,20),PosProteinCountProb decimal(23,20),NegProteinCountProb decimal(23,20),PosMissedCleavageProb decimal(23,20),NegMissedCleavageProb decimal(23,20),PosPPMErrorProb decimal(23,20),NegPPMErrorProb decimal(23,20),PosProb decimal(23,20),NegProb decimal(23,20),Log10EM decimal(23,20),EM decimal(23,20))"

def createEMFullTable(conn):
	conn.execute("create table if not exists EMFull"+TABLEFIELDS);
	conn.execute("delete from EMFull");

def readFile(infile,startLine,endLine,delimiter=' '):
	dataList=[]
	with open(infile,"r") as f:
		for line in itertools.islice(f,startLine,endLine):
			line=line.strip()
			splitLine = line.split(delimiter)
			#find EM to use for statistics test via the rule:
			#1) calculate LOG EM from EM probability field  EM probability pos 17
			#2) make new "EM" field which is {log(1-EM probability) if EM probability >= 0; OR log(EM probability) if EM probability < 0}
			#3) sort by new EM field
			Log10EM=math.log(float(splitLine[17]))
			if (splitLine[18]>-1*Log10EM):
				EM=splitLine[18]
			else:
				EM=Log10EM
			if (splitLine[7]=="Indeterminate"):
				ModSize=-1
			else:
				ModSize=splitLine[7]
			dataList.append((splitLine[0],splitLine[1],splitLine[2],splitLine[3],splitLine[4],splitLine[5],splitLine[6],ModSize,splitLine[8],splitLine[9],splitLine[10],splitLine[11],splitLine[12],splitLine[13],splitLine[14],splitLine[15],splitLine[16],splitLine[17],splitLine[18],splitLine[19],splitLine[20],splitLine[21],splitLine[22],splitLine[23],splitLine[24],splitLine[25],splitLine[26],splitLine[27],splitLine[28],splitLine[29],splitLine[30],splitLine[31],splitLine[32],splitLine[33],splitLine[34],Log10EM,EM))
	return dataList

def	populateEMTable(outputFolder,conn):
	infile=outputFolder+"/EM_Results_EMProbs_END_TOPONLY.tdv"
	lineCount=0
	f=open(infile)
	lineCount= sum(1 for line in f)
	iteration=(lineCount-1)/MAXLINE
	remainder=(lineCount-1)%MAXLINE
	if (remainder>0):
		iteration=iteration+1
	for i in range(0,iteration):
		startVal=1+i*MAXLINE
		endVal=1+(i+1)*MAXLINE
		if endVal>lineCount:
			endVal=lineCount
		dataList=readFile(infile,startVal,endVal,'	')
		loadEMFullTable(conn,dataList)
		
def loadEMFullTable(conn,dataList):
	conn.executemany("insert into EMFull(ScanF,Charge,MatchingTagLength,SpectrumProbabilityScore,Specificity,Modified,PPMError,ModSize,NumMods,UniqueSiblings,ContextModVariants,NumModOccurrences,ModClass,NumSingleModsFound,Context,ModContext,ModTuple,EMProbability,OneSubtractlg10EM,PosSpectrumProb,NegSpectrumProb,PosModProb,NegModProb,PosContextProb,NegContextProb,PosDBMatchProb,NegDBMatchProb,PosProteinCountProb,NegProteinCountProb,PosMissedCleavageProb,NegMissedCleavageProb,PosPPMErrorProb,NegPPMErrorProb,PosProb,NegProb,Log10EM,EM) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",dataList)
	conn.commit()
	
def getEMHighLowCnt(cursor):
	SQL="select count(*) from EMFull where EMProbability>=0.5"
	cursor.execute(SQL)
	highCountRow=cursor.fetchone()	
	EMHighCount=highCountRow[0]
	SQL="select count(*) from EMFull where EMProbability<0.5"
	cursor.execute(SQL)
	lowCountRow=cursor.fetchone()
	EMLowCount=lowCountRow[0]
	return EMHighCount,EMLowCount
	
def testCases(cursor):
	countsHigh=[]
	countsLow=[]	
	for case in TESTCASES:
		fieldName=case[0]
		fieldCondition=case[1]
		SQL="select count(*) from EMFull where EMProbability>=0.5 and "+fieldName+str(fieldCondition)
		cursor.execute(SQL)
		countRow=cursor.fetchone()	
		count=countRow[0]
		countsHigh.append(count)
		SQL="select count(*) from EMFull where EMProbability<0.5 and "+fieldName+str(fieldCondition)
		cursor.execute(SQL)
		countRow=cursor.fetchone()	
		count=countRow[0]
		countsLow.append(count)
	return countsHigh,countsLow

def calculation(EMHighCount,EMLowCount,countsHigh,countsLow):
# 4a) Spectrum score: frequency of high magnitude (>=3) Spectrum scores relative to total high EM or low EM populations. 
# For the A375 data set, this should be ~0.84 for high EM scores and 0.12 for the low EM population. 
# 4b) Tag Length: frequency of long (>=8) tags relative to total high EM or low EM populations. For the A375 data set, this should be ~0.81 for the high EM scores and 0.07 for the low EM scores. 
# 4c) Modification size: frequency of large (>0) relative to total high EM or low EM populations. For the A375 data set, this should be ~0.004 for the high EM scores and 0.24 for the low EM scores. 	
	highRatio_High=[]
	lowRatio_High=[]
	highRatio_Low=[]
	lowRatio_Low=[]
	for count in countsHigh:
		highRatio_High.append(float(count)/EMHighCount)
		lowRatio_High.append(float(count)/EMLowCount)
	for count in countsLow:
		highRatio_Low.append(float(count)/EMHighCount)
		lowRatio_Low.append(float(count)/EMLowCount)
	return (highRatio_High,lowRatio_High,highRatio_Low,lowRatio_Low)

def calculation_Verify(highRatio_High,lowRatio_Low):
#5) calculate x^2 statistic for each value in 4a-c: ((high EM value - low EM value)^2) / low EM value 
#For Spectrum score, this should be 4.1 
#For Tag length this should be 7.7 
#For Mod size, this should be 0.24 
	verifyVals=[]
	for i in range(len(highRatio_High)):
		val=((highRatio_High[i]-lowRatio_Low[i])**2)/lowRatio_Low[i]
		verifyVals.append(val)
	return verifyVals

def display2Decimal(val):
	newVal=float("{0:.2f}".format(val))
	return newVal
	
def finalMessage(verifyVals):
#Report three checks:
#1) Is the Spectrum score x^2 > 1? Yes: PASS. No: FAIL: Warn that match between high-scoring peptides and spectra is unusually low 
#2) Is the Tag length x^2 > 1? Yes: PASS. No: FAIL: Warn that unexpectedly few long de novo alignments with database 
#3) Is the mod size x^2 > 1? Yes: FAIL. Warn that Tag Graph had to insert far more large modifications that would be expected. No: PASS.
	result=""
	if (verifyVals[0]>1):
		result="Spectrum score test: result=%s : Pass.\n"%display2Decimal(verifyVals[0])
	else:
		result="Spectrum score test: result=%s : FAIL: match between high-scoring peptides and spectra is unusually low.\n"%display2Decimal(verifyVals[0])
	if (verifyVals[1]>1):
		result+="Tag length test: result=%s : Pass.\n"%display2Decimal(verifyVals[1])
	else:
		result+="Tag length test: result=%s : FAIL: unexpectedly few long de novo alignments with database.\n"%display2Decimal(verifyVals[1])
	if (verifyVals[2]>1):
		result+="Mod Size test: result=%s : FAIL: Tag Graph had to insert far more large modifications that would be expected.\n"% display2Decimal(verifyVals[2])
	else:
		result+="Mod Size test: result=%s : Pass.\n"%display2Decimal(verifyVals[2])
	return result
	
def verifyEM(outputFolder):
 	dbFile=outputFolder+"/results.db"
	conn=lite.connect(dbFile) 
	cursor=conn.cursor()
	createEMFullTable(conn)
	populateEMTable(outputFolder,conn)
	(EMHighCount,EMLowCount)=getEMHighLowCnt(cursor)
	(countsHigh,countsLow)=testCases(cursor)
	(highRatio_High,lowRatio_High,highRatio_Low,lowRatio_Low)=calculation(EMHighCount,EMLowCount,countsHigh,countsLow)
	verifyVals=calculation_Verify(highRatio_High,lowRatio_Low)
	result=finalMessage(verifyVals)
	return result

if __name__ == '__main__':
	if (len(sys.argv)<2):
		print " Error! \n Usage: python verifyEM.py [outputFolder]  \n"
		print " eg: python verifyEM.py /lab/samba/shared/Users/Sarah/A375 \n"
		# python verifyEM.py /lab/samba/shared/Users/Sarah/A375
		sys.exit(2)
	outputFolder=sys.argv[1]
	result=verifyEM(outputFolder)
	print result
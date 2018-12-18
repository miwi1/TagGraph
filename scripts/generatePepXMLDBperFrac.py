#slin 201707  added lib path and others.
#slin 20170921 Issue #12

import os
import sys
import glob
import pickle
from math import floor
import sqlite3 as lite
PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)

import DataFile
import ArgLib
import re
import xml.dom
import xml.etree.ElementTree as ET
from xml.etree import ElementTree
import time, datetime
from collections import defaultdict
import collections
from decimal import Decimal

AAMASS = {
'A':71.03711,
'B':114.534935,
'C':103.00919,
'D':115.02694,
'E':129.04259,
'F':147.06841,
'G':57.02146,
'H':137.05891,
'I':113.08406,
'J':0,
'K':128.09496,
'L':113.08406,
'M':131.04049,
'N':114.04293,
'O':100.076239,
'P':97.05276,
'Q':128.05858,
'R':156.10111,
'S':87.03203,
'T':101.04768,
'U':150.9536351,
'V':99.06841,
'W':186.07931,
'X':113.08406,
'Y':163.06333,
'Z':128.5505853
}

PROTONMASS=1.007823
H2O=18.01056469
PROTON=1.00727647
CConstant=57.021464
OXIDIZEDMETHIONINE=15.994915
maxModLength=30
MODIFIERCNT=4
NONEWORD=' \t\n\r'

AA = {
'Ala':'A',
'Arg':'R',
'Asn':'N',
'Asp':'D',
'Cys':'C',
'Gln':'Q',
'Glu':'E',
'Gly':'G',
'His':'H',
'Ile':'I',
'Leu':'L',
'Lys':'K',
'Met':'M',
'Phe':'F',
'Pro':'P',
'Ser':'S',
'Thr':'T',
'Trp':'W',
'Tyr':'Y',
'Val':'V',
'Sec':'?',
'Xle':'?',
'Camcys':'?',
'Pyro-glu':'?',
'Ethylaminoala':'?',
'Methylaminoala':'?',
'Aha':'?',
'Lacticacid':'?',
'Metox':'?',
'Pyruvicacid':'?',
'Glusa':'?',
'Oxoalanine':'?',
'Oxolactone':'?',
'Hydroxykynurenin':'?',
'Hpg':'?',
'Dha':'?',
'Npo':'?',
'Pyrrolidone':'?',
'Pyrrolidinone':'?',
'Orn':'?',
'Allysine':'?',
'Kynurenin':'?',
'Aminoadipicacid':'?',
'Hse':'?',
'Thiocarboxy':'?',
'Carboxy':'?',
'Hsl':'?'
}

def writeNodeToXML(rootNode,outFile):
	output_file=open(outFile,'w' )
	output_file.write('<?xml version="1.0"?>\n')
	rootNode=indent(rootNode)
	str=xml.etree.ElementTree.tostring(rootNode)
	output_file.write(str)
	output_file.close()

def writeRootNode(tag,dict):
	rootNode=ET.Element(tag)
	for item in dict:
		rootNode.set(item,dict[item])
	return rootNode

def writeInternalNode(node,tag,dict):
	internalNode=ET.SubElement(node,tag)
	for item in dict:
		internalNode.set(item,dict[item])
	return internalNode

def writeLeafNode(node,tag,dict):
	leafNode=ET.SubElement(node,tag)
	for item in dict:
		leafNode.set(item,dict[item])

def writeElement(node,elementName,listDictElement):
	for i in range(len(listDictElement)):
		dicElement = listDictElement[i]
		writeLeafNode(node,elementName,dicElement)

def insertStrAtPos(orgString,pos,subStr):
	 return orgString[:pos]+subStr+orgString[pos:] 

def writeAminoacidMod(node,aminoacidsDiffsSorted):
	if (len(aminoacidsDiffsSorted) > 0):
		for aminoacidInfo in aminoacidsDiffsSorted:
			splitAminoacidInfo=aminoacidInfo.split(":")
			aminoacid=splitAminoacidInfo[0]
			massdiff=splitAminoacidInfo[1]
			modMass=splitAminoacidInfo[2]
			description=aminoacid+splitAminoacidInfo[4]
			dictAminoacidModification={"aminoacid":aminoacid,"mass":str("{0:.6f}".format(round(float(modMass),6))),"massdiff":str("{0:.6f}".format(round(float(massdiff),6))),"variable":'Y',"description":description}
			writeLeafNode(node,"aminoacid_modification",dictAminoacidModification)

def writeCNTermsList(node,cnTermsList):
	if (len(cnTermsList) > 0):
		for CNTerms in cnTermsList:
			if CNTerms!="":
				splitCNTerms=CNTerms.split(":")
				CNTerm=splitCNTerms[0]
				massdiff=splitCNTerms[1]
				modMass=massdiff
				dictCNTerms={"terminus":CNTerm,"mass":str("{0:.6f}".format(round(float(modMass),6))),"massdiff":str("{0:.6f}".format(round(float(massdiff),6))),"variable":'Y'}
				writeLeafNode(node,"terminal_modification",dictCNTerms)

def writeAlternativeProtein(node,splitProtein):
	for i in xrange(len(splitProtein)):
		if (splitProtein[i].strip(NONEWORD)!=""):
			if (i==0):
				continue;
			altProtein=splitProtein[i].replace(" '","").replace("'","").split(" ")
			altProteinName=altProtein[0]
			altProteinDescription=splitProtein[i].replace(altProteinName+" ","").replace(" '","")
			dictAlternativeProtein={"protein":altProteinName,"protein_descr":altProteinDescription}
			writeLeafNode(node,"alternative_protein",dictAlternativeProtein)

########################
## obtained online.
## no xml format package, so use indent function to format xml better.
def indent(elem, level=0):
	i="\n" + level*"  "
	j="\n" + (level-1)*"  "
	if len(elem):
		if not elem.text or not elem.text.strip():
			elem.text=i + "  "
		if not elem.tail or not elem.tail.strip():
			elem.tail=i
		for subelem in elem:
			indent(subelem, level+1)
		if not elem.tail or not elem.tail.strip():
			elem.tail=j
	else:
		if level and (not elem.tail or not elem.tail.strip()):
			elem.tail=j
	return elem

def intLenth2Str(input):
	if (input<10):
		return "0"+str(input)
	else:
		return str(input) 

def modsArray(curMods,massArray):
	if (curMods.find(",",0)>-1):
		splitCurMods = curMods.split("), ((")
		for j in xrange(len(splitCurMods)):
			mods_ = splitCurMods[j].replace("\"", "").replace("(", "").replace(")", "").replace("[", "").replace("]", "").replace("\'", "")
			eachItem= mods_.split(",")
			curName=eachItem[0].strip().lower()
			curMass=eachItem[1].strip()
			lenCurMass=len(curMass)
			if (massArray.has_key(curName)):
				for k in xrange(len(massArray[curName])):
					if abs(float(massArray[curName][k])-float(curMass))>1:
						massArray[curName].append(curMass)
						break
					else:
						if (len(massArray[curName][k])<len(curMass)):
							massArray[curName][k]=curMass
							break
			else:
				massArray[curName].append(curMass)
	return massArray

def readTopResultFileModsArray(infile,delimiter=' '):
	mods=[]
	rowCnt=0
	header=""
	with open(infile) as f:
	   for line in f:
			line=line.strip()
			if (rowCnt > 0):
				splitLine = line.split(delimiter)
				mods.append(splitLine[16])
			if (rowCnt==0):
				header=line
				rowCnt=rowCnt+1
	massArray=defaultdict(list)
	for i in xrange(len(mods)):
		massArray=modsArray(mods[i],massArray)
	return (massArray)

def modList(splitLine,modStartPos,lineLen,modsCnt,inc):
	subMod=[]
	pad1=0
	if ((lineLen-(modStartPos+1))%MODIFIERCNT ==0):
		pad1=0
	else:
		if int((floor((lineLen-(modStartPos+1))/MODIFIERCNT) *MODIFIERCNT) + inc) >= (lineLen-modStartPos):
			pad1=0
		else:
			pad1=1 
	curModsCnt=((lineLen-(modStartPos+1))/MODIFIERCNT)+pad1
	for i in xrange(curModsCnt):
		subMod.append(splitLine[modStartPos+(i*MODIFIERCNT)+inc])
	for i in xrange(curModsCnt,modsCnt+1):
		subMod.append("")
	return subMod

def getAminoAcids(modMassArray,oneMod,cnTermsL,aminoAcidL):
	nTermMass=""
	cTermMass=""
	nTermMod=""
	cTermMod=""
	modfifyInfo=""
	strMod=re.sub('[\[\]]','',oneMod)
	mod=re.sub('\), \(\(', ')|((',str(strMod))
	splitMods=mod.split("|")
	nTerm=""
	cTerm=""
	for i in xrange(len(splitMods)):
		curMods=splitMods[i].replace(",","|").replace("(","").replace(")","").replace("'","")
		splitCurMods=curMods.split("|")
		if (len(splitCurMods)>4):
			element=splitCurMods[3].strip(NONEWORD)
			massdiff=splitCurMods[1].strip(NONEWORD)
			if ((splitCurMods[3].strip()=="C-term") or (splitCurMods[3].strip()=="N-term")):
				if (splitCurMods[3].strip()=="C-term"):
					cTermMass=massdiff
					cTermMod=splitCurMods[0]
					cTerm="c:"+str(massdiff)
				if (splitCurMods[3].strip()=="N-term"):
					nTermMass=massdiff
					nTermMod=splitCurMods[0]
					nTerm="n:"+str(massdiff)
				if (cTerm in cnTermsL):
					continue
				else:
					cnTermsL.append(cTerm)
				if (nTerm in cnTermsL):
					continue
				else: 
					cnTermsL.append(nTerm)
			if ((len(splitCurMods)==6) and (splitCurMods[3].strip()!="C-term") and (splitCurMods[3].strip()!="N-term")):
				element=splitCurMods[3].strip(NONEWORD)
				curName=splitCurMods[0].strip(NONEWORD)
				massdiff=splitCurMods[1].strip(NONEWORD)
				massUsed=massdiff
				baseMass=AAMASS[element]
				for k in xrange(len(modMassArray[curName.lower()])):
					if abs(float(modMassArray[curName.lower()][k])-float(massdiff))<0.1:
						massUsed=str(round(Decimal(modMassArray[curName.lower()][k]),6))
						break
				modMass=float(massUsed)+float(baseMass)
				aminoacids=element+":"+str(massUsed)+":"+str(modMass)+":"+str(baseMass)+":"+"(+"+str(round(float(massUsed),2))+")"
				modfifyInfo=modfifyInfo+"|"+aminoacids+":"+str(splitCurMods[5].strip(NONEWORD))
				if (aminoacids in aminoAcidL):
					continue
				else:
					aminoAcidL.append(aminoacids)
	return (nTermMass,nTermMod,cTermMass,cTermMod,cnTermsL,aminoAcidL)

def writeSpectrumQuery(node,scanL,chargeL,retentionTimeL,obsMassL,emProbL,oneSubtractlg10EML,spectrumScoreL,alignmentScoreL,compositeScoreL,contextL,modContextL,modsL,proteinsL,deNovoPeptideL,deNovoScoreL,ntermAAL,finalPeptideL,ctermAAL,recalculteMassL,newPPML,modPosL,modMassL,modNameL,modNTermMassL,modNTermModL,modCTermMassL,modCTermModL,cnTermsMassL,aminoAcidL):
	totalProteinCnt=0
	proteinName=""
	proteinDescription=""
	altProteinName=""
	altProteinDescription=""
	for i in xrange(len(scanL)):
		curID=i
		curScan=scanL[i]
		context=contextL[i]
		rank=1
		assumedCharge=chargeL[i]
		proteins=re.sub('[\[\-\]]','', proteinsL[i])
		deNovoPeptide=deNovoPeptideL[i]
		mods=modsL[i]
		modContext=modContextL[i]
		precursorNeutralMass=float(obsMassL[i])-PROTONMASS
		calcNeutralPepMass=float(recalculteMassL[i])-PROTONMASS
		massdiff=float(obsMassL[i])-float(recalculteMassL[i])
		peptide=contextL[i][2:-2]
		modPos=modPosL[i]
		modMass=modMassL[i]
		modName=modNameL[i]
		modNTermMass=modNTermMassL[i]
		modNTermMod=modNTermModL[i]
		modCTermMass=modCTermMassL[i]
		modCTermMod=modCTermModL[i]
		finalPeptide=finalPeptideL[i]
		splitProtein=proteins.split("',")
		totalProteinCnt=len(splitProtein)
		listDictScores=[]
		listDictScores.append({"name":"Alignment Score","value":str(alignmentScoreL[i])})
		listDictScores.append({"name":"Spectrum Score","value":str(spectrumScoreL[i])})
		listDictScores.append({"name":"Composite Score","value":str(compositeScoreL[i])})
		listDictScores.append({"name":"De Novo Score","value":str(deNovoScoreL[i])})
		spectrumName="spectrum"+str(curScan)
		dictspectrumquery={"spectrum":spectrumName,"start_scan":str(curScan),"end_scan":str(curScan),"precursor_neutral_mass":str(precursorNeutralMass),"assumed_charge":str(assumedCharge),"index":str(curScan)}
		spectrum_query=writeInternalNode(node, "spectrum_query",dictspectrumquery)
		for i in xrange(len(splitProtein)):
			if (splitProtein[i].strip(NONEWORD)!=""):
				if (i==0):
					protein=splitProtein[i].replace("'","").split(" ")
					proteinName=protein[0]
					proteinDescription=splitProtein[i].replace(proteinName+" ","")
					break
		search_result=ET.SubElement(spectrum_query,"search_result")
		dictSearchHit={"hit_rank":str(rank),"peptide":peptide,"denovo_peptide":deNovoPeptide,"peptide_prev_aa":context[0],"peptide_next_aa":context[-1],"protein":proteinName,"Final_Peptide":finalPeptide,"num_tot_proteins":str(totalProteinCnt),"calc_neutral_pep_mass":str(calcNeutralPepMass),"massdiff":str(massdiff),"protein_descr":proteinDescription}
		search_hit=writeInternalNode(search_result,"search_hit",dictSearchHit)
		writeElement(search_hit,"search_score",listDictScores)
		if (totalProteinCnt >1):
			writeAlternativeProtein(search_hit,splitProtein)
		writeModificationInfoFromFile(search_hit,modPos,modMass,modName,finalPeptide,modNTermMass,modNTermMod,modCTermMass,modCTermMod)

def writeModificationInfoFromFile(node,modPos,modMass,modName,finalPeptide,modNTermMass,modNTermMod,modCTermMass,modCTermMod):
	dictModifiedPeptide={"modified_peptide":finalPeptide}
	if modNTermMass!="":dictModifiedPeptide["mod_nterm_mass"]=modNTermMass
	if modCTermMass!="":dictModifiedPeptide["mod_cterm_mass"]=modCTermMass
	if (len(modName)>0):
		if modName[0]!="":
			modification_info=writeInternalNode(node,"modification_info",dictModifiedPeptide)
			writeModAminoAcidMass(modification_info,modPos,modMass,modName,modNTermMod,modCTermMod,len(re.findall('[a-zA-Z]',finalPeptide)))
		
def writeModAminoAcidMass(node,modPos,modMass,modName,modNTermMod,modCTermMod,peptideLen):
	dictModAminoAcidMass={}
	skipmass=0
	for i in range(len(modPos)):
		cnTerm=0
		if modPos[i]!="":
			if int(modPos[i])==0 and modName[i].lower()==modNTermMod.lower():
				cnTerm=1
			if int(modPos[i])==peptideLen-1 and modName[i].lower()==modCTermMod.lower():
				cnTerm=1
			if cnTerm==0:
				if (modName[i].find("->",0)>=0):
					splitLine=modName[i].split("->")
					newChar=AA[splitLine[1].capitalize()]
					if (newChar!="?"):
						skipmass+=1
						dictModAminoAcidMass={"modname":modName[i],"position":str(modPos[i])}
						writeLeafNode(node, "mod_aminoacid_mass",dictModAminoAcidMass)
					else:	
						dictModAminoAcidMass={"modname":modName[i],"position":str(modPos[i]),"mass":str(modMass[i-skipmass])}
						writeLeafNode(node, "mod_aminoacid_mass",dictModAminoAcidMass)	
				else:
					dictModAminoAcidMass={"modname":modName[i],"position":str(modPos[i]),"mass":str(modMass[i-skipmass])}
					writeLeafNode(node, "mod_aminoacid_mass",dictModAminoAcidMass)	
		else:
			break
			
def	getModMassPerFraction(conn,fractionID):
	modsList=[]
	cursor=conn.cursor()
	query="select mods from result where top_result=1 and fraction_id="+str(fractionID)+" order by scan"
	cursor.execute(query)
	scans=cursor.fetchall()	
	for scan in scans:
		modsList.append(scan[0])
	massArray=defaultdict(list)
	for i in xrange(len(modsList)):
		massArray=modsArray(modsList[i],massArray)
	return (massArray)
	
def	getScansPerFraction(conn,fractionID,modMassArray,FDRCutoff,logEMCutoff):
	idL=[]
	scanL=[]
	chargeL=[]
	obsMassL=[]
	retentionTimeL=[]
	alignmentScoreL=[]
	spectrumScoreL=[]
	compositeScoreL=[]
	emProbL=[]
	oneSubtractlg10EML=[]
	contextL=[]
	modContextL=[]
	modsL=[]
	proteinsL=[]
	deNovoPeptideL=[]
	deNovoScoreL=[]
	theomassL=[]
	ppmL=[]
	ntermAAL=[]
	finalPeptideL=[]
	ctermAAL=[]
	aminoAcidL=[]
	modNTermMassL=[]
	modNTermModL=[]
	modCTermMassL=[]
	modCTermModL=[]	
	cnTermsL=[]
	modPosL=[]
	modMassL=[]
	modNameL=[]
	modAAL=[]
	cursor=conn.cursor()
	query="select r.id,r.scan,r.charge,r.obs_mh,r.retention_time,r.alignment_score,r.spectrum_score,r.composite_score,r.em_probability,r.log_em_probability,r.context,r.mod_context,r.mods,r.proteins,r.de_novo_peptide,r.de_novo_score,f.newTheoMass,f.newPPm,f.nTerm,f.finalPeptide,f.cTerm,c.FDR "
	for i in range(0,maxModLength):
		query+=",m.modPos"+str(i)
		query+=",m.modName"+str(i)
		query+=",m.modAA"+str(i)
		query+=",m.modMass"+str(i)
	query+=" from result r left join finalPeptidePPMInfo f on r.id=f.id left join cleanDataFDR c on r.id=c.id left join modsPerScan m on r.id=m.id where r.top_result=1 and r.fraction_id="+str(fractionID)+" and r.id=f.id and r.fraction_id=f.fractionID and r.id=m.id and r.fraction_id=m.fractionID order by r.scan"
	cursor.execute(query)
	scans=cursor.fetchall()
	for curScan in scans:
		floatFDRCutoff=float(FDRCutoff)
		floatlogEMCutoff=float(logEMCutoff)
		if (float(curScan[21])<=floatFDRCutoff or float(curScan[9])>=floatlogEMCutoff):	
			idL.append(curScan[0])
			scanL.append(curScan[1])
			chargeL.append(curScan[2])
			obsMassL.append(curScan[3])
			retentionTimeL.append(curScan[4])
			alignmentScoreL.append(curScan[5])
			spectrumScoreL.append(curScan[6])
			compositeScoreL.append(curScan[7])
			emProbL.append(curScan[8])
			oneSubtractlg10EML.append(curScan[9])
			contextL.append(curScan[10])
			modContextL.append(curScan[11])
			curMods=curScan[12]
			modsL.append(curMods)
			proteinsL.append(curScan[13])
			deNovoPeptideL.append(curScan[14])
			deNovoScoreL.append(curScan[15])
			theomassL.append(curScan[16])
			ppmL.append(curScan[17])
			ntermAAL.append(curScan[18])
			finalPeptideL.append(curScan[19])
			ctermAAL.append(curScan[20])
			(nTermMass,nTermMod,cTermMass,cTermMod,cnTermsL,aminoAcidL)=getAminoAcids(modMassArray,curMods,cnTermsL,aminoAcidL)
			modNTermMassL.append(nTermMass)
			modNTermModL.append(nTermMod)
			modCTermMassL.append(cTermMass)
			modCTermModL.append(cTermMod)
			submodPosL=[]
			submodNameL=[]
			submodAAL=[]
			submodMassL=[]
			startPos=22
			for i in range(0,maxModLength):
				if (curScan[startPos+(i*4)]!=''):
					submodPosL.append(curScan[startPos+(i*4)])
					submodNameL.append(curScan[1+startPos+(i*4)])
					submodAAL.append(curScan[2+startPos+(i*4)])
					submodMassL.append(curScan[3+startPos+(i*4)])
				else:
					break
			modPosL.append(submodPosL)
			modNameL.append(submodNameL)
			modAAL.append(submodAAL)
			modMassL.append(submodMassL)
	return (scanL,chargeL,retentionTimeL,obsMassL,emProbL,oneSubtractlg10EML,spectrumScoreL,alignmentScoreL,compositeScoreL,contextL,modContextL,modsL,proteinsL
	,deNovoPeptideL,deNovoScoreL,ntermAAL,finalPeptideL,ctermAAL,theomassL,ppmL,modPosL,modMassL,modNameL,modNTermMassL,modNTermModL,modCTermMassL,modCTermModL,sorted(cnTermsL),sorted(aminoAcidL))

if __name__=='__main__':
	#eg: python generatePepXMLDBperFrac.py /home/slin/sampleInputFiles/TG/TAG_GRAPH_Tryp_CysCarbam_MetOx.ini 10 0.1 400 200 /home/slin/sampleInputFiles/FMIndices/20141209_UniHUMAN_cRAP_ILEq.fm resources/AllChargeDist_posOnlyDependence_20150808_HumanProt500000.pck resources/unimodDict_noLabels_20160724.pck /home/slin/sampleInputFiles/TG/mzXML_EM_output/ /tmp/A375mzXML_28981/ 0.01 20
	print "Start Time:", time.asctime(time.localtime(time.time()))
	arg_init=sys.argv[1]
	arg_ppmstd=sys.argv[2]
	arg_modtolerance=sys.argv[3]
	arg_maxcounts=sys.argv[4]
	arg_modmaxcounts=sys.argv[5]
	arg_fmindex=sys.argv[6]
	arg_model=sys.argv[7]
	arg_unimoddict=sys.argv[8]
	arg_output=sys.argv[9]
	arg_dtadir=sys.argv[10]
	FDRCutoff=sys.argv[11]
	logEMCutoff=sys.argv[12]
	outDir=os.path.join(arg_dtadir,'taggraph')
	peaksDir=os.path.join(arg_dtadir,'de_novo')
	dataDir=os.path.join(arg_dtadir,'data')
	summaryFileExtension=""
	summaryFileBase=""
	for file in os.listdir(arg_output): #get latest final output file.
		if file.find("results.db")==0:
			fileFound=1
			finaloutputFile=file
			break;
	if fileFound:
		fracMappingFile=open(arg_output+'/fileFractionMapping.pck','rb')
		fileFractionMapping=pickle.load(fracMappingFile)
		fractionMapping={}
		for i in range(0,len(fileFractionMapping)):
			fractionMapping[fileFractionMapping[i][0]]=fileFractionMapping[i][1]
		resultsDBFile=os.path.join(arg_output,"results.db")
		infile=os.path.join(arg_output,finaloutputFile)
		conn=lite.connect(resultsDBFile)
		conn.execute("PRAGMA max_page_count=max_page;")
		conn.execute("PRAGMA temp_store=2;")
		cursor=conn.cursor()
		listDictInputParams=[]
		listDictInputParams.append({"name":"init","value":str(arg_init)})
		listDictInputParams.append({"name":"ppmstd","value":str(arg_ppmstd)})
		listDictInputParams.append({"name":"modtolerance","value":str(arg_modtolerance)})
		listDictInputParams.append({"name":"model","value":str(arg_model)})
		listDictInputParams.append({"name":"fmindex","value":str(arg_fmindex)})
		listDictInputParams.append({"name":"unimoddict","value":str(arg_unimoddict)})
		listDictInputParams.append({"name":"maxcounts","value":str(arg_maxcounts)})
		listDictInputParams.append({"name":"modmaxcounts","value":str(arg_modmaxcounts)})
		listDictInputParams.append({"name":"dtadir","value":str(arg_dtadir)})  
		(tm_year,tm_mon,tm_mday,tm_hour,tm_min,tm_sec,tm_wday,tm_yday,tm_isdst)=time.localtime(os.path.getmtime(resultsDBFile))
		fileCreationTime=str(tm_year)+"-"+intLenth2Str(tm_mon)+"-"+intLenth2Str(tm_mday)+"T"+intLenth2Str(tm_hour)+":"+intLenth2Str(tm_min)+":"+intLenth2Str(tm_sec)
		query="select distinct(fraction_id) from result"
		cursor.execute(query)
		fractionIDs=cursor.fetchall()
		for fractID in fractionIDs:
			fractionID=fractID[0]
			mzMLmxXMLFile=fractionMapping[intLenth2Str(fractionID)]
			if mzMLmxXMLFile.lower().endswith(".mzxml"):
				summaryFileExtension=mzMLmxXMLFile[-6:]
				summaryFileBase==dataDir+"/"+mzMLmxXMLFile[0:-6]
			if mzMLmxXMLFile.lower().endswith(".mzml"):
				summaryFileExtension=mzMLmxXMLFile[-5:]
				summaryFileBase=dataDir+"/"+mzMLmxXMLFile[0:-5]
			tagGraphpepXMLName="TagGraph.Fraction"+intLenth2Str(fractionID)+".pep.xml"
			tagGraphpepXML=arg_output+"/"+tagGraphpepXMLName
			modMassArray=getModMassPerFraction(conn,fractionID)
			(scanL,chargeL,retentionTimeL,obsMassL,emProbL,oneSubtractlg10EML,spectrumScoreL,alignmentScoreL,compositeScoreL,contextL,modContextL,modsL,proteinsL,deNovoPeptideL,deNovoScoreL,ntermL,finalPeptideL,ctermL,theoMassL,ppmL,modPosL,modMassL,modNameL,modNTermMassL,modNTermModL,modCTermMassL,modCTermModL,cnTermsL,aminoAcidL)=getScansPerFraction(conn,fractionID,modMassArray,FDRCutoff,logEMCutoff)
			dictPipelineAnalysis={"date":fileCreationTime,"summary_xml":tagGraphpepXMLName}
			msms_pipeline_analysis=writeRootNode("msms_pipeline_analysis", dictPipelineAnalysis)
			dictMsmsRunSummary={"base_name":summaryFileBase,"raw_data_type":"raw","raw_data":summaryFileExtension}
			msms_run_summary=writeInternalNode(msms_pipeline_analysis, "msms_run_summary",dictMsmsRunSummary)
			dictSampleEnzyme={"name":"None","fidelity":"nonspecific"}
			sample_Enzyme=writeInternalNode(msms_run_summary, "sample_enzyme",dictSampleEnzyme)
			dictSpecificity={"cut":"ARNDCEQGHILKMFPSTWYVBJOXZ","sense":"C"}
			writeLeafNode(sample_Enzyme,"specificity",dictSpecificity)
			dictSearchSummary={"base_name":summaryFileBase,"search_engine":"SEQUEST","precursor_mass_type":"monoisotopic","fragment_mass_type":"monoisotopic","search_id":"1"}
			search_summary=writeInternalNode(msms_run_summary, "search_summary",dictSearchSummary)
			dictSearchDatabase={"type":"AA"}
			writeLeafNode(search_summary, "search_database",dictSearchDatabase)
			dictEnzymaticSearchConstraint={"enzyme":"None","max_num_internal_cleavages":"100","min_number_termini":"2"}
			enzymatic_search_constraint=writeInternalNode(search_summary, "enzymatic_search_constraint",dictEnzymaticSearchConstraint)
			writeAminoacidMod(search_summary,aminoAcidL)
			writeCNTermsList(search_summary,cnTermsL)
			writeElement(search_summary,"parameter",listDictInputParams)
			writeSpectrumQuery(msms_run_summary,scanL,chargeL,retentionTimeL,obsMassL,emProbL,oneSubtractlg10EML,spectrumScoreL,alignmentScoreL,compositeScoreL,contextL,modContextL,modsL,proteinsL,deNovoPeptideL,deNovoScoreL,ntermL,finalPeptideL,ctermL,theoMassL,ppmL,modPosL,modMassL,modNameL,modNTermMassL,modNTermModL,modCTermMassL,modCTermModL,cnTermsL,aminoAcidL)
			writeNodeToXML(msms_pipeline_analysis,tagGraphpepXML)
	else:
		print "pepxml output process failed"
	print "End Time :",time.asctime(time.localtime(time.time()))
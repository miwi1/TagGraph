# slin 201707xx modify using unimod: unimodDict.pck
#               calculate FDR for all but not just with the scan in the ppm range.
#               use .ini file defined in the parameter file to obtain the static mods and other adjustments.
#               filter out the scan when  FDR<=FDRCutoff or logEM>=logEMCutoff:
# slin 20170817 fix the staticModsRemove list postion when insertion happen
# slin 20170921 Issue #11, and #12
# slin 20171028 Issue #14
# slin 20171104 Issue #13
# slin 20171213 issue #17

import csv
import os
import re
import sys
import datetime
import copy
import string
from decimal import Decimal
from warnings import filterwarnings
from collections import defaultdict
import collections
import pickle
import math
import sqlite3 as lite
import re
PAR_DIR=os.path.abspath(os.path.join(os.path.dirname(__file__),os.path.pardir))
sys.path.insert(1,PAR_DIR)
LIB_DIR=PAR_DIR+"/lib"
sys.path.insert(1,LIB_DIR)
DATABASE_SCRIPT_DIR=os.path.join(PAR_DIR,'database')
import ArgLib
import DataFile
AAMass={
'A':71.03711,
'B':115.02694,
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
'Z':129.04259
}
AMINOACID=['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val']
AA={
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
H2O=18.0105647
PROTON=1.00727647
PPMRANGE=50.0
PPMITERATIONS=20
PRECISION=6
maxModLength=30
MaxProtenLengthEach=255
def modText(eachMod,newStr,newStrAAOnly,cpstaticModsListOne,staticModsCounter,NTERM,CTERM,adjustedMass,staticModsRemove,shiftItem,dictStaticMods):
	replaceChar=''
	adjustShift=0
	if ((len(eachMod)==6) or (len(eachMod)==5)):
		replace0=eachMod[0].strip()
		diffMass=round(Decimal(eachMod[1].strip()),2)
		if len(eachMod)==6:
			replaceChar=eachMod[3].strip()
			toBeReplaced=eachMod[4].strip()
 			replacePos=int(eachMod[5].strip())
		else:
			toBeReplaced=eachMod[3].strip()
			replacePos=int(eachMod[4].strip())
		origChar=newStrAAOnly[replacePos+2:replacePos+3]
		if shiftItem:  #adjustPos: #issue#1
			for i in range(len(shiftItem)):
				adjustShift=adjustShift+shiftItem[i][1]
		if (replace0.find("->",0)>=0):  #directReplace
			splitLine=replace0.split("->")
			newChar=AA[splitLine[1].capitalize()]
			for aa in cpstaticModsListOne:  #exit if all empty.
				if (len(cpstaticModsListOne[aa]))>0:
					if (AA[splitLine[0].capitalize()]==aa) and (newChar!="?"):
						cpstaticModsListOne[aa].remove(replacePos)
						break
			for aa in dictStaticMods:  #exit if all empty.
				if (newChar==aa):
					staticModsRemove[aa].append(replacePos+adjustShift)
					staticModsCounter[aa]=staticModsCounter[aa]+1
					break
			if (newChar=="?"): #no matching, need the original Char, and diffMass
				if (replacePos==0) and (newStr[2]=="["): #check if the new string has been modified via n-terms
					newReplacePos=newStr.find(']')
					newStr=replaceAtPos(newReplacePos+1,newReplacePos+2,newStr,origChar+"["+str(diffMass)+"]")
				else:
					newStr=replaceAtPos(replacePos+2,replacePos+3,newStr,origChar+"["+str(diffMass)+"]")
				newStrAAOnly=replaceAtPos(replacePos+2,replacePos+3,newStrAAOnly,origChar)
			else:
				if (replacePos==0) and (newStr[2]=="["): #if it got modified by n-term, really need to find the first char pos, replace it!
					newReplacePos=newStr.find(']')
					newStr=replaceAtPos(newReplacePos+1,newReplacePos+2,newStr,newChar)
				else:
					newStr=replaceAtPos(replacePos+2,replacePos+3,newStr,newChar)
				newStrAAOnly=replaceAtPos(replacePos+2,replacePos+3,newStrAAOnly,newChar)
		else:
			if (replace0=='Insertion'):
				insertLen=len(toBeReplaced)
				shiftItem.append([replacePos,insertLen])
				for aa in cpstaticModsListOne:
					curModsList=findChr(toBeReplaced,aa)
					for i in range(0,len(curModsList)):
						if (curModsList[i]):
							curModsList[i]=curModsList[i]+replacePos
					for i in range(0,len(cpstaticModsListOne[aa])):
						if (cpstaticModsListOne[aa][i]>=replacePos):
							cpstaticModsListOne[aa][i]=cpstaticModsListOne[aa][i]+insertLen+adjustShift
				for aa in dictStaticMods: #slin issue#1
					if (len(staticModsRemove[aa])>0):
						for i in range(len(staticModsRemove[aa])):
							if (staticModsRemove[aa][i]>=replacePos):
								staticModsRemove[aa][i]+=insertLen
				for aa in dictStaticMods:
					curModsList=findChr(toBeReplaced,aa)
					for i in range(len(curModsList)):
						staticModsRemove[aa].append(replacePos+curModsList[i])
				newStr=newStr[0:replacePos+2]+toBeReplaced+newStr[replacePos+2:]
				newStrAAOnly=newStrAAOnly[0:replacePos+2]+toBeReplaced+newStrAAOnly[replacePos+2:]
			elif (replace0=='Deletion'):
				delLen=len(toBeReplaced)
				shiftItem.append([replacePos,-delLen])
				for aa in cpstaticModsListOne: 
					delCList=findChr(toBeReplaced,aa)
					for i in range(len(delCList)):
						removePos=delCList[i]+replacePos
						cpstaticModsListOne[aa].remove(removePos) #delete the pos
					for i in range(len(cpstaticModsListOne[aa])): #shift the position
						if (cpstaticModsListOne[aa][i]>replacePos):
							cpstaticModsListOne[aa][i]=cpstaticModsListOne[aa][i]-delLen+adjustShift
					endReplacePos=replacePos+delLen #slin issue #1
				#print "sarah check staticModsRemove %s" % staticModsRemove
				for aa in staticModsRemove: #slin issue#1
					if (len(staticModsRemove[aa])>0):
						curModsList=findChr(toBeReplaced,aa)
						#sarah check len(staticModsRemove[aa]) 1 curModsList [8, 9, 11] toBeReplaced TGTGSGASCCPC replacePos 8 endReplacePos 20
						#print  "sarah check delLen %s len(staticModsRemove[aa]) %s curModsList %s toBeReplaced %s replacePos %s endReplacePos %s" % (delLen,len(staticModsRemove[aa]),curModsList,toBeReplaced,replacePos,endReplacePos)
						for i in range(len(staticModsRemove[aa])):
							if (staticModsRemove[aa][i]>=replacePos) and (staticModsRemove[aa][i]<endReplacePos):
								staticModsRemove[aa][i]=staticModsRemove[aa][i]-delLen
						for i in reversed(xrange(len(staticModsRemove[aa]))): #issue#1
							if (staticModsRemove[aa][i]>=replacePos and staticModsRemove[aa][i]<endReplacePos):
								del staticModsRemove[aa][i]
							if (staticModsRemove[aa][i]>endReplacePos):
								staticModsRemove[aa][i]=staticModsRemove[aa][i]-delLen+adjustShift
				delLength=len(toBeReplaced)-1
				newStr=newStr[0:replacePos+2]+newStr[replacePos+3+delLength:]
				newStrAAOnly=newStrAAOnly[0:replacePos+2]+newStrAAOnly[replacePos+3+delLength:]
			elif(replace0=='Isobaric Substitution'): #do nothing
				newStr=newStr
				newStrAAOnly=newStrAAOnly
			elif (len(replaceChar)==1):
				if (replacePos==0) and (newStr[2]=="["): #if it got modified by n-term, really need to find the first char pos, replace it!
					newReplacePos=newStr.find(']')
					newStr=replaceAtPos(newReplacePos+1,newReplacePos+2,newStr,replaceChar+"["+str(diffMass)+"]")
				else:
					newStr=replaceAtPos(replacePos+2,replacePos+3,newStr,replaceChar+"["+str(diffMass)+"]")
				newStrAAOnly=replaceAtPos(replacePos+2,replacePos+3,newStrAAOnly,replaceChar)
			else:
				if ((toBeReplaced=="N-term") and (replaceChar=="N-term")):
					diffMass+=NTERM
					newStr=replaceAtPos(replacePos+2,replacePos+3,newStr,"["+str(round(Decimal(diffMass),2))+"]"+origChar)
					newStrAAOnly=replaceAtPos(replacePos+2,replacePos+3,newStrAAOnly,origChar)
				else:
					if (replacePos==0) and (newStr[2]=="["): #if it got modified by nterm,really need to find the first char pos, replace it!
						newReplacePos=newStr.find(']')
						newStr=replaceAtPos(newReplacePos+1,newReplacePos+2,newStr,origChar+"["+str(diffMass)+"]")
					else:
						newStr=replaceAtPos(replacePos+2,replacePos+3,newStr,origChar+"["+str(diffMass)+"]")
					newStrAAOnly=replaceAtPos(replacePos+2,replacePos+3,newStrAAOnly,origChar)
	for aa in staticModsRemove:  #slin issue#1
		sorted(staticModsRemove[aa])
	return (newStr,newStrAAOnly,cpstaticModsListOne,staticModsCounter,adjustedMass,staticModsRemove,shiftItem)

def modsSubtitionPerScan(curMods,modContext,context,mPoundStr,cpstaticModsDict,NTERM,CTERM,dictDiffMods,dictDiffModsMass,dictStaticMods,dictStaticModsName,diffModsPos,diffModsStr,dictDiffModsWithOriginal,staticModsRemove,staticModsCounter):
	substituteInfo=[]
	mPoundModFinal=[]
	firstItem=[]
	lastItem=[]
	shiftItem=[]
	newStr=context
	newStrAAOnly=context
	mPoundCnt=0
	adjustedMass=0
	modNTerm=0
	modCTerm=0
	if (len(cpstaticModsDict))>=0:
		for aa in cpstaticModsDict:
			staticModsRemove[aa]=[]
			staticModsCounter[aa]=0
	if (curMods.find(",",0)==-1): #condition1: no mod exist
		contextWithAAOnly=modContext
		(modContext)=adjStaticMod(modContext,dictStaticMods,staticModsRemove)
		if (NTERM!=0):
			adjustedMass+=NTERM
			modContext=modContext[:2]+"["+str(round(Decimal(NTERM),2))+"]"+modContext[2:]
		if (CTERM!=0):
			adjustedMass+=NTERM
			modContext=modContext[:-2]+"["+str(round(Decimal(CTERM),2))+"]"+modContext[-2:]
	elif (curMods.find("), ((")==-1): #condition2: if only 1 mod exist
		fullMods=curMods.replace("\"","").replace("(","").replace(")","").replace("[","").replace("]","").replace("\'","")
		eachMod= fullMods.split(",")
		if len(eachMod)==6:
			replaceChar=eachMod[3].strip()
			if (replaceChar=="N-term"):
				modNTerm=1
			if (replaceChar=="C-term"):
				modCTerm=1
		(modContext,contextWithAAOnly,cpstaticModsDict,staticModsCounter,adjustedMass,staticModsRemove,shiftItem)=modText(eachMod,newStr,newStrAAOnly,cpstaticModsDict,staticModsCounter,NTERM,CTERM,adjustedMass,staticModsRemove,shiftItem,dictStaticMods)
		(modContext)=adjStaticMod(modContext,dictStaticMods,staticModsRemove)
		if (NTERM!=0 and modNTerm==0):
			(modContext,adjustedMass)=modNCTerm(modContext,NTERM,"N-term",adjustedMass)
		if (CTERM!=0 and modCTerm==0):
			(modContext,adjustedMass)=modNCTerm(modContext,CTERM,"C-term",adjustedMass)
	elif (curMods.find("), ((")>-1):
		splitCurMods=curMods.split("), ((")
		for j in reversed(xrange(len(splitCurMods))):  #reverse order.
			mods_=splitCurMods[j].replace("\"","").replace("(","").replace(")","").replace("[","").replace("]","").replace("\'","")
			eachMod=mods_.split(",")
			if len(eachMod)==6:
				replaceChar=eachMod[3].strip()
				if (replaceChar=="N-term"):
					modNTerm=1
				if (replaceChar=="C-term"):
					modCTerm=1
			(modContext,contextWithAAOnly,cpstaticModsDict,staticModsCounter,adjustedMass,staticModsRemove,shiftItem)=modText(eachMod,newStr,newStrAAOnly,cpstaticModsDict,staticModsCounter,NTERM,CTERM,adjustedMass,staticModsRemove,shiftItem,dictStaticMods)
			newStr=modContext
			newStrAAOnly=contextWithAAOnly
		(modContext)=adjStaticMod(modContext,dictStaticMods,staticModsRemove)
		if (NTERM!=0 and modNTerm==0):
			(modContext,adjustedMass)=modNCTerm(modContext,NTERM,"N-term",adjustedMass)
		if (CTERM!=0 and modCTerm==0):
			(modContext,adjustedMass)=modNCTerm(modContext,CTERM,"C-term",adjustedMass)
	modMPoundIDList=[]
	diffModMassChange=0
	for aa in dictDiffMods: #diffModsPos,diffModsStr
		(modMPoundIDList)=diffModPrep(modContext,contextWithAAOnly,diffModsStr[aa],aa,dictDiffMods[aa],diffModsPos[aa],dictDiffMods)
		(modContext,mPoundCnt,mPoundModFinal,diffModMassChange)=diffModProcess(modContext,contextWithAAOnly,diffModsStr[aa],modMPoundIDList,aa,dictDiffMods[aa],str(round(Decimal(dictDiffModsMass[aa],6),2)),diffModMassChange,dictStaticMods,dictDiffModsWithOriginal)
	return(modContext,contextWithAAOnly,mPoundCnt,mPoundModFinal,cpstaticModsDict,staticModsCounter,adjustedMass,diffModMassChange,staticModsRemove)

def adjStaticMod(newStr,dictStaticMods,staticModsRemove):
	subString=newStr[2:-2]
	for aa in dictStaticMods:
		modMass=str(round(float(dictStaticMods[aa]),2))
		modWithMass=aa+"["+modMass+"]"
		removeLength=len(modWithMass)-1 #subtract 1 for the AA
		subString=subString.replace(modWithMass,aa).replace(aa,modWithMass)
		for i in reversed(xrange(len(staticModsRemove[aa]))):  
			counter=-1
			cPos=staticModsRemove[aa][i]
			for j in xrange(len(subString)):
				if counter==cPos:
					subString=replaceAtPos(j,j+removeLength,subString,"")
					break
				if isAlpha(subString[j]):
					counter+=1
	finalStr=newStr[:2]+subString+newStr[-2:]
	return (finalStr)

def diffModPrep(modContext,newStrAAOnly,mPoundStr,aa,subSymbol,modMPoundIDList,dictDiffMods):
	padLeft=0
	padRight=0
	mPoundCnt=0
	mPoundStrPosList=findChr(mPoundStr,subSymbol)
	if (len(mPoundStrPosList)>1):
		for i in xrange(len(mPoundStrPosList)):
			for j in xrange(len(mPoundStrPosList)-1,i-1,-1):
				if mPoundStrPosList[i]==0:
					padLeft=0
				elif mPoundStrPosList[i]==1:
					padLeft=1
				elif mPoundStrPosList[i]==2:
					padLeft=2
				elif mPoundStrPosList[i]==3:
					padLeft=3
				else:
					padLeft=4
				if mPoundStrPosList[j]==len(mPoundStr):
					padRight=0
				elif mPoundStrPosList[j]==len(mPoundStr)+1:
					padRight=1
				elif mPoundStrPosList[j]==len(mPoundStr)+2:
					padRight=2
				elif mPoundStrPosList[j]==len(mPoundStr)+3:
					padRight=3
				else:
					padRight=4
				curStrMPound=mPoundStr[mPoundStrPosList[i]-padLeft:mPoundStrPosList[j]+padRight]
				curStr=curStrMPound.replace(subSymbol,"")
				foundMPoundPos=newStrAAOnly.find(curStr)
				if (foundMPoundPos>=0):
					for k in (i,j):
						modMPoundIDList.append(k)
	else:
		curStrMPound=mPoundStr
		curStr=curStrMPound.replace(subSymbol,"")
		foundMPoundPos=newStrAAOnly.find(curStr)
		if (foundMPoundPos>=0):
			modMPoundIDList.append(0)
	newList=list(set(modMPoundIDList))
	newList.sort()
	return (newList)

def diffModProcess(modContext,newStrAAOnly,mPoundStr,modMPoundIDList,aa,subSymbol,aaMass,diffModMassChange,dictStaticMods,dictDiffModsWithOriginal):
	#mPoundStr, the orginal greedy M# string
	mPoundStrPos=[]
	mPoundStrPos=findChr(mPoundStr,subSymbol)
	mPoundModFinal=[]
	mPoundCnt=0
	if (len(mPoundStrPos)>0):
		mStr=mPoundStr.replace(subSymbol,"")
		for l in xrange(len(mPoundStrPos)):
			mPoundStrPos[l]=mPoundStrPos[l]-l-1
		for i in xrange(len(mPoundStrPos)-1,-1,-1):
			if not i in modMPoundIDList:
				mPoundStrPos.pop(i)
		if (dictDiffModsWithOriginal[aa]=="1"):
			lookFor=aa+"["+ str(round(Decimal(dictStaticMods[aa]),2))+"]["+aaMass+"]"
		else:
			lookFor=aa+"["+aaMass+"]"
		for l in xrange(len(mPoundStrPos)):
			aaCount=0
			setbreak=0
			padLeft=0
			padshift=0
			if mPoundStrPos[l]==0:
				padLeft=0
			elif mPoundStrPos[l]==1:
				padLeft=1
			elif mPoundStrPos[l]==2:
				padLeft=2
			else:
				padLeft=3
			padRight=0
			if (mPoundStrPos[l]==len(mPoundStr)):
				padRight=0
			if (mPoundStrPos[l]+1<=len(mPoundStr)):
				padRight=1
			if (mPoundStrPos[l]+2<=len(mPoundStr)):
				padRight=2
			if (mPoundStrPos[l]+3<=len(mPoundStr)):
				padRight=3
			curSubStr=mStr[mPoundStrPos[l]-padLeft:mPoundStrPos[l]+padRight]
			foundMPoundPos=newStrAAOnly.find(curSubStr)
			for i in xrange(2,len(modContext)-1):
				if (modContext[i]>='A' and modContext[i]<='Z'):
					aaCount+=1
				if (aaCount==foundMPoundPos+padLeft-1) and (foundMPoundPos!=-1):
					if modContext[i]==aa.upper():
						subModContext=modContext[i:len(modContext)-1]
						if modContext[i+1]=="[":
							if (subModContext.find(lookFor)<>0):
								modContext=replaceAtPos(i-padshift,i+1-padshift,modContext,aa.upper()+"["+aaMass+"]")
								diffModMassChange=diffModMassChange+float(aaMass)
								mPoundCnt=mPoundCnt+1
								mPoundModFinal.append(foundMPoundPos+padLeft-2)
						else:
							modContext=replaceAtPos(i-padshift,i+1-padshift,modContext,aa.upper()+"["+aaMass+"]")
							diffModMassChange=diffModMassChange+float(aaMass)
							mPoundCnt=mPoundCnt+1
							mPoundModFinal.append(foundMPoundPos+padLeft-2)
		if (dictDiffModsWithOriginal[aa]=="1"): #remove the static mod, remove [xx], subtract mass
			toBeReplaced=aa+"["+aaMass+"]"
			cnt=modContext.count(lookFor)
			modContext=modContext.replace(lookFor,toBeReplaced)
			diffModMassChange=diffModMassChange-(float(dictStaticMods[aa])*cnt)
			lookFor=aa+"["+aaMass+"]["+str(round(Decimal(dictStaticMods[aa]),2))+"]"
			cnt=modContext.count(lookFor)
			modContext=modContext.replace(lookFor,toBeReplaced)
			diffModMassChange=diffModMassChange-(float(dictStaticMods[aa])*cnt)
	return (modContext,mPoundCnt,mPoundModFinal,diffModMassChange)

def splitModsSubtition(newDic,allMods,modContextList,contextList,maxModLength,deNovoPeptide,staticModsList,NTERM,CTERM,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictStaticMods,dictStaticModsName,dictDiffModsWithOriginal,dictStaticCNTermModsName):
	newContext=[]
	modsInfo=[]
	newTheoMassList=[]
	posOfMPound=[]
	curList=[]
	mPoundModFinalStr=[]
	massArray=defaultdict(list)
	mPoundStr=""
	mPoundModFinalList=[]
	newMassArray={}
	for i in xrange(len(allMods)):
		modsArray(allMods[i],massArray)
	newMassArray=cleanupMassArray(massArray)
	newDic=mergeToDictionary(newDic,newMassArray)
	diffModsPos={}
	diffModsStr={}
	for i in xrange(len(allMods)):
		mPoundStr=""
		theoMass=0
		origTheoMass=0
		mPoundCnt=0
		for aa in dictDiffMods:
			posOfMPound=[pos for pos,ltr in enumerate(deNovoPeptide[i]) if ltr==dictDiffMods[aa]]
			diffModsPos[aa]=posOfMPound
			padRight=0
			padLeft=0
			if (len(posOfMPound)>=1):
				if posOfMPound[0]==1:
					padLeft=1
				elif posOfMPound[0]==2:
					padLeft=2
				else:
					padLeft=3
				if (len(posOfMPound)>1):
					if posOfMPound[len(posOfMPound)-1]==len(deNovoPeptide[i]):
						padRight=0
					elif  posOfMPound[len(posOfMPound)-1]==(len(deNovoPeptide[i])-1):
						padRight=1
					elif  posOfMPound[len(posOfMPound)-1]==(len(deNovoPeptide[i])-2):
						padRight=2
					else:
						padRight=3
				else:
					if posOfMPound[0]==len(deNovoPeptide[i]):
						padRight=0
					elif  posOfMPound[0]==len(deNovoPeptide[i])-1:
						padRight=1
					elif  posOfMPound[0]==len(deNovoPeptide[i])-2:
						padRight=2
					else:
						padRight=3
			if (len(posOfMPound)>1):
				mPoundStr=deNovoPeptide[i][posOfMPound[0]-padLeft:posOfMPound[len(posOfMPound)-1]+padRight]
			if (len(posOfMPound)==1):
				mPoundStr=deNovoPeptide[i][posOfMPound[0]-padLeft:posOfMPound[0]+padRight]
			diffModsStr[aa]=mPoundStr
		cpstaticModsListOne=copy.copy(staticModsList[i])
		staticModsCounter={}
		staticModsRemove={}
		(context_Cleaned,contextWithAAOnly,mPoundCnt,mPoundModFinal,cpstaticModsListOne,staticModsCounter,adjustedMass,diffModMassChanged,staticModsRemove)=modsSubtitionPerScan(allMods[i],modContextList[i],contextList[i],mPoundStr,cpstaticModsListOne,NTERM,CTERM,dictDiffMods,dictDiffModsMass,dictStaticMods,dictStaticModsName,diffModsPos,diffModsStr,dictDiffModsWithOriginal,staticModsRemove,staticModsCounter)
		cpstaticModsListOne1=copy.copy(cpstaticModsListOne)	
		mPoundModFinalList.append(mPoundModFinal)
		newCTL=list(contextWithAAOnly[2:-2])
		for z in xrange(len(newCTL)):
			theoMass+=float(AAMass[newCTL[z]])
			for aa in dictStaticMods:
				if newCTL[z]==aa:
					theoMass+=float(dictStaticMods[aa])
					break
		subtractMass=0
		for aa in dictStaticMods:
			subtractMass=float(subtractMass)+(float(dictStaticMods[aa])*float(staticModsCounter[aa]))
		ifCTermPos=len(contextWithAAOnly[2:-2])-1
		#slin: issue 13 use newDic instead of the massArray
		(sqlVal,fileVal,staticModMassChanged)=splitMods(allMods[i],maxModLength,newDic,mPoundModFinal,cpstaticModsListOne1,NTERM,CTERM,dictStaticCNTermModsName,ifCTermPos)
		theoMass=theoMass+H2O+PROTON+staticModMassChanged+diffModMassChanged-subtractMass+adjustedMass
		#slin: issue 14
		if context_Cleaned.find("][")>-1:
			splitLine=fileVal.split("\t")
			peptideLength=0
			for j in xrange(len(context_Cleaned[2:-2])):
				if isAlpha(context_Cleaned[2:-2][j]):
					peptideLength=peptideLength+1
			findPos=[match.start() for match in re.finditer(re.escape("]["),context_Cleaned)]
			#now remove all nonalpha character before find pos to know the pos of such replacement
			for pos in findPos[::-1]:
				additionStrC=""
				newDisplay=Decimal(0)
				cur_Pos=len(re.sub("[^a-zA-Z]+","",context_Cleaned[2:pos]))-1
				#check if cur_pos the last postion of the peptide, if so, need to check if C-term
				firstHalf=len(context_Cleaned[:pos])-context_Cleaned[:pos][::-1].index("[")
				counter=0
				for j in xrange(len(context_Cleaned[pos:])-2):
					if isAlpha(context_Cleaned[pos:][j]):
						break
					else:
						counter=counter+1
				secondHalf=pos+counter-1
				for j in xrange(len(splitLine)/4):
					if str(splitLine[(j*4+1)])==str(cur_Pos):
						if (cur_Pos==(peptideLength-1)): #possible C-term
							if str(splitLine[(j*4+3)])=="C-term":
								additionStrC="]["+ str(round(Decimal(splitLine[((j+1)*4)]),2))
							else:
								newDisplay=newDisplay+Decimal(splitLine[((j+1)*4)])
						else:
							newDisplay=newDisplay+Decimal(splitLine[((j+1)*4)])
					if str(splitLine[(j*4+1)])=="":
						break;
				context_Cleaned=context_Cleaned[:firstHalf]+str(round(Decimal(newDisplay),2))+additionStrC+context_Cleaned[secondHalf:]
		newContext.append(context_Cleaned)
		newTheoMassList.append(theoMass)
		modsInfo.append(fileVal)
	return newContext,modsInfo,newTheoMassList

def cleanupMassArray(massArray):
	newMassArray={}
	for curkey,curlist in massArray.iteritems():
		sortedList=sorted(map(str_to_float_with_precision,list(set(curlist))))
		sortedStrList=map(str,sortedList)
		if (len(sortedList)>=2):
			for k in xrange(len(sortedList)-1):
				lenK=len(sortedStrList[k].split(".")[1])
				lenKplus1=len(sortedStrList[k+1].split(".")[1])
				if (lenK>2 and lenKplus1>2):
					continue
				if (abs(float(sortedList[k]-sortedList[k+1]))<0.01):
					if (lenK<=2):
						sortedList[k]=sortedList[k+1]
						sortedStrList[k]=sortedStrList[k+1]
					if (lenKplus1<=2):
						sortedList[k+1]=sortedList[k]
						sortedStrList[k+1]=sortedStrList[k]
		newSortedList=[]
		for i in xrange(len(sortedList)-1):
			if str(sortedList[i])[::-1].find('.')>=2:
				newSortedList.append(sortedList[i])
		newMassArray[curkey]=sorted(list(set(newSortedList)))
	return(newMassArray)
	
def mergeToDictionary(newDic,newMassArray):
	for curKey,curlist in newMassArray.iteritems():
		if curKey in newDic.keys():
			newList=curlist+newDic[curKey]
			newDic[curKey]=sorted(list(set(newList)))
		else:
			newDic[curKey]=sorted(list(set(curlist)))
	return(newDic)

#obtain online: https://stackoverflow.com/questions/18149400/convert-a-list-of-string-floats-to-a-list-of-floats-with-2-decimal-points
def str_to_float_with_precision(item):
	precision=7
	return float(Decimal(item, precision))

def modsArray(curMods,massArray):
	if (curMods.find(",",0)>-1):
		splitCurMods=curMods.split("), ((")
		for j in xrange(len(splitCurMods)):
			mods_=splitCurMods[j].replace("\"","").replace("(","").replace(")","").replace("[","").replace("]","").replace("\'","")
			eachMod= mods_.split(",")
			curName=eachMod[0].strip().lower()
			curMass=eachMod[1].strip()
			massArray[curName].append(curMass)

def replaceAtPos(id1,id2,originalstr,str):
	return originalstr[:id1]+str+originalstr[id2:]

def findChr(str,ch):
	return [i for i,ltr in enumerate(str) if ltr==ch]

def modVals(mergeCandModPound,lastPos,curPos,mass,name,pos,modAA,padding,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictStaticMods,dictStaticModsName):
	for curKey in sorted(mergeCandModPound):
		if (int(curKey)>=int(curPos)+padding):
			break
		else:
			#loop for 2 dictionaries, dictDiffMods/and staticMods, repalce needed
			for aa in dictDiffMods:
				if mergeCandModPound[curKey]==aa.upper():
					name=name+",\'"+dictDiffModsName[aa]+"\'"
					pos=pos+","+str(curKey)
					mass=mass+","+str(dictDiffModsMass[aa])
					modAA=modAA+","+aa.upper()
					break
			for aa in dictStaticMods:
				if mergeCandModPound[curKey]==aa.upper():
					name=name+",\'"+dictStaticModsName[aa]+"\'"
					pos=pos+","+str(curKey)
					mass=mass+","+str(dictStaticMods[aa])
					modAA=modAA+","+aa.upper()
					break
			del mergeCandModPound[curKey]
	return(mergeCandModPound,mass,name,pos,modAA)

def splitMods(curMods,maxModLength,newDic,mPoundMod,cpstaticModsListOne1,NTERM,CTERM,dictStaticCNTermModsName,ifCTermPos):
	mass=""
	name=""
	pos=""
	modAA=""
	sqlVal=""
	fileVal=""
	skipped=0
	massChanged=0
	lenCurMods=0
	mergeCandModPound={}
	#create dictionray for MPound and C
	for aa in cpstaticModsListOne1: #exit if all empty.
		if (len(cpstaticModsListOne1[aa]))>=0:
			for i in range(0,len(cpstaticModsListOne1[aa])):
				mergeCandModPound[cpstaticModsListOne1[aa][i]]=aa
	loopStart=0
	loopEnd=maxModLength
	if (curMods.find(",",0)==-1): #condition1: no mod exist
		if (NTERM!=0):
			sqlVal=sqlVal+",0,"+dictStaticCNTermModsName['N-term']+",N-term,"+str(NTERM)
			mass=","+str(NTERM)+mass
			name=","+dictStaticCNTermModsName['N-term']+name
			pos=",0"+pos
			modAA=",N-term"+modAA
			loopStart=1
			
		if (CTERM!=0):
			loopEnd=loopEnd-1
			sqlVal=sqlVal+","+str(ifCTermPos)+","+dictStaticCNTermModsName['C-term']+",C-term,"+str(CTERM)
			mass=mass+","+str(CTERM)
			name=name+","+dictStaticCNTermModsName['C-term']
			pos=pos+","+str(ifCTermPos)
			modAA=modAA+",C-term"
		if (len(mergeCandModPound)>0):
			(mergeCandModPound,mass,name,pos,modAA)=modVals(mergeCandModPound,0,10000,mass,name,pos,modAA,0,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictStaticMods,dictStaticModsName)
		(sqlVal,fileVal)=finalModsStrings(mass,name,pos,modAA,maxModLength)
	else:
		lastPos=0
		modAAupdate=0
		splitCurMods=curMods.split("), ((")
		lenCurMods=len(splitCurMods)
		modCTerm=0
		modNTerm=0
		if (lenCurMods>0):
			padding=0
			for j in xrange(len(splitCurMods)):
				massSkip=0
				getKey=splitCurMods[j].split(",")
				if (j==0):
					curKey=getKey[0].strip()[4:-1].lower()
				else:
					curKey=getKey[0].strip()[1:-1].lower()
				mods_=splitCurMods[j].replace("\"","").replace("(","").replace(")","").replace("[","").replace("]","").replace("\'","")
				eachMod=mods_.split(",")
				modAAupdate=0
				if eachMod[0].strip()!='Isobaric Substitution':
					if len(eachMod)==6:
						curPos=int(eachMod[5].strip())
					else:
						curPos=int(eachMod[4].strip())
					curMass=eachMod[1].strip()
					curName=eachMod[0].strip()
					curAA=eachMod[3].strip()
					(mergeCandModPound,mass,name,pos,modAA)=modVals(mergeCandModPound,lastPos,int(curPos),mass,name,pos,modAA,padding,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictStaticMods,dictStaticModsName)
					if (curAA=="C-term"):
						modCTerm=1
					if (curAA=="N-term"):
						modNTerm=1
					if (curName=='Insertion') or (curName=='Deletion'):
						massSkip=1
						replaceChar=eachMod[3].strip()
						if (curName=='Insertion'):
							padding=len(replaceChar)
						else:
							padding=len(replaceChar)*-1
						modAA=modAA+","
						modAAupdate=1
					else:
						pos=pos+","+str(curPos+padding)
						name=name+",\'"+curName+"\'"
					massUsed=round(Decimal(curMass),6)
					if (curName.find("->",0)>=0):
						splitLine=curName.split("->")
						newChar=AA[splitLine[1].capitalize()]
						if (newChar!="?"): #matching, no diff Mass
							massSkip=1
							mass=mass+","
							modAA=modAA+","
							modAAupdate=1
					if massSkip==0:
						for k in xrange(len(newDic[curName.lower()])):
							if abs(float(newDic[curName.lower()][k])-float(curMass))<0.1:
								massUsed=round(Decimal(newDic[curName.lower()][k]),6)
							if len(eachMod)==6:
								toBeReplaced=eachMod[4].strip()
							else:
								toBeReplaced=eachMod[3].strip()
							if (toBeReplaced=='N-term'):
								massUsed=Decimal(massUsed,6)+Decimal(NTERM,6)
							elif (toBeReplaced=='C-term'):
								massUsed=Decimal(massUsed,6)+Decimal(CTERM,6)
							break
						mass=mass+","+str(massUsed)
						massChanged+=float(massUsed)
					lastPos=int(curPos)
					if (modAAupdate==0):
						modAA=modAA+","+curAA
			(mergeCandModPound,mass,name,pos,modAA)=modVals(mergeCandModPound,lastPos,10000,mass,name,pos,modAA,padding,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictStaticMods,dictStaticModsName)
			if (NTERM!=0) and (modNTerm==0):
				mass=","+str(NTERM)+mass
				name=","+dictStaticCNTermModsName['N-term']+name
				pos=",0"+pos
				modAA=",N-term"+modAA
			if (CTERM!=0) and (modCTerm==0):
				mass=mass+","+str(CTERM)
				name=name+","+dictStaticCNTermModsName['C-term']
				pos=pos+","+str(ifCTermPos)
				modAA=modAA+",C-term"
		(sqlVal,fileVal)=finalModsStrings(mass,name,pos,modAA,maxModLength)
	return (sqlVal,fileVal,massChanged)

def finalModsStrings(mass,name,pos,modAA,maxModLength):
	sqlVal=""
	fileVal=""
	shiftPos=0
	nameL=name.split(",")
	posL=pos.split(",")
	massL=mass.split(",")
	modAAL=modAA.split(",")
	for i in range(1,len(nameL)):
		sqlVal=sqlVal+","+posL[i]+","+nameL[i]+","+modAAL[i]+","+massL[i]
	for i in range(name.count(','),maxModLength):
		sqlVal=sqlVal+",None,None,None,None"
		fileVal=sqlVal.replace(",","\t").replace("None","").replace("\"","").replace("\'","")
	return (sqlVal,fileVal)

def CalculatePPMInfo(obsMass,newTheoMass,oneSubtractlg10EM):
	tenPowerNegOneSubtractlg10EMList=[]
	newPPMList=[]
	for i in xrange(len(obsMass)):
		newPPM=(float(obsMass[i])-float(newTheoMass[i]))/float(newTheoMass[i])*1000000
		newPPMList.append(newPPM)
		tenPowerNegOneSubtractlg10EMList.append(math.pow(10,(-1*float(oneSubtractlg10EM[i]))))
	return newPPMList,tenPowerNegOneSubtractlg10EMList

def finalPeptide(peptide):
	return re.sub(r'\d+','',peptide).replace('[','').replace(']','').replace('.','')

def insertBaseDataForPPMCalculation(conn,idList,chargeList,PPMList,newPPMList,EMProbabilityList,oneSubtractlg10EMList,finalPeptideList,tenPowerNegOneSubtractlg10EMList):
	dataList=[]
	for i in xrange(len(chargeList)):
		dataList.append((idList[i],chargeList[i],PPMList[i],newPPMList[i],round(newPPMList[i],1),EMProbabilityList[i],oneSubtractlg10EMList[i],tenPowerNegOneSubtractlg10EMList[i],finalPeptideList[i][2:-2]))
	conn.executemany("INSERT INTO baseDataForPPMCalculation(id,charge,PPM,newPPM,roundedPPM,EMProbability,oneSubtractlg10EM,tenPowerNegOneSubtractlg10EM,finalPeptide) VALUES (?,?,?,?,?,?,?,?,?)",dataList)
	conn.commit()

def createPPMRangesFDRTables(conn):
	conn.execute("DROP TABLE if exists ppmRanges")
	conn.execute("DROP TABLE if exists finalPeptidePPMInfo")
	conn.commit()
	conn.execute("CREATE TABLE ppmRanges(id int,fractionID int,lowerBound decimal(5,2),upperBound decimal(5,2),mean decimal(14,10),stdev decimal(14,10))")
	conn.execute("CREATE TABLE finalPeptidePPMInfo(id int,fractionID int,scan int,newTheoMass decimal(14,10),newPPM decimal(18,10),nTerm varchar(1),finalPeptide TEXT,cTerm varchar(1));")
	conn.commit()
	mods=""
	for i in range(0,maxModLength):
		mods=mods+",modPos"+str(i)+" int"
		mods=mods+",modName"+str(i)+" varchar(255)"
		mods=mods+",modAA"+str(i)+" varchar(50)"
		mods=mods+",modMass"+str(i)+" float"
	conn.commit()
	conn.execute("DROP TABLE if exists modsPerScan")
	conn.execute("DROP TABLE if exists cleanDataFDR")
	conn.execute("DROP TABLE if exists cleanDataFDRUnique")
	conn.execute("DROP TABLE if exists FDREMMapping")
	conn.commit()
	conn.execute("CREATE TABLE modsPerScan(id int,fractionID int,scan int" + mods+")")
	conn.execute("CREATE TABLE cleanDataFDR(id int,charge int,fractionID int,PPM decimal(14,8),roundedPPM decimal(10,1),EMProbability decimal(14,10),oneSubtractlg10EM decimal(14,10),newEM decimal(8,2),finalPeptide varchar(512),FDR decimal(8,4))")
	conn.execute("CREATE TABLE cleanDataFDRUnique(id int,fractionID int,charge int,PPM decimal(14,8),roundedPPM decimal(10,1),EMProbability decimal(14,10),oneSubtractlg10EM decimal(14,10),newEM decimal(8,2),finalPeptide varchar(512),FDR decimal(8,4))")
	conn.execute("CREATE TABLE FDREMMapping(newEM decimal(8,2),FDR decimal(14,10))")
	conn.commit()

def createPPMTables(conn):
	conn.execute("DROP TABLE if exists baseDataForPPMCalculation")
	conn.execute("DROP TABLE if exists cleanDataPPM")
	conn.execute("DROP TABLE if exists ppmRangesPreFraction")
	conn.commit()
	conn.execute("CREATE TABLE baseDataForPPMCalculation(id int,charge int,PPM decimal(14,8),newPPM decimal(14,8),roundedPPM decimal(10,1),EMProbability decimal(14,10),OneSubtractlg10EM decimal(14,10),tenPowerNegOneSubtractlg10EM decimal(14,10),finalPeptide varchar(512))")
	conn.execute("CREATE TABLE cleanDataPPM (id int,charge int,PPM decimal(14,8),roundedPPM decimal(10,1),EMProbability decimal(14,10),oneSubtractlg10EM decimal(14,10),tenPowerNegOneSubtractlg10EM decimal(14,10),finalPeptide varchar(512))")
	conn.execute("CREATE TABLE ppmRangesPreFraction(id int,lowerBound decimal(5,2),upperBound decimal(5,2),mean decimal(14,10),stdev decimal(14,10))")
	conn.commit()

def insertCleanDataPPM(conn): #per fraction use.
	conn.execute("insert into cleanDataPPM(id,charge,PPM,roundedPPM,EMProbability,oneSubtractlg10EM,tenPowerNegOneSubtractlg10EM,finalPeptide) select a.id,a.charge,a.ppm,a.roundedPPM,a.EMProbability,a.oneSubtractlg10EM,tenPowerNegOneSubtractlg10EM,a.FinalPeptide from baseDataForPPMCalculation a,(select max(oneSubtractlg10EM) oneSubtractlg10EM,charge,finalPeptide from baseDataForPPMCalculation group by finalPeptide,charge) b where ((a.PPM<=50 and a.PPM>=-50) or (a.newPPM<=50 and a.newPPM>=-50)) and (a.oneSubtractlg10EM=b.oneSubtractlg10EM) and (a.charge=b.charge) and (a.finalPeptide=b.finalPeptide) order by a.oneSubtractlg10EM desc")
	conn.commit()

def populatePerfectFDRTableDummy(conn):
	FDRDefault=[]
	i=(-1*PPMRANGE)-5
	while (i<=(PPMRANGE+5)):
		FDRDefault.append((round(i,2),-1000))
		i=i+0.01
	conn.executemany("insert into FDREMMapping(newEM,FDR) values(?,?)",FDRDefault)
	conn.commit()
	conn.execute("DROP INDEX IF EXISTS FDREMMapping_EM")
	conn.commit()
	conn.execute("CREATE INDEX FDREMMapping_EM ON FDREMMapping(newEM)")
	conn.commit()

def findPPMRange(conn,fractionID):
	cursor=conn.cursor()
	cursor1=conn.cursor()
	ppmRangeUpperBound={}
	ppmRangeLowerBound={}
	ppmRangesList=[]
	ppmRangeList=[]
	query="select roundedPPM,count(*) as total,(count(*)-round(sum(TenPowerNegOneSubtractlg10EM))) as TP,round(sum(TenPowerNegOneSubtractlg10EM)) as FP,roundedPPM*(count(*)-sum(TenPowerNegOneSubtractlg10EM)) as weightedSum from cleanDataPPM group by roundedPPM order by roundedPPM"
	cursor1.execute(query)
	rows=cursor1.fetchall()
	aa_histogramTP=defaultdict(lambda:0)
	aa_histogramFP=defaultdict(lambda:0)
	aa_histogramWeightedSum=defaultdict(lambda:0)
	for row in rows:
		PPM=row[0]
		tp=row[2]
		fp=row[3]
		weightedSum=row[4]
		aa_histogramFP[PPM]=fp
		aa_histogramTP[PPM]=tp
		aa_histogramWeightedSum[PPM]=weightedSum
	weighted_sum=0 # ppm * (total - fp). This is more accurate than ppm*(max(0,total-fp)) in larger ppm ranges
	n_tp=0 #number of tp identifications
	for key,value in aa_histogramTP.iteritems():
		weighted_sum+=aa_histogramWeightedSum[key]
		n_tp+=aa_histogramTP[key]
	if (n_tp==0): # find starting mean
		prior_mean=0.0
	else:
		prior_mean=round(Decimal((weighted_sum/n_tp),PRECISION),PRECISION)
	numerator=0
	for key,value in aa_histogramTP.iteritems():
		decimalKey=Decimal((key),2)
		numerator+=Decimal(aa_histogramTP[key],PRECISION)*(decimalKey-Decimal(prior_mean,PRECISION))*(decimalKey-Decimal(prior_mean,PRECISION))
	if (n_tp<=1):
		prior_stdev=0
	else:
		prior_stdev=round(Decimal(math.sqrt(abs(Decimal(numerator,PRECISION)/(Decimal(n_tp,PRECISION)-1))),PRECISION),PRECISION)
	upperbound=round(Decimal((Decimal(prior_mean,PRECISION)+(Decimal(prior_stdev,PRECISION)*3)),PRECISION),1)
	lowerbound=round(Decimal((Decimal(prior_mean,PRECISION)-(Decimal(prior_stdev,PRECISION)*3)),PRECISION),1)
	upperbound=upperbound if upperbound<PPMRANGE else PPMRANGE
	lowerbound=lowerbound if lowerbound>(-1*PPMRANGE) else (-1*PPMRANGE)
	# initialize prior mean, stdev
	prior_prior_mean=1e6
	prior_prior_stdev=1e6
	flg_zero_stdev=0
	prior_tp=n_tp
	flg_bad_tp=0
	if printMsg:
		print "after first iteration: prior_stdev: %s upperbound: %s lowerbound: %s " %(prior_stdev,upperbound,lowerbound)
	if (prior_stdev>0) and (prior_stdev>0):
		if (upperbound<=PPMRANGE or lowerbound>=(-1*PPMRANGE)): 
			iteration=0
			while ((iteration<PPMITERATIONS) and (abs(Decimal(prior_mean,PRECISION)-Decimal(prior_prior_mean,PRECISION))/abs(Decimal(prior_prior_mean,PRECISION)) > 0.05) or (abs(Decimal(prior_stdev,PRECISION)-Decimal(prior_prior_stdev,PRECISION))/abs(Decimal(prior_prior_stdev,PRECISION))>0.05) and (flg_zero_stdev==0) and (flg_bad_tp==0)):
				iteration+=1
				n_tp=0
				weighted_sum=0
				for key, value in sorted(aa_histogramTP.iteritems()):
					if (key>=lowerbound) and (key<=upperbound):
						weighted_sum+=aa_histogramWeightedSum[key]
						n_tp+=aa_histogramTP[key]
				if (n_tp==0):
					new_mean=0
				else:
					new_mean=round(Decimal(weighted_sum,PRECISION)/Decimal(n_tp,PRECISION),PRECISION)
				numerator=0
				for key, value in sorted(aa_histogramTP.iteritems()):
					if (key>=lowerbound) and (key<=upperbound):
						decimalKey=Decimal((key),2)
						numerator+=Decimal(aa_histogramTP[key],PRECISION)*(decimalKey-Decimal(new_mean,PRECISION))*(decimalKey-Decimal(new_mean,PRECISION))
				if ((float(n_tp)-1)==0):
					new_stdev=0
				else:
					new_stdev=math.sqrt(abs(Decimal(numerator,PRECISION))/(Decimal(n_tp,PRECISION)-1))
				if (new_stdev==0): # bounds are the same don't use this iteration
					new_stdev=prior_prior_stdev
					new_mean=prior_prior_mean
					flg_zero_stdev=1
				elif ((prior_tp!=0) and (Decimal(n_tp,PRECISION)/Decimal(prior_tp,PRECISION)<0.75)): # new range, while more uniform throws out too many right answers
					new_stdev=prior_prior_stdev
					new_mean=prior_prior_mean
					flg_bad_tp=1
				upperbound=round(Decimal((Decimal(new_mean,PRECISION)+(Decimal(new_stdev,PRECISION)*3)),PRECISION),1)
				lowerbound=round(Decimal((Decimal(new_mean,PRECISION)-(Decimal(new_stdev,PRECISION)*3)),PRECISION),1)
				upperbound=upperbound if (upperbound<PPMRANGE) else PPMRANGE
				lowerbound=lowerbound if (lowerbound>(-1*PPMRANGE)) else (-1*PPMRANGE)
				if printMsg:
					print "iteration: %s prior_stdev: %s upperbound: %s lowerbound: %s " %(iteration+1,prior_stdev,upperbound,lowerbound)
				prior_tp=n_tp
				prior_prior_mean=prior_mean
				prior_prior_stdev=prior_stdev
				prior_mean=new_mean
				prior_stdev=new_stdev
		else:
			new_mean=prior_mean
			lowerbound=-1*PPMRANGE
			upperbound=PPMRANGE
			new_stdev=prior_stdev
	else:
		new_mean=prior_mean
		lowerbound=lowerbound
		upperbound=upperbound
		new_stdev=prior_stdev
	ppmRangeUpperBound=upperbound
	ppmRangeLowerBound=lowerbound
	ppmRangesList.append((fractionID,lowerbound,upperbound,new_mean,new_stdev))
	ppmRangeList.append((lowerbound,upperbound,new_mean,new_stdev)) 
	cursor.executemany("INSERT INTO ppmRanges(fractionID,lowerBound,upperBound,mean,stdev) VALUES (?,?,?,?,?)",ppmRangesList)
	cursor.executemany("INSERT INTO ppmRangesPreFraction(lowerBound,upperBound,mean,stdev) VALUES (?,?,?,?)",ppmRangeList)
	conn.commit()
	return ppmRangeUpperBound,ppmRangeLowerBound

def insertCleanDataFDR(conn,fractionID,ppmRangeUpperBound,ppmRangeLowerBound):  #per fraction now.
	query="insert into cleanDataFDR(id,fractionID,charge,PPM,roundedPPM,EMProbability,oneSubtractlg10EM,finalPeptide) select id,"+str(fractionID)+",charge,ppm,roundedPPM,EMProbability,oneSubtractlg10EM,FinalPeptide from baseDataForPPMCalculation"
	conn.execute(query)
	conn.commit()
	conn.execute("DROP INDEX IF EXISTS cleanDataFDR_id")
	conn.commit()
	conn.execute("CREATE INDEX cleanDataFDR_id ON cleanDataFDR(id)")
	conn.commit()
	newEMs=[]
	cursor=conn.cursor()
	dataSetQuery="select id,EMProbability,oneSubtractlg10EM from cleanDataFDR where fractionID="+str(fractionID)
	cursor.execute(dataSetQuery)
	dataSetRows=cursor.fetchall()
	for dataSetRow in dataSetRows:  #sqlite doesnt have built in log function, need to do it seperately.
		id=dataSetRow[0]
		EMProbability=dataSetRow[1]
		oneSubtractlg10EM=dataSetRow[2]
		if (EMProbability>0.5):
			newEM=oneSubtractlg10EM
		elif (EMProbability==0):
			newEM=0.0000000001
		else:
			newEM=math.log(EMProbability,10)
		newEMs.append((round(newEM,2),id))
	conn.executemany("UPDATE cleanDataFDR SET NewEM=? WHERE id=?",newEMs)
	conn.commit()
	conn.execute("insert into cleanDataFDRUnique(id,fractionID,charge,PPM,roundedPPM,EMProbability,oneSubtractlg10EM,NewEM,finalPeptide) select a.id, a.fractionID,a.charge,a.ppm,a.roundedPPM,a.EMProbability,a.oneSubtractlg10EM,a.NewEM,a.FinalPeptide from cleanDataFDR a,(select max(oneSubtractlg10EM) oneSubtractlg10EM,charge,finalPeptide from cleanDataFDR where fractionID="+str(fractionID)+" group by charge,finalPeptide) b where a.fractionID="+str(fractionID)+" and a.oneSubtractlg10EM=b.oneSubtractlg10EM and a.charge=b.charge and a.finalPeptide=b.finalPeptide")
	conn.commit()
	conn.execute("DROP INDEX IF EXISTS cleanDataFDRUnique_id")
	conn.commit()
	conn.execute("CREATE INDEX cleanDataFDRUnique_id ON cleanDataFDRUnique(id)")
	conn.commit()

def calFalseDiscoveryRate(conn):
	MAXCNT=10000  #we want to update table in batch to avoid the list oversize for big search
	total=0
	fp=0
	minFDR=1000
	maxFDR=-1000
	cursor=conn.cursor()
	query="select id,(1-EMProbability) pOfIncorrect,NewEM from cleanDataFDRUnique order by oneSubtractlg10EM desc"
	cursor.execute(query)
	rows=cursor.fetchall()
	localCounter=0
	FDRInfo=[]
	FDRUpdate=[]
	for row in rows:
		id=row[0]
		total+=1
		fp+=row[1]
		NewEM=row[2]
		if (minFDR>NewEM):
			minFDR=NewEM-0.5
		if (maxFDR<NewEM):
			maxFDR=NewEM+0.5
		FDR=Decimal(fp,8)/Decimal(total,8)
		FDRInfo.append((round(FDR,4),round(NewEM,2)))
		FDRUpdate.append((round(FDR,4),id))
		if (localCounter==MAXCNT):
			conn.executemany("UPDATE cleanDataFDRUnique SET FDR=? WHERE id=?",FDRUpdate)
			conn.executemany("UPDATE cleanDataFDR SET FDR=? WHERE id=?",FDRUpdate)
			conn.executemany("UPDATE FDREMMapping SET FDR=? WHERE newEM=?",FDRInfo)
			conn.commit()
			FDRInfo=[]
			FDRUpdate=[]
			localCounter=0
		localCounter=localCounter+1
	conn.executemany("UPDATE cleanDataFDRUnique SET FDR=? WHERE id=?",FDRUpdate)
	conn.executemany("UPDATE cleanDataFDR SET FDR=? WHERE id=?",FDRUpdate)
	conn.executemany("UPDATE FDREMMapping SET FDR=? WHERE newEM=?",FDRInfo)
	conn.commit()
	populatePerfectFDRTable(conn,FDRInfo)
	updateFDRonDuplicateRows(conn)

def populatePerfectFDRTable(conn,FDRInfo):
	cursor=conn.cursor()
	conn.commit()
	dataSetQuery="SELECT newEM,FDR FROM FDREMMapping WHERE FDR<>-1000 ORDER BY newEM"
	cursor.execute(dataSetQuery)
	dataSetRows=cursor.fetchall()
	newEMList=[]
	FDRList=[]
	verifyFDR=[]
	newEMListValue=[]
	for dataSetRow in dataSetRows:
		newEMList.append(dataSetRow[0])
		FDRList.append(dataSetRow[1])
	for i in range(0,len(newEMList)-1):
		strEM=newEMList[i]
		endEM=newEMList[i+1]
 		strEMFDR=FDRList[i]
		endEMFDR=FDRList[i+1]
		gap=(endEMFDR-strEMFDR)/((endEM-strEM)/0.01)
		for j in range(1,abs(int((endEM-strEM)/0.01))+1):
			newEMListValue.append((round(strEMFDR+(j*gap),4),round(strEM+(j*0.01),2)))
	conn.executemany("UPDATE FDREMMapping SET FDR=? WHERE newEM=?",newEMListValue)
	conn.commit()

def updateFDRonDuplicateRows(conn):
	updateSQL="UPDATE cleanDataFDR SET FDR=(SELECT e.FDR FROM FDREMMapping e WHERE cleanDataFDR.newEM= e.newEM) WHERE cleanDataFDR.FDR IS NULL"
	conn.execute(updateSQL)
	conn.commit()

def getScansPerFraction(conn,fractionID,dictDiffMods,dictStaticMods):
	cursor=conn.cursor()
	query="SELECT id,scan,context,mods,mod_context,de_novo_peptide,obs_mh,em_probability,log_em_probability,charge,ppm FROM result WHERE top_result=1 and fraction_id="+str(fractionID)+" ORDER BY id"
	cursor.execute(query)
	scans=cursor.fetchall()
	idList=[]
	scanList=[]
	contextList=[]
	modsList=[]
	modContextList=[]
	deNovoPeptideList=[]
	obsMassList=[]
	EMProbabilityList=[]
	oneSubtractlg10EMList=[]
	chargeList=[]
	staticModsList=[]
	PPMList=[]
	for scan in scans:
		staticModsDict={}
		idList.append(scan[0])
		scanList.append(scan[1])
		contextList.append(scan[2])
		modsList.append(scan[3])
		modContextList.append(scan[4])
		deNovoPeptideList.append(scan[5])
		obsMassList.append(scan[6])
		EMProbabilityList.append(scan[7])
		oneSubtractlg10EMList.append(scan[8])
		chargeList.append(scan[9])
		PPMList.append(scan[10])
		if (len(dictStaticMods))>0:
			for localAA in dictStaticMods:
				oneStaticModList=findChr(scan[2][2:-2],localAA)
				staticModsDict[localAA]=oneStaticModList
		staticModsList.append(staticModsDict)
	return(idList,scanList,contextList,modContextList,modsList,deNovoPeptideList,obsMassList,EMProbabilityList,oneSubtractlg10EMList,chargeList,staticModsList,PPMList)

def insertProcessedInfo(fractionID,idList,scanList,newPPMList,newTheoMassList,finalPeptideList,modsInfoList):
	finalPeptidePPMInfoList=[]
	modsList=[]
	for i in range(0,len(idList)):
		curModList=[]
		finalPeptidePPMInfoList.append((idList[i],fractionID,scanList[i],newTheoMassList[i],newPPMList[i],finalPeptideList[i][:1],finalPeptideList[i][2:-2],finalPeptideList[i][-1:]))
		curModList.append(idList[i])
		curModList.append(fractionID)
		curModList.append(scanList[i])
		splitModInfos=modsInfoList[i].split("\t")
		for j in xrange(1,(maxModLength*4)+1):
			if(len(splitModInfos)>1):
				curModList.append(splitModInfos[j])
			else:
				curModList.append('')
		modsList.append(curModList)
	conn.executemany("INSERT INTO finalPeptidePPMInfo(id,fractionID,scan,newTheoMass,newPPM,nTerm,finalPeptide,cTerm) VALUES(?,?,?,?,?,?,?,?)",finalPeptidePPMInfoList)
	conn.commit()
	query="insert into modsPerScan(id,fractionID,scan"
	for i in range(0,maxModLength):
		query=query+",modPos"+str(i)
		query=query+",modName"+str(i)
		query=query+",modAA"+str(i)
		query=query+",modMass"+str(i)
	query=query+") values (?"
	for i in range(0,maxModLength*4+2):
		query=query+",?"
	query=query+")"
	conn.executemany(query,modsList)
	conn.commit()
	
def isAlpha(x):
	return x in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

def writeHeader(fout):
	newHeader="ScanF\tCharge\tRetention Time\tObs M+H\tTheo M+H\tPPM\tPPM Upper Bound\tPPM Lower Bound\tPPM In Range\tFDR\tEM Probability\t1-lg10 EM\tSpectrum Score\tAlignment Score\tComposite Score\tUnique Siblings\tContext Mod Variants\tNum Mod Occurrences\tContext\tMod Context\tMods\tMod Ambig Edges\tMod Ranges\tProteins\tDe Novo Peptide\tDe Novo Score\tMatching Tag Length\tNum Matches\tN-terminus\tFinal Peptide\tC-terminus"
	fout.write("%s"%(newHeader))
	for i in range(0,maxModLength):
		fout.write("\tmodPos"+str(i))
		fout.write("\tmodName"+str(i))
		fout.write("\tmodAA"+str(i))
		fout.write("\tmodMass"+str(i))
	fout.write("\n")

def processProteinForPrintout(protein,DisplayProteinNum):
	splitVia="', '"
	newProtein=""
	splitProtein=protein[1:-1].split(splitVia)
	cntProtein=len(splitProtein)
	displayCnt=int(DisplayProteinNum)
	proteinCnt=displayCnt
	appendSquarBracket="]"
	if cntProtein<=displayCnt:
		proteinCnt=cntProtein
	if cntProtein>0:
		newProtein=str(cntProtein)+":: "
		for i in range(0,int(proteinCnt)):
			if i==0:
				newProteinsub="["+splitProtein[i][:MaxProtenLengthEach]
			else:
				newProteinsub=newProteinsub+", "+splitProtein[i][:MaxProtenLengthEach]
		newProtein=newProtein+newProteinsub+appendSquarBracket
	return newProtein

def createIndexes(conn):
	cursor=conn.cursor()
	cursor.execute("create index IF NOT EXISTS result_1 on result(id)")
	cursor.execute("create index IF NOT EXISTS result_2 on result(top_result,fraction_id)")
	cursor.execute("create index IF NOT EXISTS finalPeptidePPMInfo_1 on finalPeptidePPMInfo(id)")
	cursor.execute("create index IF NOT EXISTS cleanDataFDR_1 on cleanDataFDR(id,FDR)")
	cursor.execute("create index IF NOT EXISTS modsPerScan_1 on modsPerScan(id)")
	conn.commit()

def writeScans(conn,fout,outFile,fractionMapping,fractionID,outputPerFraction,FDRCutoff,logEMCutoff,DisplayProteinNum,writeFileCnt):
	writeFileCnt=writeFileCnt+1
	if (fractionID<10):
		curFractionID="0"+str(fractionID)
	else:
		curFractionID=str(fractionID)
	if (writeFileCnt==1):
		outFileFull=outFile+".txt"
		fout=open(outFileFull,'w')
		writeHeader(fout)
	if (outputPerFraction=="yes"):
		outFileNamePreFraction=outFile+"_"+fractionMapping[curFractionID]+".txt"
		foutPreFraction=open(outFileNamePreFraction,'w')
		writeHeader(foutPreFraction)
	cursorPPM=conn.cursor()
	query="select lowerBound,upperBound from ppmRanges where fractionID="+str(fractionID)
	cursorPPM.execute(query)
	ppm=cursorPPM.fetchone()
	ppmRangeLowerBound=ppm[0]
	ppmRangeUpperBound=ppm[1]
	cursor=conn.cursor()
	#slin# check SQL for error
	query="select r.*,c.FDR,c.newEM,f.newTheoMass,f.newPPm,f.nTerm,f.finalPeptide,f.cTerm, m.* from result r left join finalPeptidePPMInfo f on r.id=f.id left join cleanDataFDR c on r.id=c.id left join modsPerScan m on r.id=m.id where r.top_result=1 and r.fraction_id="+str(fractionID) #issue 12
	cursor.execute(query)
	scans=cursor.fetchall()
	#4:theo_mh,6:ppm,11:log_em_probability,32:FDR,34:newTheomMass,35:newPPM
	for scan in scans:
		floatFDRCutoff=float(FDRCutoff)
		floatlogEMCutoff=float(logEMCutoff)
		if scan[32]==None:
			checkFDR=100000
		else:
			checkFDR=scan[32]
		if (float(checkFDR)<=floatFDRCutoff or float(scan[11])>=floatlogEMCutoff):
			ppmInRange='No'
			protein=processProteinForPrintout(scan[18],DisplayProteinNum)
			theoMass=scan[34]
			ppm=scan[35]
			if (ppmRangeLowerBound<=scan[35] and scan[35]<=ppmRangeUpperBound):
				ppmInRange='Yes'
				theoMass=scan[34]
				ppm=scan[35]
			elif (ppmRangeLowerBound<=scan[6] and scan[6]<=ppmRangeUpperBound):
				ppmInRange='Yes'
				theoMass=scan[4]
				ppm=scan[6]
			if scan[32]==None:
				curFDR=''
			else:
				curFDR=scan[32]
			if scan[33]==None:
				curNewEM=''
			else:
				curNewEM=scan[33]
			#theomass:34, ppm:35 then on
			fout.write("%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(fractionMapping[curFractionID],scan[1],scan[2],scan[5],scan[3],theoMass,ppm,ppmRangeUpperBound,ppmRangeLowerBound,ppmInRange,curFDR,scan[10],scan[11],scan[8],scan[7],scan[9],scan[26],scan[27],scan[28],scan[12],scan[13],scan[15],scan[16],scan[17],protein,scan[21],scan[23],scan[19],scan[24],scan[36],scan[37],scan[38],scan[42],scan[43],scan[44],scan[45],scan[46],scan[47],scan[48],scan[49],scan[50],scan[51],scan[52],scan[53],scan[54],scan[55],scan[56],scan[57],scan[58],scan[59],scan[60],scan[61],scan[62],scan[63],scan[64],scan[65],scan[66],scan[67],scan[68],scan[69],scan[70],scan[71],scan[72],scan[73],scan[74],scan[75],scan[76],scan[77],scan[78],scan[79],scan[80],scan[81],scan[82],scan[83],scan[84],scan[85],scan[86],scan[87],scan[88],scan[89],scan[90],scan[91],scan[92],scan[93],scan[94],scan[95],scan[96],scan[97],scan[98],scan[99],scan[100],scan[101],scan[102],scan[103],scan[104],scan[105],scan[106],scan[107],scan[108],scan[109],scan[110],scan[111],scan[112],scan[113],scan[114],scan[115],scan[116],scan[117],scan[118],scan[119],scan[120],scan[121],scan[122],scan[123],scan[124],scan[125],scan[126],scan[127],scan[128],scan[129],scan[130],scan[131],scan[132],scan[133],scan[134],scan[135],scan[136],scan[137],scan[138],scan[139],scan[140],scan[141],scan[142],scan[143],scan[144],scan[145],scan[146],scan[147],scan[148],scan[149],scan[150],scan[151],scan[152],scan[153],scan[154],scan[155],scan[156],scan[157],scan[158],scan[159],scan[160],scan[161]))
			if (outputPerFraction=="yes"):
				foutPreFraction.write("%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(fractionMapping[curFractionID],scan[1],scan[2],scan[5],scan[3],theoMass,ppm,ppmRangeUpperBound,ppmRangeLowerBound,ppmInRange,curFDR,scan[10],scan[11],scan[8],scan[7],scan[9],scan[26],scan[27],scan[28],scan[12],scan[13],scan[15],scan[16],scan[17],protein,scan[21],scan[23],scan[19],scan[24],scan[36],scan[37],scan[38],scan[42],scan[43],scan[44],scan[45],scan[46],scan[47],scan[48],scan[49],scan[50],scan[51],scan[52],scan[53],scan[54],scan[55],scan[56],scan[57],scan[58],scan[59],scan[60],scan[61],scan[62],scan[63],scan[64],scan[65],scan[66],scan[67],scan[68],scan[69],scan[70],scan[71],scan[72],scan[73],scan[74],scan[75],scan[76],scan[77],scan[78],scan[79],scan[80],scan[81],scan[82],scan[83],scan[84],scan[85],scan[86],scan[87],scan[88],scan[89],scan[90],scan[91],scan[92],scan[93],scan[94],scan[95],scan[96],scan[97],scan[98],scan[99],scan[100],scan[101],scan[102],scan[103],scan[104],scan[105],scan[106],scan[107],scan[108],scan[109],scan[110],scan[111],scan[112],scan[113],scan[114],scan[115],scan[116],scan[117],scan[118],scan[119],scan[120],scan[121],scan[122],scan[123],scan[124],scan[125],scan[126],scan[127],scan[128],scan[129],scan[130],scan[131],scan[132],scan[133],scan[134],scan[135],scan[136],scan[137],scan[138],scan[139],scan[140],scan[141],scan[142],scan[143],scan[144],scan[145],scan[146],scan[147],scan[148],scan[149],scan[150],scan[151],scan[152],scan[153],scan[154],scan[155],scan[156],scan[157],scan[158],scan[159],scan[160],scan[161]))
	if (outputPerFraction=="yes"):
		foutPreFraction.close()
	return (writeFileCnt,fout)

def loadFractionInfo(fracMappingFile):
	fileFractionMapping=pickle.load(fracMappingFile)
	fractionMapping={}
	fractionList=[]
	for i in range(0,len(fileFractionMapping)):
		fractionMapping[fileFractionMapping[i][0]]=fileFractionMapping[i][1]
		fractionList.append(int(fileFractionMapping[i][0]))
	return (fractionMapping,fractionList)

def loadInitFile(initFile):
	paramsDict=DataFile.parseParams_v1(initFile)
	NTERM=0
	CTERM=0	
	dictStaticMods={}
	dictStaticModsName={}
	dictStaticCNTermModsName={}
	if (len(paramsDict['Static Mods'])>0): #the value to add to exist mods val
		for modList in paramsDict['Static Mods']:
			if modList[1]=="C-term":
				CTERM=float(paramsDict['Static Mods'][modList])
				dictStaticCNTermModsName[modList[1]]=modList[0]
			elif modList[1]=="N-term":
				NTERM=float(paramsDict['Static Mods'][modList])
				dictStaticCNTermModsName[modList[1]]=modList[0]
			else:
				dictStaticMods[modList[1]]=paramsDict['Static Mods'][modList]
				dictStaticModsName[modList[1]]=modList[0]
	dictDiffMods={}
	dictDiffModsMass={}
	dictDiffModsName={}
	dictDiffModsWithOriginal={}
	if (len(paramsDict['Diff Mods'])>0):  
		for aa in paramsDict['Diff Mods']:
			dictDiffMods[paramsDict['Diff Mods'][aa][1]]=aa
			dictDiffModsMass[paramsDict['Diff Mods'][aa][1]]=paramsDict['Diff Mods'][aa][2]
			dictDiffModsName[paramsDict['Diff Mods'][aa][1]]=paramsDict['Diff Mods'][aa][0]
			dictDiffModsWithOriginal[paramsDict['Diff Mods'][aa][1]]=paramsDict['Diff Mods'][aa][3]
	if (len(paramsDict['Amino Acids'])>0):  #replace the AA array values with the one found in .ini file
		for aa in paramsDict['Amino Acids']:
			AAMass[aa]=paramsDict['Amino Acids'][aa][2]
	return (dictStaticMods,dictStaticModsName,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictDiffModsWithOriginal,NTERM,CTERM,dictStaticCNTermModsName)

def modNCTerm(newStr,ncTermAdjustMass,ncTerm,adjustedMass):
	if (ncTerm=="N-term"):
		adjustedMass+=ncTermAdjustMass
		newStr=newStr[:2]+"["+str(round(Decimal(ncTermAdjustMass),2))+"]"+newStr[2:]
	if (ncTerm=="C-term"):
		adjustedMass+=ncTermAdjustMass
		newStr=newStr[:-2]+"["+str(round(Decimal(ncTermAdjustMass),2))+"]"+newStr[-2:]
	return (newStr,adjustedMass)

def loadUnimodDict(UnimodDictFile):
	UnimodDict=pickle.load(UnimodDictFile)
	return UnimodDict

if __name__ == '__main__':
	if (len(sys.argv)<6):
		print " Error! \n Usage: python parseTopResultDB.py [outputFolder] [initFile] [outputPerFraction] [FDRCutoff] [logEMCutoff] [DisplayProteinNum] [debug/or noDebug]  \n"
		print " eg: python parseResultsDB.py /home/slin/sampleInputFiles/TG/mzXML_EM_output /home/slin/sampleInputFiles/TG/TAG_GRAPH_Tryp_CysCarbam_MetOx.ini No 0.01 2 5 debug"
		print " eg: python parseResultsDB.py /home/slin/sampleInputFiles/TG/mzML_EM_output /home/slin/sampleInputFiles/TG/TAG_GRAPH_Tryp_CysCarbam_MetOx.ini No 0.01 2 5 debug"
		sys.exit(2)
	debug=""
	newDic={}
	tgmodFile=PAR_DIR+"/resources/tagGraphModMassPair.pck"
	originalUnimodFile=PAR_DIR+"/resources/unimodDict_noLabels.pck"
	if os.path.isfile(tgmodFile):
		UnimodDictFile=open(tgmodFile,'rb')
		newDic=loadUnimodDict(UnimodDictFile)
	else:
		UnimodDictFile=open(originalUnimodFile,'rb')
		origUnimodDict=loadUnimodDict(UnimodDictFile)
		UnimodDict={k.lower():v for k,v in origUnimodDict.iteritems()}
		sdic=sorted(UnimodDict.items())
		for k,v in sdic:
			newDic[k]=[UnimodDict[k]['mass']]
	printMsg=True
	outputFolder=sys.argv[1]
	initFile=sys.argv[2]
	outputPerFraction=sys.argv[3].lower()
	FDRCutoff=sys.argv[4]
	logEMCutoff=sys.argv[5]
	DisplayProteinNum=sys.argv[6]
	if (len(sys.argv)==8):
		debug=sys.argv[7]
	if (debug=="debug"):
		printMsg=True
	for dirName,subdirList,fileList in os.walk(outputFolder): #walk and find topResult.tdv file
		for fname in fileList:
			if fname.find("_TopResults.tdv")>0:
				outfile=fname.replace("_TopResults.tdv","_TopResults")
				break
	resultsDBFile=outputFolder+"/"+"results.db"
	outFile=outputFolder+"/"+outfile
	fracMappingFile=open(outputFolder+'/fileFractionMapping.pck','rb')
	(fractionMapping,fractionList)=loadFractionInfo(fracMappingFile)
	(dictStaticMods,dictStaticModsName,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictDiffModsWithOriginal,NTERM,CTERM,dictStaticCNTermModsName)=loadInitFile(initFile)
	conn=lite.connect(resultsDBFile)
	conn.execute("PRAGMA max_page_count=max_page;")
	conn.execute("PRAGMA temp_store=2;")
	createPPMRangesFDRTables(conn)
	conn.commit()
	cursor=conn.cursor()
	populatePerfectFDRTableDummy(conn)
	for fractionID in fractionList:
		createPPMTables(conn)
		if printMsg:
			print ("Fration # %s" % fractionID)
		(idList,scanList,contextList,modContextList,modsList,deNovoPeptideList,obsMassList,EMProbabilityList,oneSubtractlg10EMList,chargeList,staticModsList,PPMList)=getScansPerFraction(conn,fractionID,dictDiffMods,dictStaticMods)
		if printMsg:
			print "done getScansPerFraction"
		(finalPeptideList,modsInfoList,newTheoMassList)=splitModsSubtition(newDic,modsList,modContextList,contextList,maxModLength,deNovoPeptideList,staticModsList,NTERM,CTERM,dictDiffMods,dictDiffModsMass,dictDiffModsName,dictStaticMods,dictStaticModsName,dictDiffModsWithOriginal,dictStaticCNTermModsName)
		if printMsg:
			print "done splitModsSubtition"
		(newPPMList,tenPowerNegOneSubtractlg10EMList)=CalculatePPMInfo(obsMassList,newTheoMassList,oneSubtractlg10EMList)
		if printMsg:
			print "done reCalPPMInfo"
		insertProcessedInfo(fractionID,idList,scanList,newPPMList,newTheoMassList,finalPeptideList,modsInfoList)
		if printMsg:
			print "done processDataUpload"
		insertBaseDataForPPMCalculation(conn,idList,chargeList,PPMList,newPPMList,EMProbabilityList,oneSubtractlg10EMList,finalPeptideList,tenPowerNegOneSubtractlg10EMList)
		if printMsg:
			print "done insertBaseDataForPPMCalculation"
		insertCleanDataPPM(conn)
		if printMsg:
			print "done insertCleanDataPPM"
		(ppmRangeUpperBound,ppmRangeLowerBound)=findPPMRange(conn,fractionID)
		if printMsg:
			print "done findPPMRange"
		##add extra step to move clean FDR to the keeper table.
		insertCleanDataFDR(conn,fractionID,ppmRangeUpperBound,ppmRangeLowerBound)
		if printMsg:
			print "done insertCleanDataFDR"
	calFalseDiscoveryRate(conn)
	createIndexes(conn)
	fout=open(outFile,'w') #dummy file handle
	writeFileCnt=0
	for fractionID in fractionList:
		print "write File: fraction: %s @ %s"%(fractionID,datetime.datetime.now())
		(writeFileCnt,fout)=writeScans(conn,fout,outFile,fractionMapping,fractionID,outputPerFraction,FDRCutoff,logEMCutoff,DisplayProteinNum,writeFileCnt)
	fout.close()
	os.remove(outFile)
	UnimodDictFile.close()
	keyToRemove = ('undefined mass shift','insertion','deletion')
	for curKey in keyToRemove:
		newDic.pop(curKey, None)
	with open(tgmodFile,'wb') as handle:
		pickle.dump(newDic,handle,protocol=0)
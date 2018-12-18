'''
Created on July 6, 2011

@author: Arun
'''

import numpy as np
from Constants import mods, aminoacids
import Constants

import pickle
import DataFile
import time
import copy
import re
import math

class ProbNetwork:

    ionRules = {
    'b': lambda PM, m, Nmod, Cmod: m + mods['H+'] + Nmod,
    'b-H2O': lambda PM, m, Nmod, Cmod: m + mods['H+'] - mods['H2O'] + Nmod,
    'b-NH3': lambda PM, m, Nmod, Cmod: m + mods['H+'] - mods['NH3'] + Nmod,
    'b-H2O-H2O': lambda PM, m, Nmod, Cmod: m + mods['H+'] - 2 * mods['H2O'] + Nmod,
    'b-H2O-NH3': lambda PM, m, Nmod, Cmod: m + mods['H+'] - mods['H2O'] - mods['NH3'] + Nmod,
    'a': lambda PM, m, Nmod, Cmod: m + mods['H+'] - mods['CO'] + Nmod,
    'a-H2O': lambda PM, m, Nmod, Cmod: m + mods['H+'] - mods['CO'] - mods['H2O'] + Nmod,
    'a-NH3': lambda PM, m, Nmod, Cmod: m + mods['H+'] - mods['CO'] - mods['NH3'] + Nmod,
    'y': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] + mods['H2O'] + Cmod,
    'y-H2O': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] + Cmod,
    'y-NH3': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] + mods['H2O'] - mods['NH3'] + Cmod,
    'y-H2O-H2O': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] - mods['H2O'] + Cmod,
    'y-H2O-NH3': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] - mods['NH3'] + Cmod,
    'c': lambda PM, m, Nmod, Cmod: m + mods['H+'] + mods['NH3'] + Nmod,
    'z': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] + mods['H2O'] - mods['NH2'] + Cmod,
    'c-1': lambda PM, m, Nmod, Cmod: m + mods['H+'] + mods['NH3'] - mods['Hyd'] + Nmod,
    'z+1': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] + mods['Hyd'] + mods['H2O'] - mods['NH2'] + Cmod,
    'c-H2O': lambda PM, m, Nmod, Cmod: m + mods['H+'] + mods['NH3'] - mods['H2O'] + Nmod,
    'b2': lambda PM, m, Nmod, Cmod: (m + 2*mods['H+'] + Nmod)/2,
    'y2': lambda PM, m, Nmod, Cmod: (PM - m + 2*mods['H+'] + mods['H2O'] + Cmod)/2,
    'b3': lambda PM, m, Nmod, Cmod: (m + 3*mods['H+'] + Nmod)/3,
    'y3': lambda PM, m, Nmod, Cmod: (PM - m + 3*mods['H+'] + mods['H2O'] + Cmod)/3
    }

    deltaRules = {
    'b': lambda PM, m, Nmod, Cmod: m - mods['H+'] - Nmod,
    'y': lambda PM, m, Nmod, Cmod: PM - m + mods['H+'] + mods['H2O'] + Cmod,
    }

    paramRules = {
    'AAClass': lambda self, PM, nodeInfo: ProbNetwork.getAAEquivClass(nodeInfo),
    'pos': lambda self, PM, nodeInfo: ProbNetwork.getBinnedRelPos(self, PM, nodeInfo),
    'cleavage': lambda self, PM, nodeInfo: ProbNetwork.getCleavageMotif(self, nodeInfo),
    'term': lambda self, PM, nodeInfo: ProbNetwork.getTerminus(nodeInfo)
    }

    binSizes = {
    'ion': lambda self: self._config['header']['bin_sizes']['ion'],
    'pos': lambda self: 5,
    'AAClass': lambda self: 16,
    'cleavage': lambda self: len(self._config['header']['enzyme']['specificity']) + 1,
    'term': lambda self: 3

    }

    def __init__(self, configIn, modelIn=None):
        with open(configIn, 'r') as fin:
            self._config = pickle.load(fin)
            self._dists = self._config['header']['distributions']
        if modelIn:
            with open(modelIn, 'r') as fin:
                self._model = pickle.load(fin)
        else:
            self.initializeCounts()
            print 'Entering training mode.'

    def getIonIntensities(self, spec):
        nodes = nodeInfoGen(spec.seq)
        Ints = []
        for nodeInfo in nodes:
            Inds = {}
            m = nodeInfo['prm']
            for ion in self._ions:
                ionMass = self.ionRules[ion](spec.pm, m, spec.Nmod, spec.Cmod)
                Inds[ion] = spec.getIntScore(ionMass)
            Ints.extend([Inds])

        return Ints

    def getIonIntensitiesForPRM(self, spec, prm):
        Ints = []
        for ion in self._ions:
            ionMass = self.ionRules[ion](spec.pm, prm, spec.Nmod, spec.Cmod)
            Ints += [spec.getIntScore(ionMass)]

        return Ints


    def getInds(self, spec, nodeInfo, pm, dist):
        Inds = {}

        for ion in self._config[dist]['ions']:
            ionMass = self.ionRules[ion](pm, nodeInfo['prm'], 0, 0)
            Inds[ion] = spec.getIntScore(ionMass)


        for param in self._config[dist]['params']:
            paramBin = self.paramRules[param](self, pm, nodeInfo)
            Inds[param] = paramBin

        return Inds

    def getModelProb(self, Inds, dist):
        prob = 1
        #print Inds
        try:
            for attr in self._model[dist].keys():
                conds = self._config[dist]['model'][attr]
                attrLoc = (Inds[attr],) + tuple(Inds[cond] for cond in conds)
                prob *= self._model[dist][attr][attrLoc]

            return prob
        except AttributeError:
            print 'No model loaded!'
            # Return 1 if distribution isn't in model
        except KeyError:
            return 1

    def getNoiseProb(self, spec, nodeInfo, pm, dist):
        prob = 1
        m = nodeInfo['prm']
        for ion in self._config[dist]['ions']:
            ionMass = self.ionRules[ion](pm, m, 0, 0)
            prob *= spec.getNoiseProb(ionMass)

        return prob


    def initializeCounts(self):
        self._counts = {}
        print self._config
        for dist in self._dists:
            self._counts[dist] = {}
            for attr in self._config[dist]['model'].keys():
                if attr in self.ionRules:
                    shape = [self.binSizes['ion'](self)]
                elif attr in self.paramRules:
                    shape = [self.binSizes[attr](self)]
                else:
                    raise ValueError('%s is not a valid model attribute. Try one of the following: %s' % (attr, self.binSizes.keys()))

                conds = self._config[dist]['model'][attr]
                for cond in conds:
                    if cond in self.ionRules:
                        shape.extend([self.binSizes['ion'](self)])
                    else:
                        shape.extend([self.binSizes[cond](self)])

                #add uniform prior
                self._counts[dist][attr] = np.ones(shape)

    def addToCounts(self, spec, dist):
        try:
            nodes = Constants.nodeInfoGen(spec.seq, considerTerminalMods=True, addTerminalNodes=True)
            for nodeInfo in nodes:
                Inds = self.getInds(spec, nodeInfo, spec.pm, dist)
                try:
                    for attr in self._counts[dist].keys():
                        conds = self._config[dist]['model'][attr]
                        attrLoc = (Inds[attr],) + tuple(Inds[cond] for cond in conds)
                        self._counts[dist][attr][attrLoc] += 1
                except KeyError:
                    pass

        except AttributeError:
            print 'No counts matrix detected. Not in training mode!'

    def getModelFromCounts(self, modelOut=None, comment=None):
        try:
            self._model = copy.deepcopy(self._counts)
            for dist in self._dists:
                for attr in self._model[dist].keys():
                    total = 0
                    for i in range(self._model[dist][attr].shape[0]):
                        total += self._model[dist][attr][i]

                    for i in range(self._model[dist][attr].shape[0]):
                        self._model[dist][attr][i] = self._model[dist][attr][i] / total

            self._model['header'] = self._config['header']
            if comment:
                self._model['header']['info'] += ' ' + comment


            if modelOut:
                with open(modelOut, 'w') as fout:
                    pickle.dump(self._model, fout)
            return self._model
        except AttributeError:
            print 'No counts matrix detected. Not in training mode!'

    def getModelHyperParameters(self):
        return self._model['hyper_parameters']

    def getHyperParameters(self, pairConfigName):
        try:
            return self._model['hyper_parameters'][pairConfigName]
        except KeyError:
            raise KeyError("No Optimal Hyper Parameters recorded for pair configuration %s" % (pairConfigName,))

    def getCleavageMotif(self, nodeInfo):
        specificity = self._config['header']['enzyme']['specificity']

        term = ProbNetwork.getTerminus(nodeInfo)

        for i, specRule in enumerate(specificity):
            if term != 1:
                CTrue = bool(re.match(specRule[1], nodeInfo['lattAA']))
            else:
                CTrue = True

            if term != 0:
                NTrue = bool(re.match(specRule[0], nodeInfo['formAA']))
            else:
                NTrue = True

            if NTrue and CTrue:
                return i

        return i + 1


    @staticmethod
    def getBinnedRelPos(self, PM, nodeInfo):
        m = nodeInfo['prm']
        numBins = ProbNetwork.binSizes['pos'](self)
        return min(np.floor(numBins * m / PM), numBins - 1)

    @staticmethod
    def getAAEquivClass(nodeInfo):
        l = nodeInfo['lattAA']
        f = nodeInfo['formAA']
        if l == 'P':
            if f == 'P':
                return 15
            else:
                return 0
        elif f == 'P':
            return 1
        elif l == 'G':
            if f == 'G':
                return 15
            else:
                return 2
        elif f == 'G':
            return 3
        elif l == 'R' or l == 'K':
            return 4
        elif f == 'H':
            return 5
        elif l == 'H':
            return 6
        elif f == 'D' or f == 'E':
            return 7
        elif l == 'D' or l == 'E':
            return 8
        elif f == 'J' or f == 'V':
            return 9
        elif l == 'J' or l == 'V':
            return 10
        elif f == 'S' or f == 'T':
            return 11
        elif l == 'S' or l == 'T':
            return 12
        elif f == 'N':
            return 13
        elif l == 'N':
            return 14
        else:
            return 15

    @staticmethod
    def getTerminus(nodeInfo):
        if nodeInfo['formAA'] == None:
            return 0
        elif nodeInfo['lattAA'] == None:
            return 1
        else:
            return 2

    @staticmethod
    def printIons(precMass, prm, NMod=0, CMod=0):
        PM = precMass - mods['H+'] - mods['H2O']
        for ion in ProbNetwork.ionRules:
            print ion + ': ' + str(ProbNetwork.ionRules[ion](PM, prm, NMod, CMod))


class Spectrum:

    def __init__(self, probNet, mH, Nmod, Cmod, spectrum, epsilon=0.02, sequence=None, intEp=5, resolution=None, windowSize=100, useMemo=False):
        self.e = epsilon
        self.probNet = probNet
        self._mH = mH
        self.pm = mH - mods['H2O'] - mods['H+'] - Nmod - Cmod
        self.Nmod = Nmod
        self.Cmod = Cmod
        self.origSpec = spectrum
        if not resolution:
            self._res = epsilon / intEp
            self._intEp = intEp
        else:
            self._res = min(resolution, epsilon / 2)
            self._intEp = int(np.ceil(epsilon / self._res))
        self._spec = self.prepareSpectrum(spectrum)
        self._window = windowSize
        if sequence:
            self.seq = sequence

        self._useMemo = useMemo
        if useMemo:
            self._memo = {}

    def initializeNoiseModel(self):
        masses = self.origSpec[:, 0]
        masses = np.sort(masses)
        noiseVector = np.zeros(ProbNetwork.binSizes['ion'](self.probNet))
        noiseCounts = np.zeros((np.ceil(masses[-1] / self._window) + 1, ProbNetwork.binSizes['ion'](self.probNet)))

        currWind = 0
        for mass in masses:
            if mass > (currWind + 1) * self._window:
                noiseCounts[currWind, :] = noiseVector
                noiseVector = np.zeros(ProbNetwork.binSizes['ion'](self.probNet))
                currWind = np.floor(mass / self._window)

            noiseVector[self.getIntScore(mass)] += 1

        noiseCounts[currWind, :] = noiseVector
        self.noiseModel = np.zeros(noiseCounts.shape)
        for window in range(noiseCounts.shape[0]):
            for intensity in range(noiseCounts.shape[1]):
                if intensity == 0:
                    prob = 1
                else:
                    prob = noiseCounts[window][intensity] * self.e / self._window

                for j in range(intensity + 1, noiseCounts.shape[1]):
                    prob *= 1 - (noiseCounts[window][j] * self.e / self._window)

                self.noiseModel[window][intensity] = prob

        #assume prior of 1 count in all zero slots
        self.noiseModel[np.where(self.noiseModel == 0)] = self.e / self._window

    def prepareSpectrum(self, spectrum):
        maxMass = self.hashMass(self._mH)
        #keep extra row for possible peak metadata
        spec = np.zeros((maxMass + self._intEp, 2), dtype=np.int8)

        normInts = self.intensityScoreLogTIC(spectrum)
        for i, pair in enumerate(spectrum):
            masses = self.getBuckets(pair[0])
            try:
                for mass in masses:
                    spec[mass][0] = max(normInts[i], spec[mass][0])
            except IndexError:
                pass

        return spec

    def hashSpectrum(self, spectrum):
        maxMass = self.hashMass(self._mH)
        #keep extra row for possible peak metadata
        spec = np.zeros((maxMass + self._intEp, 2), dtype=np.int8)

        for pair in spectrum:
            masses = self.getBuckets(pair[0])
            try:
                for mass in masses:
                    spec[mass][0] = max(pair[1], spec[mass][0])
            except IndexError:
                pass

        return spec

    def useToTrain(self, dist):
        if hasattr(self, 'seq'):
            self.probNet.addToCounts(self, dist)
        else:
            print 'Must provide sequence to use for training!'

    def getPScore(self, revMap=None, ambigEdges=None):
        self.initializeNoiseModel()
        pScore = 0
        if hasattr(self, 'seq'):
            prmLadder = getPRMLadder(self.seq, 0, ambigEdges=ambigEdges, addEnds=False)
            for i, prm in enumerate(prmLadder):
                if revMap:
                    pScore += self.getNodeScore(prm=prm, formAA=revMap[self.seq[i]], lattAA=revMap[self.seq[i + 1]])
                else:
                    pScore += self.getNodeScore(prm=prm, formAA=self.seq[i], lattAA=self.seq[i + 1])
        else:
            print 'Must provide sequence for scoring!'
            pScore = False

        return pScore


    def getIntScore(self, mass):
        m = self.hashMass(mass)
        try:
            return self._spec[m][0]
        except IndexError:
            return 0

    def getPriorScore(self, **nodeInfo):
        Inds = self.probNet.getInds(self, nodeInfo, getIons=False)
        return np.log(self.probNet.getModelProb(Inds, 'prior'))

    def getNodeScore(self, nodeInfo, pm, dist):

        # See if node has already been scored
        if self._useMemo:
            try:
                return self._memo[self.hashMass(nodeInfo['prm'])]
            except KeyError:
                pass

        Inds = self.probNet.getInds(self, nodeInfo, pm, dist)
        CIDProb = self.probNet.getModelProb(Inds, dist)
        NoiseProb = self.probNet.getNoiseProb(self, nodeInfo, pm, dist)

#        print nodeInfo, Inds, np.log(CIDProb * priorProb / NoiseProb)
#        return np.log(CIDProb * priorProb / NoiseProb)
        nodeScore = np.log(CIDProb/NoiseProb)

        if self._useMemo:
            masses = self.getBuckets(nodeInfo['prm'])
            for mass in masses:
                self._memo[mass] = nodeScore

        return nodeScore

    def getPRMScore(self, prm):
        nodeInfo = {'prm': prm}

        Inds = self.probNet.getInds(self, prm)
        CIDProb = self.probNet.getModelProb(Inds, 'likelihood')
        NoiseProb = self.probNet.getNoiseProb(self, nodeInfo)

        nodeScore = np.log(CIDProb/NoiseProb)
        return nodeScore if nodeScore > 0 else 3*nodeScore

    def getNoiseProb(self, mass):
        window = np.floor(mass / self._window)
        intensity = self.getIntScore(mass)
        if window >= 0 and window < self.noiseModel.shape[0]:
            if self.noiseModel[window][intensity] == 0:
                print mass, intensity, window
            return self.noiseModel[window][intensity]
        else:
            if intensity == 0:
                return 1
            else:
                return 0

    """
    def recordProb(self, mass, prob):
        intMass = self.hashMass(mass)
        if intMass < self._spec.shape[0]:
            if self._spec[intMass][0] == 0:
                self._spec[intMass][1] = max(self._spec[intMass][1], prob)
            else:
                peakMasses = self.getPeak(intMass)
                for mass in peakMasses:
                    self._spec[mass][1] = max(self._spec[mass][1], prob)
    """

    # linear intensity binning based on fraction of total TIC
    def intensityScoreLinearTIC(self, spectrum):
        TIC = np.sum(spectrum[:,1])
        numBins = ProbNetwork.binSizes['ion'](self.probNet)

        normIntArr = np.ceil(spectrum[:,1]/TIC * (numBins - 1))
        return normIntArr

    # log intensity binning based on fraction of total TIC (anything above 10% of TIC goes in highest bin)
    def intensityScoreLogTIC(self, spectrum):
        TIC = np.sum(spectrum[:,1])
        numBins = ProbNetwork.binSizes['ion'](self.probNet)

        # log base 10 makes score of any peak above 10% of TIC  = numBins - 1
        normIntArr = numBins - 1 + np.ceil(np.log10(spectrum[:,1]/TIC))
        # Make sure that minimum normalized intensity is 1
        return np.maximum(np.ones(spectrum.shape[0]), normIntArr)

    # binning based on percentage (bins are equally sized percentage intervals)
    def intensityScoreLinearPercent(self, spectrum):
        intOrd = np.argsort(spectrum[:,1])
        numBins = ProbNetwork.binSizes['ion'](self.probNet)

        normIntArr = np.ceil((1 - intOrd.astype(np.float)/spectrum.shape[0]) * (numBins - 1))
        return normIntArr

    # binning based on log of percentage (top 5% of peaks are in highest bin). This binning scheme is valid to a maximum of 6 bins, otherwise second highest bin has less peaks than highest bin
    def intensityScoreLogPercent(self, spectrum):
        intOrd = np.argsort(spectrum[:,1]) + 1
        numBins = ProbNetwork.binSizes['ion'](self.probNet)

        base = math.pow(0.05, -1/(numBins - 2))
        normIntArr = -1 * np.floor(np.log(intOrd.astype(np.float)/(spectrum.shape[0] + 1))/np.log(base))
        return np.minimum(normIntArr, np.ones(spectrum.shape[0]) * (numBins - 1))

    # log intensity binning based on fold difference over baseline intensity (average of intensity of lower third of spectrum)
    def intensityScoreLogBaseline(self, spectrum):
        ordInts = np.sort(spectrum[:, 1])
        baseline = np.average(ordInts[0:np.ceil(spectrum.shape[0] / 3)])
        numBins = ProbNetwork.binSizes['ion'](self.probNet)

        normIntArr = np.ceil(np.log10(spectrum[:,1]/baseline)) + 1
        # make sure minimum intensity is 1 and maximum intensity is numBins - 1
        normIntArr = np.minimum(np.maximum(np.ones(spectrum.shape[0]), normIntArr), np.ones(spectrum.shape[0])*(numBins-1))

        return normIntArr

    def intensityScore(self, intensity, baseline):
        numBins = ProbNetwork.binSizes['ion'](self.probNet)
        normInt = np.ceil(np.log10(intensity / baseline))
        if normInt + 1 >= numBins - 1:
            return numBins - 1
        elif normInt + 1 <= 1:
            return 1
        else:
            return normInt + 1

    def getPeak(self, intMass):
        peakInt = self._spec[intMass][0]
        if peakInt > 0:
            if peakInt == self._spec[intMass + 1][0]:
                return (intMass, intMass + 1)
            elif self._spec[intMass - 1][0] == peakInt:
                return (intMass - 1, intMass)
            else:
                print "ERROR: Spectrum improperly hashed!"
        else:
            print "No peak detected at position: ", intMass

    def getBuckets(self, mass):
        hMass = self.hashMass(mass)
        return [m for m in range(hMass - self._intEp, hMass + self._intEp + 1)]

    def hashMass(self, mass):
        return int(np.round(mass / self._res))

    def unhashMass(self, hMass):
        return hMass * self._res

    def getHypothesisProb(self, sequence):
        prm = 0
        for j in range(len(sequence) - 1):
            prm += aminoacids[sequence[j]][2]
            formAA = sequence[j]
            lattAA = sequence[j + 1]
            self.probNet.getCIDProb(self, {'prm': prm, 'formAA': formAA, 'lattAA': lattAA}, recordProbs=True)

        prob = 1
        i = 1
        while i < self._spec.shape[0]:
            if self._spec[i, 1] > 0:
                prob *= self._spec[i, 1]
                #skip past rest of peak buckets
                while self._spec[i + 1, 1] == self._spec[i, 1]:
                    i += 1

            i += 1

        if prob != 1:
            return prob
        else:
            return 0


    def clearSpectrumProb(self):
        self._spec[:, 1] = 0

    def correctParentMass(self, halfWindow=0.1):
        masses = self.origSpec[:, 0]
        numEp = int(np.floor(halfWindow / self.e))

        maxOverlap = 0
        massCorrection = 0
        for i in range(-1 * numEp, numEp + 1):
            testMH = self._mH + self.e * i
            overlap = 0
            revMasses = testMH - masses + mods['H+']
            for mass in revMasses:
                if self._spec[self.hashMass(mass)][0] > 0:
                    overlap += 1

            #print i, overlap
            if overlap > maxOverlap:
                maxOverlap = overlap
                massCorrection = self.e * i
            elif overlap == maxOverlap:
                if np.abs(self.e * i) < np.abs(massCorrection):
                    massCorrection = self.e * i
            else:
                pass

        self._mH += massCorrection
        self.pm += massCorrection
        return massCorrection

def nodeInfoGen(seq, startMass=0, aaDict=aminoacids, addTerminalNodes=False):
    i = 0
    AAs = []
    aa = ''
    prm = startMass

    while i < len(seq):
        try:
            while i < len(seq):
                aamass = aaDict[aa + seq[i]][2]
                aa += seq[i]
                i += 1
        except KeyError:
            if len(AAs) > 0:
                if AAs[-1] != '':
                    yield {'prm': prm, 'lattAA': aa, 'formAA': AAs[-1]}
                else:
                    raise KeyError('%s is not a valid amino acid' % (aa + seq[i],))
            elif addTerminalNodes:
                yield {'prm': prm, 'lattAA': aa, 'formAA': None}

            AAs += [aa]
            aa = ''
            prm += aamass
            aamass = 0

    if len(AAs) > 0:
        try:
            aamass = aaDict[aa][2]
            yield {'prm': prm, 'lattAA': aa, 'formAA': AAs[-1]}
        except KeyError:
            raise KeyError('%s is not a valid amino acid' % (aa,))

    if addTerminalNodes and aa:
        if len(AAs) == 0:
            yield {'prm': prm, 'lattAA': aa, 'formAA': None}

        yield {'prm': prm + aamass, 'lattAA': None, 'formAA': aa}

def getPRMLadder(seq, startMass, ambigAA='X', ambigEdges=None, epsilon=1, addEnds=True, considerTerminalMods=True):
    nodeGen = Constants.nodeInfoGen(seq, startMass=startMass, addTerminalNodes=addEnds, considerTerminalMods=considerTerminalMods)
    seqIndex = 1
    PRMLadder = []

    while True:
        try:
            node = nodeGen.next()
            PRMLadder.extend([node['prm']])
        except KeyError:
            if seq[seqIndex - 1] == ambigAA:
                edge = ambigEdges[0]
                ambigEdges = ambigEdges[1:]
                if np.abs(edge[0] - (PRMLadder[-1] if len(PRMLadder) > 0 else startMass)) < epsilon:
                    PRMLadder.extend([edge[1]])
                    PRMLadder.extend(getPRMLadder(seq=seq[seqIndex:], startMass=edge[1], ambigAA=ambigAA, ambigEdges=ambigEdges, epsilon=epsilon, addEnds=False))
                    break
                else:
                    print 'ERROR: Ambiguous edges do not correspond to ambiguous regions of sequence for PRM =', PRMLadder[-1] if len(PRMLadder) > 0 else startMass, 'and ambiguous edges', edge
                    return False
            else:
                raise ValueError('Unknown amino acid found %s' % node['lattAA'])
        except StopIteration:
            break

    return PRMLadder

if __name__ == '__main__':
    #aminoacids['C'] = (aminoacids['C'][0], aminoacids['C'][1], aminoacids['C'][2] + mods['Carbamidomethyl'], aminoacids['C'][3])

    #PN = ProbNetwork('config.txt', 'model.txt')
    """
    scanInfo = DataFile.getScanInfo('adevabhaktuni_1310166306.csv')

    dirPath = 'C:\\Users\\Arun\\Pythonprojects\\DeNovoSequencing\\LF2_short_HCD+CID_ath001862_244\\'
    dtaNames = DataFile.getDTAFNamesInDir(dirPath)

    scansIter = iter(dtaNames)
    currScanInfo = scansIter.next()
    for dta in dtaNames:
        precMass = DataFile.getPrecMassAndCharge(dta)[0]
        spectra = DataFile.getMassIntPairs(dta)
        S = Spectrum(PN, precMass, 0.0, 0.0, spectra)
        corr = S.correctParentMass()
        if np.abs(corr) > 0.04:
            print dta, corr

    """
    paramsDict = DataFile.parseParams('/home/arun/Documents/LADS_SILAC_Trypsin.ini')
    print getPRMLadder('A', 0)
    """
    heavyPath = "C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3367.3367.1.dta"
    lightPath = "C:\\Users\\Arun\\DropBox\\SpectraCorrelation\\244.3383.3383.1.dta"
    heavyPairs = DataFile.getMassIntPairs(heavyPath)
    lightPairs = DataFile.getMassIntPairs(lightPath)
    heavyPrecMass, heavyCharge = DataFile.getPrecMassAndCharge(heavyPath)
    lightPrecMass, lightCharge = DataFile.getPrecMassAndCharge(lightPath)

    heavySpec = Spectrum(PN, heavyPrecMass, 0, mods['*'], heavyPairs)
    lightSpec = Spectrum(PN, lightPrecMass, 0, 0, lightPairs)
    heavySpec.initializeNoiseModel()
    lightSpec.initializeNoiseModel()
    print heavySpec.noiseModel
    print lightSpec.noiseModel
    t1 = time.time()
    print heavySpec.getNodeScore(prm=955.458, formAA='H', lattAA='P')
    print lightSpec.getNodeScore(prm=955.458, formAA='H', lattAA='P')
    t2 = time.time()
    print t2-t1
    """

    """
    paramsDict = DataFile.parseParams('/home/arun/Documents/LADS_SILAC_Trypsin.ini')
    PN = ProbNetwork('/home/arun/Documents/LysC_likelihood_prior_config.txt')
    PN.initializeCounts()

    dirPath = '/home/arun/Proteomics_Data/LF2_short_HCD+CID_ath001862_244/'
    dta = '244.1266.1266.1.dta'
    precMass = DataFile.getPrecMassAndCharge(dirPath+dta)[0]
    massIntPairs = DataFile.getMassIntPairs(dirPath+dta)
    epsilon=20*precMass*10**-6
    S = Spectrum(PN, precMass, 0.0, 0.0, massIntPairs, epsilon=epsilon, sequence='RVAEDDEDDDVDTK*K*')
    S.useToTrain()
    """


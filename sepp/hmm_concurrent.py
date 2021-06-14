import copy
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib
from subprocess import run

from Bio import SeqIO

from sepp import get_logger
from sepp.filemgr import get_root_temp_dir
from sepp.scheduler import JobPool
from sepp.jobs import HMMAlignJob,HMMBuildJob,HMMSearchJob


hmmSeqFile = ''
queryName = ''
trueAlignment = ''
predictionName = ''
dirName = ''
dataSetName = ''

_LOG = get_logger(__name__)

_LOG = get_logger(__name__)

def ensureFolder(fileName):
    pathlib.Path(fileName).parents[0].mkdir(parents=True, exist_ok=True)
    # def removeEnd(fileName):
    #     fileName = '/'.join(fileName.split('/')[:-1])
    #     return (fileName)
    # fileName = removeEnd(fileName)
    # if os.path.exists(fileName) == False:
    #     toMake = []
    #     while os.path.exists(fileName) == False:
    #         toMake.append(fileName)
    #         fileName = removeEnd(fileName)
    #     toMake = np.array(toMake)[-1::-1]
    #     for a in range(len(toMake)):
    #         os.mkdir(toMake[a])

def seqToArray(seqs):
    seqArray = []
    for seq in seqs:
        seqArray.append(list(seq))
    seqArray = np.array(seqArray)
    return seqArray

def arrayToSeq(array):
    seqs = []
    for a in range(len(array)):
        seq = ''.join(list(array[a]))
        seqs.append(seq)
    seqs = np.array(seqs)
    return seqs

def setAllFileNames(givenHSF, givenQN, givenTA, givenPN, givenDir, givenDSN):
    global hmmSeqFile 
    hmmSeqFile = givenHSF
    global queryName 
    queryName = givenQN
    global trueAlignment
    trueAlignment = givenTA
    global predictionName
    predictionName = givenPN
    global dirName
    dirName = givenDir
    global dataSetName
    dataSetName = givenDSN

def giveAllFileNames():
    global hmmSeqFile
    global queryName
    global trueAlignment
    global predictionName
    global dataSetName

    fileNames = [hmmSeqFile, queryName, trueAlignment, predictionName, dataSetName]
    return fileNames

def giveQueries():
    queryName = giveQueryFileName()
    file_obj = open(queryName)
    queries = []
    count1 = 0
    for line_number, i in enumerate(file_obj):
        if line_number % 2 == 0:
            queries.append(i[1:-1])
    queries = np.array(queries)
    return queries

def removeMultInvert(seqs):
    lowercase = np.array(list('abcdefghijklmnopqrstuvwxyz'))
    seqArray = seqToArray(seqs)
    seqArray = seqArray.T
    seqArray2 = []
    for a in range(0, len(seqArray)):
        seq = seqArray[a]
        argsLower = np.argwhere(np.isin(seq, lowercase))[:, 0]
        if argsLower.shape[0] == 0:
            seqArray2.append(np.copy(seqArray[a]))
        else:
            for b in range(len(argsLower)):
                seq2 = np.zeros(len(seq)).astype(str)
                seq2[:] = '-'
                seq2[argsLower[b]] = seq[argsLower[b]]
                seqArray2.append(np.copy(seq2))
    seqArray2 = np.array(seqArray2).T
    seqs = arrayToSeq(seqArray2)
    return seqs

def removeEmptyColumns(seqs):
    seqArray = []
    for seq in seqs:
        seqArray.append(copy.copy(list(seq)))
    seqArray = np.array(seqArray).T
    seqArray2 = []
    for a in range(0, len(seqArray)):
        size1 = seqArray[a][seqArray[a] != '-'].shape[0]
        if size1 != 0:
            seqArray2.append(np.copy(seqArray[a]))
    seqArray2 = np.array(seqArray2).T
    seq2 = []
    for a in range(0, len(seqArray2)):
        str1 = ''.join(list(seqArray2[a]))
        seq2.append(str1)
    seq2 = np.array(seq2)
    return seq2

def saveFasta(name, formatedData):
    ensureFolder(name)
    with open(name, 'w') as f:
        f.write('\n'.join(formatedData))

def loadFastaFormat(name):
    file_obj = open(name)
    data = []
    name = ''
    for line_number, i in enumerate(file_obj):
        name = name + i
        if line_number % 2 == 1:
            name = name[:-1]
            data.append(name)
            name = ''
    return data

def loadFastaBasic(name):
    file_obj = open(name)
    keys = []
    seqs = []
    count1 = 0
    for line_number, i in enumerate(file_obj):
        if line_number % 2 == 0:
            keys.append(i[1:-1])
        else:
            seq1 = i.replace('\n', '')
            seqs.append(seq1)
    keys, seqs = np.array(keys), np.array(seqs)
    return keys, seqs

def saveFastaBasic(name, keys, seqs):
    data = []
    for a in range(0, len(keys)):
        data.append('>' + keys[a] + '\n' + seqs[a])
    data[-1] = data[-1] + '\n'
    saveFasta(name, data)

def saveSequenceFileNames():
    fileNames = giveAllFileNames()
    src0 = fileNames[0]
    sequenceFileNames = []
    Anums = []
    namesA = os.listdir(src0)
    for nameA in namesA:
        if nameA[:2] == 'A_':
            Anum = nameA.split('_')
            Anum = int(Anum[-1])
            src1 = src0 + nameA + '/'
            names = os.listdir(src1)
            for name in names:
                if name[-5:] == 'fasta':
                    src2 = src1 + name
                    sequenceFileNames.append(src2)
                    Anums.append(Anum)
    sequenceFileNames = np.array(sequenceFileNames)
    Anums = np.array(Anums)
    removalArgs = []
    for a in range(len(sequenceFileNames)):
        for b in range(len(sequenceFileNames)):
            if b > a:
                keys1, _ = loadFastaBasic(sequenceFileNames[a])
                keys2, _ = loadFastaBasic(sequenceFileNames[b])
                if len(keys1) == len(keys2):
                    if len(np.intersect1d(keys1, keys2)) == len(keys1):
                        removalArgs.append(b)
    keepArgs = np.arange(len(sequenceFileNames))
    removalArgs = np.array(removalArgs)
    keepArgs = keepArgs[np.isin(keepArgs, removalArgs) == False]

    Anums = Anums[keepArgs]
    sequenceFileNames = sequenceFileNames[keepArgs]
    sequenceFileNames = sequenceFileNames[np.argsort(Anums)]

    inputFiles = giveAllFileNames()
    dataFileName = inputFiles[4]

    sequenceFileNames_saveName = get_root_temp_dir() + '/data/internalData/' + dataFileName + '/seqFileNames/names.npy'
    ensureFolder(sequenceFileNames_saveName)
    np.save(sequenceFileNames_saveName, sequenceFileNames)

def giveSequenceFileNames():
    inputFiles = giveAllFileNames()
    dataFileName = inputFiles[4]
    sequenceFileNames_saveName = get_root_temp_dir() + '/data/internalData/' + dataFileName + '/seqFileNames/names.npy'
    sequenceFileNames = np.load(sequenceFileNames_saveName)
    return sequenceFileNames

def giveQueryFileName():
    fileNames = giveAllFileNames()
    queryName = fileNames[1]
    return queryName

def saveInitialHMM(abstract_algorithm, minSize=1):
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        src = sequenceFileNames[a]
        inputFiles = giveAllFileNames()
        dataFileName = inputFiles[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFileName + "/initialHMM/test" + str(a) + ".hmm"
        ensureFolder(hmmName)
        addHMMBuildJob(abstract_algorithm, hmmName, src)
    JobPool().wait_for_all_jobs()

def giveHMMversion():
    hmmVersion = ''
    return hmmVersion

def addHMMSearchJob(abstract_algorithm, hmmName, queryName, scoreName):
    ensureFolder(scoreName)
    hmmsearch_job = HMMSearchJob(**vars(abstract_algorithm.options.hmmsearch))
    hmmsearch_job.setup(hmmName, queryName, scoreName, elim=abstract_algorithm.elim)
    hmmsearch_job.results_on_temp = False
    JobPool().enqueue_job(hmmsearch_job)

def addHMMBuildJob(abstract_algorithm, hmmName, seqName):
    ensureFolder(hmmName)
    hmmbuild_job = HMMBuildJob()
    hmmbuild_job.setup(seqName, hmmName, **vars(abstract_algorithm.options.hmmbuild))
    JobPool().enqueue_job(hmmbuild_job)

def addHMMAlignJob(abstract_algorithm, hmmName, queryName, predictionName):
    ensureFolder(predictionName)
    hmmalign_job = HMMAlignJob(**vars(abstract_algorithm.options.hmmalign))
    # _LOG.debug(str(abstract_algorithm.options.hmmalign))
    hmmalign_job.setup(hmmName, queryName, predictionName)
    JobPool().enqueue_job(hmmalign_job)

def runHMMbuild(abstract_algorithm, hmmName, seqName):
    # ensureFolder(hmmName)
    hmmbuild_job = HMMBuildJob()
    hmmbuild_job.setup(seqName, hmmName, **vars(abstract_algorithm.options.hmmbuild))
    hmmbuild_job.run()

def runHMMalign(abstract_algorithm, hmmName, queryName, predictionName):
    # ensureFolder(predictionName)
    hmmalign_job = HMMAlignJob()
    hmmalign_job.setup(hmmName, queryName, predictionName, **vars(abstract_algorithm.options.hmmalign))
    hmmalign_job.run()

def saveScoreSimple(abstract_algorithm, hmmNames, queryNames, scoreName):
    queryDict = {}
    queries = giveQueries()
    for queryNum in range(0, len(queries)):
        query = queries[queryNum]
        queryDict[query] = queryNum
    sequenceFileNames = giveSequenceFileNames()
    hmmNum = len(hmmNames)
    data = np.zeros((len(queries), hmmNum))
    randInt = str(np.random.randint(100000000000))
    for a in range(0, len(hmmNames)):
        hmmName = hmmNames[a]
        scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + randInt + ".txt"
        queryName = queryNames[a]
        ensureFolder(scoreNameTemp)
        addHMMAlignJob(abstract_algorithm, hmmName, queryName, scoreNameTemp)
    JobPool().wait_for_all_jobs()

    for a in range(0, len(hmmNames)):
        hmmName = hmmNames[a]
        scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + randInt + ".txt"
        queryName = queryNames[a]
        ensureFolder(scoreNameTemp)
        dataOld = []
        file_obj = open(scoreNameTemp)
        count1 = 0
        for line_number, i in enumerate(file_obj):
            ar = i.split('  ')
            listI = []
            bitscoreList = []
            if line_number == 1:
                posScore = i.find('score')
            if line_number == 2:
                argSpace = np.argwhere(np.array(list(i)) == ' ')[:, 0]
            if (len(ar) >= 30) and (line_number >= 3):
                sequenceName = i[:20].replace(' ', '')
                bitScore = float(i[posScore-1:posScore+5].replace(' ', ''))
                dataOld.append([sequenceName, bitScore])

        dataOld = np.array(dataOld)
        for b in range(0, len(dataOld)):
            queryN = queryDict[dataOld[b, 0]]
            data[queryN, a] = float(dataOld[b, 1])
    ensureFolder(scoreName)
    np.save(scoreName, data)
    os.remove(scoreNameTemp)

def saveNewScores(abstract_algorithm, strategyName):
    dataFolderName = giveAllFileNames()[4]
    queryNameO = giveQueryFileName()
    queryData = loadFastaFormat(queryNameO)
    queryNames = []
    hmmNames = []

    queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")
    numHMM = int(np.max(queryToHmm[:, 1])+1)
    for b in range(0, numHMM):
        argsHMM = queryToHmm[:, 0][queryToHmm[:, 1] == b]
        dataNow = []
        keysNow = []
        for a in argsHMM:
             dataNow.append(queryData[a])
        dataFolderName = giveAllFileNames()[4]
        queryName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/inputQuery/'  + str(b) + '.fasta'
        ensureFolder(queryName)
        saveFasta(queryName, dataNow)
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(b) + ".hmm"
        ensureFolder(hmmName)
        queryNames.append(queryName)
        hmmNames.append(hmmName)
    dataFolderName = giveAllFileNames()[4]
    scoreName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/hmmScores/scoresFull/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(abstract_algorithm, hmmNames, queryNames, scoreName)

def saveScoreFromBool(abstract_algorithm, hmmNames, queryHMM, scoreName, noSave=False):
    queryDict = {}
    queries = giveQueries()
    for queryNum in range(0, len(queries)):
        query = queries[queryNum]
        queryDict[query] = queryNum
    sequenceFileNames = giveSequenceFileNames()
    hmmNum = len(hmmNames)
    data = np.zeros((len(queries), hmmNum))
    randInt1 = str(np.random.randint(100000000000))
    randInt2 = str(np.random.randint(100000000000))
    for a in range(0, len(hmmNames)):
        argsInclude = queryHMM[:, 0][queryHMM[:, 1] == a]
        if argsInclude.shape[0] > 0:
            hmmName = hmmNames[a]
            scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + randInt1 + ".txt"
            queryName = giveQueryFileName()
            queryData = loadFastaFormat(queryName)
            queryNameTemp = get_root_temp_dir() + '/data/temporaryStorage/temp_' + randInt2 + '.fasta'
            queryNow = []
            for b in argsInclude:
                queryNow.append(queryData[b])
            ensureFolder(queryNameTemp)
            saveFasta(queryNameTemp, queryNow)
            ensureFolder(scoreNameTemp)
            addHMMAlignJob(abstract_algorithm, hmmName, queryNameTemp, scoreNameTemp)
    JobPool().wait_for_all_jobs()
    for a in range(0, len(hmmNames)):
        argsInclude = queryHMM[:, 0][queryHMM[:, 1] == a]
        if argsInclude.shape[0] > 0:
            scoreNameTemp = get_root_temp_dir() + "/data/temporaryStorage/temp_" + randInt1 + ".txt"
            dataOld = []
            file_obj = open(scoreNameTemp)
            count1 = 0
            for line_number, i in enumerate(file_obj):
                ar = i.split('  ')
                listI = []
                bitscoreList = []
                if line_number == 1:
                    posScore = i.find('score')
                if line_number == 2:
                    argSpace = np.argwhere(np.array(list(i)) == ' ')[:, 0]
                if (len(ar) >= 30) and (line_number >= 3):
                    sequenceName = i[:20].replace(' ', '')
                    bitScore = float(i[posScore-1:posScore+5].replace(' ', ''))
                    dataOld.append([sequenceName, bitScore])
            dataOld = np.array(dataOld)
            for b in range(0, len(dataOld)):
                queryN = queryDict[dataOld[b, 0]]
                data[queryN, a] = float(dataOld[b, 1])
    os.remove(scoreNameTemp)
    os.remove(queryNameTemp)
    if noSave:
        return data
    else:
        ensureFolder(scoreName)
        np.save(scoreName, data)

def compareScores(strategyName):
    scoreStrat = np.load("%s/hmmScores/scoresFull/" % dirName + strategyName + "_full.npy")
    scoreUPP = np.load("%s/hmmScores/full.npy" % dirName)
    plt.plot(np.max(scoreStrat, axis=1)-np.max(scoreUPP, axis=1))
    plt.show()

def saveScoresOriginal(abstract_algorithm):
    sequenceFileNames = giveSequenceFileNames()
    hmmNames = []
    queryNames = []

    for a in range(0, len(sequenceFileNames)):
        dataFolderName = giveAllFileNames()[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/initialHMM/test" + str(a) + ".hmm"
        hmmNames.append(hmmName)
        queryName = giveQueryFileName()
        queryNames.append(queryName)

    dataFolderName = giveAllFileNames()[4]
    scoreName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(abstract_algorithm, hmmNames, queryNames, scoreName)

def processScores():
    ''' Note: is this deprecated
    '''
    sequenceFileNames = giveSequenceFileNames()
    hmmNames = []
    queryNames = []

    for a in range(0, len(sequenceFileNames)):
        dataFolderName = giveAllFileNames()[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/initialHMM/test" + str(a) + ".hmm"
        hmmNames.append(hmmName)
        queryName = giveQueryFileName()
        queryNames.append(queryName)

    dataFolderName = giveAllFileNames()[4]
    scoreName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(hmmNames, queryNames, scoreName)

def processScoresOld():
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        data = []
        file_obj = open("%s/hmmScores/raw/" % dirName + str(a) + ".txt")
        count1 = 0
        for line_number, i in enumerate(file_obj):
            ar = i.split('  ')
            if (len(ar) >= 30) and (line_number >= 3):
                bitScore = float(i[76:83].replace(' ', ''))
                sequenceName = i[:20].replace(' ', '')
                data.append([sequenceName, bitScore])
        data = np.array(data)
        np.save("%s/hmmScores/processed/" %  dirName + str(a) + ".npy", data)

def saveScoresBySeq():
    queryDict = {}
    queries = giveQueries()
    for queryNum in range(0, len(queries)):
        query = queries[queryNum]
        queryDict[query] = queryNum
    sequenceFileNames = giveSequenceFileNames()
    hmmNum = len(sequenceFileNames)
    data = np.zeros((len(queries), hmmNum))
    for a in range(0, hmmNum):
        dataFolderName = giveAllFileNames()[4]
        dataOld = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/processed/" + str(a) + ".npy")
        for b in range(0, len(dataOld)):
            queryN = queryDict[dataOld[b, 0]]
            data[queryN, a] = float(dataOld[b, 1])

    dataFolderName = giveAllFileNames()[4]
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy", data)

def saveDecomposition():
    dataFolderName = giveAllFileNames()[4]
    sets = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy", allow_pickle=True)
    treeData = np.zeros((len(sets), len(sets)))
    for a in range(0, len(sets)):
        for b in range(0, len(sets)):
            set1, set2 = np.array(sets[a]), np.array(sets[b])
            if np.intersect1d(set1, set2).shape[0] == set2.shape[0]:
                treeData[a, b] = 1
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/treeDecomp.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/treeDecomp.npy", treeData)

def saveHMMsets():
    sequenceFileNames = giveSequenceFileNames()
    sets = []
    for a in range(0, len(sequenceFileNames)):
        src = sequenceFileNames[a]
        set1 = []
        file_obj = open(src)
        count1 = 0
        for line_number, i in enumerate(file_obj):
            if line_number % 2 == 0:
                set1.append(i[1:-1])

        sets.append(set1)



    #print ("Q")

    #It gives warning. Fine for now but maybe change at some point.

    dataFolderName = giveAllFileNames()[4]

    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy", sets)

def findDecomposition():
    dataFolderName = giveAllFileNames()[4]
    treeData = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/treeDecomp.npy")
    return treeData

def findTreeEdges():
    treeData = findDecomposition()
    treeData[np.arange(treeData.shape[0]), np.arange(treeData.shape[0])] = 0
    edges = np.zeros(treeData.shape)
    treeSum = np.sum(treeData, axis=1)
    nodes = []
    for a in range(treeData.shape[0]):
        treeSumRemain = np.copy(treeSum)
        nodesArray = np.array(nodes).astype(int)
        treeSumRemain[nodesArray] = -1

        argMax = int(np.argmax(treeSumRemain))
        nodes.append(argMax)
        treeSumIn = np.copy(treeSum) + ((1-treeData[:, argMax]) * 1000000)
        argMin = int(np.argmin(treeSumIn))

        if treeData[argMin, argMax] == 1:
            edges[argMin, argMax] = 1
    return edges

def findFatherNode():
    edges = findTreeEdges()
    fatherNode = np.argmax(edges, axis=0)
    return fatherNode

def findBrotherNode():
    edges = findTreeEdges()
    fatherNode = findFatherNode()
    kidsOfFather = edges[fatherNode]
    kidsOfFather[np.arange(edges.shape[0]), np.arange(edges.shape[0])] = 0
    brotherNode = np.argmax(kidsOfFather, axis=1)
    return brotherNode

def interweive(data):
    M = data.shape[0]
    arange1 = np.arange(data.shape[1])
    arange2 = arange1 // M
    arange3 = arange1 % M
    data4 = data[arange3, arange2]
    return data4

def reallignHMMSeq():
    sequenceFileNames = giveSequenceFileNames()
    rootName = sequenceFileNames[0]
    rootKeys, rootSeqs = loadFastaBasic(rootName)

    for a in range(0, len(sequenceFileNames)):
        name = sequenceFileNames[a]
        nowKeys, nowSeqs = loadFastaBasic(name)
        argsIn = np.argwhere(np.isin(rootKeys, nowKeys))[:, 0]
        nowKeys, nowSeqs = rootKeys[argsIn], rootSeqs[argsIn]
        dataFolderName = giveAllFileNames()[4]
        ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSeqAlign/" + str(a) + '.fasta')
        saveFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSeqAlign/" + str(a) + '.fasta', nowKeys, nowSeqs)

def compareHMM():
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        HMM1 = "%s/initialHMM/testStandardAlign" % dirName + str(a) + ".hmm"
        HMM2 = "%s/initialHMM/test" % dirName + str(a) + ".hmm"
        cmd = 'diff ' + HMM1 + ' ' + HMM2
        print (cmd)

def saveAdjustedScore():
    dataFolderName = giveAllFileNames()[4]
    sets = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmSets.npy", allow_pickle=True)
    sizes1 = []
    for a in range(len(sets)):
        sizes1.append(len(sets[a]))
    sizes1 = np.array(sizes1)
    scores = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy")
    heightShift = np.log2(sizes1) * 1
    for a in range(len(scores[0])):
        scores[:, a] = scores[:, a] + heightShift[a]
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy")
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy", scores)

def saveInitialSteps(abstract_algorithm):
    saveSequenceFileNames()
    saveInitialHMM(abstract_algorithm)
    saveScoresOriginal(abstract_algorithm)
    saveHMMsets()
    saveDecomposition()
    reallignHMMSeq()
    saveAdjustedScore()

def doPoisonRemoval(maxHMM, scores, queryNum):
    treeData = findDecomposition()
    treeSum = np.sum(treeData, axis=1)
    brotherNode = findBrotherNode()
    fatherNode = findFatherNode()
    treeDataMinusSame = np.copy(treeData)
    treeDataMinusSame[np.arange(treeData.shape[0]), np.arange(treeData.shape[0])] = 0
    nodeToCheck = treeDataMinusSame[maxHMM]
    argsCheck = np.argwhere(nodeToCheck == 1)

    queryCheck = queryNum[argsCheck[:, 0]]
    nodeCheck = argsCheck[:, 1]

    nodeFather = fatherNode[nodeCheck]
    nodeBrother = brotherNode[nodeCheck]

    scoreNode = scores[queryCheck, nodeCheck]
    scoreFather = scores[queryCheck, nodeFather]
    scoreBrother = scores[queryCheck, nodeBrother]

    ep1 = 1e-3
    diff1 = ((scoreFather + ep1) / (scoreNode + ep1)) - 1.1
    diff2 = ((scoreBrother + ep1) / (scoreFather + ep1)) - 1.1

    argsPoison = np.argwhere(np.logical_and(  diff1 > 0, diff2 > 0 ))[:, 0]

    queryPoison = queryCheck[argsPoison]
    nodePoison = nodeCheck[argsPoison]

    keepMask = treeData[maxHMM]
    keepMask[queryPoison, nodePoison] = 0

    return keepMask

def saveFromBooleanInclude(hmmLeafs, strategyName):
    dataFolderName = giveAllFileNames()[4]
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(hmmLeafs)):
        newFasta = []
        for b in range(0, len(hmmLeafs[0])):
            if hmmLeafs[a, b] == 1:
                name = get_root_temp_dir() + '/alignData/hmmSeqAlign/' + str(b) + '.fasta'
                data = loadFastaFormat(name)
                keysUsed, _ = loadFastaBasic(name)
                newFasta = newFasta + data
        saveFasta(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", newFasta)

def scoresToHMMSeq(strategyName):
    dataFolderName = giveAllFileNames()[4]
    scores = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy")
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")
    ensureFolder(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/")
    if strategyName in ['stefan_UPP', 'stefan_UPPadjusted']:
        if strategyName == 'stefan_UPP':
            scores = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/full.npy")
        scores = scores / np.max(scores)
        treeData = findDecomposition()
        treeSum = np.sum(treeData, axis=1)
        maxHMM = np.argmax(scores, axis=1).astype(int)
        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)
        HMMinverse = np.array([np.arange(HMMinverse.shape[0]), HMMinverse]).T
        print(f"saving to {get_root_temp_dir()}/data/internalData/{dataFolderName}/{strategyName}/queryToHmm/original.npy")
        np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)
        for a in range(0, len(HMMunique)):
            HMMnum = HMMunique[a]
            name = get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(a) + '.fasta'
            keys, seqs = loadFastaBasic(name)
            saveFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys, seqs)
    if strategyName in ['stefan_fastUPP', 'stefan_fastUPPexception', 'stefan_fastUPPearly']:
        if strategyName == 'stefan_fastUPPearly':
            scoresFull = np.load('%s/ensembleData/Searcher/scoreFiles/Early_score.npy' % dirName)
        else:
            scoresFull = np.load('%s/ensembleData/Searcher/scoreFiles/score.npy' % dirName)
        maxHMM = np.argmax(scoresFull, axis=1).astype(int)
        if strategyName == 'stefan_fastUPPexception':
            scores_original = np.load("%s/hmmScores/full.npy" % dirName)
            treeData = findDecomposition()
            treeSum = np.sum(treeData, axis=1)
            treeTop = np.unique(treeSum)[-3:]
            maxHMM2 = np.argmax(scores_original, axis=1).astype(int)
            argsReplace = np.argwhere(np.isin(treeSum[maxHMM], treeTop))[:, 0]
            maxHMM[argsReplace] = maxHMM2[argsReplace]
        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)
        HMMinverse = np.array([np.arange(HMMinverse.shape[0]), HMMinverse]).T
        np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)
        for a in range(0, len(HMMunique)):
            HMMnum = HMMunique[a]
            name = get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            keys, seqs = loadFastaBasic(name)
            saveFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys, seqs)

def generateNewHMM(abstract_algorithm, strategyName):
    dataFolderName = giveAllFileNames()[4]
    numHMM = int(np.max(np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")[:, 1])+1)
    for a in range(0, numHMM):
        src = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta"
        keys, seqs = loadFastaBasic(src)
        seqMatrix = []
        for seq in seqs:
            seqMatrix.append(list(seq))
        seqMatrix = np.array(seqMatrix)
        seqMatrix2 = seqMatrix.T
        usedCols = []
        for b in range(seqMatrix2.shape[0]):
            list1 = seqMatrix2[b]
            sizeReal = list1[list1 != '-'].shape[0]
            if sizeReal != 0:
                usedCols.append(b)
        usedCols = np.array(usedCols)
        usedCols_filename = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/columnSets/" + str(a) + ".npy"
        ensureFolder(usedCols_filename)
        np.save(usedCols_filename, usedCols)
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(a) + ".hmm"
        ensureFolder(hmmName)
        addHMMBuildJob(abstract_algorithm, hmmName, src)
    JobPool().wait_for_all_jobs()

def resortToUPP(strategyName, doResort=True):
    dataFolderName = giveAllFileNames()[4]
    scoreStrat_file = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/hmmScores/scoresFull/full.npy"
    ensureFolder(scoreStrat_file)
    scoreStrat = np.load(scoreStrat_file)
    scoreUPP_file = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy"
    ensureFolder(scoreUPP_file)
    scoreUPP = np.load(scoreUPP_file)
    if not doResort:
        scoreUPP[:] = -1000
    diff = np.max(scoreStrat, axis=1) - np.max(scoreUPP, axis=1)
    uppChoice = np.argmax(scoreUPP, axis=1)
    strategyChoice = np.argmax(scoreStrat, axis=1)
    argsUppBetter = np.argwhere(diff < 0)[:, 0]
    uppToUse = uppChoice[argsUppBetter]
    _, uppToUseInverse = np.unique(uppToUse, return_inverse=True)
    HMMinverse_old = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")
    max1 = np.max(HMMinverse_old[:, 1])
    HMMinverse = np.copy(strategyChoice)
    HMMinverse[argsUppBetter] = max1 + 1 + uppToUseInverse
    #TODO what about HMM which are now removed
    uppPosition = np.zeros(HMMinverse.shape[0]) - 1
    uppPosition[argsUppBetter] = uppChoice[argsUppBetter]
    HMMinverse_file = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy"
    ensureFolder(HMMinverse_file)
    np.save(HMMinverse_file, HMMinverse)
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/UPPused.npy", uppPosition)
    uppToUseUnique, index1 = np.unique(uppToUse, return_index=True)
    newNumsFull = []
    for a in range(len(uppToUseUnique)):
        uppNum = uppToUseUnique[a]
        newNum = HMMinverse[argsUppBetter[index1[a]]]
        newNumsFull.append(newNum)
        keys, seqs = loadFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(uppNum) + '.fasta')
        seqMatrix = []
        for seq in seqs:
            seqMatrix.append(list(seq))
        seqMatrix = np.array(seqMatrix)
        seqMatrix2 = seqMatrix.T
        usedCols = []
        for b in range(seqMatrix2.shape[0]):
            list1 = seqMatrix2[b]
            sizeReal = list1[list1 != '-'].shape[0]
            if sizeReal != 0:
                usedCols.append(b)
        usedCols = np.array(usedCols)
        usedCols_file = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/newHMM/columnSets/' + str(newNum) + ".npy"
        ensureFolder(usedCols_file)
        np.save(usedCols_file, usedCols)
    newNumsFull = np.array(newNumsFull)
    theoryNewNums = HMMinverse[HMMinverse > max1]
    diffIssue = np.unique(newNumsFull) - np.unique(theoryNewNums)
    diffIssue = np.sum(np.abs(diffIssue))

    assert diffIssue == 0

def stockholmToFasta(stockholmName, fastaName):
    fakeFastaName = '%s/temporaryFileSave/fakeFasta_1.fasta' % dirName
    records = SeqIO.parse(stockholmName, "stockholm")
    count = SeqIO.write(records, fakeFastaName, "fasta")
    dataAlign = loadFastaFormat(fakeFastaName)

    newData = ''.join(dataAlign)
    newData = newData.replace('\n', '')
    newData = newData.split('>')[1:]

    keys = []
    seqs = []
    for a in range(0, len(newData)):
        keys.append(newData[a][:10])
        seqs.append(newData[a][10:])

    saveFastaBasic(fastaName, keys, seqs)

def loadStockholm(stockholmName):
    data = []
    started = False
    second = False
    startedLine = 0
    secondLine = 10000000000
    file_obj = open(stockholmName)
    count1 = 0
    for line_number, i in enumerate(file_obj):
        if started:
            position = (line_number - startedLine - 1) % (secondLine - startedLine)
            if position >= len(data):
                data.append([])
            data[position].append(i)

        if (i == '\n') and ((started == True) and (second == False)):
            second = True
            secondLine = line_number

        if (i == '\n') and (started == False):
            started = True
            startedLine = line_number

    data2 = []
    keys = []

    for a in range(0, len(data) - 1):
        str1 = ''
        key =  data[a][0]
        keyList = np.array(list(key))
        argFirst = np.argwhere(keyList == ' ')[-1, 0] + 1
        key = key[:argFirst]
        keys.append(key)

        for b in range(0, len(data[a])):
            str2 = data[a][b][argFirst:-1]
            str1 += str2
        data2.append(str1)

    data2 = np.array(data2)
    keys = np.array(keys)
    return keys, data2

def loadStockholmOnlySeqs(stockholmName):
    keys, seqs = loadStockholm(stockholmName)
    seqNum = (len(seqs) - 2) // 2
    seqsNew = []
    keysNew = []
    for a in range(0, seqNum):
        seqsNew.append(seqs[a * 2])
        key = keys[a*2]
        keyList = np.array(list(key))
        argFirst = np.argwhere(keyList == ' ')[0, 0]
        key = key[:argFirst]
        keysNew.append(key)

    seqsNew = np.array(seqsNew)
    keysNew = np.array(keysNew)

    return keysNew, seqsNew

def txtToFasta(txtName, fastaName):
    file_obj = open(txtName)
    dataAlign = []
    count1 = 0
    for line_number, i in enumerate(file_obj):
        dataAlign.append(i)
    newData = ''.join(dataAlign)
    newDataList = np.array(list(newData))
    allbreaks = np.argwhere(newDataList == '\n')[0::2, 0]
    allCarrot = np.argwhere(newDataList == '>')[:, 0]
    keyLengths = (allbreaks - allCarrot - 1).astype(int)

    newData = newData.replace('\n', '')
    newData = newData.split('>')[1:]

    keys = []
    seqs = []
    for a in range(0, len(newData)):
        keys.append(newData[a][:keyLengths[a]])
        seqs.append(newData[a][keyLengths[a]:])
    keys, seqs = np.array(keys), np.array(seqs)
    saveFastaBasic(fastaName, keys, seqs)

def alignQueries(abstract_algorithm, strategyName):
    dataFolderName = giveAllFileNames()[4]
    queryName = giveQueryFileName()
    queryData = loadFastaFormat(queryName)
    queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy")
    uppHMM = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/UPPused.npy")

    predictionDataFull = []
    for b in range(0, int(np.max(queryToHmm)) + 1):
        argsHMM = np.argwhere(queryToHmm == b)[:, 0]
        if argsHMM.shape[0] > 0:
            dataNow = []
            keysNow = []
            for a in range(0, len(queryData)):
                if queryToHmm[a] == b:
                    dataNow.append(queryData[a])
            assert len(dataNow) == argsHMM.shape[0]
            queryName1 = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/inputQuery/' + str(b) + '.fasta'
            ensureFolder(queryName1)
            saveFasta(queryName1, dataNow)
            predictionName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/predictedQuery/' + str(b) + '.sto'
            ensureFolder(predictionName)
            if uppHMM[argsHMM[0]] == -1:
                hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(b) + ".hmm"
            else:
                uppNum = int(uppHMM[argsHMM[0]])
                hmmName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/initialHMM/test" + str(uppNum) + ".hmm"
            ensureFolder(hmmName)
            addHMMAlignJob(abstract_algorithm, hmmName, queryName1, predictionName)
    JobPool().wait_for_all_jobs()

def mergeAlignments(strategyName, overlapLowercase=True):
    dataFolderName = giveAllFileNames()[4]
    queryToHmm = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy")
    predictionDataFull = []
    backboneKeys, backboneSeqs = loadFastaBasic(get_root_temp_dir() + "/data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(0) + '.fasta')
    backBoneChoice = np.zeros((queryToHmm.shape[0], len(backboneSeqs[0]))).astype(str)
    backBoneChoice[:] = '-'
    backBoneChoiceBool = np.zeros(backBoneChoice.shape)

    fullInsertions = []
    for a in range(len(backboneSeqs[0])+1):
        fullInsertions.append([])
    fullInsertionsIndex = []
    for a in range(len(backboneSeqs[0])+1):
        fullInsertionsIndex.append([])

    fullInsertionsNumber = np.zeros(len(backboneSeqs[0])+1).astype(int)
    queryNames = np.zeros(queryToHmm.shape[0]).astype(str)
    
    for a in range(0, int(np.max(queryToHmm)) + 1):
        argsHMM = np.argwhere(queryToHmm == a)[:, 0]
        if argsHMM.shape[0] > 0:
            predictionName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/predictedQuery/'  + str(a) + '.sto'
            _, insertions = loadStockholm(predictionName)
            insertions = insertions[-1]
            insertions = np.array(list(insertions))
            insertions[insertions == '.'] = 1
            insertions[insertions == 'x'] = 0
            insertions = insertions.astype(int)

            keys, seqs = loadStockholmOnlySeqs(predictionName)
            queryNames[argsHMM] = keys
            seqsArray = []
            for seq in seqs:
                seqsArray.append(list(seq))
            seqsArray = np.array(seqsArray)
            queryFileName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/inputQuery/'  + str(a) + '.fasta'
            queryKey, querySeq = loadFastaBasic(queryFileName)

            assert len(keys) == len(queryKey)
            for b in range(len(queryKey)):
                assert queryKey[b] == keys[b]

            usedCols = np.load(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/newHMM/columnSets/' + str(a) + ".npy")
            matchPosition = np.argwhere(insertions == 0)[:, 0]

            try:
                assert len(usedCols) == len(matchPosition)
            except:
                print ("This issue means that the number of match states from the HMM is")
                print ("not equal to the number of non-empty columns on the sequences")
                print ("the hmm was trained on.")
                assert len(usedCols) == len(matchPosition)

            backBoneChoice[np.ix_(argsHMM, usedCols)] = seqsArray[:, matchPosition]
            backBoneChoiceBool[np.ix_(argsHMM, usedCols)] = 1

            insertionsNew = np.concatenate((np.zeros(1), insertions))
            insertionsNew = np.concatenate((insertionsNew, np.zeros(1)))
            insertionsNew = insertionsNew.astype(int)

            insertions_cumsum = np.cumsum(insertionsNew)
            matchPosition = np.argwhere(insertionsNew == 0)[:, 0]
            insertNumber = insertions_cumsum[matchPosition[1:]] - insertions_cumsum[matchPosition[:-1]]
            usedColsExtra = np.concatenate((np.zeros(1), usedCols+1)).astype(int)

            if True:
                usedColsExtra[-1] = len(backboneSeqs[0])

            for b in range(len(insertNumber)):
                if insertNumber[b] > 0:
                    toInsert = seqsArray[:, matchPosition[b]:matchPosition[b]+insertNumber[b]]
                    insertPosition = usedColsExtra[b]
                    fullInsertions[insertPosition].append(np.copy(toInsert))
                    fullInsertionsIndex[insertPosition].append(np.copy(argsHMM))
                    if overlapLowercase:
                        fullInsertionsNumber[insertPosition] = max(insertNumber[b], fullInsertionsNumber[insertPosition])
                    else:
                        fullInsertionsNumber[insertPosition] += insertNumber[b]

    #newAlignment = np.zeros((queryToHmm.shape[0], int(fullInsertionsNumber.shape[0] + np.sum(fullInsertionsNumber))  )).astype(str)
    newAlignment = np.zeros((queryToHmm.shape[0] + len(backboneKeys), int(fullInsertionsNumber.shape[0] + np.sum(fullInsertionsNumber))  )).astype(str)

    newAlignment[:] = '-'
    newAlignmentBool = np.zeros(newAlignment.shape)
    colIndex = np.cumsum(fullInsertionsNumber+1).astype(int)
    colIndex = np.concatenate((np.zeros(1), colIndex[:-1])).astype(int)
    #newAlignment[:, colIndex[1:]] = backBoneChoice #[1:] removes fake "before everything" match state.
    #newAlignmentBool[:, colIndex[1:]] = backBoneChoiceBool * 2
    newAlignment[:queryToHmm.shape[0], colIndex[1:]] = backBoneChoice #[1:] removes fake "before everything" match state.
    newAlignmentBool[:queryToHmm.shape[0], colIndex[1:]] = backBoneChoiceBool * 2

    for a in range(0, len(fullInsertions)):
        start1 = colIndex[a] + 1
        start2 = start1
        for b in range(0, len(fullInsertions[a])):
            toInsert = fullInsertions[a][b]
            insertQuery = fullInsertionsIndex[a][b]
            size1 = toInsert.shape[1]
            insertPosition = np.arange(size1) + start2
            newAlignment[np.ix_(insertQuery, insertPosition)] = np.copy(toInsert)
            newAlignmentBool[np.ix_(insertQuery, insertPosition)] = 1

            if not overlapLowercase:
                start2 += size1
    newAlignment[newAlignment == '.'] = '-'

    #newAlignmentBool
    queryNames = list(queryNames)
    for a in range(len(backboneKeys)):
        #backboneKeys, backboneSeqs
        seq1 = backboneSeqs[a]
        key1 = backboneKeys[a]

        newAlignment[a + queryToHmm.shape[0], colIndex[1:]] = np.copy(np.array(list(seq1)))
        queryNames.append(key1)
    queryNames = np.array(queryNames)

    newAlignment = newAlignment[:, 1:]
    newAlignmentBool = newAlignmentBool[:, 1:]
    colIndexTrue = (colIndex[1:] - 1).astype(int)

    newAlignmentString = []
    for a in range(newAlignment.shape[0]):
        str1 = ''.join(list(newAlignment[a]))
        newAlignmentString.append(copy.copy(str1))
    newAlignmentString = np.array(newAlignmentString)

    fastaName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignmentFasta.fasta'
    ensureFolder(fastaName)
    saveFastaBasic(fastaName, queryNames, newAlignmentString)

    fastaNameQuery = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/query_alignmentFasta.fasta'
    saveFastaBasic(fastaNameQuery, queryNames[:queryToHmm.shape[0]], newAlignmentString[:queryToHmm.shape[0]])

    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignment.npy', newAlignment)
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignmentBool.npy', newAlignmentBool)
    np.save(get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/columnIndex.npy', colIndexTrue)

def InputMergeAlignments(queryNamesFull, alignments, allInsertions, columnSets, Ncolumns, overlapLowercase=True):
    Nquery = 0
    for alignment in alignments:
        Nquery += len(alignment)

    backBoneChoice = np.zeros((Nquery, Ncolumns)).astype(str)
    backBoneChoice[:] = '-'
    backBoneChoiceBool = np.zeros(backBoneChoice.shape)

    fullInsertions = []
    for a in range(Ncolumns+1):
        fullInsertions.append([])

    fullInsertionsIndex = []
    for a in range(Ncolumns+1):
        fullInsertionsIndex.append([])

    fullInsertionsNumber = np.zeros(Ncolumns+1).astype(int)
    queryNames = np.zeros(Nquery).astype(str)
    queryCount = 0
    for a in range(0, len(alignments)):
        sizeAlign = len(alignments[a])
        argsHMM = np.arange(sizeAlign) + queryCount
        queryCount += sizeAlign
        if argsHMM.shape[0] > 0:
            insertions = allInsertions[a]
            insertions = np.array(list(insertions))
            insertions[insertions == '.'] = 1
            insertions[insertions == 'x'] = 0
            insertions = insertions.astype(int)
            seqs = alignments[a]
            queryNames[argsHMM] = queryNamesFull[argsHMM]
            seqsArray = []
            for seq in seqs:
                seqsArray.append(list(seq))
            seqsArray = np.array(seqsArray)
            usedCols = columnSets[a]
            matchPosition = np.argwhere(insertions == 0)[:, 0]
            backBoneChoice[np.ix_(argsHMM, usedCols)] = seqsArray[:, matchPosition]
            backBoneChoiceBool[np.ix_(argsHMM, usedCols)] = 1

            insertionsNew = np.concatenate((np.zeros(1), insertions))
            insertionsNew = np.concatenate((insertionsNew, np.zeros(1)))
            insertionsNew = insertionsNew.astype(int)

            insertions_cumsum = np.cumsum(insertionsNew)
            matchPosition = np.argwhere(insertionsNew == 0)[:, 0]
            insertNumber = insertions_cumsum[matchPosition[1:]] - insertions_cumsum[matchPosition[:-1]]

            usedColsExtra = np.concatenate((np.zeros(1), usedCols+1)).astype(int)

            if True:
                usedColsExtra[-1] = Ncolumns

            for b in range(len(insertNumber)):
                if insertNumber[b] > 0:
                    toInsert = seqsArray[:, matchPosition[b]:matchPosition[b]+insertNumber[b]]
                    insertPosition = usedColsExtra[b]
                    fullInsertions[insertPosition].append(np.copy(toInsert))
                    fullInsertionsIndex[insertPosition].append(np.copy(argsHMM))
                    if overlapLowercase:
                        fullInsertionsNumber[insertPosition] = max(insertNumber[b], fullInsertionsNumber[insertPosition])
                    else:
                        fullInsertionsNumber[insertPosition] += insertNumber[b]

    newAlignment = np.zeros((Nquery, int(fullInsertionsNumber.shape[0] + np.sum(fullInsertionsNumber))  )).astype(str)
    newAlignment[:] = '-'
    newAlignmentBool = np.zeros(newAlignment.shape)
    colIndex = np.cumsum(fullInsertionsNumber+1).astype(int)
    colIndex = np.concatenate((np.zeros(1), colIndex[:-1])).astype(int)

    newAlignment[:, colIndex[1:]] = backBoneChoice #[1:] removes fake "before everything" match state.
    newAlignmentBool[:, colIndex[1:]] = backBoneChoiceBool * 2

    for a in range(0, len(fullInsertions)):
        start1 = colIndex[a] + 1
        start2 = start1
        for b in range(0, len(fullInsertions[a])):
            toInsert = fullInsertions[a][b]
            insertQuery = fullInsertionsIndex[a][b]
            size1 = toInsert.shape[1]
            insertPosition = np.arange(size1) + start2
            newAlignment[np.ix_(insertQuery, insertPosition)] = np.copy(toInsert)
            newAlignmentBool[np.ix_(insertQuery, insertPosition)] = 1
            if not overlapLowercase:
                start2 += size1

    newAlignment[newAlignment == '.'] = '-'
    newAlignment = newAlignment[:, 1:]
    newAlignmentBool = newAlignmentBool[:, 1:]
    colIndexTrue = (colIndex[1:] - 1).astype(int)
    newAlignmentString = []
    for a in range(newAlignment.shape[0]):
        str1 = ''.join(list(newAlignment[a]))
        newAlignmentString.append(copy.copy(str1))
    newAlignmentString = np.array(newAlignmentString)
    return newAlignmentString

def easyInputMerge(alignments, columnSets, overlapLowercase=True):
    Ncolumns = 0
    for set1 in columnSets:
        Ncolumns = max(Ncolumns, np.max(set1))
    Ncolumns = int(Ncolumns) + 1
    Nquery = 0
    for align in alignments:
        for seq in align:
            Nquery += 1
    lowercase = np.array(list('abcdefghijklmnopqrstuvwxyz'))
    queryNamesFull = np.arange(Nquery).astype(str)
    allInsertions = []
    for a in range(len(alignments)):
        seqArray = []
        for seq in alignments[a]:
            seqArray.append(list(seq))
        seqArray = np.array(seqArray).T
        insertions1 = ''
        for b in range(len(seqArray)):
            seq1 = seqArray[b]
            if np.intersect1d(seq1, lowercase).shape[0] == 0:
                insertions1 = insertions1 + 'x'
            else:
                insertions1 = insertions1 + '.'
        allInsertions.append(insertions1)
    seqsNew = InputMergeAlignments(queryNamesFull, alignments, allInsertions, columnSets, Ncolumns, overlapLowercase=overlapLowercase)
    return seqsNew

def checkAlignmentsValid(strategyName):
    newAlignment = np.load('%s/hmmQueryList/merged/' % dirName + strategyName + '_alignment.npy')
    newAlignmentBool = np.load('%s/hmmQueryList/merged/' % dirName + strategyName + '_alignmentBool.npy')
    queryToHmm = np.load("%s/queryToHmm/" % dirName + strategyName + ".npy")
    for a in range(0, int(np.max(queryToHmm)) + 1):
        argsHMM = np.argwhere(queryToHmm == a)[:, 0]
        if argsHMM.shape[0] > 0:
            predictionName = get_root_temp_dir() + 'alignData/hmmQueryList/predictedQuery/' + strategyName + '_'  + str(a) + '.sto'
            keys, seqs = loadStockholmOnlySeqs(predictionName)
            seqsArray = []
            for b in range(0, len(seqs)):
                seq = seqs[b]
                indexInNew = argsHMM[b]
                seq2 = np.array(list(seq))
                newSeq = newAlignment[indexInNew][newAlignmentBool[indexInNew] != 0]
                newSeqBool = newAlignmentBool[indexInNew][newAlignmentBool[indexInNew] != 0]
                for c in range(0, len(seq2)):
                    if seq2[c] != newSeq[c]:
                        print ("Issue")
                        quit()

def buildAlignMerge(abstract_algorithm, strategyName, doResort=True):
    generateNewHMM(abstract_algorithm, strategyName)

    saveNewScores(abstract_algorithm, strategyName)
    resortToUPP(strategyName, doResort=doResort)

    alignQueries(abstract_algorithm, strategyName)
    mergeAlignments(strategyName)

def saveTrueUPPSubset():
    fileNames = giveAllFileNames()
    dataFolderName = giveAllFileNames()[4]
    predictionName = fileNames[3]

    predictionNameStefan = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/' + 'stefan_UPP' '/hmmQueryList/merged/alignmentFasta.fasta'
    predictionNameNew = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'

    ensureFolder(predictionNameStefan)
    ensureFolder(predictionNameNew)
    trueKey, trueSeq = loadFastaBasic(predictionName)
    predKey, predSeq = loadFastaBasic(predictionNameStefan)

    trueKey, trueSeq = trueKey[np.isin(trueKey, predKey)], trueSeq[np.isin(trueKey, predKey)]
    saveFastaBasic(predictionNameNew, trueKey, trueSeq)

def saveTrueAlignFasta():
    fileNames = giveAllFileNames()
    trueAlignment = fileNames[2]
    dataFolderName = giveAllFileNames()[4]
    trueFasta = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/trueAlignment/true_align_fragged_fasta.fasta'
    ensureFolder(trueFasta)
    txtToFasta(trueAlignment, trueFasta)

def doAllSteps(strategyName):
    buildAlignMerge(strategyName)

def scoreAlignment(strategyName):
    dataFolderName = giveAllFileNames()[4]
    saveTrueUPPSubset()
    saveTrueAlignFasta()
    trueFasta = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/trueAlignment/true_align_fragged_fasta.fasta'
    predictionName = get_root_temp_dir() + "/data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/query_alignmentFasta.fasta'
    ensureFolder(predictionName)

    UPPName = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'
    ensureFolder(UPPName)
    trueFastaNew = get_root_temp_dir() + '/data/internalData/' + dataFolderName + '/trueAlignment/subset/true_fasta.fasta'
    ensureFolder(trueFastaNew)

    queryName = giveQueryFileName()
    queryKey, _ = loadFastaBasic(queryName)
    trueKey, trueSeq = loadFastaBasic(trueFasta)
    predKey, predSeq = loadFastaBasic(predictionName)
    uppKey, uppSeq = loadFastaBasic(UPPName)

    argsIn1 = []
    for a in range(len(predKey)):
        argUse1 = np.argwhere(trueKey == predKey[a])[0, 0]
        argsIn1.append(argUse1)
    argsIn1 = np.array(argsIn1)
    trueKey, trueSeq = trueKey[argsIn1], trueSeq[argsIn1]

    saveFastaBasic(trueFastaNew, trueKey, trueSeq)

    cmd = "java -jar ./FastSP.jar -ml -r " + trueFastaNew + " -e " + predictionName

    scoreInfo  = os.popen(cmd).readlines()

    scoreInfo2 = []
    for a in range(0, len(scoreInfo)):
        string1 = scoreInfo[a].split(' ')[1].replace('\n', '')
        scoreInfo2.append(float(string1))

    print (strategyName)
    print (scoreInfo)

def compareToUPP():
    keys1, seqs1 = loadFastaBasic(predictionStefan)
    keys2, seqs2 = loadFastaBasic(predictionUPP)

    keys2, seqs2 = keys2[np.isin(keys2, keys1)], seqs2[np.isin(keys2, keys1)]

    print (keys1[0])
    print (keys2[0])
    print (seqs1[0])
    print ('')
    print (seqs2[0])


def hierchySearch(abstract_algorithm, fakeSimulate=False):

    earlyStop = False
    observeNephew = False
    dataFolderName = giveAllFileNames()[4]
    # np.load("./data/internalData/" + dataFolderName + "/hmmScores/full.npy")
    scores_original = np.load(get_root_temp_dir() + "/data/internalData/%s/hmmScores/full.npy" % (dataFolderName))

    hmmNames = []
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        inputFiles = giveAllFileNames()
        dataFileName = inputFiles[4]
        hmmName = get_root_temp_dir() + "/data/internalData/" + dataFileName + "/initialHMM/test" + str(a) + ".hmm"
        hmmNames.append(hmmName)

    queryName = giveQueryFileName()
    queryData = loadFastaFormat(queryName)
    Nquery = int(len(queryData))
    Nhmm = int(len(hmmNames))

    edges = findTreeEdges()
    children = np.argsort(edges, axis=1)[:, -2:]
    childNum = np.sum(edges, axis=1)
    brotherNode = findBrotherNode()

    treeData = findDecomposition()
    treeSum = np.sum(treeData, axis=1)
    sortedSize = np.argsort(treeSum)
    scoresFull = np.zeros((Nquery, Nhmm)) - 1000
    scoreAdjustment = np.zeros(scoresFull.shape)
    initialSearchSet = np.array([0])

    queryPart = np.arange(Nquery * len(initialSearchSet)) % Nquery
    hmmPart = np.array(initialSearchSet).repeat(Nquery)

    queryHMM = np.array([queryPart, hmmPart]).T
    scoreName = ''
    if fakeSimulate:
        scores = np.copy(scores_original)
    else:
        scores = saveScoreFromBool(abstract_algorithm, hmmNames, queryHMM, scoreName, noSave=True)
    scoresFull[queryHMM[:, 0], queryHMM[:, 1]] = scores[queryHMM[:, 0], queryHMM[:, 1]]

    queryHMM_original = np.copy(queryHMM)
    maxHMM = np.argmax(scoresFull, axis=1)
    boolDone = np.zeros((Nquery, Nhmm))
    boolDone[:, 0] = 1
    maxHMMlist = []
    scoresRecord = []

    a = 1
    done = False
    while not done:
        maxNotLeaf = np.argwhere(childNum[maxHMM] != 0)[:, 0]
        if observeNephew:
            maxHMMChild_original = children[maxHMM[maxNotLeaf]]
            maxHMMChild_original = np.concatenate((maxHMMChild_original, children[brotherNode[maxHMM[maxNotLeaf]]]), axis=1)
        else:
            maxHMMChild_original = children[maxHMM[maxNotLeaf]]

        maxHMMChild = np.copy(maxHMMChild_original).reshape((maxHMMChild_original.size,))
        queryCorrespond = maxNotLeaf[np.arange(maxHMMChild.size) // (maxHMMChild_original.shape[1]) ]
        queryToHMM_try = np.array([queryCorrespond, maxHMMChild]).T
        queryHMM = queryToHMM_try[boolDone[queryToHMM_try[:, 0], queryToHMM_try[:, 1] ] == 0]
        boolDone[queryHMM[:, 0], queryHMM[:, 1] ] = 1

        if queryHMM.shape[0] != 0:
            scoreName = ''
            if fakeSimulate:
                scores = np.copy(scores_original)
            else:
                scores = saveScoreFromBool(abstract_algorithm, hmmNames, queryHMM, scoreName, noSave=True)

            scoresFull[queryHMM[:, 0], queryHMM[:, 1]] = np.copy(scores[queryHMM[:, 0], queryHMM[:, 1]])
            newScores1 = scoresFull[maxNotLeaf, maxHMMChild_original[:, 0]]
            newScores2 = scoresFull[maxNotLeaf, maxHMMChild_original[:, 1]]
            newScoresMax = np.argmax(np.array([newScores1, newScores2]), axis=0)

            newChild = maxHMMChild_original[np.arange(newScores1.shape[0]), newScoresMax]
            if earlyStop:
                maxHMM = np.argmax(scoresFull, axis=1)
            else:
                maxHMM[maxNotLeaf] = newChild

            scoreStep = scoresFull[np.arange(maxHMM.shape[0]), maxHMM]
            scoresRecord.append(np.copy(scoreStep))
            maxHMMlist.append(np.copy(maxHMM[:10]))
        else:
            done = True
        a += 1

    scoresRecord = np.array(scoresRecord)

    np.save('%s/ensembleData/Searcher/scoreFiles/score.npy' % dirName, scoresFull)
    np.save('%s/ensembleData/Searcher/scoreFiles/Bool.npy' % dirName, boolDone)

    maxHMM_guess = np.argmax(scoresFull, axis=1)
    maxHMM_original = np.argmax(scores_original, axis=1)

    argDiff = np.argwhere((maxHMM_original-maxHMM_guess) != 0)[:, 0]
    argDiffTop = argDiff[treeSum[maxHMM_guess[argDiff]] == 999]
    print (np.mean(scores_original[argDiff, maxHMM_original[argDiff]] - scores_original[argDiff, maxHMM_guess[argDiff]]))
    print (argDiff.shape)
    print (np.unique(treeSum[maxHMM_guess[argDiff]], return_counts=True)[1][-3:])
    print (np.unique(treeSum[maxHMM_guess], return_counts=True)[1][-3:])
    bitscore_guess = scores_original[np.arange(Nquery), maxHMM_guess]
    bitscore_original = scores_original[np.arange(Nquery), maxHMM_original]

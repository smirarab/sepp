import numpy as np
import matplotlib.pyplot as plt
import copy
import os

#os.system("hmmbuild ./stefanData/globins4.hmm ./stefanData/globins4.sto")
hmmSeqFile = ''
queryName = ''
trueAlignment = ''
predictionName = ''

def ensureFolder(fileName):

    def removeEnd(fileName):
        fileName = '/'.join(fileName.split('/')[:-1])
        return (fileName)


    fileName = removeEnd(fileName)

    if os.path.exists(fileName) == False:

        toMake = []
        while os.path.exists(fileName) == False:
            toMake.append(fileName)
            fileName = removeEnd(fileName)

        toMake = np.array(toMake)[-1::-1]

        for a in range(len(toMake)):
            os.mkdir(toMake[a])

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

def setAllFileNames(givenHSF, givenQN, givenTA, givenPN):
    global hmmSeqFile 
    hmmSeqFile = givenHSF
    global queryName 
    queryName = givenQN
    global trueAlignment
    trueAlignment = givenTA
    global predictionName
    predictionName = givenPN

def giveAllFileNames():
    
    global hmmSeqFile
    global queryName
    global trueAlignment
    global predictionName

    if False:
        hmmSeqFile = './alignData/tmpfiles/R13/output.4za7fub8/root/P_0/'
        queryName = './alignData/tmpfiles/R13/output.4za7fub8/fragment_chunks/fragment_chunk_0ndvaj14s.fasta'
        trueAlignment = './alignData/UnalignFragTree/low_frag/1000M1/R13/true_align_fragged.txt'
        predictionName = './alignData/UPPoutput/R13_output_alignment.fasta'

    if False:
        hmmSeqFile = './alignData/tmpfiles/R14/output.1bkwwmuv/root/P_0/'
        queryName = './alignData/tmpfiles/R14/output.1bkwwmuv/fragment_chunks/fragment_chunk_0451q1ctg.fasta'
        trueAlignment = './alignData/UnalignFragTree/high_frag/1000M1/R14/true_align_fragged.txt'
        predictionName = './alignData/UPPoutput/R14_output_alignment.fasta'

    if False:
        hmmSeqFile = './alignData/tmpfiles/R14_decomp5/output.gmue63ev/root/P_0/'
        queryName = './alignData/tmpfiles/R14_decomp5/output.gmue63ev/fragment_chunks/fragment_chunk_00twg8q8d.fasta'
        trueAlignment = './alignData/UnalignFragTree/high_frag/1000M1/R14/true_align_fragged.txt'
        predictionName = './alignData/UPPoutput/R14_output_alignment.fasta'

    if False:
        hmmSeqFile = './alignData/tmpfiles/R14_decomp1/output.vlkvwo6g/root/P_0/'
        queryName = './alignData/tmpfiles/R14_decomp1/output.vlkvwo6g/fragment_chunks/fragment_chunk_0vxzt57yz.fasta'
        trueAlignment = './alignData/UnalignFragTree/high_frag/1000M1/R14/true_align_fragged.txt'
        predictionName = './alignData/UPPoutput/R14_output_alignment.fasta'

    if False: # test query sequences
        hmmSeqFile = './alignData/tmpfiles/M4/R0/Rep0_rtol_1e-1max/output.t5g238nm/root/P_0/'
        queryName = './alignData/tmpfiles/M4/R0/Rep0_rtol_1e-1max/output.t5g238nm/fragment_chunks/fragment_chunk_049mh1i47.fasta'
        #trueAlignment = './alignData/UnalignFragTree/high_frag/1000M1/R0/true_align_fragged.txt'
        trueAlignment = "./alignData/subsetTrueAln/M4R0_truealign_frag.txt"
        predictionName = './alignData/UPPoutput/M4Rep0_rtol_1e-1max_output_alignment.fasta'

    if False: # test subset decomposition
        hmmSeqFile = './alignData/tmpfiles/M4/R0/30/output.q22vqcip/root/P_0/'
        queryName = './alignData/tmpfiles/M4/R0/30/output.q22vqcip/fragment_chunks/fragment_chunk_0p__pj929.fasta'
        trueAlignment = "./alignData/UnalignFragTree/high_frag/1000M4/R0/true_align_fragged.txt"
        predictionName = './alignData/UPPoutput/M4R0decomp30_output_alignment.fasta'

    fileNames = [hmmSeqFile, queryName, trueAlignment, predictionName]

    return fileNames




def giveQueries():
    #file_obj = open('./alignData/unaligned.txt')
    queryName = giveQueryFileName()
    file_obj = open(queryName)
    #file_obj = open('./alignData/tmpfiles/output.ped7oo42/backbone/queryqqpulswm.fas')
    queries = []
    count1 = 0
    for line_number, i in enumerate(file_obj):
        #print ([i])
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

    #print (seqArray2[0])
    #quit()
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
    with open(name, 'w') as f:
        f.write('\n'.join(formatedData))

def loadFastaFormat(name):
    file_obj = open(name)
    data = []
    name = ''
    for line_number, i in enumerate(file_obj):


        #if line_number < 2:
        #if line_number >= 1498:
        #    print ([i])
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
        #print ([i])
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
    #src0 = './alignData/tmpfiles/output.ped7oo42/root/P_0/'
    #src0 = './alignData/tmpfiles/R13/output.4za7fub8/root/P_0/'
    #src0 = './alignData/tmpfiles/R14_high/output.zn48i3r_/root/P_0/'
    #src0 = './alignData/tmpfiles/R14/output.zn48i3r_/root/P_0/'
    #src0 = './alignData/tmpfiles/R14/output.1bkwwmuv/root/P_0/'

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

    sequenceFileNames_saveName = './data/internalData/' + dataFileName + '/seqFileNames/names.npy'
    ensureFolder(sequenceFileNames_saveName)
    np.save(sequenceFileNames_saveName, sequenceFileNames)

#saveSequenceFileNames()
#quit()




def giveSequenceFileNames():
    '''
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
    sequenceFileNames = sequenceFileNames[np.argsort(Anums)]
    '''

    inputFiles = giveAllFileNames()
    dataFileName = inputFiles[4]

    sequenceFileNames_saveName = './data/internalData/' + dataFileName + '/seqFileNames/names.npy'

    sequenceFileNames = np.load(sequenceFileNames_saveName)


    return sequenceFileNames

def giveQueryFileName():

    #queryName = './alignData/tmpfiles/output.ped7oo42/backbone/queryqqpulswm.fas'
    #queryName = './alignData/tmpfiles/R13/output.4za7fub8/fragment_chunks/fragment_chunk_0ndvaj14s.fasta'
    #queryName = './alignData/tmpfiles/R14/output.zn48i3r_/fragment_chunks/fragment_chunk_0s4gnzx9j.fasta'
    #queryName = './alignData/tmpfiles/R14/output.1bkwwmuv/fragment_chunks/fragment_chunk_0451q1ctg.fasta'

    fileNames = giveAllFileNames()
    queryName = fileNames[1]

    return queryName


def saveInitialHMM(minSize=1):
    sequenceFileNames = giveSequenceFileNames()

    for a in range(0, len(sequenceFileNames)):
        src = sequenceFileNames[a]
        #os.system("hmmbuild ./alignData/initialHMM/test" + str(a) + ".hmm " + src)
        #os.system("./hmmer-3.0/src/hmmbuild  --symfrac 0.0 --informat afa  ./alignData/initialHMM/test" + str(a) + ".hmm " + src)

        #hmmName = "./alignData/initialHMM/test" + str(a) + ".hmm"

        inputFiles = giveAllFileNames()
        dataFileName = inputFiles[4]

        hmmName = "./data/internalData/" + dataFileName + "/initialHMM/test" + str(a) + ".hmm"
        ensureFolder(hmmName)
        runHMMbuild(hmmName, src)


#saveInitialHMM(minSize=1)
#quit()


def giveHMMversion():

    #hmmVersion = './hmmer-3.0/src/'
    hmmVersion = ''

    return hmmVersion



def saveScore(hmmName, queryName, scoreName):

    ensureFolder(scoreName)

    hmmVersion = giveHMMversion()
    os.system(hmmVersion + "hmmsearch --noali --tblout " + scoreName + " --cpu 1 -E 99999999 --max " + hmmName + " " + queryName)

def runHMMbuild(hmmName, seqName):

    ensureFolder(hmmName)

    hmmVersion = giveHMMversion()
    os.system(hmmVersion + "hmmbuild  --symfrac 0.0 --informat afa  " + hmmName + " " + seqName)




def runHMMalign(hmmName, queryName, predictionName):

    ensureFolder(predictionName)

    hmmVersion = giveHMMversion()
    if hmmVersion == '':
        os.system(hmmVersion + "hmmalign --dna -o " + predictionName + " " + hmmName + ' ' + queryName + " ")
    else:
        os.system(hmmVersion + "hmmalign --allcol  --dna -o " + predictionName + " " + hmmName + ' ' + queryName + " ")






def saveScoreSimple(hmmNames, queryNames, scoreName):

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
        #scoreNameTemp = "./alignData/hmmScores/temporaryStorage/temp1.txt"

        scoreNameTemp = "./data/temporaryStorage/temp_" + randInt + ".txt"
        queryName = queryNames[a]


        ensureFolder(scoreNameTemp)
        saveScore(hmmName, queryName, scoreNameTemp)

        #print ('HI')

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
                #print ([i[argSpace[4]:argSpace[5]]])
                #bitScore = float(i[argSpace[4]:argSpace[5]].replace(' ', ''))

                bitScore = float(i[posScore-1:posScore+5].replace(' ', ''))




                dataOld.append([sequenceName, bitScore])


                #print ('bitScore ', bitScore)

            #quit()

        dataOld = np.array(dataOld)

        for b in range(0, len(dataOld)):
            queryN = queryDict[dataOld[b, 0]]
            #print (queryN, a, b)
            data[queryN, a] = float(dataOld[b, 1])

    ensureFolder(scoreName)
    np.save(scoreName, data)

    os.remove(scoreNameTemp)



def saveNewScores(strategyName):

    dataFolderName = giveAllFileNames()[4]

    queryNameO = giveQueryFileName()
    queryData = loadFastaFormat(queryNameO)

    queryNames = []
    hmmNames = []


    queryToHmm = np.load("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")

    numHMM = int(np.max(queryToHmm[:, 1])+1)
    for b in range(0, numHMM):


        #argsHMM = np.argwhere(queryToHmm == b)[:, 0]
        argsHMM = queryToHmm[:, 0][queryToHmm[:, 1] == b]
        dataNow = []
        keysNow = []
        #for a in range(0, len(queryData)):
        #    if queryToHmm[a] == b:
        #        dataNow.append(queryData[a])
        for a in argsHMM:
             dataNow.append(queryData[a])
        #queryName = './alignData/hmmQueryList/inputQuery/' + strategyName + '_'  + str(b) + '.fasta'

        dataFolderName = giveAllFileNames()[4]


        queryName = "./data/internalData/" + dataFolderName + "/" + strategyName + '/inputQuery/'  + str(b) + '.fasta'

        ensureFolder(queryName)
        saveFasta(queryName, dataNow)

        hmmName = "./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(b) + ".hmm "
        ensureFolder(hmmName)
        #queryNames.append(queryNameO)
        queryNames.append(queryName)
        hmmNames.append(hmmName)


    dataFolderName = giveAllFileNames()[4]

    scoreName = "./data/internalData/" + dataFolderName + "/" + strategyName + "/hmmScores/scoresFull/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(hmmNames, queryNames, scoreName)




def saveScoreFromBool(hmmNames, queryHMM, scoreName, noSave=False):

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



            scoreNameTemp = "./data/temporaryStorage/temp_" + randInt1 + ".txt"
            queryName = giveQueryFileName()
            queryData = loadFastaFormat(queryName)

            queryNameTemp = './data/temporaryStorage/temp_' + randInt2 + '.fasta'

            #argsInclude = np.argwhere(queryBool[a] == 1)[:, 0]

            queryNow = []
            for b in argsInclude:
                queryNow.append(queryData[b])

            ensureFolder(queryNameTemp)
            saveFasta(queryNameTemp, queryNow)

            ensureFolder(scoreNameTemp)
            saveScore(hmmName, queryNameTemp, scoreNameTemp)

            #print ('HI')

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
                    #print (i)
                    sequenceName = i[:20].replace(' ', '')
                    #print ([i[argSpace[4]:argSpace[5]]])
                    #bitScore = float(i[argSpace[4]:argSpace[5]].replace(' ', ''))

                    bitScore = float(i[posScore-1:posScore+5].replace(' ', ''))



                    dataOld.append([sequenceName, bitScore])

                    #if 0 in argsInclude:
                    #    print ('bitScore ', bitScore)
                    #    quit()

                #quit()

            dataOld = np.array(dataOld)

            for b in range(0, len(dataOld)):
                queryN = queryDict[dataOld[b, 0]]
                #print (queryN, a, b)
                data[queryN, a] = float(dataOld[b, 1])

    #print (scoreName)

    #print (np.max(data[0]))

    os.remove(scoreNameTemp)
    os.remove(queryNameTemp)

    if noSave:
        return data
    else:
        ensureFolder(scoreName)
        np.save(scoreName, data)


#saveNewScores('stefan_removeLeafExtra')
#saveNewScores('stefan_addUnpoison')
#saveNewScores('stefan_removePoisonGrandpa')
#quit()


def compareScores(strategyName):

    scoreStrat = np.load("./alignData/hmmScores/scoresFull/" + strategyName + "_full.npy")
    #scoreStrat = np.load("./alignData/hmmScores/fullReAlign.npy")
    #scoreStrat = np.load("./alignData/hmmScores/fullRoot.npy")
    #scoreStrat = np.load("./alignData/hmmScores/full.npy")
    scoreUPP = np.load("./alignData/hmmScores/full.npy")
    #scoreUPP = np.load("./alignData/hmmScores/testStandardAlign.npy")
    #scoreUPP = np.load("./alignData/hmmScores/fullReAlign.npy")


    #HMMinverse = np.load("./alignData/queryToHmm/" + strategyName + ".npy")
    #HMMunique = np.load("./alignData/queryToHmm/Unique_" + strategyName + ".npy")


    #img = np.zeros(scoreStrat.shape)
    #img[np.arange(scoreUPP.shape[0]), np.argmax(scoreUPP, axis=1)] = 1
    #img[np.arange(scoreUPP.shape[0]), np.argmax(scoreStrat, axis=1)] = 1

    #argBad = np.argmax(np.abs(np.max(scoreStrat, axis=1) - np.max(scoreUPP, axis=1) ))

    #HMMunique = np.unique(HMMinverse).astype(int)


    #scoreStrat2 = np.zeros(scoreStrat.shape)
    #scoreStrat2[:, HMMunique] = scoreStrat[:, :HMMunique.shape[0]]

    #plt.plot(np.max(scoreStrat, axis=1))
    #plt.plot(np.max(scoreUPP, axis=1))
    plt.plot(np.max(scoreStrat, axis=1)-np.max(scoreUPP, axis=1))
    plt.show()



#compareScores('stefan_removePoison')
#compareScores('stefan_removePoisonLeaf')
#compareScores('stefan_removePoisonInternal')
#compareScores('stefan_addUnpoison')
#compareScores('stefan_removePoisonGrandpa')
#compareScores('stefan_brotherMixture')
#quit()
########compareScores('stefan_removeLeaf')



def saveScoresOriginal():
    sequenceFileNames = giveSequenceFileNames()

    hmmNames = []
    queryNames = []

    for a in range(0, len(sequenceFileNames)):
        #hmmName = "./alignData/initialHMM/test" + str(a) + ".hmm"

        dataFolderName = giveAllFileNames()[4]

        hmmName = "./data/internalData/" + dataFolderName + "/initialHMM/test" + str(a) + ".hmm"

        hmmNames.append(hmmName)
        queryName = giveQueryFileName()
        queryNames.append(queryName)


    dataFolderName = giveAllFileNames()[4]
    scoreName = "./data/internalData/" + dataFolderName + "/hmmScores/full.npy"
    ensureFolder(scoreName)
    saveScoreSimple(hmmNames, queryNames, scoreName)

def processScores():
    sequenceFileNames = giveSequenceFileNames()

    hmmNames = []
    queryNames = []

    for a in range(0, len(sequenceFileNames)):
        dataFolderName = giveAllFileNames()[4]

        hmmName = "./data/internalData/" + dataFolderName + "/initialHMM/test" + str(a) + ".hmm"

        #hmmName = "./alignData/initialHMM/test" + str(a) + ".hmm"
        #hmmName = "./alignData/initialHMM/testRoot" + str(a) + ".hmm"
        #hmmName = "./alignData/initialHMM/testReAlign" + str(a) + ".hmm"
        #hmmName = "./alignData/initialHMM/testStandardAlign" + str(a) + ".hmm"

        hmmNames.append(hmmName)
        queryName = giveQueryFileName()
        queryNames.append(queryName)

    dataFolderName = giveAllFileNames()[4]
    scoreName = "./data/internalData/" + dataFolderName + "/hmmScores/full.npy"
    #scoreName = "./alignData/hmmScores/full.npy"
    #scoreName = "./alignData/hmmScores/fullRoot.npy"
    #scoreName = "./alignData/hmmScores/fullReAlign.npy"
    #scoreName = "./alignData/hmmScores/testStandardAlign.npy"
    ensureFolder(scoreName)
    saveScoreSimple(hmmNames, queryNames, scoreName)


#processScores()
#quit()

def processScoresOld():
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        #data = np.loadtxt("./alignData/hmmScores/" + str(a) + ".txt", dtype=str)
        #a = 51
        data = []
        file_obj = open("./alignData/hmmScores/raw/" + str(a) + ".txt")
        count1 = 0
        for line_number, i in enumerate(file_obj):
            #print (line_number)
            #print ([i])

            ar = i.split('  ')

            #print (i)
            #if line_number > 20:
            #    quit()

            if (len(ar) >= 30) and (line_number >= 3):
                bitScore = float(i[76:83].replace(' ', ''))
                sequenceName = i[:20].replace(' ', '')

                #print (bitScore)
                #quit()

                #print ([sequenceName, bitScore])
                data.append([sequenceName, bitScore])

        data = np.array(data)
        np.save("./alignData/hmmScores/processed/" + str(a) + ".npy", data)


#processScores()
#quit()

#queries = giveQueries()
#print (len(queries))
#quit()

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

        dataOld = np.load("./data/internalData/" + dataFolderName + "/hmmScores/processed/" + str(a) + ".npy")
        for b in range(0, len(dataOld)):
            #dataOld[b, 0]
            #queryDict[dataOld[b, 0]]
            queryN = queryDict[dataOld[b, 0]]
            data[queryN, a] = float(dataOld[b, 1])

    dataFolderName = giveAllFileNames()[4]
    ensureFolder("./data/internalData/" + dataFolderName + "/hmmScores/full.npy")
    np.save("./data/internalData/" + dataFolderName + "/hmmScores/full.npy", data)


#saveScoresBySeq()
#quit()

def saveDecomposition():
    dataFolderName = giveAllFileNames()[4]

    sets = np.load("./data/internalData/" + dataFolderName + "/hmmSets.npy", allow_pickle=True)

    treeData = np.zeros((len(sets), len(sets)))
    for a in range(0, len(sets)):
        for b in range(0, len(sets)):
            set1, set2 = np.array(sets[a]), np.array(sets[b])
            if np.intersect1d(set1, set2).shape[0] == set2.shape[0]:
                treeData[a, b] = 1

    ensureFolder("./data/internalData/" + dataFolderName + "/treeDecomp.npy")
    np.save("./data/internalData/" + dataFolderName + "/treeDecomp.npy", treeData)

#saveDecomposition()
#quit()

def saveHMMsets():
    sequenceFileNames = giveSequenceFileNames()
    sets = []

    for a in range(0, len(sequenceFileNames)):
        src = sequenceFileNames[a]

        #print (src)

        #print ("F")
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

    ensureFolder("./data/internalData/" + dataFolderName + "/hmmSets.npy")
    np.save("./data/internalData/" + dataFolderName + "/hmmSets.npy", sets)

#saveHMMsets()
#quit()


def findDecomposition():

    dataFolderName = giveAllFileNames()[4]
    treeData = np.load("./data/internalData/" + dataFolderName + "/treeDecomp.npy")

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

        #treeSumIn = np.zeros(treeSum.shape) + 10000000 #treeSumIn things containing argmax
        #treeSumIn[nodesArray] = treeSum[nodesArray]
        #treeSumIn = np.copy(treeSum)
        treeSumIn = np.copy(treeSum) + ((1-treeData[:, argMax]) * 1000000)
        argMin = int(np.argmin(treeSumIn))

        if treeData[argMin, argMax] == 1:
            edges[argMin, argMax] = 1
        #edges[argMax, argMin] = 1



    #plt.imshow(edges)
    #plt.show()

    return edges

#findTreeEdges()
#quit()

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

    #print (data.shape)

    data4 = data[arange3, arange2]

    #print (data4.shape)

    return data4






def reallignHMMSeq():

    sequenceFileNames = giveSequenceFileNames()

    rootName = sequenceFileNames[0]

    rootKeys, rootSeqs = loadFastaBasic(rootName)
    #quit()

    for a in range(0, len(sequenceFileNames)):

        name = sequenceFileNames[a]
        nowKeys, nowSeqs = loadFastaBasic(name)

        argsIn = np.argwhere(np.isin(rootKeys, nowKeys))[:, 0]
        nowKeys, nowSeqs = rootKeys[argsIn], rootSeqs[argsIn]

        #nowSeqs = removeEmptyColumns(nowSeqs)

        dataFolderName = giveAllFileNames()[4]

        #if False:
        ensureFolder("./data/internalData/" + dataFolderName + "/hmmSeqAlign/" + str(a) + '.fasta')
        saveFastaBasic("./data/internalData/" + dataFolderName + "/hmmSeqAlign/" + str(a) + '.fasta', nowKeys, nowSeqs)


#reallignHMMSeq()
#quit()



def compareHMM():

    sequenceFileNames = giveSequenceFileNames()

    for a in range(0, len(sequenceFileNames)):

        HMM1 = "./alignData/initialHMM/testStandardAlign" + str(a) + ".hmm "
        HMM2 = "./alignData/initialHMM/test" + str(a) + ".hmm "

        #print (nameTemp, src)
        cmd = 'diff ' + HMM1 + ' ' + HMM2

        print (cmd)
        quit()

        from subprocess import run
        scoreInfo  = os.popen(cmd).readlines()

        print (scoreInfo)

#compareHMM()
#quit()

def saveAdjustedScore():

    dataFolderName = giveAllFileNames()[4]
    sets = np.load("./data/internalData/" + dataFolderName + "/hmmSets.npy", allow_pickle=True)

    sizes1 = []
    for a in range(len(sets)):
        sizes1.append(len(sets[a]))

    sizes1 = np.array(sizes1)


    scores = np.load("./data/internalData/" + dataFolderName + "/hmmScores/full.npy")

    '''
    #scores = scores #/ np.max(scores)

    treeData = findDecomposition()

    treeSum = np.sum(treeData, axis=1)

    #print (np.unique(treeSum, return_counts=True))

    argsInternal = np.argwhere(treeSum > 2)[:, 0]
    treeDataLeaf = np.copy(treeData)
    treeDataLeaf[:, argsInternal] = 0
    treeSumLeaf = np.sum(treeDataLeaf, axis=1)

    #print (np.max(scores))
    #quit()

    #plt.plot(treeSumLeaf)
    #plt.show()
    #quit()
    '''

    #heightShift = np.log2(sizes1) * 0.01 * np.max(scores)
    #heightShift = np.log2(sizes1) * 0.005 * np.max(scores)
    heightShift = np.log2(sizes1) * 1

    for a in range(len(scores[0])):
        scores[:, a] = scores[:, a] + heightShift[a]

    ensureFolder("./data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy")
    np.save("./data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy", scores)

#saveAdjustedScore()
#quit()


def saveInitialSteps():
    saveSequenceFileNames()
    saveInitialHMM()
    saveScoresOriginal()
    saveHMMsets()
    saveDecomposition()
    reallignHMMSeq()
    saveAdjustedScore()

#saveInitialSteps()
#quit()


def doPoisonRemoval(maxHMM, scores, queryNum):


    #print ("A")
    treeData = findDecomposition()
    treeSum = np.sum(treeData, axis=1)
    brotherNode = findBrotherNode()
    fatherNode = findFatherNode()

    treeDataMinusSame = np.copy(treeData)
    treeDataMinusSame[np.arange(treeData.shape[0]), np.arange(treeData.shape[0])] = 0

    nodeToCheck = treeDataMinusSame[maxHMM]
    #print (nodeToCheck[0])
    argsCheck = np.argwhere(nodeToCheck == 1)

    queryCheck = queryNum[argsCheck[:, 0]]
    nodeCheck = argsCheck[:, 1]

    #print (queryCheck[:5])
    #print (nodeCheck[:5])
    #quit()

    nodeFather = fatherNode[nodeCheck]
    nodeBrother = brotherNode[nodeCheck]

    scoreNode = scores[queryCheck, nodeCheck]
    scoreFather = scores[queryCheck, nodeFather]
    scoreBrother = scores[queryCheck, nodeBrother]

    #print (scoreNode[:5])
    #print (scoreFather[:5])
    #print (scoreBrother[:5])

    #diff1 = scoreFather - scoreNode
    #diff2 = scoreBrother - scoreFather
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

    #keysEverUsed = np.array([])
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(hmmLeafs)):
        newFasta = []
        for b in range(0, len(hmmLeafs[0])):
            if hmmLeafs[a, b] == 1:
                #name = sequenceFileNames[b]
                name = './alignData/hmmSeqAlign/' + str(b) + '.fasta'
                data = loadFastaFormat(name)
                keysUsed, _ = loadFastaBasic(name)
                #keysEverUsed = np.concatenate((keysEverUsed, keysUsed))
                newFasta = newFasta + data


        #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", newFasta)
        saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", newFasta)





def scoresToHMMSeq(strategyName):

    #scores is a matrix of bit scores where the 0 axis is query, and the 1 axis is hmm
    #scores = np.load("./alignData/hmmScores/full.npy")

    dataFolderName = giveAllFileNames()[4]


    scores = np.load("./data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy")


    ensureFolder("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")
    ensureFolder("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/")



    if strategyName in ['stefan_basic', 'stefan_removeLeaf', 'stefan_onlyOneAbove', 'stefan_useAboveLeaf']:

        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)


        #unique1 = np.unique(treeSum)

        oneMask = np.zeros(scores.shape)
        oneMask[:, treeSum==1] = 1
        #oneMask[:, treeSum<=2] = 1
        maxHMM = np.argmax(scores, axis=1)

        if strategyName == 'stefan_useAboveLeaf':
            oneMask[:, treeSum==3] = 1

        usedMask = treeData[maxHMM, :]

        usedOneMask = usedMask * oneMask

        meanUsed = np.sum(scores * usedOneMask, axis=1) / np.sum(usedOneMask, axis=1)
        meanUsed2 = meanUsed.repeat(scores.shape[1]).reshape(scores.shape)


        aboveMask = np.zeros(scores.shape)
        aboveMask[((scores + 0.01)/(meanUsed2 + 0.01)) > (2/3)] = 1

        aboveMaskOne = aboveMask * oneMask

        if strategyName == 'stefan_removeLeaf':
            aboveMaskOne = aboveMaskOne * usedMask

        if strategyName == 'stefan_onlyOneAbove':

            containsMax = treeData[:, maxHMM].T
            containsMax[np.arange(maxHMM.shape[0]), maxHMM] = 0
            treeSumRepeat = treeSum.repeat(maxHMM.shape[0]).reshape((treeSum.shape[0], maxHMM.shape[0])).T
            treeSumRepeat[containsMax == 0] = 1000000
            directAbove = np.argmin(treeSumRepeat, axis=1).astype(int)
            directAboveLeafs = treeData[directAbove]
            aboveMaskOne = directAboveLeafs * usedMask


        #aboveScoreInterweiv = interweive(np.array([scores, usedOneMask, aboveMaskOne]))
        #plt.imshow(aboveScoreInterweiv[60:100, treeSum==1])
        ##plt.imshow(aboveScoreInterweiv[0:30])
        #plt.show()
        #quit()

        aboveMaskOneString = []
        for a in range(0, len(aboveMaskOne)):
            string1 = ''
            for b in range(0, len(aboveMaskOne[a])):
                string1 = string1 + str(int(aboveMaskOne[a, b])) + ':'
            aboveMaskOneString.append(string1)
        aboveMaskOneString = np.array(aboveMaskOneString)

        _, hmmUniqueIndex, queryToHmm = np.unique(aboveMaskOneString, return_index=True, return_inverse=True)

        hmmLeafs = aboveMaskOne[hmmUniqueIndex]

        queryToHmm = np.array([np.arange(queryToHmm.shape[0]), queryToHmm]).T



        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", queryToHmm)

        #keysEverUsed = np.array([])
        sequenceFileNames = giveSequenceFileNames()
        for a in range(0, len(hmmLeafs)):
            newFasta = []
            for b in range(0, len(hmmLeafs[0])):
                if hmmLeafs[a, b] == 1:
                    #name = sequenceFileNames[b]
                    #name = './alignData/hmmSeqAlign/' + str(b) + '.fasta'
                    name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(b) + '.fasta'
                    data = loadFastaFormat(name)
                    keysUsed, _ = loadFastaBasic(name)
                    #keysEverUsed = np.concatenate((keysEverUsed, keysUsed))
                    newFasta = newFasta + data


            #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", newFasta)
            saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", newFasta)

        #keysEverUsed = np.unique(keysEverUsed)
        #keysAll, seqsAll = loadFastaBasic('./alignData/hmmSeqAlign/' + str(0) + '.fasta')
        #seqNotUsed = seqsAll[np.isin(keysAll, keysEverUsed) == False]
        #keyNotUsed = keysAll[np.isin(keysAll, keysEverUsed) == False]
        #if keyNotUsed.shape[0] != 0:
        #    saveFastaBasic('./alignData/newHMM/unusedSeq/' + strategyName + '.fasta', keyNotUsed, seqNotUsed)

    if strategyName in ['stefan_UPP', 'stefan_UPPadjusted']:


        if strategyName == 'stefan_UPP':
            scores = np.load("./data/internalData/" + dataFolderName + "/hmmScores/full.npy")

        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        #unique1 = np.unique(treeSum)

        #oneMask = np.zeros(scores.shape)
        #oneMask[:, treeSum==1] = 1
        maxHMM = np.argmax(scores, axis=1).astype(int)

        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)

        HMMinverse = np.array([np.arange(HMMinverse.shape[0]), HMMinverse]).T

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", HMMinverse)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)
        #np.save("./alignData/queryToHmm/Unique_" + strategyName + ".npy", HMMunique)

        #keysEverUsed = np.array([])
        #sequenceFileNames = giveSequenceFileNames()
        for a in range(0, len(HMMunique)):

            HMMnum = HMMunique[a]
            name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(a) + '.fasta'

            #data = loadFastaFormat(name)
            #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", data)


            keys, seqs = loadFastaBasic(name)
            #seqs = removeEmptyColumns(seqs)
            saveFastaBasic("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys, seqs)


    if strategyName == 'stefan_allwaysRoot':

        scores = scores / np.max(scores)

        HMMroot = np.zeros(scores.shape[0]).astype(int)

        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMroot)



        #name = "./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq.fasta"
        name = "./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(0) + ".fasta"
        data = loadFastaFormat(name)
        #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(0) + ".fasta", data)

    if strategyName in ['stefan_removePoison', 'stefan_removePoisonLeaf', 'stefan_removePoisonInternal', 'stefan_removePoisonGrandpa']:
        #quit()
        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        #unique1 = np.unique(treeSum)

        oneMask = np.zeros(scores.shape)
        oneMask[:, treeSum==1] = 1
        maxHMM = np.argmax(scores, axis=1)

        fatherNode = findFatherNode()

        usedOneMask3 = np.copy(treeData[maxHMM])

        if strategyName == 'stefan_removePoisonGrandpa':
            maxHMM2 = fatherNode[maxHMM]
            maxHMM2[maxHMM == 0] = 0

            maxHMM3 = fatherNode[maxHMM2]
            maxHMM3[maxHMM2 == 0] = 0

            maxHMM = np.copy(maxHMM3)


        usedOneMask2 = doPoisonRemoval(maxHMM, scores, np.arange(scores.shape[0]))

        usedOneMask2 = usedOneMask2 + usedOneMask3
        usedOneMask2[usedOneMask2>1] = 1


        #usedOneMask2[queryArgs, poisonNode] = 0
        aboveMaskOne = np.copy(usedOneMask2)

        #print (usedOneMask2[244])

        #plt.plot(usedOneMask2[244]+1)
        #plt.show()

        #print (np.sum(usedOneMask2, axis=1))

        #21
        #quit()

        aboveMaskOneString = []
        for a in range(0, len(aboveMaskOne)):
            string1 = ''
            for b in range(0, len(aboveMaskOne[a])):
                string1 = string1 + str(int(aboveMaskOne[a, b])) + ':'
            aboveMaskOneString.append(string1)
        aboveMaskOneString = np.array(aboveMaskOneString)

        _, hmmUniqueIndex, queryToHmm = np.unique(aboveMaskOneString, return_index=True, return_inverse=True)

        hmmLeafs = aboveMaskOne[hmmUniqueIndex]

        queryToHmm = np.array([np.arange(queryToHmm.shape[0]), queryToHmm]).T


        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", queryToHmm)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", queryToHmm)

        saveFromBooleanInclude(hmmLeafs, strategyName)

    if strategyName == 'stefan_SEPP':
        #quit()
        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        #unique1 = np.unique(treeSum)

        argsLeaf = np.argwhere(treeSum == 1)[:, 0]

        maxHMM = np.argmax(scores[:, argsLeaf], axis=1).astype(int)
        maxHMM = argsLeaf[maxHMM]

        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", HMMinverse)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)


        for a in range(0, len(HMMunique)):

            HMMnum = HMMunique[a]
            #name = './alignData/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            data = loadFastaFormat(name)

            #saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq.fasta", data)
            saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", data)

            #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", data)


    if strategyName in ['stefan_removeLeafExtra', 'stefan_removeLeafExtra2']:
        #quit()
        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        #unique1 = np.unique(treeSum)

        oneMask = np.zeros(scores.shape)
        oneMask[:, treeSum==1] = 1
        maxHMM = np.argmax(scores, axis=1)

        edges = findTreeEdges()

        parentNode = (edges.T)[maxHMM]
        parentNode[:, 0] = 0.5
        parentNode = np.argmax(parentNode, axis=1)

        leafOfParent = treeData[parentNode, :]

        usedMask = treeData[maxHMM, :]

        usedOneMask = usedMask * oneMask

        meanUsed = np.sum(scores * usedOneMask, axis=1) / np.sum(usedOneMask, axis=1)
        meanUsed2 = meanUsed.repeat(scores.shape[1]).reshape(scores.shape)


        aboveMask = np.zeros(scores.shape)
        #aboveMask[((scores + 0.01)/(meanUsed2 + 0.01)) > (2/3)] = 1
        if strategyName == 'stefan_removeLeafExtra':
            aboveMask[((scores + 0.01)/(meanUsed2 + 0.01)) > (1.0)] = 1


        if strategyName == 'stefan_removeLeafExtra2':
            aboveMask2 = np.copy(aboveMask)
            print ("Hi")
            aboveMask[((scores + 0.01)/(meanUsed2 + 0.01)) > (1.0)] = 1
            aboveMask2[((scores + 0.01)/(meanUsed2 + 0.01)) > (0.8)] = 1


        #aboveMaskOne = aboveMask * oneMask
        aboveMaskOne = aboveMask * oneMask * leafOfParent

        if strategyName == 'stefan_removeLeafExtra2':

            aboveMaskOne = aboveMaskOne + (usedMask * oneMask * aboveMask2)
        else:
            aboveMaskOne = aboveMaskOne + (usedMask * oneMask)
        aboveMaskOne[aboveMaskOne > 1] = 1
        aboveMaskOne = aboveMaskOne.astype(int)




        aboveMaskOneString = []
        for a in range(0, len(aboveMaskOne)):
            string1 = ''
            for b in range(0, len(aboveMaskOne[a])):
                string1 = string1 + str(int(aboveMaskOne[a, b])) + ':'
            aboveMaskOneString.append(string1)
        aboveMaskOneString = np.array(aboveMaskOneString)

        _, hmmUniqueIndex, queryToHmm = np.unique(aboveMaskOneString, return_index=True, return_inverse=True)

        hmmLeafs = aboveMaskOne[hmmUniqueIndex]

        queryToHmm = np.array([np.arange(queryToHmm.shape[0]), queryToHmm]).T

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", queryToHmm)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", queryToHmm)

        #keysEverUsed = np.array([])
        sequenceFileNames = giveSequenceFileNames()
        for a in range(0, len(hmmLeafs)):
            newFasta = []
            for b in range(0, len(hmmLeafs[0])):
                if hmmLeafs[a, b] == 1:
                    #name = sequenceFileNames[b]
                    #name = './alignData/hmmSeqAlign/' + str(b) + '.fasta'
                    name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(b) + '.fasta'
                    data = loadFastaFormat(name)
                    keysUsed, _ = loadFastaBasic(name)
                    #keysEverUsed = np.concatenate((keysEverUsed, keysUsed))
                    newFasta = newFasta + data

            #saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq.fasta", newFasta)
            saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", newFasta)
            #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", newFasta)

        #keysEverUsed = np.unique(keysEverUsed)
        #keysAll, seqsAll = loadFastaBasic('./alignData/hmmSeqAlign/' + str(0) + '.fasta')
        #seqNotUsed = seqsAll[np.isin(keysAll, keysEverUsed) == False]
        #keyNotUsed = keysAll[np.isin(keysAll, keysEverUsed) == False]
        #if keyNotUsed.shape[0] != 0:
        #    saveFastaBasic('./alignData/newHMM/unusedSeq/' + strategyName + '.fasta', keyNotUsed, seqNotUsed)



    if strategyName == 'stefan_brotherMixture':
        #quit()
        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        #unique1 = np.unique(treeSum)

        #oneMask = np.zeros(scores.shape)
        #oneMask[:, treeSum==1] = 1
        maxHMM = np.argmax(scores, axis=1).astype(int)

        brotherNode = findBrotherNode()

        maxHMM_brother = brotherNode[maxHMM]
        maxHMM_brother[maxHMM == 0] = 0




        HMMunique, HMMindex, HMMinverse = np.unique(maxHMM, return_index=True, return_inverse=True)

        HMMunique_brother = maxHMM_brother[HMMindex]

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", HMMinverse)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)

        #keysEverUsed = np.array([])
        #sequenceFileNames = giveSequenceFileNames()
        for a in range(0, len(HMMunique)):

            HMMnum = HMMunique[a]
            HMMnum_brother = HMMunique_brother[a]
            #name = './alignData/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum_brother) + '.fasta'
            #name_brother = './alignData/hmmSeqAlign/' + str(HMMnum_brother) + '.fasta'

            keys, seqs = loadFastaBasic(name)
            keys_brother, seqs_brother = loadFastaBasic(name)

            size1 = len(keys_brother)
            size2 = int(size1 * 0.5)
            randomSubset = np.random.choice(size1, size=size2, replace=False)

            keys_brother, seqs_brother = keys_brother[randomSubset], seqs_brother[randomSubset]

            keys_all = np.concatenate((keys, keys_brother))
            seqs_all = np.concatenate((seqs, seqs_brother))


            #seqs = removeEmptyColumns(seqs)
            #saveFastaBasic("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", keys_all, seqs_all)

            saveFastaBasic("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys_all, seqs_all)



    if strategyName in ['stefan_addUnpoison', 'stefan_addUnpoisonUncle']:

        #quit()
        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        #unique1 = np.unique(treeSum)

        oneMask = np.zeros(scores.shape)
        oneMask[:, treeSum==1] = 1
        maxHMM_original = np.argmax(scores, axis=1)

        usedMask_original = treeData[maxHMM_original]

        brotherNode = findBrotherNode()
        fatherNode = findFatherNode()

        if strategyName == 'stefan_addUnpoisonUncle':
            maxHMM1 = fatherNode[maxHMM_original]
            maxHMM1[maxHMM_original == 0] = 0

            maxHMM = brotherNode[maxHMM1]
            maxHMM[maxHMM1 == 0] = 0



        if strategyName == 'stefan_addUnpoison':
            maxHMM = brotherNode[maxHMM_original]
            maxHMM[maxHMM_original == 0] = 0


        toAdd = doPoisonRemoval(maxHMM, scores, np.arange(scores.shape[0]))


        aboveMaskOne = np.copy(usedMask_original * oneMask)
        aboveMaskOne = aboveMaskOne + toAdd
        aboveMaskOne[aboveMaskOne > 1] = 1


        aboveMaskOneString = []
        for a in range(0, len(aboveMaskOne)):
            string1 = ''
            for b in range(0, len(aboveMaskOne[a])):
                string1 = string1 + str(int(aboveMaskOne[a, b])) + ':'
            aboveMaskOneString.append(string1)
        aboveMaskOneString = np.array(aboveMaskOneString)

        _, hmmUniqueIndex, queryToHmm = np.unique(aboveMaskOneString, return_index=True, return_inverse=True)

        hmmLeafs = aboveMaskOne[hmmUniqueIndex]

        queryToHmm = np.array([np.arange(queryToHmm.shape[0]), queryToHmm]).T

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", queryToHmm)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", queryToHmm)

        #keysEverUsed = np.array([])
        sequenceFileNames = giveSequenceFileNames()
        for a in range(0, len(hmmLeafs)):
            newFasta = []
            for b in range(0, len(hmmLeafs[0])):
                if hmmLeafs[a, b] == 1:
                    #name = sequenceFileNames[b]
                    #name = './alignData/hmmSeqAlign/' + str(b) + '.fasta'
                    name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(b) + '.fasta'
                    data = loadFastaFormat(name)
                    keysUsed, _ = loadFastaBasic(name)
                    #keysEverUsed = np.concatenate((keysEverUsed, keysUsed))
                    newFasta = newFasta + data


            #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", newFasta)
            saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", newFasta)

        #keysEverUsed = np.unique(keysEverUsed)
        #keysAll, seqsAll = loadFastaBasic('./alignData/hmmSeqAlign/' + str(0) + '.fasta')
        #seqNotUsed = seqsAll[np.isin(keysAll, keysEverUsed) == False]
        #keyNotUsed = keysAll[np.isin(keysAll, keysEverUsed) == False]
        #if keyNotUsed.shape[0] != 0:
        #    saveFastaBasic('./alignData/newHMM/unusedSeq/' + strategyName + '.fasta', keyNotUsed, seqNotUsed)


    if strategyName == 'stefan_removePoisonGrandpa3':

        #quit()
        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        #unique1 = np.unique(treeSum)

        oneMask = np.zeros(scores.shape)
        oneMask[:, treeSum==1] = 1
        maxHMM_original = np.argmax(scores, axis=1)

        usedMask_original = treeData[maxHMM_original]

        #brotherNode = findBrotherNode()
        fatherNode = findFatherNode()

        maxHMM1 = fatherNode[maxHMM_original]
        maxHMM1[maxHMM_original == 0] = 0

        maxHMM2 = fatherNode[maxHMM1]
        maxHMM2[maxHMM1 == 0] = 0

        toAdd1 = doPoisonRemoval(maxHMM_original, scores)
        toAdd2 = doPoisonRemoval(maxHMM1, scores)


        aboveMaskOne = np.copy(usedMask_original * oneMask)

        #aboveMaskOne1 = aboveMaskOne + toAdd1
        #aboveMaskOne1[aboveMaskOne1 > 1] = 1

        aboveMaskOne1 = np.copy(toAdd1)

        aboveMaskOne2 = aboveMaskOne + toAdd2
        aboveMaskOne2[aboveMaskOne2 > 1] = 1

        aboveMaskOne = np.concatenate((aboveMaskOne1, aboveMaskOne2), axis=0)


        aboveMaskOneString = []
        for a in range(0, len(aboveMaskOne)):
            string1 = ''
            for b in range(0, len(aboveMaskOne[a])):
                string1 = string1 + str(int(aboveMaskOne[a, b])) + ':'
            aboveMaskOneString.append(string1)
        aboveMaskOneString = np.array(aboveMaskOneString)

        queryList = np.concatenate((np.arange(scores.shape[0]), np.arange(scores.shape[0])))


        _, hmmUniqueIndex, queryToHmm = np.unique(aboveMaskOneString, return_index=True, return_inverse=True)

        hmmLeafs = aboveMaskOne[hmmUniqueIndex]

        queryToHmm = np.array([queryList, queryToHmm]).T

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", queryToHmm)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", queryToHmm)

        #keysEverUsed = np.array([])
        sequenceFileNames = giveSequenceFileNames()
        for a in range(0, len(hmmLeafs)):
            newFasta = []
            for b in range(0, len(hmmLeafs[0])):
                if hmmLeafs[a, b] == 1:
                    #name = sequenceFileNames[b]
                    #name = './alignData/hmmSeqAlign/' + str(b) + '.fasta'
                    name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(b) + '.fasta'
                    data = loadFastaFormat(name)
                    keysUsed, _ = loadFastaBasic(name)
                    #keysEverUsed = np.concatenate((keysEverUsed, keysUsed))
                    newFasta = newFasta + data


            #saveFasta("./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta", newFasta)
            #saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq.fasta", newFasta)
            saveFasta("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", newFasta)

        #keysEverUsed = np.unique(keysEverUsed)
        #keysAll, seqsAll = loadFastaBasic('./alignData/hmmSeqAlign/' + str(0) + '.fasta')
        #seqNotUsed = seqsAll[np.isin(keysAll, keysEverUsed) == False]
        #keyNotUsed = keysAll[np.isin(keysAll, keysEverUsed) == False]
        #if keyNotUsed.shape[0] != 0:
        #    saveFastaBasic('./alignData/newHMM/unusedSeq/' + strategyName + '.fasta', keyNotUsed, seqNotUsed)


    if strategyName in ['stefan_highBiasedUPP', 'stefan_highBiasedUPP2']:

        #quit()
        scores = scores / np.max(scores)

        treeData = findDecomposition()

        treeSum = np.sum(treeData, axis=1)

        argsInternal = np.argwhere(treeSum != 1)[:, 0]
        treeDataLeaf = np.copy(treeData)
        treeDataLeaf[:, argsInternal] = 0
        treeSumLeaf = np.sum(treeDataLeaf, axis=1)

        if strategyName == 'stefan_highBiasedUPP':
            heightShift = np.log2(treeSumLeaf) * 0.01
        if strategyName == 'stefan_highBiasedUPP2':
            heightShift = np.log2(treeSumLeaf - 0.75) * 0.005


        for a in range(len(scores[0])):
            scores[:, a] = scores[:, a] + heightShift[a]

        #unique1 = np.unique(treeSum)

        #oneMask = np.zeros(scores.shape)
        #oneMask[:, treeSum==1] = 1
        maxHMM = np.argmax(scores, axis=1).astype(int)

        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)

        HMMinverse = np.array([np.arange(HMMinverse.shape[0]), HMMinverse]).T

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", HMMinverse)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)
        #np.save("./alignData/queryToHmm/Unique_" + strategyName + ".npy", HMMunique)

        #keysEverUsed = np.array([])
        #sequenceFileNames = giveSequenceFileNames()
        for a in range(0, len(HMMunique)):

            HMMnum = HMMunique[a]
            #name = './alignData/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum) + '.fasta'

            keys, seqs = loadFastaBasic(name)
            #seqs = removeEmptyColumns(seqs)
            saveFastaBasic("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys, seqs)


    if strategyName in ['stefan_fastUPP', 'stefan_fastUPPexception', 'stefan_fastUPPearly']:

        #quit()
        if strategyName == 'stefan_fastUPPearly':
            scoresFull = np.load('./alignData/Searcher/scoreFiles/Early_score.npy')
        else:
            scoresFull = np.load('./alignData/Searcher/scoreFiles/score.npy')

        maxHMM = np.argmax(scoresFull, axis=1).astype(int)

        if strategyName == 'stefan_fastUPPexception':
            scores_original = np.load("./alignData/hmmScores/full.npy")

            treeData = findDecomposition()
            treeSum = np.sum(treeData, axis=1)
            treeTop = np.unique(treeSum)[-3:]

            maxHMM2 = np.argmax(scores_original, axis=1).astype(int)
            argsReplace = np.argwhere(np.isin(treeSum[maxHMM], treeTop))[:, 0]
            maxHMM[argsReplace] = maxHMM2[argsReplace]



        HMMunique, HMMinverse = np.unique(maxHMM, return_inverse=True)

        HMMinverse = np.array([np.arange(HMMinverse.shape[0]), HMMinverse]).T

        #np.save("./alignData/queryToHmm/original/" + strategyName + ".npy", HMMinverse)
        np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy", HMMinverse)

        for a in range(0, len(HMMunique)):

            HMMnum = HMMunique[a]
            #name = './alignData/hmmSeqAlign/' + str(HMMnum) + '.fasta'
            name = "./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(HMMnum) + '.fasta'

            keys, seqs = loadFastaBasic(name)
            #seqs = removeEmptyColumns(seqs)
            saveFastaBasic("./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta", keys, seqs)




    if strategyName == 'stefan_printTesting':

        #quit()
        scores = scores - np.min(scores)
        scores = scores / np.max(scores)

        #plt.imshow(scores[:50])
        #plt.show()
        #quit()

        #scores = scores[:10]

        treeData = findDecomposition()

        print ("B")

        treeSum = np.sum(treeData, axis=1)

        '''
        argsSmall = np.argwhere(treeSum == 1)[:, 0]
        for argSmall in argsSmall:
            name = './alignData/hmmSeqAlign/' + str(argSmall) + '.fasta'
            data = loadFastaFormat(name)
            print (len(data))
        quit()
        '''

        #unique1 = np.unique(treeSum)

        maxHMM = np.argmax(scores, axis=1)

        queryNum = np.arange(maxHMM.shape[0])
        nonPoison = doPoisonRemoval(maxHMM, scores, queryNum)

        #print (maxHMM.shape)
        #print (np.unique(treeSum[maxHMM], return_counts=True))
        #quit()

        oneMask = np.zeros(scores.shape)
        oneMask[:, treeSum==1] = 1
        maxHMM = np.argmax(scores, axis=1)

        usedMask = treeData[maxHMM, :]

        for a in range(0, len(usedMask)):
            nonPoisonNow = nonPoison[a]
            usedNow = np.argwhere(usedMask[a] == 1)[:, 0]

            treeNow = treeData[usedNow, :][:, usedNow]

            nonPoisonNow = nonPoisonNow[usedNow].repeat(usedNow.shape[0]).reshape((usedNow.shape[0], usedNow.shape[0]))


            #usedOther = np.argwhere(np.sum(treeNow, axis=0) != 0)[:, 0]
            #treeNow = treeNow[:, usedOther]

            scores2 = scores[a, usedNow].repeat(usedNow.shape[0]).reshape((usedNow.shape[0], usedNow.shape[0]))


            plt.imshow(treeNow * scores2 * nonPoisonNow)
            plt.show()

        #img = interweive(np.array([usedMask, scores]))

        #plt.imshow(img[:50])
        #plt.show()
        quit()





#scoresToHMMSeq('stefan_basic')
#quit()

#scoresToHMMSeq('stefan_printTesting')
#scoresToHMMSeq('stefan_UPP')
#scoresToHMMSeq('stefan_UPPadjusted')

#scoresToHMMSeq('stefan_allwaysRoot')
#scoresToHMMSeq('stefan_removeLeaf')
#scoresToHMMSeq('stefan_onlyOneAbove')
#scoresToHMMSeq('stefan_useAboveLeaf')
#scoresToHMMSeq('stefan_removePoison')
#scoresToHMMSeq('stefan_removePoisonLeaf')
#scoresToHMMSeq('stefan_removePoisonInternal')


#scoresToHMMSeq('stefan_SEPP')
#scoresToHMMSeq('stefan_removeLeafExtra')
#scoresToHMMSeq('stefan_removeLeafExtra2')
#scoresToHMMSeq('stefan_addUnpoison')
#scoresToHMMSeq('stefan_removePoisonGrandpa')

#scoresToHMMSeq('stefan_addUnpoisonUncle')
#scoresToHMMSeq('stefan_brotherMixture')

#scoresToHMMSeq('stefan_removePoisonGrandpaTwo')
#scoresToHMMSeq('stefan_removePoisonGrandpa3')

#scoresToHMMSeq('stefan_highBiasedUPP')
#scoresToHMMSeq('stefan_highBiasedUPP2')

#scoresToHMMSeq('stefan_fastUPP')

#scoresToHMMSeq('stefan_fastUPPexception')
#scoresToHMMSeq('stefan_fastUPPearly')
#quit()





def generateNewHMM(strategyName):

    dataFolderName = giveAllFileNames()[4]



    #numHMM = int(np.max(np.load("./alignData/queryToHmm/original/" + strategyName + ".npy")[:, 1])+1)
    numHMM = int(np.max(np.load("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")[:, 1])+1)

    #print (numHMM)
    #quit()

    for a in range(0, numHMM):
        #print ("A: ", a)
        #src = "./alignData/newHMM/newHMMseq/" + strategyName + "_" + str(a) + ".fasta"
        src = "./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/newHMMseq/" + str(a) + ".fasta"
        #data = loadFastaFormat(src)

        #print ("A")
        #keys, seqs = loadFastaBasic(src)
        #print (len(seqs[0]))
        #seqs = removeEmptyColumns(seqs)
        #print (len(seqs[0]))

        #'''
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
                #size1 += 1
                usedCols.append(b)

        usedCols = np.array(usedCols)
        usedCols_filename = "./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/columnSets/" + str(a) + ".npy"
        ensureFolder(usedCols_filename)
        np.save(usedCols_filename, usedCols)
        #'''

        #os.system("hmmbuild ./alignData/newHMM/hmm/" + strategyName + "_" + str(a) + ".hmm " + src)
        #os.system("./hmmer-3.0/src/hmmbuild --symfrac 0.0 --informat afa ./alignData/newHMM/hmm/" + strategyName + "_" + str(a) + ".hmm " + src)

        #hmmName = "./alignData/newHMM/hmm/" + strategyName + "_" + str(a) + ".hmm"
        hmmName = "./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(a) + ".hmm"

        print ("HI")
        ensureFolder(hmmName)

        runHMMbuild(hmmName, src)


#generateNewHMM('stefan_UPP')
#quit()



def resortToUPP(strategyName, doResort=True):

    dataFolderName = giveAllFileNames()[4]

    scoreStrat_file = "./data/internalData/" + dataFolderName + "/" + strategyName + "/hmmScores/scoresFull/full.npy"
    ensureFolder(scoreStrat_file)
    scoreStrat = np.load(scoreStrat_file)
    #scoreUPP = np.load("./alignData/hmmScores/full.npy")
    #scoreUPP = np.load("./alignData/hmmScores/fullAdjusted.npy")
    scoreUPP_file = "./data/internalData/" + dataFolderName + "/hmmScores/fullAdjusted.npy"
    ensureFolder(scoreUPP_file)
    scoreUPP = np.load(scoreUPP_file)
    if not doResort:
        scoreUPP[:] = -1000

    diff = np.max(scoreStrat, axis=1) - np.max(scoreUPP, axis=1)

    uppChoice = np.argmax(scoreUPP, axis=1)
    strategyChoice = np.argmax(scoreStrat, axis=1)

    argsUppBetter = np.argwhere(diff < 0)[:, 0]

    uppToUse = uppChoice[argsUppBetter]

    #np.save('./alignData/hmmScores/useUPP/queryNums.npy')

    _, uppToUseInverse = np.unique(uppToUse, return_inverse=True)

    #HMMinverse_old = np.load("./alignData/queryToHmm/original/" + strategyName + ".npy")
    HMMinverse_old = np.load("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/original.npy")
    max1 = np.max(HMMinverse_old[:, 1])

    HMMinverse = np.copy(strategyChoice)
    #HMMinverse[argsUppBetter] = np.max(HMMinverse) + 1 + uppToUseInverse
    HMMinverse[argsUppBetter] = max1 + 1 + uppToUseInverse



    #TODO what about HMM which are now removed


    uppPosition = np.zeros(HMMinverse.shape[0]) - 1
    uppPosition[argsUppBetter] = uppChoice[argsUppBetter]

    HMMinverse_file = "./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy"
    ensureFolder(HMMinverse_file)

    np.save(HMMinverse_file, HMMinverse)
    #"./alignData/queryToHmm/withUPP/UPPused_" + strategyName + ".npy"
    np.save("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/UPPused.npy", uppPosition)


    uppToUseUnique, index1 = np.unique(uppToUse, return_index=True)

    newNumsFull = []

    for a in range(len(uppToUseUnique)):

        uppNum = uppToUseUnique[a]


        newNum = HMMinverse[argsUppBetter[index1[a]]]

        newNumsFull.append(newNum)


        keys, seqs = loadFastaBasic("./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(uppNum) + '.fasta')
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
                #size1 += 1
                usedCols.append(b)

        usedCols = np.array(usedCols)
        #usedCols_file = './alignData/newHMM/columnSets/'  + strategyName + "_" + str(newNum) + ".npy"
        usedCols_file = "./data/internalData/" + dataFolderName + "/" + strategyName + '/newHMM/columnSets/' + str(newNum) + ".npy"
        ensureFolder(usedCols_file)
        #np.save('./alignData/newHMM/columnSets/'  + strategyName + "_" + str(newNum) + ".npy", usedCols)
        np.save(usedCols_file, usedCols)

    newNumsFull = np.array(newNumsFull)
    theoryNewNums = HMMinverse[HMMinverse > max1]

    diffIssue = np.unique(newNumsFull) - np.unique(theoryNewNums)
    diffIssue = np.sum(np.abs(diffIssue))

    assert diffIssue == 0



#resortToUPP('stefan_removePoison')
#quit()


def stockholmToFasta(stockholmName, fastaName):

    from Bio import SeqIO

    fakeFastaName = './alignData/temporaryFileSave/fakeFasta_1.fasta'
    records = SeqIO.parse(stockholmName, "stockholm")
    count = SeqIO.write(records, fakeFastaName, "fasta")

    dataAlign = loadFastaFormat(fakeFastaName)

    newData = ''.join(dataAlign)
    newData = newData.replace('\n', '')

    newData = newData.split('>')[1:]

    #print (len(newData[0]))

    keys = []
    seqs = []
    for a in range(0, len(newData)):
        keys.append(newData[a][:10])
        seqs.append(newData[a][10:])

    #print (len(seqs[0]))

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

        #print ([i])

        #print (i)

        if started:

            position = (line_number - startedLine - 1) % (secondLine - startedLine)
            if position >= len(data):
                data.append([])
            #print (position)
            data[position].append(i)

        if (i == '\n') and ((started == True) and (second == False)):
            second = True
            secondLine = line_number
            #print ('second ', line_number)

        if (i == '\n') and (started == False):
            started = True
            startedLine = line_number
            #print ('start ', line_number)


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

    #lowercase = np.array(list('abcdefghijklmnopqrstuvwxyz'))

    #print (originalLength)
    #quit()

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
        #key = key.replace(' ', '')
        #print ([key])
        #quit()

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

    #print (keyLengths[:5])
    #print (newData[:1000])
    #quit()

    newData = newData.replace('\n', '')

    newData = newData.split('>')[1:]

    #print (newData[0])
    #quit()

    #print (len('AAAAAAAAZY'))
    #quit()

    keys = []
    seqs = []
    for a in range(0, len(newData)):
        #if a < 5:
        #    print (newData[a])
        #    print (keyLengths[a])
        keys.append(newData[a][:keyLengths[a]])
        seqs.append(newData[a][keyLengths[a]:])
        #keys.append(newData[a][:10])
        #seqs.append(newData[a][10:])

    #print (keys[:10])
    #quit()

    keys, seqs = np.array(keys), np.array(seqs)

    saveFastaBasic(fastaName, keys, seqs)

#txtToFasta('./alignData/trueAlign/true_align.txt', './alignData/trueAlign/true_align_fasta.fasta')
#txtToFasta('./alignData/trueAlign/true_align_fragged.txt', './alignData/trueAlign/true_align_fragged_fasta.fasta')
#quit()



def alignQueries(strategyName):

    import copy

    dataFolderName = giveAllFileNames()[4]

    queryName = giveQueryFileName()
    queryData = loadFastaFormat(queryName)
    #queryData = loadFastaFormat('./alignData/tmpfiles/output.ped7oo42/backbone/queryqqpulswm.fas')


    #queryToHmm = np.load("./alignData/queryToHmm/" + strategyName + ".npy")

    #queryToHmm = np.load("./alignData/queryToHmm/withUPP/HMMused_" + strategyName + ".npy")
    #uppHMM = np.load("./alignData/queryToHmm/withUPP/UPPused_" + strategyName + ".npy")

    queryToHmm = np.load("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy")
    uppHMM = np.load("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/UPPused.npy")



    predictionDataFull = []

    for b in range(0, int(np.max(queryToHmm)) + 1):
        #print (b)

        argsHMM = np.argwhere(queryToHmm == b)[:, 0]

        if argsHMM.shape[0] > 0:

            #argsHMM = np.concatenate((np.array([0]), argsHMM))

            #print (argsHMM.shape)

            dataNow = []
            keysNow = []

            for a in range(0, len(queryData)):
                if queryToHmm[a] == b:
                    dataNow.append(queryData[a])

            #print (dataNow[0])

            assert len(dataNow) == argsHMM.shape[0]


            queryName1 = "./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/inputQuery/' + str(b) + '.fasta'
            ensureFolder(queryName1)
            saveFasta(queryName1, dataNow)

            predictionName = "./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/predictedQuery/' + str(b) + '.sto'
            ensureFolder(predictionName)
            if uppHMM[argsHMM[0]] == -1:
                hmmName = "./data/internalData/" + dataFolderName + "/" + strategyName + "/newHMM/hmm/" + str(b) + ".hmm"
            else:
                uppNum = int(uppHMM[argsHMM[0]])
                hmmName = "./data/internalData/" + dataFolderName + "/initialHMM/test" + str(uppNum) + ".hmm "
                #hmmName = "./alignData/newHMM/hmm/" + strategyName + "_" + str(uppNum) + ".hmm"

            ensureFolder(hmmName)

            #os.system("hmmalign -o " + predictionName + " " + hmmName + " " + name)
            #os.system("hmmalign --allcol --dna " + predictionName + " " + hmmName + " " + name)
            #os.system("hmmalign --allcol " + hmmName + " " + name + ' ' + hmmName + " ")
            #os.system("./hmmer-3.0/src/hmmalign --allcol  --dna " + predictionName + " " + hmmName + ' ' + name + " ")
            #os.system("./hmmer-3.0/src/hmmalign --allcol  --dna " + hmmName + ' ' + name + " ")

            #os.system("./hmmer-3.0/src/hmmalign --allcol  --dna -o " + predictionName + " " + hmmName + ' ' + name + " ")
            runHMMalign(hmmName, queryName1, predictionName)




#alignQueries('stefan_UPP')
#quit()

def mergeAlignments(strategyName, overlapLowercase=True):

    dataFolderName = giveAllFileNames()[4]

    #queryToHmm = np.load("./alignData/queryToHmm/" + strategyName + ".npy")
    queryToHmm = np.load("./data/internalData/" + dataFolderName + "/" + strategyName + "/queryToHmm/withUPP/HMMused.npy")

    predictionDataFull = []

    backboneKeys, backboneSeqs = loadFastaBasic("./data/internalData/" + dataFolderName + '/hmmSeqAlign/' + str(0) + '.fasta')

    #colInserts = np.zeros(len(backboneSeqs[0])+1)

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
            #a = 1
            predictionName = "./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/predictedQuery/'  + str(a) + '.sto'
            _, insertions = loadStockholm(predictionName)
            insertions = insertions[-1]

            #print (insertions)
            #quit()

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


            queryFileName = "./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/inputQuery/'  + str(a) + '.fasta'
            #queryFileName = giveQueryFileName()
            queryKey, querySeq = loadFastaBasic(queryFileName)

            assert len(keys) == len(queryKey)
            for b in range(len(queryKey)):
                assert queryKey[b] == keys[b]


            #print (seqs[0])
            #quit()

            usedCols = np.load("./data/internalData/" + dataFolderName + "/" + strategyName + '/newHMM/columnSets/' + str(a) + ".npy")
            matchPosition = np.argwhere(insertions == 0)[:, 0]

            #argsHMM2 = argsHMM.repeat(usedCols.shape[0])
            #usedCols2 = usedCols.repeat(argsHMM.shape[0]).reshape((usedCols.shape[0], argsHMM.shape[0])).T.reshape((usedCols.shape[0] * argsHMM.shape[0],))

            try:
                assert len(usedCols) == len(matchPosition)
            except:
                print ("This issue means that the number of match states from the HMM is")
                print ("not equal to the number of non-empty columns on the sequences")
                print ("the hmm was trained on.")
                assert len(usedCols) == len(matchPosition)

            backBoneChoice[np.ix_(argsHMM, usedCols)] = seqsArray[:, matchPosition]


            backBoneChoiceBool[np.ix_(argsHMM, usedCols)] = 1
            #backBoneChoice[argsHMM2, usedCols2] = seqsArray[:, matchPosition].reshape((seqsArray.shape[0] * matchPosition.shape[0],))


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

                    #print (argsHMM[0])
                    #print (toInsert[0])

                    insertPosition = usedColsExtra[b]

                    fullInsertions[insertPosition].append(np.copy(toInsert))

                    fullInsertionsIndex[insertPosition].append(np.copy(argsHMM))

                    if overlapLowercase:
                        fullInsertionsNumber[insertPosition] = max(insertNumber[b], fullInsertionsNumber[insertPosition])
                    else:
                        fullInsertionsNumber[insertPosition] += insertNumber[b]

            #quit()


    #print (str(list(backBoneChoice[0])))


    newAlignment = np.zeros((queryToHmm.shape[0], int(fullInsertionsNumber.shape[0] + np.sum(fullInsertionsNumber))  )).astype(str)
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

            #if 189 in insertQuery:
            #    print ("insertQuery")
            #    print (toInsert[0])

            newAlignment[np.ix_(insertQuery, insertPosition)] = np.copy(toInsert)
            newAlignmentBool[np.ix_(insertQuery, insertPosition)] = 1

            if not overlapLowercase:
                start2 += size1

    #print (list(newAlignment[189]))
    #print (list(newAlignment[189][newAlignmentBool[189] == 1]))

    #lowercase = np.array(list('abcdefghijklmnopqrstuvwxyz'))
    #uppercase = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))

    #for a in range(len(lowercase)):
    #    newAlignment[newAlignment == lowercase[a]] = uppercase[a]

    newAlignment[newAlignment == '.'] = '-'
    #print ('a' in newAlignment)


    newAlignment = newAlignment[:, 1:]
    newAlignmentBool = newAlignmentBool[:, 1:]
    colIndexTrue = (colIndex[1:] - 1).astype(int)

    newAlignmentString = []
    for a in range(newAlignment.shape[0]):
        str1 = ''.join(list(newAlignment[a]))
        newAlignmentString.append(copy.copy(str1))
    newAlignmentString = np.array(newAlignmentString)

    #print (newAlignmentString[2])

    #newAlignmentString = removeMultInvert(newAlignmentString)

    #print (newAlignmentString[2])
    #quit()

    fastaName = "./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignmentFasta.fasta'
    ensureFolder(fastaName)
    saveFastaBasic(fastaName, queryNames, newAlignmentString)

    np.save("./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignment.npy', newAlignment)
    np.save("./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignmentBool.npy', newAlignmentBool)
    np.save("./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/columnIndex.npy', colIndexTrue)



#mergeAlignments('stefan_removePoison')
#quit()



def InputMergeAlignments(queryNamesFull, alignments, allInsertions, columnSets, Ncolumns, overlapLowercase=True):

    #queryToHmm = np.load("./alignData/queryToHmm/" + strategyName + ".npy")
    #predictionDataFull = []
    #backboneKeys, backboneSeqs = loadFastaBasic('./alignData/hmmSeqAlign/' + str(0) + '.fasta')

    #New
    Nquery = 0
    for alignment in alignments:
        Nquery += len(alignment)
    #End New


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

        #argsHMM = np.argwhere(queryToHmm == a)[:, 0]
        sizeAlign = len(alignments[a])
        argsHMM = np.arange(sizeAlign) + queryCount
        queryCount += sizeAlign
        if argsHMM.shape[0] > 0:
            #a = 1
            #predictionName = './alignData/hmmQueryList/predictedQuery/' + strategyName + '_'  + str(a) + '.sto'
            #_, insertions = loadStockholm(predictionName)
            #insertions = insertions[-1]

            insertions = allInsertions[a]

            #print (insertions)
            #quit()

            insertions = np.array(list(insertions))
            insertions[insertions == '.'] = 1
            insertions[insertions == 'x'] = 0
            insertions = insertions.astype(int)

            #keys, seqs = loadStockholmOnlySeqs(predictionName)
            seqs = alignments[a]
            queryNames[argsHMM] = queryNamesFull[argsHMM]
            seqsArray = []
            for seq in seqs:
                seqsArray.append(list(seq))
            seqsArray = np.array(seqsArray)



            #queryFileName = './alignData/hmmQueryList/inputQuery/' + strategyName + '_'  + str(a) + '.fasta'
            #queryKey, querySeq = loadFastaBasic(queryFileName)
            #assert len(keys) == len(queryKey)
            #for b in range(len(queryKey)):
            #    assert queryKey[b] == keys[b]


            #print (seqs[0])
            #quit()

            #usedCols = np.load('./alignData/newHMM/columnSets/'  + strategyName + "_" + str(a) + ".npy")
            usedCols = columnSets[a]
            matchPosition = np.argwhere(insertions == 0)[:, 0]

            #argsHMM2 = argsHMM.repeat(usedCols.shape[0])
            #usedCols2 = usedCols.repeat(argsHMM.shape[0]).reshape((usedCols.shape[0], argsHMM.shape[0])).T.reshape((usedCols.shape[0] * argsHMM.shape[0],))

            backBoneChoice[np.ix_(argsHMM, usedCols)] = seqsArray[:, matchPosition]
            backBoneChoiceBool[np.ix_(argsHMM, usedCols)] = 1
            #backBoneChoice[argsHMM2, usedCols2] = seqsArray[:, matchPosition].reshape((seqsArray.shape[0] * matchPosition.shape[0],))


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

                    #print (argsHMM[0])
                    #print (toInsert[0])

                    insertPosition = usedColsExtra[b]

                    fullInsertions[insertPosition].append(np.copy(toInsert))

                    fullInsertionsIndex[insertPosition].append(np.copy(argsHMM))

                    if overlapLowercase:
                        fullInsertionsNumber[insertPosition] = max(insertNumber[b], fullInsertionsNumber[insertPosition])
                    else:
                        fullInsertionsNumber[insertPosition] += insertNumber[b]

            #quit()


    #print (str(list(backBoneChoice[0])))


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

            #if 189 in insertQuery:
            #    print ("insertQuery")
            #    print (toInsert[0])

            newAlignment[np.ix_(insertQuery, insertPosition)] = np.copy(toInsert)
            newAlignmentBool[np.ix_(insertQuery, insertPosition)] = 1

            if not overlapLowercase:
                start2 += size1

    #print (list(newAlignment[189]))
    #print (list(newAlignment[189][newAlignmentBool[189] == 1]))

    #lowercase = np.array(list('abcdefghijklmnopqrstuvwxyz'))
    #uppercase = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))

    #for a in range(len(lowercase)):
    #    newAlignment[newAlignment == lowercase[a]] = uppercase[a]

    newAlignment[newAlignment == '.'] = '-'
    #print ('a' in newAlignment)


    newAlignment = newAlignment[:, 1:]
    newAlignmentBool = newAlignmentBool[:, 1:]
    colIndexTrue = (colIndex[1:] - 1).astype(int)

    newAlignmentString = []
    for a in range(newAlignment.shape[0]):
        str1 = ''.join(list(newAlignment[a]))
        newAlignmentString.append(copy.copy(str1))
    newAlignmentString = np.array(newAlignmentString)

    #print (newAlignmentString[2])

    #newAlignmentString = removeMultInvert(newAlignmentString)

    #print (newAlignmentString[2])
    #quit()

    #print (newAlignmentString)

    #fastaName = './alignData/hmmQueryList/merged/' + strategyName + '_alignmentFasta.fasta'
    #saveFastaBasic(fastaName, queryNames, newAlignmentString)

    #np.save('./alignData/hmmQueryList/merged/' + strategyName + '_alignment.npy', newAlignment)
    #np.save('./alignData/hmmQueryList/merged/' + strategyName + '_alignmentBool.npy', newAlignmentBool)
    #np.save('./alignData/hmmQueryList/merged/' + strategyName + '_columnIndex.npy', colIndexTrue)

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

    #alignments = [np.array(['tgca-T--acg']), np.array(['tgAATA-tc'])]
    #allInsertions = ['....xxxx...', '..xxxxx..']
    #columnSets = [np.array([0, 2, 3, 4]), np.array([0, 1, 2, 3, 4])]
    #Ncolumns = 5
    #queryNamesFull = np.array(["X", "Y"])

    seqsNew = InputMergeAlignments(queryNamesFull, alignments, allInsertions, columnSets, Ncolumns, overlapLowercase=overlapLowercase)

    return seqsNew




def checkAlignmentsValid(strategyName):

    newAlignment = np.load('./alignData/hmmQueryList/merged/' + strategyName + '_alignment.npy')
    newAlignmentBool = np.load('./alignData/hmmQueryList/merged/' + strategyName + '_alignmentBool.npy')


    queryToHmm = np.load("./alignData/queryToHmm/" + strategyName + ".npy")


    for a in range(0, int(np.max(queryToHmm)) + 1):

        argsHMM = np.argwhere(queryToHmm == a)[:, 0]
        if argsHMM.shape[0] > 0:
            #a = 1
            predictionName = './alignData/hmmQueryList/predictedQuery/' + strategyName + '_'  + str(a) + '.sto'
            #_, insertions = loadStockholm(predictionName)
            #insertions = insertions[-1]
            #insertions = np.array(list(insertions))

            keys, seqs = loadStockholmOnlySeqs(predictionName)
            seqsArray = []
            for b in range(0, len(seqs)):
                seq = seqs[b]
                indexInNew = argsHMM[b]
                #print (seq)
                seq2 = np.array(list(seq))
                newSeq = newAlignment[indexInNew][newAlignmentBool[indexInNew] != 0]
                newSeqBool = newAlignmentBool[indexInNew][newAlignmentBool[indexInNew] != 0]
                #print (list(newAlignment[indexInNew]))

                #newSeq = np.sort(newSeq)
                #seq2 = np.sort(seq2)

                #print (len(seq2))

                for c in range(0, len(seq2)):
                    if seq2[c] != newSeq[c]:
                        print ("Issue")
                        quit()
                #quit()
                #print (list(newSeq))
                #quit()

    print ("All Good")

#checkAlignmentsValid('stefan_UPP')
#quit()


def buildAlignMerge(strategyName, doResort=True):
    generateNewHMM(strategyName)

    saveNewScores(strategyName)
    resortToUPP(strategyName, doResort=doResort)

    alignQueries(strategyName)
    mergeAlignments(strategyName)



#buildAlignMerge('stefan_UPP', doResort=False)
#buildAlignMerge('stefan_UPPadjusted')
#quit()

#buildAlignMerge('stefan_basic')
#buildAlignMerge('stefan_removeLeaf')
#buildAlignMerge('stefan_removePoison')
#buildAlignMerge('stefan_removePoisonLeaf')
#buildAlignMerge('stefan_removePoisonInternal')

#buildAlignMerge('stefan_SEPP')
#buildAlignMerge('stefan_removeLeafExtra')
#buildAlignMerge('stefan_removeLeafExtra2')
#buildAlignMerge('stefan_addUnpoison')
#buildAlignMerge('stefan_removePoisonGrandpa')
#buildAlignMerge('stefan_addUnpoisonUncle')
#buildAlignMerge('stefan_removePoisonGrandpa3')

#buildAlignMerge('stefan_highBiasedUPP', doResort=False)
#buildAlignMerge('stefan_highBiasedUPP2', doResort=False)


#buildAlignMerge('stefan_fastUPP', doResort=False)
#buildAlignMerge('stefan_fastUPPexception', doResort=False)
#buildAlignMerge('stefan_fastUPPearly', doResort=False)

#quit()




def saveTrueUPPSubset():

    fileNames = giveAllFileNames()
    dataFolderName = giveAllFileNames()[4]

    predictionName = fileNames[3]

    #predictionName = './alignData/UPPoutput/output_alignment.fasta'
    #predictionName = './alignData/UPPoutput/R' + str(13) + '_output_alignment.fasta'
    #predictionName = './alignData/UPPoutput/R' + str(14) + '_output_alignment.fasta'


    #predictionNameStefan = './alignData/hmmQueryList/merged/' + 'stefan_UPP' + '_alignmentFasta.fasta'
    #predictionNameNew = './alignData/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'

    predictionNameStefan = './data/internalData/' + dataFolderName + '/' + 'stefan_UPP' '/hmmQueryList/merged/alignmentFasta.fasta'
    predictionNameNew = './data/internalData/' + dataFolderName + '/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'

    ensureFolder(predictionNameStefan)
    ensureFolder(predictionNameNew)
    trueKey, trueSeq = loadFastaBasic(predictionName)
    predKey, predSeq = loadFastaBasic(predictionNameStefan)

    trueKey, trueSeq = trueKey[np.isin(trueKey, predKey)], trueSeq[np.isin(trueKey, predKey)]

    ##########trueSeq = removeMultInvert(trueSeq)


    #trueSeq = removeEmptyColumns(trueSeq) #Does not affect score at all.

    saveFastaBasic(predictionNameNew, trueKey, trueSeq)

#saveTrueUPPSubset()
#quit()

#to R19
def saveTrueAlignFasta():

    #for a in [13, 14]:#range(20):
    #    if a != 16:
    fileNames = giveAllFileNames()
    trueAlignment = fileNames[2]
    dataFolderName = giveAllFileNames()[4]

    #trueAlignment = './alignData/UnalignFragTree/low_frag/1000M1/R' + str(13) + '/true_align_fragged.txt'

    #trueFasta = './alignData/trueAlignment/original/R' + str(13) + '_true_align_fragged_fasta.fasta'
    #trueFasta = './alignData/trueAlignment/original/true_align_fragged_fasta.fasta'
    trueFasta = './data/internalData/' + dataFolderName + '/trueAlignment/true_align_fragged_fasta.fasta'
    ensureFolder(trueFasta)
    #txtToFasta(trueAlignment, trueFasta)

    #trueAlignment = './alignData/UnalignFragTree/high_frag/1000M1/R14/true_align_fragged.txt'
    #trueFasta = './alignData/trueAlignment/original/R14_true_align_fragged_fasta.fasta'
    txtToFasta(trueAlignment, trueFasta)

#saveTrueAlignFasta()
#quit()

def doAllSteps(strategyName):
    #saveInitialSteps()
    #scoresToHMMSeq(strategyName)
    buildAlignMerge(strategyName)

#doAllSteps('stefan_UPP')
#quit()


def scoreAlignment(strategyName):

    dataFolderName = giveAllFileNames()[4]

    saveTrueUPPSubset()
    saveTrueAlignFasta()



    #trueFasta = './data/internalData/' + dataFolderName + '/trueAlignment/original/true_align_fragged_fasta.fasta'
    trueFasta = './data/internalData/' + dataFolderName + '/trueAlignment/true_align_fragged_fasta.fasta'

    #predictionName = './data/internalData/' + dataFolderName + '/hmmQueryList/merged/' + strategyName + '_alignmentFasta.fasta'
    predictionName = "./data/internalData/" + dataFolderName + "/" + strategyName + '/hmmQueryList/merged/alignmentFasta.fasta'
    #predictionName = './alignData/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'
    ensureFolder(predictionName)

    UPPName = './data/internalData/' + dataFolderName + '/hmmQueryList/merged/' + 'stefan_TrueUPP' + '_alignmentFasta.fasta'
    ensureFolder(UPPName)

    trueFastaNew = './data/internalData/' + dataFolderName + '/trueAlignment/subset/true_fasta.fasta'
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

    #print (len(trueKey))




    saveFastaBasic(trueFastaNew, trueKey, trueSeq)

    #cmd = "java -jar ./FastSP.jar -r " + trueFastaNew + " -e " + predictionName
    cmd = "java -jar ./FastSP.jar -ml -r " + trueFastaNew + " -e " + predictionName

    from subprocess import run
    scoreInfo  = os.popen(cmd).readlines()
    #Example: scoreInfo = ['SP-Score 0.0\n', 'Modeler 0.0\n', 'SPFN 1.0\n', 'SPFP 1.0\n', 'Compression 0.12854722575180008\n', 'TC 0.0\n']

    scoreInfo2 = []
    for a in range(0, len(scoreInfo)):
        string1 = scoreInfo[a].split(' ')[1].replace('\n', '')
        scoreInfo2.append(float(string1))

    print (strategyName)
    print (scoreInfo)



#scoreAlignment('stefan_basic')
#quit()
# scoreAlignment('stefan_UPP')
# scoreAlignment('stefan_UPPadjusted')
#quit()
#scoreAlignment('stefan_allwaysRoot')
#scoreAlignment('stefan_removeLeaf')
#scoreAlignment('stefan_onlyOneAbove')
#scoreAlignment('stefan_useAboveLeaf')
#scoreAlignment('stefan_TrueUPP')
#scoreAlignment('stefan_removePoison')
#scoreAlignment('stefan_removePoisonLeaf')
#scoreAlignment('stefan_removePoisonInternal')
#scoreAlignment('stefan_SEPP')
#scoreAlignment('stefan_removeLeafExtra')
#scoreAlignment('stefan_removeLeafExtra2')

#scoreAlignment('stefan_removePoisonGrandpa')
#scoreAlignment('stefan_addUnpoisonUncle')
#scoreAlignment('stefan_brotherMixture')
#scoreAlignment('stefan_addUnpoison')
#scoreAlignment('stefan_removePoisonGrandpa3')



#scoreAlignment('stefan_highBiasedUPP')
#scoreAlignment('stefan_highBiasedUPP2')

#scoreAlignment('stefan_fastUPP')
#scoreAlignment('stefan_fastUPPexception')
#scoreAlignment('stefan_fastUPPearly')







#quit()

def compareToUPP():

    #UnalignFragTree

    #predictionStefan = './alignData/hmmQueryList/merged/' + 'stefan_UPP' + '_alignmentFasta.fasta'
    #predictionUPP = './alignData/UPPoutput/output_alignment.fasta'


    #predictionStefan = './alignData/hmmQueryList/merged/' + 'stefan_UPP' + '_alignmentFasta.fasta'
    #predictionUPP = './alignData/UPPoutput/output_alignment.fasta'

    keys1, seqs1 = loadFastaBasic(predictionStefan)
    keys2, seqs2 = loadFastaBasic(predictionUPP)

    keys2, seqs2 = keys2[np.isin(keys2, keys1)], seqs2[np.isin(keys2, keys1)]

    print (keys1[0])
    print (keys2[0])
    print (seqs1[0])
    print ('')
    print (seqs2[0])



#compareToUPP()



#quit()

#print ("Align Length: ", len(dataAlign1), ', ', len(dataAlign2))

#os.system("java -jar ./FastSP.jar -r " + trueAlignName + " -e " + predictionName_fasta)




#quit()

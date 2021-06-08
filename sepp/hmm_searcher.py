import numpy as np
import matplotlib.pyplot as plt
from sepp.hmm import *

def hierchySearch(fakeSimulate=False):

    earlyStop = False
    observeNephew = False

    scores_original = np.load("./alignData/hmmScores/full.npy")

    hmmNames = []
    sequenceFileNames = giveSequenceFileNames()
    for a in range(0, len(sequenceFileNames)):
        hmmName = "./alignData/initialHMM/test" + str(a) + ".hmm"
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
        scores = saveScoreFromBool(hmmNames, queryHMM, scoreName, noSave=True)
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
                scores = saveScoreFromBool(hmmNames, queryHMM, scoreName, noSave=True)

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

    np.save('./alignData/Searcher/scoreFiles/score.npy', scoresFull)
    np.save('./alignData/Searcher/scoreFiles/Bool.npy', boolDone)

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


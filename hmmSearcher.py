import numpy as np
import matplotlib.pyplot as plt
from stefanHMM import *

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
    if False:
        argsInternal = np.argwhere(treeSum != 1)[:, 0]
        treeDataLeaf = np.copy(treeData)
        treeDataLeaf[:, argsInternal] = 0
        treeSumLeaf = np.sum(treeDataLeaf, axis=1)
        heightShift = np.log2(treeSumLeaf) * 0.01 * 364.2
        for a in range(len(scoreAdjustment[0])):
            scoreAdjustment[:, a] = heightShift[a]


    #initialSearchSet = [0, children[0, 0], children[0, 1], children[children[0, 0], 0], children[children[0, 1], 0], children[children[0, 0], 1], children[children[0, 1], 1]]
    #initialSearchSet = sortedSize[-31:]
    #initialSearchSet = sortedSize[-7:]
    initialSearchSet = np.array([0])

    queryPart = np.arange(Nquery * len(initialSearchSet)) % Nquery
    hmmPart = np.array(initialSearchSet).repeat(Nquery)

    queryHMM = np.array([queryPart, hmmPart]).T
    #queryHMM = np.array([np.arange(Nquery), np.zeros(Nquery)]).T.astype(int)
    #scoreName = './alignData/Searcher/scoreFiles/0.npy'
    scoreName = ''
    if fakeSimulate:
        scores = np.copy(scores_original)
    else:
        scores = saveScoreFromBool(hmmNames, queryHMM, scoreName, noSave=True)
    scoresFull[queryHMM[:, 0], queryHMM[:, 1]] = scores[queryHMM[:, 0], queryHMM[:, 1]]

    queryHMM_original = np.copy(queryHMM)

    #plt.show()

    #maxHMM = np.zeros(Nquery).astype(int)
    maxHMM = np.argmax(scoresFull, axis=1)
    boolDone = np.zeros((Nquery, Nhmm))
    boolDone[:, 0] = 1
    #boolDone[queryHMM[:, 0], queryHMM[:, 1]] = 1

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


        #queryCorrespond = np.array([maxNotLeaf, maxNotLeaf]).T

        #print (np.max(maxNotLeaf))

        maxHMMChild = np.copy(maxHMMChild_original).reshape((maxHMMChild_original.size,))
        #queryCorrespond = queryCorrespond.reshape((maxHMMChild.size,))
        queryCorrespond = maxNotLeaf[np.arange(maxHMMChild.size) // (maxHMMChild_original.shape[1]) ]

        queryToHMM_try = np.array([queryCorrespond, maxHMMChild]).T


        queryHMM = queryToHMM_try[boolDone[queryToHMM_try[:, 0], queryToHMM_try[:, 1] ] == 0]
        #assert queryToHMM_try.shape[0] == queryHMM.shape[0]

        boolDone[queryHMM[:, 0], queryHMM[:, 1] ] = 1


        if queryHMM.shape[0] != 0:

            scoreName = ''

            if fakeSimulate:
                scores = np.copy(scores_original)
            else:
                scores = saveScoreFromBool(hmmNames, queryHMM, scoreName, noSave=True)

                #scores = saveScoreFromBool(hmmNames, np.concatenate((queryHMM, queryHMM_original), axis=0), scoreName, noSave=True)


            #scoresFull[queryHMM[:, 0], queryHMM[:, 1]] = scores[queryHMM[:, 0], queryHMM[:, 1]]

            scoresFull[queryHMM[:, 0], queryHMM[:, 1]] = np.copy(scores[queryHMM[:, 0], queryHMM[:, 1]])




            newScores1 = scoresFull[maxNotLeaf, maxHMMChild_original[:, 0]]
            newScores2 = scoresFull[maxNotLeaf, maxHMMChild_original[:, 1]]
            newScoresMax = np.argmax(np.array([newScores1, newScores2]), axis=0)

            #if (np.min(newScores1) < -100) or (np.min(newScores2) < -100):
            #    print (np.min(newScores1))
            #    print (np.min(newScores2))
            #    quit()

            newChild = maxHMMChild_original[np.arange(newScores1.shape[0]), newScoresMax]

            if earlyStop:
                maxHMM = np.argmax(scoresFull, axis=1)
            else:
                maxHMM[maxNotLeaf] = newChild

            scoreStep = scoresFull[np.arange(maxHMM.shape[0]), maxHMM]
            scoresRecord.append(np.copy(scoreStep))

            maxHMMlist.append(np.copy(maxHMM[:10]))

            #print (np.mean(treeSum[maxHMM]))

        else:
            done = True



        a += 1



    scoresRecord = np.array(scoresRecord)


    #np.save('./alignData/Searcher/scoreFiles/Early_score.npy', scoresFull)
    #np.save('./alignData/Searcher/scoreFiles/Early_Bool.npy', boolDone)
    np.save('./alignData/Searcher/scoreFiles/score.npy', scoresFull)
    np.save('./alignData/Searcher/scoreFiles/Bool.npy', boolDone)
    #treeData = np.load('./alignData/treeDecomp.npy')

    maxHMM_guess = np.argmax(scoresFull, axis=1)
    maxHMM_original = np.argmax(scores_original, axis=1)

    argDiff = np.argwhere((maxHMM_original-maxHMM_guess) != 0)[:, 0]

    argDiffTop = argDiff[treeSum[maxHMM_guess[argDiff]] == 999]

    print (np.mean(scores_original[argDiff, maxHMM_original[argDiff]] - scores_original[argDiff, maxHMM_guess[argDiff]]))

    #plt.plot(np.mean(scoresRecord[:, :], axis=1))
    #plt.show()

    #plt.plot(np.mean(scoresRecord[:, argDiffTop], axis=1))
    #plt.show()

    print (argDiff.shape)
    print (np.unique(treeSum[maxHMM_guess[argDiff]], return_counts=True)[1][-3:])
    print (np.unique(treeSum[maxHMM_guess], return_counts=True)[1][-3:])

    #print (np.mean(scores_original[argDiff, maxHMM_original[argDiff]]))
    #print (np.mean(scores_original[np.arange(maxHMM_original.shape[0]), maxHMM_original]))

    #plt.hist(scores_original[np.arange(Nquery), maxHMM_original], bins=100)
    #plt.show()

    #plt.hist(scores_original[argDiff, maxHMM_original[argDiff]], bins=100)
    #plt.show()

    bitscore_guess = scores_original[np.arange(Nquery), maxHMM_guess]
    bitscore_original = scores_original[np.arange(Nquery), maxHMM_original]

    #plt.hist((bitscore_original - bitscore_guess)[argDiff], bins=100)
    #plt.show()


    #queryTest = np.arange(Nquery * 3) // Nquery
    #HMMList = np.concatenate((np.zeros(Nquery), np.zeros(Nquery) + children[0, 0], np.zeros(Nquery) + children[0, 1]))


# hierchySearch(fakeSimulate=True)
# quit()



'''
scoresFull = np.load('./alignData/Searcher/scoreFiles/score.npy')
scores = np.load("./alignData/hmmScores/fullAdjusted.npy")
maxHMM_guess = np.argmax(scoresFull, axis=1)
maxHMM_original = np.argmax(scores, axis=1)

treeData = findDecomposition()
treeSum = np.sum(treeData, axis=1)
treeSumRepeat = treeSum.repeat(treeSum.shape[0]).reshape(treeData.shape)
treeDataWithSum = treeData * (10000 - treeSumRepeat)
treeDataWithSum_T = treeDataWithSum.T

multSum = treeDataWithSum_T[maxHMM_guess] * treeDataWithSum_T[maxHMM_original]
LCA = np.argmax(multSum, axis=1)

#plt.plot(treeSum[LCA])

argsDiff = np.argwhere((maxHMM_guess - maxHMM_original) != 0)[:, 0]

plt.plot(treeSum[maxHMM_guess[argsDiff]])
plt.plot(treeSum[maxHMM_original[argsDiff]])
plt.show()





quit()

plt.plot(maxHMM)
plt.plot(maxHMMGuess)
plt.show()
'''





#saveScoreFromBool

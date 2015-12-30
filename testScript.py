import numpy as np
import os

import math
import pickle

from sklearn import linear_model,cross_validation
from sklearn.cross_validation import KFold
from sklearn import grid_search

from sklearn.linear_model import SGDClassifier

from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import f1_score

from sklearn.metrics import confusion_matrix

supportVectors = [False,False,False,True,False,False,True,False,False,False,False,True,False,False,False,False,False,False,False,False,False,True,False,True,False,True,True,False,True,False,False,True,True,False,False,False,False,False,False,False,True,False,False,False,False,True,True,False,False,True,True,True,False,False,False,False,True,True,False,True,True,True,False,False,False,True,False,True,False,True,True,True,False,False,False,False,False,False,False,True,True,True,True,True,True,True,True,False,True,True,True,True,False,False,False,True,False,True,True,True,True,True,False,False,False,True,True,False,False,False,False,True,False,True,False,True,True,False,True,True,True,False,False,False,False,False,False,True,False,True,True,True,False,False,False,False,False,True,True,True,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,True,True,False,False,False,False,True,False,False,False,False,False,True,False,False,False,False,False,False,True,True,False,True,True,True,True,False,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,False,True,False,False,False,False,False,False,False,False,False,False,False,True ]

def sigmoidfn(x):
    return 1/(1+math.exp(-x))

#
foldFileName = 'prot_list'

foldFile = open(foldFileName, 'r+')

f = open('lr_classifier.pickle')
classifier = pickle.load(f)
f.close()

# For trying with the SVM Classifier
# f = open('svm_classifier.pickle')
# svm_classifier = pickle.load(f)
# f.close()

f = open('lrfpr.pickle')
lrfpr = pickle.load(f)
f.close()


f = open('selectfpr.pickle')
selectfpr = pickle.load(f)
f.close()

i11 = 0       # do for 2-3 proteins initially to check

accuracyArr = []
mccArr = []
fScoreArr = []

for proteinName in foldFile.readlines():
    
    # Flag to account for any processing errors while building the feature vectors
    flag = 0
    
    # Load PRSA file for each protein 
    prsaFileName = './PRSA/'+proteinName[:-1]+'.prsa'
    prsaFileArr = np.loadtxt(prsaFileName, usecols=(3,))
    # Load ss2 file
    ss2FileName = './psipred/'+proteinName[:-1]+'.ss2'
    ss2FileArr = np.loadtxt(ss2FileName, skiprows=2, usecols=(0,3,4,5), ndmin=2)
    # Load pssm file
    pssmFileName = './PSSM/'+proteinName[:-1]+'.pssm'
    pssmFile = open(pssmFileName, 'r+')
    pssmList = []
    
    i = 0
    for lines in pssmFile.readlines():
        if i<3:
            pass
        else:
            try:
                pssmList.append([int(lines[i:i+3]) for i in range(9,len(lines)-94,3)])     #162-94
            except:
                #print proteinName, "Check line ", i
                flag = 1
                break

        i+=1

    # Just accounting for the errors in processing in pssm file
    if flag == 1:
        continue
    else:
        pass

    # Load PPI-site File for final comparison after prediction
    # Only ppiFileArr[4] through to ppiFileArr[-4] is used for each protein because negihbouring 4 acids' values are taken per amino acid
    ppiFileName = './PPI-site/'+proteinName[:-1]+'.int'
    ppiFileArr = np.loadtxt(ppiFileName,usecols=(0,))

    yActualArr = ppiFileArr[4:-4]         #Will be used for final testing

    X = []  # Set of feature vectors to be fed into the classifer at one go for one set of values


    j = 4       # Start with the 4th position for each protein and move on to the L-4th
                # j is the position of each amino acid within the protein

    while j<len(ppiFileArr)-4:
        # 1. Get 9 PSSM value sets - each is a set of 20 values => 9*20 = 180
        
        x = []  # Per-amino-acid feature vector

        i = j-4
        while i<j+5:
            try:
                for n in pssmList[i]:
                    x.append(sigmoidfn(n))

            except:
                print "some error while getting pssm values"
            i+=1

        # 2. Get 9 ss2 value-sets (each is a set of 3 values)      => 27 values
        i = j-4
        while i<j+5:
            x.append(ss2FileArr[i][1])
            x.append(ss2FileArr[i][2])
            x.append(ss2FileArr[i][3])

            i+=1 

        # 3. get 9 prsa values
        i = j-4
        while i<j+5:
            x.append(prsaFileArr[i])
            i+=1

        X.append(x)         # Readying each feature vector
        j+=1
   

    del pssmList[:]
    np.delete(ss2FileArr, [i for i in range(len(ss2FileArr))])
    np.delete(prsaFileArr, [i for i in range(len(prsaFileArr))])
    np.delete(ppiFileArr, [i for i in range(len(ppiFileArr))])
    pssmFile.close()
    
    i11 += 1

    # Trying lrfpr w/ selectfpr
    xnew2 = []
    # supportVectors = selectfpr.get_support()
    # print len(supportVectors)
    for jq in range(len(X)):
        xq = []
        for iq in range(len(supportVectors)):
            if supportVectors[iq] == True:
                xq.append(X[jq][iq])
        xnew2.append(xq)
    # print len(xnew2), len(xnew2[0])

    yPredArrFPR = lrfpr.predict(xnew2[:])
    print proteinName
    print "Accuracy: ", lrfpr.score(xnew2[:],yActualArr)
    print "MCC: ", matthews_corrcoef(yActualArr,yPredArrFPR)
    print "f-score: ",f1_score(yActualArr,yPredArrFPR)
    print confusion_matrix(yActualArr,yPredArrFPR)

    # yPredArr = []
    # yPredArr = classifier.predict(X[:])
    # print proteinName,
    # accuracyArr.append(classifier.score(X[:],yActualArr) )
    # mccArr.append(matthews_corrcoef(yActualArr,yPredArr))
    # fScoreArr.append(f1_score(yActualArr,yPredArr))
    # print "Accuracy:",accuracyArr[-1]
    # print "\tMCC: ",mccArr[-1]
    # print "\tF-measure: ",fScoreArr[-1]
    # print

    del X[:]
    # np.delete(yPredArr, [i for i in range(len(yPredArr))])
    np.delete(yActualArr, [i for i in range(len(yActualArr))])

# print "Total Average Accuracy for test-72: ", float(sum(accuracyArr))/len(accuracyArr)
# print "Total Average MCC Score for test-72: ", float(sum(mccArr))/len(mccArr)
# print "Total Average f-score for test 72: ", float(sum(fScoreArr))/len(fScoreArr)
# print 
# print "Max Accuracy obtained for a chain in test-72: ", max(accuracyArr)
# print "Max MCC obtained for a chain in test-72: ", max(mccArr)
# print "Max f-score obtained for a chain in test-72: ", max(fScoreArr)
# print 
# print "Min Accuracy obtained for a chain in test-72: ", min(accuracyArr)
# print "Min MCC obtained for a chain in test-72: ", min(mccArr)
# print "Min f-score obtained for a chain in test-72: ", min(fScoreArr)    

# posmcc = 0
# posmccSum = 0
# negmcc = 0
# for i in range(len(mccArr)):
#     if mccArr[i]>=0:
#         posmcc+=1
#         posmccSum+=mccArr[i]
#     else:
#         negmcc+=1

# print "Number of positive mcc values: ", posmcc
# print "Avg of only pos mcc values: ", float(posmccSum)/posmcc
# print "Number of negative mcc values: ", negmcc

    # Remove
    # if i11>5:   
    #     break

    
"""
Final Comments:

Unable to process the PSSM file for 2 protein chains out of the 60 proteins 1pxv_A and 1pxv_C.


A high number of negative examples in test set is giving - bad results


"""


"""
# xnew1 = SelectFdr(f_classif).fit_transform(x,y)
# xnew2 = SelectFpr(f_classif).fit_transform(x,y)
# print len(xnew1[0]), len(xnew2[0])


# ypredFdr = linear_model.LogisticRegression().fit(xnew1[:8000],y[:8000]).predict(xnew1[8000:])
# # print ypredFdr
# print linear_model.LogisticRegression().fit(xnew1[:8000],y[:8000]).score(xnew1[8000:],y[8000:])
# print "MCC: ", matthews_corrcoef(y[8000:],ypredFdr)
# print "f1_score: ", f1_score(y[8000:], ypredFdr)

# ypredFpr = linear_model.LogisticRegression().fit(xnew2[:8000],y[:8000]).predict(xnew2[8000:])
# print linear_model.LogisticRegression().fit(xnew2[:8000],y[:8000]).score(xnew2[8000:],y[8000:])
# print "MCC: ", matthews_corrcoef(y[8000:],ypredFpr)
# print "f1_score: ", f1_score(y[8000:], ypredFpr)
"""
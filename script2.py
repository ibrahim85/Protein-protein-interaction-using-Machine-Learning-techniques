import numpy as np
import os

#
foldFileName = 'IDfold10.INTlist'

foldFile = open(foldFileName, 'r+')

# for fil in os.listdir('.'):
#     if fil[:6] == 'IDfold':
#         foldFileName = fil
#         foldFile = open(foldFileName,'r+')

#         #These lists contain the entire set of protein chains and their corresponding positively interacting
#         #residues and negatively interacting residues
proteinList = []
pos_residuesList = []
neg_residuesList = []

for lines in foldFile.readlines():
    protein, pos_residues, neg_residues = lines.split()
    proteinList.append(protein)
    pos_residuesList.append(map(int,pos_residues.split(',')))
    neg_residuesList.append(map(int,neg_residues.split(',')))



for k in range(len(proteinList)):

    # print pos_residuesList[k]
    # print neg_residuesList[k]
    # print proteinList[k]

    ss2FileName = './psipred/'+proteinList[k]+'.ss2'
    pssmFileName = './PSSM/'+proteinList[k]+'.pssm'
    backbonePredFileName = './DynaMine/'+proteinList[k]+'_backbone.pred'

    ss2File = open(ss2FileName,'r+')
    pssmFile = open(pssmFileName, 'r+')
    backbonePredFile = open(backbonePredFileName,'r+')

    pssmList = []

    i = 0
    for lines in pssmFile.readlines():
        if i<3:
            pass
        else:
            try:
                pssmList.append([int(lines[i:i+3]) for i in range(9,len(lines)-94,3)])     #162-94
            except:
                print "check this line"
                print proteinList[k]
                print i
        i+=1


    ss2FileArr = np.loadtxt(ss2FileName, skiprows=2, usecols=(0,3,4,5), ndmin=2)

    
    #print ss2FileArr[105][1]

    backbonePredArr = np.loadtxt(backbonePredFileName,usecols=(1,))


    #Open file to write to
    fvFile = foldFileName[2:7]+'0.fv'       #0.fv only for 10

    #Write the feature vector to the file

    fold1fv = open(fvFile,'a')
    

    for j in pos_residuesList[k]:
        fold1fv.write('1\t')

        i = j-4
        while i<j+5:
            try:

                for n in pssmList[i]:
                    fold1fv.write('%d' %n)
                    fold1fv.write('\t')

            except:
                print "check inside for loop"
                print 'k value: ', k
                print proteinList[k]
                print len(pos_residuesList[k])
            i+=1

        i = j-4
        while i<j+5:
            fold1fv.write('%f' %ss2FileArr[i][1])
            fold1fv.write('\t')
            fold1fv.write('%f' %ss2FileArr[i][2])
            fold1fv.write('\t')
            fold1fv.write('%f' %ss2FileArr[i][3])
            fold1fv.write('\t')
            i+=1

        i = j-4
        while i<j+5:
            fold1fv.write('%f' %backbonePredArr[i])
            fold1fv.write('\t')
            i+=1

        fold1fv.write('\n')

    for j in neg_residuesList[k]:
        fold1fv.write('-1\t')

        i = j-4

        while i<j+5:
            try:
                for n in pssmList[i]:
                    fold1fv.write('%d' %n)
                    fold1fv.write('\t')
            except:
                print "check inside for loop"
                #print pssmList[i]
                print 'k value: ', k
                print proteinList[k]
                print len(pos_residuesList[k])
            i+=1

        i = j-4
        while i<j+5:
            fold1fv.write('%f' %ss2FileArr[i][1])
            fold1fv.write('\t')
            fold1fv.write('%f' %ss2FileArr[i][2])
            fold1fv.write('\t')
            fold1fv.write('%f' %ss2FileArr[i][3])
            fold1fv.write('\t')
            i+=1

        i = j-4
        while i<j+5:
            fold1fv.write('%f' %backbonePredArr[i])
            fold1fv.write('\t')
            i+=1

        fold1fv.write('\n')

    del pssmList[:]
    np.delete(ss2FileArr, [i for i in range(len(ss2FileArr))])
    np.delete(backbonePredArr, [i for i in range(len(backbonePredArr))])

    fold1fv.close()

"""

i = 0
for lines in ss2File.readlines():
    if i<2:
        pass
    else:
        pos_1_indexed, alpha1,alpha2, v1,v2,v3= lines.split()

        pos_1_indexedList.append(pos_1_indexed)
        v1List.append(v1)
        v2List.append(v2)
        v3List.append(v3)
    i+=1
print pos_1_indexedList[pos_residuesList[0][0]], v1List[pos_residuesList[0][0]], v2List[pos_residuesList[0][0]],v3List[pos_residuesList[0][0]]

"""



"""
ALTERNATE using numpy:
str = pssmFile.readline()
[str[i:i+3] for i in range(0,len(str),3)]

['   ', ' 4 ', 'P  ', ' -3', ' -6', ' -6', ' -6', ' -7', ' -6', ' -5', ' -7', ' -7', ' -7', ' -7', ' -3', ' -7', ' -8', '  9', ' -3', ' -5', ' -8', ' -7', ' -7', '   ', ' 2 ', '  0', '   ', '0  ', ' 0 ', '  0', '   ', '0  ', ' 0 ', '  0', '   ', '0  ', ' 0 ', '  0', '   ', '2  ', ' 0 ', '  0', '  9', '5  ', ' 2 ', '  0', '   ', '0  ', ' 0 ', '  0', '  3', '.16', ' 0.', '94\n']


str = pssmFile.readline()       #Ignore first 3 lines
[str[i:i+3] for i in range(9,len(str)-94,3)]

"""

"""
# Doing only for the first protein



ss2FileName = './psipred/'+proteinList[0]+'.ss2'
pssmFileName = './PSSM/'+proteinList[0]+'.pssm'
backbonePredFileName = './DynaMine/'+proteinList[0]+'_backbone.pred'

ss2File = open(ss2FileName,'r+')
pssmFile = open(pssmFileName, 'r+')
backbonePredFile = open(backbonePredFileName,'r+')

pos_1_indexedList = []
v1List = []
v2List = []
v3List = []

pssmList = []

i = 0
for lines in pssmFile.readlines():
    if i<3:
        pass
    else:
        #str1=pssmFile.readline()
        #print str1
        pssmList.append([int(lines[i:i+3]) for i in range(9,len(lines)-94,3)])     #162-94=
    i+=1


print pssmList[0]
print pssmList[32]
print pssmList[2]

ss2FileArr = np.loadtxt(ss2FileName, skiprows=2, usecols=(0,3,4,5), ndmin=2)

#the three required values at position 32
print ss2FileArr[32][1]
print ss2FileArr[32][2]
print ss2FileArr[32][3]


backbonePredArr = np.loadtxt(backbonePredFileName,usecols=(1,))

print backbonePredArr[32]



#Open file to write to
#Write the feature vector to the file

fold1fv = open('fold1.fv','a')

for j in pos_residuesList[0]:
    fold1fv.write('1\t')

    i = int(j)
    i = i-4
    while i<j+5:
        for n in pssmList[i]:
            fold1fv.write('%d' %n)
            fold1fv.write('\t')
        i+=1

    i = j-4
    while i<j+5:
        fold1fv.write('%f' %ss2FileArr[i][1])
        fold1fv.write('\t')
        fold1fv.write('%f' %ss2FileArr[i][2])
        fold1fv.write('\t')
        fold1fv.write('%f' %ss2FileArr[i][3])
        fold1fv.write('\t')
        i+=1

    i = j-4
    while i<j+5:
        fold1fv.write('%f' %backbonePredArr[i])
        fold1fv.write('\t')
        i+=1

    fold1fv.write('\n')

for j in neg_residuesList[0]:
    fold1fv.write('-1\t')

    i = int(j)
    i = i-4

    while i<j+5:
        for n in pssmList[i]:
            fold1fv.write('%d' %n)
            fold1fv.write('\t')
        i+=1

    i = j-4
    while i<j+5:
        fold1fv.write('%f' %ss2FileArr[i][1])
        fold1fv.write('\t')
        fold1fv.write('%f' %ss2FileArr[i][2])
        fold1fv.write('\t')
        fold1fv.write('%f' %ss2FileArr[i][3])
        fold1fv.write('\t')
        i+=1

    i = j-4
    while i<j+5:
        fold1fv.write('%f' %backbonePredArr[i])
        fold1fv.write('\t')
        i+=1

    fold1fv.write('\n')

fold1fv.close()

"""
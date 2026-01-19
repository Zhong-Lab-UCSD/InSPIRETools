#read in mapped read pairs n check protein-coding genes
from collections import defaultdict
import sys
import glob
import random

bcDirPath=sys.argv[1]
geneDicPath=sys.argv[2]
outDirPath=sys.argv[3]
jobId=sys.argv[4]

#read in libOrder
libIdList=set()
fileList=glob.glob('%s/*'%(bcDirPath))
dicId_order={}
libOrderList=[]
for file in fileList:
    libOrder=file[-6:-4]
    libOrderList.append(libOrder)
    with open(file,'r') as f:
        next(f)
        for line in f:
            splitLine=line.strip().split('\t')
            libId=splitLine[-1]
            libIdList.add(libId)
            dicId_order[libId]=libOrder
            break
            

#read in refseq dic
dicIdGeneName={}
dicIdGeneType={}
with open('%s'%(geneDicPath),'r') as f:
    for line in f:
        splitLine=line.strip().split(',')
        dicIdGeneName[splitLine[0]]=splitLine[1]
        dicIdGeneType[splitLine[0]]=splitLine[2]
        
#read readId-molecule info
dicDrugInfo_mole={}
for libOrder in libOrderList:
    with open('%s/BB-Codon Map_HGP0001-OpenDEL00%s.txt'%(bcDirPath,libOrder),'r') as f:
        next(f)
        next(f)
        for line in f:
            splitLine=line.strip().split('\t')
            cycle=splitLine[1]
            drugMole=splitLine[2]
            drugOrder=splitLine[0]
            keyInfo=';'.join([libOrder,cycle,drugOrder])
            dicDrugInfo_mole[keyInfo]=drugMole

#read id umi
dicId_umi={}
fileList=glob.glob('%s/%sreadId_UMI/*'%(outDirPath,jobId))
for file in fileList:
    with open(file,'r') as f:
        for line in f:
            splitLine=line.strip().split(',')
            while len(splitLine[1])<8:
                ha=random.sample(['A','T','G','C'],1)
                splitLine[1]+=ha[0]
            dicId_umi[splitLine[0]]=splitLine[1]
            
#read end info
comboSet=set()
dicId_durgInfo={}
with open('%s/%sintermediateFiles/endInfo.csv'%(outDirPath,jobId),'r') as f:
    next(f)
    for line in f:
        splitLine=line.strip().split(',')
        readId,lib,combo=splitLine[0],splitLine[-2],splitLine[-1]
        umi=dicId_umi[readId]
        finalCombo=';'.join([lib]+[combo]+[umi])
        dicId_durgInfo[readId]=finalCombo
        
#general stats
statsFile=open('%s/generalStats.csv'%(outDirPath),'w')
libReadCount=0
with open('%s/libIds_stats.csv'%(outDirPath),'r') as f:
    for line in f:
        splitLine=line.strip().split(',')
        libReadCount+=int(splitLine[1])
        
count1=0
count2=0
count3=0
with open('%s/drugIds_stats.csv'%(outDirPath),'r') as f:
    for line in f:
        splitLine=line.strip().split(',')
        count1+=int(splitLine[1])
        count2+=int(splitLine[2])
        count3+=int(splitLine[3])
statsFile.write('Reads with libID,%d\n'%(libReadCount))
statsFile.write('Reads with C1ID,%d\n'%(count1))
statsFile.write('Reads with C2ID,%d\n'%(count2))
statsFile.write('Reads with C3ID,%d\n'%(count3))
            
finalInfoSet=set()
validCount=0
roughSet1=set()
readId_set1=set()
dupCount=0

targetFile=open('%sSmoProteinAssociations.csv'%(outDirPath),'w')
targetFile.write(
    'readId,drugLibrary,cycle1_index,cycle2_index,cycle3_index,cycle1_molecule,cycle2_molecule,cycle3_molecule,protein\n')

with open('%s/%salignment/mapped.sorted.bed'%(outDirPath,jobId),'r') as f:
    for line in f:
        splitLine=line.strip().split('\t')
        readId=splitLine[3]
        txId=splitLine[0]
        cigar=splitLine[6]
        start,end=splitLine[1],splitLine[2]
        roughSet1.add(readId)
        if dicIdGeneType[txId]=='protein_coding' and readId not in readId_set1:
            readId_set1.add(readId)
            mapInfo=';'.join([txId,start,end])
            drugInfo=dicId_durgInfo[readId]
            finalInfo=drugInfo+mapInfo
            if finalInfo not in finalInfoSet:
                validCount+=1
                protein=dicIdGeneName[txId]
                [lib,combo,umi]=drugInfo.split(';')
                [c1Index,c2Index,c3Index]=combo.split('-')
                key1=';'.join([lib,'1',c1Index])
                key2=';'.join([lib,'2',c2Index])
                key3=';'.join([lib,'3',c3Index])
                c1Mole,c2Mole,c3Mole=dicDrugInfo_mole[key1],dicDrugInfo_mole[key2],dicDrugInfo_mole[key3]
                infoList=','.join([readId,lib,c1Index,c2Index,c3Index,c1Mole,c2Mole,c3Mole,protein])
                targetFile.write(infoList)
                targetFile.write('\n')
            else:
                dupCount+=1
            finalInfoSet.add(finalInfo)
targetFile.close()

statsFile.write('Reads mapped to tx,%d\n'%(len(roughSet1)))
statsFile.write('Reads mapped to protein,%d\n'%(len(readId_set1)))
statsFile.write('dedpued DPI rads,%d\n'%(validCount))
statsFile.close()
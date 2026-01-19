import sys
import glob
from collections import defaultdict
import numpy as np

#help functions
def neighbors(pattern, d):
    if d == 0:
        return(pattern)
    if len(pattern) == 1:
        return(['A','T','C','G'])
    nlist = []
    suffixn = neighbors(SUFFIX(pattern),d)
    for string in suffixn:
        if hamming(string, SUFFIX(pattern)) < d:
            for x in ['A','T','C','G']:
                nlist.append(x+string)
        else:
            nlist.append(pattern[0]+string)
    return(nlist)

def SUFFIX(pattern):
    return(pattern[1:])

def hamming(str1, str2):
    count = 0 
    assert (len(str1)==len(str2)), "length doesn't match"
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            count += 1
    return(count)

def n_mis_neighbors(words_dict, d):
  '''
  words dict with keys the reference barcode and values the barcode name
  '''
  hash_ = {}
  words = list(words_dict.keys())
  for i in words:
    word = i
    name = words_dict[i]
    # remaining bcd string
    bcd_nei = neighbors(word, d)
    for b in bcd_nei:
      hash_[b] = name
  return hash_

bcDirPath=sys.argv[1]
r1Path=sys.argv[2]
r2Path=sys.argv[3]
outDirPath=sys.argv[4]
jobId=sys.argv[5]

#read in library barcodes
fileList=glob.glob('%s/*'%(bcDirPath))
dicBC_id={}
libIdSet=set()
for file in fileList:
    libId=file[-6:-4]
    with open(file,'r') as f:
        next(f)
        for line in f:
            splitLine=line.strip().split('\t')
            libBC=splitLine[-1]
            dicBC_id[libBC]=libId
            libIdSet.add(libId)
            break
dicLibBC_hash = n_mis_neighbors(dicBC_id, 1)

#check library barcodes
i=0
temp=[]
flag=False
diclibId_fastq=defaultdict(list)
diclibId_umi=defaultdict(list)
diclibId_readSet=defaultdict(set)
#R1
with open('%s'%(r1Path),'r') as f:
        for line in f:
            temp.append(line)
            i+=1
            if i==1:
                splitLine=line.strip().split()
                readId=splitLine[0][1:]
            if i==2:
                seq=line[:-1]
                for j in range(39,len(seq)-12-9):
                    tempSeq=seq[j:j+9]
                    libId=dicLibBC_hash.get(tempSeq, 'NA')
                    if libId!='NA':
                        newSeq=seq[:j]
                        umi=seq[j+9:j+len(libId)+12]
                        diclibId_umi[libId].append([readId,umi])
                        diclibId_readSet[libId].add(readId)
                        temp[1]=newSeq+'\n'
                        flag=True
                        break
            if i==4:
                if flag:
                    diclibId_fastq[libId].append(temp)
                i=0
                temp=[]
                flag=False
#write into file
for libId in libIdSet:
    targetFile=open('%s/%sprocessedFastq/R1.lib%s.fastq'%(outDirPath,jobId,libId),'w')
    ha=diclibId_fastq[libId]
    for chunk in ha:
        for line in chunk:
            targetFile.write(line)
    targetFile.close()
    
#R2
i=0
temp=[]
flag=False
diclibId_fastq=defaultdict(list)
with open('%s'%(r2Path),'r') as f:
        for line in f:
            temp.append(line)
            i+=1
            if i==1:
                splitLine=line.strip().split()
                readId=splitLine[0][1:]
            if i==2:
                seq=line[:-1]
                for j in range(39,len(seq)-12-9):
                    tempSeq=seq[j:j+9]
                    libId=dicLibBC_hash.get(tempSeq, 'NA')
                    if libId!='NA':
                        newSeq=seq[:j]
                        temp[1]=newSeq+'\n'
                        umi=seq[j+9:j+len(libId)+12]
                        diclibId_umi[libId].append([readId,umi])
                        diclibId_readSet[libId].add(readId)
                        flag=True
                        break
            if i==4:
                if flag:
                    diclibId_fastq[libId].append(temp)
                i=0
                temp=[]
                flag=False
#write into file
statsFile=open('%s/libIds_stats.csv'%(outDirPath),'w')
for libId in libIdSet:
    targetFile=open('%s/%sprocessedFastq/R2.lib%s.fastq'%(outDirPath,jobId,libId),'w')
    ha=diclibId_fastq[libId]
    for chunk in ha:
        for line in chunk:
            targetFile.write(line)
    targetFile.close()
    #write umis
    umiFile=open('%s/%sreadId_UMI/lib%s.csv'%(outDirPath,jobId,libId),'w')
    for ha in diclibId_umi[libId]:
        [readId,umi]=ha
        umiFile.write('%s,%s\n'%(readId,umi))
    umiFile.close()
    readCount=len(diclibId_readSet[libId])
    statsFile.write('%s,%d\n'%(libId,readCount))    
statsFile.close()    

#check drug barcodes
statsFile=open('%s/drugIds_stats.csv'%(outDirPath),'w')
infoFile=open('%s/%sintermediateFiles/endInfo.csv'%(outDirPath,jobId),'w')
infoFile.write('readId,drugEnd,proteinEnd,libIndex,drugIndexCombo\n')
for libId in libIdSet:
    #read in drug barcodes by cycle
    drugBC_id_c1={}
    drugBC_id_c2={}
    drugBC_id_c3={}
    with open('%s/BB-Codon Map_HGP0001-OpenDEL00%s.txt'%(bcDirPath,libId),'r') as f:
        next(f)
        next(f)
        for line in f:
            splitLine=line.strip().split('\t')
            cycle=int(splitLine[1])
            drugBC=splitLine[3]
            drugId=splitLine[0]
            if cycle==1:
                drugBC_id_c1[drugBC]=drugId
            if cycle==2:
                drugBC_id_c2[drugBC]=drugId
            if cycle==3:
                drugBC_id_c3[drugBC]=drugId
    dicDrugBC_hash_c1 = n_mis_neighbors(drugBC_id_c1, 1)
    dicDrugBC_hash_c2 = n_mis_neighbors(drugBC_id_c2, 1)
    dicDrugBC_hash_c3 = n_mis_neighbors(drugBC_id_c3, 1)
    
    #read in R1 fastq file
    i=0
    count1=0
    count2=0
    count3=0
    with open('%s/%sprocessedFastq/R1.lib%s.fastq'%(outDirPath,jobId,libId),'r') as f:
        for line in f:
            i+=1
            if i==1:
                splitLine=line.strip().split()
                readId=splitLine[0][1:]
            if i==2:
                seq=line[:-1]
                #check cycle 1 drug
                for j in range(len(seq)-33):
                    tempSeq=seq[j:j+11]
                    drugId1=dicDrugBC_hash_c1.get(tempSeq, 'NA')
                    if drugId1!='NA':
                        count1+=1
                        newSeq=seq[j+11:]
                        #check cycle2 drug
                        for k in range(len(newSeq)-22):
                            tempSeq=newSeq[k:k+11]
                            drugId2=dicDrugBC_hash_c2.get(tempSeq, 'NA')
                            if drugId2!='NA':
                                count2+=1
                                newSeq2=newSeq[k+11:]
                                #check cycle3 drug
                                for l in range(len(newSeq2)-11):
                                    tempSeq=newSeq2[l:l+11]
                                    drugId3=dicDrugBC_hash_c3.get(tempSeq, 'NA')
                                    if drugId3!='NA':
                                        count3+=1                                    
                                        infoFile.write('%s,%d,%s,%s-%s-%s\n'%(readId,1,libId,drugId1,drugId2,drugId3))
                                        break                    
            if i==4:
                i=0
    
    #read in R2 fastq file
    i=0
    with open('%s/%sprocessedFastq/R2.lib%s.fastq'%(outDirPath,jobId,libId),'r') as f:
        for line in f:
            i+=1
            if i==1:
                splitLine=line.strip().split()
                readId=splitLine[0][1:]
            if i==2:
                seq=line[:-1]
                #check cycle 1 drug
                for j in range(len(seq)-33):
                    tempSeq=seq[j:j+11]
                    drugId1=dicDrugBC_hash_c1.get(tempSeq, 'NA')
                    if drugId1!='NA':
                        count1+=1
                        newSeq=seq[j+11:]
                        #check cycle2 drug
                        for k in range(len(newSeq)-22):
                            tempSeq=newSeq[k:k+11]
                            drugId2=dicDrugBC_hash_c2.get(tempSeq, 'NA')
                            if drugId2!='NA':
                                count2+=1
                                newSeq2=newSeq[k+11:]
                                #check cycle3 drug
                                for l in range(len(newSeq2)-11):
                                    tempSeq=newSeq2[l:l+11]
                                    drugId3=dicDrugBC_hash_c3.get(tempSeq, 'NA')
                                    if drugId3!='NA':
                                        count3+=1                                    
                                        infoFile.write('%s,%d,%s,%s-%s-%s\n'%(readId,2,libId,drugId1,drugId2,drugId3))
                                        break                    
            if i==4:
                i=0

    statsFile.write('%s,%d,%d,%d\n'%(libId,count1,count2,count3))
infoFile.close()
statsFile.close()
print ('haha')


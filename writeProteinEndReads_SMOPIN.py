import sys

outDirPath=sys.argv[1]
jobId=sys.argv[2]
r1Path=sys.argv[3]
r2Path=sys.argv[4]

#assign reads to protein end
dicId_r1r2={}
readIdSet=set()
with open('%s/%sintermediateFiles/endInfo.csv'%(outDirPath,jobId),'r') as f:
    next(f)
    for line in f:
        splitLine=line.strip().split(',')
        readId,drugEnd=splitLine[0],splitLine[1]
        readIdSet.add(readId)
        dicId_r1r2[readId]=drugEnd
        
#write fastq files
targetFile1=open('%s/%sprocessedFastq/proteinEnd.fastq'%(outDirPath,jobId),'w')
#R1
temp=[]
flag1=False 
i=0
with open('%s'%(r1Path),'r') as f:
        for line in f:
            temp.append(line)
            i+=1
            if i==1:
                splitLine=line.strip().split()
                readId=splitLine[0][1:]
                if readId in readIdSet:
                    #drug end is read2 = protein end is read1
                    if dicId_r1r2[readId]=='2':
                        flag1=True
            if i==4:
                if flag1:
                    for ha in temp:
                        targetFile1.write(ha)
                i=0
                temp=[]
                flag1=False
#R2                
temp=[]
flag1=False
i=0
with open('%s'%(r2Path),'r') as f:
        for line in f:
            temp.append(line)
            i+=1
            if i==1:
                splitLine=line.strip().split()
                readId=splitLine[0][1:]
                if readId in readIdSet:
                    #drug end is read1 = protein end is read2
                    if dicId_r1r2[readId]=='1':
                        flag1=True
            if i==4:
                if flag1:
                    for ha in temp:
                        targetFile1.write(ha)
                i=0
                temp=[]
                flag1=False 
            
targetFile1.close()
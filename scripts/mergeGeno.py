import argparse
import os
import re
import subprocess
from collections import defaultdict
from subprocess import PIPE


def processSnplist(fileHandles):
    nFiles = len(fileHandles)
    snpDict = defaultdict(int)
    for fileC, snpFile in enumerate(fileHandles):
        for snp in snpFile:
            snpParse = snp.strip()
            snpDict[snpParse] += 1
    return(snpDict)

def processOut(out, snpDict, nFiles):
    snpCounter = 0
    with open(out, 'w') as outFile:
        for snp, val in snpDict.items():
            if val == nFiles:
                snpCounter += 1
                outFile.write(snp+'\n')
    if snpCounter > 0:
        print('MATCHED {} SNPS ACROSS {} BEDFILES'.format(snpCounter, nFiles))
    else:
        print('UNABLE TO FIND CONSENSUS AMONG THE SNPLISTS')

def checkPlink():
    plinkC = subprocess.Popen('which plink', stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = plinkC.communicate()
    if not stderr:
        print('FOUND PLINK {}'.format(stdout.decode()))
        return True
    else:
        raise FileNotFoundError('PLINK IS NOT IN PATH ??')

def makeConsensusBed(bedName, commonSnpList, exclude, ntry):
    if ntry > 0: # first time
        outPath = bedName+'_Try{}'.format(ntry)
    else: #other times
        outPath=re.sub('Try[\d+]', 'Try{}'.format(ntry+1),bedName)
    if os.path.exists(bedName+'.bed'):
        if not exclude:
            command = 'plink --bfile {} --extract {} --make-bed --out {}'.format(bedName, commonSnpList, outPath)
        else:
            command = 'plink --bfile {} --exclude {} --make-bed --out {}'.format(bedName, commonSnpList, outPath)
        extractSnp = subprocess.Popen(command,  stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = extractSnp.communicate()
        if 'done' in stdout.decode():
            print('EXTRACT SUCCESSFUL FROM BASE {} TO NEW BED {}'.format(bedName, outPath))
            return outPath
        else:
            print(stderr.decode())
            raise FileNotFoundError('EXTRACT FAILED CHECK PATH OF {}'.format(bedName))
    else:
        raise FileNotFoundError('CHECK {} '.format(bedName))
    
def processBedextract(bedList, commonSnpList, exclude, ntry):
    newBeds = []
    for bed in bedList:
        print(bed)
        out=makeConsensusBed(bedName=bed, commonSnpList=commonSnpList, exclude=exclude, ntry=ntry)
        if out:
            newBeds.append(out)
        else:
            raise ValueError('CHECK BED SUBSET {}'.format(bed))
    if newBeds:
        return newBeds
    else:
        return False

def bedMerge(newBeds, ntry, bedOut):
    outList = 'mergeList.txt'
    mergedBed = '{}_Try{}'.format(bedOut, ntry)
    outF=open(outList, 'w')
    for n, bed in enumerate(newBeds):
        if n > 0:
            outF.write(bed+'_Try{}\n'.format(ntry))
        else:
            baseBed = bed+'_Try{}'.format(ntry)
    outF.close()
    command = "plink --bfile {} --merge-list {} --make-bed --out {}".format(baseBed, outList, mergedBed)
    print(command)
    mergeCommand = subprocess.Popen(command,  stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = mergeCommand.communicate()
    print(stderr)
    if 'missnp' in stderr.decode(): ## trigger consensus bed
        print('CHECK MISNP OUTPUT')
        return False
    else:
        return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=argparse.FileType('r'), nargs='+', help='input a list of snplists from different GWAS batches')
    parser.add_argument('-out',required=True, help='outfile a stringname with a suffix .txt')
    parser.add_argument('-bedList', required=True, help='merged output of the binary files')
    parser.add_argument('-bedOut', required=True, help='merged output of the binary files')
    args = parser.parse_args()
    fileHandles = args.file
    out = args.out
    bedOut = args.bedOut
    bedList = [i.strip() for  i in open(args.bedList)]
    nFiles = len(fileHandles)
    print('FOUND {} SNPLISTS'.format(nFiles))
    snpDict=processSnplist(fileHandles)
    processOut(out, snpDict, nFiles) # will write out consensus list
    if checkPlink(): ## check whether plink is present
        exclude=False
        ntry = 1
        consensusBedlist=processBedextract(bedList, out, exclude, ntry)
    if consensusBedlist:
        print('PROCEEDING TO BEDMERGE')
        mergeCall=bedMerge(newBeds=bedList, ntry=ntry, bedOut=bedOut)
        if not mergeCall: ## retry but passing misnp error
            misnpFile = '{}_Try{}-merge.missnp'.format(bedOut, ntry)
            ntry += 1
            print('3+ ALLELEIC VARIANTS DOING MISNP')
            exclude = True
            consensusBedlist = processBedextract(bedList, misnpFile, exclude, ntry)
            if consensusBedlist:
                print('PROCEEDING TO BEDMERGE')
                mergeCall=bedMerge(newBeds=bedList, ntry=ntry, bedOut=bedOut)
                if mergeCall:
                    print('MERGE SUCESSFUL')

if __name__ == '__main__':main()

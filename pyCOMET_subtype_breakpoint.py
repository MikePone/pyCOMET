#!/usr/bin/env python3

'''
    pyCOMET_subtype.py

    a python3 based program to use PPMD models for subtype classification

    The program takes as input the following:
        - a model description file in msgpack format
        - a fixed context size

    This model gives log likelihoods of the probablities of seeing a residue at
    a particular position given the context

    These likelihoods then can be used to predict the subtype
    
    if a query is predicted as potential recombinant, the following is done:
        - average likelihood value is calculated for each nucleotide position within a flanking window of size 15
        - This is done for all subtypes returned by COMET's decision tree
        - a graph is created with these likelihood values across the sequence length for all reported subtypes
        - This graph can aid the identification of potential recombination breakpoints  
'''

import math, os, sys, argparse, time
import umsgpack
import json
from Bio import SeqIO
import numpy as np

import subprocess

#import matplotlib.pyplot as plt

##********************************************************##
def contextThreshold(x):
    try:
        x = int(x)
    except ValueError as e:
        sys.exit(e)

    if x <= 0:
        raise argparse.ArgumentTypeError('%r must be a positive integer' % (x,))

    return x
#*******************************************************

#*******************************************************
def getNonRedundantListOrder(lst):
    seen = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]
#*******************************************************    
    
#*******************************************************    
def pLike(c,ctx,cLen,sType,pScore=0.0):
    
    conList = [] # holds exclusion kmers
    kSize = cLen # current value for k
    for z in range(cLen+1):
        cCtx = ctx[z:] 
        
        if cCtx not in ppmd_models[sType][kSize]:
            pass
        
        else: # k-mer present in the model
            exclusion = 0
            for kmer in conList:
                if kmer in ppmd_models[sType][kSize][cCtx]:
                    exclusion += ppmd_models[sType][kSize][cCtx][kmer]

            
            sumCount = sum(ppmd_models[sType][kSize][cCtx].values())
            
            if c in ppmd_models[sType][kSize][cCtx]:
                try:
                    prob = float(ppmd_models[sType][kSize][cCtx][c]) / (sumCount - exclusion)
                    pScore += math.log10(prob)
                    #print(prob,pScore)
                    return pScore
                except (ZeroDivisionError, ValueError) as e:
                    sys.exit(e)
            else:
                conList = list(ppmd_models[sType][kSize][cCtx].keys())
                conList.remove('esc')
                
                escape = ppmd_models[sType][kSize][cCtx]['esc']
                
                try:
                    prob = float(escape) / (sumCount - exclusion)
                    pScore += math.log10(prob)
                except (ZeroDivisionError, ValueError) as e:
                    sys.exit(e)
                
        kSize -= 1
    
    return pScore

#*******************************************************

#*******************************************************    
def pProb(c,ctx,cLen,sType,pScore=1.0):
    
    conList = [] # holds exclusion kmers
    kSize = cLen # current value for k
    for z in range(cLen+1):
        cCtx = ctx[z:] 
        
        if cCtx not in ppmd_models[sType][kSize]:
            pass
        
        else: # k-mer present in the model
            exclusion = 0
            for kmer in conList:
                if kmer in ppmd_models[sType][kSize][cCtx]:
                    exclusion += ppmd_models[sType][kSize][cCtx][kmer]

            
            sumCount = sum(ppmd_models[sType][kSize][cCtx].values())
            
            if c in ppmd_models[sType][kSize][cCtx]:
                try:
                    prob = float(ppmd_models[sType][kSize][cCtx][c]) / (sumCount - exclusion)
                    #pScore += math.log10(prob)
                    pScore *= prob
                    #print(prob,pScore)
                    return pScore
                except (ZeroDivisionError, ValueError) as e:
                    sys.exit(e)
            else:
                conList = list(ppmd_models[sType][kSize][cCtx].keys())
                conList.remove('esc')
                
                escape = ppmd_models[sType][kSize][cCtx]['esc']
                
                try:
                    prob = float(escape) / (sumCount - exclusion)
                    #pScore += math.log10(prob)
                    pScore *= prob
                except (ZeroDivisionError, ValueError) as e:
                    sys.exit(e)
                
        kSize -= 1
    
    return pScore

#*******************************************************

#*******************************************************
def challenge(sLike,target,nSubtypes,start,end,thr):
    '''This function computes the sum of likelihoods for each subtypes 
    in the given window and finds the most likely subtype 
    '''
    
    # get the sum of likelihoods for each subtypes
    #sumLL = [sum(sLike[i][start:end])for i in range(nSubtypes)]
    sumLL = np.sum(sLike[:,start:end],axis=1)
    
    # find the maximum likelihood and index
    #maxLL = max(sumLL)
    maxLL = np.amax(sumLL)
    
    #maxIndex = sumLL.index(maxLL)
    maxIndex = np.argmax(sumLL)

    # return best matching subtype according to COMET's decision tree
    if (maxLL - sumLL[target]) <= thr:
        return target
    else:
        return maxIndex
#*******************************************************

#*******************************************************
def detectBreakpoints(assignedSubtypes,subtypes,query,wSize,args,zName):

    
    # convert sequences into upper case and remove gaps
    qSeq = str(query.seq).upper().replace('-','')   
    
    qLen = len(qSeq)
    
    # get number of assigned subtypes
    nParents = len(assignedSubtypes)

    # create a matrix to hold the probability values   
    pMatrix = np.zeros((nParents,qLen))
    
    # create a numpy matrix to hold the average probability  values
    aMatrix = np.zeros((nParents,qLen))
    
    cSize = args.context
    
    #print(assignedSubtypes)
    
    #print(cSize,nParents,qLen)
    
    # generate probability for each site for each subtype
    for r in range(nParents):
       probs = np.zeros(qLen)
       
       # get the probabilities
       for j in range(qLen):
           c = qSeq[j]
           start = j - cSize if j > cSize else 0
           ctx = qSeq[start:j]
           
           lctx = j - start
           
           if set(ctx).issubset('ACGT'):
               pPos = pProb(c,ctx,lctx,subtypes[assignedSubtypes[r]])
               probs[j] = pPos
           
       pMatrix[r] = probs
           
    #print(query.id,'\n',pMatrix[:,10:15])
         
    
    # calculate average scores at each position for each assigned subtypes
    # for insufficient values on both side, average will be calculated at 
    # positions 7..(qLen-7)
    hwSize = wSize // 2
    lastPos = qLen - hwSize
    
    sName = zName + '/' + query.id + '.prob.txt'
    pName = zName + '/' + query.id + '.prob.png'
    
    lName = zName + '/' + 'breakpoints.log'
    lh = open(lName,'a')
    
    # create the probability distribution file
    fh = open(sName,'w')
    fh.write('Positions')
    
    for sub in assignedSubtypes:
        sType = subtypes[sub]
        fh.write('\t{}'.format(sType))
    #fh.write('\n')
    
    for i in range(qLen):
        fh.write('\n{}'.format(str(i)))
        
        # for all subtypes in assignedSubtypes
        for k in range(nParents):
            #print(k)
            #if i < hwSize:
                #start = 0
            #else:
                #start = i - hwSize
 
            start = 0 if i < hwSize else (i-hwSize)

            end = i + hwSize + 1
            
            tLL = pMatrix[k][start:end]
            aMatrix[k][i] = np.average(tLL)
                        
            #if i < hwSize or i >= lastPos:
                #aMatrix[k][i] = pMatrix[k][i]
            #else:
                #start = i - hwSize
                #end = i + hwSize + 1
                #tLL =  pMatrix[k][start:end]
                #aMatrix[k][i] = np.average(tLL)
            
            fh.write('\t{}'.format(aMatrix[k][i]))
    
        #fh.write('\n')
    
    # create the probability distribution plot using R script
    cl = ['Rscript','analyseBreakpoints.R',sName,pName, query.id]
    try:
        subprocess.check_call(cl,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
        pass
        
    fh.close()
    
    
    
    #plt.plot(aMatrix[0])
    #plt.xlabel('Nucleotide positions')
    #plt.ylabel('Probability')
    #plt.show()
    return
    
#*******************************************************

#****************************************************
def check_subtype(sLike,query,subtypes,nSubtypes,qLen,pIndex,args,zName):
    '''
        Uses COMET's desition tree to call subtypes for a query
            sLike: likelihood matrix
            seqID: sequence identifier
            subtypes: list of subtype names
            nSubtypes: number of reference subtypes
            qLen: length of the query sequence
            args: command line arguments
        
    '''
    # get the sequence ID
    seqId = query.id
    
    thr = args.thr
    # get the sum of likelihoods for all subtypes
    #sumLL = [sum(sLike[i]) for i in range(nSubtypes)]
    sumLL = np.sum(sLike,axis=1)
    #print(sumLL)

    # find the index of most likely subtype PURE/CRF
    #maxS = max(sumLL)
    maxS = np.amax(sumLL)
    
    #S = sumLL.index(maxS)
    S = np.argmax(sumLL)
    
    # find the most likely PURE subtype
    #maxPS = max(sumLL[pIndex:])
    maxPS = np.amax(sumLL[pIndex:])
    #PS = (sumLL[pIndex:].index(maxPS)) + pIndex
    PS = np.argmax(sumLL[pIndex:]) + pIndex
    
    #print(seqId, subtypes[S],maxS,subtypes[PS],maxPS)
    #with open(args.outFile,'a') as fh:
        #fh.write('{}\t{}\t{}\t{}\t{}\n'.format(seqId, subtypes[S],maxS,subtypes[PS],maxPS)) 
    #return
    
    wSize = args.wSize # set window size
    bSize = args.bSize # set the step size

    # Pre-compute number of windows
    numOfWindows = int((qLen-wSize)/bSize)+1
    #print('numOfWindows',numOfWindows)
    
    # check the PURE subtype first
    # create a list of subtype assignment for each window
    subAssignment = [PS]*numOfWindows
    iWindow = 0    

    # get the most likely subtype for each window
    for i in range(0,numOfWindows*bSize,bSize):
        start = i
        end = i + wSize
        #print(iWindow)
        subAssignment[iWindow] = challenge(sLike,PS,nSubtypes,start,end,thr)
        iWindow += 1
        
    
    if S == PS:
        assignedSubtypes = getNonRedundantListOrder(subAssignment)
        if len(assignedSubtypes) == 1:
            msg = seqId + '\t' + subtypes[PS] + '\t(PURE)'
            #print(msg)
            #detectBreakpoints(assignedSubtypes,subtypes,query,101,args)
            return msg
        else:
            msg = seqId + '\t' + 'unassigned_1\t' 
            for asub in assignedSubtypes:
                msg += subtypes[asub] + ' '
            #print(msg)
            detectBreakpoints(assignedSubtypes,subtypes,query,101,args,zName)
            return msg
    else: # S != PS
        assignedSubtypes = getNonRedundantListOrder(subAssignment)
        if len(assignedSubtypes) == 1:
            msg = seqId + '\t' + subtypes[PS] + '\t(Check ' + subtypes[S] + ')'
            #print(msg)
            #detectBreakpoints([S,PS],subtypes,query,101,args)
            return msg
        else: # needs checking CRF
            subAssignment = [S]*numOfWindows
            iWindow = 0    
            # get the most likely subtype for each window
            for i in range(0,numOfWindows*bSize,bSize):
                start = i
                end = i + wSize
                #print(iWindow)
                subAssignment[iWindow] = challenge(sLike,S,nSubtypes,start,end,thr)
                iWindow += 1
            
            assignedSubtypes = getNonRedundantListOrder(subAssignment)
            if len(assignedSubtypes) == 1:
                msg = seqId + '\t' + subtypes[S] + '\t(CRF)'
                #print(msg)
                #detectBreakpoints(assignedSubtypes,subtypes,query,101,args)
                return msg
            else:
                msg = seqId + '\t' + 'unassigned_2\t' 
                for asub in assignedSubtypes:
                    msg += subtypes[asub] + ' '
                #print(msg)
                detectBreakpoints(assignedSubtypes,subtypes,query,101,args,zName)
                return msg
#****************************************************            

#def calculateProbability(query,)

#****************************************************
def calculateLogLikelihood(query,subtypes,nSubtypes,pIndex,args,zName):
    '''
        Calculates likelihood values for a sequence
    '''
    # convert sequences into upper case and remove gaps
    qSeq = str(query.seq).upper().replace('-','')
    
    
    qLen = len(qSeq)

    # Create a list to hold likelihoods for all the nucleotide positions
    # each row represents likelihoods generated by each of the subtype PPMD models

    #lMatrix = [[]]*nSubtypes
    lMatrix = np.zeros((nSubtypes,qLen))
    
    # get the context size
    
    cSize = int(args.context)

    # generate likelihoods based on each subtype models
    for r in range(nSubtypes):
        # 'bits' holds likelihood values of all the nucleotide positions
        #bits = [0.0]*qLen
        bits = np.zeros(qLen)
        # get likelihood for all the residue positions
        for j in range(qLen):  # qLen
            # the nucleotide at position 'j'
            c = qSeq[j]
            
            # do not calculate for ambiguous characters
            #if c not in dna:
            #    continue

            # extract the current context
            start = j - cSize if j > cSize else 0
            ctx = qSeq[start:j]
 
            lctx = j-start

            # calculate likelihood for only 'ACGT'
            if set(ctx).issubset('ACGT'):
                lPos = pLike(c,ctx,lctx,subtypes[r])
                bits[j] = lPos
            #lPos = pLike(c,ctx,lctx,subtypes[r])
            #if lPos == None:
                #print(j,subtypes[r],c,ctx)
            
            #bits[j] = lPos
            
            #lPos = predictLikelihood(c,ctx,lctx,[],subtypes[r],0.0)
            #bits[j] = sum(bitList)

        #break
        #print(bits)
        #print(loaded_subtypes[r],sum(bits))

        ## create a list of likelihood values
        ## for each query based on each reference
        ## for each residue position
        lMatrix[r] = bits

    #print(sLike[0][0:10])
    #print(len(sLike))

    return check_subtype(lMatrix,query,subtypes,nSubtypes,qLen,pIndex,args,zName)

#****************************************************

##**************************************************

def getArguments():
    '''
        Parse all the command line arguments from the user
    '''

    parser = argparse.ArgumentParser(description='Predicts subtype of a sequence based on PPMD models trained using reference sequences', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-q','--query',required=True,help='Query sequence file in FASTA format for subtype classification')
    parser.add_argument('-c','--context',type=contextThreshold,default=8,help='Context size for PPMD model (default: 8)')
    parser.add_argument('-w','--wSize',type=contextThreshold,default=100,help='Window size for COMET decision tree (default: 100)')
    parser.add_argument('-b','--bSize',type=contextThreshold,default=3,help='Step size of the windows for COMET decision tree (default: 3)')
    parser.add_argument('-m','--modelFile',required=True,help='PPMD model file for reference sequences')
    parser.add_argument('-o','--outFile',required=True,help='Output file for subtype prediction results')
    parser.add_argument('-t','--thr',type=contextThreshold,default=28,help='Threshold difference used in decision tree (default: 28)')


    args = parser.parse_args()

    return args

##********************************************************##

#***********************************************************

if __name__=="__main__":
    args = getArguments()

    ## get context size
    cSize = int(args.context)

    ## Create dna alphabet
    #dna = ['A','C','G','T']

    ## define the digits list
    #digits = list('0123456789')

    ## read in the model file and populate the ppmd list
    #fh = open(args.modelFile,'rb')
    #ppmd_models = umsgpack.load(fh)
    try:
        with open(args.modelFile,'r') as fh:
            ppmd_models = json.load(fh)
    except FileNotFoundError as e:
        eMsg = '\nThe pyCOMET model file <{}> could not be found.'.format(args.modelFile)
        eMsg += ' Please try again with correct model file name.\n'
        print(eMsg)
        sys.exit()

    ## get names of subtype present in the loaded PPMD models
    loaded_subtypes = sorted(list(ppmd_models.keys()))

    ## Remove 'CPZ' from subtype list
    try:
        loaded_subtypes.remove('CPZ')
    except:
        pass

    #print(ppmd_models['A1'][8])
    #sys.exit(0)

    nSubtypes = len(loaded_subtypes)
    
    # get the index of 'A1'; this marks the start index of PURE subtypes
    pIndex = loaded_subtypes.index('A1')
       
    ## read in the query sequences
    try:
        qSeqs = list(SeqIO.parse(args.query,'fasta'))
    except FileNotFoundError as e:
        eMsg = '\nThe query sequence file <{}> could not be found.'.format(args.query)
        eMsg += ' Please try again with correct sequence file name.\n'
        print(eMsg)
        sys.exit()

    ## check if the sequences were read properly
    if len(qSeqs) == 0:
        msg = 'Query sequences were not read properly'
        msg += '\nPlease run again with valid FASTA sequence file with at least one sequence\n'
        sys.exit(msg)
        
        
    # Create a directory to stote the output of the pyCOMET run
    timeNow = time.strftime('%Y-%m-%d-%H%M%S')
    zName = 'pyCOMET.' + timeNow 
    print(zName)
    try:
        os.mkdir(zName)
        
    except OSError as e:
        sys.exit(e)
                  
    oName = zName + '/' + args.outFile
    
    #sys.exit()
    
    ## open output file for storing predicted subtypes
    #fh = open(args.outFile,'w')

    # calls calculateLogLikelihood() to generate subtypes
    for query in qSeqs:  #len(seqs)
        #print(query.id,len(query.seq))
        # calculate likelihood matrix
        tMsg = calculateLogLikelihood(query,loaded_subtypes,nSubtypes,pIndex,args,zName)
        
        print('{}'.format(tMsg))
        
        with open(oName,'a') as fh:
            fh.write('{}\n'.format(tMsg))
        
        
        
    #fh.close()
#***********************************************************

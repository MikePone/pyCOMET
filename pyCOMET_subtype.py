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
'''

import math, os, sys, argparse, time
import umsgpack
import json
from Bio import SeqIO
import numpy as np

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
    if (maxLL - sumLL[target]) <= 28:
        return target
    else:
        return maxIndex
#*******************************************************

#****************************************************
def check_subtype(sLike,seqId,subtypes,nSubtypes,qLen,pIndex,args):
    '''
        Uses COMET's desition tree to call subtypes for a query
            sLike: likelihood matrix
            seqID: sequence identifier
            subtypes: list of subtype names
            nSubtypes: number of reference subtypes
            qLen: length of the query sequence
            args: command line arguments
        
    '''
    thr = 28
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
    
    #print(subtypes[S],maxS,subtypes[PS],maxPS)
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
            return msg
        else:
            msg = seqId + '\t' + 'unassigned_1\t' 
            for asub in assignedSubtypes:
                msg += subtypes[asub] + ' '
            #print(msg)
            return msg
    else: # S != PS
        assignedSubtypes = getNonRedundantListOrder(subAssignment)
        if len(assignedSubtypes) == 1:
            msg = seqId + '\t' + subtypes[PS] + '\t(Check ' + subtypes[S] + ')'
            #print(msg)
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
            
            if len(assignedSubtypes) == 1:
                msg = seqId + '\t' + subtypes[S] + '\t(CRF)'
                #print(msg)
                return msg
            else:
                msg = seqId + '\t' + 'unassigned_2\t' 
                for asub in assignedSubtypes:
                    msg += subtypes[asub] + ' '
                #print(msg)
                return msg
#****************************************************            

#****************************************************
def calculateLogLikelihood(query,subtypes,nSubtypes,pIndex,args):
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

            #if len(set(dna) | set(ctx)) == 4: # no ambiguous characters
                # holds the likelihood values for a nucleotide site
                # unsuccessful match with a longer context results in searching
                # with a smaller context hence multiple likelihood values
                # employs the idea: log(a*b*c) = log(a)+log(b)+log(c)
                #bitList = list()

                # calls the predict function
                # returns the bitList with the likelihood values
                #predictLikelihood(c,ctx,lctx,[],loaded_subtypes[r],bitList)

                # combines the probability from each context into a single likelihood
                #bits[j] = sum(bitList)
            #print(c,ctx,lctx,subtypes[r])
            
            
            lPos = pLike(c,ctx,lctx,subtypes[r])
            #if lPos == None:
                #print(j,subtypes[r],c,ctx)
            bits[j] = lPos
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

    return check_subtype(lMatrix,query.id,subtypes,nSubtypes,qLen,pIndex,args)

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

    ## open output file for storing predicted subtypes
    #fh = open(args.outFile,'w')

    # calls calculateLogLikelihood() to generate subtypes
    for query in qSeqs:  #len(seqs)
        #print(query.id,len(query.seq))
        # calculate likelihood matrix
        tMsg = calculateLogLikelihood(query,loaded_subtypes,nSubtypes,pIndex,args)
        
        print('{}'.format(tMsg))
        
        with open(args.outFile,'a') as fh:
            fh.write('{}\n'.format(tMsg))
        
        
        
    #fh.close()
#***********************************************************

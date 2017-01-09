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
from Bio import SeqIO

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
# This is a function implementing recursion
def predictLikelihood(c,ctx,cLen,conList,sType,bitList):
    '''
    predictLikelihood() returns the log likelihood values of predicting a nucleotide given a context
    argument list:
        c - nucleotide
        ctx - current context
        cLen - length of the context
        conList - list of nucleotides matched for the previous longer context
        sType - subtype in 'loaded_subtype' for which likelihood is calculated
        bitList - list of log likelihood values for all the context sizes
    '''
    # cLen == -1 means no match was found in any context
    # returns with likelihood values
    if cLen == -1:
        return
    
    # if the context is not present in the PPMD model for the subtype
    # recursively call predict with the smaller context
    # does not assign any likelihood
    if ctx not in ppmd_models[sType][cLen]:
        predictLikelihood(c,ctx[1:],cLen-1,[],sType,bitList)
    
    # if context present in the PPMD model for sType
    else:
        # copy the dictionary of the current context for the subtype
        cContext = ppmd_models[sType][cLen][ctx].copy()
        
        # create variable to hold the count for nucleotides matched 
        # in the previous context for exclusion
        exclusion = 0
        
        
        # find total count of the nucleotides followed by the previous context
        # conList contains the nucleotides present in the previous context
        for key in conList:
            if key in cContext:
                exclusion += cContext[key]
        
        # we need to calculate the escape count
        # in PPMD escape count is = number of nucleotides / 2
        escape = len(cContext.keys()) / 2
        
        # We need the total count of nucleotides for probability calculation
        cSum = sum(cContext.values())
        
        # calculate probability if nucleutide follows the current context
        if c in cContext:
            # prob = count_of_c / (sum_of_count + escape - exclusion)
            try:
                prob = float(cContext[c]) / (cSum + escape - exclusion)
            except ZeroDivisionError as e:
                sys.exit(e)
            
            # calculate the log likelihood and append to bitList
            bitList.append(math.log(prob))
            # return to the calling function
            return
        
        else:
            # if 'c' does not follow the current context
            # call predictLikelihood() with the smaller context
            # we need to find the exclusion characters to be passed on to the smaller context
            exc = list(cContext.keys())
        
            # log likelihood of 'escape' will be added to bitList
            try:
                prob = float(escape) / float(cSum + escape - exclusion)
            except ValueError as e:
                sys.exit(e)
            bitList.append(math.log(prob))
        
            # call predictLikelihood()
            predictLikelihood(c,ctx[1:],cLen-1,exc,sType,bitList)

#***************************************************

#****************************************************
def challangeType(sLike,target,start,end):
    '''
        calculates sum of log likelihood for each reference in the window
        if other - PS > 28 for all others returns 'True', else returns 'False'
    '''
    
    psLike = sLike[target][start:end]
    
    for i in range(len(sLike)):
        #if ppmd[i]['type'] == ppmd[target]['type']:
            #continue
        oLike = sLike[i][start:end]
                
        diff = sum(oLike) - sum(psLike)
        
        if diff > 28:
            #print('challenge: false', i, target, diff)
            return False
        
    return True
#****************************************************

#****************************************************
def check_subtype(sLike,seqId,args):
    '''
        Uses COMET's desition tree to call subtypes for a query
    '''
    # get the number of reference subtypes
    # and the size of the likelihood array
    numRef = len(sLike)
    numSites = len(sLike[0])
    
    S = None # holds the index for most likely subtype
    PS = None # holds the index for most likely PURE subtype
    
    sll = None # holds the minimum likelihood for S
    psll = None # holds the minimum likelihood for PS 
    
    sFlag = False # S and PS are not the same
    
    ## decide S and PS
    for i in range(numRef):
        lsum = sum(sLike[i])
        #print(lsum)
        
        if not sll: # sll yet to define
            sll = lsum
            S = i
        elif sll < lsum: # better match found
            sll = lsum
            S = i     
        
        if loaded_subtypes[i][0] not in digits: # reference is a PURE subtype
            if not psll: # psll yet to define
                psll = lsum
                PS = i
            elif psll < lsum:
                psll = lsum
                PS = i
    
    #print(S,sll,PS,psll) 
    
    wSize = args.wSize # set window size
    bSize = args.bSize # set the step size
    
    # Pre-compute number of windows
    numOfWindows = int((numSites-wSize)/bSize)+1
    #print('numOfWindows',numOfWindows)

    if S == PS:
        ## check whether PURE subtype is best match 
        ## in each of the windows
        #print('Checking S', S)
        bMatch = True
        
        for i in range(0,numOfWindows*bSize,bSize):
            start = i
            end = i + wSize
            bMatch = challangeType(sLike,S,start,end)
        
            if not bMatch: # challengeType reurned False
                return 'UNASSIGNED, true'
        
        if bMatch:
            rTxt = loaded_subtypes[S] + ' (PURE)'
            return rTxt    
    
    else: # S != PS
        bMatch = True
        for i in range(0,numOfWindows*bSize,bSize):
            start = i
            end = i + wSize
            bMatch = challangeType(sLike,PS,start,end)
            
            if not bMatch: # challengeType returned False
                break
        
        if bMatch: # 'True' found but needs to check CRF - PS(check S)
            rTxt = loaded_subtypes[PS] + ' (check ' + loaded_subtypes[S] + ')'
            return rTxt
        
        else: # needs to check CRF
            bMatch = True
            for i in range(0,numOfWindows*bSize,bSize):
                start = i
                end = i + wSize
                bMatch = challangeType(sLike,S,start,end)
                
                if not bMatch: # returns False
                    return 'UNASSIGNED, crf'
            
            if bMatch:
                rTxt = loaded_subtypes[S] + ' (CRF)'
                return rTxt
#****************************************************

#****************************************************
def calculateLogLikelihood(qSeq,seqId,args):
    '''
        Calculates likelihood values for a sequence
    '''
    # Create a list to hold likelihoods for all the nucleotide positions
    # each row represents likelihoods generated by each of the subtype PPMD models
    sLike = list()
    
    dna = ['A','T','C','G']
    
    # generate likelihoods based on each subtype models
    for r in range(len(loaded_subtypes)): 
        # 'bits' holds likelihood values of all the nucleotide positions
        bits = list() 
        # get likelihood for all the residue positions
        for j in range(len(qSeq)):
            # the nucleotide at position 'j'
            c = qSeq[j]
            # do not calculate for ambiguous characters
            if c not in 'ATCG':
                continue
            
            # extract the current context
            start = j - cSize if j > cSize else 0
            ctx = qSeq[start:j]
            
            if len(set(dna) | set(ctx)) == 4: # no ambiguous characters
                # holds the likelihood values for a nucleotide site
                # unsuccessful match with a longer context results in searching 
                # with a smaller context hence multiple likelihood values
                # employs the idea: log(a*b*c) = log(a)+log(b)+log(c)
                bitList = list()
                
                # calls the predict function 
                # returns the bitList with the likelihood values
                predictLikelihood(c,ctx,len(ctx),[],loaded_subtypes[r],bitList)
                
                # combines the probability from each context into a single likelihood
                bits.append(sum(bitList)) 
            
        #print(bits)
        #print(loaded_subtypes[r],sum(bits))
            
        ## create a list of likelihhod values
        ## for each query based on each reference
        ## for each residue position
        sLike.append(bits)
            
    #print(sLike[0][0:10])    
    #print(len(sLike))    
    
    return check_subtype(sLike,seqId,args)

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


    args = parser.parse_args()
    
    return args
    
##********************************************************##

#***********************************************************

if __name__=="__main__":
    args = getArguments()

    ## get context size
    cSize = int(args.context)
        
    ## Create dna alphabet
    dna = ['A','T','C','G']
    
    ## define the digits list
    digits = list('0123456789') 
    
    ## read in the model file and populate the ppmd list
    fh = open(args.modelFile,'rb')
    ppmd_models = umsgpack.load(fh)
    fh.close()

    ## get names of subtype present in the loaded PPMD models
    loaded_subtypes = list(ppmd_models.keys())
    
    ## read in the query sequences
    seqs = list(SeqIO.parse(args.query,'fasta'))
 
    ## check if query file is empty
    if len(seqs) == 0:
        msg = 'Query file does not have valid FASTA sequences.'
        msg += '\nPlease run again with valid FASTA sequence file\n' 
        sys.exit(msg)
    
    ## open output file for storing predicted subtypes
    fh = open(args.outFile,'w')
     
    # calls calculateLogLikelihood() to generate subtypes
    for i in range(len(seqs)):  #len(seqs)
        # convert sequences into upper case and remove gaps
        qSeq = str(seqs[i].seq).upper().replace('-','')
        
        # call likelihood function that returns a string
        tMsg = calculateLogLikelihood(qSeq,seqs[i].id,args)
        fh.write('{}\t{}\n'.format(seqs[i].id,tMsg))
        #print(seqs[i].id, tMsg)

    fh.close()
#***********************************************************

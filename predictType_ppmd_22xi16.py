#!/usr/bin/env python3

'''
    predictType_ppmd.py
    
    a python3 based program to use PPMD models for subtype classification
     
    The program takes as input the following:
        - a model description file in plain text format
        - a fixed context size
    
    This model gives log likelihoods of the probablities of seeing a residue at 
    a particular position given the context
    
    These likelihoods then can be used to predict the subtype
'''

import math, os, sys, argparse, time

from Bio import SeqIO

#bitList = list() # likelihood values for each context for each residue


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

##*****************************************************************************
## predict gives a log likelihood of finding a character 
## at a particular site given the context
##*****************************************************************************
def predict(c, ctx, s, conList,pTable,bitList):
    #print('***')
    #print(c,ctx,s,conList)
    #print(pTable[0])
    

    if s == -1:
        #bitList.append(-math.log(1/len(dna),2))
        return
        
    if ctx in pTable[s]:
        cp = pTable[s][ctx].copy()
        #print(ctx,c,cp)
        exclusion = 0
        for key in conList:
            if key in cp and key in 'ATCG':
                #cp.pop(key)
                #ch = cp.pop(key)
                #print('ch',key)
                exclusion += cp[key] # only contains count of ATCG
                #print('exclusion',exclusion)
        # escape will only count presence of ATCG and divide by 2
        escape = 0.5 * len(set(cp.keys()) & set('ATCG'))
        #escape = len(cp.keys())/float(2)
        #print('escape',escape)
       
        # csum will have the count of 'ATCG' residues
        csum = 0
        if 'A' in cp.keys():
            csum += cp['A']  
        if 'T' in cp.keys():
            csum += cp['T']  
        if 'C' in cp.keys():
            csum += cp['C']  
        if 'G' in cp.keys():
            csum += cp['G']  

        #csum = sum(cp.values())
        #print('escape',escape,'csum',csum, 'exclusion',exclusion)

        if c in cp:
            try:
                prob = float(cp[c]) / float(escape + csum - exclusion)
                #print('prob',prob,c,cp)
            except ZeroDivisionError as e:
                #print(e)
                sys.exit()
            bitList.append(math.log(prob,2))
            #probList.append(prob)
            #print output % ("' '" if c == ' ' else c, prob)
            #o1 = "' '" if c == ' ' else c
            #print(output % ("' '" if c == ' ' else c,prob))
        else:
            if csum > 0:
                try:
                    prob = float(escape) / float(escape + csum - exclusion)
                except ValueError as e:
                    sys.exit('ValueError')
                bitList.append(math.log(prob,2))
                #probList.append(prob)
                #print(output % (escape, prob))
                #print(c,ctx[1:],s-1,impossible + list(cp.keys()))
                predict(c, ctx[1:], s-1, list(cp.keys()),pTable,bitList)
    else:
        predict(c, ctx[1:], s-1, conList, pTable,bitList)
#******************************************************************************

#******************************************************************************
def loadModels(mFile):
    '''
        Loads the model definitions into the ppmd list
    '''
    nSeq = -1
    contextNum = None
    cChar = None
    res = None

    
    with open(mFile,'r') as rFile:
        for line in rFile:
            words = line.strip().split('\t')
            if 'ref' in words[0]:
                eD = dict()
                ppmd.append(eD)
                nSeq += 1
                
                ppmd[nSeq]['name'] = words[1]
                ppmd[nSeq]['type'] = words[2]
                ppmd[nSeq]['PC'] = words[3]
                ppmd[nSeq]['table'] = list()
                
                #print(ppmd)
                
            elif 'C' in words[0]:
                contextNum = int(words[1])
                tDict = dict() 
                ppmd[nSeq]['table'].append(tDict)
                
                #print(ppmd[0])
            
            elif 'K' in words[0]:
                if len(words) == 1:
                    cChar = ''
                else:
                    cChar = words[1]
                
                ppmd[nSeq]['table'][contextNum][cChar] = dict()
                
                #print(ppmd[0])
            
            elif 'R' in words[0]:
                res = words[1]
                count = float(words[2])
                
                ppmd[nSeq]['table'][contextNum][cChar][res] = count   
                
                #print(ppmd[0]['table'][1]) 
    
    #print(ppmd[0]['table'][1])            
#****************************************************

##***************************************************
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
            print('challenge: false', i, target, diff)
            return False
        
    return True
#***********************************************        
    
##**********************************************
def check_subtype(sLike,seqId):
    '''
        Uses COMET's desition tree to call subtypes for a query
    '''
    
    #print(sLike[0])
    #print(len(sLike),len(sLike[1]))        
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
        
        if ppmd[i]['PC'] == 'P': # reference is a PURE subtype
            if not psll: # psll yet to define
                psll = lsum
                PS = i
            elif psll < lsum:
                psll = lsum
                PS = i
    
    print(S,sll,PS,psll) 
    
    wSize = 100 # set window size
    bSize = 3 # set the step size
    
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
            rTxt = ppmd[S]['type'] + ' (PURE)'
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
            rTxt = ppmd[PS]['type'] + ' (check ' + ppmd[S]['type'] + ')'
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
                rTxt = ppmd[S]['type'] + ' (CRF)'
                return rTxt
            
                   
##***************************************************  

##***************************************************    
def calculateLogLikelihood(msg,seqId):
    '''
        Calculates likelihood values for a sequence
    '''
    sLike = list()
    
    dna = ['A','T','C','G']
        
    for r in range(len(ppmd)): #len(ppmd)
        bits = list() # holds likelihood values for each sites in a query
        for j in range(len(msg)):
            c = msg[j]
            if c not in 'ATCG':
                continue
            
            start = j - cSize if j > cSize else 0
            ctx = msg[start:j]
            
            if len(set(dna) | set(ctx)) <= 4: # no ambiguous characters
                bitList = list() # likelihood values for each context for each residue
                
                predict(c,ctx,len(ctx),[],ppmd[r]['table'],bitList)
                
                # combines the probability from each context into a single likelihood
                bits.append(sum(bitList)) 
            
        #print(bits)
        #print(r,sum(bits))
            
        ## create a list of likelihhod values
        ## for each query based on each reference
        ## for each residue position
        sLike.append(bits)
            
    #print(sLike[1][0:10])    
        
    return check_subtype(sLike,seqId)
    
##******************************************************

##***************************************************    
def subtype_calling(query):
    '''
        Reads in each of the query sequences
        Calculates the likelihoods for each positions for each reference
        predicts the subtype as: 
            - Pure or CRF
            - Pure with potential CRF
            - Unassigned
          
    '''
    seqs = list(SeqIO.parse(query,'fasta'))
    
    for i in range(len(seqs)):  #len(seqs)
        #print('id:',seqs[i].id)
        
        msg = str(seqs[i].seq).upper()
        
        tMsg = calculateLogLikelihood(msg,seqs[i].id)
        print(seqs[i].id, tMsg)
        
##******************************************    

#*******************************************
def deAlign(iFile, dFile):
    '''
      - Removes gaps (if any) from the input sequence file
    ''' 
  
    #print("\nRemoving the gap characters '-'/'.' from the sequences")
  
    # Reading sequence file in fasta format
    seqs = list(SeqIO.parse(iFile,'fasta'))
  
    if len(seqs) >= 1: # at least one sequence present | file is in FASTA format
  
        st = ''
  
        for seq in seqs:
          st += '>' + seq.id + '\n' + str(seq.seq).replace('-','').replace('.','') + '\n'
  
        fh = open(dFile,'w')
        fh.write(st)
        fh.close()
  
        msg = '[' + time.strftime('%d %b %H:%M:%S') + ']'
        msg += ' Gapless query sequence file written in <%s>\n' % dFile
        print(msg)
  
    else: # no sequence present or wrong format
        msg = '\n\nError: Could not read the input sequence file.'
        msg += '\n       Make sure the file is in FASTA format'
        msg += '\n       and at least one sequnce present in file\n'
        sys.exit(msg) 

#******************************************************************************

#******************************************************************************
def getArguments():
    '''
        Parse all the command line arguments from the user
    '''
    
    parser = argparse.ArgumentParser(description='Predicts subtype of a sequence based on PPMD models trained using reference sequences', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-q','--query',required=True,help='Query sequence file in FASTA format for subtype classification')
    parser.add_argument('-c','--context',type=contextThreshold,default=8,help='Context size for PPMD model (default: 8)')
    parser.add_argument('-m','--modelFile',required=True,help='PPMD model file for reference sequences')


    args = parser.parse_args()
    
    return args
    
##********************************************************##

##********************************************************##
if __name__=="__main__":
    args = getArguments()

    ## get context size
    cSize = int(args.context)
        
    dna = ['A','T','C','G']
    ## creates an empty list to load PPMD models
    ppmd = list()
    
    ## read in the model file and populate the ppmd list
    loadModels(args.modelFile)
    #print(ppmd[0]['table'][0:3])
    
    ## Dealign the query sequence
    deAlign(args.query,'query.dealign.fas')
    ## subtype calling
    subtype_calling('query.dealign.fas')
    
##***************************************************

#!/usr/bin/env python3

'''
    predictType_ppmd_train.py
    
    a python3 based program to create a PPMD model for subtype classification
     
    The program takes as input the following:
        - a reference sequence file in FASTA format
        - a fixed context size
        - an output file where the PPMD model will be stored
    
    This model can then be used to get log likelihoods of the probablities of 
    seeing a residue at a particular position given the context
    
    These likelihoods then can be used to predict the subtype
'''

import math, os, sys, argparse, time

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

##****************************************************
## update function creates the PPMD model for each
## of the training sequnce and add to the dictionary
##****************************************************
def update(c, ctx):
    for i in range(0,len(ctx)+1):
        pre = ctx[i:]
        s = len(pre)
        if pre not in table[s]:
            table[s][pre] = dict()
        if c not in table[s][pre]:
            table[s][pre][c] = 0.5
        else:
            table[s][pre][c] += 1

#*******************************************************

##********************************************************##
def getArguments():
    '''
        Parse all the command line arguments from the user
    '''
    
    parser = argparse.ArgumentParser(description='Creates a PPMD model from a set of reference sequences given a fixed context size for subtype classification', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r','--reference',required=True,help='Reference sequence file in FASTA format for creating PPMD model')
    parser.add_argument('-c','--context',required=True,type=contextThreshold,default=8,help='Context size for PPMD model (default: 8)')
    parser.add_argument('-m','--modelFile',required=True,help='Output file to save the PPMD model')


    args = parser.parse_args()
    
    return args
    
##********************************************************##

if __name__=="__main__":
    args = getArguments()
    
    #print(args.reference,args.context,args.modelFile)    
    
    ## get the context size
    cSize = int(args.context)
    
    ## create the alphabet
    dna = ['A','T','C','G']
    
    ## define the digits list
    digits = list('0123456789')
    
    ## read the reference sequence file
    seqs = list(SeqIO.parse(args.reference,'fasta'))
    
    ## create a dictionary for each of the reference sequences
    #ppmc = [dict() for i in range(len(seqs))]

    ## start writing the output model file
    fh = open(args.modelFile,'w')
    #fh.write('Models 2\n')
        
    ##**********************
    ## for each sequence create a PPMD model and update the dictionary
    ##***********************
    for i in range(len(seqs)): #len(seqs)
        sType = seqs[i].id.split('.')[0]
        
        #ppmd[i]['name'] = seqs[i].id
        #ppmd[i]['PC'] = 'C' if sType[0] in digits else 'P'
        #ppmd[i]['type'] = sType
        
        name = seqs[i].id
        PC = 'C' if sType[0] in digits else 'P'
        
        
        ## PPMD model for that reference sequence
        table = [dict() for j in range(cSize + 1)] # {0,1,...,k}
        
        ## convert sequence into a string
        sequence = str(seqs[i].seq).upper()
         
        for s in range(len(sequence)):
            c = sequence[s]
            start = s - cSize if s > cSize else 0
            ctx = sequence[start:s]
            if len(set(dna) | set(ctx)) <= 4: # ctx does not contain ambiguous characters
                update(c,ctx)
            
        ## add the PPMD model to the dictionary
        #ppmd[i]['table'] = table 
        
        
        fh.write('ref\t%s\t%s\t%s\n' % (name,sType,PC))
        
        
        for t in range(cSize + 1):
            fh.write('C\t%s\n' % t)
            #print('contextSize\t%d' % t)
            
            keys = list(table[t].keys())
            for key in keys:
                fh.write('K\t%s\n' % key)
                #print('contex\t%s' % key)
                
                chars = list(table[t][key].keys())
                for ch in chars:
                    fh.write('R\t%s\t%.1f\n' % (ch,table[t][key][ch]))
                    #print('char\t%s\t%d' % (ch,table[t][key][ch]))
           
         
    fh.close()
        
    #print(ppmc[0:2])
      
    
    

#!/usr/bin/env python3

'''
pyCOMET_train.py

A python3 based program to create PPMD models for subtype classification

The program takes as input the following:
    - reference sequence file in FASTA format
    - context size (default 8)
    - name of the output file where PPMD model files are written
'''

from Bio import SeqIO
from collections import Counter
import json

import math, os, sys, argparse, time

##********************************************************##
def contextThreshold(x):
    '''
    The context must be a positive integer
    '''
    try:
        x = int(x)
    except ValueError as e:
        sys.exit(e)

    if x <= 0:
        raise argparse.ArgumentTypeError('%r must be a positive integer' % (x,))

    return x
#*******************************************************

#******************************************************
def make_subtype_dict(refs,subtypes):
    '''
    This function creates a dictionary os the subtypes
        - keys of the dictionary are the subtype labels
        - Dictionary values are the reference sequences for that subtype
    '''
    # get the names of the distinct subtypes
    unique_subtypes=list(set(subtypes))
    # add subtypes entries in the dictionary
    refdict={s:[str(refs[i].seq) for i in range(len(refs)) if subtypes[i]==s] for s in unique_subtypes}
    # return the dictionary and the list of the unique subtypes
    return(refdict,unique_subtypes)
#******************************************************

#******************************************************
def make_refdict(fn,fmt="fasta",sep=".",field=0):
    '''
    This function calls make_subtype_dict() to create subtype dictionary
    '''
    # Reads in the reference sequence file
    refs=list(SeqIO.parse(fn,fmt))
    # extracts subtype labels from reference sequences
    subtypes=[r.description.split(sep)[field] for r in refs]
    # calls make_subtype_dict() to create dictionary
    refdict,unique_subtypes=make_subtype_dict(refs,subtypes)
    # returns the subtype dictionary to the calling function
    return(refs,refdict,unique_subtypes)
#******************************************************

#******************************************************
def ppmdTrain(c, ctx,sType):
    '''
    This function creates/updates the PPMD model for a given subtype from reference sequences

    Function arguments:
        - c: character
        - ctx: context of size 0..k
        - sType: subtype of the reference context
    This function builds a list of frequencies for each nucleotide given a context

    '''
    #for each context of size n get all the lower contexts
    # example context: abc, all contexts: abc, bc, c
    for i in range(0,len(ctx)+1):
        # get current context
        pre = ctx[i:]
        s = len(pre)
        # if the context is new; add to the dictionary
        # for that subtype
        if pre not in ppmd[sType][s]:
            ppmd[sType][s][pre] = dict()
            ppmd[sType][s][pre]['esc'] = 0
        # if character encountered first for a given context
        if c not in ppmd[sType][s][pre]:
            ppmd[sType][s][pre][c] = 0.5
            ppmd[sType][s][pre]['esc'] += 0.5
        # if character is seen before for a given context
        else:
            ppmd[sType][s][pre][c] += 1
#******************************************************

##********************************************************##
def getArguments():
    '''
        Parse all the command line arguments from the user
    '''

    parser = argparse.ArgumentParser(description='Creates PPMD models from a set of reference sequences given a fixed context size for subtype classification', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-r','--reference',required=True,help='Reference sequence file in FASTA format for creating PPMD models')
    parser.add_argument('-c','--context',type=contextThreshold,default=8,help='Context size for PPMD models (default: 8)')
    parser.add_argument('-m','--modelFile',required=True,help='Output file to save the PPMD model')
    parser.add_argument('-s','--sep',default='.',help='separator for header name (default=".")')
    parser.add_argument('-f','--field',type=int,default=0,help='Position of subtype name in the header (default = 0')


    args = parser.parse_args()

    return args

##********************************************************##

#***********************************************************

if __name__=="__main__":
    args = getArguments()

    ## get the context size
    cSize = int(args.context)

    ## create the alphabet
    dna = 'ACGT'

    ## define the digits list
    #digits = list('0123456789')

    ## get the name of the reference file
    refFileName = args.reference

    ## create the dictionary using reference file
    refSeqs,refdict,subtypes=make_refdict(refFileName,sep=args.sep,field=int(args.field))

    ## number of unique subtypes
    nSubtypes = len(subtypes)
    #print('Reference sequence file has the following subtypes:\n {}'.format(unique_subtypes))
    
    ## create an empty dictionary to hold PPMD models
    #ppmd = dict()
    ppmd = {s: [{} for k in range(cSize+1)] for s in subtypes}
    
    # Update the PPMD model one sequence at a time
    for seq in refSeqs:
        # extract the subtype of the reference sequence
        sType = seq.description.split(args.sep)[args.field]
    
        # convert the sequence into a string and remove gaps
        sequence = str(seq.seq).upper().replace('-','')
        # get the length of the sequences
        sLen = len(sequence)
    
        # Update the PPMD for each nucleotide sites
        for j in range(sLen):
            # get the nucleotide
            ch = sequence[j]
            if ch not in 'ACGT':
                continue
        
            # get the context
            start = j - cSize if j > cSize else 0
            context = sequence[start:j]
            # check if context has ambiguous residue
            if set(context).issubset('ACGT'):
                ppmdTrain(ch,context,sType)
    
    with open('pycomet.model.test.json','w') as fh:
        json.dump(ppmd,fh)    
#*************************************************************

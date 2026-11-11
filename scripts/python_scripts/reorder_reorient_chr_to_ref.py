#!/usr/bin/env python3

import os
import argparse
import logging
import pandas as pd
import json
import sys
from PowerhouseOfTheCell import *

# This script will take in the output files from gfastats and busco
# Input
# sample name
# fasta file for sample
# coords.idx file from DotPrep.py script
# number of chromosomes of the ref and query are assumed to be the same
# order of records in fasta file is assumed to be in ascending order. i.e the first chromosme is Chromosome 1, no matter how its ordered 

# Output
# fasta file with reoriented chromosoems
# excel file with remappings sense check

vars = argparse.ArgumentParser(description='Rename and Reorient Scaffolds In A Fasta File')

vars.add_argument('-s',dest='sample', type=str, required=True, help='sample(s) to process stats for')

vars.add_argument('--sample_fasta',dest='sampleFasta', type=str, required=True, help='fasta file for sample')

vars.add_argument('--ref_fasta',dest='refFasta', type=str, required=True, help='fasta file for reference')

vars.add_argument('-nhap',dest='hapNum', type=str, required=True, help='haplome number')

vars.add_argument('-c',dest='coordsIdxFile', type=str, required=True, help='path to coors idx file')

vars.add_argument('-o','--output_prefix',dest='outputPrefix', type=str, required=True, help='prefix for the output excel file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=True, help='name of directory for output. will be stored in current working directory if none given')

args = vars.parse_args()

# load in fasta and get ref dicts
sampleFastaRefDict = load_fasta(args.sampleFasta)
refFastaRefDict = load_fasta(args.refFasta)

# get number of chromosomes
nchr = len(sampleFastaRefDict)

# create dict of chromosome in reference, and its number in the chromosome ordering
refChrom_chromNumber = {j.split()[0]:i for i,j in enumerate(list(refFastaRefDict.keys()),start=1) }

# load in coors idx file at different row intervals in the file

# start index for the query df is after the #query line, and that comes after the #ref line. Since it is 0 based we start at 0 which is # ref then add the nchr, then add 3. 1 for the header line in # ref, 1 for the actual line # ref. and the other for the header # query line
queryDfStartIndex = nchr + 3

# we take nchr times 2 because dtaat in ref and query df contains nch. Then we add another 5. 1 for # ref line, 1 for header line in # ref, one for # query line, one for the # query header line. one for # overvirew
overviewDfStartIndex = (nchr*2) + 5

# load in dfs
queryDf = pd.read_csv(args.coordsIdxFile,sep=',',skiprows=queryDfStartIndex,nrows=nchr,index_col='query')
overivewDf = pd.read_csv(args.coordsIdxFile,sep=',',skiprows=overviewDfStartIndex)

# create dict for renamed samples
renamedSampleFastaRefDict = {}

# create sense check df
senseCheckRenaming = {'reference_Chr_name':[],
                      'query_Chr_original_name':[],
                      'query_Chr_new_name':[],
                      'performed_reverse_complement':[]}

# go through each header in our fasta sample, and check if it has to be reverse complemented
for header,record in sampleFastaRefDict.items():

    # get row of data
    row = queryDf.loc[header]

    # check orientation relative to ref
    if row['orientation'] == '-':

        # reverse complement the sequence
        record.sequence = find_reverse_complement_fasta(sequence=record.sequence)

        senseCheckRenaming['performed_reverse_complement'].append('no')

    elif row['orientation'] == '+':
        senseCheckRenaming['performed_reverse_complement'].append('yes')

    # next check what chromosome it has been renamed to 

    # for a given query, get only the unique alignments found
    unqiueAlignmentsOnlyDf = overivewDf[(overivewDf['query'] == header) & (overivewDf['tag'] == 'unique')]

    # sometimes there can be more than one ref chromosome with a unique alignment to the query so we want to get the ref chr with the most alignments to our query
    refChr_unqiueHitsCount = unqiueAlignmentsOnlyDf['ref'].value_counts().reset_index().set_index('ref')
    
    # get the ref with the most unique hits to the query
    referenceMatch = refChr_unqiueHitsCount['count'].idxmax()

    # get new name for the query
    newQueryChrName = f"{args.sample}_H{args.hapNum}_Chr{refChrom_chromNumber[referenceMatch]}"

    # set new output header for the record
    record.output_header = f">{newQueryChrName}\n"

    # set new output sequence
    record.output_sequence = create_multiline_fasta(sequence=record.sequence,characters_per_line=60)

    # build new outdict
    renamedSampleFastaRefDict[newQueryChrName] = record

    # build sensecheck
    senseCheckRenaming['reference_Chr_name'].append(referenceMatch)
    senseCheckRenaming['query_Chr_original_name'].append(header)
    senseCheckRenaming['query_Chr_new_name'].append(newQueryChrName)

# write output fasta
output_fasta_file(FastaRefDict=renamedSampleFastaRefDict,output_directory=args.output_directory,file_name=f'{args.sample}_Haplome{args.hapNum}.fa')

# create sense check df
senseCheckDf = pd.DataFrame(senseCheckRenaming)
write_to_tsv_and_excel(senseCheckDf,outputDir=args.output_directory,csv=False,sheetName='Chromosome Remappings',outputFilePrefix=args.outputPrefix)



    

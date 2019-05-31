# coding: utf-8
import numpy as np
import pandas as pd
import re
import os

import Bio
from Bio import motifs
from Bio import SeqIO

def read_fasta(fasta_fp):
    fasta_lines = open(fasta_fp).readlines()
    sequences_dct = dict()
    sequences_array = []
    if re.search(">",fasta_lines[0]) is None:
        print("This is not a valid fasta format")
    for line in fasta_lines:
        if line.strip():
            if re.search(">", line):
                sp = line.replace(">","").replace("\n","").strip()
            else:
                seq = line.strip("\n")
                sequences_array.append(seq)
                sequences_dct[sp] = seq
    sequences_array = np.array(sequences_array)
    return sequences_dct


def write_fasta(out_fp, alignments_dct):
    '''
    :param out_fp: output file path
    :param alignments_dct: A dictionary containing alignments. Keys are species names, and values are sequences
    :return: nothing
    '''
    string = ""
    for seq in alignments_dct:
        string += ">" + seq + "\n"
        string += alignments_dct[seq] + "\n"
    with open(out_fp, 'w') as g:
        g.write(string)
    return


def convert_sequences_to_array(sequences):
    '''
    inputs: sequence of nucleotides represented as a string composed of A, C, G, T
    outputs: a list of numpy array representations of a sequence with:
             A = [1, 0, 0, 0]
             C = [0, 1, 0, 0]
             G = [0, 0, 1, 0]
             T = [0, 0, 0, 1]
             N = [0.25,0.25,0.25,0.25]
             
    '''

    nucleotide_array_dict = {'A': [1, 0, 0, 0],
                             'C': [0, 1, 0, 0],
                             'G': [0, 0, 1, 0],
                             'T': [0, 0, 0, 1],
                             'N': [0.25,0.25,0.25,0.25]}
    
#     nucleotide_array_dict = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 0}
    sequence_array_list = []
    for seq in sequences:
        seq_array = []
        for nuc in seq:
            seq_array.append(nucleotide_array_dict[nuc])
        seq_array = np.array(seq_array)
        sequence_array_list.append(seq_array)
    return sequence_array_list


def pad_sequence_arrays(seq_arrays, pad_array, num, direction):
    '''
    inputs:
        seq_arrays: an array of sequences in form of matrix: 
            A = [1, 0, 0, 0]
            C = [0, 1, 0, 0]
            G = [0, 0, 1, 0]
            T = [0, 0, 0, 1]
            N = [0.25,0.25,0.25,0.25]
        pad_array: array to be padded to @seq_arrays
        num: number of times that @pad_array is padded
        direction: 
            'both' -- both ends
            'pre' -- pad to the front only
            'post' -- pad to the back only
    outputs: a list of numpy array representations of a sequence
             
    '''
    
    padded_seq = []
    for i in seq_arrays:
        if direction == 'both':
            padded_seq.append([pad_array]*num + i.tolist() + [pad_array]*num)
        elif direction == 'pre':
            padded_seq.append([pad_array]*num + i.tolist())
        elif direction == 'post':
            padded_seq.append(i.tolist() + [pad_array]*num)
        else:
            print('Wrong Parameter!!')
    return np.array(padded_seq)


def peakFile2bed(peakFile, OnlyRegion=True):
    file = open(peakFile, 'r')
    lines = file.readlines()
    file.close()
    for l in range(len(lines)):
        if lines[l][0] != '#':
            start_line = l-1
            break
    
    with open('tmp0.txt', 'w') as tmpf:
        for line in lines[start_line:]:
            tmpf.write(line)
    
    df = pd.read_csv('tmp0.txt', sep='\t')
    bed_file = peakFile.split('.')[0]+'.bed'
    print('BED file:', bed_file)
    if OnlyRegion:
        df.iloc[:,1:4].to_csv(bed_file, sep='\t', index=None, header=None)
    else:
        df.iloc[:,1:5].to_csv(bed_file, sep='\t', index=None, header=None)
    os.remove('tmp0.txt')

    
def load_motifs(motif_dir='/home/zes017/curated_motifs_jasparFormat/', key_id=0):
    motif_dict = {}
    nuc = ['A', 'C', 'G', 'T']
    for mf in os.listdir(motif_dir):
        with open(motif_dir + mf) as f:
            m = motifs.read(f, 'jaspar')
            counts = np.array([m.counts[n] for n in nuc])
            m.pseudocounts = np.mean(counts.sum(axis=0))//10
            m.background = None
            if key_id == 0:
                motif_dict[m.base_id] = m
            elif key_id == 1:
                motif_dict[m.name] = m
            
    return motif_dict
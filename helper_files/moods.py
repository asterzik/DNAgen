#!/usr/bin/env python
import MOODS.scan
import MOODS.tools
import MOODS.parsers

import sys
import time
import os
#file_name = 'data/chr1.fa'
#with open(file_name, "r") as file_handle:
#    first = file_handle.next().strip()
#    if first[0] == '>':
#        first = ''
#    seq = "".join([first] + [line.strip() for line in file_handle])

def motif_match(seq_array, matrix_path='./helper_files/motif_matrix/'):
    '''Takes sequence as a string and matrix directory path'''
        
    # ---- parameters ----
    pseudocount = 0.00001
    pvalue = 0.0001
    # ---- end parameters ----
    
    seq=''
    for i in seq_array:
        seq += i


    # ---- compute bg ----
    bg = MOODS.tools.bg_from_sequence_dna(seq,1)

    # ---- process matrices ----
    matrix_names = [filename for filename in os.listdir(matrix_path) if filename[-4:] == '.pfm']
    matrices = []
    for filename in matrix_names:
        matrices.append(MOODS.parsers.pfm_to_log_odds(matrix_path + filename, bg, pseudocount))
    # reverse complements
    #matrices = matrices + [MOODS.tools.reverse_complement(m) for m in matrices]
    rev_matrices=[]
    count = 0
    for m in matrices:
        count+=1
        rev_matrices.append(MOODS.tools.reverse_complement(m))
    matrices.extend(rev_matrices)
    thresholds = []
    for m in matrices:
        thresholds.append(MOODS.tools.threshold_from_p(m, bg, pvalue))

    # ---- scanning ----
    results = MOODS.scan.scan_dna(seq, matrices, bg, thresholds, 7)

    # ---- process results ----
    fr = results[:count]
    rr = results[count:]

    # mix the results together, use + and - to indicate strand
    results = [ [(r.pos, r.score, '+') for r in fr[i]] + [(r.pos, r.score, '-') for r in rr[i]] for i in range(count)]
    for (matrix,matrix_name,result) in zip(matrices, matrix_names, results):
        l = len(matrix[0])
        for r in sorted(result, key=lambda r: r[0]):
            strand = r[2]
            pos = r[0]
            hitseq = seq[pos:pos+l]
            # print matrix_name + '|' + str(pos) + '|' + strand + '|' + hitseq + '|'  + str(r[1])
    total = sum([len(r) for r in results])
    hits = []
    for i, name, m in zip(range(len(results)), matrix_names, results):
        hits.append(len(m))
    return hits, total


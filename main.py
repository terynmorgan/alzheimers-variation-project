#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 15:54:31 2021
AP
@author: terynmorgan
"""
#import modules
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq

seq_record = SeqIO.parse('seq.fasta.fa', 'fasta')
seq_list = []

# Extracts sequence from seq_Record handle
for record in seq_record:
    seq = record.seq
    
# Puts sequence into a list so it can be mutated later
for i in seq:
    seq_list.append(i)

# Predicted TF binding sites from PROMO
# Factor name, position
instances = [
    # MZF-1 [T00529] (0-6)
    Seq('AGGGGG')  
    # COE1 [T01112] (181-187)
    , Seq('TCCCTT')
    # HNF-3 [T02277]  (100-105)
    , Seq('TGTTTC')  
    # c-Ets-1 [T00112] (81-86)
    , Seq('CTGGAA')
    # c-Ets-2 [T01397] (191-197)
    , Seq('GCAGTT')
    # c-Myb [T00137] (190-196)
    , Seq('GGCAGT')
    # MED8 [T05391] (26-32)
    , Seq('ATTTCA')            
    # TFs in range of SNV
    # TMF [T00835] (123-128)
    , Seq('AAGTAT')
    # #IPF1 [T04362] (120-125)
    , Seq('TGTACG')
    # TGA1a [T00829] (121-126)
    , Seq('GTACGT')
]

# Creates motif object
motif_obj = motifs.create(instances)  
# Creates position weight matrix (PWM) from TFs
pwm = motif_obj.counts.normalize()    
# Converts PWM into PSSM
pssm = pwm.log_odds()         

# Mutates G at postion 125 to C
seq_list[125] = 'C'
# Put this strand into PROMO to obtain TFs
seq_list = ''.join([str(elem) for elem in seq_list])
# The TF sites will be the same for areas that don't include the 125 mutated position
mutated_instances = [
    # MZF-1 [T00529] (0-6)
    Seq('AGGGGG')  
    # COE1 [T01112] (181-187)
    , Seq('TCCCTT')
    # HNF-3 [T02277]  (100-105)
    , Seq('TGTTTC')  
    # c-Ets-1 [T00112] (81-86)
    , Seq('CTGGAA')
    # c-Ets-2 [T01397] (191-197)
    , Seq('GCAGTT')
    # c-Myb [T00137] (190-196)
    , Seq('GGCAGT')
    # MED8 [T05391] (26-32)
    , Seq('ATTTCA')            
    # TFs in range of SNP with mutation (125--> C)
    # TMF [T00835] (123-128)
    , Seq('AACTAT') 
    # IPF1 [T04362] (120-125)
    , Seq('TGTACC')
    # TGA1a [T00829] (121-126)
    , Seq('GTACCT')
]

# Creates motif object for mutated strand TFs
mutated_motif_obj = motifs.create(mutated_instances) 
# Position weight matrix from  TFs
mutated_pwm = m_mut.counts.normalize()       
# PWM to PSSM
mutated_pwm = mutated_pwm.log_odds()            

# Background value set to None resets to a uniform distribution
# Motif PSSM mean
mean = pssm.mean()                                    
mutuated_mean = pssm_mut.mean()     

# Motif PSSM standard deviation
std = pssm.std()
mutuated_std = pssm_mut.std()

# m.background to a single value will be interpreted as the GC content
motif_obj.background = 0.8
mean_background = pssm.mean(motif_obj.background)             #0.5842829802829755
mut_mean_background = pssm_mut.mean(motif_obj.background)     #0.5842829802829755

# Standard deviation with background value
std_background = pssm.std(motif_obj.background)               #1.2280435787449955
mutuated_std_background = pssm_mut.std(motif_obj.background)       #1.2280435787449955

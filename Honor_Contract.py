#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 15:54:31 2021
AP
@author: terynmorgan
"""
#Three main genes are considered the main risk factors for EOAD (Early Onset Alzheimer's Disease): 
#APP on chromosome 21q21.3, PSEN1 on 14q24.3, and PSEN2 on 1q42.13. 

#For this project, I identified the Missense Variant SNV on the PSEN1 gene chr14.73217225 (GRCh38.p12) through the NCBI gene software. 
#The alleles associated with this SNV are G>A/G>T. 
#Looking at rs661 in the Gene Variation Viewer on NCBI, I extracted the 100bp surrounding the SNV in the fasta file: 'seq.fasta.fa'. 


#import Bio
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
seq_Record=SeqIO.parse('seq.fasta.fa', 'fasta')
seq_list=[]
#gets sequence from seq_Record handlde
for record in seq_Record:
    seq=record.seq
#puts the sequnece into a list so it can be mutated later
#Seq object don't support item assignment
for i in seq:
    seq_list.append(i)

#predicted TF binding sites from PROMO
#JASPAR and TRANSFAC did not have PSEN1 in their database
                            #factor name, position
instances= [Seq('AGGGGG'),  #MZF-1 [T00529] (0-6)
            Seq('TCCCTT'),  #COE1 [T01112] (181-187)
            Seq('TGTTTC'),  #HNF-3 [T02277]  (100-105)
            Seq('CTGGAA'),  #c-Ets-1 [T00112] (81-86)
            Seq('GCAGTT'),  #c-Ets-2 [T01397] (191-197)
            Seq('GGCAGT'),  #c-Myb [T00137] (190-196)
            Seq('ATTTCA'),  #MED8 [T05391] (26-32)
#TFs in range of SNV
            Seq('AAGTAT'),  #TMF [T00835] (123-128)
            Seq('TGTACG'),  #IPF1 [T04362] (120-125)
            Seq('GTACGT')]  #TGA1a [T00829] (121-126)
m=motifs.create(instances)  #created motif object
pwm=m.counts.normalize()    #position weight matrix from TFs
pssm=pwm.log_odds()         #converts pwm into pssm

#mutates G at postion 125 to C
#position 125 supposed to be A/T, with SNV missense variation goes to G
#want to see score if goes to alternative of G
seq_list[125]='C'
#put this strand into PROMO to obtain TFs
listToStr= ''.join([str(elem) for elem in seq_list])
#the TF sites will be the same for areas that don't include the 125 mutated position
                                    #factor name, position
mutated_instances= [Seq('AGGGGG'),  #MZF-1 [T00529] (0-6)]
                    Seq('TCCCTT'),  #COE1 [T01112] (181-187)
                    Seq('TGTTTC'),  #HNF-3 [T02277]  (100-105)
                    Seq('CTGGAA'),  #c-Ets-1 [T00112] (81-86)
                    Seq('GCAGTT'),  #c-Ets-2 [T01397] (191-197)
                    Seq('GGCAGT'),  #c-Myb [T00137] (190-196)
                    Seq('ATTTCA'),  #MED8 [T05391] (26-32)
#TFs in range of SNP with mutation (125--> C)
                    Seq('AACTAT'),  #TMF [T00835] (123-128)
                    Seq('TGTACC'),  #IPF1 [T04362] (120-125)
                    Seq('GTACCT')]  #TGA1a [T00829] (121-126)
m_mut=motifs.create(mutated_instances) #created motif object for mutated strand TFs
pwm_mut=m_mut.counts.normalize()       #position weight matrix from  TFs
pssm_mut=pwm_mut.log_odds()            #pwm to pssm

#background value set to None resets to a uniform distribution
#motif PSSM mean
mean=pssm.mean()                                    
mut_mean=pssm_mut.mean()                           
#motif PSSM standard deviation
std=pssm.std()
mut_std=pssm_mut.std()
#m.background to a single value will be interpreted as the GC content
m.background=0.8
mean_background=pssm.mean(m.background)             #0.5842829802829755
mut_mean_background=pssm_mut.mean(m.background)     #0.5842829802829755
#standard deviation with background value
std_background=pssm.std(m.background)               #1.2280435787449955
mut_std_background=pssm_mut.std(m.background)       #1.2280435787449955

#Conclusion:
    #since not many TFs bind to the SNV itself AND only using 10 TFs,
    #when comparing before and after the mutation,
    #the slight variations in the pssm and pwm didn't account enough to make a significant diff in mean/std when using a background value
    #however, these variations can be seen looking at the motifs for mean/std w/o background value
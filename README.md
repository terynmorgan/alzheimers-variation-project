# Analyzing-Genetic-Variation-in-Alzheimers-Disease
Independent project done through INFO-B 210: Information Infrastructure outside of required coursework through Honors Contract with IUPUI Honors College.

Fall 2021
Programming Language: Python
Description: 
- Python script utilizes Biopython to compare and calculate the transcription binding factors of a single nucleotide polymorphism (SNP) in an Alzheimer's Disease gene before and after mutation

Background:
Three main genes are considered the main risk facotrs for Early Onset Alzheimer's Disease (EOAD): APP on chromosome 21q21.3, PSEN1 on 14q24.3, and PSEN2 on 1q42.13. For this project, I identified the Missense Variation single nucleotide polymorphism (SNP) on PSEN1 gene chr14.73217225 (GRCh38.p12) through NCBI gene software. A SNP is a genetic variation in single position within an individual's DNA sequence. A SNP can either occur as an replacement, insertion, or deletion of one of the four nucleotide bases (A, C, G, T). In the case of the PSEN1 gene, the missense variation mutate the A/T allele to G at position 125. 

The transcription binding score for a transcription sites comes from a Position Specific Scoring Matrix (PSSM). Transcription factor binding sites (TFBS) defined in the TRANSFAC and JASPAR database are used to construct specific binding site weight matrices for TFBS prediction for you. Since the PSEN1 gene wasn't in either of the aforementioned databases, I used PROMO to identity TFBS from an inputted DNA sequence. Using the top 10 TFBS motifs from PROMO, I converted each into a a Position Weight Matrix (PWM) and then into a PSSM. A sequence motif is a nucleotide or amino-acid sequence pattern of the same length. 

Steps of development:
1. Using NCBI SNP, identified the location of SNP and which gene it belongs to
2. Based on the gene information from NCBI RefSeq, identified which intron or exon the SNP was located in 
3. Extracted the 100 base pair sequence surrounding the SNP from RefSeq
4. Compared the extracted sequence with the transcription factor binding site to analyze the potential binding site
5. Calculated the binding score and identified the top 10 transcription factors that can potentially bind to the 100 base pair sequence 
6. Mutated the SNP (replaced G allele with C) and repeated step 4
7. For the top 10 transcription factors identified in step 5, calculated the binding score again and visualize binding score before and after mutation

Required Files:
seq.fasta.fa
INFO210_Honors_Contract.py

Conclusion:
Since not many TFs bind to the SNV itself and only using 10 TFs, when comparing before and after the mutation, the slight variations in the PSSM and PWM didn't account enough to make a significant different in mean and standard deviation when using a background value. However, these variations can be seen looking at the motifs for mean and standard deviation without background value.

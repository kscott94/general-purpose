""" Sere W. wants to know the histidine rations of all T. kodakarensis coding sequences"""

import pandas as pd
from Bio import SeqIO

fasta_dict = {}

with open('./data/tk_genes.fa', 'r') as fa:
    for record in SeqIO.parse(fa, "fasta"):
        gene_id = record.name
        sequence = ''.join(record.seq)
        fasta_dict[gene_id] = sequence

def translate(seq):
    genecode = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'stop', 'TAG': 'stop',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'stop', 'TGG': 'W'}
    seq=seq.upper()
    codons = []
    aa_seq = []
    reading_frame = 0
    for frame in range(int(len(seq)/3)):
        codons.append(seq[reading_frame:reading_frame+3])
        reading_frame += 3
    for codon in codons:
        aa_seq.append(genecode[codon])
    return(''.join(aa_seq[:-1]))

aa_dict = {}

for record in fasta_dict:
    aa_dict[record] = translate(fasta_dict[record])

aa_len = {}
aa_H_count = {}
aa_H_ratio = {}

for record in aa_dict:
    count = 0
    protein_len = len(aa_dict[record])
    aa_len[record] = protein_len
    for i in aa_dict[record]:
        if i == "H":
            count += 1
        aa_H_count[record] = count
    aa_H_ratio[record] = round(count/protein_len, 3)

dictionaries = [aa_dict, aa_len, aa_H_count, aa_H_ratio]

keys = []
sequences = []
lens = []
H_counts = []
H_ratios = []

for dict in dictionaries:
    for record in dict:
            if record not in keys:
                keys.append(record)
            else:
                continue

for key in keys:
    sequences.append(aa_dict[key])
for key in keys:
    lens.append(aa_len[key])
for key in keys:
    H_counts.append(aa_H_count[key])
for key in keys:
    H_ratios.append(aa_H_ratio[key])

complete_df = pd.DataFrame(list(zip(keys, sequences,lens,H_counts,H_ratios)), columns=['gene_id', "aa_seq", "aa_length", "H_count","H_ratio"]).sort_values("gene_id")
print(complete_df)
complete_df.to_csv("./data/Histidine_ratios", sep = '\t', index=False)

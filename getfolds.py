""" find patterns in RNA folding propensities amount all Thermococcus kodakarensis RNAs. Based on minimum free energy """


import pandas as pd

file = "./TS559_dotbracket.seq"

genes = []
sequences = []
struc = []

with open(file, "r") as file:
    file = file.readlines()
    name_counter = 0
    seq_counter = 1
    struc_counter = 2

    length = int(len(file)/3)
    while struc_counter <=length:
        for x in range(1,length):
            gene_id = file[name_counter].replace('\n', ""). replace(">","")
            sequence_fa = file[seq_counter].replace('\n', "")
            structure_dot = file[struc_counter].replace('\n', "")

            genes.append(gene_id)
            sequences.append(sequence_fa)
            struc.append(structure_dot)

            name_counter += 3
            seq_counter += 3
            struc_counter += 3

output = list(zip(genes, struc))
outfile = pd.DataFrame(output, columns=["gene_ID", "struc"])
outfile[["structure", "MFA"]] = outfile.struc.str.split(" ", n=1, expand=True)
outfile['MFA'] = outfile.MFA.str.replace("(", "")
outfile['MFA'] = outfile.MFA.str.replace(")", "")
outfile['length'] = outfile.structure.str.len()
outfile = outfile[['gene_ID','structure', 'MFA', 'length']]


outfile.to_csv("./output",index=False, sep="\t")

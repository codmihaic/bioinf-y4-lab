from Bio.Seq import Seq

dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print("DNA:", dna_seq)

rna_seq = dna_seq.transcribe()
print("RNA:", rna_seq)

protein_seq = rna_seq.translate()
print("Protein:", protein_seq)

rev_comp_seq = dna_seq.reverse_complement()
print("Reverse complement:", rev_comp_seq)

gc = (dna_seq.count("G") + dna_seq.count("C")) / len(dna_seq)
print("GC fraction:", round(gc, 3))

motif = "ATG"
positions = [
    i for i in range(len(dna_seq) - len(motif) + 1)
    if str(dna_seq[i:i + 3]) == motif
]
print("ATG positions:", positions)

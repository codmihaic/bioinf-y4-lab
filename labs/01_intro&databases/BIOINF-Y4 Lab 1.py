#Lab Activity 1:Accessing Biological Databases
#Set Up Entrez:

from Bio import Entrez
Entrez.email = "your.email@example.com"  # Replace with your email
#Note: NCBI requires an email address for usage tracking.

#Search the Database:
#Let's search for the BRCA1 gene in Homo sapiens.

handle = Entrez.esearch(db="nucleotide", term="BRCA1[Gene] AND Homo sapiens[Organism]")
record = Entrez.read(handle)
id_list = record["IdList"]
print(f"Found {len(id_list)} records.")

#Fetch a record
handle = Entrez.efetch(db="nucleotide", id=id_list[0], rettype="gb", retmode="text")
from Bio import SeqIO
record = SeqIO.read(handle, "genbank")
print(f"Sequence ID: {record.id}")
print(f"Description: {record.description}")


#Lab Activity 2: Working with DNA Sequences

#Importing Biopython Modules
#Open your Python environment and import the necessary modules:
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Creating a DNA Sequence Object
#Let's create a DNA sequence:

dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
print("DNA Sequence:", dna_seq)

#Transcription
#Transcribe the DNA to RNA:
rna_seq = dna_seq.transcribe()
print("RNA Sequence:", rna_seq)

#Translation
#Translate the RNA to a protein sequence:
protein_seq = rna_seq.translate()
print("Protein Sequence:", protein_seq)

#Reverse Complement
#Get the reverse complement of the DNA sequence:
rev_comp_seq = dna_seq.reverse_complement()
print("Reverse Complement:", rev_comp_seq)

#Analyzing the Sequence
#Calculate the GC content:
from Bio.SeqUtils import GC
gc_content = GC(dna_seq)
print("GC Content:", gc_content)
#Finding Motifs
motif = "ATG"  # Start codon
positions = [pos for pos in range(len(record.seq)) if record.seq.startswith(motif, pos)]
print(f"Motif '{motif}' found at positions: {positions}")

#Lab Activity 3: Exploring Genetic Variation
#Let's investigate known mutations in the BRCA1 gene.
#Accessing SNP Data
#We'll use dbSNP to retrieve information about single nucleotide polymorphisms.
    #Fetching SNP Data
handle = Entrez.esearch(db="snp", term="BRCA1[gene] AND Homo sapiens[organism]")
record = Entrez.read(handle)
snp_ids = record["IdList"]
print(f"Found {len(snp_ids)} SNPs associated with BRCA1.")
    #Fetching SNP Details
snp_record <- entrez_fetch(db="snp", id=snp_ids[1], rettype="docset", retmode="text")
cat(snp_record)


#Lab activity 4: Working with multiple sequences
#Download a multi-FASTA file containing mitochondrial sequences from different species.
#Parse the file and calculate GC content for each sequence.
#Compare the GC content across species and discuss any patterns observed.
records = list(SeqIO.parse("mitochondrial_sequences.fasta", "fasta"))
for record in records:
    gc_content = GC(record.seq)
    print(f"{record.id}: GC Content = {gc_content:.2f}%")




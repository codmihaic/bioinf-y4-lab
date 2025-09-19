from Bio import Entrez

# Always provide your email when using Entrez
Entrez.email = "your.email@example.com"

# Search for BRCA1 gene in the human genome
handle = Entrez.esearch(db="nucleotide", term="BRCA1[Gene] AND Homo sapiens[Organism]")
record = Entrez.read(handle)

# Get list of IDs returned by the search
id_list = record["IdList"]
print(f"Found {len(id_list)} records for BRCA1.")
print(f"First record ID: {id_list[0]}")

# Fetch the sequence using the ID
handle = Entrez.efetch(db="nucleotide", id=id_list[0], rettype="fasta", retmode="text")
fasta_data = handle.read()

# Print the sequence in FASTA format
print(fasta_data)

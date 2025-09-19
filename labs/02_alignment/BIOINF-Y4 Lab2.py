Introduction to Sequence Alignment
Types of Sequence Alignment
Global Alignment: Needleman-Wunsch Algorithm
Local Alignment: Smith-Waterman Algorithm
Substitution Matrices: PAM and BLOSUM
Applications of Sequence Alignment
Lab activities: Practical Sequence Alignment
    Using BLAST for Local Alignment
    Performing Global and Multiple Sequence Alignment with Clustal Omega
    Retrieving and Analyzing Real Genomic Data
    Perform assignments 2.1 -> 2.10, compare your results with Blast and Clustal Omega
Next week: phylogenetics and phylogenic trees

Lab Activity 1: Using BLAST for Local Alignment
Step 1: Accessing NCBI BLAST

Navigate to the NCBI BLAST website.: https://blast.ncbi.nlm.nih.gov/Blast.cgi
Choose the nucleotide BLAST option.
Step 2: Inputting a Query Sequence

In the text box, input the following sequence (which is a portion of the human hemoglobin gene):

ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG

Click BLAST to initiate the search.

Step 3: Analyzing the Results

Once BLAST completes, you’ll see a list of similar sequences from various organisms. BLAST highlights the regions of alignment, showing where your query sequence aligns with the database sequences.
Discussion:
What does the alignment tell you about the similarity between these sequences?
How significant are the results (look at the E-values)?


Lab Activity 2: Performing Global and Multiple Sequence Alignment with Clustal Omega
Step 1: Installing Clustal Omega

If you don’t already have Clustal Omega installed, you can install it via the command line:

sudo apt-get install clustalo

Step 2: Prepare Your Input Sequences

Let’s work with the following three sequences, which are related protein sequences:

>Sequence 1
MAIVMGRWKGAR...
>Sequence 2
MTIVMGRTWNGAA...
>Sequence 3
MAIVMGRWGGH...
Save these sequences to a file named sequences.fasta.

Step 3: Running Clustal Omega

In your terminal, run the following command to perform a multiple sequence alignment:

clustalo -i sequences.fasta -o aligned.fasta --force

This will produce an aligned sequence file named aligned.fasta.

Step 4: Analyzing the Results

Open the aligned file and observe the conserved regions. Clustal Omega marks regions where the sequences are identical.
Discussion:
What regions are conserved across all sequences?
Why do you think certain regions show variation while others are conserved?

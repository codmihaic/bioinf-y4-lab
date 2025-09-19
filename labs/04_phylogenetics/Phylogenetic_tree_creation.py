#We’ll be working with the sequence.fasta file from last week, building a phylogenetic tree to analyze evolutionary relationships.

#Step 1: Loading DNA Sequence Data

from Bio import SeqIO
sequences = list(SeqIO.parse('sequence.fasta', 'fasta'))
for seq in sequences:
    print(seq.id, seq.seq)
#This code displays each sequence ID and the sequence itself, giving us a sense of what we’re working with.

#Step 2: Sequence Alignment and Distance Matrix

from Bio.Phylo.TreeConstruction import DistanceCalculator
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(sequences)
print(distance_matrix)
#The distance matrix shows how similar or different each sequence is. Higher values indicate higher similarity.

#Step 3: Building a Neighbor-Joining Tree
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)
#This code creates an unrooted Neighbor-Joining tree, focusing on relationships without an assumed origin point.

#Step 4: Tree Visualization
from Bio import Phylo
import matplotlib.pyplot as plt

Phylo.draw(tree)
plt.show()
#The resulting tree shows how each sequence is related. Closer branches indicate greater similarity, which can suggest recent evolutionary divergence

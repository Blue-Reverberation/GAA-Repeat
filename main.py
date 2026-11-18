import re
import Bio
from Bio import SeqIO

test_data = SeqIO.read("/home/nerd/Documents/College/Python/GAA-Repeat/test-data/pure-run.fasta", "fasta") # Import data
sequence = str(test_data.seq) # Convert relevant object to string
#print(sequence) # debug
#print(type(sequence)) # debug


def indexed_match(motif): # Function to find motifs & locations
    for m in re.finditer(motif, sequence): # For each match of a given motif in the sequence,
        print('%02d-%02d: %s' % (m.start(), m.end(), m.group(0))) # Print the location of the match alongside the matched string




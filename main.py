import re
import Bio
from Bio import SeqIO

## Todo
# 1. 
#
#

## Flowchart
# Take input sequence, find pure runs and list coordinates
# Append coordinates if distance between them < max_gap
# 


test_data = SeqIO.read("/home/nerd/Documents/College/Python/GAA-Repeat/test-data/10-gap.fasta", "fasta") # Import data
test_sequence = str(test_data.seq) # Convert relevant object to string


def indexed_match(motif, sequence): # Function to find motifs & locations
    xlist = [] # Initialise empty list
    for m in re.finditer(motif, sequence): # For each match of a given motif in the sequence,
       #print('%02d-%02d: %s' % (m.start(), m.end(), m.group(0))) # Print the location of the match alongside the matched string
        xlist.append((m.start(), m.end())) # Append the start and end coordinates of each found pure run
    return xlist

match_list = indexed_match("(?:GAA)+", test_sequence)

#print(match_list) # [(0,15), (25,40)]
longer_sample = [(5,10), (15,20), (30,40), (50,60), (70,80), (90,100), (110,120), (130,140), (150,160), (170,180), (800,900)]


def typeshi(runs, max_gap):
    tracts = []
    if not runs: # check for empty lit
        return []
    
    curr_start, curr_end = runs[0]
    print("[0] start, end:", curr_start, curr_end)

    for start, end in runs[1:]:
        print("[>0] start, end:", start, end)
        gap = start - curr_end
        
        if gap <= max_gap:
            curr_end = end
            print("current end:", curr_end)
            print("current gap:", gap)
        else:
            tracts.append(((curr_start),(curr_end)))
            print("Current tract:",tracts)
            curr_end = end
            curr_start = start
            print("current start/end after gap:", curr_start, curr_end)

    tracts.append(((curr_start),(curr_end)))
    return tracts

        
        

test = typeshi(longer_sample, 20)

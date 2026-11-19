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

## Todo
# 1. Record span as well as start and end coordinates
# 2. Record interruption coordinates and save to a separate list
# 3. Count number of units in each tract
# 4. Record longest running pure run in each tract
# 5. Create dictionary storing this data and return that rather than tracts list
def build_tract_dict(runs, sequence):
    


def join_units(runs, max_gap):
    tracts = []
    tract_dict = [
    ]

    if not runs: # Check for empty list
        return []
    
    curr_start, curr_end = runs[0] # For first entry, set current start and end coordinates to the first entry in the list
    print("[0] start, end:", curr_start, curr_end) # Test print to check that the first entry is being read correctly

    for start, end in runs[1:]: # For each subsequent entry in the list, starting from the second entry
        print("[>0] start, end:", start, end) 
        gap = start - curr_end # Calculate the gap between the gap between the previous entry and the current entry
        
        if gap <= max_gap: # If the gap is less than or equal to the maximum gap, append the current end coordinate to the current tract
            
            curr_end = end # Update the current end coordinate to the end coordinate of the current entry
            print("current end:", curr_end) # Test print to check that the current end coordinate is being updated correctly
            print("current gap:", gap) # Test print to check that the gap is being calculated correctly
        else:
            tracts.append(((curr_start),(curr_end))) # Otherwise, append the current run to the list of tracts
            print("Current tract:",tracts) # Test print to check that the current tract is being appended correctly
            curr_end = end # And reset the current start and end coordinates to the current entry
            curr_start = start
            print("current start/end after gap:", curr_start, curr_end) #Test print to check that the current start and end coordinates are being reset correctly
    tracts.append(((curr_start),(curr_end))) # Append the final run to the list of tracts
    return tracts # Return the list of tracts

        
        

test = join_units(longer_sample, 10)
print(test) 

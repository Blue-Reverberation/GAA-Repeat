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
longer_sample = [(5,10), (15,20), (30,40)]

def coord_distance(listx):  
    count = 1 
    distance_list = []
    if len(listx) > 1:
        for i in listx: # For each index in matchlist, take n+1[0] - n[1]
            distances = (int(listx[count][0]) - listx[count-1][1])
            distance_list.append(distances)     
            if count < len(listx) - 1: # iterate
                count += 1
            else:
                return distance_list
    else:
        return(listx)
                

coordinate_distances = coord_distance(longer_sample)

#print(coordinate_distances)

def listjoiner(match_list, distance_list, max_gap):
    
    n = 0
    tracts = []
    is_same = []
    if len(match_list) > 1:
        for i in distance_list: # For each item in the distances
            if n == 0:
                if distance_list[n] < max_gap:#print(match_list[n][0],match_list[n+1][1])  # Combine the matches together
                    tracts.append((match_list[n][0],match_list[n+1][1]))
                    n += 1
                    print(tracts)
            elif distance_list[n] < max_gap:# If the distance is below the threshold
                    #print(match_list[n][1],match_list[n+1][1])  # Combine the matches together
                    tracts.append((match_list[n][1],match_list[n+1][1]))
                    n += 1
                    print(tracts)
            else:
                n+1
        if tracts == is_same:
            return tracts
        else:
            print("repeat. tracts is:", tracts)
            print("repeat. is_same is:", is_same)
            is_same = listjoiner(tracts, coord_distance(tracts), max_gap)
            return listjoiner(tracts, coord_distance(tracts), max_gap)
    else:
        print("returning match list", match_list)
        return match_list

test = listjoiner(longer_sample, coordinate_distances, 20)


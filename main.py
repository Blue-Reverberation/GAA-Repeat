import re
import Bio
import pathlib
from Bio import SeqIO


def indexed_match(motif, sequence): # Function to find motifs & locations
    xlist = [] # Initialise empty list
    for m in re.finditer(motif, sequence): # For each match of a given motif in the sequence,
       #print('%02d-%02d: %s' % (m.start(), m.end(), m.group(0))) # Print the location of the match alongside the matched string
        xlist.append((m.start(), m.end())) # Append the start and end coordinates of each found pure run
    return xlist

def tract_builder(runs, sequence): # Builds tracts into a dictionary from list of runs
    interruptions = [] #Initialise empty list
    for (_, previous_end), (next_start, _) in zip(runs, runs[1:]): # Unpacking relevant start/end data in neighbouring runs

        length = next_start - previous_end # Calculate interruption length

        # Defining what should go inside the interruptions dictionary
        interruptions.append({
            "start": previous_end,
            "end": next_start,
            "length": length,
            "sequence": sequence[previous_end:next_start],
            "multiple_of_3": length % 3 == 0,
        })

        start, end = runs[0][0], runs[-1][1] # Start and end of tract
        pure_units = sum((e-s)//3 for s, e in runs) # Calculate the number of units

    return {
            "start": start,
            "end": end,
            "span": end - start,
            "units": pure_units, # uninterrupted
            "longest_run": max((e-s)//3 for s,e in runs), # The longest uninterrupted run
            "interruptions": interruptions,
        }


def join_units(runs, max_gap, min_units, sequence): 
    # Function which joins pure runs if within parameters, then calls tract_builder() to create a list of dictionaries

    ## Todo
    # 1. Rigorously check edge cases

    if not runs:
        return []

    tracts_prefilter = [] # Initialise empty list for dictionaries
    tracts_postfilter = [] # Initialise empty list for filtered dictionaries containing only entries with units > min_units
    # Seeding the first tract with the first entry, as it has no previous run to compare 
    # current is a list holding runs that are being built into a tract
    current = [runs[0]] 
    
    for run in runs[1:]:
        gap = run[0] - current[-1][1] # Calculate the distance between the candidate pure run to the most recent addition to the tract being built
        if gap <= max_gap: # If the gap is within threshold,
            current.append(run) # Append the run to the tract being built
        else: # If the candidate pure runs distance exceeds the threshold,
            tracts_prefilter.append(tract_builder(current, sequence)) # Finish the tract and build into dictionary with tract_builder()
            current = [run] # Overwrite the old, now built track with the current run

        
    tracts_prefilter.append(tract_builder(current, sequence)) # Append the last track

    for entry in tracts_prefilter: # For each entry in the compiled dictionaries,
        if entry["units"] >= min_units: # if the units entry is greater than the minimum threshold
            tracts_postfilter.append(entry) # append the entry to the filtered list. Non appended items get discarded.
            
    return tracts_postfilter


if __name__ == "__main__": # Guarding

    test_data = SeqIO.read("test-data/multi-interrupt-test1.fasta", "fasta") # Import data
    test_sequence = str(test_data.seq) # Convert relevant object to string

    match_list = indexed_match("(?:GAA)+", test_sequence)

    test = join_units(match_list, 100, 3, test_sequence)
    print(test)
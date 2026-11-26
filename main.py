import re
import Bio
import pathlib
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def indexed_match(motif, sequence): # Function to find motifs & locations
    '''Returns list of (start, end) coordinates, greedy algorithm '''
    #print("index motif:", motif)
    reglist = [] # Initialise empty list
    for m in re.finditer(motif, sequence): # For each match of a given motif in the sequence,
       #print('%02d-%02d: %s' % (m.start(), m.end(), m.group(0))) # Print the location of the match alongside the matched string
        reglist.append((m.start(), m.end())) # Append the start and end coordinates of each found pure run

    return reglist

def tract_builder(runs, sequence, is_complement): # Builds tracts into a dictionary from list of runs
    '''Builds tract information dictionary from a single tract, including a nested interruption dictionary for a tract if it includes one'''
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
            "units": pure_units, 
            "longest_run": max((e-s)//3 for s,e in runs), # The longest uninterrupted run
            "interruptions": interruptions,
            "reverse_complement": is_complement
        }


def join_units(runs, max_gap, min_units, sequence, is_complement): 
    '''Joins pure runs if within max_gap parameters. Calls tract_builder() to create a tract information dictionary. Returns a list of each tract dictionary.'''
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
            tracts_prefilter.append(tract_builder(current, sequence, is_complement)) # Finish the tract and build into dictionary with tract_builder()
            current = [run] # Overwrite the old, now built track with the current run

    tracts_prefilter.append(tract_builder(current, sequence, is_complement)) # Append the last track

    for entry in tracts_prefilter: # For each entry in the compiled dictionaries,
        if entry["units"] >= min_units: # if the units entry is greater than the minimum threshold
            tracts_postfilter.append(entry) # append the entry to the filtered list. Non appended items get discarded.
            
    return tracts_postfilter

def parse_args(argv=None): # Argument parser 
        parser = argparse.ArgumentParser(description="Find repeats of motifs in a sequence")
        parser.add_argument("fasta", help="Path to FASTA file with a single sequence")
        parser.add_argument("--motif", default="GAA", help="repeat motif (default = GAA)")
        parser.add_argument("--max-gap", type=int, default=9, help="Largest gap (in bases) which is still accepted as joined motifs (default = 9)")
        parser.add_argument("--min-unit", type=int, default=3, help="Minimum number of units in a tract to be reported (default = 3)")
        args = parser.parse_args(argv)

        args.motif = args.motif.upper() # Make inputs upper case

        if not args.motif or set(args.motif) - set("ACTG"): # Using set logic (empty set = falsy) to determine if input contains non-base characters
            parser.error("--motif must only contain 'A', 'C', 'T' or 'G'") 

        if args.max_gap < 0: # Check for negative
            parser.error("--max-gap must be at least 0")

        if args.min_unit < 1: # Check for 0 or negative
            parser.error("--min-unit must be at least 1")

        return args

def clean_print(list_of_tracts):
    '''Prints information from a list of tract dictionaries cleanly in terminal'''
    span_data = []
    unit_data = []
    counter = 1
    
    for tract in list_of_tracts:
        print("TRACT #" + str(counter))
        counter += 1

        curr_start = tract["start"]
        print("start:" ,curr_start)

        curr_end = tract["end"]
        print("end:" ,curr_end)

        curr_span = tract["span"]
        print("span:" ,curr_span)
        span_data.append(curr_span)

        curr_units = tract["units"]
        print("units:" ,curr_units)
        unit_data.append(curr_units)

        curr_longest_run = tract["longest_run"]
        print("longest pure run:" ,curr_longest_run)

        curr_reverse_complement = tract["reverse_complement"]
        print("reverse complement?:" ,curr_reverse_complement, "\n")

        if not tract["interruptions"] == []:
            interruptions = tract["interruptions"]
            print("Joined interruptions present: \n")

            for interruption in interruptions:
                print("Interruption start:" ,interruption["start"])
                print("Interruption end:", interruption["end"])
                print("Interruption span:", interruption["length"])
                print("sequence:", interruption["sequence"])
                print("multiple of 3?:", interruption["multiple_of_3"], "\n")
            print("\n")
            #print(tract["interruptions"], "\n")
        else:
            print("No joined interruptions", "\n \n")

    #print("span data:", span_data)
    #print("unit data:", unit_data)
    return span_data, unit_data


if __name__ == "__main__": # Guarding

    args = parse_args()

    test_data = SeqIO.read(args.fasta, "fasta") # Import data
    test_sequence = str(test_data.seq)          # Convert relevant information to string
    #print("sequence:", test_sequence)


    reg_pure_runs = indexed_match("(?:"+args.motif+")+", test_sequence)   # Create list of pure runs
    #print("runs:", pure_runs)
    #print("args.motif:", args.motif)

    reg_tracts = join_units(reg_pure_runs, args.max_gap, args.min_unit, test_sequence, False) # Create list of dictionaries
    #print("reg_tracts:", reg_tracts)

    # Repeat steps for complementary motif
    reverse_comp_motif = str(Seq(args.motif).reverse_complement())    
    #print("complementary motif:", reverse_comp_motif)
    #print("equal?:", reverse_comp_motif == args.motif)

    if not reverse_comp_motif == args.motif: # Checking for motifs which are their own reverse complement
        revcomp_pure_runs = indexed_match("(?:"+reverse_comp_motif+")+", test_sequence)
        #print("revcomp_pure_runs:", revcomp_pure_runs)

        revcomp_tracts = join_units(revcomp_pure_runs, args.max_gap, args.min_unit, test_sequence, True)

    #print("revcomp_tracts:", revcomp_tracts) 


    reg_span_data = [] # Initialise empty lists to hold data from dictionary – Will be used for histograms
    reg_unit_data = []

    reg_span_data, reg_unit_data = clean_print(reg_tracts) # Using clean_print to append data
    # print("assigned span data:", reg_span_data)
    # print("assigned unit data;", reg_unit_data)

    # Repeat for reverse complement
    revcomp_span_data = []
    revcomp_unit_data = []
    
    revcomp_span_data, revcomp_unit_data = clean_print(revcomp_tracts)
    # print("assigned revcomp span data:", revcomp_span_data)
    # print("assigned revcomp unit data:", revcomp_unit_data)

 

## Todo
# 1. Distinguish between revcomp and reg in dictionary   # Done
# 2. Testing!   # Almost done
# 3. README.md  # WIP
# 4. Histogram for data     # Unfinished
    # i) Pure run Length
    # ii) Interruption length
# 5. Printing data cleanly  # Done!
# 6. Automatic pytest   # In progress

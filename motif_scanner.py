### COPYRIGHT ###
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.








import re
import Bio
import pathlib
import argparse
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq



def indexed_match(motif, sequence): # Function to find motifs & locations
    '''Returns list of (start, end) coordinates, greedy algorithm '''
    #print("index motif:", motif)
    reglist = [] # Initialise empty list
    for m in re.finditer(motif, sequence): # For each match of a given motif in the sequence,
        reglist.append((m.start(), m.end())) # Append the start and end coordinates of each found pure run

    return reglist

def tract_builder(runs, sequence, is_complement, motif): # Builds tracts into a dictionary from list of runs
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
    pure_units = sum((e-s)//len(motif) for s, e in runs) # Calculate the number of units

    return {
            "start": start,
            "end": end,
            "span": end - start,
            "units": pure_units, 
            "longest_run": max((e-s)//len(motif) for s,e in runs), # The longest uninterrupted run
            "interruptions": interruptions,
            "reverse_complement": is_complement
        }


def join_units(runs, max_gap, min_units, sequence, is_complement, motif): 
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
            tracts_prefilter.append(tract_builder(current, sequence, is_complement, motif)) # Finish the tract and build into dictionary with tract_builder()
            current = [run] # Overwrite the old, now built track with the current run

    tracts_prefilter.append(tract_builder(current, sequence, is_complement, motif)) # Append the last track

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
        parser.add_argument("--show-interruptions", action="store_true", help="Show full list of all joined interruptions")
        parser.add_argument("--plot", action="store_true", help="Plot histogram of pure rune length & interruption length")
        args = parser.parse_args(argv)

        args.motif = args.motif.upper() # Make inputs upper case

        if not args.motif or set(args.motif) - set("ACTG"): # Using set logic (empty set = falsy) to determine if input contains non-base characters
            parser.error("--motif must only contain 'A', 'C', 'T' or 'G'") 

        if args.max_gap < 0: # Check for negative
            parser.error("--max-gap must be at least 0")

        if args.min_unit < 1: # Check for 0 or negative
            parser.error("--min-unit must be at least 1")

        return args

def clean_print(list_of_tracts, show_interruptions):
    '''Prints information from a list of tract dictionaries cleanly in terminal'''
    span_data = []
    unit_data = []
    pure_run_length_data = []
    interruption_length_data = []
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
            print("-Joined interruptions present- \n")
            if show_interruptions == True:
                interruptions = tract["interruptions"]

                for interruption in interruptions:
                    print("Interruption start:" ,interruption["start"])
                    print("Interruption end:", interruption["end"])
                    print("Interruption span:", interruption["length"])
                    interruption_length_data.append(interruption["length"])
                    print("sequence:", interruption["sequence"])
                    print("multiple of 3?:", interruption["multiple_of_3"], "\n")
                print("\n")
                #print(tract["interruptions"], "\n")
        else:
            print("No joined interruptions", "\n \n")

    #print("span data:", span_data)
    #print("unit data:", unit_data)
    return span_data, unit_data, interruption_length_data


def plot_results(args_plot, reg_units_per_tract, revcomp_units_per_tract, reg_interruption_length_data, revcomp_interruption_length_data):
    if args_plot == True:
        import matplotlib
        matplotlib.use("Agg")

        repeat_units_per_tract = reg_units_per_tract + revcomp_units_per_tract
        if not repeat_units_per_tract == []:
            plt.hist(repeat_units_per_tract,
                    bins=[x - 0.5 for x in range(min(repeat_units_per_tract), max(repeat_units_per_tract) + 2)])
            plt.xlabel("Repeat units per tract")
            plt.ylabel("Number of runs")
            plt.savefig("run_lengths.png", dpi=150)
            plt.close()
            print("Saved run_lengths.png to project directory")

        interruption_length_data = reg_interruption_length_data + revcomp_interruption_length_data
        if not interruption_length_data == []:
            plt.hist(interruption_length_data,
                    bins=[x - 0.5 for x in range(min(interruption_length_data), max(interruption_length_data) + 2)])
            plt.xlabel("Interruption length in bases")
            plt.ylabel("Number of interruptions")
            plt.savefig("interruption_lengths.png", dpi=150)
            plt.close()
            print("Saved interruption_lengths.png to project directory")


def import_fasta(path):
    file_data = SeqIO.read(path, "fasta") # Import data
    sequence = str(file_data.seq).upper()         # Convert relevant information to string
    return sequence

def scan_print_reg_motif(sequence, motif):
    reg_pure_runs = indexed_match("(?:"+motif+")+", sequence)   # Create list of pure runs
    reg_tracts = join_units(reg_pure_runs, args.max_gap, args.min_unit, sequence, False, motif) # Create list of tract dictionaries

    reg_span_data = [] # Initialise empty lists to hold data from dictionary – Will be used for histograms
    reg_units_per_tract = []

    reg_span_data, reg_units_per_tract, reg_interruption_length_data = clean_print(reg_tracts, args.show_interruptions) # Using clean_print to append data
    return reg_units_per_tract, reg_interruption_length_data


def scan_print_revcomp_motif(sequence, motif):
    reverse_comp_motif = str(Seq(motif).reverse_complement().upper())  
    revcomp_pure_run_length_data = [] # Initialise empty lists
    revcomp_span_data = []
    revcomp_units_per_tract = []
    
    if not reverse_comp_motif == motif: # Checking for motifs which are their own reverse complement
        revcomp_pure_runs = indexed_match("(?:"+reverse_comp_motif+")+", sequence)
        revcomp_tracts = join_units(revcomp_pure_runs, args.max_gap, args.min_unit, sequence, True, reverse_comp_motif)

        revcomp_span_data, revcomp_units_per_tract, revcomp_interruption_length_data = clean_print(revcomp_tracts, args.show_interruptions)
    return revcomp_units_per_tract, revcomp_interruption_length_data

def main():
    sequence = import_fasta(args.fasta)
    reg_units_per_tract, reg_interruption_length_data = scan_print_reg_motif(sequence, args.motif)
    revcomp_units_per_tract, revcomp_interruption_length_data = scan_print_revcomp_motif(sequence, args.motif)

    plot_results(args.plot, reg_units_per_tract, revcomp_units_per_tract, reg_interruption_length_data, revcomp_interruption_length_data)

if __name__ == "__main__": # Guarding

    args = parse_args()
    main()
    


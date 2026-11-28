from motif_scanner import indexed_match, join_units, import_fasta, scan_print_reg_motif, scan_print_revcomp_motif
import Bio
from Bio import SeqIO
from Bio.Seq import Seq


## Empty Template - join_units()
# "start": ,
# "end": ,
# "span": ,
# "units": , 
# "longest_run": , 
# "interruptions": 
#     [
#     {'start': , 
#     'end': ,
#     'length': , 
#     'sequence': , 
#     'multiple_of_3': }
#     ],
# "reverse_complement": 
# "reverse_complement": False
def test_indexed_match_finds_two_runs():
    assert indexed_match("(?:GAA)+", "GAAGAATTGAA") == [(0, 6), (8, 11)]

def test_tract_kept_when_units_equal_min_units():
    tracts = join_units([(0, 9)], 9, 3, "GAAGAAGAA", False, "GAA")
    assert len(tracts) == 1
    assert tracts[0]["units"] == 3

def test_no_match():
    assert indexed_match("(?:GAA)+", "TTTTTTTTTT") == []

def test_empty_runs():
    empty = indexed_match("(?:GAA)+", "")
    assert join_units(empty, 9, 3, "", False, "GAA") == []

def test_join_at_max_gap():
    # "GAAGAATTTGAA"
    runs = indexed_match("(?:GAA)+", "GAAGAATTTGAA")
    assert runs == [(0,6), (9,12)]
    assert join_units(runs, 3, 0, "GAAGAATTTGAA", False, "GAA") == [
    {
        "start": 0,
        "end": 12,
        "span": 12,
        "units": 3, 
        "longest_run": 2, 
        "interruptions": 
            [
            {'start': 6, 
            'end': 9,
            'length': 3, 
            'sequence': 'TTT', 
            'multiple_of_3': True}
            ],
        "reverse_complement": False
    }
    ] 

    assert join_units(runs, 2, 0, "GAAGAATTTGAA", False, "GAA") == [
    {
        "start": 0,
        "end": 6,
        "span": 6,
        "units": 2, 
        "longest_run": 2, 
        "interruptions": [],
        "reverse_complement": False 
    },
    {
        "start": 9,
        "end": 12,
        "span": 3,
        "units": 1, 
        "longest_run": 1, 
        "interruptions": [],
        "reverse_complement": False 
    }
    ]

def test_trailing_bases_no_interruption():
    runs = indexed_match("(?:GAA)+", "GAAGAAGAATTTTT")
    assert join_units(runs, 9, 3, "GAAGAAGAATTTTT", False, "GAA") == [
        {
            "start": 0,
            "end": 9,
            "span": 9,
            "units": 3, 
            "longest_run": 3, 
            "interruptions": [],
            "reverse_complement": False 
        }
        ]

def test_frameshift_flagged():
    runs = indexed_match("(?:GAA)+", "GAAGAATTTGAAGAA")
    assert runs == [(0,6), (9,15)]
    assert join_units(runs, 3, 0, "GAAGAATTTGAAGAA", False, "GAA") == [
        {
        "start": 0,
        "end": 15,
        "span": 15,
        "units": 4, 
        "longest_run": 2, 
        "interruptions": 
            [
            {'start': 6, 
            'end': 9,
            'length': 3, 
            'sequence': 'TTT', 
            'multiple_of_3': True}
            ],
        "reverse_complement": False
        }
    ]

def test_own_reverse_complement():
    sequence = "ATATATTATATAT"
    motif = "AT"
    runs = indexed_match("(?:AT)+", sequence)
    assert runs == [(0,6), (7,13)]
    assert join_units(runs, 3, 0, sequence, False, motif) == [{
        "start": 0,
        "end": 13,
        "span": 13,
        "units": 6, 
        "longest_run": 3, 
        "interruptions": 
            [
            {'start': 6, 
            'end': 7,
            'length': 1, 
            'sequence': 'T', 
            'multiple_of_3': False}
            ],
        "reverse_complement": False
    }]
    revcomp_units_per_tract, revcomp_interruption_length_data = scan_print_revcomp_motif("ATATATTATATAT", "AT", 9, 3, True)
    assert revcomp_units_per_tract == []
    assert revcomp_interruption_length_data == []



def test_check_lower_case_sequence(tmp_path):
    fasta = tmp_path/"lower.fasta"
    fasta.write_text(">lower\ngaagaagaa\n")
    sequence = import_fasta(fasta)
    assert sequence == "GAAGAAGAA"







                        
        


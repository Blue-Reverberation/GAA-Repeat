# Motif Repeat Scanner

## Definitions
**Motif** - A short sequence which will be scanned for run-on repeats
**Unit** – A single unit is the *motif* being scanned for (E.g, GAA)  
**Pure run** – A pure run refers to one or more *units* joined without interruption  
**max_gap** – The largest number of nucleotides which is accepted for *tract* joining.   
**Tract** – A tract is made up of one or more pure runs joined, each of which is at most *max_gap* nucleotides from the neighbouring pure run.   
**min_units** – A tract is discarded if the number of units in a given tract is less than *min_units*  


## Usage
The only requisite for usage is a standard FASTA file containing genetic code.

In a terminal run:

>python main.py [-h] [--motif] [--max-gap [--min-unit] [--plot]  path/to/file.fasta 

  -h, --help           Shows help message and exit  
  --motif MOTIF        Repeat motif of interest (default = GAA)  
  --max-gap MAX_GAP    Largest gap between motifs which is still accepted and joined   
  --min-unit MIN_UNIT   Minimum number of units in a tract to be reported  
  --plot    Plots histogram of pure run length & interruption length

The scanner will **automatically scan for the reverse complement of the motif.**

## Output
The output will include:

1. Every tract found greater than *min_units* units long  

Each tract outputs the following information:
>start:   (Start coordinates of tract)  
end:    (End coordinates of tract)  
span:   (Length of tract in bases)  
units:   (Number of motifs found)  
longest pure run:   (Highest number of motifs found with 0 interruptions)  
reverse complement?:  (Whether the tract is made up of the original motif or of its reverse complement)
>  
>**Joined interruptions present: [...]**  


2. Every joined interruption within a tract

Each joined interruption outputs the following information:
>"Interruption start":    (Start coordinates of interruption)     
"Interruption end"   (End coordinates of iterruption)  
"Interruption span":   (Length of interruption in bases)  
"sequence":     (The interruption itself)  
"multiple of 3?":    (If the interruption length is divisible by 3)  
  
3. Histogram summarising data

Two separate histograms will be output, summarising the distribution of:
- Pure-run length
- Interruption length

PNG files of the histograms will be saved to the project folder.


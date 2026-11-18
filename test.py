from Bio import SeqIO

for record in SeqIO.parse("home/nerd/Documents/College/Python/GAA-Repeat/test-data/pure-run.fasta", "fasta"):
    print(record.id)
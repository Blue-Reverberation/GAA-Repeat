# Motif Repeat Scanner

## Definitions
**Unit** – A single unit is the motif being scanned for (E.g, GAA)
**Pure run** – A pure run refers to one or more units joined without interruption
**max_gap** – The largest number of nucleotides which is accepted for *tract* joining
**Tract** – A tract is made up of one or more pure runs joined, each of which is at most *max_gap* nucleotides from the neighbouring pure run
**min_units** – A tract is discarded if the number of units in a given tract is less than *min_units*


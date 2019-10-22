# motif-search example

Comparison of Branch-and-Bound and Greedy motif finding. Branch and Bound will choke on test1.fasta but Greedy will work.

Command overview:
MotifSearch sequence_filename.txt w op
Description: Program for finding motifs in sequences
using branch and bound method and greedy search.

sequence_file: file containing sequences
w: length of desired motif
op: selected method, bb for branch and bound, greedy
for greedy search

Building:
A makefile is included.

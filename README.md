# Ramachandran Plots

Repository for Bootcamp Summer Programming Project

This script generates Ramachandran plots from a PBD structure file. Structure
files can be provided as paths to local .pdb files, or PDBIDs searchable on the RCSB
Protein Data Bank (http://www.rcsb.org). 

If a directory is provided as the first argument to the script, then the directory is 
searched for any .pdb files that can be used. 

If multiple structure files are provided, in either .pdb or PDBID form, the generated 
Ramachandran plot will contain all residues present in all structures. 

If the "-i" flag is provided, separate Ramachandran plots for glycine, pre-proline
and all other residues will be generated. Plots for 'glycine', 'pre-proline', or 'others'
can be indicated to generate that specific plot. 

This script requires matplotlib (>2.0) for plotting functionality, and
standard library modules urllib, math, os, argparse and f-strings (> Python 3.5) for functionality. 
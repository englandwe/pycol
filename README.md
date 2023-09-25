# pycol
Coloring of residues in PyMol using discrete properties

Custom coloring of residues by discrete identifiers

Takes a list of custom discrete identifiers and colors, generates numerical identifiers and writes as b-factors; then makes pretty colors.

Expects whitespace-separated hex color data in the form 'identifier #XXXXXX' or 'identifier 0xXXXXXX', one entry per line.

Acceptable examples: 
    
data1 #FF00FF
    
data2 0x00FF00
    
Comments (lines starting with #) and blank lines are ignored.

Data file can be in one of three formats (all whitespace-separated, one entry per line):

four columns: chain residueNumber residueName data

A 105 TYR data1

three columns: chain residueNumber data

A 104 data2

two columns acceptable if there's only one chain: residueNumber data

103 data3

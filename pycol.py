'''
Custom coloring of residues by discrete identifiers

Takes a list of custom discrete identifiers and colors, generates numerical identifiers and writes as b-factors
Then makes pretty colors

Expects whitespace-separated hex color data in the form 'identifier #XXXXXX' or 'identifier 0xXXXXXX', one entry per line:
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
'''

import sys
import re
from pymol import stored

comblank=re.compile('^\s*$|^\s*#')

#this function inspired by and largely stolen from data2bfactor.py (http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/)
def importResData(datafile):
    bdat = {}
    chain = ''
    count = 0
    dts = set()
    with open(datafile, 'rU') as infile:
        for line in infile:
            if not comblank.match(line):
                words = line.split()
                count += 1
                if len(words) >= 4:
                    chain = words[0]
                    resi = words[1]
                    resn = words[2]
                    if chain == '-':
                        chain = ''
                    data = words[3]
                #assume we skipped the residue name, not the chain
                elif len(words) == 3:
                    chain = words[0]
                    resi = words[1]
                    resn = ''
                    data = words[2]
                elif len(words) == 2:
                    resi = words[0]
                    data = words[1]
                    resn = ''
                else:
                    sys.stderr.write("\nError in reading data files -- check number of columns.\n")
                    sys.stderr.write("Number of columns: %d  on line number %d\n" % (len(words),count))
                    sys.exit(1)

                bdat.setdefault(chain, {})[resi] = (data, resn)
                dts.add(data)

    return bdat,dts

def importColors(colorfile,dts):
    cols={}
    speclist=['gray']
    ctr=1
    with open(colorfile) as infile:
        for line in infile:
            if not comblank.match(line):
                tmpline=line.strip().split()
                if tmpline[0] in dts:
                    proper_hex=re.sub('^#','0x',tmpline[1])
                    cols[tmpline[0]]=float(ctr)
                    ctr+=1
                    speclist.append(proper_hex)
    return cols,speclist

#also derived from data2bfactor.py (http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/)
def pycol(mol='', datafile='', colorfile='', property='b', quiet=0):
    b_dict,dts = importResData(datafile)
    cols,spectrum_cols=importColors(colorfile,dts)
    quiet = int(quiet) == 1

    def b_lookup(chain, resi, name, b):
        try:
            if chain in b_dict:
                b_id = b_dict[chain][resi][0]
            else:
                b_id = b_dict[''][resi][0]
            b=cols[b_id]
            if not quiet: print('///%s/%s/%s new: %f' % (chain, resi, name, b))
        except KeyError:
            if not quiet: print('///%s/%s/%s keeping: %f' % (chain, resi, name, b))
        return b
    stored.b = b_lookup
    cmd.alter(mol, 'b=0.0')
    cmd.alter(mol, '%s=stored.b(chain, resi, name, %s)' % (property, property))
    cmd.spectrum('b',' '.join(spectrum_cols))
    cmd.rebuild()

cmd.extend('pycol',pycol)


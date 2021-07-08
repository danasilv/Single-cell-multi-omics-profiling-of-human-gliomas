#!/usr/bin/env python

from __future__ import (division,
                        print_function,
                        absolute_import,
                        unicode_literals)

import sys
import numpy as np


def main():
    if "-h" in sys.argv or "--help" in sys.argv:
        print("Prints out the file name, number of discordant reads, number of concordant reads for each file passed in", file=sys.stderr)
        print("compute_pdr.py\t[CpG_context1.txt CpG_context2.txt ...]", file=sys.stderr)
        return

    filenames = sys.argv[1:] 

    #we need to pass in OB, OT, Run1, Run2
    #read_id, n_meth, n_umeth 
    include_pairs = False
    header = ["FileName", "DiscordantReads", "ConcordantReads"]
    print(" ".join(header))
    for filename in filenames:
        read_ids = []
        current_read = "" 
        current_meth = 0
        current_umeth = 0
        with open(filename) as fd:
            for i, line in enumerate( fd ):
                if i % 10000 == 0:
                    print("Proccesed lines", i, "\r", file=sys.stderr)
                if i == 0: continue
                line_arr = line.strip().split()
                
                if line_arr[0] != current_read:
                    read_ids.append([current_read, current_meth, current_umeth])
                    current_read = line_arr[0]
                    current_meth, current_umeth = 0,0
                    start_pos = int(line_arr[3])
                elif not include_pairs:
                    chrm, pos = line_arr[2], int(line_arr[3])
                    if abs(pos - start_pos) >= 90:
                    continue
                    

                if line_arr[-1] == 'z':
                    current_umeth += 1
                else:
                    current_meth += 1 
    
        print(filename, " ".join(list(map(str, compute_pdr(read_ids)))))
    
def compute_pdr(reads):
    #read_id, umeth, meth
    #pdr is defined as n-discordant reads / total reads
    #ignore reads where meth + umeth < 4 or > 6

    """
    @returns n_discordant, n_concordant. Add them to get the total number.
    """
    
    reads = np.asarray(reads) 
    read_ids = reads[:, 0]
    meth_umeth = np.asarray(reads[:, 1:], dtype=int)

    meth_umeth_qualified = meth_umeth[np.logical_and(np.sum(meth_umeth, axis = 1) >= 4, np.sum(meth_umeth, axis = 1) <= 6)] 
    
    n_concordant = 0
    for row in meth_umeth_qualified:
        if row[0] < 1 or row[1] < 1:
            n_concordant += 1
    
    return len(meth_umeth_qualified) - n_concordant, n_concordant
    
        

if __name__ == "__main__":
    main()

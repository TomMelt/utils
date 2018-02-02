#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys
import getopt
import re

HARTREE_TO_EV=27.2113839
pd.options.display.float_format = "{:>10.3f}".format

def main(argv):

    for arg in argv:
        filename=arg
        getTargetStatesMolPro(filename)

def getTargetStatesMolPro(inFilename):
    with open(inFilename) as infile:
        num,sym,energy=[],[],[]
        for line in infile:
            if "!MCSCF STATE" in line and "Energy" in line:
                numSym=line.split()[2]
                num.append(int(numSym.split(".")[0]))
                sym.append(int(numSym.split(".")[1])-1)
                energy.append(float(line.split()[-1]))

        symLookup=pd.DataFrame({"sym": [0,1,2,3], "key": ["a1","b1","b2","a2"]}, columns=["sym","key"])
        df=pd.DataFrame({"num": num, "sym": sym, "energy": energy}, columns=["num","sym","energy"])
        df=pd.merge(df, symLookup, how='left', on='sym')

        df["energy"]=df["energy"].apply(lambda x: (x)*HARTREE_TO_EV)

        print(df[["num","key","energy","sym"]].sort_values(by=["energy"]).head(30).to_string(index=False))
    return


if __name__ == "__main__": main(sys.argv[1:])

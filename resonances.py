#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys

GOODNESS_CUTOFF=1.0E-1
WIDTH_CUTOFF=2.0

def getResonanceInformation(infile):
    currentSym=""
    symmetry, position, width, goodness=[],[],[],[]
    for line in infile:
        if "Scattering state symmetry is" in line:
            currentSym=list(filter(None,line.replace("\n","").split(" ")))[-1]
        if "Positions" in line:
            symmetry.append(currentSym)
            position.append(float(list(filter(None,line.replace("\n","").replace("D","E").split(" ")))[-1]))
        if "Widths" in line:
            width.append(float(list(filter(None,line.replace("\n","").replace("D","E").split(" ")))[-1]))
        if "Goodness" in line:
            goodness.append(float(list(filter(None,line.replace("\n","").replace("D","E").split(" ")))[-1]))
    return symmetry, position, width, goodness

#==============================
#==== read info from files ====
#==============================

pd.options.display.float_format = '{:>10.3f}'.format

if len(sys.argv) > 1:
    try:
        print("Exclusion Criteria:")
        print("GOODNESS_CUTOFF=",GOODNESS_CUTOFF)
        print("WIDTH_CUTOFF=",WIDTH_CUTOFF)
        for filename in sys.argv[1:]:
            print("filename::",filename)
            with open(filename) as infile:
                symmetry, position, width, goodness=getResonanceInformation(infile)
                df=pd.DataFrame({"symmetry": symmetry, "position": position, "width": width, "goodness": goodness}, columns=["symmetry", "position", "width", "goodness"])
                df=df.ix[df["goodness"]<GOODNESS_CUTOFF]
                df=df.ix[df["width"]<WIDTH_CUTOFF]
                print(df.to_string(index=False))
    except (IOError, OSError) as e:
        print("ERROR:: file \""+e.filename+"\" does not exist")
        quit()
else:
    print("ERROR:: Please provide atleast one path to reson_message file")
    print("        e.g. ./resonances.py /home/USER/Documents/reson_message")

quit()


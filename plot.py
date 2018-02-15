#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import getopt
import re

HARTREE_TO_EV=27.2113839
pd.options.display.float_format = "{:>10.3f}".format
USAGE='''Usage :: ./scat_input.py [-h] [-rs] -p <path> -c <frozen> <occupied> <virtual>
-h help
-r create rmat?.data files
-s create ci?np1.data files
-p path to folder which contains citargmod.prop and rmat?.data
-c frozen, occupied and virtual MO's in order of symmetry (A1,B1,B2,A2)
   e.g., -c 4,0,0,0 6,2,2,0 4,2,2,0
   would correspond to:
       4 A1, 0 B1, 0 B2, 0 A2 frozen
       6 A1, 2 B1, 2 B2, 0 A2 occupied
       4 A1, 2 B1, 2 B2, 0 A2 virtual'''
OPTIONS=["r","s"]

def main(argv):

#    path,frozen,occupied,virtual="","","",""
#    processed=[]
#    for i, arg in enumerate(argv):
#        if arg == "-h":
#            print(USAGE)
#            quit()
#        if arg == "-p":
#            path=argv[i+1]
#        if arg == "-c":
#            frozen,occupied,virtual=argv[i+1:i+4]
#
#    argv=[i for i in argv if i not in ["-p","-c",path,frozen,occupied,virtual]]
#    argv="".join(argv).replace("-","")
#
#    for arg in argv:
#        if arg not in OPTIONS:
#            print("Error :: Unknown option \'"+str(arg)+"\'")
#            print(USAGE)
#            quit()

    print("arguments are ",argv)

    if argv[0] is not None:
        try:
            df=pd.DataFrame()
            rvalues=[]
            print("Reading molpro output files...")
            for i,filename in enumerate(argv):
                r,sym,energy=getMCSCF(filename)
                rvalues.append(r)
                if filename is argv[0]:
                    df=pd.DataFrame(columns=sym)
                    df.loc[i]=energy
                else:
                    df.loc[i]=energy

            df=df.rename({i:r for i,r in enumerate(rvalues)}, axis="index")
            df=df.rename({"1":"A1","2":"B1","3":"B2","4":"A2"}, axis="columns")
            df=df.apply(pd.to_numeric)
            groundEnergy=df.values.min()
            print("Min. Ground State Energy (eV) ::",groundEnergy*HARTREE_TO_EV)
            df=df.applymap(lambda x: (x-groundEnergy)*HARTREE_TO_EV)

            print(df)

            plot(df)

        except (IOError, OSError) as e:
            print("ERROR:: one path is invalid")
            quit()
    else:
        print("ERROR:: Please provide at least one path")
        print("        e.g. ./plot.py molpro.out")

def getMCSCF(filename):

    sym,energy=[],[]
    rindex=0
    r=""
    with open(filename) as infile:
        for i,line in enumerate(infile):
            if "Bond lengths in Bohr" in line:
                rindex=i+3
            if all(word in line for word in ["!MCSCF STATE","Energy"]):
                sym.append(line.split()[2].split(".")[-1])
                energy.append(line.split()[-1])
            if rindex > 0 and i==rindex:
                r=line.replace("(","").replace(")","").split()[0]
    r=float(r)
    return r,sym,energy

def plot(df):
    df.drop(["B2"],axis=1,inplace=True)
    print(df)
    df.plot()
    plt.legend(loc='best')
    plt.show()

if __name__ == "__main__": main(sys.argv[1:])

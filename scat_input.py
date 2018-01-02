#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys

HARTREE_TO_EV=27.2113839
pd.options.display.float_format = '{:>10.3f}'.format


def main():
    if len(sys.argv) > 1:
        try:
            for filename in sys.argv[1:]:
                with open(filename) as infile:
                    df=getTargetInformation(infile)
                    #print(df.to_string(index=False))

                    # number of target states per spin-space sym
                    numTarg=(df["mult"]+df["sym"]).value_counts(sort=False)

                    df=df.sort_values(by=["mult","sym"])
                    df["energy"].replace(r"D",r"E",regex=True,inplace=True)
                    df=df.apply(pd.to_numeric)
                    groundEnergy=min(df["energy"])
                    df["energy"]=df["energy"].apply(lambda x: (x-groundEnergy)*HARTREE_TO_EV)

                    # get idtarg (used in rmat files)
                    idtarg=fortstring(df["num"])
                    nvo=fortstring(np.zeros(len(df["num"]), dtype=np.int))

                    print(df.sort_values(by=["energy"]).to_string(index=False))

                    print("add following lines to rmat?.data")
                    print("\tidtarg="+idtarg)
                    print("\tnvo="+nvo)

        except (IOError, OSError) as e:
            print("ERROR:: file \""+e.filename+"\" does not exist")
            quit()
    else:
        print("ERROR:: Please provide atleast one path to reson_message file")
        print("        e.g. ./resonances.py /home/USER/Documents/reson_message")

def getTargetInformation(infile):

    #=================================================
    #==== Get order of target states from denprop ====
    #=================================================

    num,mult,sym,energy=[],[],[],[]

    for line in infile:
        cols=line.split()
        if cols[0]=="5":
            num.append(cols[1])
            sym.append(cols[4])
            mult.append(cols[5])
            energy.append(cols[8])

    return pd.DataFrame({"num": num, "sym": sym, "mult": mult, "energy": energy}, columns=["num","sym","mult","energy"])

def convertToeV(string):
    return float(string)*HARTREE_TO_EV

def fortstring(array):
    string=""
    for elem in array:
        string+=str(elem)+","
    return string





if __name__ == "__main__": main()

##================================
##==== Produce CI Input Files ====
##================================
#
## number of continuum orbitals per spatial sym
#numCont=np.array([53,34,34,22])
#
## form a space-spin symmetry table (8x8 matrix)
#symTable=np.array([0,1,2,3,1,0,3,2,2,3,0,1,3,2,1,0])
#symTable=symTable.reshape(4,4)
#symTable=np.concatenate((symTable, symTable), axis=0)
#symTable=np.concatenate((symTable, symTable), axis=1)
#
#print("\nSymmetry Table c2v:\n", symTable)
#
#
#for scatSym in range(0,8):
#
#    f = open('ci'+str(scatSym+1)+'np1.data', 'w')
#
#    if scatSym < 4:
#        numSym=4
#    else:
#        numSym=8
#
#    notgt=""
#    mcont=""
#    numtgt=""
#    for targSym in range(0,numSym):
#        notgt+=str(numCont[symTable[scatSym,targSym]])+","
#        mcont+=str(symTable[scatSym,targSym])+","
#        numtgt+=str(numTarg[targSym])+","
#
#    fileText=["&INPUT", "\tiexpc=1,", "\ticitg=1,", "\tntgsym="+str(numSym)+",", "\tidiag=2,", "\tnftg=0,", "\tnotgt="+notgt, "\tmcont="+mcont, "\tnumtgt="+numtgt, "&END", "&CINORN", "name = '(N+1) electron CI',", "nciset="+str(scatSym+1), "/"]
#
#    for lines in fileText:
#        f.write(lines+"\n")
#
#
#print("\ndone")
#

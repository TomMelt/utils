#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys
import getopt
import re

HARTREE_TO_EV=27.2113839
DENPROP="citargmod.prop"
NP1SWEDMOS="../output/swedmos.np1.out"
NOELEC=15
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

    path,frozen,occupied,virtual="","","",""
    processed=[]
    for i, arg in enumerate(argv):
        if arg == "-h":
            print(USAGE)
            quit()
        if arg == "-p":
            path=argv[i+1]
        if arg == "-c":
            frozen,occupied,virtual=argv[i+1:i+4]

    argv=[i for i in argv if i not in ["-p","-c",path,frozen,occupied,virtual]]
    argv="".join(argv).replace("-","")

    for arg in argv:
        if arg not in OPTIONS:
            print("Error :: Unknown option \'"+str(arg)+"\'")
            print(USAGE)
            quit()

    if path is not None:
        try:
            print("Extracting info from Denprop...")
            filename=path+"/"+DENPROP
            numTarg,ntarg,idtarg,nvo=getTargetInformation(filename)

            print("Extracting Target Energies from molpro.out...")
            getTargetStatesMolPro(path+"/molpro.out")

            if occupied is not "":
                createCongenFiles(path,frozen,occupied,virtual)

            if "s" in argv:
                print("Extracting No. of Scatt. Orbs. from swedmos.np1...")
                filename=path+"/"+NP1SWEDMOS
                numContOrbitals=getNumContOrb(filename)
                print("Creating ci?np1.data files...")
                createScatCInp1files(numTarg,numContOrbitals)

            if "r" in argv:
                print("Creating rmat?.data files...")
                for rmatFile in FilenameGenerator("rmat",".data",1,8):
                    inFilename=path+"/"+rmatFile
                    outFilename=rmatFile
                    createRmatFiles(inFilename,outFilename,ntarg,idtarg,nvo)

        except (IOError, OSError) as e:
            print("ERROR:: path \""+path+"\" is not a valid path")
            quit()
    else:
        print("ERROR:: Please provide at least one path")
        print("        e.g. ./resonances.py /home/USER/Documents/input")

def getTargetInformation(filename):

    with open(filename) as infile:

        num,mult,sym,energy=[],[],[],[]

        for line in infile:
            cols=line.split()
            if cols[0]=="5":
                num.append(cols[1])
                sym.append(cols[4])
                mult.append(cols[5])
                energy.append(cols[8])

        df=pd.DataFrame({"num": num, "sym": sym, "mult": mult, "energy": energy}, columns=["num","sym","mult","energy"])

        df=df.sort_values(by=["mult","sym"])
        # number of target states per spin-space sym
        numTarg=df.groupby(["mult","sym"]).count().num.tolist()

        df["energy"].replace(r"D",r"E",regex=True,inplace=True)
        df=df.apply(pd.to_numeric)
        groundEnergy=min(df["energy"])
        print("Ground State Energy (eV) ::",groundEnergy*HARTREE_TO_EV)
        df["energy"]=df["energy"].apply(lambda x: (x-groundEnergy)*HARTREE_TO_EV)

        # get idtarg (used in rmat files)
        ntarg=len(df["num"])
        idtarg=fortstring(df["num"])
        nvo=fortstring(np.zeros(ntarg,dtype=int))

        print(df.sort_values(by=["energy"]).to_string(index=False))

    return numTarg,ntarg,idtarg,nvo

def fortstring(array):
    string=""
    for elem in array:
        string+=str(elem)+","
    return string

def FilenameGenerator(part1,part2,start,finish):
    for i in range(start,finish+1):
        yield part1+str(i)+part2

def createRmatFiles(inFilename,outFilename,ntarg,idtarg,nvo):
    with open(inFilename) as infile:
        outfile=open(outFilename,"w")
        for line in infile:
            if "ntarg=" in line:
                outfile.write("\tntarg="+str(ntarg)+",\n")
            elif "idtarg=" in line:
                outfile.write("\tidtarg="+idtarg+"\n")
            elif "nvo=" in line:
                outfile.write("\tnvo="+nvo+"\n")
            else:
                outfile.write(line)
    return


def createScatCInp1files(numTarg,numContOrbitals):
    symTable=symmetryTable()
    #print("\nSymmetry Table c2v:\n", symTable)

    for scatSym in range(0,8):

        outfile = open("ci"+str(scatSym+1)+"np1.data", "w")

        if scatSym < 4:
            numSym=4
        else:
            numSym=8

        notgt=""
        mcont=""
        numtgt=""
        for targSym in range(0,numSym):
            notgt+=str(numContOrbitals[symTable[scatSym,targSym]])+","
            mcont+=str(symTable[scatSym,targSym])+","
            numtgt+=str(numTarg[targSym])+","

        fileText=["&INPUT", "\tiexpc=1,", "\ticitg=1,", "\tntgsym="+str(numSym)+",", "\tidiag=2,", "\tnftg=0,", "\tnotgt="+notgt, "\tmcont="+mcont, "\tnumtgt="+numtgt, "&END", "&CINORN", "name = '(N+1) electron CI',", "nciset="+str(scatSym+1), "/"]

        for lines in fileText:
            outfile.write(lines+"\n")

    return

def getNumContOrb(filename):
    numContOrbitals=[]
    with open(filename) as infile:
        for line in infile:
            if "No. of scattering orbitals" in line:
                numContOrbitals.append(line.split()[-1])
    return numContOrbitals

def createCongenFiles(path,frozen,occupied,virtual):

    getMOOrder(path+"/molpro.out")

    frozen=np.array(frozen.split(","),dtype=int)
    occupied=np.array(occupied.split(","),dtype=int)
    virtual=np.array(virtual.split(","),dtype=int)
    continuum=np.array([int(i) for i in getNumContOrb(path+"/"+NP1SWEDMOS)])

    numFrozen=frozen.sum()*2
    numCAS=NOELEC-numFrozen

    #==================================
    #           congen?.data
    #==================================
    nelecp=np.array([numFrozen,numCAS])
    pqnStr=""
    mshl,nshlp=[],[]
    for MOType in [frozen,occupied]:
        nshlpCounter=0
        for sym,i in enumerate(MOType):
            if i > 0:
                nshlpCounter=nshlpCounter+1
                mshl.append(sym)
                if MOType is frozen:
                    pqnStr+="0,1,"+str(i)+", "
                if MOType is occupied:
                    pqnStr+="0,"+str(frozen[sym]+1)+","+str(i)+", "
        nshlp.append(nshlpCounter)

    for congenFile in FilenameGenerator("congen",".data",1,8):
        inFilename=path+"/"+congenFile
        outFilename=congenFile
        with open(inFilename) as infile:
            outfile=open(outFilename,"w")
            for line in infile:
                if "nob=" in line:
                    outfile.write("\tnob="+",".join(map(str,occupied))+",\n")
                elif "nob0=" in line:
                    outfile.write("\tnob0="+",".join(map(str,occupied))+",\n")
                elif "nelecp=" in line:
                    outfile.write("\tnelecp="+",".join(map(str,nelecp))+",\n")
                elif "nshlp=" in line:
                    outfile.write("\tnshlp="+",".join(map(str,nshlp))+",\n")
                elif "mshl=" in line:
                    outfile.write("\tmshl="+",".join(map(str,mshl))+",\n")
                elif "pqn=" in line:
                    outfile.write("\tpqn="+pqnStr+"\n")
                else:
                    outfile.write(line)

    #==================================
    #           congen?np1.data
    #==================================
    symTable=symmetryTable()
    print("\nSymmetry Table c2v:\n", symTable)
    for i,congenFile in enumerate(FilenameGenerator("congen","np1.data",1,8)):
        inFilename=path+"/"+congenFile
        outFilename=congenFile
        numGroups=4 if i < 4 else 8
        groupName=""
        groupCount=0
        with open(inFilename) as infile:
            outfile=open(outFilename,"w")
            for line in infile:
                if "nob=" in line:
                    outfile.write("\tnob="+",".join(map(str,occupied+virtual+continuum))+",\n")
                elif "nob0=" in line:
                    outfile.write("\tnob0="+",".join(map(str,occupied))+",\n")
                elif "gname=" in line:
                    groupName=line
                    outfile.write(line)
                elif "continuum" in groupName:
                    if "nelecp=" in line:
                        outfile.write("\tnelecp="+",".join(map(str,nelecp))+",1,\n")
                    elif "nshlp=" in line:
                        outfile.write("\tnshlp="+",".join(map(str,nshlp))+",1,\n")
                    elif "mshl=" in line:
                        mshlCont=symTable[i,groupCount]
                        outfile.write("\tmshl="+", ".join(map(str,mshl))+", "+str(mshlCont)+",\n")
                        groupCount=groupCount+1
                    elif "pqn=" in line:
                        mshlCont=symTable[i,groupCount]
                        orbital=occupied[mshlCont]+virtual[mshlCont]+1
                        pqnStrCont="0,"+str(orbital)+","+str(orbital+1)+", "
                        outfile.write("\tpqn="+pqnStr+pqnStrCont+"\n")
                    else:
                        outfile.write(line)
                elif "virtual" in groupName:
                    if "nelecp=" in line:
                        outfile.write("\tnelecp="+",".join(map(str,nelecp))+",1,\n")
                    elif "nshlp=" in line:
                        nshlpVirt=sum(1 for j in virtual if j > 0)
                        outfile.write("\tnshlp="+",".join(map(str,nshlp))+","+str(nshlpVirt)+",\n")
                    elif "mshl=" in line:
                        mshlVirt=[k for k,j in enumerate(virtual) if j > 0]
                        outfile.write("\tmshl="+", ".join(map(str,mshl))+", "+", ".join(map(str,mshlVirt))+",\n")
                        groupCount=groupCount+1
                    elif "pqn=" in line:
                        pqnStrVirt=""
                        for k,j in enumerate(virtual):
                            if j > 0:
                                 pqnStrVirt=pqnStrVirt+"0,"+str(occupied[k]+1)+","+str(occupied[k]+virtual[k])+", "
                        outfile.write("\tpqn="+pqnStr+pqnStrVirt+"\n")
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)


    return

def getMOOrder(inFilename):
    with open(inFilename) as infile:
        j=-1
        num,sym,occ,energy=[],[],[],[]
        for i,line in enumerate(infile):
            if "Orb     Occ        Energy       Coefficients" in line:
                j=i
            if j>0 and i>j:
                tempStr=line[0:30]
                if bool(re.search(r'\d', tempStr)):
                    numSym=tempStr.split()[0]
                    num.append(int(numSym.split(".")[0]))
                    sym.append(int(numSym.split(".")[1])-1)
                    occ.append(-1.0*float(tempStr.split()[1]))
                    energy.append(float(tempStr.split()[-1])*HARTREE_TO_EV)
            if "Natural orbital dump" in line:
                j=-1

        symLookup=pd.DataFrame({"sym": [0,1,2,3], "key": ["a1","b1","b2","a2"]}, columns=["sym","key"])
        df=pd.DataFrame({"num": num, "sym": sym, "occ": occ, "energy": energy}, columns=["num","sym","occ","energy"])
        df=pd.merge(df, symLookup, how='left', on='sym')
        print("First 30 MO's")
        print(df[["num","key","energy","occ","sym"]].sort_values(by=["occ","energy"]).head(30).to_string(index=False))
    return

def getTargetStatesMolPro(inFilename):
    with open(inFilename) as infile:
        num,sym,energy=[],[],[]
        for line in infile:
            if "!MCSCF STATE" in line and "Energy" in line:
                print(line.split())
                numSym=line.split()[2]
                num.append(int(numSym.split(".")[0]))
                sym.append(int(numSym.split(".")[1])-1)
                energy.append(float(line.split()[-1]))

        symLookup=pd.DataFrame({"sym": [0,1,2,3], "key": ["a1","b1","b2","a2"]}, columns=["sym","key"])
        df=pd.DataFrame({"num": num, "sym": sym, "energy": energy}, columns=["num","sym","energy"])
        df=pd.merge(df, symLookup, how='left', on='sym')

        groundEnergy=min(df["energy"])
        print("Ground State Energy (eV) ::",groundEnergy*HARTREE_TO_EV)
        df["energy"]=df["energy"].apply(lambda x: (x-groundEnergy)*HARTREE_TO_EV)

        print(df[["num","key","energy","sym"]].sort_values(by=["energy"]).head(30).to_string(index=False))
    return

def symmetryTable():
    # form a space-spin symmetry table (8x8 matrix)
    symTable=np.array([0,1,2,3,1,0,3,2,2,3,0,1,3,2,1,0])
    symTable=symTable.reshape(4,4)
    symTable=np.concatenate((symTable, symTable), axis=0)
    symTable=np.concatenate((symTable, symTable), axis=1)
    return symTable


if __name__ == "__main__": main(sys.argv[1:])

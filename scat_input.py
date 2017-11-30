#!/usr/bin/python3
import numpy as np

def fortstring(array):
    string=""
    for i in range(len(array)):
        string+=str(array[i])+","
    return string

#=================================================
#==== Get order of target states from denprop ====
#=================================================

filename = 'citargmod.prop'

# get no of target states (used in rmat files)
ntarg=np.genfromtxt(filename, dtype=int, usecols=(4), max_rows=1)

# store target states to np array
data=np.genfromtxt(filename, dtype=None, skip_header=3, usecols=(1,4,5,8), names=('num','sym','mult','energy'), max_rows=ntarg)

# sort states by the multiplicity then the symmetry
data=np.sort(data, order=('mult','sym'))

# get idtarg (used in rmat files)
idtarg=fortstring(data['num'])
nvo=fortstring(np.zeros((ntarg), dtype=np.int))

# number of target states per spin-space sym
numTarg=data['mult']*10+data['sym']
unique, numTarg= np.unique(numTarg, return_counts=True)

print("idtarg="+idtarg)
print("nvo="+nvo)

#================================
#==== Produce CI Input Files ====
#================================

# number of continuum orbitals per spatial sym
numCont=np.array([53,34,34,22])

# form a space-spin symmetry table (8x8 matrix)
symTable=np.array([0,1,2,3,1,0,3,2,2,3,0,1,3,2,1,0])
symTable=symTable.reshape(4,4)
symTable=np.concatenate((symTable, symTable), axis=0)
symTable=np.concatenate((symTable, symTable), axis=1)

print("\nSymmetry Table c2v:\n", symTable)


for scatSym in range(0,8):

    f = open('ci'+str(scatSym+1)+'np1.data', 'w')

    if scatSym < 4:
        numSym=4
    else:
        numSym=8

    notgt=""
    mcont=""
    numtgt=""
    for targSym in range(0,numSym):
        notgt+=str(numCont[symTable[scatSym,targSym]])+","
        mcont+=str(symTable[scatSym,targSym])+","
        numtgt+=str(numTarg[targSym])+","

    fileText=["&INPUT", "\tiexpc=1,", "\ticitg=1,", "\tntgsym="+str(numSym)+",", "\tidiag=2,", "\tnftg=0,", "\tnotgt="+notgt, "\tmcont="+mcont, "\tnumtgt="+numtgt, "&END", "&CINORN", "name = '(N+1) electron CI',", "nciset="+str(scatSym+1), "/"]

    for lines in fileText:
        f.write(lines+"\n")


print("\ndone")


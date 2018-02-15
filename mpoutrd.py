#!/usr/bin/python3
import numpy as np
import pandas as pd

def fortstring(array):
    string=""
    for i in range(len(array)):
        string+=str(array[i])+", "
    return string

def stringFixedLeft(outfile,line):
    outfile.write("!  ")
    for num, part in enumerate(line):
        outfile.write("%-10s"%part)
        if (num+1)%10 == 0:
            outfile.write("\n!  ")
    outfile.write("\n")
    return

def stringFixedRight(outfile,line):
    for num, part in enumerate(line):
        outfile.write(" %10s "%part)
        if (num+1)%10 == 0:
            outfile.write("\n")
    outfile.write("\n")
    return

#=============================================
#==== RMat default order of coefficients  ====
#=============================================
rmatCoeff=['1s','2px','2py','2pz','3d0','3d2-','3d1+','3d2+','3d1-','4f1+','4f1-','4f0','4f3+','4f2-','4px','4f3-','4f2+','4py','4pz','5g0','5g2-','5g1+','5g4+','5g3-','5g2+','5g4-','5g3+','5d2-','5d1+','5s','5g1-','5d2+','5d1-','5d0']

rmatLookup = pd.DataFrame({"coeff": rmatCoeff, "rmatOrder": range(0,len(rmatCoeff))})
print("Only the following coefficients are supported")
print(rmatLookup)

#=============================================
#==== import useful info from molpro.out  ====
#=============================================
MOLPRO_COEFF_WIDTH=10
MOLPRO_COEFF_START_COL=32
START_LINE=0
LINES=0

numContractions=np.zeros((4,), dtype=int)
linesPerContraction=np.zeros((4,), dtype=int)
filename = 'molpro.out'
data=[]

with open(filename) as infile:
    for num, line in enumerate(infile, 1):
        # find the number of contractions for each space symmetry (A1,B1,B2,A2)
        if "NUMBER OF CONTRACTIONS" in line:
            parts = line.split(" ")
            for part in parts:
                if "A1" in part:
                    numContractions[0]=int(str.replace(part,"A1",""))
                if "B1" in part:
                    numContractions[1]=int(str.replace(part,"B1",""))
                if "B2" in part:
                    numContractions[2]=int(str.replace(part,"B2",""))
                if "A2" in part:
                    numContractions[3]=int(str.replace(part,"A2",""))

        # find the line number where the coefficients begin
        if "Orb     Occ        Energy       Coefficients" in line:
            START_LINE=num+2
            print(START_LINE)
            for i, n in enumerate(numContractions):
                linesPerContraction[i] = (int((n-1)/MOLPRO_COEFF_WIDTH)+2)
                print(linesPerContraction[i],n,MOLPRO_COEFF_WIDTH,int(n/MOLPRO_COEFF_WIDTH))
                LINES += linesPerContraction[i]*(n+1)
                print(START_LINE+LINES,linesPerContraction)

        # read in all required data
        if num >= START_LINE and num < START_LINE+LINES:
            data.append(line)

print("numContractions=",numContractions)
#print(data)

coeffinfo=[]

A1df=pd.DataFrame()
B1df=pd.DataFrame()
B2df=pd.DataFrame()
A2df=pd.DataFrame()

# read data into pandas framework
for i, df in enumerate([A1df, B1df, B2df, A2df]):

    START_LINE=np.sum((numContractions[0:i]+1)*linesPerContraction[0:i])

    temp=[]

    for line in data[START_LINE:START_LINE+linesPerContraction[i]]:
        temp+= list(filter(None, line[MOLPRO_COEFF_START_COL:].replace("\n","").split(" ")))

    site, coeff=temp[::2], temp[1::2]

    tempdf=pd.DataFrame({'site': site,'coeff': coeff})
    df = pd.concat([df, tempdf], axis=1)

    for eachContraction in range(0,numContractions[i]):

        START_LINE+=linesPerContraction[i]
        coeffValues=[]

        for num, line in enumerate(data[START_LINE:START_LINE+linesPerContraction[i]]):

            if num == 0:
                tempInfo=list(filter(None, line[:MOLPRO_COEFF_START_COL].replace("\n","").split(" ")))
                # tempInfo is changed to match the format of MPOUTRD, this step is not necessary
                tempInfo[0]+=";"
                tempInfo[1]="Energy="
                tempInfo[2]+=" H"
                coeffinfo.append(tempInfo)

            coeffValues+= list(filter(None, line[MOLPRO_COEFF_START_COL:].replace("\n","").split(" ")))

        tempdf=pd.DataFrame({'MO '+str(eachContraction+1): coeffValues})
        df = pd.concat([df, tempdf], axis=1)

    df=pd.merge(df, rmatLookup, how='left', on='coeff')
    df=df.sort_values(by=['site','rmatOrder'])

    if i==0:
        A1df = df
    if i==1:
        B1df = df
    if i==2:
        B2df = df
    if i==3:
        A2df = df

# Write the information to swedmos.inp file
filename="swedmos.inp"

topString="&INPUT\
        \nTITLE=\'\',\
        \nIVCS=0,\
        \nNSYMT=4,\
        \nNBFT= "+fortstring(numContractions)+"\
        \nNOBT= "+fortstring(numContractions)+"\
        \nIND=  1,  1,  1,  "+str(numContractions[0])+",    1,  2,  1,  "+str(numContractions[1])+",    1,  3,  1,  "+str(numContractions[2])+",    1,  4,  1,  "+str(numContractions[3])+",\
        \nVMORB=\n"

outfile=open(filename, 'w')

outfile.write(topString)

for i, df in enumerate([A1df, B1df, B2df, A2df]):
    stringFixedLeft(outfile, df['coeff'])
    outfile.write("\n")
    for eachContraction in range(0,numContractions[i]):
        stringFixedRight(outfile, df['MO '+str(eachContraction+1)])
        outfile.write("\n")
outfile.write("\n&END")


quit()


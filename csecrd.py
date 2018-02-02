#!/usr/bin/python3
import numpy as np
import pandas as pd
import sys

def FilenameGenerator(part1,part2,start,finish):
    for i in range(start,finish+1):
        yield part1+str(i)+part2

def getCrossSection(infile):
    energy,crosssection=[],[]
    isRead=False
    rows,cols=0,0
    for line in infile:
        if "T-Matrices" in line:
            isRead=False
        if isRead:
            x,y=line.replace("D","E").split()[1:3]
            energy.append(x)
            crosssection.append(y)
        if "TOTAL" in line:
            isRead=True
            cols=cols+1
    rows=int(len(crosssection)/cols)
    for i in energy, crosssection:
        i = np.array(i,dtype=float)
        i = i.reshape(-1,cols)
    energy=np.array(energy,dtype=float)
    crosssection=np.array(crosssection,dtype=float)
    energy=energy.reshape(cols,-1)
    crosssection=crosssection.reshape(cols,-1)
    df=pd.DataFrame({'energy':energy[1,:].tolist(), 'crosssection':crosssection.sum(axis=0).tolist()},columns=["energy","crosssection"])
    return df

#==============================
#==== read info from files ====
#==============================

pd.options.display.float_format = '{:>10.5f}'.format

if len(sys.argv) > 1:
    try:
        for i,filename in enumerate(FilenameGenerator(sys.argv[1]+"/fort.3","",0,7)):
            print("filename::",filename)
            with open(filename) as infile:
                df=getCrossSection(infile)
                df.to_csv(sys.argv[1]+"/plot3"+str(i)+".dat",sep=" ",header=False,index=False)
    except (IOError, OSError) as e:
        print("ERROR:: file \""+e.filename+"\" does not exist")
        quit()
else:
    print("ERROR:: Please provide atleast one path to fort.3? file")
    print("        e.g. ./csecrd /home/USER/Documents/input/fort.3?")

quit()


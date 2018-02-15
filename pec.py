#!/usr/bin/python3
import numpy as np
import os
import pandas as pd
import sys

ANGSTROM_TO_AU=1.8897261243

def main(args):

    # read submission script
    filename="sub.job"
    subline=[]
    with open(filename) as infile:
        for i,line in enumerate(infile):
            subline.append(line)

    # read input file
    filename="molpro.inp"
    outline=[]
    changelines=[]
    with open(filename) as infile:
        for i,line in enumerate(infile):
            outline.append(line)
            if "GEOMETRY" in line.upper():
                changelines=[i+1,i+2]

    # mass of N and O atoms
    massN=14.00670
    massO=15.99940
    massTotal=massN+massO

    # internuclear separation range
    rmin=1.05
    rmax=1.40
    dr=0.01
    n=int(round((rmax-rmin)/dr,1))+1

    # create output file for each distance
    for i in range(n):
        r=rmin+dr*i
        dispN= massO/massTotal*r
        dispO=-massN/massTotal*r

        outline[changelines[0]]="N, 0.000000, 0.000000, {0:0.6f}\n".format(dispN)
        outline[changelines[1]]="O, 0.000000, 0.000000, {0:0.6f}\n".format(dispO)

        filename="r{0:0.2f}".format(r).replace(".","")
        subname="rsub{0:0.2f}".format(r).replace(".","")
        outfile=open(filename,"w")
        for line in outline:
            outfile.write(line)
        subfile=open(subname,"w")
        for line in subline:
            if "job_name" in line:
                subfile.write("#$ -N j{0:0.2f}\n".format(r).replace(".",""))
            elif "-wd" in line:
                folders=os.getcwd().split("/")
                for folder in folders[:]:
                    folders.remove(folder)
                    if folder == "Scratch":
                        break
                subfile.write("#$ -wd /home/ucaptme/Scratch/NO/"+"/".join(folders)+"\n")
            elif "molpro.inp" in line:
                subfile.write("/shared/ucl/apps/molpro/2015.1.3/bin/molpro --no-xml-output {}\n".format(filename))
            else:
                subfile.write(line)


if __name__ == "__main__": main(sys.argv[1:])

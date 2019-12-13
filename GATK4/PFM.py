#!/usr/bin/env python3

####### Author
# Pauline Fourgoux p.f.m.fourgoux@qmul.ac.uk ; pauline.fourgoux@gmail.com
# alobley@gmail.com
####### Packages call
import pandas as pd
import sys
import os
####### Funtions
def usage():
    print("""
    computing_RNAseq_AF.py takes the ID of the sample to open a vcf file and calculates the real AF for RNA-seq mutation calling.

    Example of command line:
    =======================

    python computing_RNAseq_AF.py --in-vcf myin.vcf --out-vcf myout.vcf

    obligatory:
    ==========

    --in-vcf  -> VCF
    --out-vcf -> VCF
    """)

####### Main
# Get arguments
infile=""
outfile=""
#---------------------------------------------------
try:
    infile = sys.argv[sys.argv.index("--in-vcf")+1]
    print("sample to treat:", sample_ID)
except:
    usage()
    print("ERROR : Please enter the input vcf", infile)
    sys.exit()
try:
    outfile = sys.argv[sys.argv.index("--out-vcf")+1]
except:
    usage()
    print("ERROR : please enter the path to out vcf file", outfile)
    sys.exit()


vcf = pd.read_csv(infile,skiprows=56,sep="\t")

for row in vcf:

        fields=vcf.split(":")
        for i in range(1,len(fields)):
                if re.match(fields[i],"GT:DP:AF"):
                        data =  fields[i].split(':')
                        DP=data.index("DP")
                        AD=data.index("AD")
                #print(AD)
                # The comma ',' separates the AD of the different allele.
                # We want to take the last one that is the one of the mutation allele.
                        j=i+1
                        field=fields[j]
                        data=field.split(":")
                        AD = float(data[AD].split(","))
                        DP = float(data[DP])
                        AF = round(AD/DP,4)
                        #---- Replacing AF in the INFO column
                        INFO = fields[i]+":AF"
                        field= field+(":").join(INFO)
                        fields[i]=field
                        i=len(fields)
        row=("\t").join(fields)
        write_csv(row,output)

print("Output written")
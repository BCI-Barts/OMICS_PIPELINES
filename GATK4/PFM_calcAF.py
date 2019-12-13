#!/usr/bin/env python3

####### Author
# Pauline Fourgoux p.f.m.fourgoux@qmul.ac.uk ; pauline.fourgoux@gmail.com

####### Packages call
import pandas as pd
import sys
import os

####### Funtions
def usage():
    print("""
    computing_RNAseq_AF.py takes the ID of the sample to open a vcf file and calculates the real AF for RNA-seq mutation calling.

    exemple of command line:
    =======================

    python computing_RNAseq_AF.py --id  1360_T --path /path/to/vcf/

    obligatory:
    ==========

    --id -> sample ID

    --path -> path to the vcf file

    """)

####### Main

# Get arguments

try:
    sample_ID = sys.argv[sys.argv.index("--id")+1]
    print("sample to treat:", sample_ID)
except:
    usage()
    print("ERROR : please enter the sample ID")
    sys.exit()
try:
    path = sys.argv[sys.argv.index("--path")+1]
    print("path to file:", path)
except:
    usage()
    print("ERROR : please enter the path to vcf file")
    sys.exit()


inputfile = "output_f_" + sample_ID + ".vcf"
vcf_path = os.path.join(path,inputfile)

vcf = pd.read_table(vcf_path,skiprows=56)
rows = vcf.shape[0]

for r in list(range(rows)):
        #print(vcf.loc[[r]])
        if vcf.at[r,'FORMAT'] == "GT:AD:DP:GQ:PL":
                # Retrieving AD and DP to calculate real AF
                #print(vcf.loc[[r]])
                sample = vcf.at[r,'sample']
                sample = sample.split(':')
                DP = float(sample[2])
                #print(DP)
                if DP>0:
                        AD = sample[1]
                        #print(AD)
                        # The comma ',' separates the AD of the different allele.
                        # We want to take the last one that is the one of the mutation allele.
                        if "," in AD:
                                AD = AD.split(',')
                                AD = AD[-1]
                        AD = float(AD)
                        AF = round(AD/DP,4)
                        # Replacing AF in the INFO column
                        INFO = vcf.at[r,'INFO']
                        INFO = INFO.split(';')
                        INFO[1] = "AF=" + str(AF)
                        #print(AF)
                        #print(INFO[1])
                        INFO = ';'.join(INFO)
                        vcf.at[r,'INFO'] = INFO
                        #print(vcf.loc[[r]])

outputfile = "output_AF_" + sample_ID + ".vcf"
output_vcf_path = os.path.join(path,outputfile)
print("Writing output file:",output_vcf_path)
vcf.to_csv(output_vcf_path, '\t', header=True, index=None)
print("Output written")



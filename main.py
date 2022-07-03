import sys,os
import pandas as pd
from csv import reader, writer

# #Reading csv file
# data = pd.read_csv('CSQ.csv')

# #result file open
# result = open('result.fasta',"w")
# results = writer(open('result.csv',"w"),delimiter=',')

# df = pd.DataFrame(data)

# #filtering data
# filter_data = df[df["Type_of_Variant"]=="missense_variant"]
# NG = filter_data[filter_data["Gene"] == "BF2"]
# #printing the filtered data
# print(NG)
# NG.to_csv('ng.csv',index=False)

# #converting csv to FASTA

# # for items in NG:
# #     result.writelines(">"+items[0]+items[1]+"\n")
# # result.close()
# with result as outfile:
#     for item in NG:
#         outfile.writelines(">"+item+"\n"+''+"\n")

# outfile.close()

import vcf
import allel

# data = vcf.Reader(open('am.vcf','r'))
# for record in data:
#     print(record)
callset = allel.read_vcf('am.vcf')
print(sorted(callset.keys()))
# print(callset['samples'])

gt = allel.GenotypeArray(callset['calldata/GT'])
gt.count_het(axis=1)
ac = gt.count_alleles()
print(ac)
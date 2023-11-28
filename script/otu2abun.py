import os, argparse
import pandas as pd
import re
#set the input file
parser = argparse.ArgumentParser(description="Simulate the abundance file from ASV or OTU file")
parser.add_argument('-i', '--inputFilePath', metavar='fasta file', required=True, help='Input ASV or OTU file with fasta formate.')
parser.add_argument('-o', '--outputFilePath', metavar='simulate abundance file', required=True, help='Output file with table file.')
#get the user type informations.
args = parser.parse_args()

#print(args.inputFilePath)
#Open the input file, fasta file.
print("Open file: ",args.inputFilePath)
file1=open(args.inputFilePath,'r')
out=open(args.outputFilePath,"w")
#Read the content of file.
Lines=file1.readlines()
#Set the inital array.
faname=[]

#Impelete the file by line.
for line in Lines:
    #remove \n content.
    #print(line.strip())
    #Using regular expression, match the >name.
    matchObj=re.match(r'^>(.*)',line.strip())
    if matchObj:
        #print (matchObj.group(1))
        faname.append(matchObj.group(1))

file1.close()
#Now we get fasta names.

#print(faname)
faname_len=len(faname)
print("The total number of elements of array is ",faname_len,"Get the name of fasta file",faname[0],faname[1],"...")
title=['otu']
numlab=list(range(1,faname_len+1))
for i, x in enumerate(numlab):
    numlab[i]='s'+str(x)
numlab=['otu']+numlab
title2="\t".join(str(x) for x in numlab)
title2=title2+"\n"
out.write(title2)

for i, name in enumerate(faname):
    #print(i,name)
    ahead=[0]*i
    atail=[0]*(faname_len-i-1)
    #print(ahead)
    #print(1)
    #print(atail)
    simulate=[name]
    simulate=simulate+ahead+[1]+atail
    #print(simulate)
    #print(len(simulate))
    #print("\t".join(str(x) for x in simulate))
    info="\t".join(str(x) for x in simulate)
    info=info+"\n"
    out.write(info)
print("The process is successfully!!!")
out.close()





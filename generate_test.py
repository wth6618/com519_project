import pandas as pd
from Bio import SeqIO
"""
used to generatae testing dataset
"""
def generate_test():
    arrayA = []
    arrayB = []
    with open("Challenge_9934185_scoring/test_validation_set.txt",'r') as file:
        line = file.readline().split('\t')
        while line[0] != '':
            #print(line)
            if line[0] == '###\n':
                line = file.readline().split('\t')
                continue
            if line[2][10] == 'A':
                arrayA.append(line[2])
            else:
                arrayB.append(line[2])
            line = line = file.readline().split('\t')
    return arrayA, arrayB

A,B = generate_test()

print(len(A))

for a in A:
    print(a)
transcript_fileA = SeqIO.parse(open('Challenge_9934185_A.chromosomes/A.transcripts.fasta'),'fasta')
transcript_fileB = SeqIO.parse(open('Challenge_9934185_C.chromosomes/C.transcripts.fasta'),'fasta')

file = open('testA_trans.fasta','w+')
fileB = open('testC_trans.fasta','w+')
for trans in transcript_fileA:
    name, seq = trans.id, str(trans.seq)
    if name in A:
        file.write(">" + name + "\n" + seq + "\n")

for trans in transcript_fileB:
    name, seq = trans.id, str(trans.seq)
    if name in B:
        fileB.write(">" + name + "\n" + seq + "\n")



import pandas as pd
from Bio import SeqIO
#df = pd.read_csv("Challenge_9934185_scoring/dev_validation_set.tsv",delimiter='\t')
def generate_test():
    arrayA = []
    arrayB = []
    with open("Challenge_9934185_scoring/dev_validation_set.txt",'r') as file:
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
transcript_fileB = SeqIO.parse(open('Challenge_9934185_B.chromosomes/B.transcripts.fasta'),'fasta')

file = open('testA_trans.fasta','w+')
fileB = open('testB_trans.fasta','w+')
for trans in transcript_fileA:
    name, seq = trans.id, str(trans.seq)
    if name in A:
        file.write(">" + name + "\n" + seq + "\n")

for trans in transcript_fileB:
    name, seq = trans.id, str(trans.seq)
    if name in B:
        fileB.write(">" + name + "\n" + seq + "\n")
#with open('testB_trans.fasta','w+') as file:
#    for name_B in B:
#        for trans in transcript_fileB:
#            name, seq = trans.id, str(trans.seq)
#            if name == name_B:
#                file.write(">" + name_B + "\n" +seq + "\n")



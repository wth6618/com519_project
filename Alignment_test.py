import mappy as mp
import timeit
from Bio import SeqIO
import gffutils
import subprocess
import os
FeatureDB_A = gffutils.FeatureDB('saved_data/A_GFFDB')
FeatureDB_B = gffutils.FeatureDB('saved_data/B_GFFDB')
FeatureDB_C= gffutils.FeatureDB('saved_data/C_GFFDB')
transcript_file = SeqIO.parse(open('Challenge_9934185_A.chromosomes/A.chromosomes.fasta'),'fasta')
for trans in transcript_file:
        name, seq = trans.id, str(trans.seq)
        print(name)
        if name == 'chromosomeA1':
                input = seq
                break



#print(input)
a = mp.Aligner(seq= input,preset='map-ont')
#a = mp.Aligner("tryB_transcripts.fasta",preset='map-ont')  # load or build index
if not a: raise Exception("ERROR: failed to load/build index")
#s = a.seq("MT_human", 100, 200)     # retrieve a subsequence from the index
#mp.revcomp(s)
start = timeit.default_timer()

for name, seq, qual in mp.fastx_read("tryB_transcripts.fasta"): # read a fasta/q sequence
        for hit in a.map(seq): # traverse alignments

                print(name)
                print('source len {}'.format(len(seq)))
                print('target alignment info')
                print("{}\n{}\n{}\n{}".format(hit.ctg, hit.r_st, hit.r_en,hit.cigar_str))
                length = len(seq)
                print('target length {}'.format(length))
                #print(hit.ctg_len, hit.r_en-hit.r_st,len(seq))
                #print('blen {}'.format(hit.blen))
                #print('mlen {}'.format(hit.mlen))
                print('NM {}\n'.format(hit.NM))

                accuracy1NM = 1.0 - float(hit.NM/ length)
                #accuracyb = float(hit.blen/length)
                #accuracym = float(hit.mlen/length)
                print(accuracy1NM,hit.mapq)
                start, end = hit.r_st-10, hit.r_en+10
                seqid = 'chromosomeA1'
                if accuracy1NM > 0.89:
                        print(start,end)
                        region = (seqid,start,end)
                        features = FeatureDB_A.region(region,featuretype='mRNA')
                        count = 0
                        for f in features:
                                count+=1
                                print(f)
                        if count == 1: print('unique_transcript')








stop = timeit.default_timer()
print('Time: ', stop - start)
import mappy as mp
import timeit
import subprocess
import os
from My_Class import *
from util import *
import numpy as np
from Bio import SeqIO
from Pre_process import preprocessing
import distance
from alignment_TG import Alignment_TG
title = '#SourceA	SourceB	SourceA_Transcript_ID	SourceB_Transcript_ID	Call	Score	SourceA_Gene	SourceB_Gene	Category\n'



FeatureDB_A,FeatureDB_B,FeatureDB_C = preprocessing()
print('finish preprocess')

"""
Session: A & B
start Mapping B to A"""
print('Session: A & B\n start Mapping B to A\n')
start = timeit.default_timer()
file_handle = open('Outputs/output_BA.txt','w+')
file_handle.write(title)


output = Alignment_TG(FeatureDB_B,FeatureDB_A,'B','A',file_handle,source_file='Challenge_9934185_B.chromosomes/B.transcripts.fasta',
                                                      target_file='Challenge_9934185_A.chromosomes/A.chromosomes.fasta')
np.save('saved_data/BA.npy',output,allow_pickle=True)

print('start Mapping A to B')
output = Alignment_TG(FeatureDB_A,FeatureDB_B,'A','B',file_handle,source_file='Challenge_9934185_A.chromosomes/A.transcripts.fasta',
                                                      target_file='Challenge_9934185_B.chromosomes/B.chromosomes.fasta')
np.save('saved_data/AB.npy',output,allow_pickle=True)
file_handle.close()
print('done')

stop = timeit.default_timer()
print('Time: ', stop - start)
"""
Session END


Session: A & C
start Mapping C to A"""

print('Session: A & C\n start Mapping C to A\n')
start = timeit.default_timer()


file_handle = open('Outputs/output_AC.txt','w+')
file_handle.write(title)
output = Alignment_TG(FeatureDB_C,FeatureDB_A,'C','A',file_handle,source_file='Challenge_9934185_C.chromosomes/C.transcripts.fasta',
                                                      target_file='Challenge_9934185_A.chromosomes/A.chromosomes.fasta')
np.save('saved_data/CA.npy',output,allow_pickle=True)
print('start Mapping A to C\n')
output = Alignment_TG(FeatureDB_A,FeatureDB_C,'A','C',file_handle,source_file='Challenge_9934185_A.chromosomes/A.transcripts.fasta',
                                                      target_file='Challenge_9934185_C.chromosomes/C.chromosomes.fasta')
np.save('saved_data/AC.npy',output,allow_pickle=True)
file_handle.close()
print('done')

stop = timeit.default_timer()
print('Time: ', stop - start)

"""
Session END


Session: B & C
start Mapping B to C"""
print('Session: B & C\n start Mapping B to C\n')
start = timeit.default_timer()

file_handle = open('Outputs/output_BC.txt','w+')
file_handle.write(title)
output = Alignment_TG(FeatureDB_B,FeatureDB_C,'B','C',file_handle,source_file='Challenge_9934185_B.chromosomes/B.transcripts.fasta',
                                                      target_file='Challenge_9934185_C.chromosomes/C.chromosomes.fasta')
np.save('saved_data/BC.npy',output,allow_pickle=True)
print('start Mapping C to B\n')
output = Alignment_TG(FeatureDB_C,FeatureDB_B,'C','B',file_handle,source_file='Challenge_9934185_C.chromosomes/C.transcripts.fasta',
                                                      target_file='Challenge_9934185_B.chromosomes/B.chromosomes.fasta')
np.save('saved_data/CB.npy',output,allow_pickle=True)
file_handle.close()
print('done')

stop = timeit.default_timer()
print('Time: ', stop - start)
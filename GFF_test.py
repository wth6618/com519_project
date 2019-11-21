import os
cwd = os.getcwd()
cwd2 = cwd+'/GFFs'
cwd +='/gffcompare'


#os.chdir(cwd)
os.system('pwd')
#os.system('export PATH=$PATH:/home/wth6618/PycharmProjects/innocentive_Comp519-416/gffcompare')
#os.chdir(cwd2)
#os.system('pwd')
os.system('gffcompare -R -r ./Challenge_9934185_B.chromosomes/B.gff3 -o cuffcmp ./Challenge_9934185_C.chromosomes/C.gff3')

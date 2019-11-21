from My_Class import Exons
from My_Class import Transcript
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def write_title(file_handle):
    file_handle.write(
        '#SourceA	SourceB	SourceA_Transcript_ID	SourceB_Transcript_ID	Call	Score	SourceA_Gene	SourceB_Gene	Category\n')
    return

def find_Transcript(Chrom,trans_name):
    if Chrom == 'A':
        file = 'Challenge_9934185_A.chromosomes/A.transcripts.fasta'
    elif Chrom == 'B':
        file = 'Challenge_9934185_B.chromosomes/B.transcripts.fasta'
    else:
        file = 'Challenge_9934185_C.chromosomes/C.transcripts.fasta'
    transcript_file = SeqIO.parse(open(file),'fasta')
    for trans in transcript_file:
        name,seq = trans.id, str(trans.seq)
        print(name)
        if name == trans_name:
            return seq
def find_gene(Chrom,chromo_name,g_start,g_end):
    if Chrom == 'A':
        file = 'Challenge_9934185_A.chromosomes/A.chromosomes.fasta'
    elif Chrom == 'B':
        file = 'Challenge_9934185_B.chromosomes/B.chromosomes.fasta'
    else:
        file = 'Challenge_9934185_C.chromosomes/C.chromosomes.fasta'
    chromo_file = SeqIO.parse(open(file),'fasta')
    for chrom in chromo_file:
        name,seq = chrom.id, str(chrom.seq)

        if name == chromo_name:
            gene = seq[g_start:g_end]
            print(name)
            return gene
def get_chromo(trans):
    return trans[10]

def get_allexon_seq(trans):
    t_name = trans.name
    #t_seq = find_Transcript(get_chromo(t_name),t_name)
    t_seq = find_gene(get_chromo(t_name),trans.chrom,trans.start,trans.end)
    print(len(t_seq))
    exon_seqs = []
    if len(trans.exons) !=0:
        for exon in trans.exons:
            print(exon.start,exon.end)
            exon_seqs.append(t_seq[exon.start: exon.end+1])

    return exon_seqs


def Add_new_trans(trans_name,All_trans,file):
    GFF_handle = open(file,'r')
    line = GFF_handle.readline().split()
    line = GFF_handle.readline().split()
    trans_find = False
    while not trans_find:
        if line[2] == 'gene':
            chromosome = line[0]
            gene = line[-1].split(';')[0][3:]
            line = GFF_handle.readline().split()
            while line[2] !='gene':
                if line[2] == 'mRNA':
                    transcript = line[-1].split(';')[0][3:]
                    #print(transcript)
                    if transcript == trans_name:
                        c_start = int(line[3])
                        c_end = int(line[4])
                        #print(chromosome, gene, transcript)
                        print('add')
                        print(line)
                        line = GFF_handle.readline().split()
                        print(line)
                        new_transcript = extractGFF_info(chromosome,gene,c_start,c_end,GFF_handle,transcript)
                        All_trans[transcript] = new_transcript
                        trans_find = True
                    else:
                        line = GFF_handle.readline().split()
                else:
                    line = GFF_handle.readline().split()

        else:
            line = GFF_handle.readline().split()

    return




def extractGFF_info(chromosome,gene, c_start,c_end,handler,transcript,findseq = False):


    line = handler.readline().split()
    exons = []
    print('extract')
    print(line)
    while line[2] == 'intron':
        line = handler.readline().split()

    while line[2] == 'exon':
        start = int(line[3])-c_start
        end = int(line[4]) - c_start
        name = line[-1].split(';')[0][3:]
        parent = line[-1].split(';')[-1][7:]
        new_exon = Exons(name, start,end,chromosome,transcript,gene)
        #new_exon.print_info()
        exons.append(new_exon)

        line = handler.readline().split()

    if findseq:
        seq = find_Transcript('B',transcript)
        new_transcript = Transcript(chromosome, transcript, c_start, c_end, exons,seq)
        print(seq)
    else:
        new_transcript = Transcript(chromosome, transcript, c_start, c_end, exons)

    return new_transcript


def Compare_exons(trans1, trans2):
    # check number of exons
    print(len(trans1.exons),len(trans2.exons))
    if len(trans1.exons) != len(trans2.exons):
        print(trans1.exons)

        print('different length')
        category = 'changed_exons'
        return category,False
    else: return 'exact_match',True
    # check content by alignment


"""    exonset1 = get_allexon_seq(trans1)
    exonset2 = get_allexon_seq(trans2)
    for i in range(len(trans1.exons)):
        exon1 = exonset1[i]
        exon2 = exonset2[i]
        alignment_result = pairwise2.align.globalms(exon1, exon2, 1.5, -1, -1, -0.1,score_only = True)
        score = min(len(exon2),len(exon1))
        if alignment_result > score:
            score = alignment_result
            print(score)
            return True
        else:
            print('not similar')
            return False"""

        #for a in alignment_result:
        #    print(format_alignment(*a))



#All_trans = dict()
#Add_new_trans('transcriptB3',All_trans)
#Add_new_trans('transcriptB1',All_trans)
#Add_new_trans('transcriptB5',All_trans)
#Add_new_trans('transcriptB2',All_trans)
#print(All_trans)
#print(All_trans['transcriptB2'].name)

#print(get_allexon_seq(All_trans['transcriptB2']))
#Compare_exons(All_trans['transcriptB1'],All_trans['transcriptB3'])
#GFF_handle = open('exon_try.txt','r')
#line = GFF_handle.readline().split()
#print(line)
def GFF_get_gene(C,transcript_name):
    if C =='A':
        file = 'Challenge_9934185_A.chromosomes/A.gff3'
    elif C =='B':
        file = 'Challenge_9934185_B.chromosomes/B.gff3'
    else:
        file = 'Challenge_9934185_C.chromosomes/C.gff3'
    #print('name {}'.format(transcript_name))
    with open(file,'r') as GFF:
        find = False
        line = GFF.readline()
        #print(line)

        line = GFF.readline().split()
        #print(line)
        while not find:
            #print(line)
            if line[2] == 'gene':

                gene = line[-1].split(';')[0][3:]
                line = GFF.readline().split()
                while line[2] != 'gene':
                    if line[2] == 'mRNA':

                        transcript = line[-1].split(';')[0][3:]
                        if transcript == transcript_name:
                            find = True
                            return gene
                        else:line = GFF.readline().split()
                    else:
                        line = line = GFF.readline().split()
            else:
                line =line = GFF.readline().split()

#print(GFF_get_gene('C','transcriptC43470'))

"""        exons = []
        while line[2] == 'exon':
            start = int(line[3])-c_start
            end = int(line[4]) - c_start
            name = line[-1].split(';')[0][3:]
            parent = line[-1].split(';')[-1][7:]
            new_exon = Exons(name, start,end,chromosome,parent)
            new_exon.print_info()
            exons.append(new_exon)

            line = GFF.readline().split()
        new_transcript = Transcript()
        print(GFF.readline().split())
        print(GFF.readline().split())"""
class Exons(object):
    def __init__(self, name, start,end, chromo,parent,gene):
        self.name = name
        self.gene = gene
        self.start = start
        self.end = end
        self.chromo = chromo
        self.parent = parent
    def print_info(self):
        print('exon_name:{},(start,end): ({},{}),chromo: {}, parent: {}'.format(
            self.name,self.start, self.end,self.chromo,self.parent))

    def get_sequence(self):

        return

class Transcript(object):
    def __init__(self,Chromosome,name,start,end,exons,gene, Seq = None):
        self.chrom = Chromosome
        self.gene = gene
        self.name = name
        self.start = start
        self.end = end
        self.seq = Seq
        self.exons = exons

"""SourceA	SourceB	SourceA_Transcript_ID	SourceB_Transcript_ID	Call	Score	SourceA_Gene	SourceB_Gene	Categ"""
class Match(object):
    def __init__(self,SourceA,SourceB,SourceA_Transcript_ID,SourceB_Transcript_ID,Call,Score,SourceA_Gene,SourceB_Gene,Categ):
        self.sourceA = SourceA
        self.sourceB = SourceB
        self.SourceA_Transcript_ID = SourceA_Transcript_ID
        self.SourceB_Transcript_ID = SourceB_Transcript_ID
        self.call = Call
        self.score = Score
        self.SourceA_Gene = SourceA_Gene
        self.SourceB_Gene = SourceB_Gene
        self.category = Categ
    def print_line_result(self,hanlde):
        output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.sourceA,self.sourceB,self.SourceA_Transcript_ID,
                                                               self.SourceB_Transcript_ID,self.call,self.score,
                                                               self.SourceA_Gene,self.SourceB_Gene,self.category)
        print(output)
        hanlde.write(output)
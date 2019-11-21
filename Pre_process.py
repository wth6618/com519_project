import gffutils
from pathlib import Path
from Bio import SeqIO
import mappy as mp

#db = gffutils.create_db('Challenge_9934185_A.chromosomes/A.gff3', 'saved_data/A_GFFDB', keep_order=True,merge_strategy='merge', sort_attribute_values=True,id_spec={"gene": ["ID","Name"], "mRNA": ["ID", "transcript_id"],"exon":["ID","Name"],"CDS":["ID","Name"]})

#db = gffutils.create_db('Challenge_9934185_A.chromosomes/A.gff3', 'saved_data/A_GFFDB',force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
#FeatureDB_A = gffutils.FeatureDB('saved_data/A_GFFDB')
#FeatureDB_B = gffutils.FeatureDB('saved_data/B_GFFDB')
#Feature_db = gffutils.FeatureDB('saved_data/B_GFFDB')
#print('db finished')
#db1 = gffutils.create_db('Challenge_9934185_B.chromosomes/B.gff3', 'saved_data/B_GFFDB', keep_order=True,merge_strategy='merge', sort_attribute_values=True,id_spec={"gene": ["ID","Name"], "mRNA": ["ID", "transcript_id"],"exon":["ID","Name"]})
#print('db1 finished')

#db2 = gffutils.create_db('Challenge_9934185_C.chromosomes/C.gff3', 'saved_data/C_GFFDB', keep_order=True,merge_strategy='merge', sort_attribute_values=True,id_spec={"gene": ["ID","Name"], "mRNA": ["ID", "transcript_id"],"exon":["ID","Name"],"CDS":["ID","Name"],"intron":["ID","Name"]})


#print('db2 finished')
def preprocessing():
    A_db = Path("saved_data/A_GFFDB")
    if A_db.is_file():
        FeatureDB_A = gffutils.FeatureDB('saved_data/A_GFFDB')
    else:
        print('generateing data')
        db = gffutils.create_db('Challenge_9934185_A.chromosomes/A.gff3', 'saved_data/A_GFFDB', keep_order=True,
                                merge_strategy='merge', sort_attribute_values=True,
                                id_spec={"gene": ["ID", "Name"], "mRNA": ["ID", "transcript_id"],
                                         "exon": ["ID", "Name"], "CDS": ["ID", "Name"]})
        print('done')
        FeatureDB_A = gffutils.FeatureDB('saved_data/A_GFFDB')
    B_db = Path("saved_data/B_GFFDB")
    if B_db.is_file():
        FeatureDB_B = gffutils.FeatureDB('saved_data/B_GFFDB')
    else:
        db = gffutils.create_db('Challenge_9934185_B.chromosomes/B.gff3', 'saved_data/B_GFFDB', keep_order=True,
                                merge_strategy='merge', sort_attribute_values=True,
                                id_spec={"gene": ["ID", "Name"], "mRNA": ["ID", "transcript_id"],
                                         "exon": ["ID", "Name"], "CDS": ["ID", "Name"]})
        FeatureDB_B = gffutils.FeatureDB('saved_data/B_GFFDB')
    C_db = Path("saved_data/C_GFFDB")
    if B_db.is_file():
        FeatureDB_C = gffutils.FeatureDB('saved_data/C_GFFDB')
    else:
        db = gffutils.create_db('Challenge_9934185_C.chromosomes/C.gff3', 'saved_data/C_GFFDB', keep_order=True,
                                merge_strategy='merge', sort_attribute_values=True,
                                id_spec={"gene": ["ID", "Name"], "mRNA": ["ID", "transcript_id"],
                                         "exon": ["ID", "Name"], "CDS": ["ID", "Name"]})
        FeatureDB_C = gffutils.FeatureDB('saved_data/C_GFFDB')
    return (FeatureDB_A,FeatureDB_B,FeatureDB_C)


'''children = Feature_db.children('geneB2',featuretype=('exon','mRNA'))
parents = Feature_db.parents('transcriptB5.exon.2',featuretype='gene')
print(Feature_db['transcriptB5.exon.2'])
#schema = Feature_db.schema()
#print(schema)
for child in children:
    print(child)

parents = list(parents)
print(parents[0].start)
for p in parents:
    print('hi')
    print(p)

for feature in Feature_db.featuretypes():
    print(feature)'''

'''exon = Feature_db['transcriptB5.exon.2']
#print(Feature_db.schema())
print(exon.seqid)
seq = exon.sequence('Challenge_9934185_B.chromosomes/B.chromosomes.fasta')
print(len(seq))'''

'''FeatureDB_C = gffutils.FeatureDB('saved_data/C_GFFDB')
start, end = 16324190, 16325199
seqid = 'chromosomeC1'
region = (seqid, start, end)
features = FeatureDB_C.region(region)
for f in features:
    print(f)'''

'''featureA = FeatureDB_A['transcriptA80928']
featureB = FeatureDB_B['transcriptB22040']
print('A: start {}, end{}'.format(featureA.start,featureB.end))
for child in FeatureDB_A.children('transcriptA45689',featuretype='exon'):
    print(child)
print('B: start {}, end{}'.format(featureB.start,featureB.end))

for child in FeatureDB_B.children('transcriptB21774',featuretype='exon'):
    print(child)'''

"""db = gffutils.create_db('Challenge_9934185_B.chromosomes/B.gff3', 'saved_data/B_GFFDB', keep_order=True,
                        merge_strategy='merge', sort_attribute_values=True,
                        id_spec={"gene": ["ID", "Name"], "mRNA": ["ID", "transcript_id"],
                                 "exon": ["ID", "Name"], "CDS": ["ID", "Name"]})
print('done')"""
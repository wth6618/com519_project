import numpy as np
import distance


outputs = np.load('saved_data/BA.npy',allow_pickle=True).item()
write_file = open('Outputs/BA_solution.tsv','w+')

print(outputs)
u_count = 0
g_count = 0
for key in outputs.keys():
    for M in outputs[key]:
        if M[0] == 'unique_transcript':
            u_count +=1
        if M[0] == 'gene_fusion':
            g_count +=1

    if u_count ==0:
        if g_count >0:
            # write gene_fusion

            output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.sourceA, self.sourceB,
                                                                   self.SourceA_Transcript_ID,
                                                                   self.SourceB_Transcript_ID, self.call, self.score,
                                                                   self.SourceA_Gene, self.SourceB_Gene, self.category)





"""def check_features(name_S,features_T,low_boundary, upper_boundary,FeatureDB_T,FeatureDB_S,count):
    f_arrary = []
    for f in features_T:
        f_arrary.append(f)
        print(f)

    if len(f_arrary) == 1:
        if f_arrary[0].start> low_boundary and f_arrary[0].end < upper_boundary:
            print(f_arrary[0].id)
            #parent = FeatureDB.parents(f_arrary[0].id,featuretype = 'gene')
            return 'unique_transcript', 'exact_match',f_arrary[0].id
        children_T = list(get_children(f_arrary[0].id, FeatureDB_T, 'exon'))
        children_S = list(get_children(name_S, FeatureDB_S, 'exon'))
        # check if they have same number of exons
        if len(list(children_T)) != len(list(children_S)) :
            if count > 0:
                return '','',''
            else:
                return 'absent_transcript','',''
        # check exon boundary
        else:

            if (children_T[0].end <low_boundary) or (children_T[-1].start >upper_boundary):
                if count > 0:
                    return '', '', ''
                else:
                    return 'absent_transcript', '', ''
            else:
                print(f_arrary[0].id)

                return 'unique_transcript', '', f_arrary[0].id


    else:
        print('length not 1')
        return '','','' 
"""
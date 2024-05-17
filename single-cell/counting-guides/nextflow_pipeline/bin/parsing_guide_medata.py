#!/usr/bin/env python

import pandas as pd
import sys

#Simulating a chromossome with each guides
#Creating a bed files and pseudo isoform files to generate a gtf to be used to create a guide star referece
df = pd.read_excel(sys.argv[1]) # This is the script input guide_metadata from nextflow
df['protospacer'] = df['sgRNA_sequences'].apply(lambda x:   x)  # Probably this will be not needed case we use the correct reference
df_metadata_raw = df.copy()
dict_proto_to_guide_id = dict(zip(df_metadata_raw['protospacer'].values, df_metadata_raw['sgRNA_ID'].values))
bed_pd = pd.DataFrame()
bed_pd['seqid'] =  [f"chr_{dict_proto_to_guide_id[x]}_peseudo_chr_{x}"   for x in  df_metadata_raw['protospacer'] ] 
bed_pd['star'] = 0
bed_pd['end'] = df_metadata_raw['protospacer'].apply(len) 
bed_pd['source'] =  df_metadata_raw['protospacer']
bed_pd['type'] = df_metadata_raw['protospacer']
bed_pd['strand'] = '+'
bed_pd['thickStart'] = 1
bed_pd['thickEnd'] = bed_pd['end'] 
bed_pd['color'] = "0,0,255"
bed_pd['blockCount'] = 1
bed_pd['blockSizes'] = bed_pd['end'] 
bed_pd['blockStarts'] = 0
bed_pd.to_csv('guides.bed', header=None, index=None, sep='\t')
bed_pd[['seqid', 'source']].to_csv('isoforms.txt', header=None, index=None, sep='\t')
#creating the pseud_genome from the guide sequences
with open('pseudo_genome.fa', 'w') as save_genome:
    save_genome.write('\n'.join(bed_pd.apply(lambda x : f'>{x["seqid"]}\n{x["source"]}' ,axis=1 ).tolist()))

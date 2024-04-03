#!/usr/bin/env python
import pandas as pd
import sys
import anndata
import os
import subprocess

# def create_ann_data(dir_files, type_file):
#     dir_files
#     detecting_mtx_cmd = f"find  {dir_files} |grep matrix.mtx | grep {type_file}"
#     print (detecting_mtx_cmd)
#     detecting_mtx =!$detecting_mtx_cmd
#     print (detecting_mtx)
#     bar_to_extract =  detecting_mtx[0]
#     print (f'reading matrix...{bar_to_extract}')
#     test_ann = anndata.read_mtx(bar_to_extract)
#     bar_names_cmd = f"find  {dir_files}      | grep barcodes.tsv |grep {type_file}"
#     print (f'{bar_names_cmd}')

#     bar_names = !$bar_names_cmd
#     bar_names_cat = bar_names[0]
#     bars =  !cat $bar_names_cat
#     features_file_cmd = f"find  {dir_files}  | grep features.tsv | grep {type_file}"
#     print (f'{features_file_cmd}')

#     features_file = !$features_file_cmd
#     feature_selected = features_file[0]
#     feature = !cat $feature_selected
    
#     print (pd.DataFrame(feature).head())
#     test_ann.obs.index = [ x for x in  pd.DataFrame([i.split('\t') for i in feature])[1]]

#     #print (bars)
#     test_ann.var.index = [ x.split('-')[0] for x in  pd.DataFrame(bars)[0]]
#     return bars, test_ann



def run_system_command(command):
    """
    Execute a system command and return the output.
    """
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stderr:
        print(f"Error executing command '{command}': {result.stderr}")
    return result.stdout.strip().split('\n')

def create_ann_data(dir_files, type_file):
    detecting_mtx_cmd = f"find {dir_files} -name 'matrix.mtx' -and -path '*{type_file}*' "
    print(detecting_mtx_cmd)
    detecting_mtx = run_system_command(detecting_mtx_cmd)
    print(detecting_mtx)
    bar_to_extract = detecting_mtx[0]
    print(f'reading matrix...{bar_to_extract}')
    test_ann = anndata.read_mtx(bar_to_extract)

    bar_names_cmd = f"find {dir_files} -name 'barcodes.tsv' -and -path '*{type_file}*' "
    print(bar_names_cmd)
    bar_names = run_system_command(bar_names_cmd)
    bar_names_cat = bar_names[0]
    bars = run_system_command(f"cat {bar_names_cat}")

    features_file_cmd = f"find {dir_files} -name 'features.tsv' -and -path '*{type_file}*' "
    print(features_file_cmd)
    features_file = run_system_command(features_file_cmd)
    feature_selected = features_file[0]
    feature = run_system_command(f"cat {feature_selected}")
    
    print(pd.DataFrame(feature).head())
    test_ann.obs.index = [x for x in pd.DataFrame([i.split('\t') for i in feature])[1]]

    test_ann.var.index = [x.split('-')[0] for x in pd.DataFrame(bars)[0]]
    return bars, test_ann



bar, ann_out = create_ann_data(sys.argv[1], 'raw')
df = pd.read_excel(sys.argv[2]) # This is the script input guide_metadata from nextflow
df['protospacer'] = df['sgRNA_sequences'].apply(lambda x:   x)  # Probably this will be not needed case we use the correct reference
seq_to_original_id = dict(zip(df['sgRNA_sequences'].values, df['sgRNA_ID'].values)) # i think it should be 
ann_out.obs.index = [seq_to_original_id[x.split('_')[-1]] for x in ann_out.obs.index]
ann_out.write_h5ad('guides_star_solo.h5ad')

#this should be a plot module
# df_guides_total = pd.DataFrame( { 'guides_id': ann_out.obs.index, 'total_UMIs': ann_out.X.sum(1).flatten().tolist()[0] })
# df_guides_total['target'] =  df_guides_total['guides_id'].apply(lambda x: x.split('_')[0]  )
# df_guides_total   = df_guides_total.sort_values('guides_id')
# df_guides_total.index = df_guides_total['guides_id'].values
# df_guides_total.query('target == "BCL11A" ').plot(kind='bar')
#plt.savefig('test.png')

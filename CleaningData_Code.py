#Importing Libraries 
import pandas as pd
import numpy as np
import os

#File Paths
input_file_path = r"C:\Users\sauri\OneDrive\Desktop\Fem-1 Project\GSE174646_flp21_wt_early_embryo_rawcounts_GEO.txt"
output_directory = r"C:\Users\sauri\OneDrive\Desktop\Fem-1 Project"
output_file_name = "cleaned_GSE174646.txt"
output_file_path = os.path.join(output_directory, output_file_name)

#Importing Files for Analayis
wt_columns = [
    'WT_11_mRNA_trimmed.bam',
    'WT_22_mRNA_trimmed.bam',
    'WT_30_mRNA_trimmed.bam'
]
gene_id_column = 'Geneid'

#Using Pandas to delete the columns not needed in the file
df = pd.read_csv(input_file_path, sep='\t', skiprows=2, header=0, comment='#')

df_selected = df[[gene_id_column] + wt_columns].copy()

def calculate_positional_stats(row):
    wt_values_str = [str(row[col]) for col in wt_columns]

    wt_lists_numeric = []
    max_len = 0
    for val_str in wt_values_str:
        parts = val_str.split(';')
        numeric_parts = []
        for p in parts:
            try:
                numeric_parts.append(float(p))
            except ValueError:
                numeric_parts.append(np.nan)
        wt_lists_numeric.append(numeric_parts)
        if len(numeric_parts) > max_len:
            max_len = len(numeric_parts)

    for i in range(len(wt_lists_numeric)):
        if len(wt_lists_numeric[i]) < max_len:
            wt_lists_numeric[i].extend([np.nan] * (max_len - len(wt_lists_numeric[i])))

    np_wt_values_2d = np.array(wt_lists_numeric, dtype=object)

    positional_means = []
    positional_stds = []

    for pos_idx in range(max_len):
        current_position_values = np_wt_values_2d[:, pos_idx]

        numeric_pos_values = current_position_values[~pd.isna(current_position_values)].astype(float)

        if len(numeric_pos_values) > 0:
            positional_means.append(np.nanmean(numeric_pos_values))
            positional_stds.append(np.nanstd(numeric_pos_values))
        else:
            positional_means.append(np.nan)
            positional_stds.append(np.nan)

#using numpy to calculate mean and standard deviation of wildstrains
    mean_parts = [f"{x:.4f}" if not np.isnan(x) else 'NaN' for x in positional_means]
    std_parts = [f"{x:.4f}" if not np.isnan(x) else 'NaN' for x in positional_stds]

    mean_str = ';'.join(mean_parts)
    std_str = ';'.join(std_parts)

    if max_len > 0 and not mean_str.endswith(';'):
        mean_str += ';'
    if max_len > 0 and not std_str.endswith(';'):
        std_str += ';'

    return pd.Series([mean_str, std_str])

df_selected[[f'Mean of {gene_id_column} Strains', f'Standard Deviation of {gene_id_column} Strains']] = df_selected.apply(
    calculate_positional_stats, axis=1
)
#outputting final file
df_output = df_selected[[gene_id_column, f'Mean of {gene_id_column} Strains', f'Standard Deviation of {gene_id_column} Strains']].copy()

df_output.to_csv(output_file_path, sep='\t', index=False)
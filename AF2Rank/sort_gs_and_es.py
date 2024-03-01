import pandas as pd
import sys

output_data_xlsx = "/Users/chenea/work/VSCode/af2rank_0923/af2rank_output_data.xlsx"

output_df = pd.read_excel(output_data_xlsx)
column_names = [column for column in output_df if column != 'Unnamed: 0']
columns = [output_df[column] for column in output_df if column != 'Unnamed: 0']
ground_states = columns[-1] #remove last column

gs_data = []
for i in range(len(ground_states)):
    new_line = [column[i] for column in columns]
    if ground_states[i] == columns[7][i]:
        new_line = new_line[7:14]+new_line[:7]+[new_line[14]]
    gs_data.append(new_line)

excel_writer = gs_data
output_df = pd.DataFrame(excel_writer,
                        columns=["pdb_gs", "fsr_tm_gs", "fsr_plddt_gs", "tm_gs", "plddt_gs", "ptm_gs", "composite_gs",
                        "pdb_es", "fsr_tm_es", "fsr_plddt_es", "tm_es", "plddt_es", "ptm_es", "composite_es", "ground_states"])
                        
output_df.to_excel("gs_es.xlsx")
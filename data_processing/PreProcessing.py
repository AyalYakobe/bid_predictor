import pandas as pd
import string
import numpy as np

dfFinal=pd.read_csv("LiuLloydGevers1_AbundanceData_WSpecies_WMetadata_IBD.csv", header=None, low_memory=False)


section = dfFinal.iloc[1:1169, 7:]
section = section.astype(float)
dfFinal.iloc[1:1169, 7:]=section
dfFinal=dfFinal.iloc[0:1169,:]
row_sums = dfFinal.iloc[1:,7:].sum(axis=1)
count=1
for entry in row_sums:
    if(int(entry)==0):
        dfFinal.drop(count,inplace=True)
    count+=1
row_sums = dfFinal.iloc[1:,7:].sum(axis=1)
# Normalize each row
dfFinal.iloc[1:,7:] = dfFinal.iloc[1:, 7:].div(row_sums, axis=0)
Study=[]
Study.append("Study")
for i in range(3,85):
    Study.append("Liu")
for i in range(85,261):
    Study.append("Lloyd")
for i in range(261,1167):
    Study.append("Gevers1")
dfFinal.insert(loc=4, column="4",value=Study)

liu_row_index = dfFinal[dfFinal.iloc[:, 4] == "Liu"].index[0]
Lloyd_row_index = dfFinal[dfFinal.iloc[:, 4] == "Lloyd"].index[0]
Gevers1_row_index = dfFinal[dfFinal.iloc[:, 4] == "Gevers1"].index[0]

# Calculate the sum for each column starting from the row where "Liu" is found
column_sums_Liu= dfFinal.iloc[1:83, 8:].sum()
column_sums_Lloyd= dfFinal.iloc[83:259, 8:].sum()
column_sums_Gevers1= dfFinal.iloc[259:, 8:].sum()

# Iterate over each column and sub NaNs if not found in that study
for col in dfFinal.columns[8:]:
    if column_sums_Liu[col] == 0.0:
        dfFinal.iloc[1:83, col] = np.nan
    if column_sums_Lloyd[col] == 0.0:
        dfFinal.iloc[83:259, col] = np.nan
    if column_sums_Gevers1[col] == 0.0:
        dfFinal.iloc[259:, col] = np.nan 


dfFinal.to_csv("LiuLloydGevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized_WStudy_WNA.csv", index=False, na_rep='NaN')
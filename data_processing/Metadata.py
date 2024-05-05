import pandas as pd
import string
import numpy as np

df = pd.read_csv('Liu/SraRunTable_Liu_Metadata.txt', header=None)
df=df.iloc[1:,:]

Metadata=[]

for entry in range(len(df)):
    if(df.iloc[entry,14]=="Crohn's disease"):
        status="IBD"
    elif(df.iloc[entry,14]=="Not IBD"):
        status="Control"
    else:
        status="N/A"
    if("Biopsy" in df.iloc[entry,26]):
        type= "Biopsy Ileum"
    else:
        type="Stool"
    if(df.iloc[entry,24]=="male"):
        sex="Male"
    else:
        sex="Female"
    subject=df.iloc[entry,-3]
    collectionWeek=0
    Metadata.append([df.iloc[entry,0],subject, collectionWeek, status, sex, df.iloc[entry,22], type])



df2 = pd.read_csv('Lloyd-Price(HMP2)/SraRunTable_Lloyd-Price_Metadata.txt', header=None)

df4=df2.iloc[1:,:]
filtered_df2 = df2[df2.iloc[:, 0].astype(str).str.startswith('SRR67')]

df5 = pd.read_csv('HMP2/hmp2_metadata_2018-08-20_Lloyd-Price_Metadata.csv', header=None, low_memory=False)
filtered_df5 = df5[df5.iloc[:, 4].astype(str).str.startswith('biopsy')]

for entry in range(len(filtered_df2)):
    row = filtered_df5[filtered_df5.iloc[:,1] == filtered_df2.iloc[entry,40]]
    if(filtered_df2.iloc[entry,38]=="Crohn''s disease" or filtered_df2.iloc[entry,38]=="CD"):
        status="IBD"
    elif(filtered_df2.iloc[entry,38]=="nonIBD"):
        status="Control"
    elif(filtered_df2.iloc[entry,38]=="Ulcerative Colitis" or filtered_df2.iloc[entry,38]=="UC"):
        status="IBD"
    else:
        status="N/A"
    type = "Biopsy "+ str(row.iloc[0,40])
    if(filtered_df2.iloc[entry,37]=="male"):
        sex="Male"
    else:
        sex="Female"
    age= row.iloc[0,69]
    subject=row.iloc[0,26]
    collectionWeek=row.iloc[0,5]
    Metadata.append([filtered_df2.iloc[entry,0],subject, collectionWeek, status, sex,age, type ])


df7 = pd.read_csv('Gevers1/1-s2.0-S1931312814000638-mmc2_Gevers1_Metadata.csv', header=None)

df8 = pd.read_csv('Gevers1/SraRunTable_Gevers1_Metadata.txt', header=None)
print(df8.iloc[:,:7])
df9 = pd.read_csv('Gevers1/1939_20230206-083256_Gevers1_Metadata.txt',sep='\t', header=None)

df10 = pd.read_csv('Gevers1/1998_20230206-084333_Gevers1_Metadata.txt',sep='\t', header=None)


UsefulAccessions={}
for entry in range(1,len(df7)):
    if(df7.iloc[entry,6]=="stool"):
        type="Stool"
    elif("Ileum" in df7.iloc[entry,6]):
        type= "Biopsy Ileum"
    elif("Rectum" in df7.iloc[entry,6]):
        type= "Biopsy Rectum"
    else:
        type="Biopsy nan"
    if(df7.iloc[entry,3]=="CD"):
        status="IBD"
    elif(df7.iloc[entry,3]=="Not IBD"):
        status="Control"
    elif(df7.iloc[entry,3]=="UC"):
        status="IBD"
    else:
        status="N/A"
    if(str(df7.iloc[entry,0])[0]=="S" and str(df7.iloc[entry,0])[2]!="B"):
        sample=str(df7.iloc[entry,0])[0:2]+"B"+str(df7.iloc[entry,0])[2:]
    else:
        sample=str(df7.iloc[entry,0])

    UsefulAccessions[sample]=[status,type,df7.iloc[entry,1]]


for entry in range(1,len(df9)):
    name=df9.iloc[entry,1]
    if(name[-1]=="S"):
        name=name[0:2]+"-"+name[3:]
    elif(name[0]=="S"):
        name = name.translate(str.maketrans('', '', string.punctuation))
    if(df9.iloc[entry,55]=="male"):
        sex="Male"
    else:
        sex="Female"
    if name in UsefulAccessions:
        if(len(UsefulAccessions[name])==3):
            UsefulAccessions[name].append(df9.iloc[entry,30])
            UsefulAccessions[name].append(sex)


for entry in range(1,len(df10)):
    name=str(df10.iloc[entry,1])
    if(name[-1]=="S"):
        name=name[0:2]+"-"+name[3:]
    elif(name[0]=="S" and name[2]!="B"):
        name = name[0:2]+"B"+name[2:].translate(str.maketrans('', '', string.punctuation))
    elif(name[0]=="S" and name[2]=="B" and name[-4].isnumeric()==False):
        name.translate(str.maketrans('', '', string.punctuation))
        name = name[0:-3]+"0"+name[-3:].translate(str.maketrans('', '', string.punctuation))
    elif(name[0]=="S"):
        name = name.translate(str.maketrans('', '', string.punctuation))
    if(df10.iloc[entry,53]=="male"):
        sex="Male"
    else:
        sex="Female"
    if name in UsefulAccessions:
        if(len(UsefulAccessions[name])==3):
            UsefulAccessions[name].append(df10.iloc[entry,28])
            UsefulAccessions[name].append(sex)


for entry in UsefulAccessions:
    if len(UsefulAccessions[entry])==3:
        UsefulAccessions[entry].append("N/A")
        UsefulAccessions[entry].append("N/A")
    UsefulAccessions[entry].append([])
    
SRAs350=set()
SRAs175=set()
for entry in range(1,len(df8)):
    if(str(df8.iloc[entry,3])=="AMPLICON"):
        name=str(df8.iloc[entry,44])
        if(name[0]=="P"):
            name=name.split(".",1)[1]
        if(name[-1]=="S"):
            name=name[0:6]+".S"
        elif(name[0]=="S"):
            name = name.translate(str.maketrans('', '', string.punctuation))
        if(name[0]=="S" and name[2]!="B"):
            name = name[0:2]+"B"+name[2:].translate(str.maketrans('', '', string.punctuation))
        elif(name[0]=="S" and name[2]=="B" and name[-4].isnumeric()==False):
            name = name[0:-3]+"0"+name[-3:].translate(str.maketrans('', '', string.punctuation))
        if name in UsefulAccessions:
            UsefulAccessions[name][5].append(df8.iloc[entry,0])
            if(int(df8.iloc[entry,4])==350):
                SRAs350.add(str(df8.iloc[entry,0]))
            elif(int(df8.iloc[entry,4])==175):
                SRAs175.add(str(df8.iloc[entry,0]))


addToMetadata=[]
countNA=0
for entry in UsefulAccessions:
    if(len(UsefulAccessions[entry][5])>0):
        for SRA in UsefulAccessions[entry][5]:
            if(UsefulAccessions[entry][4]!="N/A"):
                addToMetadata.append([SRA, UsefulAccessions[entry][2], 0, UsefulAccessions[entry][0],UsefulAccessions[entry][4],UsefulAccessions[entry][3],UsefulAccessions[entry][1]])
            else:
                countNA+=1
                SRAs350.discard(str(SRA))
                SRAs175.discard(str(SRA))


            

#SRAs = sorted(SRAs, key=lambda x: x[0])
SRAs350=sorted(SRAs350)
SRAs175=sorted(SRAs175)
addToMetadata=sorted(addToMetadata, key=lambda x: x[0])
Metadata.extend(addToMetadata)

filename = "SRAsUsed.txt"


# Writing the list elements to the text file
with open(filename, 'w') as f:
    for item in SRAs350:
        f.write("%s\n" % (item+"_1.fastq.gz"))
        f.write("%s\n" % (item+"_2.fastq.gz"))

df3 = pd.DataFrame(Metadata, columns=["Sample", "Subject", "Collection Week","Diagnosis", "Sex", "Age", "Sample Type"])
df3.to_csv("LiuLloydGeversMetadataIBD.csv", index=False)

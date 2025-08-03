# Usage: generateColor2.py coolerPath cpg.bed output.tsv
# ls clean3_pairs | sed -e "s/.pairs.gz//g" | parallel -j 30 --progress "cooler cload pairs -c1 2 -c2 4 -p1 3 -p2 5 /share/Data/public/ref_genome/mouse_ref/M23/raw_data/chr.len:1000000 clean3_pairs/{}.pairs.gz cools/{}.1m.cool"

import cooler
import numpy as np
import pandas as pd
import sys

clr = cooler.Cooler(sys.argv[1])
linearCpGData = pd.read_csv(sys.argv[2],sep="\t",header=None)

# fill in missing rows to 0 for linearCpGData
new_df = pd.DataFrame()
for chromosome, group in linearCpGData.groupby(0):
    group = group.sort_values(by=1)
    last_end = 0
    new_rows = []

    interval = group.iloc[1][1] - group.iloc[0][1] if len(group) > 1 else 0

    for index, row in group.iterrows():
        start = row[1]
        end = row[2]
        value = row[3]

        while last_end < start:
            new_rows.append([chromosome, last_end, last_end+interval, 0])
            last_end += interval

        new_rows.append([chromosome, start, end, value])
        last_end = end

    new_group = pd.DataFrame(new_rows)
    new_df = pd.concat([new_df, new_group])
linearCpGData = new_df
linearCpGData.columns=["chrom","start","end","cpg"]
chromlist = clr.chroms()[:][["name"]].values.T.tolist()[0]

res = []
for chrom in chromlist:
    matrix = clr.matrix(balance=False).fetch(chrom).astype("int")
    matrix = np.diag(np.ones(matrix.shape[0])) + matrix
    matrix[matrix>0] = 1
    linearCpGData_bychrom = linearCpGData.query("chrom == @chrom")
    if matrix.shape[0] != linearCpGData_bychrom.shape[0]:
        matrix = matrix[:linearCpGData_bychrom.shape[0],:linearCpGData_bychrom.shape[0]]
    linear_cpg_vector = linearCpGData_bychrom["cpg"].to_numpy()
    location = linearCpGData_bychrom[["chrom","start","end"]]
    location["scAB"] = np.dot(matrix,linearCpGData.query('chrom == @chrom')["cpg"].values) / (np.sum(matrix,axis=0)+1)
    res.append(location)
    
pd.concat(res).to_csv(sys.argv[3],sep="\t",header=None,index=None)
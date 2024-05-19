import cooler
import numpy as np
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor
import tqdm

def calc_scab_cell(clr,linearCpGData):
    """
    Calculate scA/B for single cell data
    clr: cooler object
    linearCpGData: pandas dataframe with CpG information
    """
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

    return pd.concat(res).set_index(["chrom","start","end"])



# example usage
# cellnames = metadata["cellname"].values
# linearCpGData = pd.read_csv("/share/Data/public/ref_genome/mouse_ref/M23/CpG/normal_cpg/M23.CpG.200000.txt",sep="\t",header=None)
# linearCpGData.columns = ["chrom","start","end","cpg"]

# linearCpGData = linearCpGData.query('chrom in @chroms')
# linearCpGData.columns = [0,1,2,3]

# def task_func(params):
#     cooler_obj, linearCpG = params
#     return calc_scab_cell(cooler_obj, linearCpG)

# with ProcessPoolExecutor(max_workers=40) as executor:
#     tasks = [(cooler.Cooler(f"/home/zliu/shareb/zliu/analysis/mouse_brain/Tan2020/cools/{cellname}.200k.cool"), linearCpGData) for cellname in cellnames]
#     res = list(tqdm.tqdm(executor.map(task_func, tasks), total=len(cellnames)))

# scAB = pd.concat(res, axis=1)
# scAB.columns = cellnames
# temp = scAB.index
# temp = [i[0] + "-" + str(i[1]) + "-" + str(i[2]) for i in temp]
# scAB.index = temp
# scAB = scAB.rank(axis=0, method="average", pct=True) 
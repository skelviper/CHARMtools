import pandas as pd
from pybedtools import BedTool
import numpy as np
import argparse
import matplotlib.pyplot as plt

def calculate_enrichment(tss_df, fragment_df, bin_size=100, flank=2000, chrom_length="/share/Data/public/ref_genome/human_ref/GRCh37d5/raw_data/hg19.chr.len",return_data = False):

    tss_bedtool = BedTool.from_dataframe(tss_df)
    fragment_bedtool = BedTool.from_dataframe(fragment_df)
    tss_flanking_bedtool = tss_bedtool.slop(b=(flank + bin_size/2), g=chrom_length)
    fragments_intersect_tss = fragment_bedtool.intersect(tss_flanking_bedtool, wa=True, wb=True)
    intersect_df = fragments_intersect_tss.to_dataframe(names=['chrom', 'start', 'end', 'tss_chrom', 'tss_start', 'tss_end', 'name', 'score', 'strand'])
    matrix_shape = (len(intersect_df), int(2 * flank + bin_size))
    binary_matrix = np.zeros(matrix_shape, dtype=int)

    for idx, row in intersect_df.iterrows():
        tss_start = row['tss_start']
        read_start = row['start']
        read_end = row['end']

        rel_start = read_start - tss_start
        rel_end = read_end - tss_start

        if rel_start < 0:
            rel_start = 0
        if rel_end > matrix_shape[1]:
            rel_end = matrix_shape[1]
        binary_matrix[idx, rel_start:rel_end] = 1

    summed_vector = binary_matrix.sum(axis=0)
    binned_vector = [summed_vector[i:i + bin_size] for i in range(0, len(summed_vector), bin_size)]

    mean_first_bin = binned_vector[0].mean()
    mean_middle_bin = binned_vector[len(binned_vector) // 2].mean()
    mean_last_bin = binned_vector[-1].mean()

    enrichment = mean_middle_bin / ((mean_first_bin + mean_last_bin) * 0.5)
    
    if return_data:
        return binary_matrix,enrichment
    else:
        return enrichment
    

def plot_enrichment(enrichment_vector, enrichment_score):
    # Set the figure size
    plt.figure(figsize=(4, 4))

    
    
    # Plot the enrichment vector
    plt.plot(enrichment_vector.sum(axis=0))

    # Add text to the plot
    plt.text(0.05, 0.9, "Enrichment score: %.3f" % enrichment_score, 
        horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes,
        fontsize=12)
    
    # Set labels and ticks
    plt.xlabel("Distance from TSS")
    plt.ylabel("Reads")
#    plt.ylim(top = 1.1*max(enrichment_vector.sum(axis=0)))

    vec_len = len(enrichment_vector.sum(axis=0))
    plt.xticks([0, vec_len/4, vec_len / 2, vec_len / 4 * 3, vec_len],
                ["-2k", "-1k", "0", "+1k", "+2k"])

    # Show the plot
    plt.show()

# Example usage:
# enrichment_vec, enrichment_score = tss_enrichment.calculate_enrichment(tss_df,human_fragments.iloc[:,0:3],return_data=True)
# plot_enrichment(enrichment_vec, enrichment_score)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate TSS enrichment')
    parser.add_argument('--tss', type=str, help='Path to TSS file')
    parser.add_argument('--fragments', type=str, help='Path to fragments file')
    parser.add_argument('--bin_size', type=int, default=100, help='Bin size for TSS enrichment calculation')
    parser.add_argument('--flank', type=int, default=2000, help='Flank size for TSS enrichment calculation')
    parser.add_argument('--chrom_length', type=str, default="/share/Data/public/ref_genome/human_ref/GRCh37d5/raw_data/hg19.chr.len", help='Path to chromosome length file')
    parser.add_argument('--output', type=str, default="tss_enrichment.txt", help='Output file name')
    args = parser.parse_args()

    tss_df = pd.read_csv(args.tss, sep="\t", header=None, names=['chrom', 'start', 'end'])
    tss_df = tss_df.iloc[:,0:3]

    fragment_df = pd.read_csv(args.fragments, sep="\t", header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    fragment_df = fragment_df.iloc[:,0:3]

    enrichment = calculate_enrichment(tss_df, fragment_df, bin_size=args.bin_size, flank=args.flank, chrom_length=args.chrom_length)
    print("Enrichment of {} on {} is {}".format(args.fragments, args.tss, enrichment))



# Usage in jupyter notebook
# import pandas as pd
# import multiprocessing as mp
# from tqdm import tqdm
# import warnings
# warnings.filterwarnings("ignore", category=UserWarning, message=".*inconsistent naming convention.*")

# def calculate_enrichment_wrapper(args):
#     cell_id, group, tss_df = args
#     group_subset = group.iloc[:,0:3]
#     _, enrichment_score = tss_enrichment.calculate_enrichment(tss_df, group_subset, return_data=True)
#     return [cell_id, len(group), enrichment_score]

# # Assuming you have the tss_enrichment module and tss_df loaded
# # and your human_fragments dataframe has columns 'cell_id', 'start', 'end'

# # Group human_fragments by cell_id
# grouped = human_fragments.groupby('cell_id')

# # Prepare arguments for parallel processing
# args_list = [(cell_id, group, tss_df) for cell_id, group in grouped]

# # Use a multiprocessing Pool to process the data in parallel
# with mp.Pool(processes=10) as pool:
#     results = list(tqdm(pool.imap_unordered(calculate_enrichment_wrapper, args_list), total=len(args_list)))

# # Convert the results list to a pandas DataFrame
# enrichment_df = pd.DataFrame(results, columns=['cell_id', 'fragments_count', 'enrichment_score'])

# # Display the enrichment_df
# #print(enrichment_df)

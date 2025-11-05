import pandas as pd
from scipy.spatial.distance import pdist, squareform
import numpy as np
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.distance import permdisp  # used for PERMDISP proxy
import numpy as np
from statsmodels.stats.multitest import multipletests
import os, argparse

# module add hpc-env/13.1
# module add Biopython/1.81-foss-2023a


def parse_arguments():
    parser = argparse.ArgumentParser(usage = "python3 permanova.py -i STR -g STR -o STR (-s | --no-s) -n STR [-h]\n\nPython 3 required!!\n\nThis program calculates Dxy as according to Irwin et al. (2016), Molecular Ecology 26(18), pp4488-4507. This between-group nucleotide differentiation (Dxy)  is determed as p1(1-p2) + p2(1-p1), where p1 is the grequecy of a given allele in the first population and p2 is the frequency of that allele in the second population and ranges from 0 to 1. This was originally made for both gppfst and minotaur. It also calculates pi for both populations and the difference between them. A .dxy output file is also printed to a desired directory.\n\nAuthor: Tane Kafle\n\n")
    parser.add_argument("-i", "--input_gtmat", required=True, help="input gtmat file pathway.")
    parser.add_argument("-g", "--input_groupings", required=True, help="input groupings file pathway.")
    parser.add_argument("-o", "--output_dir", required=True, help="output directory pathway. Please include final forward slash.")
    parser.add_argument("-s", "--spc_mode", help="True if comparing one species from different locations, False if comparing many species from a single location.", action=argparse.BooleanOptionalAction)
    parser.add_argument("-n", "--sample_name", required=True, help="Name of sample for output files.", type=str)

    return parser.parse_args()


def clean_indices(indices):
    """
    Cleans a list of indices by removing the duplicated part after an underscore.

    Args:
        indices (list): A list of strings, where some strings may have a
                        duplicated part separated by an underscore.

    Returns:
        list: A new list with the cleaned indices.
    """
    cleaned_indices = []
    for index in indices:
        # Check if the string contains an underscore
        if '_' in index:
            # Split the string at the first underscore
            parts = index.split('_', 1)
            # Check if the two parts are identical
            if parts[0] == parts[1]:
                cleaned_indices.append(parts[0])
            else:
                # If they're not identical, keep the original string
                cleaned_indices.append(index)
        else:
            # If there's no underscore, keep the original string
            cleaned_indices.append(index)
    return cleaned_indices

# python3 -m venv bioenv
# source bioenv/bin/activate
# pip install scikit-bio scipy pandas matplotlib seaborn
# python3 -m IPython

# Load the .traw file
#traw = pd.read_csv("/fs/dss/work/doau0129/chapter3/tane/allDATA/bySPC/pue_ld_pruned_gtmat.subsampled.traw", delim_whitespace=True)


def main():
    args = parse_arguments()

    os.makedirs(args.output_dir, exist_ok=True)

    traw = pd.read_csv(args.input_gtmat, sep=r"\s+")

    # Remove first 6 metadata columns (CHR, SNP, BP, A1, TEST, NMISS)

    snp_data = traw.iloc[:, 6:]
    # Transpose so samples are rows, SNPs are columns
    snp_matrix = snp_data.T
    snp_matrix = snp_matrix.fillna(0)  # or .dropna()


    # Optionally, assign sample IDs if available

    new_indices = clean_indices(snp_matrix.index)
    snp_matrix.index = new_indices

    # Load metadata
    #metadata = pd.read_csv("metadata.csv")

    # scp doau0129@hpcl001.hpc.uni-oldenburg.de:/fs/dss/work/doau0129/chapter3/tane/allDATA/bySPC/pue.txt .
    #metadata_df = pd.read_csv("/fs/dss/work/doau0129/chapter3/tane/allDATA/bySPC/pue.txt")
    metadata_df = pd.read_csv(args.input_groupings, header=None, names = ['sample'])

    #is_spc_mode_str = args.spc_mode
    #is_spc_mode = is_spc_mode_str.lower() in ('true', 't', '1')
    midfix=''
    print(args.spc_mode)
    if args.spc_mode:
        metadata_df["group"] = metadata_df["sample"].str[-3:]
        midfix='sm'
    else:
        print("going location mode")
        metadata_df["group"] = metadata_df["sample"].str[-6:-3]
        midfix='lm'

    #metadata_df = metadata_df.set_index("sample")


    # Optionally, assign sample IDs if available
    snp_matrix.index = metadata_df['sample'].to_list()

    #metadata = metadata.set_index("sample")

    #grouping = metadata["group"]


    # @to/do I need to remove all samples with sample size one.
    numsamples = metadata_df['group'].value_counts()


    groups2keep = numsamples[numsamples > 1].index.to_list()

    #metadata_filt_df = metadata_df[metadata_df['group'].isin(groups2keep)]
    #metadata_filt_df.index = metadata_filt_df['sample'].to_list()

    metadata_filt_df = metadata_df[metadata_df['group'].isin(groups2keep)].copy()



    # Ensure the sample order matches the SNP matrix
    #snp_matrix = snp_matrix.loc[metadata_filt_df['sample'].to_list()]
    samples2keep = metadata_filt_df['sample'].to_list()
    snp_matrix_filt = snp_matrix.loc[samples2keep]

    # Use Euclidean distance; you can choose others (e.g., 'cityblock', 'cosine')
    distance_array = pdist(snp_matrix_filt, metric='euclidean')
    distance_matrix = squareform(distance_array)


    dm = DistanceMatrix(distance_matrix, ids=snp_matrix.index)

    
    grouping = metadata_filt_df['group'].to_list()
    #permanova_result = permanova(dm, grouping, permutations=999)
    #permdisp_results = permdisp(dm, grouping, permutations=999)
    #print(permanova_result)
    #print(permdisp_results)

    #results = []
    #for i in [permanova_result, permdisp_results]:
    #    res_tmp = []
    #    for j in i:
    #        res_tmp.append(j)
    #    results.append(res_tmp)
        
    #results_df = pd.DataFrame(results, columns = ['test', 'test_statistic_name', 'sample_size', 'number_of_groups', 'test_statistic', 'p-value', 'number_of_permutations' ]) # why I didn't build df row by row: https://stackoverflow.com/questions/13784192/creating-an-empty-pandas-dataframe-and-then-filling-it


    #results_df.to_csv(os.path.join(args.output_dir, f'{args.sample_name}.csv') )

    #@to/do save this to file

    # i want a dataframe with:
    # group1 group2 num_in_group1 num_in_group2 test_statistic p_value


    groups = metadata_filt_df['group'].unique()

    #results = {}
    permanova_pairwise = []
    for i, g1 in enumerate(groups):
        for g2 in groups[i+1:]:
            # subset samples in these two groups
            print(f"comparing {g1} and {g2}")
            idx = metadata_filt_df[metadata_filt_df['group'].isin([g1, g2])].index
            subset_matrix = snp_matrix.iloc[idx]
            subset_metadata = metadata_filt_df.loc[idx]
            
            # compute distance matrix (e.g., Euclidean)
            dist_array = pdist(subset_matrix, metric='euclidean')
            dist_matrix = squareform(dist_array)
            
            # create DistanceMatrix object for scikit-bio
            dm = DistanceMatrix(dist_matrix, ids=subset_matrix.index)
            
            # run PERMANOVA
            result_pa = permanova(dm, subset_metadata['group'].to_list(), permutations=999)
            # run PERMdisp
            result_pd = permdisp(dm, subset_metadata['group'].to_list(), permutations=999)
        
            #results[(g1, g2)] = result
            permanova_pairwise.append([g1, g2, numsamples[g1], numsamples[g2], result_pa['test statistic'], result_pa['p-value'],  result_pd['test statistic'], result_pd['p-value']])
            


    permanova_df = pd.DataFrame(permanova_pairwise) # why I didn't build df row by row: https://stackoverflow.com/questions/13784192/creating-an-empty-pandas-dataframe-and-then-filling-it

    #permanova_df.rename(columns = ['spc1', 'spc2', 'n_spc1', 'n_spc2', 'permanova_teststat', 'permanova_pval', 'permadisp_teststat', 'permadisp_pval'])
    permanova_df.rename(columns = {0:'spc1', 1:'spc2', 2:'n_spc1', 3:'n_spc2', 4:'permanova_teststat', 
                                5: 'permanova_pval', 6: 'permadisp_teststat', 7: 'permadisp_pval'}, inplace=True)


    permanova_df['permanova_corr_pval'] = multipletests(permanova_df['permanova_pval'], alpha=0.05, method='fdr_bh')[1]
    permanova_df['permadisp_corr_pval'] = multipletests(permanova_df['permadisp_pval'], alpha=0.05, method='fdr_bh')[1]

    permanova_df.to_csv(os.path.join(args.output_dir, f'{args.sample_name}.{midfix}.pairwise.csv') )


    # # Extract distances per sample to group centroid
    # group_dists = []
    # for i, sample_id in enumerate(dm.ids):
    #     group = grouping[sample_id]
    #     dist = dm[i, :].mean()
    #     group_dists.append((sample_id, group, dist))

    # df_disp = pd.DataFrame(group_dists, columns=["sample", "group", "mean_distance"])

    # sns.boxplot(x="group", y="mean_distance", data=df_disp)
    # plt.title("PERMDISP-like Visualization")
    # plt.show()


if __name__ == "__main__":
    main()




# BASEDIR="/fs/dss/work/doau0129/chapter3/tane/allDATA/bySPC/"
# samplecode="pue"
# request_path = os.path.join(BASEDIR, f'{samplecode}.txt')
#'/Users/tkafle/Documents/hamlet_snp/analysis_files/permanova.py'

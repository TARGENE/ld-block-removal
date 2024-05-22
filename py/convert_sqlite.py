# Import packages 
import pandas as pd
import numpy as np
import sys
import dask.dataframe as dd
import matplotlib.pyplot as plt
import seaborn as sns
import sqlite3
from argparse import ArgumentParser

pd.set_option('display.max_columns', 500)

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--input", type = str, help="Input csv file with QTLs"
    )

    return parser.parse_args()

def get_LD_limits(rsid):
    rsid_filename = rsid.replace(":", "_")
    sqlite_file = rsid_filename+".sqlite"
    dat = sqlite3.connect(sqlite_file) # Connect to database
    query = dat.execute("SELECT * FROM LDView") # Select view of database. 
    cols = [column[0] for column in query.description] # Selcting columns from database
    qctools= pd.DataFrame.from_records(data = query.fetchall(), columns = cols) # Converting to pandas database
    if qctools.shape[0]==0:
        sys.stdout.write(f"No LD information available for {rsid}") 
        return pd.NA, pd.NA, pd.NA 
    else:
        # Check we are only looking at comparisons with the rsid of interest
        qctools= qctools[qctools["variant1_rsid"]==rsid] 
        # Check we are only considering every other SNP in LD block once
        qctools = qctools.drop_duplicates(subset=['variant2_rsid'],keep= 'last')
        qctools = qctools.drop_duplicates(subset=['variant2_rsid'],keep= 'last')
        qctools["distance"] = (qctools["variant2_position"]-qctools["variant1_position"].iloc[0])/1000
        outlier_down = qctools[qctools["dosage_r2"]>=0.1].sort_values(by="distance")["distance"].iloc[0]
        outlier_up = qctools[qctools["dosage_r2"]>=0.1].sort_values(by="distance",ascending=False)["distance"].iloc[0]   
        LDblock_length = outlier_up-outlier_down
        return outlier_down, outlier_up, LDblock_length

def plot_histogram(snps, outdir):
    # Plot the distribution of LD block lengths
    g, ax = plt.subplots(figsize=(14, 7))
    ax.set_facecolor('whitesmoke')
    sns.set_theme()
    sns.set_color_codes()
    counts, bins, bars = plt.hist(snps['LDblock_length'], color = "skyblue", lw=1)
    ax.set_title("LD block lengths")
    ax.set_ylabel("Frequency")
    ax.set_xlabel("LD block length (kb)")
    outname = outdir+'/LD_block_length_histogram.png'
    g.savefig(outname)

if __name__ == '__main__':
    input_file = get_input().input
    # Read in SNP file
    SNPs = pd.read_csv(input_file, sep = ',')

    # Apply function to all SNPs:
    SNPs[['LDblock_lower','LDblock_upper','LDblock_length']] = SNPs.apply(lambda row: get_LD_limits(row['RSID']), axis=1, result_type="expand")

    # Drop missing values
    SNPs = SNPs.dropna()
    # Have lower and upper bound
    SNPs["lower_bound"]=SNPs["POS"]+SNPs["LDblock_lower"]*1000
    SNPs["upper_bound"]=SNPs["POS"]+SNPs["LDblock_upper"]*1000
    SNPs.sort_values(by="LDblock_length")

    # Save to .csv file 
    basename = input_file.replace(".csv","")
    outname = basename+"_with_LD_blocks.csv"
    SNPs.to_csv(outname, index=False, header=True, sep=",")

    plot_histogram(SNPs, '.')
   

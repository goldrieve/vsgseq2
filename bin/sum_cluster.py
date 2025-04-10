import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def summarise_cluster_expression(tpm_file, cluster_summary_output, highest_tpm_vsg_output, fasta_file, output_fasta):
    """
    Calculates the sum of expression for each cluster and sample in a TPM file,
    identifies the vsg_id with the highest mean TPM value for each cluster,
    saves both results to separate CSV files, and extracts the corresponding
    sequences from a FASTA file, renaming them with the cluster number.

    Args:
        tpm_file (str): Path to the input TPM file (CSV format).
        cluster_summary_output (str): Path to the output CSV file for cluster summary.
        highest_tpm_vsg_output (str): Path to the output CSV file for highest TPM VSG IDs.
        fasta_file (str): Path to the input multi-FASTA file.
        output_fasta (str): Path to the output FASTA file containing extracted sequences.

    Returns:
        None
    """

    # Load the TPM file
    try:
        tpm_df = pd.read_csv(tpm_file)
    except FileNotFoundError:
        print(f"Error: TPM file not found at {tpm_file}")
        return
    except pd.errors.EmptyDataError:
        print(f"Error: TPM file is empty: {tpm_file}")
        return
    except pd.errors.ParserError:
        print(f"Error: Could not parse TPM file.  Check the format: {tpm_file}")
        return

    # Ensure 'cluster' and 'vsg_id' columns exist
    if 'cluster' not in tpm_df.columns:
        print("Error: 'cluster' column not found in the TPM file.")
        return
    if 'vsg_id' not in tpm_df.columns:
        print("Error: 'vsg_id' column not found in the TPM file.")
        return

    # Identify sample columns (all columns except 'vsg_id' and 'cluster')
    sample_cols = [col for col in tpm_df.columns if col not in ['vsg_id', 'cluster']]

    # Ensure there are sample columns
    if not sample_cols:
        print("Error: No sample columns found in the TPM file (expecting all columns except 'vsg_id' and 'cluster' to be samples).")
        return

    # Calculate the mean TPM value for each vsg_id across all samples
    try:
        tpm_df['mean_tpm'] = tpm_df[sample_cols].mean(axis=1)
    except KeyError as e:
        print(f"Error: Column not found during mean calculation: {e}")
        return
    except Exception as e:
        print(f"An unexpected error occurred during mean calculation: {e}")
        return

    # Group by 'cluster' and find the vsg_id with the highest mean TPM value
    try:
        highest_tpm_vsg = tpm_df.groupby('cluster').apply(lambda x: x.loc[x['mean_tpm'].idxmax(), 'vsg_id']).reset_index(name='highest_tpm_vsg_id')
    except Exception as e:
        print(f"An unexpected error occurred during grouping and max identification: {e}")
        return

    # Group by 'cluster' and sum the expression values for each sample
    try:
        cluster_summary = tpm_df.groupby('cluster')[sample_cols].sum()
    except KeyError as e:
        print(f"Error: Column not found during aggregation: {e}")
        return
    except Exception as e:
        print(f"An unexpected error occurred during aggregation: {e}")
        return

    # Save the cluster summary to a CSV file
    try:
        cluster_summary.to_csv(cluster_summary_output)
        print(f"Cluster summary saved to {cluster_summary_output}")
    except Exception as e:
        print(f"Error: Could not save cluster summary to CSV file: {e}")

    # Save the highest TPM VSG IDs to a CSV file
    try:
        highest_tpm_vsg.to_csv(highest_tpm_vsg_output, index=False)
        print(f"Highest mean TPM VSG IDs saved to {highest_tpm_vsg_output}")
    except Exception as e:
        print(f"Error: Could not save highest mean TPM VSG IDs to CSV file: {e}")

    # Extract sequences from FASTA file
    try:
        # Load the highest TPM VSG IDs
        highest_tpm_vsg_df = pd.read_csv(highest_tpm_vsg_output)
        highest_tpm_vsg_ids = highest_tpm_vsg_df['highest_tpm_vsg_id'].tolist()
        cluster_dict = pd.Series(highest_tpm_vsg_df.cluster.values,index=highest_tpm_vsg_df.highest_tpm_vsg_id).to_dict()

        # Load the FASTA file
        records = list(SeqIO.parse(fasta_file, "fasta"))

        # Extract matching sequences and rename them
        extracted_records = []
        for record in records:
            if record.id in highest_tpm_vsg_ids:
                cluster_number = cluster_dict[record.id]
                record.id = f"cluster_{cluster_number}"  # Changed the header
                record.description = ""  # Remove description
                extracted_records.append(record)

        # Write the extracted sequences to a new FASTA file
        SeqIO.write(extracted_records, output_fasta, "fasta")
        print(f"Extracted sequences saved to {output_fasta}")

    except FileNotFoundError:
        print(f"Error: FASTA file not found at {fasta_file}")
        return
    except Exception as e:
        print(f"Error: Could not extract sequences from FASTA file: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <tpm_file> <cluster_summary_output> <highest_tpm_vsg_output> <fasta_file> <output_fasta>")
        sys.exit(1)

    tpm_file = sys.argv[1]
    cluster_summary_output = sys.argv[2]
    highest_tpm_vsg_output = sys.argv[3]
    fasta_file = sys.argv[4]
    output_fasta = sys.argv[5]

    # Find the vsg_id with the highest mean TPM value for each cluster
    summarise_cluster_expression(tpm_file, cluster_summary_output, highest_tpm_vsg_output, fasta_file, output_fasta)
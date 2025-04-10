import pandas as pd
import argparse

def add_length_to_tpm(tpm_file, length_file, output_file):
    """
    Adds the length column from the VSG length count file to the TPM file based on vsg_id,
    placing the length column in the second position.

    Args:
        tpm_file (str): Path to the input TPM file.
        length_file (str): Path to the VSG length count file.
        output_file (str): Path to the output TPM file with lengths added.
    """
    # Load the TPM file and length file into DataFrames
    tpm_df = pd.read_csv(tpm_file)
    length_df = pd.read_csv(length_file)

    # Ensure the 'vsg_id' column exists in both files
    if 'vsg_id' not in tpm_df.columns:
        print("Error: 'vsg_id' column not found in the TPM file.")
        return
    if 'vsg_id' not in length_df.columns or 'length' not in length_df.columns:
        print("Error: 'vsg_id' or 'length' column not found in the length file.")
        return

    # Merge the TPM DataFrame with the length DataFrame on 'vsg_id'
    merged_df = pd.merge(tpm_df, length_df, on='vsg_id', how='left')

    # Get the list of columns
    cols = list(merged_df.columns)

    # Move 'sequence_length' to the second position
    if 'length' in cols:
        cols.remove('length')
        cols.insert(1, 'length')  # Insert at index 1 (second position)
        merged_df = merged_df[cols]

    # Save the updated DataFrame to the output file
    merged_df.to_csv(output_file, index=False)
    print(f"Updated TPM file with lengths saved to {output_file}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Add VSG sequence lengths to a TPM file based on vsg_id, placing length in the second column.")
    parser.add_argument("tpm_file", help="Path to the input TPM file.")
    parser.add_argument("length_file", help="Path to the VSG length count file.")
    parser.add_argument("output_file", help="Path to the output TPM file with lengths added.")

    # Parse command line arguments
    args = parser.parse_args()

    # Call the function with the provided arguments
    add_length_to_tpm(args.tpm_file, args.length_file, args.output_file)
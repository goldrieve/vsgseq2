import re
import pandas as pd

def add_cluster_to_lines(cluster_file, processed_cluster_file):
    """
    Process the cluster file to add cluster numbers to each line.
    """
    with open(cluster_file, 'r') as infile, open(processed_cluster_file, 'w') as outfile:
        current_cluster = None
        for line in infile:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = line.split()[1]  # Extract cluster number
                outfile.write(line + '\n')  # Write the cluster header
            elif line:
                match = re.search(r'>([^\s]+)', line)
                if match:
                    transcript_id = match.group(1)
                    outfile.write(f"{line}\t{current_cluster}\n")  # Add cluster to the line
                else:
                    outfile.write(line + '\n')  # Write other lines as is

def find_matching_transcripts(tpm_file, processed_cluster_file, output_file):
    """
    Search for every entry in the 'vsg_id' column of the tpm file
    in the processed cluster file, add a new column 'cluster' to the tpm file,
    and save the updated file with 'cluster' as the second column.
    """
    # Load the tpm file
    tpm_df = pd.read_csv(tpm_file)
    
    # Ensure the 'vsg_id' column exists in the tpm file
    if 'vsg_id' not in tpm_df.columns:
        print("Error: 'vsg_id' column not found in the tpm file.")
        return
    
    # Read the processed cluster file
    with open(processed_cluster_file, 'r') as f:
        cluster_lines = f.readlines()
    
    # Create a list to store cluster numbers
    cluster_column = []
    
    # Iterate through each vsg_id in the tpm file
    for vsg_id in tpm_df['vsg_id']:
        # Search for the vsg_id in the processed cluster file
        found = False
        for line in cluster_lines:
            if vsg_id in line:
                # Extract the cluster number (last column of the line)
                cluster_number = line.strip().split()[-1]
                cluster_column.append(cluster_number)
                found = True
                break
        if not found:
            cluster_column.append('NA')
    
    # Add the cluster column to the tpm DataFrame
    tpm_df['cluster'] = cluster_column
    
    # Reorder columns to make 'cluster' the second column
    columns = list(tpm_df.columns)
    reordered_columns = [columns[0], 'cluster'] + columns[1:-1]
    tpm_df = tpm_df[reordered_columns]
    
    # Save the updated DataFrame to a new file
    tpm_df.to_csv(output_file, index=False)
    print(f"Updated TPM file saved to {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script.py <tpm_file> <cluster_file> <output_file>")
        sys.exit(1)
    
    # Input and output file names
    tpm_file = sys.argv[1]
    cluster_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Temporary file to store the processed cluster file
    processed_cluster_file = "processed_cluster_file.txt"
    
    # Step 1: Process the cluster file to add cluster numbers
    add_cluster_to_lines(cluster_file, processed_cluster_file)
    
    # Step 2: Update the TPM file with the cluster column
    find_matching_transcripts(tpm_file, processed_cluster_file, output_file)
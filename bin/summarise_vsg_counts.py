import pandas as pd
import argparse

def calculate_cluster_stats(input_file, output_file):
    """
    Calculate max, mean, median, count of clusters with max, mean, and median > 100,
    and count of clusters with max, mean, and median > 10,000 for each cluster in the input CSV file.
    Filter out clusters with a max value below 100.

    Args:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to the output CSV file with calculated stats.
    """
    # Load the CSV file into a DataFrame
    df = pd.read_csv(input_file, sep=',')  # Adjust delimiter if necessary (e.g., sep='\t' for tab-separated files)

    # Ensure the 'cluster' column exists
    if 'cluster' not in df.columns:
        raise KeyError("The 'cluster' column is not found in the input file. Check the column names.")

    # Calculate max, mean, median for each cluster
    stats_df = df.set_index('cluster').agg(['max', 'mean', 'median'], axis=1)

    # Count clusters with max, mean, and median > 100
    clusters_above_100 = {
        'max_above_100': (stats_df['max'] > 100).sum(),
        'mean_above_100': (stats_df['mean'] > 100).sum(),
        'median_above_100': (stats_df['median'] > 100).sum()
    }

    # Count clusters with max, mean, and median > 10,000
    clusters_above_10000 = {
        'max_above_10000': (stats_df['max'] > 10000).sum(),
        'mean_above_10000': (stats_df['mean'] > 10000).sum(),
        'median_above_10000': (stats_df['median'] > 10000).sum()
    }

    # Print the counts    print(f"\nFiltered cluster statistics saved to ")
    print(str({output_file}) + "100:" + str(clusters_above_100))
    print(str({output_file}) + "1000:" + str(clusters_above_10000))

    # Filter out clusters with max < 100
    filtered_stats_df = stats_df[(stats_df['max'] >= 100)]

    # Save the filtered results to a new CSV file
    filtered_stats_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate max, mean, median, and count clusters with values > 100 and > 10,000 for each cluster in a CSV file. Filter out clusters with a max value below 100.")
    parser.add_argument("input_file", help="Path to the input CSV file.")
    parser.add_argument("output_file", help="Path to the output CSV file.")

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the function with the provided arguments
    calculate_cluster_stats(args.input_file, args.output_file)

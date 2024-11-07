import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compute Jaccard Similarity between gene sets from two CSV files.')
    parser.add_argument('-t', '--path-ground-truth', type=str, help='Path to the ground truth CSV file')
    parser.add_argument('-p', '--path-prediction', type=str, help='Path to the CSV file from the GRN tool')
    args = parser.parse_args()

    file_truth = args.path_ground_truth
    file_pred = args.path_prediction

    # Read CSV files
    df1 = pd.read_csv(file_truth, header=None)
    df2 = pd.read_csv(file_pred, header=None)

    # Clean and extract gene sets
    genes1 = set(df1[0].dropna().str.strip().str.upper())
    genes2 = set(df2[0].dropna().str.strip().str.upper())

    # Compute Jaccard Similarity
    intersection = genes1.intersection(genes2)
    union = genes1.union(genes2)
    jaccard_similarity = len(intersection) / len(union) if union else 0.0

    return jaccard_similarity

if __name__ == '__main__':
    main()

import requests
from bs4 import BeautifulSoup
import re
import pandas as pd
import argparse
import sys

def find_dataset_ids_for_gene(gene_name):
    """Find dataset IDs for a given gene name."""
    # Construct the URL with the specified gene name
    url = f"http://www.licpathway.net/KnockTFv2/search/search_tf_result.php?tf_name={gene_name}"
    
    # Send a request to the URL
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"Failed to retrieve the webpage: {e}")
        return []
    
    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Search for the Dataset ID pattern
    dataset_ids = set()
    for link in soup.find_all('a', href=True):
        matches = re.findall(r'DataSet_\d+_\d+', link['href'])
        for match in matches:
            dataset_ids.add(match)

    if dataset_ids:
        dataset_ids = list(dataset_ids)
        print("Dataset IDs found:", dataset_ids)
        return dataset_ids
    else:
        print("Dataset IDs not found.")
        return []

def fetch_and_process_data(dataset_id, gene_name, species, row_num, sort_meth):
    """Fetch and process data for a specific dataset ID."""
    params = {
        'species': species,
        'sample_tf_name': gene_name,
        'sample_id': dataset_id,
        'sel_row_num': row_num,  # Number of top genes
        'sort_meth': sort_meth   # Sorting method ('abs', 'up', 'down')
    }
    # Define the URL and parameters
    url = 'http://www.licpathway.net/KnockTFv2/search/sample_result_figure_rank/figure_rank.php'
    
    # Send the GET request
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
    except requests.RequestException as e:
        print(f"Error fetching data for dataset {dataset_id}: {e}")
        return pd.DataFrame()
    except ValueError as e:
        print(f"Error decoding JSON data for dataset {dataset_id}: {e}")
        return pd.DataFrame()

    # Process the data
    data_network = []

    # Extract the TF name
    try:
        tf_name = data[0]['figure_rank_50_name'][0]['name']
    except (KeyError, IndexError) as e:
        print(f"Error extracting TF name for dataset {dataset_id}: {e}")
        return pd.DataFrame()
    
    # Loop through the target genes
    for element in data[0]['figure_rank_50_name'][1:]:
        row = {
            'T(co)F': tf_name,
            'Target': element['name'],
            'log2FC': element['log2FC']
        }
        data_network.append(row)
    
    # Convert to DataFrame
    df = pd.DataFrame(data_network)
    return df

def main():
    """Main function to orchestrate data fetching and processing."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Fetch and process gene data.')
    parser.add_argument('-g', '--gene_name', required=True, help='Gene name')
    parser.add_argument('-n', '--row_num', type=int, default=100, help='Number of top genes to retrieve')
    parser.add_argument('-s', '--species', default='Mus musculus', choices=['Mus musculus', 'Homo sapiens'], help='Species name')
    parser.add_argument('-o', '--output', default=None, help='Output CSV file name')
    parser.add_argument('--sort_method', choices=['abs', 'up', 'down'], default='abs', help='Sorting method (abs, up, down)')
    args = parser.parse_args()
    
    gene_name = args.gene_name
    row_num = args.row_num
    species = args.species
    sort_meth = args.sort_method
    output_file = args.output or f'{gene_name}_data.csv'
    
    # Find dataset IDs
    dataset_ids = find_dataset_ids_for_gene(gene_name)
    if not dataset_ids:
        print("No dataset IDs found.")
        sys.exit(1)
    
    # Process each dataset ID
    combined_df = pd.DataFrame()
    for dataset_id in dataset_ids:
        print(f"Processing dataset {dataset_id}...")
        df = fetch_and_process_data(dataset_id, gene_name, species, row_num, sort_meth)
        if not df.empty:
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    if not combined_df.empty:
        # Write combined data to CSV
        combined_df.to_csv(output_file, index=False)
        print(f"Data successfully written to {output_file}")
    else:
        print("No data to write.")
        sys.exit(1)

if __name__ == '__main__':
    main()

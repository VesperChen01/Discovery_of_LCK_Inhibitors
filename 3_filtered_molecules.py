import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

def create_gradient_background(x, y):
    X, Y = np.meshgrid(x, y)
    Z = (X - x.min()) / (x.max() - x.min())  
    cmap = LinearSegmentedColormap.from_list('my_custom_cmap', [(1, 0, 0), (0, 1, 0)], N=256)
    return X, Y, Z, cmap

def plot_with_gradient(dataframe, x_column, y_column, dpi=300):
    x = dataframe[x_column]
    y = dataframe[y_column]
    x_range = np.linspace(x.min(), x.max(), 500)
    y_range = np.linspace(y.min(), y.max(), 500)
    
    X, Y, Z, cmap = create_gradient_background(x_range, y_range)
    
    plt.figure(figsize=(10, 6))
    plt.imshow(Z, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower', aspect='auto', cmap=cmap, alpha=0.5)
    sc = plt.scatter(x, y, color='blue', alpha=0.8)

    plt.axvline(x=-7, color='r', linestyle='--')
        
    
    plt.xlabel(f'{x_column}')
    plt.ylabel(y_column)
    plt.title('Glide_SP vs PLEANT')
   # plt.colorbar(sc, label=x_column)
   # plt.grid(True)
    plt.savefig('/data/crf/test/image.png', dpi=300)
    plt.close()
    
def filter_data(dataframe, pleant_threshold, Glide_SP):
    filtered_df = dataframe[(dataframe['PLEANT'] > pleant_threshold) & (dataframe['Glide_SP'] < Glide_SP)]
    return filtered_df[['CAS', 'PLEANT', 'Glide_SP']]  # Assuming 'CAS' is the molecule name column

# Main script
if __name__ == "__main__":
    # Load data
    data_file_path = '/data/crf/test/molecular_scores.csv'  # Replace with your data file path
    data = pd.read_csv(data_file_path)

    # Filter data
    pleant_threshold = 7
    Glide_SP = -7
    filtered_data = filter_data(data, pleant_threshold, Glide_SP)

    # Save filtered data to CSV
    filtered_data_file_path = 'filtered_molecules.csv'  # Replace with your desired output file path
    filtered_data.to_csv(filtered_data_file_path, index=False)
    print(f"Filtered data saved to {filtered_data_file_path}")

    # Plot data and save the figure
    plot_with_gradient(data, 'Glide_SP', 'PLEANT', plt.savefig)

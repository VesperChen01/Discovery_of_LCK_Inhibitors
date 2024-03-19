import matplotlib.pyplot as plt
import numpy as np

def draw_protein_donuts(data, protein_id_fontsize=12, percent_text_fontsize=10, filename="protein_donuts.png"):
    """
    Draws and saves a donut plot for the given protein data as a PNG file with 300 DPI.
    
    Parameters:
    data (dict): A dictionary with protein IDs as keys and a list of values representing
                 the distribution of amino acid residues in different regions.
    protein_id_fontsize (int): Font size for the protein IDs in the center of the donuts.
    percent_text_fontsize (int): Font size for the percentage texts on the donuts.
    filename (str): The name of the file to save the plot.
    """
    
    # Define labels and colors for the plot
    labels = ['Most favoured regions [A,B,L]', 'Additional allowed regions [a,b,l,p]', 
              'Generously allowed regions[~a~b,~l~p]', 'Disallowed regions']
    colors = ['#FF9999', '#66B2FF', '#99FF99', '#FF0000']
    
    # Create subplots for each protein structure in a 3x3 grid
    nrows, ncols = 3, 3
    fig, axs = plt.subplots(nrows, ncols, figsize=(15, 15))
    axs = axs.flatten()  # Flatten to 1D array for easy iteration

    # Define autopct to increase the font size of the percentage text
    def autopct_format(pct):
        return f"{pct:.1f}%" if pct >= 5 else ""

    # Plot each protein's data
    for i, (protein, values) in enumerate(data.items()):
        ax = axs[i]
        wedges, _, autotexts = ax.pie(values, labels=None, colors=colors,
                                      autopct=autopct_format, pctdistance=0.85, startangle=90)
        
        # Draw a circle at the center of pie to make it look like a donut
        centre_circle = plt.Circle((0, 0), 0.70, fc='white')
        ax.add_artist(centre_circle)

        # Change the font size of the autotexts (percentage labels)
        for autotext in autotexts:
            autotext.set_fontsize(percent_text_fontsize)

        # Add protein ID at the center of the donut with the specified font size
        ax.text(0, 0, protein, ha='center', va='center', fontsize=protein_id_fontsize)

        # Equal aspect ratio ensures that pie is drawn as a circle.
        ax.axis('equal')

    # Hide any unused subplots if there are fewer than 9 proteins
    for i in range(len(data), 9):
        axs[i].set_visible(False)

    # Adjust legend
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=label,
                                  markerfacecolor=color, markersize=10) for label, color in zip(labels, colors)]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.01))

    # Adjust overall spacing to prevent clipping of legend and subplots
    fig.subplots_adjust(bottom=0.2)

    # Save the figure as a PNG file with 300 DPI
    plt.savefig(filename, format='png', dpi=300)
    plt.close(fig)  # Close the figure to prevent it from displaying

    return filename

# Data for each protein
protein_data = {
    "1QP7": [87.9, 10.8, 0.7, 0.6],
    "2OF2": [93.3, 5.9, 0.8, 0.0],
    "2OFV": [89.9, 9.2, 0.9, 0.0],
    "3AC1": [93.3, 6.2, 0.4, 0.1],
    "3AD5": [93.8, 5.8, 0.0, 0.4],
    "3BYM": [93.3, 6.7, 0.0, 0.0],
    "3BYO": [91.2, 8.8, 0.0, 0.0],
    "3LCK": [94.6, 5.4, 0.0, 0.0],
    "6PDJ": [92.9, 6.2, 0.8, 0.1]
}

# Draw and save the donuts in a 3x3 grid
filename = draw_protein_donuts(protein_data, protein_id_fontsize=12, percent_text_fontsize=10)
print(f"The donut plot has been saved as: {filename}")


# Draw and save the donuts
draw_protein_donuts(protein_data, protein_id_fontsize=35, percent_text_fontsize=18, filename="protein_donuts.png")

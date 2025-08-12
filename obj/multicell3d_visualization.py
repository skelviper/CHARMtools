import matplotlib.pyplot as plt

class MultiCell3DVisualization:
    """
    Visualization methods for MultiCell3D class.
    """
    
    def plot_diff(self, diff, chrom_plot):
        """
        Plot the difference between two groups of cells.
        """
        df = diff.query('chrom == @chrom_plot')
        plt.figure(figsize=(6, 4))
        plt.plot(df['pos'], df['mean_group1'], label='Group 1', color='blue', alpha=1)
        plt.plot(df['pos'], df['mean_group2'], label='Group 2', color='red', alpha=1)

        significant_points = df[df['p_value_adj'] < 0.05]
        for index, row in significant_points.iterrows():
            plt.axvspan(row['pos'], row['pos'] + 1000000, color='grey', alpha=0.2)

        plt.ylim(0.7, 1.3)
        plt.xlabel(chrom_plot)
        plt.ylabel('Radial position')
        plt.legend()

        plt.show()
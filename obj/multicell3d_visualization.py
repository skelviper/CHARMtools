import matplotlib.pyplot as plt
import seaborn as sns

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
    
    def plot_radial_position(self, groupby, chrom=None, pos_range=None, figsize=(10, 6), 
                           alpha=0.7, s=20, title=None, save_path=None):
        """
        Plot radial position as scatter plot grouped by metadata column.
        
        Parameters:
        -----------
        groupby : str
            Column name in metadata to group cells by
        chrom : str, optional
            Specific chromosome to plot. If None, plot all chromosomes
        pos_range : tuple, optional
            Position range to plot (start, end). If None, plot all positions
        figsize : tuple, default (10, 6)
            Figure size
        alpha : float, default 0.7
            Transparency of points
        s : int, default 20
            Size of scatter points
        title : str, optional
            Plot title
        save_path : str, optional
            Path to save the figure
        
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure
        """
        # Check if radial_position matrix exists
        if "radial_position" not in self.matrices:
            raise ValueError("Radial position matrix not found. Please calculate it first.")
        
        # Get radial position data
        radial_df = self.matrices["radial_position"]
        
        # Check if groupby column exists in metadata
        if groupby not in self.metadata.columns:
            raise ValueError(f"Column '{groupby}' not found in metadata")
        
        # Filter by chromosome if specified
        if chrom is not None:
            # Extract chromosome from index (assuming format like 'chr1:1000000')
            chrom_mask = radial_df.index.str.startswith(f"{chrom}:")
            radial_df = radial_df.loc[chrom_mask]
        
        # Filter by position range if specified
        if pos_range is not None:
            start_pos, end_pos = pos_range
            # Extract position from index
            positions = radial_df.index.str.split(':').str[1].astype(int)
            pos_mask = (positions >= start_pos) & (positions <= end_pos)
            radial_df = radial_df.loc[pos_mask]
        
        if radial_df.empty:
            raise ValueError("No data available for the specified filters")
        
        # Group cells by metadata column and calculate mean
        grouped_data = []
        
        for group_name, group_cells in self.metadata.groupby(groupby):
            # Get cells that exist in both metadata and radial_df columns
            common_cells = list(set(group_cells.index) & set(radial_df.columns))
            
            if len(common_cells) == 0:
                continue
                
            # Calculate mean radial position for this group
            group_mean = radial_df[common_cells].mean(axis=1)
            
            # Extract chromosome and position for plotting
            for idx, mean_val in group_mean.items():
                if pd.notna(mean_val):
                    chrom_pos = idx.split(':')
                    if len(chrom_pos) == 2:
                        chrom_name, pos = chrom_pos
                        grouped_data.append({
                            'chromosome': chrom_name,
                            'position': int(pos),
                            'radial_position': mean_val,
                            'group': str(group_name)
                        })
        
        if not grouped_data:
            raise ValueError("No valid data found for plotting")
        
        # Convert to DataFrame
        plot_df = pd.DataFrame(grouped_data)
        
        # Create the plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get unique groups and assign colors
        unique_groups = plot_df['group'].unique()
        colors = sns.color_palette("husl", len(unique_groups))
        
        # Plot each group
        for i, group in enumerate(unique_groups):
            group_data = plot_df[plot_df['group'] == group]
            ax.scatter(group_data['position'], group_data['radial_position'], 
                      c=[colors[i]], label=group, alpha=alpha, s=s)
        
        # Customize plot
        ax.set_xlabel('Genomic Position')
        ax.set_ylabel('Radial Position')
        
        if title is None:
            if chrom is not None:
                title = f'Radial Position by {groupby} - {chrom}'
            else:
                title = f'Radial Position by {groupby}'
        ax.set_title(title)
        
        # Add legend
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Format x-axis for better readability
        ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
        
        # Adjust layout
        plt.tight_layout()
        
        # Save if path provided
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
# Cell3D Visualization Module - 3D plotting and visualization
import warnings
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

class Cell3DVisualization:
    """Visualization methods for Cell3D objects"""
    
    def plot3D(self, feature=None, genome_coord=None, smooth=False, smooth_sigma=1.0, 
               color_map='viridis', point_size=50, alpha=0.7, figsize=(10, 8), 
               title=None, save_path=None, show_axes=True, background_color='white'):
        """Create 3D scatter plot of the chromatin structure"""
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
        except ImportError:
            raise ImportError("matplotlib is required for plotting")
        
        if self.on_disk:
            self.to_memory()
        
        # Get data to plot
        if genome_coord is None:
            plot_data = self.tdg.copy()
        else:
            plot_data = self.tdg.query(f"chrom == '{genome_coord}'").copy()
        
        if len(plot_data) == 0:
            raise ValueError("No data to plot")
        
        # Sort by chromosome and position for consistent plotting
        plot_data = plot_data.sort_values(['chrom', 'pos'])
        
        # Extract coordinates
        x = plot_data['x'].values
        y = plot_data['y'].values
        z = plot_data['z'].values
        
        # Determine colors
        if feature is None:
            # Color by chromosome
            unique_chroms = plot_data['chrom'].unique()
            chrom_colors = {chrom: i for i, chrom in enumerate(unique_chroms)}
            colors = [chrom_colors[chrom] for chrom in plot_data['chrom']]
            color_label = 'Chromosome'
        else:
            if feature not in plot_data.columns:
                raise ValueError(f"Feature '{feature}' not found in data")
            colors = plot_data[feature].values
            color_label = feature
        
        # Apply smoothing if requested
        if smooth and feature is not None:
            # Create a grid for interpolation
            xi = np.linspace(x.min(), x.max(), 50)
            yi = np.linspace(y.min(), y.max(), 50)
            zi = np.linspace(z.min(), z.max(), 50)
            
            # Interpolate colors to grid
            points = np.column_stack((x, y, z))
            xi_grid, yi_grid, zi_grid = np.meshgrid(xi, yi, zi, indexing='ij')
            grid_points = np.column_stack((xi_grid.ravel(), yi_grid.ravel(), zi_grid.ravel()))
            
            try:
                interpolated_colors = griddata(points, colors, grid_points, method='linear')
                interpolated_colors = interpolated_colors.reshape(xi_grid.shape)
                
                # Apply Gaussian smoothing
                interpolated_colors = gaussian_filter(interpolated_colors, sigma=smooth_sigma)
                
                # Map back to original points
                colors = griddata(grid_points, interpolated_colors.ravel(), points, method='nearest')
            except Exception as e:
                warnings.warn(f"Smoothing failed: {e}. Using original colors.")
        
        # Create the plot
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        # Set background color
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        
        if background_color != 'white':
            ax.xaxis.pane.set_edgecolor(background_color)
            ax.yaxis.pane.set_edgecolor(background_color)
            ax.zaxis.pane.set_edgecolor(background_color)
        
        # Create scatter plot
        scatter = ax.scatter(x, y, z, c=colors, cmap=color_map, s=point_size, alpha=alpha)
        
        # Add colorbar
        if feature is not None:
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=20)
            cbar.set_label(color_label)
        
        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        if title is None:
            if genome_coord:
                title = f'3D Structure - {genome_coord}'
            else:
                title = '3D Chromatin Structure'
            if feature:
                title += f' (colored by {feature})'
        
        ax.set_title(title)
        
        # Hide axes if requested
        if not show_axes:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_zlabel('')
        
        # Save or show
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_distance_matrix(self, genome_coord=None, method='euclidean', figsize=(8, 6), 
                           color_map='viridis', title=None, save_path=None):
        """Plot distance matrix heatmap"""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            raise ImportError("matplotlib and seaborn are required for plotting")
        
        if self.on_disk:
            self.to_memory()
        
        # Get data
        if genome_coord is None:
            plot_data = self.tdg.copy()
        else:
            plot_data = self.tdg.query(f"chrom == '{genome_coord}'").copy()
        
        if len(plot_data) == 0:
            raise ValueError("No data to plot")
        
        # Calculate distance matrix
        coords = plot_data[['x', 'y', 'z']].values
        from scipy.spatial.distance import pdist, squareform
        
        if method == 'euclidean':
            distances = squareform(pdist(coords, metric='euclidean'))
        else:
            distances = squareform(pdist(coords, metric=method))
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot heatmap
        im = ax.imshow(distances, cmap=color_map, aspect='auto')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Distance')
        
        # Set title
        if title is None:
            if genome_coord:
                title = f'Distance Matrix - {genome_coord}'
            else:
                title = 'Distance Matrix'
        
        ax.set_title(title)
        ax.set_xlabel('Genomic Position Index')
        ax.set_ylabel('Genomic Position Index')
        
        # Save or show
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_feature_distribution(self, feature_name, bins=50, figsize=(8, 6), 
                                color='skyblue', title=None, save_path=None):
        """Plot distribution of a feature"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")
        
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found in data")
        
        # Get feature values
        feature_values = self.tdg[feature_name].values
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot histogram
        ax.hist(feature_values, bins=bins, color=color, alpha=0.7, edgecolor='black')
        
        # Set labels and title
        ax.set_xlabel(feature_name)
        ax.set_ylabel('Frequency')
        
        if title is None:
            title = f'Distribution of {feature_name}'
        
        ax.set_title(title)
        
        # Add statistics
        mean_val = np.mean(feature_values)
        std_val = np.std(feature_values)
        ax.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
        ax.axvline(mean_val + std_val, color='orange', linestyle='--', alpha=0.7, label=f'+1 SD: {mean_val + std_val:.2f}')
        ax.axvline(mean_val - std_val, color='orange', linestyle='--', alpha=0.7, label=f'-1 SD: {mean_val - std_val:.2f}')
        
        ax.legend()
        
        # Save or show
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_radial_distribution(self, feature_name, n_bins=10, figsize=(8, 6), 
                               color='green', title=None, save_path=None):
        """Plot radial distribution of a feature"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")
        
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found in data")
        
        # Calculate radial distribution
        distribution, bin_edges = self.feature_radial_distribution(feature_name, n_bins=n_bins)
        
        # Calculate bin centers
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot line
        ax.plot(bin_centers, distribution, marker='o', color=color, linewidth=2, markersize=6)
        
        # Set labels and title
        ax.set_xlabel('Radial Distance from Center')
        ax.set_ylabel(f'Mean {feature_name}')
        
        if title is None:
            title = f'Radial Distribution of {feature_name}'
        
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        
        # Save or show
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.tight_layout()
        return fig, ax
    
    def plot_chromosome_comparison(self, feature_name, figsize=(12, 8), 
                                 color_map='Set3', title=None, save_path=None):
        """Plot feature comparison across chromosomes"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")
        
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found in data")
        
        # Group by chromosome
        chrom_data = self.tdg.groupby('chrom')[feature_name].agg(['mean', 'std', 'count']).reset_index()
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot bars with error bars
        x_pos = np.arange(len(chrom_data))
        bars = ax.bar(x_pos, chrom_data['mean'], yerr=chrom_data['std'], 
                     capsize=5, alpha=0.7, color=plt.cm.get_cmap(color_map)(np.linspace(0, 1, len(chrom_data))))
        
        # Set labels and title
        ax.set_xlabel('Chromosome')
        ax.set_ylabel(f'Mean {feature_name}')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(chrom_data['chrom'], rotation=45)
        
        if title is None:
            title = f'{feature_name} by Chromosome'
        
        ax.set_title(title)
        
        # Add value labels on bars
        for i, (bar, mean_val, count) in enumerate(zip(bars, chrom_data['mean'], chrom_data['count'])):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + chrom_data['std'].iloc[i],
                   f'{mean_val:.2f}\n(n={count})', ha='center', va='bottom', fontsize=8)
        
        # Save or show
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.tight_layout()
        return fig, ax
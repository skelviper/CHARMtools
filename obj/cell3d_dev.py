# Cell3D Development Module - Development and experimental features
import warnings
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.optimize import minimize

class Cell3DDev:
    """Development and experimental features for Cell3D objects"""
    
    def _fdg_read_pairs(self, pairs_file, max_distance=2000000):
        """Read pairs file for FDG algorithm"""
        try:
            pairs_df = pd.read_csv(pairs_file, sep='\t', header=None)
            
            # Assume standard pairs format: chr1, pos1, chr2, pos2, [additional columns]
            if pairs_df.shape[1] >= 4:
                pairs_df.columns = ['chr1', 'pos1', 'chr2', 'pos2'] + [f'col_{i}' for i in range(4, pairs_df.shape[1])]
            else:
                raise ValueError("Pairs file must have at least 4 columns")
            
            # Filter by distance if on same chromosome
            if max_distance is not None:
                same_chr_mask = pairs_df['chr1'] == pairs_df['chr2']
                distance_mask = np.abs(pairs_df['pos1'] - pairs_df['pos2']) <= max_distance
                pairs_df = pairs_df[~same_chr_mask | distance_mask]
            
            return pairs_df
            
        except Exception as e:
            raise ValueError(f"Error reading pairs file: {str(e)}")
    
    def _calc_forces(self, coords, pairs_data, k_spring=1.0, k_repulsion=1.0, 
                    equilibrium_distance=1.0, repulsion_distance=2.0):
        """Calculate forces for FDG algorithm"""
        n_points = len(coords)
        forces = np.zeros_like(coords)
        
        # Spring forces from pairs
        for _, pair in pairs_data.iterrows():
            try:
                # Find indices for the pair positions
                idx1 = self._position_to_index(pair['chr1'], pair['pos1'])
                idx2 = self._position_to_index(pair['chr2'], pair['pos2'])
                
                if idx1 is not None and idx2 is not None and idx1 != idx2:
                    # Calculate spring force
                    diff = coords[idx2] - coords[idx1]
                    distance = np.linalg.norm(diff)
                    
                    if distance > 0:
                        direction = diff / distance
                        force_magnitude = k_spring * (distance - equilibrium_distance)
                        force = force_magnitude * direction
                        
                        forces[idx1] += force
                        forces[idx2] -= force
            except:
                continue
        
        # Repulsion forces between all points
        for i in range(n_points):
            for j in range(i + 1, n_points):
                diff = coords[j] - coords[i]
                distance = np.linalg.norm(diff)
                
                if distance > 0 and distance < repulsion_distance:
                    direction = diff / distance
                    force_magnitude = k_repulsion / (distance**2 + 1e-10)
                    force = force_magnitude * direction
                    
                    forces[i] -= force
                    forces[j] += force
        
        return forces
    
    def _position_to_index(self, chrom, pos):
        """Convert genomic position to index in tdg dataframe"""
        if self.on_disk:
            self.to_memory()
        
        try:
            mask = (self.tdg['chrom'] == chrom) & (self.tdg['pos'] == pos)
            indices = self.tdg.index[mask].tolist()
            return indices[0] if indices else None
        except:
            return None
    
    def _fdg(self, pairs_data, n_iterations=1000, learning_rate=0.01, 
            k_spring=1.0, k_repulsion=1.0, convergence_threshold=1e-6):
        """Force-directed graph layout algorithm"""
        if self.on_disk:
            self.to_memory()
        
        # Initialize coordinates randomly if not present
        if 'x' not in self.tdg.columns or 'y' not in self.tdg.columns or 'z' not in self.tdg.columns:
            n_points = len(self.tdg)
            self.tdg['x'] = np.random.randn(n_points)
            self.tdg['y'] = np.random.randn(n_points)
            self.tdg['z'] = np.random.randn(n_points)
        
        coords = self.tdg[['x', 'y', 'z']].values.copy()
        
        prev_energy = float('inf')
        
        for iteration in range(n_iterations):
            # Calculate forces
            forces = self._calc_forces(coords, pairs_data, k_spring, k_repulsion)
            
            # Update coordinates
            coords += learning_rate * forces
            
            # Calculate energy for convergence check
            if iteration % 100 == 0:
                energy = np.sum(forces**2)
                
                if abs(prev_energy - energy) < convergence_threshold:
                    print(f"Converged at iteration {iteration}")
                    break
                
                prev_energy = energy
                print(f"Iteration {iteration}, Energy: {energy:.6f}")
        
        # Update coordinates in dataframe
        self.tdg['x'] = coords[:, 0]
        self.tdg['y'] = coords[:, 1]
        self.tdg['z'] = coords[:, 2]
        
        return coords
    
    def _fdg_from_pairs(self, pairs_file, max_distance=2000000, n_iterations=1000, 
                       learning_rate=0.01, k_spring=1.0, k_repulsion=1.0):
        """Run FDG algorithm from pairs file"""
        # Read pairs data
        pairs_data = self._fdg_read_pairs(pairs_file, max_distance)
        
        # Run FDG algorithm
        final_coords = self._fdg(pairs_data, n_iterations, learning_rate, k_spring, k_repulsion)
        
        return final_coords
    
    def optimize_structure(self, method='gradient_descent', **kwargs):
        """Optimize 3D structure using various methods"""
        if self.on_disk:
            self.to_memory()
        
        coords = self.tdg[['x', 'y', 'z']].values
        
        if method == 'gradient_descent':
            optimized_coords = self._gradient_descent_optimization(coords, **kwargs)
        elif method == 'simulated_annealing':
            optimized_coords = self._simulated_annealing_optimization(coords, **kwargs)
        elif method == 'genetic_algorithm':
            optimized_coords = self._genetic_algorithm_optimization(coords, **kwargs)
        else:
            raise ValueError(f"Unsupported optimization method: {method}")
        
        # Update coordinates
        self.tdg['x'] = optimized_coords[:, 0]
        self.tdg['y'] = optimized_coords[:, 1]
        self.tdg['z'] = optimized_coords[:, 2]
        
        return optimized_coords
    
    def _gradient_descent_optimization(self, coords, learning_rate=0.01, n_iterations=1000, **kwargs):
        """Gradient descent optimization"""
        current_coords = coords.copy()
        
        for i in range(n_iterations):
            # Calculate gradient (simplified)
            gradient = self._calculate_gradient(current_coords)
            
            # Update coordinates
            current_coords -= learning_rate * gradient
            
            if i % 100 == 0:
                energy = self._calculate_energy(current_coords)
                print(f"Iteration {i}, Energy: {energy:.6f}")
        
        return current_coords
    
    def _simulated_annealing_optimization(self, coords, initial_temp=1000, cooling_rate=0.95, 
                                        min_temp=1, n_iterations=1000, **kwargs):
        """Simulated annealing optimization"""
        current_coords = coords.copy()
        best_coords = coords.copy()
        
        current_energy = self._calculate_energy(current_coords)
        best_energy = current_energy
        
        temperature = initial_temp
        
        for i in range(n_iterations):
            # Generate neighbor solution
            neighbor_coords = current_coords + np.random.normal(0, 0.1, current_coords.shape)
            neighbor_energy = self._calculate_energy(neighbor_coords)
            
            # Accept or reject
            delta_energy = neighbor_energy - current_energy
            
            if delta_energy < 0 or np.random.random() < np.exp(-delta_energy / temperature):
                current_coords = neighbor_coords
                current_energy = neighbor_energy
                
                if current_energy < best_energy:
                    best_coords = current_coords.copy()
                    best_energy = current_energy
            
            # Cool down
            temperature = max(min_temp, temperature * cooling_rate)
            
            if i % 100 == 0:
                print(f"Iteration {i}, Temperature: {temperature:.2f}, Best Energy: {best_energy:.6f}")
        
        return best_coords
    
    def _genetic_algorithm_optimization(self, coords, population_size=50, n_generations=100, 
                                      mutation_rate=0.1, crossover_rate=0.8, **kwargs):
        """Genetic algorithm optimization"""
        # Initialize population
        population = []
        for _ in range(population_size):
            individual = coords + np.random.normal(0, 0.5, coords.shape)
            population.append(individual)
        
        for generation in range(n_generations):
            # Evaluate fitness
            fitness_scores = [1.0 / (self._calculate_energy(individual) + 1e-10) for individual in population]
            
            # Selection
            new_population = []
            for _ in range(population_size):
                # Tournament selection
                tournament_indices = np.random.choice(population_size, 3, replace=False)
                tournament_fitness = [fitness_scores[i] for i in tournament_indices]
                winner_idx = tournament_indices[np.argmax(tournament_fitness)]
                new_population.append(population[winner_idx].copy())
            
            # Crossover and mutation
            for i in range(0, population_size - 1, 2):
                if np.random.random() < crossover_rate:
                    # Single-point crossover
                    crossover_point = np.random.randint(1, len(coords))
                    
                    child1 = new_population[i].copy()
                    child2 = new_population[i + 1].copy()
                    
                    child1[crossover_point:] = new_population[i + 1][crossover_point:]
                    child2[crossover_point:] = new_population[i][crossover_point:]
                    
                    new_population[i] = child1
                    new_population[i + 1] = child2
                
                # Mutation
                if np.random.random() < mutation_rate:
                    mutation_strength = 0.1
                    new_population[i] += np.random.normal(0, mutation_strength, coords.shape)
                
                if np.random.random() < mutation_rate:
                    mutation_strength = 0.1
                    new_population[i + 1] += np.random.normal(0, mutation_strength, coords.shape)
            
            population = new_population
            
            if generation % 10 == 0:
                best_fitness = max(fitness_scores)
                print(f"Generation {generation}, Best Fitness: {best_fitness:.6f}")
        
        # Return best individual
        final_fitness = [1.0 / (self._calculate_energy(individual) + 1e-10) for individual in population]
        best_idx = np.argmax(final_fitness)
        
        return population[best_idx]
    
    def _calculate_energy(self, coords):
        """Calculate energy of the current configuration"""
        # Simple energy function based on pairwise distances
        n_points = len(coords)
        energy = 0.0
        
        for i in range(n_points):
            for j in range(i + 1, n_points):
                distance = np.linalg.norm(coords[i] - coords[j])
                # Lennard-Jones-like potential
                energy += 1.0 / (distance**12 + 1e-10) - 1.0 / (distance**6 + 1e-10)
        
        return energy
    
    def _calculate_gradient(self, coords):
        """Calculate gradient of energy function"""
        n_points = len(coords)
        gradient = np.zeros_like(coords)
        
        for i in range(n_points):
            for j in range(n_points):
                if i != j:
                    diff = coords[i] - coords[j]
                    distance = np.linalg.norm(diff)
                    
                    if distance > 0:
                        direction = diff / distance
                        # Derivative of Lennard-Jones-like potential
                        force_magnitude = -12.0 / (distance**13 + 1e-10) + 6.0 / (distance**7 + 1e-10)
                        gradient[i] += force_magnitude * direction
        
        return gradient
    
    def experimental_feature_analysis(self, feature_name, method='pca'):
        """Experimental feature analysis methods"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found")
        
        feature_values = self.tdg[feature_name].values
        coords = self.tdg[['x', 'y', 'z']].values
        
        if method == 'pca':
            from sklearn.decomposition import PCA
            
            # Combine coordinates and feature for PCA
            combined_data = np.column_stack([coords, feature_values])
            
            pca = PCA(n_components=3)
            transformed_data = pca.fit_transform(combined_data)
            
            return {
                'transformed_coordinates': transformed_data,
                'explained_variance_ratio': pca.explained_variance_ratio_,
                'components': pca.components_
            }
        
        elif method == 'tsne':
            from sklearn.manifold import TSNE
            
            # Use t-SNE for dimensionality reduction
            combined_data = np.column_stack([coords, feature_values])
            
            tsne = TSNE(n_components=3, random_state=42)
            transformed_data = tsne.fit_transform(combined_data)
            
            return {
                'transformed_coordinates': transformed_data,
                'kl_divergence': tsne.kl_divergence_
            }
        
        else:
            raise ValueError(f"Unsupported analysis method: {method}")
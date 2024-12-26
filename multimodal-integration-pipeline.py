import os
import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse

class MultiModalIntegration:
    def __init__(self, data_dir, rna_file="rna.h5ad", atac_file="atac.h5ad"):
        """
        Initialize the pipeline with data directory and file names
        """
        self.data_dir = data_dir
        self.rna_file = rna_file
        self.atac_file = atac_file
        self.rna = None
        self.atac = None
        self.combined = None
        self.model = None
        
    def load_data(self):
        """
        Load RNA and ATAC data from h5ad files
        """
        try:
            os.chdir(self.data_dir)
            self.rna = sc.read_h5ad(self.rna_file)
            self.atac = sc.read_h5ad(self.atac_file)
            
            # Make variable names unique
            self.rna.var_names_make_unique()
            self.atac.var_names_make_unique()
            
            print(f"Loaded RNA data: {self.rna.shape}")
            print(f"Loaded ATAC data: {self.atac.shape}")
            
        except Exception as e:
            raise Exception(f"Error loading data: {str(e)}")
            
    def prepare_rna_data(self):
        """
        Prepare RNA data for integration
        """
        rna_var = pd.DataFrame({
            'ID': self.rna.var_names.tolist(),
            'modality': ['RNA'] * len(self.rna.var_names),
            'feature_types': ['Gene'] * len(self.rna.var_names)
        })
        
        rna_obs = self.rna.obs[['Mutation']].copy()
        
        self.rna_adata = sc.AnnData(
            X=self.rna.X.copy(),
            obs=rna_obs,
            var=rna_var
        )
        
    def prepare_atac_data(self):
        """
        Prepare ATAC data for integration
        """
        atac_var = pd.DataFrame({
            'ID': self.atac.var_names.tolist(),
            'modality': ['ATAC'] * len(self.atac.var_names),
            'feature_types': ['Peak'] * len(self.atac.var_names),
            'chr': self.atac.var_names.str.split('-').str[0],
            'start': self.atac.var_names.str.split('-').str[1],
            'end': self.atac.var_names.str.split('-').str[2]
        })
        
        atac_obs = self.atac.obs[['Mutation']].copy()
        
        self.atac_adata = sc.AnnData(
            X=self.atac.X.copy(),
            obs=atac_obs,
            var=atac_var
        )
        
    def combine_modalities(self):
        """
        Combine RNA and ATAC data
        """
        # Ensure obs indices are strings
        self.rna_adata.obs_names = self.rna_adata.obs_names.astype(str)
        self.atac_adata.obs_names = self.atac_adata.obs_names.astype(str)
        
        # Concatenate
        self.combined = sc.concat(
            [self.rna_adata, self.atac_adata],
            axis=1,
            join='outer',
            merge='unique',
            label='modality',
            keys=['RNA', 'ATAC']
        )
        
        # Add modality information
        self.combined.obs['modality'] = pd.Categorical(
            ['RNA' if idx in self.rna_adata.obs_names else 'ATAC' 
             for idx in self.combined.obs_names]
        )
        
        # Add metadata
        self.combined.uns['modality_counts'] = {
            'RNA': sum(self.combined.var['modality'] == 'RNA'),
            'ATAC': sum(self.combined.var['modality'] == 'ATAC')
        }
        
        # Filter features
        sc.pp.filter_genes(self.combined, min_cells=int(self.combined.shape[0] * 0.01))
        print(f"Combined data shape after filtering: {self.combined.shape}")
        
    def train_multivi(self, save_dir="multivi_model"):
        """
        Train MultiVI model
        """
        scvi.model.MULTIVI.setup_anndata(
            self.combined,
            batch_key="modality"
        )
        
        self.model = scvi.model.MULTIVI(
            self.combined,
            n_genes=(self.combined.var["modality"] == "RNA").sum(),
            n_regions=(self.combined.var["modality"] == "ATAC").sum(),
        )
        
        self.model.train()
        
        # Save model
        model_dir = os.path.join(self.data_dir, save_dir)
        self.model.save(model_dir, overwrite=True)
        
    def process_latent_space(self, min_dist=0.2):
        """
        Process and visualize latent space
        """
        MULTIVI_LATENT_KEY = "X_multivi"
        self.combined.obsm[MULTIVI_LATENT_KEY] = self.model.get_latent_representation()
        
        # Compute neighbors and UMAP
        sc.pp.neighbors(self.combined, use_rep=MULTIVI_LATENT_KEY)
        sc.tl.umap(self.combined, min_dist=min_dist)
        
    def visualize_results(self, imputed_expression=None):
        """
        Create visualizations of the results
        """
        # UMAP plot
        sc.pl.umap(self.combined, color="modality")
        
        if imputed_expression is not None:
            fig, axes = plt.subplots(2, 2, figsize=(15, 15))
            
            # Expression distribution
            sns.histplot(np.mean(imputed_expression.values, axis=0), bins=50, ax=axes[0,0])
            axes[0,0].set_title('Mean Expression Distribution')
            axes[0,0].set_xlabel('Mean Expression')
            
            # Feature variance
            feature_var = np.var(imputed_expression.values, axis=0)
            sns.histplot(feature_var, bins=50, ax=axes[0,1])
            axes[0,1].set_title('Feature Variance Distribution')
            axes[0,1].set_xlabel('Variance')
            
            # Heatmap
            n_cells, n_features = 100, 50
            top_var_idx = np.argsort(feature_var)[-n_features:]
            random_cell_idx = np.random.choice(imputed_expression.shape[0], n_cells, replace=False)
            heatmap_data = imputed_expression.values[random_cell_idx][:, top_var_idx]
            
            sns.heatmap(heatmap_data, 
                       yticklabels=False, 
                       xticklabels=False, 
                       cmap='viridis', 
                       ax=axes[1,0])
            axes[1,0].set_title('Expression Heatmap\n(Random cells, top variable features)')
            
            # Correlation plot
            n_features_corr = 20
            top_var_idx_corr = np.argsort(feature_var)[-n_features_corr:]
            correlation_matrix = np.corrcoef(imputed_expression.values[:, top_var_idx_corr].T)
            
            sns.heatmap(correlation_matrix,
                       xticklabels=False,
                       yticklabels=False,
                       cmap='coolwarm',
                       center=0,
                       ax=axes[1,1])
            axes[1,1].set_title('Feature Correlation Heatmap\n(Top variable features)')
            
            plt.tight_layout()
            plt.show()

def main():
    # Initialize and run pipeline
    pipeline = MultiModalIntegration("/data/MOFA/")
    
    # Run pipeline steps
    pipeline.load_data()
    pipeline.prepare_rna_data()
    pipeline.prepare_atac_data()
    pipeline.combine_modalities()
    pipeline.train_multivi()
    pipeline.process_latent_space()
    pipeline.visualize_results()

if __name__ == "__main__":
    main()

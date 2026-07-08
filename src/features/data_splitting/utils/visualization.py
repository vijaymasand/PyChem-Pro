import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from .descriptors import preprocess_descriptors_numpy

def generate_pca_plot(descriptor_df, train_indices, test_indices, save_path):
    """
    Generates a 2D PCA projection of the chemical space.
    Highlights the Train (blue) and Test (orange) subsets.
    Saves the plot as an image.
    """
    try:
        X = preprocess_descriptors_numpy(descriptor_df)
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X)
        
        plt.figure(figsize=(6, 4.5))
        
        # Create boolean mask
        train_mask = np.zeros(X.shape[0], dtype=bool)
        train_mask[train_indices] = True
        
        # Plot Train
        plt.scatter(X_pca[train_mask, 0], X_pca[train_mask, 1], 
                    c='#2b5c8f', label=f'Train Set (N={len(train_indices)})', 
                    alpha=0.6, edgecolors='none', s=25)
        # Plot Test
        plt.scatter(X_pca[~train_mask, 0], X_pca[~train_mask, 1], 
                    c='#d95f02', label=f'Test Set (N={len(test_indices)})', 
                    alpha=0.8, edgecolors='none', s=25)
        
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
        plt.title('Dataset Split in Chemical Space (PCA)')
        plt.legend(loc='best')
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(save_path, dpi=150)
        plt.close()
        return True
    except Exception as e:
        print(f"Error generating PCA plot: {e}")
        return False

def compute_split_stats(df, train_indices, test_indices, target_col=None):
    """
    Computes summary statistics comparing the Train and Test subsets.
    If a target_col is provided and numeric, evaluates distribution alignment.
    """
    n_train = len(train_indices)
    n_test = len(test_indices)
    n_total = n_train + n_test
    
    stats = {
        "total_count": n_total,
        "train_count": n_train,
        "test_count": n_test,
        "train_pct": (n_train / n_total) * 100 if n_total > 0 else 0,
        "test_pct": (n_test / n_total) * 100 if n_total > 0 else 0,
    }
    
    if target_col and target_col in df.columns:
        # Check if target is numeric
        target_series = pd.to_numeric(df[target_col], errors='coerce')
        if not target_series.isna().all():
            train_targets = target_series.iloc[train_indices].dropna().values
            test_targets = target_series.iloc[test_indices].dropna().values
            
            if len(train_targets) > 0 and len(test_targets) > 0:
                stats["target_name"] = target_col
                stats["train_mean"] = np.mean(train_targets)
                stats["train_std"] = np.std(train_targets)
                stats["test_mean"] = np.mean(test_targets)
                stats["test_std"] = np.std(test_targets)
                
                # Perform Kolmogorov-Smirnov test to check distribution similarity
                ks_stat, ks_pval = ks_2samp(train_targets, test_targets)
                stats["ks_pvalue"] = ks_pval
                stats["ks_stat"] = ks_stat
                
    return stats

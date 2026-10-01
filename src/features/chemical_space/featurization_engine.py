"""
Chemical Space Visualization - Featurization Engine
Uses PyChem-Pro's native descriptor calculator and SMILES parser.
"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler

from src.features.smiles_parser.services.parser import parse_smiles
from src.features.descriptor_calculator.descriptor_engine import DescriptorEngine
from src.features.descriptor_calculator.descriptor_types import DescriptorCategory

class FeaturizationEngine:
    def __init__(self):
        # Disable cache to prevent unbounded memory growth during large batch processing
        self.engine = DescriptorEngine(enable_cache=False)

    def featurize_dataset(self, df: pd.DataFrame, smiles_col: str, use_fingerprints: bool, use_descriptors: bool,
                          normalize: bool, progress_callback=None):
        """
        Featurizes the DataFrame natively using PyChem-Pro.
        Returns a numpy array of features X, and a filtered df containing only successfully parsed molecules.
        """
        features_list = []
        valid_indices = []

        total = len(df)
        
        categories = []
        if use_fingerprints:
            categories.append(DescriptorCategory.FINGERPRINTS)
        if use_descriptors:
            categories.extend([
                DescriptorCategory.CONSTITUTIONAL,
                DescriptorCategory.TOPOLOGICAL
            ])

        for idx, row in df.iterrows():
            smiles = row[smiles_col]
            if progress_callback:
                progress_callback(int((idx / total) * 100), f"Featurizing {idx+1}/{total}")

            try:
                # 1. Parse SMILES to Molecule natively
                molecule = parse_smiles(smiles)
                
                # 2. Calculate Descriptors
                results = self.engine.calculate_all(molecule, categories=categories)
                
                # Flatten the results into a dict
                feat_dict = {}
                import hashlib
                for desc_name, res in results.items():
                    val = res.value
                    
                    # Proactively handle 2D Morgan Fingerprint string by hashing it into a numerical vector
                    if desc_name == "MorganFingerprint" and isinstance(val, str):
                        fp_array = [0.0] * 64  # Fold into 64 dimensions to keep PCA fast
                        if val and val != "MACCS_PLACEHOLDER":
                            features = val.split("|")
                            for feature in features:
                                h_idx = int(hashlib.md5(feature.encode('utf-8')).hexdigest()[:8], 16) % 64
                                fp_array[h_idx] += 1.0
                        for i, count in enumerate(fp_array):
                            feat_dict[f"MorganFP_{i}"] = count
                        continue
                        
                    if val is None or pd.isna(val):
                        feat_dict[desc_name] = 0.0
                    elif isinstance(val, (int, float)):
                        if np.isinf(val):
                            feat_dict[desc_name] = 0.0
                        else:
                            feat_dict[desc_name] = float(val)
                    else:
                        try:
                            feat_dict[desc_name] = float(val)
                        except ValueError:
                            # Skip non-numeric placeholder descriptors
                            pass
                
                if feat_dict:
                    features_list.append(feat_dict)
                    valid_indices.append(idx)
            except Exception as e:
                # Silently skip failed molecules, could log it
                pass

        if not features_list:
            return None, df.iloc[0:0] # Return empty
        
        # Convert to DataFrame
        feat_df = pd.DataFrame(features_list)
        
        # Ensure we only keep the rows in original df that successfully parsed
        df_valid = df.iloc[valid_indices].reset_index(drop=True)

        # Impute any remaining NaNs with 0
        feat_df = feat_df.fillna(0.0)

        X = feat_df.values

        # Normalize
        if normalize and X.shape[1] > 0:
            scaler = StandardScaler()
            X = scaler.fit_transform(X)

        return X, df_valid

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from .base import BaseSplitter
from .cadex import CADEXSplitter
from ..utils.descriptors import preprocess_descriptors_numpy


class BorutaSplitter(BaseSplitter):
    """Performs Boruta feature selection using Random Forest and shadow
    features to find target-relevant descriptors, and then performs a
    Kennard-Stone split within the selected feature subspace."""

    def __init__(self, target_ratio: float, n_jobs: int = 1,
                 metric: str = 'euclidean', **kwargs):
        super().__init__(target_ratio, n_jobs, **kwargs)
        self.metric = metric

    def split(self, df, smiles_col, target_col=None, descriptor_df=None):
        if descriptor_df is None or descriptor_df.empty:
            raise ValueError("Boruta splitting requires calculated molecular descriptors.")
        if target_col is None or target_col not in df.columns:
            raise ValueError("Boruta method requires a target (activity) column.")

        # Extract target variable
        y = df[target_col].values
        valid_y_mask = ~pd.isna(df[target_col])
        if not np.any(valid_y_mask):
            raise ValueError(f"Target column '{target_col}' has no valid numeric values.")

        # Preprocess descriptors
        X = preprocess_descriptors_numpy(descriptor_df)

        # Filter for valid target values
        X_valid = X[valid_y_mask]
        y_valid = y[valid_y_mask]

        # Detect classification vs regression based on target characteristics
        unique_y = np.unique(y_valid)
        is_classification = len(unique_y) <= 5 or np.issubdtype(y_valid.dtype, np.bool_)

        n_features = X_valid.shape[1]
        hits = np.zeros(n_features)
        n_trials = 10

        # Run Boruta shadow feature iterations
        for trial in range(n_trials):
            X_shadow = np.copy(X_valid)
            # Shuffle each descriptor column independently to break correlations with target
            for col in range(n_features):
                np.random.shuffle(X_shadow[:, col])

            X_combined = np.hstack([X_valid, X_shadow])

            if is_classification:
                rf = RandomForestClassifier(
                    n_estimators=50, max_depth=5,
                    random_state=trial, n_jobs=self.n_jobs)
            else:
                rf = RandomForestRegressor(
                    n_estimators=50, max_depth=5,
                    random_state=trial, n_jobs=self.n_jobs)

            rf.fit(X_combined, y_valid)
            importances = rf.feature_importances_

            orig_importances = importances[:n_features]
            shadow_importances = importances[n_features:]

            max_shadow = np.max(shadow_importances) if len(shadow_importances) > 0 else 0.0
            hits += (orig_importances > max_shadow)

        # Select original descriptors that outperform shadow features in at least 50% of the runs
        selected_features = np.where(hits >= 5)[0]

        # Fallback: if too few features are selected, take the top 5 by default feature importance
        if len(selected_features) < 2:
            if is_classification:
                rf = RandomForestClassifier(
                    n_estimators=50, max_depth=5,
                    random_state=42, n_jobs=self.n_jobs)
            else:
                rf = RandomForestRegressor(
                    n_estimators=50, max_depth=5,
                    random_state=42, n_jobs=self.n_jobs)
            rf.fit(X_valid, y_valid)
            top_features = np.argsort(rf.feature_importances_)[::-1][:5]
            selected_features = top_features

        # Partition using Kennard-Stone in the Boruta-selected feature space
        selected_desc_df = descriptor_df.iloc[:, selected_features]
        ks_splitter = CADEXSplitter(
            target_ratio=self.target_ratio,
            n_jobs=self.n_jobs, metric=self.metric)
        train_idx, test_idx = ks_splitter.split(df, smiles_col, target_col, selected_desc_df)

        return train_idx, test_idx

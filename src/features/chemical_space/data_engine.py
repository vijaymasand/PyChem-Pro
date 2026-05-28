"""
Chemical Space Visualization - Data Engine
Handles CSV ingestion, data cleansing, and basic tracking.
"""
import pandas as pd
import logging

class DataEngine:
    def __init__(self):
        self.df = None
        self.total_rows = 0
        self.valid_compounds = 0
        self.skipped_rows = 0

    def load_csv(self, file_path: str):
        """Loads a CSV file and calculates basic stats."""
        try:
            self.df = pd.read_csv(file_path)
            self.total_rows = len(self.df)
            return True, "Success"
        except Exception as e:
            # Try a different encoding if it fails
            try:
                self.df = pd.read_csv(file_path, encoding='latin1')
                self.total_rows = len(self.df)
                return True, "Success (Latin1 encoding used)"
            except Exception as e2:
                logging.error(f"Failed to load CSV: {e2}")
                return False, str(e2)

    def get_columns(self):
        """Returns the list of columns."""
        if self.df is not None:
            return list(self.df.columns)
        return []

    def get_numeric_columns(self):
        """Returns columns that are numeric for color mapping."""
        if self.df is not None:
            numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
            return list(self.df.select_dtypes(include=numerics).columns)
        return []

    def cleanse_data(self, smiles_col: str):
        """Cleanses the dataframe: removes NaNs in SMILES, trims whitespace, removes duplicates."""
        if self.df is None or smiles_col not in self.df.columns:
            return False, "Invalid DataFrame or SMILES column."

        initial_len = len(self.df)
        
        # Drop rows where SMILES is NaN
        self.df = self.df.dropna(subset=[smiles_col])
        
        # Strip whitespace from SMILES
        self.df[smiles_col] = self.df[smiles_col].astype(str).str.strip()
        
        # Remove empty strings
        self.df = self.df[self.df[smiles_col] != ""]
        
        # Remove duplicates based on SMILES
        self.df = self.df.drop_duplicates(subset=[smiles_col])

        self.valid_compounds = len(self.df)
        self.skipped_rows = initial_len - self.valid_compounds

        # Reset index
        self.df = self.df.reset_index(drop=True)

        return True, "Data cleansed successfully."

    def get_preview(self, num_rows=5):
        """Returns the first N rows for preview."""
        if self.df is not None:
            return self.df.head(num_rows)
        return None

"""
QSAR Modeler Pro Plugin
A comprehensive QSAR modeling suite featuring:
- OLS & XGBoost Regression
- Genetic Algorithm (GA) Feature Selection
- QSARINS-style Validation (LOO, Ext, Leverage)
- SHAP Explainability & Plot Exports
##1. Split Column (0/1/2)
##In QSAR modeling, it is standard practice to split your dataset into a Training
##Set (to build the model) and a Test Set (to prove the model can predict new,
##unseen molecules).
##The Split Column is a feature that allows you to manually pre-define exactly
##which molecules go into which set using a specific column in your CSV file.
##What the numbers mean:
##When the script reads this column, it looks for specific integers:
##1 (Train): The molecule is placed in the Training Set.
##2 (Test): The molecule is placed in the Test Set (External Validation).
##0 (or blank): The molecule is completely ignored and excluded from the analysis.
##How to use it:
##In your CSV file: Create a new column (you can name it something like Set,
##Split, or Subset). For every row, type a 1 if you want it to train the model, or
##a 2 if you want it held back for testing.
##In the GUI: Under step "2. Variable Setup", click the "Split Column (0/1/2)"
##drop-down menu and select the name of that column.
##If you don't want to use it: Simply leave the drop-down set to (none).
##The engine will automatically place 100% of your valid molecules into the
##Training Set.
"""

import math
import itertools
import threading
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression, Lasso, Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import LeaveOneOut, KFold, RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import make_scorer


try:
    import xgboost as xgb
except ImportError:
    xgb = None

try:
    import shap
except ImportError:
    shap = None

# Strictly using ONLY the imports allowed by the host application's qt_compat
from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QTextEdit, QProgressBar,
    QFileDialog, QMessageBox, QComboBox, Qt, QThread, Signal
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType

# ===================================================================
# QSAR MATH UTILITIES
# ===================================================================
def adjusted_r2(r2: float, n: int, p: int) -> float:
    if n <= p + 1: return float('nan')
    return 1 - (1 - r2) * (n - 1) / (n - p - 1)

def rmse(y_true, y_pred):
    return float(np.sqrt(np.mean((np.asarray(y_true) - np.asarray(y_pred)) ** 2)))

def mae(y_true, y_pred):
    return float(np.mean(np.abs(np.asarray(y_true) - np.asarray(y_pred))))

def concordance_cc(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    mu_t, mu_p = np.mean(y_true), np.mean(y_pred)
    var_t, var_p = np.var(y_true, ddof=1), np.var(y_pred, ddof=1)
    cov_tp = np.cov(y_true, y_pred, ddof=1)[0, 1]
    return float((2 * cov_tp) / (var_t + var_p + (mu_t - mu_p) ** 2))

def loo_q2(y_true, y_pred_loo):
    y_true, yhat = np.asarray(y_true), np.asarray(y_pred_loo)
    ss_res = np.sum((y_true - yhat) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float('nan')

def leverage_hat(X: np.ndarray) -> np.ndarray:
    try:
        XtX_inv = np.linalg.pinv(X.T @ X)
        H = X @ XtX_inv @ X.T
        return np.clip(np.diag(H), 0, 1)
    except Exception:
        return np.zeros(X.shape[0])

def standardize_residuals(residuals: np.ndarray) -> np.ndarray:
    s = np.std(residuals, ddof=1) if len(residuals) > 1 else 1.0
    return residuals / s if s != 0 else np.zeros_like(residuals)

def parse_list_string(s):
    if not isinstance(s, str):
        return s
    s = s.strip()
    if s.startswith('[') and s.endswith(']'):
        content = s[1:-1].strip()
        if not content:
            return []
        parts = content.split(',')
        result = []
        for p in parts:
            p = p.strip()
            if not p:
                continue
            if (p.startswith("'") and p.endswith("'")) or (p.startswith('"') and p.endswith('"')):
                result.append(p[1:-1])
            else:
                try:
                    if '.' in p or 'e' in p.lower():
                        result.append(float(p))
                    else:
                        result.append(int(p))
                except ValueError:
                    result.append(p)
        return result
    return s

def expand_sequence_columns(df: pd.DataFrame) -> pd.DataFrame:
    new_cols = {}
    cols_to_drop = []
    for col in df.columns:
        non_nulls = df[col].dropna()
        if non_nulls.empty: continue
        first_val = non_nulls.iloc[0]
        
        if isinstance(first_val, str) and first_val.startswith('[') and first_val.endswith(']'):
            try:
                parsed = parse_list_string(first_val)
                if isinstance(parsed, (list, tuple)):
                    df[col] = df[col].apply(lambda val: parse_list_string(val) if isinstance(val, str) else val)
                    first_val = parsed
            except:
                pass
                
        if isinstance(first_val, (list, tuple, np.ndarray)):
            cols_to_drop.append(col)
            expanded = pd.DataFrame(df[col].tolist(), index=df.index)
            expanded.columns = [f"{col}_{i}" for i in range(expanded.shape[1])]
            for c in expanded.columns:
                new_cols[c] = expanded[c]
                
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)
        for c, s in new_cols.items():
            df[c] = s
    return df

# ===================================================================
# DATA CONTAINERS
# ===================================================================
@dataclass
class ModelResult:
    size: int
    descriptors: List[str]
    model_type: str
    coef: Optional[np.ndarray]
    intercept: Optional[float]
    r2_tr: float
    r2_adj: float
    rmse_tr: float
    mae_tr: float
    q2_loo: float
    rmse_cv: float
    yhat_tr: np.ndarray
    yhat_loo: np.ndarray
    metrics_ext: Dict[str, float] = field(default_factory=dict)
    hat_diag: Optional[np.ndarray] = None
    std_resid_fit: Optional[np.ndarray] = None
    xgb_best_params_: Optional[Dict[str, Any]] = None
    xgb_model_obj: Any = None
    X_tr_scaled: Optional[np.ndarray] = None

# ===================================================================
# CORE QSAR ENGINE
# ===================================================================
class QSAREngine:
    def __init__(self):
        self.df: Optional[pd.DataFrame] = None
        self.df_test: Optional[pd.DataFrame] = None
        self.split_col: Optional[str] = None
        self.y_col: Optional[str] = None
        self.x_cols: List[str] = []
        self.tr_idx: Optional[np.ndarray] = None
        self.te_idx: Optional[np.ndarray] = None

    def load_dataframe(self, df: pd.DataFrame, missing_val: Optional[float] = -999.0):
        self.df = expand_sequence_columns(df.copy())
        if missing_val is not None:
            self.df = self.df.replace(missing_val, np.nan)

    def load_test_dataframe(self, df: pd.DataFrame, missing_val: Optional[float] = -999.0):
        self.df_test = expand_sequence_columns(df.copy())
        if missing_val is not None:
            self.df_test = self.df_test.replace(missing_val, np.nan)

    def set_setup(self, y_col: str, x_cols: List[str], split_col: Optional[str] = None):
        self.y_col, self.x_cols, self.split_col = y_col, x_cols, split_col
        if self.df is None: raise ValueError("No dataset loaded")

        if self.df_test is not None:
            self.tr_idx = np.where(self.df[self.y_col].notna().values)[0]
            self.te_idx = np.where(self.df_test[self.y_col].notna().values)[0]
        elif split_col and split_col in self.df.columns:
            tr_indices = []
            te_indices = []
            for idx, x in enumerate(self.df[split_col].values):
                val = str(x).strip().lower() if pd.notna(x) else ""
                if any(val.startswith(pfx) for pfx in ("train", "tr", "1")):
                    tr_indices.append(idx)
                elif any(val.startswith(pfx) for pfx in ("test", "te", "ext", "2")):
                    te_indices.append(idx)
            self.tr_idx = np.array(tr_indices, dtype=int)
            self.te_idx = np.array(te_indices, dtype=int)
        else:
            self.tr_idx = np.where(self.df[self.y_col].notna().values)[0]
            self.te_idx = np.array([], dtype=int)

    def _prepare_xy(self, cols: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, StandardScaler]:
        df_tr = self.df.iloc[self.tr_idx][[self.y_col] + cols].dropna()
        y_tr = df_tr[self.y_col].values.astype(float)
        X_tr = df_tr[cols].values.astype(float)

        if self.df_test is not None:
            df_te = self.df_test.iloc[self.te_idx][[self.y_col] + cols].dropna()
            y_te = df_te[self.y_col].values.astype(float)
            X_te = df_te[cols].values.astype(float)
        elif len(self.te_idx) > 0:
            df_te = self.df.iloc[self.te_idx][[self.y_col] + cols].dropna()
            y_te = df_te[self.y_col].values.astype(float)
            X_te = df_te[cols].values.astype(float)
        else:
            y_te, X_te = np.array([], dtype=float), np.zeros((0, len(cols)), dtype=float)

        scaler = StandardScaler()
        return scaler.fit_transform(X_tr), y_tr, scaler.transform(X_te) if len(X_te) else X_te, y_te, scaler

    def train_model(self, cols: List[str], model_type: str, config: Optional[dict] = None) -> Optional[ModelResult]:
        if not cols: return None
        if config is None: config = {}
        X_tr, y_tr, X_te, y_te, scaler = self._prepare_xy(cols)
        n, p = X_tr.shape
        if n <= p + 1: return None

        model_obj = None
        if model_type == 'OLS':
            model_obj = LinearRegression()
        elif model_type == 'Ridge':
            alpha = float(config.get('ridge_alpha', 1.0))
            model_obj = Ridge(alpha=alpha)
        elif model_type == 'Lasso':
            alpha = float(config.get('lasso_alpha', 1.0))
            model_obj = Lasso(alpha=alpha)
        elif model_type == 'PLS':
            n_comp = int(config.get('pls_components', 2))
            model_obj = PLSRegression(n_components=min(n_comp, p))
        elif model_type == 'RF':
            n_est = int(config.get('rf_estimators', 100))
            model_obj = RandomForestRegressor(n_estimators=n_est, random_state=42)
        elif model_type == 'SVR':
            C = float(config.get('svr_C', 1.0))
            eps = float(config.get('svr_epsilon', 0.1))
            model_obj = SVR(C=C, epsilon=eps)
        elif model_type == 'XGB':
            if xgb is None: raise RuntimeError("XGBoost not installed.")
            def parse_range(val, is_float=False):
                parts = [v.strip() for v in str(val).split(',')]
                if len(parts) == 1: return [float(parts[0]) if is_float else int(parts[0])]
                return np.linspace(float(parts[0]), float(parts[1]), 10) if is_float else np.arange(int(parts[0]), int(parts[1]), max(1, (int(parts[1])-int(parts[0]))//5))
            
            xgbr = xgb.XGBRegressor(objective='reg:squarederror', random_state=42)
            param_dist = {
                'n_estimators': parse_range(config.get('n_estimators', '100,500'), False),
                'max_depth': parse_range(config.get('max_depth', '2,6'), False),
                'learning_rate': parse_range(config.get('learning_rate', '0.01,0.2'), True)
            }
            cv_folds = int(config.get('cv_folds', 5))
            n_iter = int(config.get('n_iter', 20))
            rs = RandomizedSearchCV(xgbr, param_distributions=param_dist, n_iter=n_iter,
                                    scoring=make_scorer(lambda y, yhat: -rmse(y, yhat)),
                                    cv=KFold(n_splits=min(cv_folds, len(y_tr)), shuffle=True, random_state=42), random_state=42)
            rs.fit(X_tr, y_tr)
            model_obj = rs.best_estimator_

        if model_type != 'XGB':
            model_obj.fit(X_tr, y_tr)
        
        yhat_tr = np.ravel(model_obj.predict(X_tr))

        yhat_loo = np.zeros(n)
        if n <= 50:
            cv_strategy = LeaveOneOut()
        else:
            cv_strategy = KFold(n_splits=5, shuffle=True, random_state=42)
            
        for tr, te in cv_strategy.split(X_tr):
            if model_type == 'OLS':
                fold_model = LinearRegression()
            elif model_type == 'Ridge':
                fold_model = Ridge(alpha=float(config.get('ridge_alpha', 1.0)))
            elif model_type == 'Lasso':
                fold_model = Lasso(alpha=float(config.get('lasso_alpha', 1.0)))
            elif model_type == 'PLS':
                fold_model = PLSRegression(n_components=min(int(config.get('pls_components', 2)), p))
            elif model_type == 'RF':
                fold_model = RandomForestRegressor(n_estimators=int(config.get('rf_estimators', 100)), random_state=42)
            elif model_type == 'SVR':
                fold_model = SVR(C=float(config.get('svr_C', 1.0)), epsilon=float(config.get('svr_epsilon', 0.1)))
            elif model_type == 'XGB':
                fold_model = xgb.XGBRegressor(**model_obj.get_params())
            
            fold_model.fit(X_tr[tr], y_tr[tr])
            yhat_loo[te] = np.ravel(fold_model.predict(X_tr[te]))

        metrics_ext = {}
        if len(y_te) > 0:
            yhat_te = np.ravel(model_obj.predict(X_te))
            ss_res, ss_tot = np.sum((y_te - yhat_te) ** 2), np.sum((y_te - np.mean(y_te)) ** 2)
            r2ext = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
            
            den_f1 = np.sum((y_te - np.mean(y_tr)) ** 2)
            q2_f1 = 1 - (ss_res / den_f1) if den_f1 > 0 else float('nan')
            q2_f2 = r2ext
            den_f3 = np.sum((y_tr - np.mean(y_tr)) ** 2) / n
            q2_f3 = 1 - ((ss_res / len(y_te)) / den_f3) if den_f3 > 0 else float('nan')
            
            try:
                k_val = np.sum(y_te * yhat_te) / np.sum(y_te ** 2) if np.sum(y_te ** 2) > 0 else 1.0
                k_prime_val = np.sum(y_te * yhat_te) / np.sum(yhat_te ** 2) if np.sum(yhat_te ** 2) > 0 else 1.0
            except:
                k_val, k_prime_val = 1.0, 1.0
                
            r2_zero = 1 - (np.sum((y_te - k_val * y_te)**2) / np.sum(y_te**2)) if np.sum(y_te**2) != 0 else 0
            r_prime2_zero = 1 - (np.sum((y_te - k_prime_val * yhat_te)**2) / np.sum(yhat_te**2)) if np.sum(yhat_te**2) != 0 else 0
            r2_zero_diff = abs(r2_zero - r_prime2_zero)
            
            gt_passed = (
                (r2ext > 0.5 if not np.isnan(r2ext) else False) and
                r2_zero_diff < 0.3 and
                0.85 <= k_val <= 1.15 and
                0.85 <= k_prime_val <= 1.15
            )
            
            metrics_ext = {
                "R2ext": r2ext, "RMSEext": rmse(y_te, yhat_te), "MAEext": mae(y_te, yhat_te), "CCCext": concordance_cc(y_te, yhat_te),
                "Q2_F1": q2_f1, "Q2_F2": q2_f2, "Q2_F3": q2_f3,
                "GT_k": k_val, "GT_k_prime": k_prime_val, "GT_diff": r2_zero_diff, "GT_passed": gt_passed
            }

        coef = None
        intercept = None
        if model_type in ('OLS', 'Ridge', 'Lasso'):
            coef = model_obj.coef_
            intercept = float(model_obj.intercept_)
        elif model_type == 'PLS':
            coef = model_obj.coef_.ravel()
            intercept = float(model_obj._y_mean[0])

        r2_tr = 1 - np.sum((y_tr - yhat_tr)**2) / np.sum((y_tr - np.mean(y_tr))**2)
        
        X_tr_1 = np.c_[np.ones((n, 1)), X_tr]
        hat_diag = leverage_hat(X_tr_1)
        std_resid_fit = standardize_residuals(y_tr - yhat_tr)
        xgb_params = model_obj.get_params() if model_type == 'XGB' else None

        return ModelResult(
            size=len(cols), descriptors=cols, model_type=model_type,
            coef=coef, intercept=intercept,
            r2_tr=r2_tr, r2_adj=adjusted_r2(r2_tr, n, p),
            rmse_tr=rmse(y_tr, yhat_tr), mae_tr=mae(y_tr, yhat_tr),
            q2_loo=loo_q2(y_tr, yhat_loo), rmse_cv=rmse(y_tr, yhat_loo),
            yhat_tr=yhat_tr, yhat_loo=yhat_loo, metrics_ext=metrics_ext,
            hat_diag=hat_diag, std_resid_fit=std_resid_fit,
            xgb_best_params_=xgb_params, xgb_model_obj=model_obj, X_tr_scaled=X_tr
        )

    def ga_feature_selection(self, model_type: str, min_f: int, max_f: int, pop_size: int, gen: int, config: dict, logger_callback) -> ModelResult:
        rng = np.random.RandomState(42)
        d = len(self.x_cols)

        def fitness(mask):
            cols = [self.x_cols[i] for i, b in enumerate(mask) if b]
            try:
                res = self.train_model(cols, model_type, config)
                if res is None or np.isnan(res.q2_loo): return -np.inf, None
                return (res.q2_loo, -res.rmse_cv, -len(cols)), res
            except Exception: return -np.inf, None

        population = []
        logger_callback("Initializing GA Population...")
        while len(population) < pop_size:
            m = np.zeros(d, dtype=bool); m[rng.choice(d, rng.randint(min_f, max_f+1), replace=False)] = True
            s, res = fitness(m)
            if res: population.append((s, m, res))

        best_tuple = max(population, key=lambda t: t[0])

        for g in range(gen):
            new_pop = []
            while len(new_pop) < pop_size:
                p1 = max(rng.choice(population, 3, replace=False), key=lambda t: t[0])[1]
                p2 = max(rng.choice(population, 3, replace=False), key=lambda t: t[0])[1]
                c1, c2 = p1.copy(), p2.copy()

                if rng.rand() < 0.9:
                    cx = rng.rand(d) < 0.5
                    c1, c2 = np.where(cx, p1, p2), np.where(cx, p2, p1)

                for c in (c1, c2):
                    c[rng.rand(d) < 0.08] = ~c[rng.rand(d) < 0.08]
                    k = c.sum()
                    if k < min_f: c[rng.choice(np.where(~c)[0], min_f - k, replace=False)] = True
                    elif k > max_f: c[rng.choice(np.where(c)[0], k - max_f, replace=False)] = False

                    s, res = fitness(c)
                    if res:
                        new_pop.append((s, c.copy(), res))
                        if s > best_tuple[0]: best_tuple = (s, c.copy(), res)

            population = sorted(new_pop, key=lambda t: t[0], reverse=True)[:pop_size]
            if g % 5 == 0: logger_callback(f"GA Gen {g}/{gen}: Best Q2 = {best_tuple[0][0]:.4f}")

        return best_tuple[2]

# ===================================================================
# BACKGROUND WORKER
# ===================================================================
class ModelingWorker(QThread):
    log_signal = Signal(str)
    result_signal = Signal(object)
    error_signal = Signal(str)

    def __init__(self, engine, use_ga, model_type, config):
        super().__init__()
        self.engine = engine
        self.use_ga = use_ga
        self.model_type = model_type
        self.config = config

    def run(self):
        try:
            if self.use_ga:
                self.log_signal.emit("Starting Genetic Algorithm Feature Selection...")
                res = self.engine.ga_feature_selection(
                    self.model_type,
                    int(self.config['ga_min']), int(self.config['ga_max']),
                    int(self.config['ga_pop']), int(self.config['ga_gen']),
                    self.config,
                    lambda msg: self.log_signal.emit(msg)
                )
            else:
                self.log_signal.emit(f"Fitting single {self.model_type} model...")
                res = self.engine.train_model(self.engine.x_cols, self.model_type, self.config)

            if res:
                self.result_signal.emit(res)
            else:
                self.error_signal.emit("Model failed to fit (Insufficient data or collinearity).")
        except Exception as e:
            self.error_signal.emit(str(e))

# ===================================================================
# QSAR UI WIDGET
# ===================================================================
class QsarModelerWidget(PluginWidget):
    def __init__(self, plugin):
        super().__init__(plugin)
        self.engine = QSAREngine()
        self.current_model = None
        self.worker = None
        self.setup_ui()

    def setup_ui(self):
        self.widget = QWidget()
        main_layout = QVBoxLayout(self.widget)

        # --- 1 & 2. Data Loading & Variable Setup (Horizontal Layout) ---
        top_layout = QHBoxLayout()

        # Step 1: Data Source
        data_widget = QWidget()
        data_layout = QHBoxLayout(data_widget)
        data_layout.setContentsMargins(0, 0, 0, 0)
        data_layout.addWidget(QLabel("<b>1. Data Source:</b>"))
        self.btn_load = QPushButton("Load Train CSV")
        self.btn_load.clicked.connect(self.load_data)
        data_layout.addWidget(self.btn_load)
        
        self.btn_load_test = QPushButton("Load Test CSV (Optional)")
        self.btn_load_test.clicked.connect(self.load_test_data)
        data_layout.addWidget(self.btn_load_test)
        
        self.lbl_file = QLabel("No files loaded")
        data_layout.addWidget(self.lbl_file)
        top_layout.addWidget(data_widget)

        top_layout.addSpacing(20)

        # Step 2: Variable Setup
        setup_widget = QWidget()
        setup_layout = QHBoxLayout(setup_widget)
        setup_layout.setContentsMargins(0, 0, 0, 0)
        setup_layout.addWidget(QLabel("<b>2. Variable Setup:</b>"))
        setup_layout.addWidget(QLabel("Target (Y):"))
        self.cmb_y = QComboBox()
        setup_layout.addWidget(self.cmb_y)
        setup_layout.addWidget(QLabel("Split Column:"))
        self.cmb_split = QComboBox()
        setup_layout.addWidget(self.cmb_split)

        setup_layout.addWidget(QLabel("Split Algo:"))
        from src.features.data_splitting.split_engine import DataSplitEngine
        self.cmb_split_algo = QComboBox()
        self.cmb_split_algo.addItems(DataSplitEngine.available_algorithms())
        setup_layout.addWidget(self.cmb_split_algo)

        setup_layout.addWidget(QLabel("SMILES Col:"))
        self.cmb_smiles_col = QComboBox()
        setup_layout.addWidget(self.cmb_smiles_col)

        setup_layout.addWidget(QLabel("Train Ratio:"))
        self.txt_train_ratio = QComboBox()
        self.txt_train_ratio.addItems(["0.8", "0.7", "0.75", "0.9", "0.5"])
        self.txt_train_ratio.setEditable(True)
        setup_layout.addWidget(self.txt_train_ratio)

        self.btn_run_split = QPushButton("Split")
        self.btn_run_split.clicked.connect(self.run_inline_split)
        setup_layout.addWidget(self.btn_run_split)

        top_layout.addWidget(setup_widget)
        top_layout.addStretch()
        main_layout.addLayout(top_layout)

        # --- 3. Descriptor Selection (Dual Table UX) ---
        main_layout.addWidget(QLabel("<b>3. Descriptor Selection:</b>"))
        desc_layout = QHBoxLayout()

        self.tbl_avail = QTableWidget()
        self.tbl_avail.setColumnCount(1)
        self.tbl_avail.setHorizontalHeaderLabels(["Available"])
        self.tbl_avail.horizontalHeader().setStretchLastSection(True)
        desc_layout.addWidget(self.tbl_avail)

        btn_move_layout = QVBoxLayout()
        self.btn_add = QPushButton(">>")
        self.btn_add.clicked.connect(self.add_descriptors)
        self.btn_rem = QPushButton("<<")
        self.btn_rem.clicked.connect(self.rem_descriptors)
        btn_move_layout.addStretch()
        btn_move_layout.addWidget(self.btn_add)
        btn_move_layout.addWidget(self.btn_rem)
        btn_move_layout.addStretch()
        desc_layout.addLayout(btn_move_layout)

        self.tbl_sel = QTableWidget()
        self.tbl_sel.setColumnCount(1)
        self.tbl_sel.setHorizontalHeaderLabels(["Selected"])
        self.tbl_sel.horizontalHeader().setStretchLastSection(True)
        desc_layout.addWidget(self.tbl_sel)

        main_layout.addLayout(desc_layout)

        # --- 4. Model Settings (Algorithms & Hyperparameters) ---
        main_layout.addWidget(QLabel("<b>4. Algorithm & Hyperparameters:</b>"))
        algo_layout = QHBoxLayout()
        algo_layout.addWidget(QLabel("Model:"))
        self.cmb_model = QComboBox()
        self.cmb_model.addItems(["OLS", "Ridge", "Lasso", "PLS", "RF", "SVR", "XGB"])
        algo_layout.addWidget(self.cmb_model)

        algo_layout.addWidget(QLabel("Feature Selection:"))
        self.cmb_ga = QComboBox()
        self.cmb_ga.addItems(["Use All Selected", "Run Genetic Algorithm (GA)"])
        algo_layout.addWidget(self.cmb_ga)
        algo_layout.addStretch()
        main_layout.addLayout(algo_layout)

        self.settings_tbl = QTableWidget()
        self.settings_tbl.setColumnCount(2)
        self.settings_tbl.setHorizontalHeaderLabels(["Setting", "Value"])
        self.settings_tbl.horizontalHeader().setStretchLastSection(True)
        self.settings_tbl.setMaximumHeight(150)

        default_settings = [
            ("ga_min", "1"), ("ga_max", "8"), ("ga_pop", "30"), ("ga_gen", "40"),
            ("n_estimators", "100, 500"), ("max_depth", "2, 6"),
            ("learning_rate", "0.01, 0.2"), ("n_iter", "20"), ("cv_folds", "5"),
            ("ridge_alpha", "1.0"), ("lasso_alpha", "1.0"), ("pls_components", "2"),
            ("rf_estimators", "100"), ("svr_C", "1.0"), ("svr_epsilon", "0.1")
        ]
        self.settings_tbl.setRowCount(len(default_settings))
        for i, (k, v) in enumerate(default_settings):
            k_item = QTableWidgetItem(k)
            k_item.setFlags(k_item.flags() & ~Qt.ItemIsEditable)
            self.settings_tbl.setItem(i, 0, k_item)
            self.settings_tbl.setItem(i, 1, QTableWidgetItem(v))
        main_layout.addWidget(self.settings_tbl)

        # --- 5. Run & Logs ---
        run_layout = QHBoxLayout()
        self.btn_apply_setup = QPushButton("Apply Setup")
        self.btn_apply_setup.clicked.connect(self.apply_setup)
        self.btn_apply_setup.setStyleSheet("font-weight: bold; padding: 8px;")
        run_layout.addWidget(self.btn_apply_setup)

        self.btn_run = QPushButton("🚀 Train Model")
        self.btn_run.setStyleSheet("background-color: #2980b9; color: white; font-weight: bold; padding: 8px;")
        self.btn_run.clicked.connect(self.run_model)
        run_layout.addWidget(self.btn_run)
        main_layout.addLayout(run_layout)

        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        main_layout.addWidget(self.txt_log)

        # --- 6. Export ---
        export_layout = QHBoxLayout()
        self.btn_export_csv = QPushButton("Export Predictions (CSV)")
        self.btn_export_csv.clicked.connect(self.export_csv)
        self.btn_export_csv.setEnabled(False)
        export_layout.addWidget(self.btn_export_csv)

        self.btn_export_plots = QPushButton("Generate & Save Plots (PNG)")
        self.btn_export_plots.clicked.connect(self.export_plots)
        self.btn_export_plots.setEnabled(False)
        export_layout.addWidget(self.btn_export_plots)

        main_layout.addLayout(export_layout)

    def load_data(self):
        path, _ = QFileDialog.getOpenFileName(self.widget, "Select Dataset", "", "CSV/TSV Files (*.csv *.tsv *.txt)")
        if not path: return

        try:
            df = pd.read_csv(path, sep=None, engine='python')
            self.engine.load_dataframe(df)
            self.engine.df_test = None

            clean_name = path.replace("\\", "/").split("/")[-1]
            self.lbl_file.setText(f"Train: {clean_name}")

            cols = list(self.engine.df.columns)
            self.cmb_y.clear()
            self.cmb_y.addItems(cols)
            self.cmb_split.clear()
            self.cmb_split.addItems(["(none)"] + cols)

            self.cmb_smiles_col.clear()
            self.cmb_smiles_col.addItems(cols)
            for col in cols:
                if "smiles" in col.lower():
                    self.cmb_smiles_col.setCurrentText(col)
                    break

            self.tbl_avail.setRowCount(len(cols))
            for i, c in enumerate(cols):
                self.tbl_avail.setItem(i, 0, QTableWidgetItem(c))
            self.tbl_sel.setRowCount(0)

            self.txt_log.append(f"Loaded Train CSV: {len(self.engine.df)} rows and {len(cols)} columns.")
        except Exception as e:
            QMessageBox.critical(self.widget, "Load Error", str(e))

    def load_test_data(self):
        path, _ = QFileDialog.getOpenFileName(self.widget, "Select Test Dataset", "", "CSV/TSV Files (*.csv *.tsv *.txt)")
        if not path: return

        try:
            df = pd.read_csv(path, sep=None, engine='python')
            self.engine.load_test_dataframe(df)

            clean_name = path.replace("\\", "/").split("/")[-1]
            train_name = self.lbl_file.text().split(" | ")[0]
            self.lbl_file.setText(f"{train_name} | Test: {clean_name}")

            self.txt_log.append(f"Loaded Test CSV: {len(df)} rows.")
        except Exception as e:
            QMessageBox.critical(self.widget, "Load Error", str(e))

    def run_inline_split(self):
        if self.engine.df is None:
            QMessageBox.warning(self.widget, "Split Error", "Please load a Train CSV dataset first.")
            return

        algo = self.cmb_split_algo.currentText()
        smiles = self.cmb_smiles_col.currentText()
        try:
            ratio = float(self.txt_train_ratio.currentText())
        except ValueError:
            QMessageBox.warning(self.widget, "Split Error", "Invalid Train Ratio. Please enter a float (e.g. 0.8)")
            return

        if not smiles:
            QMessageBox.warning(self.widget, "Split Error", "Please select the SMILES column.")
            return

        y_col = self.cmb_y.currentText()
        target_col = y_col if y_col else None

        try:
            split_engine = DataSplitEngine()
            split_engine.load_dataframe(self.engine.df)
            res = split_engine.split(
                algorithm=algo,
                target_ratio=ratio,
                smiles_col=smiles,
                target_col=target_col,
                desc_mode=DescriptorMode.FINGERPRINTS,
                n_jobs=1
            )
            self.engine.load_dataframe(res.annotated_df)

            cols = list(self.engine.df.columns)
            self.cmb_split.clear()
            self.cmb_split.addItems(["(none)"] + cols)
            self.cmb_split.setCurrentText("Split_Status")

            self.tbl_avail.setRowCount(len(cols))
            for i, c in enumerate(cols):
                self.tbl_avail.setItem(i, 0, QTableWidgetItem(c))
            self.tbl_sel.setRowCount(0)

            self.txt_log.append(f"Inline split succeeded using '{algo}' (Train: {len(res.train_indices)}, Test: {len(res.test_indices)}). Created and selected 'Split_Status' column.")
            self.apply_setup()
        except Exception as e:
            QMessageBox.critical(self.widget, "Split Error", str(e))

    def add_descriptors(self):
        for item in self.tbl_avail.selectedItems():
            row = self.tbl_sel.rowCount()
            self.tbl_sel.insertRow(row)
            self.tbl_sel.setItem(row, 0, QTableWidgetItem(item.text()))
            self.tbl_avail.removeRow(item.row())

    def rem_descriptors(self):
        for item in self.tbl_sel.selectedItems():
            row = self.tbl_avail.rowCount()
            self.tbl_avail.insertRow(row)
            self.tbl_avail.setItem(row, 0, QTableWidgetItem(item.text()))
            self.tbl_sel.removeRow(item.row())

    def apply_setup(self):
        y = self.cmb_y.currentText()
        split = self.cmb_split.currentText()
        split_col = None if split == "(none)" else split

        x_cols = []
        for i in range(self.tbl_sel.rowCount()):
            x_cols.append(self.tbl_sel.item(i, 0).text())

        if not y or not x_cols:
            QMessageBox.warning(self.widget, "Setup Error", "Select a Target (Y) and at least one descriptor.")
            return

        try:
            self.engine.set_setup(y, x_cols, split_col)
            self.txt_log.append(f"Setup Applied: Target='{y}', Descriptors={len(x_cols)}, Train N={len(self.engine.tr_idx)}")
        except Exception as e:
            QMessageBox.critical(self.widget, "Setup Error", str(e))

    def run_model(self):
        self.apply_setup()
        if not self.engine.x_cols:
            return

        self.btn_run.setEnabled(False)
        self.txt_log.clear()

        config = {}
        for i in range(self.settings_tbl.rowCount()):
            config[self.settings_tbl.item(i, 0).text()] = self.settings_tbl.item(i, 1).text()

        use_ga = self.cmb_ga.currentText() == "Run Genetic Algorithm (GA)"
        model_type = self.cmb_model.currentText()

        self.worker = ModelingWorker(self.engine, use_ga, model_type, config)
        self.worker.log_signal.connect(self.txt_log.append)
        self.worker.result_signal.connect(self.handle_result)
        self.worker.error_signal.connect(self.handle_error)
        self.worker.start()

    def handle_result(self, res: ModelResult):
        self.current_model = res
        self.txt_log.append("\n--- MODEL SUMMARY ---")
        self.txt_log.append(f"Algorithm: {res.model_type}")
        self.txt_log.append(f"Features ({len(res.descriptors)}): {', '.join(res.descriptors)}")
        self.txt_log.append(f"Train R2: {res.r2_tr:.3f} | Train RMSE: {res.rmse_tr:.3f} | Train MAE: {res.mae_tr:.3f}")
        self.txt_log.append(f"CV Q2: {res.q2_loo:.3f} | CV RMSE: {res.rmse_cv:.3f}")

        if res.metrics_ext:
            self.txt_log.append(f"Test R2 (F2): {res.metrics_ext.get('R2ext'):.3f} | Test RMSE: {res.metrics_ext.get('RMSEext'):.3f} | Test MAE: {res.metrics_ext.get('MAEext'):.3f}")
            self.txt_log.append(f"Test CCC: {res.metrics_ext.get('CCCext'):.3f}")
            self.txt_log.append(f"Test Q2 (F1): {res.metrics_ext.get('Q2_F1'):.3f} | Test Q2 (F3): {res.metrics_ext.get('Q2_F3'):.3f}")
            self.txt_log.append(f"Golbraikh-Tropsha k: {res.metrics_ext.get('GT_k'):.3f} | k': {res.metrics_ext.get('GT_k_prime'):.3f} | R2_diff: {res.metrics_ext.get('GT_diff'):.3f}")
            gt_status = "PASSED" if res.metrics_ext.get('GT_passed') else "FAILED"
            self.txt_log.append(f"Golbraikh-Tropsha Validation Status: {gt_status}")

        if res.model_type == 'XGB' and res.xgb_best_params_:
            self.txt_log.append(f"\nBest XGB Params: {res.xgb_best_params_}")

        self.btn_run.setEnabled(True)
        self.btn_export_csv.setEnabled(True)
        self.btn_export_plots.setEnabled(True)
        QMessageBox.information(self.widget, "Success", "Modeling completed successfully!")

    def handle_error(self, err_msg):
        self.txt_log.append(f"\n<span style='color:red'>ERROR: {err_msg}</span>")
        self.btn_run.setEnabled(True)
        QMessageBox.critical(self.widget, "Error", err_msg)

    def export_csv(self):
        if not self.current_model: return
        path, _ = QFileDialog.getSaveFileName(self.widget, "Save Predictions", "QSAR_Predictions.csv", "CSV Files (*.csv)")
        if not path: return

        y_tr = self.engine.df.iloc[self.engine.tr_idx][self.engine.y_col].values
        df_exp = pd.DataFrame({"Subset": "Train", "Observed": y_tr, "Predicted_Fit": self.current_model.yhat_tr, "Predicted_CV": self.current_model.yhat_loo})

        if len(self.engine.te_idx) > 0:
            if self.engine.df_test is not None:
                y_te = self.engine.df_test.iloc[self.engine.te_idx][self.engine.y_col].values
            else:
                y_te = self.engine.df.iloc[self.engine.te_idx][self.engine.y_col].values
                
            _, _, X_te_s, _, _ = self.engine._prepare_xy(self.current_model.descriptors)
            yhat_te = np.ravel(self.current_model.xgb_model_obj.predict(X_te_s))

            df_te = pd.DataFrame({"Subset": "Test", "Observed": y_te, "Predicted_Fit": yhat_te, "Predicted_CV": np.nan})
            df_exp = pd.concat([df_exp, df_te])

        df_exp.to_csv(path, index=False)
        self.txt_log.append(f"Saved CSV to {path}")

    def export_plots(self):
        if not self.current_model: return
        directory = QFileDialog.getExistingDirectory(self.widget, "Select Directory to Save Plots")
        if not directory: return

        directory = directory.replace("\\", "/")
        res = self.current_model
        y_tr = self.engine.df.iloc[self.engine.tr_idx][self.engine.y_col].values

        # 1. Fit Plot
        plt.figure(figsize=(6, 5))
        plt.scatter(y_tr, res.yhat_tr, alpha=0.7, edgecolors='k')
        lims = [min(y_tr.min(), res.yhat_tr.min()), max(y_tr.max(), res.yhat_tr.max())]
        plt.plot(lims, lims, 'r--')
        plt.xlabel("Observed")
        plt.ylabel("Predicted")
        plt.title(f"Predicted vs Observed ({res.model_type})")
        plt.savefig(f"{directory}/fit_plot.png", dpi=300)
        plt.close()

        # 2. Williams Plot (Applicability Domain)
        if res.hat_diag is not None:
            plt.figure(figsize=(6, 5))
            plt.scatter(res.hat_diag, res.std_resid_fit, alpha=0.7, edgecolors='k')
            plt.axhline(3, color='r', linestyle='--'); plt.axhline(-3, color='r', linestyle='--')
            h_star = 3 * (len(res.descriptors) + 1) / len(res.hat_diag)
            plt.axvline(h_star, color='r', linestyle='--')
            plt.xlabel("Leverage (h)")
            plt.ylabel("Standardized Residuals")
            plt.title(f"Williams Plot ({res.model_type})")
            plt.savefig(f"{directory}/williams_plot.png", dpi=300)
            plt.close()

        # 3. SHAP Plot (If XGBoost)
        if res.model_type == 'XGB' and shap is not None:
            try:
                plt.figure(figsize=(8, 6))
                explainer = shap.TreeExplainer(res.xgb_model_obj)
                shap_values = explainer.shap_values(res.X_tr_scaled)
                shap.summary_plot(shap_values, res.X_tr_scaled, feature_names=res.descriptors, show=False)
                plt.savefig(f"{directory}/shap_summary.png", dpi=300, bbox_inches='tight')
                plt.close()
            except Exception as e:
                self.txt_log.append(f"Could not generate SHAP plot: {e}")

        self.txt_log.append(f"Plots saved to {directory}")
        QMessageBox.information(self.widget, "Export Successful", f"Plots successfully saved to:\n{directory}")


class QsarModelerPlugin(BasePlugin):
    """
    QSAR Modeler Plugin Wrapper.
    """
    def __init__(self):
        super().__init__(PluginInfo(
            name="QSAR Modeler Pro",
            version="1.0.0",
            description="Advanced QSAR modeling with GA feature selection and XGBoost",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.widget = None

    def get_info(self) -> PluginInfo:
        return self.info

    def create_widget(self) -> 'QsarModelerWidget':
        if self.widget is None:
            self.widget = QsarModelerWidget(self)
        return self.widget

    def initialize(self):
        self.logger.info("QSAR Modeler plugin initialized successfully")
        return True

    def cleanup(self):
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("QSAR Modeler plugin cleaned up")
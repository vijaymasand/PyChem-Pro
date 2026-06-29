"""
Hybrid RFA + MARS QSAR Plugin

Features:
- Red Fox Algorithm (RFA) Feature Selection
- SplineTransformer (MARS-like) Non-linear Modeling
- Multiple ML models (MLR, Ridge, Lasso, PLS, RF, SVR, XGBoost)
- Golbraikh-Tropsha Criteria Validation
- Extended QSAR Validation Metrics
- Applicability Domain (Williams Plot)
- Interactive Graphics
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import KFold, cross_val_score, cross_val_predict, LeaveOneOut
from sklearn.linear_model import LinearRegression, Lasso, Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import SplineTransformer, StandardScaler
from sklearn.pipeline import Pipeline
from scipy.linalg import pinv

try:
    from xgboost import XGBRegressor
    XGB_AVAILABLE = True
except ImportError:
    XGB_AVAILABLE = False

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Strictly using ONLY the allowed imports from qt_compat
from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QTextEdit, QProgressBar,
    QFileDialog, QMessageBox, QComboBox, Qt, QThread, Signal,
    QTabWidget, QCheckBox, QSpinBox, QSplitter, QGroupBox, QListWidget, QListWidgetItem
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType

# ===================================================================
# QSAR & AD HELPER FUNCTIONS
# ===================================================================

def calculate_ccc(y_true, y_pred):
    if len(y_true) < 2: return 0.0
    cor = np.corrcoef(y_true, y_pred)[0][1]
    mean_true, mean_pred = np.mean(y_true), np.mean(y_pred)
    var_true, var_pred = np.var(y_true), np.var(y_pred)
    sd_true, sd_pred = np.std(y_true), np.std(y_pred)
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2
    return numerator / denominator if denominator != 0 else 0

def calculate_q2(y_true, y_pred):
    press = np.sum((y_true - y_pred) ** 2)
    tss = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1 - (press / tss) if tss != 0 else 0

def calculate_rmse(y_true, y_pred):
    return np.sqrt(mean_squared_error(y_true, y_pred))

def calculate_q2_f1(y_true_ex, y_pred_ex, y_train):
    num = np.sum((y_true_ex - y_pred_ex) ** 2)
    den = np.sum((y_true_ex - np.mean(y_train)) ** 2)
    return 1 - (num / den) if den != 0 else 0

def calculate_q2_f2(y_true_ex, y_pred_ex):
    # Mathematically equivalent to standard R2 on external set
    num = np.sum((y_true_ex - y_pred_ex) ** 2)
    den = np.sum((y_true_ex - np.mean(y_true_ex)) ** 2)
    return 1 - (num / den) if den != 0 else 0

def calculate_q2_f3(y_true_ex, y_pred_ex, y_train):
    n_ex = len(y_true_ex)
    n_tr = len(y_train)
    num = np.sum((y_true_ex - y_pred_ex) ** 2) / n_ex
    den = np.sum((y_train - np.mean(y_train)) ** 2) / n_tr
    return 1 - (num / den) if den != 0 else 0

def calculate_lof(y_true, y_pred, p):
    n = len(y_true)
    sse = np.sum((y_true - y_pred) ** 2)
    try:
        lof = (sse / n) / ((1 - (p + 1) / n) ** 2)
    except ZeroDivisionError:
        lof = np.nan
    return lof

def calculate_leverage(X_train, X_pred=None):
    X_train_1 = np.column_stack([np.ones(len(X_train)), X_train])
    try:
        pinv_mat = pinv(X_train_1.T @ X_train_1)
        hat_train = np.diag(X_train_1 @ pinv_mat @ X_train_1.T)
        if X_pred is not None:
            if len(X_pred) == 0:
                return hat_train, np.array([])
            X_pred_1 = np.column_stack([np.ones(len(X_pred)), X_pred])
            hat_pred = np.diag(X_pred_1 @ pinv_mat @ X_pred_1.T)
            return hat_train, hat_pred
        return hat_train
    except np.linalg.LinAlgError:
        nan_train = np.full(X_train.shape[0], np.nan)
        if X_pred is not None:
            nan_pred = np.full(X_pred.shape[0], np.nan) if len(X_pred) > 0 else np.array([])
            return nan_train, nan_pred
        return nan_train

def calculate_qsar_metrics(y_true, y_pred):
    metrics = {}
    metrics['R2_ext'] = r2_score(y_true, y_pred)
    metrics['RMSE_ext'] = float(np.sqrt(mean_squared_error(y_true, y_pred)))

    mean_true, var_true = np.mean(y_true), np.var(y_true)
    mean_pred, var_pred = np.mean(y_pred), np.var(y_pred)

    covar = np.mean((y_true - mean_true) * (y_pred - mean_pred))
    metrics['CCC'] = float((2 * covar) / (var_true + var_pred + (mean_true - mean_pred)**2)) if (var_true + var_pred + (mean_true - mean_pred)**2) != 0 else 0

    try:
        k, _, _, _ = np.linalg.lstsq(np.vstack([y_true, np.ones(len(y_true))]).T, y_pred, rcond=None)
        k_prime, _, _, _ = np.linalg.lstsq(np.vstack([y_pred, np.ones(len(y_pred))]).T, y_true, rcond=None)
        metrics['k'], metrics['k_prime'] = float(k[0]), float(k_prime[0])
    except:
        metrics['k'], metrics['k_prime'] = 1.0, 1.0

    r2_zero = 1 - (np.sum((y_true - y_pred)**2) / np.sum(y_true**2)) if np.sum(y_true**2) != 0 else 0
    r_prime2_zero = 1 - (np.sum((y_pred - y_true)**2) / np.sum(y_pred**2)) if np.sum(y_pred**2) != 0 else 0
    metrics['R2_zero_diff'] = float(abs(r2_zero - r_prime2_zero))

    return metrics

def check_golbraikh_tropsha(metrics):
    criteria = {
        "Q2_ext > 0.5": metrics.get('R2_ext', 0) > 0.5,
        "|R2_zero - R'2_zero| < 0.3": metrics.get('R2_zero_diff', 1) < 0.3,
        "0.85 <= k <= 1.15": 0.85 <= metrics.get('k', 0) <= 1.15,
        "0.85 <= k' <= 1.15": 0.85 <= metrics.get('k_prime', 0) <= 1.15
    }
    return all(criteria.values()), criteria


# ===================================================================
# INTERACTIVE CANVAS
# ===================================================================
class InteractiveCanvas(FigureCanvas):
    def __init__(self, width=6, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.points = []
        self.mpl_connect('button_press_event', self.on_double_click)

    def set_interactive_data(self, x_coords, y_coords, labels):
        self.points = [{'x': x, 'y': y, 'label': l} for x, y, l in zip(x_coords, y_coords, labels)]

    def clear_interactive_data(self):
        self.points = []

    def on_double_click(self, event):
        if event.dblclick and event.inaxes == self.ax and self.points:
            display_coords = self.ax.transData.transform([(pt['x'], pt['y']) for pt in self.points])
            event_coord = np.array([event.x, event.y])

            distances = np.sum((display_coords - event_coord)**2, axis=1)
            closest_idx = np.argmin(distances)

            if distances[closest_idx] < 400:
                pt = self.points[closest_idx]

                ann = self.ax.annotate(
                    f"{pt['label']}",
                    (pt['x'], pt['y']),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=plt.rcParams['font.size'] - 10,
                    bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.8),
                    arrowprops=dict(arrowstyle="->", connectionstyle="arc3")
                )

                ann.draggable(True)
                self.draw_idle()


# ===================================================================
# BACKGROUND WORKER THREAD
# ===================================================================
class RfaMarsWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    max_progress_signal = Signal(int)
    finished_signal = Signal(bool, str, object)

    def __init__(self, df, train_ids, test_ids, config):
        super().__init__()
        self.df = df
        self.train_ids = train_ids
        self.test_ids = test_ids
        self.config = config
        self.seed = 35

    def rfoa_feature_selection(self, X, y, pop_size, dimension, iterations):
        np.random.seed(self.seed)
        all_features = list(X.columns)
        dim = min(dimension, len(all_features))

        foxes = [list(np.random.choice(all_features, size=dim, replace=False)) for _ in range(pop_size)]

        def fitness(f):
            return r2_score(y, LinearRegression().fit(X[f], y).predict(X[f]))

        fit_scores = [fitness(f) for f in foxes]

        for t in range(iterations):
            sorted_idx = np.argsort(fit_scores)[::-1]
            foxes = [foxes[i] for i in sorted_idx]
            fit_scores = [fit_scores[i] for i in sorted_idx]

            for i in range(pop_size // 2, pop_size):
                new_f = foxes[i][:]
                idx = np.random.randint(dim)
                possible = [f for f in all_features if f not in new_f]
                if possible:
                    new_f[idx] = np.random.choice(possible)
                    new_fit = fitness(new_f)
                    if new_fit > fit_scores[i]:
                        foxes[i], fit_scores[i] = new_f, new_fit

            if t % max(1, (iterations // 5)) == 0 or t == iterations - 1:
                self.log_signal.emit(f"RFA Iter {t+1}/{iterations}: Best R2 = {fit_scores[0]:.4f}")

        return foxes[0]
        
    def get_model(self, model_name, alpha):
        if model_name == "MLR": return LinearRegression()
        if model_name == "Ridge": return Ridge(alpha=alpha, random_state=self.seed)
        if model_name == "Lasso": return Lasso(alpha=alpha, max_iter=10000, random_state=self.seed)
        if model_name == "PLS": return PLSRegression(n_components=2)
        if model_name == "Random Forest": return RandomForestRegressor(n_estimators=100, random_state=self.seed)
        if model_name == "SVR": return SVR(kernel='rbf')
        if model_name == "XGBoost" and XGB_AVAILABLE: return XGBRegressor(n_estimators=100, random_state=self.seed)
        return Lasso(alpha=alpha, max_iter=10000, random_state=self.seed)

    def run(self):
        try:
            # Parse Config
            pop_size = int(self.config.get('Pop Size', 30))
            iterations = int(self.config.get('Iterations', 30))
            dimension = int(self.config.get('Dimension', 5))
            n_knots = int(self.config.get('MARS Knots', 4))
            alpha = float(self.config.get('Lasso Alpha', 0.01))
            cv_folds = int(self.config.get('CV Folds', 5))
            model_type = self.config.get('Model Type', 'Lasso')
            use_mars = self.config.get('Use MARS', False)
            scale_desc = self.config.get('Scale Descriptors', True)
            runs = int(self.config.get('Y-Rand Runs', 50))

            # Prepare Data
            numeric_cols = self.df.select_dtypes(include=[np.number]).columns
            self.df[numeric_cols] = self.df[numeric_cols].fillna(self.df[numeric_cols].mean())

            train_indices = self.df[self.df.iloc[:, 0].isin(self.train_ids)].index
            test_indices = self.df[self.df.iloc[:, 0].isin(self.test_ids)].index

            Z = self.df.iloc[:, 0]
            Y = self.df.iloc[:, 1]
            X = self.df.iloc[:, 2:].select_dtypes(include=[np.number])

            X_train, X_test = X.loc[train_indices], X.loc[test_indices]
            y_train, y_test = Y.loc[train_indices], Y.loc[test_indices]
            z_train, z_test = Z.loc[train_indices], Z.loc[test_indices]

            self.log_signal.emit(f"Train shape: {X_train.shape}, Test shape: {X_test.shape}")

            # STAGE 1: RFA Selection
            self.log_signal.emit("\n--- STAGE 1: RFA Feature Selection ---")
            selected_features = self.rfoa_feature_selection(X_train, y_train, pop_size, dimension, iterations)
            self.log_signal.emit(f"> Best Descriptors: {list(selected_features)}")
            X_train_sel, X_test_sel = X_train[selected_features], X_test[selected_features]

            # STAGE 2: Transformations (MARS/Scaling)
            self.log_signal.emit(f"\n--- STAGE 2: Formatting Features ---")
            
            # The user requested that selecting a specific non-linear model overrides MARS. 
            # We will interpret "Use MARS" as valid primarily for Lasso/Ridge/MLR, but we respect the flag if set.
            # However, if it's overridden by the user, we skip it.
            if use_mars and model_type in ["Lasso", "Ridge", "MLR"]:
                spline = SplineTransformer(n_knots=n_knots, degree=1, include_bias=False).fit(X_train_sel)
                X_train_final = np.hstack([X_train_sel, spline.transform(X_train_sel)])
                X_test_final = np.hstack([X_test_sel, spline.transform(X_test_sel)])
                spline_names = spline.get_feature_names_out(selected_features)
                final_names = list(selected_features) + list(spline_names)
                self.log_signal.emit(f"Added {len(spline_names)} Spline/Knot features for MARS.")
            else:
                X_train_final, X_test_final = X_train_sel.to_numpy(), X_test_sel.to_numpy()
                final_names = list(selected_features)
                if use_mars:
                    self.log_signal.emit(f"MARS disabled for model '{model_type}'.")
            
            p = X_train_final.shape[1]
            n_tr = len(y_train)

            base_model = self.get_model(model_type, alpha)
            if scale_desc:
                pipeline = Pipeline([('scaler', StandardScaler()), ('model', base_model)])
                self.log_signal.emit("Applying Standard Scaling to features.")
            else:
                pipeline = Pipeline([('model', base_model)])

            # STAGE 3: Modeling & CV
            self.log_signal.emit(f"\n--- STAGE 3: Model Training ({model_type}) ---")
            
            pipeline.fit(X_train_final, y_train)
            y_train_pred = np.ravel(pipeline.predict(X_train_final))
            y_test_pred = np.ravel(pipeline.predict(X_test_final)) if len(y_test) > 0 else []

            # CV
            loo = LeaveOneOut()
            y_loo = np.ravel(cross_val_predict(pipeline, X_train_final, y_train, cv=loo))
            
            kfold = KFold(n_splits=cv_folds, shuffle=True, random_state=self.seed)
            y_lmo = np.ravel(cross_val_predict(pipeline, X_train_final, y_train, cv=kfold))
            
            r2_tr = r2_score(y_train, y_train_pred)
            r2_adj = 1 - (1 - r2_tr) * (n_tr - 1) / (n_tr - p - 1) if n_tr > p + 1 else np.nan
            rmse_tr = calculate_rmse(y_train, y_train_pred)
            ccc_tr = calculate_ccc(y_train, y_train_pred)
            lof = calculate_lof(y_train, y_train_pred, p)
            f_stat = (r2_tr / p) / ((1 - r2_tr) / (n_tr - p - 1)) if (n_tr > p + 1 and r2_tr != 1) else np.nan

            q2_loo = calculate_q2(y_train, y_loo)
            rmse_cv = calculate_rmse(y_train, y_loo)
            ccc_cv = calculate_ccc(y_train, y_loo)
            q2_lmo = calculate_q2(y_train, y_lmo)

            # STAGE 4: External Validation
            self.log_signal.emit("\n--- STAGE 4: External Validation ---")
            if len(y_test) > 0:
                r2_ex = r2_score(y_test, y_test_pred)
                rmse_ex = calculate_rmse(y_test, y_test_pred)
                ccc_ex = calculate_ccc(y_test, y_test_pred)
                q2_f1 = calculate_q2_f1(y_test, y_test_pred, y_train)
                q2_f2 = calculate_q2_f2(y_test, y_test_pred)
                q2_f3 = calculate_q2_f3(y_test, y_test_pred, y_train)
                gt_metrics = calculate_qsar_metrics(y_test, y_test_pred)
                passes, crits = check_golbraikh_tropsha(gt_metrics)
            else:
                r2_ex = rmse_ex = ccc_ex = q2_f1 = q2_f2 = q2_f3 = np.nan
                gt_metrics = {}
                passes, crits = False, {}

            self.log_signal.emit(f"Q2_ext: {r2_ex:.4f} | RMSE_ext: {rmse_ex:.4f} | CCC_ext: {ccc_ex:.4f}")

            if crits:
                self.log_signal.emit("-- Golbraikh-Tropsha Criteria --")
                for crit, val in crits.items():
                    self.log_signal.emit(f"{crit}: {'PASS' if val else 'FAIL'}")

            # Applicability Domain
            hat_train, hat_pred = calculate_leverage(X_train_final, X_test_final)
            warning_leverage = (3 * (p + 1)) / n_tr if n_tr > 0 else 0

            residuals_train = y_train - y_train_pred
            residuals_pred = y_test - y_test_pred if len(y_test) > 0 else []
            
            std_residuals_train = residuals_train / np.std(residuals_train) if np.std(residuals_train) != 0 else residuals_train
            std_residuals_pred = residuals_pred / np.std(residuals_train) if len(residuals_pred) > 0 and np.std(residuals_train) != 0 else residuals_pred

            # Y-Randomization
            self.log_signal.emit(f"\n--- STAGE 5: Y-Randomization ({runs} runs) ---")
            self.max_progress_signal.emit(runs)
            r2_list, q2_list, corr_list = [], [], []

            for i in range(runs):
                y_rand = np.random.permutation(y_train)
                pipeline.fit(X_train_final, y_rand)
                r2_rand = r2_score(y_rand, pipeline.predict(X_train_final))
                q2_rand = calculate_q2(y_rand, cross_val_predict(pipeline, X_train_final, y_rand, cv=loo))
                corr = np.abs(np.corrcoef(y_train, y_rand)[0, 1])

                r2_list.append(r2_rand)
                q2_list.append(q2_rand)
                corr_list.append(corr)
                self.progress_signal.emit(i + 1)

            r2_yscr = np.mean(r2_list)

            metrics_dict = {
                'R2tr': r2_tr, 'R2adj': r2_adj, 'LOF': lof, 'RMSEtr': rmse_tr, 'CCCtr': ccc_tr, 'F': f_stat,
                'R2cv': q2_loo, 'RMSEcv': rmse_cv, 'CCCcv': ccc_cv, 'Q2LMO': q2_lmo, 'R2Yscr': r2_yscr,
                'RMSEex': rmse_ex, 'R2ex': r2_ex, 'Q2-F1': q2_f1, 'Q2-F2': q2_f2, 'Q2-F3': q2_f3, 'CCCex': ccc_ex,
                'Warning_Lev': warning_leverage
            }

            train_df = pd.DataFrame({'ID': z_train, 'Set': 'Train', 'Actual': y_train, 'Predicted': y_train_pred, 'Residuals': residuals_train, 'Std_Residuals': std_residuals_train, 'Leverage': hat_train})
            for col in selected_features:
                train_df[col] = X_train_sel[col].values
                
            test_df = pd.DataFrame()
            if len(y_test) > 0:
                test_df = pd.DataFrame({'ID': z_test, 'Set': 'Test', 'Actual': y_test, 'Predicted': y_test_pred, 'Residuals': residuals_pred, 'Std_Residuals': std_residuals_pred, 'Leverage': hat_pred})
                for col in selected_features:
                    test_df[col] = X_test_sel[col].values
                    
            full_results = pd.concat([train_df, test_df], ignore_index=True)

            result_data = {
                'metrics': metrics_dict,
                'results_df': full_results,
                'plot_data': {
                    'y_train': y_train.values, 'y_train_pred': y_train_pred, 'sn_train': z_train.values,
                    'y_test': y_test.values if len(y_test) > 0 else [], 'y_test_pred': y_test_pred, 'sn_pred': z_test.values if len(y_test) > 0 else [],
                    'residuals_train': residuals_train.values, 'residuals_pred': residuals_pred.values if hasattr(residuals_pred, 'values') else residuals_pred,
                    'std_residuals_train': std_residuals_train.values if hasattr(std_residuals_train, 'values') else std_residuals_train, 
                    'std_residuals_pred': std_residuals_pred.values if hasattr(std_residuals_pred, 'values') else std_residuals_pred,
                    'hat_train': hat_train, 'hat_pred': hat_pred, 'warning_leverage': warning_leverage,
                    'yrand': (r2_list, q2_list, corr_list, r2_tr, q2_loo),
                    'heatmap_data': X_train_sel
                }
            }

            self.finished_signal.emit(True, "Analysis Complete", result_data)

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.finished_signal.emit(False, str(e), None)


# ===================================================================
# MAIN WIDGET UI
# ===================================================================
class QsarRfaWidget(PluginWidget):
    def __init__(self, plugin):
        super().__init__(plugin)
        self.dataset = None
        self.worker = None
        self.result_cache = None
        self.current_plot_data = None
        plt.rcParams.update({'font.size': 10})
        self.setup_ui()

    def setup_ui(self):
        self.widget = QWidget()
        main_layout = QVBoxLayout(self.widget)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # LEFT PANEL
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)

        # 1. Data Source
        left_layout.addWidget(QLabel("<b>1. Data Source:</b>"))
        data_layout = QHBoxLayout()
        self.btn_load = QPushButton("Load Dataset CSV")
        self.btn_load.clicked.connect(self.load_data)
        data_layout.addWidget(self.btn_load)
        self.lbl_file = QLabel("No file loaded")
        data_layout.addWidget(self.lbl_file)
        data_layout.addStretch()
        left_layout.addLayout(data_layout)

        # 2. Train/Test Split Dual-Table
        left_layout.addWidget(QLabel("<b>2. Define Train/Test Split:</b>"))
        split_layout = QHBoxLayout()
        self.tbl_test = QTableWidget()
        self.tbl_test.setColumnCount(1)
        self.tbl_test.setHorizontalHeaderLabels(["Test Set IDs"])
        self.tbl_test.horizontalHeader().setStretchLastSection(True)
        split_layout.addWidget(self.tbl_test)

        btn_move_layout = QVBoxLayout()
        self.btn_to_train = QPushButton(">> To Train >>")
        self.btn_to_train.clicked.connect(self.move_to_train)
        self.btn_to_test = QPushButton("<< To Test <<")
        self.btn_to_test.clicked.connect(self.move_to_test)
        btn_move_layout.addStretch()
        btn_move_layout.addWidget(self.btn_to_train)
        btn_move_layout.addWidget(self.btn_to_test)
        btn_move_layout.addStretch()
        split_layout.addLayout(btn_move_layout)

        self.tbl_train = QTableWidget()
        self.tbl_train.setColumnCount(1)
        self.tbl_train.setHorizontalHeaderLabels(["Training Set IDs"])
        self.tbl_train.horizontalHeader().setStretchLastSection(True)
        split_layout.addWidget(self.tbl_train)
        left_layout.addLayout(split_layout)

        auto_layout = QHBoxLayout()
        
        auto_layout.addWidget(QLabel("Train %:"))
        self.split_spin = QSpinBox()
        self.split_spin.setRange(1, 99)
        self.split_spin.setValue(80)
        auto_layout.addWidget(self.split_spin)

        self.btn_auto = QPushButton("Auto Split")
        self.btn_auto.clicked.connect(self.auto_split)
        auto_layout.addWidget(self.btn_auto)
        auto_layout.addStretch()
        left_layout.addLayout(auto_layout)

        # 3. Settings Grid
        left_layout.addWidget(QLabel("<b>3. Model & Selection Hyperparameters:</b>"))
        
        settings_layout = QHBoxLayout()
        
        # Standard Settings Table
        self.settings_tbl = QTableWidget()
        self.settings_tbl.setColumnCount(2)
        self.settings_tbl.setHorizontalHeaderLabels(["Setting", "Value"])
        self.settings_tbl.horizontalHeader().setStretchLastSection(True)
        self.settings_tbl.setMaximumHeight(150)

        default_settings = [
            ("Dimension", "5"), ("Pop Size", "30"), ("Iterations", "30"),
            ("MARS Knots", "4"), ("Lasso Alpha", "0.01"), ("CV Folds", "5")
        ]
        self.settings_tbl.setRowCount(len(default_settings))
        for i, (k, v) in enumerate(default_settings):
            k_item = QTableWidgetItem(k)
            k_item.setFlags(k_item.flags() & ~Qt.ItemIsEditable)
            self.settings_tbl.setItem(i, 0, k_item)
            self.settings_tbl.setItem(i, 1, QTableWidgetItem(v))
        
        settings_layout.addWidget(self.settings_tbl)
        
        # New Settings Group
        new_settings_group = QGroupBox("Model Options")
        new_settings_layout = QVBoxLayout()
        
        self.model_combo = QComboBox()
        models = ["Lasso", "MLR", "Ridge", "PLS", "Random Forest", "SVR"]
        if XGB_AVAILABLE: models.append("XGBoost")
        self.model_combo.addItems(models)
        
        new_settings_layout.addWidget(QLabel("Regressor:"))
        new_settings_layout.addWidget(self.model_combo)
        
        self.use_mars_check = QCheckBox("Apply MARS (Spline)")
        self.use_mars_check.setChecked(False)
        new_settings_layout.addWidget(self.use_mars_check)
        
        self.scale_check = QCheckBox("Scale Descriptors")
        self.scale_check.setChecked(True)
        new_settings_layout.addWidget(self.scale_check)
        
        self.yrand_spin = QSpinBox()
        self.yrand_spin.setRange(5, 500)
        self.yrand_spin.setValue(50)
        yrand_layout = QHBoxLayout()
        yrand_layout.addWidget(QLabel("Y-Rand Runs:"))
        yrand_layout.addWidget(self.yrand_spin)
        new_settings_layout.addLayout(yrand_layout)
        
        new_settings_group.setLayout(new_settings_layout)
        settings_layout.addWidget(new_settings_group)
        
        left_layout.addLayout(settings_layout)

        # 4. Execute & Logs
        self.btn_run = QPushButton("🚀 Run QSAR Analysis")
        self.btn_run.setStyleSheet("background-color: #2980b9; color: white; font-weight: bold; padding: 8px;")
        self.btn_run.clicked.connect(self.run_analysis)
        left_layout.addWidget(self.btn_run)
        
        self.progress = QProgressBar()
        self.progress.setVisible(False)
        left_layout.addWidget(self.progress)

        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        left_layout.addWidget(self.txt_log)

        export_layout = QHBoxLayout()
        self.btn_export = QPushButton("Export Results (CSV)")
        self.btn_export.clicked.connect(self.export_results)
        self.btn_export.setEnabled(False)
        export_layout.addWidget(self.btn_export)
        
        self.btn_export_plot = QPushButton("Export Current Plot")
        self.btn_export_plot.clicked.connect(self.export_plot)
        export_layout.addWidget(self.btn_export_plot)
        
        self.dpi_spin = QSpinBox()
        self.dpi_spin.setRange(72, 1200)
        self.dpi_spin.setValue(300)
        export_layout.addWidget(QLabel("DPI:"))
        export_layout.addWidget(self.dpi_spin)
        
        left_layout.addLayout(export_layout)
        
        splitter.addWidget(left_panel)
        
        # RIGHT PANEL (Graphs)
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)
        
        self.tabs = QTabWidget()
        
        self.canvas_exp_pred = InteractiveCanvas()
        self.canvas_residual = InteractiveCanvas()
        self.canvas_williams = InteractiveCanvas()
        self.canvas_yrand = InteractiveCanvas()
        self.canvas_heatmap = InteractiveCanvas()

        self.tabs.addTab(self.canvas_exp_pred, "Exp vs Pred")
        self.tabs.addTab(self.canvas_residual, "Residuals")
        self.tabs.addTab(self.canvas_williams, "Williams Plot")
        self.tabs.addTab(self.canvas_yrand, "Y-Randomization")
        self.tabs.addTab(self.canvas_heatmap, "Correlation Heatmap")
        
        font_layout = QHBoxLayout()
        font_layout.addWidget(QLabel("Global Font Size:"))
        self.font_spin = QSpinBox()
        self.font_spin.setRange(6, 24)
        self.font_spin.setValue(10)
        self.font_spin.valueChanged.connect(self.update_font_size)
        font_layout.addWidget(self.font_spin)
        font_layout.addStretch()
        
        right_layout.addLayout(font_layout)
        right_layout.addWidget(self.tabs)
        
        splitter.addWidget(right_panel)
        splitter.setSizes([450, 750])
        
        main_layout.addWidget(splitter)

    def load_data(self):
        path, _ = QFileDialog.getOpenFileName(self.widget, "Select Dataset", "", "CSV Files (*.csv)")
        if not path: return
        try:
            self.dataset = pd.read_csv(path, sep=None, engine='python')
            clean_name = path.replace("\\", "/").split("/")[-1]
            self.lbl_file.setText(clean_name)

            ids = sorted(list(self.dataset.iloc[:, 0]))
            self.tbl_test.setRowCount(len(ids))
            self.tbl_train.setRowCount(0)
            for i, val in enumerate(ids):
                self.tbl_test.setItem(i, 0, QTableWidgetItem(str(val)))

            self.txt_log.append(f"Loaded {len(ids)} compounds. Assuming Col 1=ID, Col 2=Activity.")
        except Exception as e:
            QMessageBox.critical(self.widget, "Load Error", str(e))

    def move_to_train(self):
        for item in self.tbl_test.selectedItems():
            row = self.tbl_train.rowCount()
            self.tbl_train.insertRow(row)
            self.tbl_train.setItem(row, 0, QTableWidgetItem(item.text()))
            self.tbl_test.removeRow(item.row())

    def move_to_test(self):
        for item in self.tbl_train.selectedItems():
            row = self.tbl_test.rowCount()
            self.tbl_test.insertRow(row)
            self.tbl_test.setItem(row, 0, QTableWidgetItem(item.text()))
            self.tbl_train.removeRow(item.row())

    def auto_split(self):
        if self.dataset is None: return
        all_ids = list(self.dataset.iloc[:, 0])
        np.random.shuffle(all_ids)
        
        train_pct = self.split_spin.value() / 100.0
        split_point = int(train_pct * len(all_ids))

        train_ids, test_ids = sorted(all_ids[:split_point]), sorted(all_ids[split_point:])

        self.tbl_train.setRowCount(len(train_ids))
        for i, val in enumerate(train_ids): self.tbl_train.setItem(i, 0, QTableWidgetItem(str(val)))

        self.tbl_test.setRowCount(len(test_ids))
        for i, val in enumerate(test_ids): self.tbl_test.setItem(i, 0, QTableWidgetItem(str(val)))

    def update_font_size(self, size):
        plt.rcParams.update({'font.size': size, 'axes.titlesize': size + 2, 'axes.labelsize': size})
        if self.current_plot_data:
            self.render_all_plots()

    def run_analysis(self):
        if self.dataset is None or self.tbl_train.rowCount() == 0:
            QMessageBox.warning(self.widget, "Warning", "Load data and assign training set first.")
            return

        self.btn_run.setEnabled(False)
        self.btn_export.setEnabled(False)
        self.txt_log.clear()

        config = {}
        for i in range(self.settings_tbl.rowCount()):
            config[self.settings_tbl.item(i, 0).text()] = self.settings_tbl.item(i, 1).text()
            
        config['Model Type'] = self.model_combo.currentText()
        config['Use MARS'] = self.use_mars_check.isChecked()
        config['Scale Descriptors'] = self.scale_check.isChecked()
        config['Y-Rand Runs'] = self.yrand_spin.value()

        train_ids = [self.tbl_train.item(i, 0).text() for i in range(self.tbl_train.rowCount())]
        test_ids = [self.tbl_test.item(i, 0).text() for i in range(self.tbl_test.rowCount())]

        self.progress.setValue(0)
        self.progress.setVisible(True)

        self.worker = RfaMarsWorker(self.dataset.copy(), train_ids, test_ids, config)
        self.worker.log_signal.connect(self.txt_log.append)
        self.worker.progress_signal.connect(self.progress.setValue)
        self.worker.max_progress_signal.connect(self.progress.setMaximum)
        self.worker.finished_signal.connect(self.handle_finish)
        self.worker.start()

    def handle_finish(self, success, msg, result_data):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        if success:
            self.result_cache = result_data['results_df']
            self.current_plot_data = result_data['plot_data']
            m = result_data['metrics']
            self.btn_export.setEnabled(True)
            
            metrics_str = f"""\n==============================
QSAR VALIDATION METRICS
==============================
--- Training set ---
R²tr       : {m['R2tr']:.4f}
R²adj      : {m['R2adj']:.4f}
RMSEtr     : {m['RMSEtr']:.4f}
CCCtr      : {m['CCCtr']:.4f}
LOF        : {m['LOF']:.4f}
F-stat     : {m['F']:.4f}

--- Cross-Validation ---
R²cv (Q²loo): {m['R2cv']:.4f}
RMSEcv     : {m['RMSEcv']:.4f}
CCCcv      : {m['CCCcv']:.4f}
Q²LMO      : {m['Q2LMO']:.4f}
R²Yscr     : {m['R2Yscr']:.4f}

--- External Prediction ---
R²ex       : {m['R2ex']:.4f}
RMSEex     : {m['RMSEex']:.4f}
CCCex      : {m['CCCex']:.4f}
Q²-F1      : {m['Q2-F1']:.4f}
Q²-F2      : {m['Q2-F2']:.4f}
Q²-F3      : {m['Q2-F3']:.4f}

Warning h* : {m['Warning_Lev']:.4f}
=============================="""
            self.txt_log.append(metrics_str)
            self.render_all_plots()
            
            QMessageBox.information(self.widget, "Success", "Analysis complete. View metrics and graphs.")
        else:
            QMessageBox.critical(self.widget, "Error", msg)

    def render_all_plots(self):
        pd_data = self.current_plot_data
        if not pd_data: return

        self.plot_exp_vs_pred(pd_data)
        self.plot_residuals(pd_data)
        self.plot_williams(pd_data)
        self.plot_yrandomization(*pd_data['yrand'])
        self.plot_heatmap(pd_data['heatmap_data'])

    def plot_exp_vs_pred(self, pd_data):
        ax = self.canvas_exp_pred.ax
        ax.clear()
        ax.scatter(pd_data['y_train'], pd_data['y_train_pred'], label='Training')
        if len(pd_data['y_test']) > 0:
            ax.scatter(pd_data['y_test'], pd_data['y_test_pred'], label='Prediction')
        minv = min(np.min(pd_data['y_train']), np.min(pd_data['y_test']) if len(pd_data['y_test']) > 0 else np.min(pd_data['y_train']))
        maxv = max(np.max(pd_data['y_train']), np.max(pd_data['y_test']) if len(pd_data['y_test']) > 0 else np.max(pd_data['y_train']))
        ax.plot([minv, maxv], [minv, maxv], 'k--')
        ax.set(xlabel="Experimental", ylabel="Predicted", title="Experimental vs Predicted")
        ax.legend()
        ax.grid(True)

        x_all = np.concatenate([pd_data['y_train'], pd_data['y_test']]) if len(pd_data['y_test']) > 0 else pd_data['y_train']
        y_all = np.concatenate([pd_data['y_train_pred'], pd_data['y_test_pred']]) if len(pd_data['y_test_pred']) > 0 else pd_data['y_train_pred']
        labels_all = np.concatenate([pd_data['sn_train'], pd_data['sn_pred']]) if len(pd_data['sn_pred']) > 0 else pd_data['sn_train']
        self.canvas_exp_pred.set_interactive_data(x_all, y_all, labels_all)
        self.canvas_exp_pred.draw()

    def plot_residuals(self, pd_data):
        ax = self.canvas_residual.ax
        ax.clear()
        ax.scatter(pd_data['y_train_pred'], pd_data['residuals_train'], label='Training')
        if len(pd_data['residuals_pred']) > 0:
            ax.scatter(pd_data['y_test_pred'], pd_data['residuals_pred'], label='Prediction')
        ax.axhline(0, color='k', linestyle='--')
        ax.set(xlabel="Predicted", ylabel="Residuals", title="Residual Plot")
        ax.legend()
        ax.grid(True)

        x_all = np.concatenate([pd_data['y_train_pred'], pd_data['y_test_pred']]) if len(pd_data['y_test_pred']) > 0 else pd_data['y_train_pred']
        y_all = np.concatenate([pd_data['residuals_train'], pd_data['residuals_pred']]) if len(pd_data['residuals_pred']) > 0 else pd_data['residuals_train']
        labels_all = np.concatenate([pd_data['sn_train'], pd_data['sn_pred']]) if len(pd_data['sn_pred']) > 0 else pd_data['sn_train']
        self.canvas_residual.set_interactive_data(x_all, y_all, labels_all)
        self.canvas_residual.draw()

    def plot_williams(self, pd_data):
        ax = self.canvas_williams.ax
        ax.clear()
        ax.scatter(pd_data['hat_train'], pd_data['std_residuals_train'], label='Training')
        if len(pd_data['hat_pred']) > 0:
            ax.scatter(pd_data['hat_pred'], pd_data['std_residuals_pred'], label='Prediction')

        ax.axhline(3, color='r', linestyle='--')
        ax.axhline(-3, color='r', linestyle='--')
        ax.axvline(pd_data['warning_leverage'], color='r', linestyle='--')
        ax.set(xlabel="Leverage", ylabel="Standardized Residual", title="Williams Plot")
        ax.legend()
        ax.grid(True)

        x_all = np.concatenate([pd_data['hat_train'], pd_data['hat_pred']]) if len(pd_data['hat_pred']) > 0 else pd_data['hat_train']
        y_all = np.concatenate([pd_data['std_residuals_train'], pd_data['std_residuals_pred']]) if len(pd_data['std_residuals_pred']) > 0 else pd_data['std_residuals_train']
        labels_all = np.concatenate([pd_data['sn_train'], pd_data['sn_pred']]) if len(pd_data['sn_pred']) > 0 else pd_data['sn_train']
        self.canvas_williams.set_interactive_data(x_all, y_all, labels_all)
        self.canvas_williams.draw()

    def plot_yrandomization(self, r2_list, q2_list, corr_list, orig_r2, orig_q2):
        ax = self.canvas_yrand.ax
        ax.clear()
        ax.scatter(corr_list, r2_list, label='Random R²', alpha=0.6)
        ax.scatter(corr_list, q2_list, label='Random Q²', alpha=0.6)
        ax.scatter([1.0], [orig_r2], color='blue', label='Original R²', s=100, marker='*')
        ax.scatter([1.0], [orig_q2], color='orange', label='Original Q²', s=100, marker='*')

        if len(corr_list) > 1:
            try:
                x_vals = np.linspace(0, 1, 50)
                ax.plot(x_vals, np.poly1d(np.polyfit(corr_list, r2_list, 1))(x_vals), 'b--', alpha=0.5)
                ax.plot(x_vals, np.poly1d(np.polyfit(corr_list, q2_list, 1))(x_vals), 'r--', alpha=0.5)
            except np.linalg.LinAlgError: pass

        ax.set(xlabel="Correlation coefficient (r)", ylabel="R² / Q²", title="Y-Randomization")
        ax.legend()
        ax.grid(True)
        self.canvas_yrand.clear_interactive_data()
        self.canvas_yrand.draw()

    def plot_heatmap(self, df_descriptors):
        self.canvas_heatmap.fig.clf()
        ax = self.canvas_heatmap.fig.add_subplot(111)
        self.canvas_heatmap.ax = ax
        
        corr = df_descriptors.corr()
        cax = ax.imshow(corr, cmap="coolwarm", aspect='auto', vmin=-1, vmax=1)
        self.canvas_heatmap.fig.colorbar(cax, ax=ax)
        
        ax.set_xticks(np.arange(len(corr.columns)))
        ax.set_yticks(np.arange(len(corr.columns)))
        ax.set_xticklabels(corr.columns, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(corr.columns, fontsize=8)
        
        ax.set_title("Descriptor Correlation Heatmap")
        self.canvas_heatmap.fig.tight_layout()
        self.canvas_heatmap.clear_interactive_data()
        self.canvas_heatmap.draw()

    def export_results(self):
        if self.result_cache is None: return
        path, _ = QFileDialog.getSaveFileName(self.widget, "Save CSV Results", "QSAR_Results.csv", "CSV (*.csv)")
        if path:
            try:
                self.result_cache.to_csv(path, index=False)
                QMessageBox.information(self.widget, "Export Success", f"Results successfully saved to {path}!")
            except Exception as e:
                QMessageBox.critical(self.widget, "Export Error", str(e))

    def export_plot(self):
        current_widget = self.tabs.currentWidget()
        if not isinstance(current_widget, FigureCanvas): return

        path, _ = QFileDialog.getSaveFileName(self.widget, "Save Plot", "plot.png", "Images (*.png *.jpg *.jpeg)")
        if path:
            dpi_val = self.dpi_spin.value()
            try:
                current_widget.fig.savefig(path, dpi=dpi_val, bbox_inches='tight')
                QMessageBox.information(self.widget, "Success", f"Plot saved successfully at {dpi_val} DPI.")
            except Exception as e:
                QMessageBox.critical(self.widget, "Export Error", f"Could not save image:\n{str(e)}")

# ===================================================================
# PLUGIN WRAPPER
# ===================================================================
class QsarRfaPlugin(BasePlugin):
    def __init__(self):
        super().__init__(PluginInfo(
            name="RFA + MARS QSAR",
            version="1.1.0",
            description="Hybrid Red Fox Algorithm feature selection and non-linear modeling with comprehensive validation",
            author="SMILES Team",
            plugin_type=PluginType.ANALYSIS,
            dependencies=[]
        ))
        self.widget = None

    def get_info(self) -> PluginInfo:
        return self.info

    def create_widget(self) -> 'QsarRfaWidget':
        if self.widget is None:
            self.widget = QsarRfaWidget(self)
        return self.widget

    def initialize(self):
        self.logger.info("RFA+MARS QSAR plugin initialized")
        return True

    def cleanup(self):
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info("RFA+MARS QSAR plugin cleaned up")
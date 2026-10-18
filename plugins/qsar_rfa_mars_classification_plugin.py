"""
Hybrid RFA + MARS Classification QSAR Plugin

Classification-focused adaptation of the regression RFA+MARS plugin.

Features
--------
- User-selectable Binary / Ternary / Multiclass analysis modes
- Red Fox Algorithm (RFA) descriptor selection using stratified CV fitness
- Optional SplineTransformer (MARS-like) feature expansion
- Logistic Regression, Ridge Classifier, Random Forest, SVC and optional XGBoost
- Stratified train/test auto-splitting
- Stratified K-fold cross-validation
- Accuracy, Balanced Accuracy, MCC, Cohen's kappa, Precision, Recall, F1
- ROC-AUC and PR-AUC when probability/decision scores are available
- Confusion matrix and ROC/PR visualizations
- Classification applicability-domain leverage analysis
- Class-label Y-randomization / permutation validation
- Per-compound predicted class, confidence, probabilities and leverage export

No dependency beyond packages already used by the original plugin is introduced.
"""

import math
import os
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_val_predict, train_test_split
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, matthews_corrcoef,
    cohen_kappa_score, precision_score, recall_score, f1_score,
    confusion_matrix, roc_auc_score, average_precision_score,
    roc_curve, precision_recall_curve
)
from sklearn.preprocessing import SplineTransformer, StandardScaler, LabelEncoder, label_binarize
from sklearn.pipeline import Pipeline
from scipy.linalg import pinv

try:
    from xgboost import XGBClassifier
    XGB_AVAILABLE = True
except ImportError:
    XGB_AVAILABLE = False

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# Strictly using ONLY the allowed imports from qt_compat
from src.shared.qt_compat import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QTextEdit, QProgressBar,
    QFileDialog, QMessageBox, QComboBox, Qt, QThread, Signal,
    QTabWidget, QCheckBox, QSpinBox, QSplitter, QGroupBox,
    QDialog
)
from src.plugins.base_plugin import BasePlugin, PluginWidget
from src.plugins.plugin_types import PluginInfo, PluginType


# ===================================================================
# CLASSIFICATION & AD HELPER FUNCTIONS
# ===================================================================

def calculate_leverage(X_train, X_pred=None):
    """Classical descriptor-space leverage; useful as an X-space AD diagnostic."""
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


def safe_mcc(y_true, y_pred):
    try:
        return float(matthews_corrcoef(y_true, y_pred))
    except Exception:
        return np.nan


def classification_metrics(y_true, y_pred, y_score=None, n_classes=2):
    """Return robust metrics for binary, ternary and general multiclass tasks."""
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    metrics = {
        'Accuracy': float(accuracy_score(y_true, y_pred)),
        'Balanced_Accuracy': float(balanced_accuracy_score(y_true, y_pred)),
        'MCC': safe_mcc(y_true, y_pred),
        'Kappa': float(cohen_kappa_score(y_true, y_pred)),
        'Precision_macro': float(precision_score(y_true, y_pred, average='macro', zero_division=0)),
        'Recall_macro': float(recall_score(y_true, y_pred, average='macro', zero_division=0)),
        'F1_macro': float(f1_score(y_true, y_pred, average='macro', zero_division=0)),
        'F1_weighted': float(f1_score(y_true, y_pred, average='weighted', zero_division=0)),
    }

    metrics['ROC_AUC'] = np.nan
    metrics['PR_AUC'] = np.nan

    if y_score is not None:
        try:
            score = np.asarray(y_score)
            if n_classes == 2:
                pos_score = score[:, 1] if score.ndim == 2 else score
                metrics['ROC_AUC'] = float(roc_auc_score(y_true, pos_score))
                metrics['PR_AUC'] = float(average_precision_score(y_true, pos_score))
            else:
                y_bin = label_binarize(y_true, classes=np.arange(n_classes))
                metrics['ROC_AUC'] = float(roc_auc_score(
                    y_true, score, multi_class='ovr', average='macro'
                ))
                metrics['PR_AUC'] = float(average_precision_score(
                    y_bin, score, average='macro'
                ))
        except Exception:
            pass

    return metrics


def get_score_matrix(estimator, X, n_classes):
    """Get probability-like class scores without requiring calibration dependencies."""
    if hasattr(estimator, 'predict_proba'):
        try:
            return np.asarray(estimator.predict_proba(X), dtype=float), 'probability'
        except Exception:
            pass

    if hasattr(estimator, 'decision_function'):
        try:
            raw = np.asarray(estimator.decision_function(X), dtype=float)
            if n_classes == 2:
                # Monotonic logistic mapping for convenient 0..1 display only.
                pos = 1.0 / (1.0 + np.exp(-np.clip(raw, -50, 50)))
                return np.column_stack([1.0 - pos, pos]), 'decision-derived'
            raw = raw - np.max(raw, axis=1, keepdims=True)
            exp = np.exp(np.clip(raw, -50, 50))
            probs = exp / np.sum(exp, axis=1, keepdims=True)
            return probs, 'decision-derived'
        except Exception:
            pass

    return None, 'unavailable'


def minimum_class_count(y):
    _, counts = np.unique(np.asarray(y), return_counts=True)
    return int(counts.min()) if len(counts) else 0


def validate_class_mode(n_classes, mode):
    if mode == 'Binary' and n_classes != 2:
        raise ValueError(f"Binary analysis requires exactly 2 classes; found {n_classes}.")
    if mode == 'Ternary' and n_classes != 3:
        raise ValueError(f"Ternary analysis requires exactly 3 classes; found {n_classes}.")
    if mode == 'Multiclass' and n_classes < 3:
        raise ValueError(f"Multiclass analysis requires at least 3 classes; found {n_classes}. Use Binary for 2 classes.")


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
        self.points = [
            {'x': x, 'y': y, 'label': l}
            for x, y, l in zip(x_coords, y_coords, labels)
        ]

    def clear_interactive_data(self):
        self.points = []

    def on_double_click(self, event):
        if event.dblclick and event.inaxes == self.ax and self.points:
            display_coords = self.ax.transData.transform(
                [(pt['x'], pt['y']) for pt in self.points]
            )
            event_coord = np.array([event.x, event.y])
            distances = np.sum((display_coords - event_coord) ** 2, axis=1)
            closest_idx = np.argmin(distances)

            if distances[closest_idx] < 400:
                pt = self.points[closest_idx]
                ann = self.ax.annotate(
                    f"{pt['label']}", (pt['x'], pt['y']),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=max(6, plt.rcParams['font.size'] - 2),
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
        self.train_ids = [str(v) for v in train_ids]
        self.test_ids = [str(v) for v in test_ids]
        self.config = config
        self.seed = int(config.get('Random Seed', 35))

    def get_model(self, model_name, class_weight=True):
        cw = 'balanced' if class_weight else None
        if model_name == 'Logistic Regression':
            return LogisticRegression(
                max_iter=5000, class_weight=cw, random_state=self.seed
            )
        if model_name == 'Ridge Classifier':
            return RidgeClassifier(class_weight=cw)
        if model_name == 'Random Forest':
            return RandomForestClassifier(
                n_estimators=300, class_weight=cw, random_state=self.seed,
                n_jobs=-1
            )
        if model_name == 'SVC':
            return SVC(
                kernel='rbf', probability=True, class_weight=cw,
                random_state=self.seed
            )
        if model_name == 'XGBoost' and XGB_AVAILABLE:
            # sample_weight is not introduced here; multiclass is inferred from y.
            return XGBClassifier(
                n_estimators=300, random_state=self.seed,
                eval_metric='logloss', n_jobs=-1
            )
        return LogisticRegression(
            max_iter=5000, class_weight=cw, random_state=self.seed
        )

    def make_pipeline(self, model_name, scale_desc, class_weight=True):
        model = self.get_model(model_name, class_weight)
        if scale_desc:
            return Pipeline([('scaler', StandardScaler()), ('model', model)])
        return Pipeline([('model', model)])

    def _fitness_cv_splits(self, y):
        min_count = minimum_class_count(y)
        n_splits = min(3, min_count)
        return n_splits if n_splits >= 2 else 0

    def rfoa_feature_selection(self, X, y, pop_size, dimension, iterations, model_name, scale_desc, class_weight, log_updates=True):
        """RFA fitness uses stratified CV balanced accuracy, not resubstitution accuracy."""
        rng = np.random.RandomState(self.seed)
        all_features = list(X.columns)
        dim = min(max(1, dimension), len(all_features))
        if dim == 0:
            raise ValueError('No numeric descriptor columns were found.')

        foxes = [list(rng.choice(all_features, size=dim, replace=False)) for _ in range(pop_size)]
        n_splits = self._fitness_cv_splits(y)
        if n_splits < 2:
            raise ValueError('Each training class must contain at least 2 samples for feature-selection validation.')
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=self.seed)

        def fitness(features):
            pipe = self.make_pipeline(model_name, scale_desc, class_weight)
            scores = []
            for tr_idx, va_idx in cv.split(X[features], y):
                try:
                    pipe.fit(X.iloc[tr_idx][features], np.asarray(y)[tr_idx])
                    pred = pipe.predict(X.iloc[va_idx][features])
                    scores.append(balanced_accuracy_score(np.asarray(y)[va_idx], pred))
                except Exception:
                    return -np.inf
            return float(np.mean(scores)) if scores else -np.inf

        fit_scores = [fitness(f) for f in foxes]

        for t in range(iterations):
            sorted_idx = np.argsort(fit_scores)[::-1]
            foxes = [foxes[i] for i in sorted_idx]
            fit_scores = [fit_scores[i] for i in sorted_idx]

            for i in range(pop_size // 2, pop_size):
                new_f = foxes[i][:]
                idx = rng.randint(dim)
                possible = [f for f in all_features if f not in new_f]
                if possible:
                    new_f[idx] = rng.choice(possible)
                    new_fit = fitness(new_f)
                    if new_fit > fit_scores[i]:
                        foxes[i], fit_scores[i] = new_f, new_fit

            if log_updates and (t % max(1, iterations // 5) == 0 or t == iterations - 1):
                self.log_signal.emit(
                    f"RFA Iter {t + 1}/{iterations}: Best CV Balanced Accuracy = {fit_scores[0]:.4f}"
                )

        return foxes[0], float(fit_scores[0])

    def run(self):
        try:
            pop_size = int(self.config.get('Pop Size', 30))
            iterations = int(self.config.get('Iterations', 30))
            dimension = int(self.config.get('Dimension', 5))
            n_knots = int(self.config.get('MARS Knots', 4))
            cv_folds = int(self.config.get('CV Folds', 5))
            model_type = self.config.get('Model Type', 'Logistic Regression')
            class_mode = self.config.get('Class Mode', 'Binary')
            use_mars = bool(self.config.get('Use MARS', False))
            scale_desc = bool(self.config.get('Scale Descriptors', True))
            class_weight = bool(self.config.get('Class Weight', True))
            nested_cv = bool(self.config.get('Nested CV', False))
            runs = int(self.config.get('Y-Rand Runs', 50))

            # ---------------------------
            # Prepare and validate data
            # ---------------------------
            if self.df.shape[1] < 3:
                raise ValueError('Dataset must contain ID, class target, and at least one descriptor column.')

            # Preserve target values; impute descriptors only.
            Z = self.df.iloc[:, 0].astype(str)
            Y_raw = self.df.iloc[:, 1]
            if Y_raw.isna().any():
                raise ValueError('Target/class column contains missing values. Please resolve them before analysis.')

            X = self.df.iloc[:, 2:].select_dtypes(include=[np.number]).copy()
            if X.shape[1] == 0:
                raise ValueError('No numeric descriptor columns found from column 3 onward.')
            X = X.replace([np.inf, -np.inf], np.nan)

            train_mask = Z.isin(self.train_ids)
            test_mask = Z.isin(self.test_ids)
            if train_mask.sum() == 0:
                raise ValueError('No training IDs matched the first dataset column.')
            if np.any(train_mask & test_mask):
                raise ValueError('Some IDs are assigned to both training and test sets.')

            # Training-derived imputation avoids information leakage.
            train_means = X.loc[train_mask].mean()
            X = X.fillna(train_means).fillna(0.0)

            encoder = LabelEncoder()
            y_all = encoder.fit_transform(Y_raw.astype(str))
            class_labels = list(encoder.classes_)
            n_classes = len(class_labels)
            validate_class_mode(n_classes, class_mode)

            y_series = pd.Series(y_all, index=self.df.index)
            X_train, X_test = X.loc[train_mask], X.loc[test_mask]
            y_train = y_series.loc[train_mask].to_numpy()
            y_test = y_series.loc[test_mask].to_numpy()
            z_train, z_test = Z.loc[train_mask], Z.loc[test_mask]

            missing_train_classes = sorted(set(range(n_classes)) - set(np.unique(y_train)))
            if missing_train_classes:
                missing_names = [class_labels[i] for i in missing_train_classes]
                raise ValueError(f'Training set is missing class(es): {missing_names}.')

            if len(y_test) > 0:
                missing_test_classes = sorted(set(range(n_classes)) - set(np.unique(y_test)))
                if missing_test_classes:
                    self.log_signal.emit(
                        'Warning: external test set does not contain every class; some external macro/AUC metrics may be unavailable.'
                    )

            class_counts = pd.Series(y_train).value_counts().sort_index()
            self.log_signal.emit(f"Analysis mode: {class_mode}; classes = {class_labels}")
            self.log_signal.emit(f"Train shape: {X_train.shape}, Test shape: {X_test.shape}")
            self.log_signal.emit(
                'Training class counts: ' + ', '.join(
                    f"{class_labels[i]}={int(class_counts.get(i, 0))}" for i in range(n_classes)
                )
            )

            min_count = minimum_class_count(y_train)
            if min_count < 2:
                raise ValueError('Every training class needs at least 2 samples for stratified validation.')
            actual_cv_folds = min(cv_folds, min_count)
            if actual_cv_folds < cv_folds:
                self.log_signal.emit(
                    f"CV folds reduced from {cv_folds} to {actual_cv_folds} because the smallest class has {min_count} samples."
                )

            # ---------------------------
            # STAGE 1: RFA selection
            # ---------------------------
            self.log_signal.emit('\n--- STAGE 1: RFA Feature Selection ---')
            selected_features, rfa_fitness = self.rfoa_feature_selection(
                X_train, y_train, pop_size, dimension, iterations,
                model_type, scale_desc, class_weight
            )
            self.log_signal.emit(f"> Best Descriptors: {list(selected_features)}")
            self.log_signal.emit(f"> RFA CV Balanced Accuracy: {rfa_fitness:.4f}")
            X_train_sel = X_train[selected_features]
            X_test_sel = X_test[selected_features]

            # ---------------------------
            # STAGE 2: MARS-like expansion
            # ---------------------------
            self.log_signal.emit('\n--- STAGE 2: Formatting Features ---')
            if use_mars:
                spline = SplineTransformer(
                    n_knots=n_knots, degree=1, include_bias=False
                ).fit(X_train_sel)
                X_train_final = np.hstack([X_train_sel, spline.transform(X_train_sel)])
                X_test_final = np.hstack([X_test_sel, spline.transform(X_test_sel)]) if len(X_test_sel) else np.empty((0, X_train_sel.shape[1] + spline.transform(X_train_sel.iloc[:1]).shape[1]))
                spline_names = spline.get_feature_names_out(selected_features)
                final_names = list(selected_features) + list(spline_names)
                self.log_signal.emit(f"Added {len(spline_names)} spline features.")
            else:
                X_train_final = X_train_sel.to_numpy()
                X_test_final = X_test_sel.to_numpy()
                final_names = list(selected_features)

            p = X_train_final.shape[1]
            n_tr = len(y_train)
            pipeline = self.make_pipeline(model_type, scale_desc, class_weight)
            if scale_desc:
                self.log_signal.emit('Applying Standard Scaling to model features.')

            # ---------------------------
            # STAGE 3: Fit + stratified CV
            # ---------------------------
            self.log_signal.emit(f"\n--- STAGE 3: Model Training ({model_type}) ---")
            pipeline.fit(X_train_final, y_train)
            y_train_pred = np.asarray(pipeline.predict(X_train_final), dtype=int)
            y_test_pred = np.asarray(pipeline.predict(X_test_final), dtype=int) if len(y_test) else np.array([], dtype=int)

            score_train, score_kind = get_score_matrix(pipeline, X_train_final, n_classes)
            score_test, _ = get_score_matrix(pipeline, X_test_final, n_classes) if len(y_test) else (None, score_kind)
            self.log_signal.emit(f"Class-score source: {score_kind}")

            skf = StratifiedKFold(
                n_splits=actual_cv_folds, shuffle=True, random_state=self.seed
            )

            # Optional true nested CV: RFA feature selection is repeated inside every
            # outer fold, preventing descriptor-selection leakage into CV estimates.
            if nested_cv:
                self.log_signal.emit(
                    'Nested CV enabled: descriptor selection is repeated inside each outer fold (slower, less biased).'
                )
                y_cv = np.empty(len(y_train), dtype=int)
                cv_score = np.full((len(y_train), n_classes), np.nan, dtype=float)
                nested_score_available = True
                nested_score_kinds = []

                for fold_no, (outer_tr, outer_va) in enumerate(skf.split(X_train, y_train), start=1):
                    self.log_signal.emit(f'Nested CV outer fold {fold_no}/{actual_cv_folds}...')
                    X_ot = X_train.iloc[outer_tr]
                    X_ov = X_train.iloc[outer_va]
                    y_ot = y_train[outer_tr]

                    fold_features, _ = self.rfoa_feature_selection(
                        X_ot, y_ot, pop_size, dimension, iterations,
                        model_type, scale_desc, class_weight, log_updates=False
                    )
                    X_ot_sel = X_ot[fold_features]
                    X_ov_sel = X_ov[fold_features]

                    if use_mars:
                        fold_spline = SplineTransformer(
                            n_knots=n_knots, degree=1, include_bias=False
                        ).fit(X_ot_sel)
                        X_ot_final = np.hstack([X_ot_sel, fold_spline.transform(X_ot_sel)])
                        X_ov_final = np.hstack([X_ov_sel, fold_spline.transform(X_ov_sel)])
                    else:
                        X_ot_final = X_ot_sel.to_numpy()
                        X_ov_final = X_ov_sel.to_numpy()

                    fold_pipe = self.make_pipeline(model_type, scale_desc, class_weight)
                    fold_pipe.fit(X_ot_final, y_ot)
                    y_cv[outer_va] = np.asarray(fold_pipe.predict(X_ov_final), dtype=int)
                    fold_score, fold_kind = get_score_matrix(fold_pipe, X_ov_final, n_classes)
                    nested_score_kinds.append(fold_kind)
                    if fold_score is None or np.asarray(fold_score).shape[1] != n_classes:
                        nested_score_available = False
                    else:
                        cv_score[outer_va, :] = np.asarray(fold_score)

                if not nested_score_available or np.isnan(cv_score).any():
                    cv_score = None
                    cv_score_kind = 'unavailable'
                else:
                    cv_score_kind = 'probability' if all(k == 'probability' for k in nested_score_kinds) else 'decision-derived'
            else:
                self.log_signal.emit(
                    'Post-selection stratified CV enabled. For a less biased estimate, enable Nested CV (slower).'
                )
                y_cv = np.asarray(cross_val_predict(
                    pipeline, X_train_final, y_train, cv=skf, method='predict'
                ), dtype=int)

                # Out-of-fold class probabilities/scores when supported.
                cv_score = None
                cv_score_kind = 'unavailable'
                if hasattr(pipeline, 'predict_proba'):
                    try:
                        cv_score = np.asarray(cross_val_predict(
                            pipeline, X_train_final, y_train, cv=skf, method='predict_proba'
                        ), dtype=float)
                        cv_score_kind = 'probability'
                    except Exception:
                        cv_score = None
                if cv_score is None and hasattr(pipeline, 'decision_function'):
                    try:
                        raw = np.asarray(cross_val_predict(
                            pipeline, X_train_final, y_train, cv=skf, method='decision_function'
                        ), dtype=float)
                        if n_classes == 2:
                            pos = 1.0 / (1.0 + np.exp(-np.clip(raw, -50, 50)))
                            cv_score = np.column_stack([1 - pos, pos])
                        else:
                            raw = raw - np.max(raw, axis=1, keepdims=True)
                            exp = np.exp(np.clip(raw, -50, 50))
                            cv_score = exp / np.sum(exp, axis=1, keepdims=True)
                        cv_score_kind = 'decision-derived'
                    except Exception:
                        cv_score = None

            train_metrics = classification_metrics(y_train, y_train_pred, score_train, n_classes)
            cv_metrics = classification_metrics(y_train, y_cv, cv_score, n_classes)
            test_metrics = classification_metrics(y_test, y_test_pred, score_test, n_classes) if len(y_test) else {
                k: np.nan for k in train_metrics.keys()
            }

            self.log_signal.emit(
                f"Training: BalAcc={train_metrics['Balanced_Accuracy']:.4f}, MCC={train_metrics['MCC']:.4f}, F1macro={train_metrics['F1_macro']:.4f}"
            )
            self.log_signal.emit(
                f"CV:       BalAcc={cv_metrics['Balanced_Accuracy']:.4f}, MCC={cv_metrics['MCC']:.4f}, F1macro={cv_metrics['F1_macro']:.4f}"
            )
            if len(y_test):
                self.log_signal.emit(
                    f"External: BalAcc={test_metrics['Balanced_Accuracy']:.4f}, MCC={test_metrics['MCC']:.4f}, F1macro={test_metrics['F1_macro']:.4f}"
                )

            # ---------------------------
            # STAGE 4: Applicability Domain
            # ---------------------------
            self.log_signal.emit('\n--- STAGE 4: Applicability Domain ---')
            # AD is calculated in the selected original-descriptor space.
            # This keeps h* interpretable and avoids spline-basis inflation.
            X_ad_train = X_train_sel.to_numpy()
            X_ad_test = X_test_sel.to_numpy()
            hat_train, hat_test = calculate_leverage(X_ad_train, X_ad_test)
            ad_p = X_ad_train.shape[1]
            warning_leverage = (3.0 * (ad_p + 1)) / n_tr if n_tr else np.nan
            train_ad = hat_train <= warning_leverage
            test_ad = hat_test <= warning_leverage if len(hat_test) else np.array([], dtype=bool)
            self.log_signal.emit(
                f"Leverage warning h*={warning_leverage:.4f}; train outside AD={int(np.sum(~train_ad))}"
            )
            if len(test_ad):
                self.log_signal.emit(f"External test outside AD={int(np.sum(~test_ad))}")

            # ---------------------------
            # STAGE 5: Y-randomization
            # ---------------------------
            self.log_signal.emit(f"\n--- STAGE 5: Class-label Y-Randomization ({runs} runs) ---")
            self.max_progress_signal.emit(runs)
            rng = np.random.RandomState(self.seed)
            yrand_balacc = []
            yrand_mcc = []

            # Use same stratified fold count, recreated for each shuffled target.
            for i in range(runs):
                y_rand = rng.permutation(y_train)
                # Guard against pathological folds by stratifying the permuted labels.
                cv_rand = StratifiedKFold(
                    n_splits=min(actual_cv_folds, minimum_class_count(y_rand)),
                    shuffle=True, random_state=self.seed + i + 1
                )
                rand_pred = cross_val_predict(
                    pipeline, X_train_final, y_rand, cv=cv_rand, method='predict'
                )
                yrand_balacc.append(float(balanced_accuracy_score(y_rand, rand_pred)))
                yrand_mcc.append(safe_mcc(y_rand, rand_pred))
                self.progress_signal.emit(i + 1)

            yrand_mean_balacc = float(np.nanmean(yrand_balacc)) if yrand_balacc else np.nan
            yrand_mean_mcc = float(np.nanmean(yrand_mcc)) if yrand_mcc else np.nan
            yrand_p_balacc = float((1 + np.sum(np.asarray(yrand_balacc) >= cv_metrics['Balanced_Accuracy'])) / (runs + 1)) if runs else np.nan
            yrand_p_mcc = float((1 + np.sum(np.asarray(yrand_mcc) >= cv_metrics['MCC'])) / (runs + 1)) if runs else np.nan

            # ---------------------------
            # Result tables
            # ---------------------------
            def build_result_frame(ids, y_true, y_pred, scores, leverage, in_ad, set_name):
                if len(y_true) == 0:
                    return pd.DataFrame()
                actual_labels = encoder.inverse_transform(np.asarray(y_true, dtype=int))
                pred_labels = encoder.inverse_transform(np.asarray(y_pred, dtype=int))
                out = pd.DataFrame({
                    'ID': np.asarray(ids),
                    'Set': set_name,
                    'Actual': actual_labels,
                    'Predicted': pred_labels,
                    'Correct': np.asarray(y_true) == np.asarray(y_pred),
                    'Leverage': leverage,
                    'In_AD': in_ad
                })
                if scores is not None:
                    scores_arr = np.asarray(scores)
                    out['Confidence'] = np.max(scores_arr, axis=1)
                    for ci, cname in enumerate(class_labels):
                        if ci < scores_arr.shape[1]:
                            out[f'P_{cname}'] = scores_arr[:, ci]
                for col in selected_features:
                    if set_name == 'Train':
                        out[col] = X_train_sel[col].values
                    else:
                        out[col] = X_test_sel[col].values
                return out

            train_df = build_result_frame(
                z_train.values, y_train, y_train_pred, score_train,
                hat_train, train_ad, 'Train'
            )
            test_df = build_result_frame(
                z_test.values, y_test, y_test_pred, score_test,
                hat_test, test_ad, 'Test'
            )
            full_results = pd.concat([train_df, test_df], ignore_index=True)

            # Confusion matrices use fixed class order.
            cm_train = confusion_matrix(y_train, y_train_pred, labels=np.arange(n_classes))
            cm_cv = confusion_matrix(y_train, y_cv, labels=np.arange(n_classes))
            cm_test = confusion_matrix(y_test, y_test_pred, labels=np.arange(n_classes)) if len(y_test) else None

            metrics_dict = {
                'mode': class_mode,
                'model': model_type,
                'n_classes': n_classes,
                'class_labels': class_labels,
                'selected_features': list(selected_features),
                'final_feature_count': len(final_names),
                'cv_folds': actual_cv_folds,
                'nested_cv': nested_cv,
                'rfa_cv_balacc': rfa_fitness,
                'train': train_metrics,
                'cv': cv_metrics,
                'test': test_metrics,
                'warning_leverage': warning_leverage,
                'train_outside_ad': int(np.sum(~train_ad)),
                'test_outside_ad': int(np.sum(~test_ad)) if len(test_ad) else 0,
                'yrand_mean_balacc': yrand_mean_balacc,
                'yrand_mean_mcc': yrand_mean_mcc,
                'yrand_p_balacc': yrand_p_balacc,
                'yrand_p_mcc': yrand_p_mcc,
                'score_kind': score_kind,
                'cv_score_kind': cv_score_kind,
            }

            plot_data = {
                'class_labels': class_labels,
                'n_classes': n_classes,
                'y_train': y_train,
                'y_train_pred': y_train_pred,
                'y_cv': y_cv,
                'y_test': y_test,
                'y_test_pred': y_test_pred,
                'score_train': score_train,
                'score_cv': cv_score,
                'score_test': score_test,
                'cm_train': cm_train,
                'cm_cv': cm_cv,
                'cm_test': cm_test,
                'hat_train': hat_train,
                'hat_test': hat_test,
                'warning_leverage': warning_leverage,
                'sn_train': z_train.values,
                'sn_test': z_test.values,
                'yrand_balacc': yrand_balacc,
                'yrand_mcc': yrand_mcc,
                'orig_cv_balacc': cv_metrics['Balanced_Accuracy'],
                'orig_cv_mcc': cv_metrics['MCC'],
                'heatmap_data': X_train_sel,
            }

            self.finished_signal.emit(True, 'Analysis Complete', {
                'metrics': metrics_dict,
                'results_df': full_results,
                'plot_data': plot_data,
            })

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
        self.current_metrics = None
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
        left_layout.addWidget(QLabel('<b>1. Data Source:</b>'))
        data_layout = QHBoxLayout()
        self.btn_load = QPushButton('Load Dataset CSV')
        self.btn_load.clicked.connect(self.load_data)
        data_layout.addWidget(self.btn_load)
        self.lbl_file = QLabel('No file loaded')
        data_layout.addWidget(self.lbl_file)
        data_layout.addStretch()
        left_layout.addLayout(data_layout)
        left_layout.addWidget(QLabel('Expected columns: ID | Class | numeric descriptors...'))

        # 2. Classification mode
        mode_layout = QHBoxLayout()
        mode_layout.addWidget(QLabel('<b>2. Analysis Type:</b>'))
        self.class_mode_combo = QComboBox()
        self.class_mode_combo.addItems(['Binary', 'Ternary', 'Multiclass'])
        mode_layout.addWidget(self.class_mode_combo)
        mode_layout.addStretch()
        left_layout.addLayout(mode_layout)

        # 3. Train/Test split
        left_layout.addWidget(QLabel('<b>3. Define Train/Test Split:</b>'))
        split_layout = QHBoxLayout()
        self.tbl_test = QTableWidget()
        self.tbl_test.setColumnCount(1)
        self.tbl_test.setHorizontalHeaderLabels(['Test Set IDs'])
        self.tbl_test.horizontalHeader().setStretchLastSection(True)
        split_layout.addWidget(self.tbl_test)

        btn_move_layout = QVBoxLayout()
        self.btn_to_train = QPushButton('>> To Train >>')
        self.btn_to_train.clicked.connect(self.move_to_train)
        self.btn_to_test = QPushButton('<< To Test <<')
        self.btn_to_test.clicked.connect(self.move_to_test)
        btn_move_layout.addStretch()
        btn_move_layout.addWidget(self.btn_to_train)
        btn_move_layout.addWidget(self.btn_to_test)
        btn_move_layout.addStretch()
        split_layout.addLayout(btn_move_layout)

        self.tbl_train = QTableWidget()
        self.tbl_train.setColumnCount(1)
        self.tbl_train.setHorizontalHeaderLabels(['Training Set IDs'])
        self.tbl_train.horizontalHeader().setStretchLastSection(True)
        split_layout.addWidget(self.tbl_train)
        left_layout.addLayout(split_layout)

        load_split_layout = QHBoxLayout()
        self.btn_load_train = QPushButton('Load Train Split CSV')
        self.btn_load_train.clicked.connect(self.load_train_split)
        self.btn_load_test = QPushButton('Load Test Split CSV')
        self.btn_load_test.clicked.connect(self.load_test_split)
        load_split_layout.addWidget(self.btn_load_train)
        load_split_layout.addWidget(self.btn_load_test)
        load_split_layout.addStretch()
        left_layout.addLayout(load_split_layout)

        auto_layout = QHBoxLayout()
        auto_layout.addWidget(QLabel('Train %:'))
        self.split_spin = QSpinBox()
        self.split_spin.setRange(50, 95)
        self.split_spin.setValue(80)
        auto_layout.addWidget(self.split_spin)
        self.btn_auto = QPushButton('Stratified Auto Split')
        self.btn_auto.clicked.connect(self.auto_split)
        auto_layout.addWidget(self.btn_auto)
        auto_layout.addStretch()
        left_layout.addLayout(auto_layout)

        # 4. Settings
        left_layout.addWidget(QLabel('<b>4. Model & Selection Hyperparameters:</b>'))
        settings_layout = QHBoxLayout()

        self.settings_tbl = QTableWidget()
        self.settings_tbl.setColumnCount(2)
        self.settings_tbl.setHorizontalHeaderLabels(['Setting', 'Value'])
        self.settings_tbl.horizontalHeader().setStretchLastSection(True)
        self.settings_tbl.setMaximumHeight(180)
        default_settings = [
            ('Dimension', '5'), ('Pop Size', '30'), ('Iterations', '30'),
            ('MARS Knots', '4'), ('CV Folds', '5'), ('Random Seed', '35')
        ]
        self.settings_tbl.setRowCount(len(default_settings))
        for i, (k, v) in enumerate(default_settings):
            k_item = QTableWidgetItem(k)
            k_item.setFlags(k_item.flags() & ~Qt.ItemIsEditable)
            self.settings_tbl.setItem(i, 0, k_item)
            self.settings_tbl.setItem(i, 1, QTableWidgetItem(v))
        settings_layout.addWidget(self.settings_tbl)

        model_group = QGroupBox('Model Options')
        model_layout = QVBoxLayout()
        self.model_combo = QComboBox()
        models = ['Logistic Regression', 'Ridge Classifier', 'Random Forest', 'SVC']
        if XGB_AVAILABLE:
            models.append('XGBoost')
        self.model_combo.addItems(models)
        model_layout.addWidget(QLabel('Classifier:'))
        model_layout.addWidget(self.model_combo)

        self.use_mars_check = QCheckBox('Apply MARS-like Spline Expansion')
        self.use_mars_check.setChecked(False)
        model_layout.addWidget(self.use_mars_check)

        self.scale_check = QCheckBox('Scale Descriptors')
        self.scale_check.setChecked(True)
        model_layout.addWidget(self.scale_check)

        self.class_weight_check = QCheckBox('Use balanced class weights')
        self.class_weight_check.setChecked(True)
        model_layout.addWidget(self.class_weight_check)

        self.nested_cv_check = QCheckBox('Nested CV validation (slow)')
        self.nested_cv_check.setChecked(False)
        model_layout.addWidget(self.nested_cv_check)

        self.yrand_spin = QSpinBox()
        self.yrand_spin.setRange(5, 500)
        self.yrand_spin.setValue(50)
        yr_layout = QHBoxLayout()
        yr_layout.addWidget(QLabel('Y-Rand Runs:'))
        yr_layout.addWidget(self.yrand_spin)
        model_layout.addLayout(yr_layout)

        model_group.setLayout(model_layout)
        settings_layout.addWidget(model_group)
        left_layout.addLayout(settings_layout)

        # Execute/logs
        self.btn_run = QPushButton('🚀 Run Classification QSAR Analysis')
        self.btn_run.setStyleSheet('background-color: #2980b9; color: white; font-weight: bold; padding: 8px;')
        self.btn_run.clicked.connect(self.run_analysis)
        left_layout.addWidget(self.btn_run)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        left_layout.addWidget(self.progress)

        self.txt_log = QTextEdit()
        self.txt_log.setReadOnly(True)
        left_layout.addWidget(self.txt_log)

        export_layout = QHBoxLayout()
        self.btn_export = QPushButton('Export Predictions (CSV)')
        self.btn_export.clicked.connect(self.export_results)
        self.btn_export.setEnabled(False)
        export_layout.addWidget(self.btn_export)
        self.btn_export_plot = QPushButton('Export Current Plot')
        self.btn_export_plot.clicked.connect(self.export_plot)
        export_layout.addWidget(self.btn_export_plot)
        self.dpi_spin = QSpinBox()
        self.dpi_spin.setRange(72, 1200)
        self.dpi_spin.setValue(300)
        export_layout.addWidget(QLabel('DPI:'))
        export_layout.addWidget(self.dpi_spin)
        left_layout.addLayout(export_layout)

        splitter.addWidget(left_panel)

        # RIGHT PANEL / plots
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)
        self.tabs = QTabWidget()

        self.canvas_confusion = InteractiveCanvas()
        self.canvas_roc = InteractiveCanvas()
        self.canvas_pr = InteractiveCanvas()
        self.canvas_ad = InteractiveCanvas()
        self.canvas_yrand = InteractiveCanvas()
        self.canvas_heatmap = InteractiveCanvas()

        self.tabs.addTab(self.canvas_confusion, 'Confusion Matrix')
        self.tabs.addTab(self.canvas_roc, 'ROC Curve')
        self.tabs.addTab(self.canvas_pr, 'Precision-Recall')
        self.tabs.addTab(self.canvas_ad, 'Applicability Domain')
        self.tabs.addTab(self.canvas_yrand, 'Y-Randomization')
        self.tabs.addTab(self.canvas_heatmap, 'Correlation Heatmap')

        font_layout = QHBoxLayout()
        font_layout.addWidget(QLabel('Global Font Size:'))
        self.font_spin = QSpinBox()
        self.font_spin.setRange(6, 24)
        self.font_spin.setValue(10)
        self.font_spin.valueChanged.connect(self.update_font_size)
        font_layout.addWidget(self.font_spin)
        font_layout.addStretch()
        right_layout.addLayout(font_layout)
        right_layout.addWidget(self.tabs)

        splitter.addWidget(right_panel)
        splitter.setSizes([500, 750])
        self.splitter = splitter
        self._active_popout = None

        top_bar = QHBoxLayout()
        top_bar.addStretch()
        self.btn_maximize = QPushButton('Maximize')
        self.btn_maximize.clicked.connect(self.toggle_maximize)
        top_bar.addWidget(self.btn_maximize)
        self.btn_fullscreen = QPushButton('Full Screen')
        self.btn_fullscreen.clicked.connect(lambda: self._open_popout(fullscreen=True))
        top_bar.addWidget(self.btn_fullscreen)
        self.btn_close = QPushButton('Close')
        self.btn_close.clicked.connect(self.widget.close)
        top_bar.addWidget(self.btn_close)
        main_layout.addLayout(top_bar)
        main_layout.addWidget(splitter)

    # ---------------------------
    # Data / split handling
    # ---------------------------
    def load_data(self):
        parent = self._get_dialog_parent()
        path, _ = QFileDialog.getOpenFileName(
            parent, 'Select Dataset CSV', '', 'CSV Files (*.csv)'
        )
        parent.raise_()
        parent.activateWindow()
        if not path:
            return
        try:
            df = pd.read_csv(path, sep=None, engine='python')
            if df.shape[1] < 3:
                raise ValueError('CSV needs at least 3 columns: ID, Class, descriptor(s).')
            self.dataset = df
            self.lbl_file.setText(os.path.basename(path))
            classes = sorted(df.iloc[:, 1].dropna().astype(str).unique().tolist())
            self.txt_log.append(
                f"Loaded dataset with {df.shape[0]} rows, {df.shape[1]-2} candidate descriptor columns."
            )
            self.txt_log.append(f"Detected {len(classes)} class(es): {classes}")
        except Exception as e:
            QMessageBox.critical(parent, 'Load Error', str(e))

    def load_train_split(self):
        parent = self._get_dialog_parent()
        path, _ = QFileDialog.getOpenFileName(
            parent, 'Select Train Split CSV', '', 'CSV Files (*.csv)'
        )
        parent.raise_()
        parent.activateWindow()
        if not path:
            return
        try:
            df = pd.read_csv(path, sep=None, engine='python')
            ids = sorted(df.iloc[:, 0].astype(str).tolist())
            self.tbl_train.setRowCount(len(ids))
            for i, val in enumerate(ids):
                self.tbl_train.setItem(i, 0, QTableWidgetItem(str(val)))
            self.txt_log.append(f"Loaded {len(ids)} training IDs from {os.path.basename(path)}")
        except Exception as e:
            QMessageBox.critical(parent, 'Load Error', str(e))

    def load_test_split(self):
        parent = self._get_dialog_parent()
        path, _ = QFileDialog.getOpenFileName(
            parent, 'Select Test Split CSV', '', 'CSV Files (*.csv)'
        )
        parent.raise_()
        parent.activateWindow()
        if not path:
            return
        try:
            df = pd.read_csv(path, sep=None, engine='python')
            ids = sorted(df.iloc[:, 0].astype(str).tolist())
            self.tbl_test.setRowCount(len(ids))
            for i, val in enumerate(ids):
                self.tbl_test.setItem(i, 0, QTableWidgetItem(str(val)))
            self.txt_log.append(f"Loaded {len(ids)} test IDs from {os.path.basename(path)}")
        except Exception as e:
            QMessageBox.critical(parent, 'Load Error', str(e))

    def auto_split(self):
        if self.dataset is None:
            return
        try:
            ids = self.dataset.iloc[:, 0].astype(str).to_numpy()
            y = self.dataset.iloc[:, 1].astype(str).to_numpy()
            train_pct = self.split_spin.value() / 100.0
            seed = 35
            for i in range(self.settings_tbl.rowCount()):
                if self.settings_tbl.item(i, 0).text() == 'Random Seed':
                    seed = int(self.settings_tbl.item(i, 1).text())
                    break

            if minimum_class_count(y) < 2:
                raise ValueError('Each class needs at least 2 samples for stratified splitting.')

            train_ids, test_ids = train_test_split(
                ids, train_size=train_pct, random_state=seed, stratify=y
            )
            train_ids = sorted(train_ids.tolist())
            test_ids = sorted(test_ids.tolist())

            self.tbl_train.setRowCount(len(train_ids))
            for i, val in enumerate(train_ids):
                self.tbl_train.setItem(i, 0, QTableWidgetItem(str(val)))
            self.tbl_test.setRowCount(len(test_ids))
            for i, val in enumerate(test_ids):
                self.tbl_test.setItem(i, 0, QTableWidgetItem(str(val)))
            self.txt_log.append(
                f"Created stratified split: train={len(train_ids)}, test={len(test_ids)}."
            )
        except Exception as e:
            QMessageBox.critical(self.widget, 'Split Error', str(e))

    def move_to_train(self):
        selected_items = self.tbl_test.selectedItems()
        if not selected_items:
            return
        rows = sorted({item.row() for item in selected_items}, reverse=True)
        for row in rows:
            item = self.tbl_test.takeItem(row, 0)
            self.tbl_train.insertRow(self.tbl_train.rowCount())
            self.tbl_train.setItem(self.tbl_train.rowCount() - 1, 0, item)

    def move_to_test(self):
        selected_items = self.tbl_train.selectedItems()
        if not selected_items:
            return
        rows = sorted({item.row() for item in selected_items}, reverse=True)
        for row in rows:
            item = self.tbl_train.takeItem(row, 0)
            self.tbl_test.insertRow(self.tbl_test.rowCount())
            self.tbl_test.setItem(self.tbl_test.rowCount() - 1, 0, item)

    # ---------------------------
    # Execute / output
    # ---------------------------
    def run_analysis(self):
        if self.dataset is None or self.tbl_train.rowCount() == 0:
            QMessageBox.warning(self.widget, 'Warning', 'Load data and assign a training set first.')
            return

        self.btn_run.setEnabled(False)
        self.btn_export.setEnabled(False)
        self.txt_log.clear()

        config = {}
        try:
            for i in range(self.settings_tbl.rowCount()):
                config[self.settings_tbl.item(i, 0).text()] = self.settings_tbl.item(i, 1).text()
            config['Model Type'] = self.model_combo.currentText()
            config['Class Mode'] = self.class_mode_combo.currentText()
            config['Use MARS'] = self.use_mars_check.isChecked()
            config['Scale Descriptors'] = self.scale_check.isChecked()
            config['Class Weight'] = self.class_weight_check.isChecked()
            config['Nested CV'] = self.nested_cv_check.isChecked()
            config['Y-Rand Runs'] = self.yrand_spin.value()

            # Early mode validation gives immediate UI feedback.
            n_classes = self.dataset.iloc[:, 1].dropna().astype(str).nunique()
            validate_class_mode(n_classes, config['Class Mode'])
        except Exception as e:
            self.btn_run.setEnabled(True)
            QMessageBox.warning(self.widget, 'Configuration Error', str(e))
            return

        train_ids = [
            self.tbl_train.item(i, 0).text()
            for i in range(self.tbl_train.rowCount())
            if self.tbl_train.item(i, 0) is not None
        ]
        test_ids = [
            self.tbl_test.item(i, 0).text()
            for i in range(self.tbl_test.rowCount())
            if self.tbl_test.item(i, 0) is not None
        ]

        if set(train_ids) & set(test_ids):
            self.btn_run.setEnabled(True)
            QMessageBox.warning(self.widget, 'Split Error', 'Training and test ID lists overlap.')
            return

        self.progress.setValue(0)
        self.progress.setVisible(True)
        self.worker = RfaMarsWorker(
            self.dataset.copy(), train_ids, test_ids, config
        )
        self.worker.log_signal.connect(self.txt_log.append)
        self.worker.progress_signal.connect(self.progress.setValue)
        self.worker.max_progress_signal.connect(self.progress.setMaximum)
        self.worker.finished_signal.connect(self.handle_finish)
        self.worker.start()

    def handle_finish(self, success, msg, result_data):
        self.btn_run.setEnabled(True)
        self.progress.setVisible(False)
        if not success:
            QMessageBox.critical(self.widget, 'Error', msg)
            return

        self.result_cache = result_data['results_df']
        self.current_plot_data = result_data['plot_data']
        self.current_metrics = result_data['metrics']
        self.btn_export.setEnabled(True)
        m = self.current_metrics

        def fmt(v):
            return 'NA' if v is None or (isinstance(v, float) and np.isnan(v)) else f'{v:.4f}'

        tr, cv, te = m['train'], m['cv'], m['test']
        metrics_str = f"""
========================================
CLASSIFICATION QSAR VALIDATION METRICS
========================================
Mode          : {m['mode']}
Model         : {m['model']}
Classes       : {m['class_labels']}
Selected      : {m['selected_features']}
CV folds      : {m['cv_folds']}
Nested CV     : {'Yes' if m['nested_cv'] else 'No'}
RFA CV BalAcc : {m['rfa_cv_balacc']:.4f}

--- Training set ---
Accuracy      : {fmt(tr['Accuracy'])}
Balanced Acc. : {fmt(tr['Balanced_Accuracy'])}
MCC           : {fmt(tr['MCC'])}
Cohen kappa   : {fmt(tr['Kappa'])}
Precision mac.: {fmt(tr['Precision_macro'])}
Recall macro  : {fmt(tr['Recall_macro'])}
F1 macro      : {fmt(tr['F1_macro'])}
F1 weighted   : {fmt(tr['F1_weighted'])}
ROC-AUC       : {fmt(tr['ROC_AUC'])}
PR-AUC        : {fmt(tr['PR_AUC'])}

--- Stratified cross-validation ---
Accuracy      : {fmt(cv['Accuracy'])}
Balanced Acc. : {fmt(cv['Balanced_Accuracy'])}
MCC           : {fmt(cv['MCC'])}
Cohen kappa   : {fmt(cv['Kappa'])}
Precision mac.: {fmt(cv['Precision_macro'])}
Recall macro  : {fmt(cv['Recall_macro'])}
F1 macro      : {fmt(cv['F1_macro'])}
F1 weighted   : {fmt(cv['F1_weighted'])}
ROC-AUC       : {fmt(cv['ROC_AUC'])}
PR-AUC        : {fmt(cv['PR_AUC'])}

--- External prediction ---
Accuracy      : {fmt(te['Accuracy'])}
Balanced Acc. : {fmt(te['Balanced_Accuracy'])}
MCC           : {fmt(te['MCC'])}
Cohen kappa   : {fmt(te['Kappa'])}
Precision mac.: {fmt(te['Precision_macro'])}
Recall macro  : {fmt(te['Recall_macro'])}
F1 macro      : {fmt(te['F1_macro'])}
F1 weighted   : {fmt(te['F1_weighted'])}
ROC-AUC       : {fmt(te['ROC_AUC'])}
PR-AUC        : {fmt(te['PR_AUC'])}

--- Applicability domain ---
Warning h*    : {m['warning_leverage']:.4f}
Train outside : {m['train_outside_ad']}
Test outside  : {m['test_outside_ad']}

--- Y-randomization ---
Mean BalAcc   : {fmt(m['yrand_mean_balacc'])}
Mean MCC      : {fmt(m['yrand_mean_mcc'])}
Perm. p BalAcc: {fmt(m['yrand_p_balacc'])}
Perm. p MCC   : {fmt(m['yrand_p_mcc'])}
========================================
"""
        self.txt_log.append(metrics_str)
        self.render_all_plots()
        QMessageBox.information(
            self.widget, 'Success',
            'Classification analysis complete. View metrics, plots and exportable predictions.'
        )

    # ---------------------------
    # Plotting
    # ---------------------------
    def update_font_size(self, size):
        plt.rcParams.update({
            'font.size': size, 'axes.titlesize': size + 2,
            'axes.labelsize': size
        })
        if self.current_plot_data:
            self.render_all_plots()

    def render_all_plots(self):
        if not self.current_plot_data:
            return
        self.plot_confusion(self.current_plot_data)
        self.plot_roc(self.current_plot_data)
        self.plot_pr(self.current_plot_data)
        self.plot_ad(self.current_plot_data)
        self.plot_yrandomization(self.current_plot_data)
        self.plot_heatmap(self.current_plot_data['heatmap_data'])

    def plot_confusion(self, pd_data):
        ax = self.canvas_confusion.ax
        ax.clear()
        cm = pd_data['cm_test'] if pd_data['cm_test'] is not None else pd_data['cm_cv']
        title_suffix = 'External Test' if pd_data['cm_test'] is not None else 'Cross-Validation'
        im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
        self.canvas_confusion.fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        labels = pd_data['class_labels']
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_yticklabels(labels)
        ax.set_xlabel('Predicted class')
        ax.set_ylabel('Actual class')
        ax.set_title(f'Confusion Matrix - {title_suffix}')
        threshold = cm.max() / 2.0 if cm.size else 0
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, str(cm[i, j]), ha='center', va='center',
                        color='white' if cm[i, j] > threshold else 'black')
        self.canvas_confusion.fig.tight_layout()
        self.canvas_confusion.clear_interactive_data()
        self.canvas_confusion.draw()

    def _score_source_for_plot(self, pd_data):
        if pd_data['score_test'] is not None and len(pd_data['y_test']):
            return pd_data['y_test'], pd_data['score_test'], 'External Test'
        if pd_data['score_cv'] is not None:
            return pd_data['y_train'], pd_data['score_cv'], 'Cross-Validation'
        return None, None, None

    def plot_roc(self, pd_data):
        ax = self.canvas_roc.ax
        ax.clear()
        y_true, scores, source = self._score_source_for_plot(pd_data)
        if y_true is None:
            ax.text(0.5, 0.5, 'ROC unavailable for this classifier', ha='center', va='center')
            ax.set_axis_off()
            self.canvas_roc.draw()
            return

        n_classes = pd_data['n_classes']
        labels = pd_data['class_labels']
        try:
            if n_classes == 2:
                fpr, tpr, _ = roc_curve(y_true, np.asarray(scores)[:, 1])
                auc = roc_auc_score(y_true, np.asarray(scores)[:, 1])
                ax.plot(fpr, tpr, label=f'ROC AUC={auc:.3f}')
            else:
                y_bin = label_binarize(y_true, classes=np.arange(n_classes))
                for i, label in enumerate(labels):
                    if len(np.unique(y_bin[:, i])) < 2:
                        continue
                    fpr, tpr, _ = roc_curve(y_bin[:, i], np.asarray(scores)[:, i])
                    auc = roc_auc_score(y_bin[:, i], np.asarray(scores)[:, i])
                    ax.plot(fpr, tpr, label=f'{label} AUC={auc:.3f}')
            ax.plot([0, 1], [0, 1], 'k--')
            ax.set(xlabel='False Positive Rate', ylabel='True Positive Rate',
                   title=f'ROC Curve - {source}')
            ax.legend(fontsize=max(6, plt.rcParams['font.size'] - 1))
            ax.grid(True)
        except Exception as e:
            ax.text(0.5, 0.5, f'ROC unavailable:\n{e}', ha='center', va='center')
            ax.set_axis_off()
        self.canvas_roc.clear_interactive_data()
        self.canvas_roc.draw()

    def plot_pr(self, pd_data):
        ax = self.canvas_pr.ax
        ax.clear()
        y_true, scores, source = self._score_source_for_plot(pd_data)
        if y_true is None:
            ax.text(0.5, 0.5, 'Precision-recall curve unavailable', ha='center', va='center')
            ax.set_axis_off()
            self.canvas_pr.draw()
            return

        n_classes = pd_data['n_classes']
        labels = pd_data['class_labels']
        try:
            if n_classes == 2:
                precision, recall, _ = precision_recall_curve(y_true, np.asarray(scores)[:, 1])
                ap = average_precision_score(y_true, np.asarray(scores)[:, 1])
                ax.plot(recall, precision, label=f'AP={ap:.3f}')
            else:
                y_bin = label_binarize(y_true, classes=np.arange(n_classes))
                for i, label in enumerate(labels):
                    if np.sum(y_bin[:, i]) == 0:
                        continue
                    precision, recall, _ = precision_recall_curve(y_bin[:, i], np.asarray(scores)[:, i])
                    ap = average_precision_score(y_bin[:, i], np.asarray(scores)[:, i])
                    ax.plot(recall, precision, label=f'{label} AP={ap:.3f}')
            ax.set(xlabel='Recall', ylabel='Precision', title=f'Precision-Recall - {source}')
            ax.legend(fontsize=max(6, plt.rcParams['font.size'] - 1))
            ax.grid(True)
        except Exception as e:
            ax.text(0.5, 0.5, f'PR curve unavailable:\n{e}', ha='center', va='center')
            ax.set_axis_off()
        self.canvas_pr.clear_interactive_data()
        self.canvas_pr.draw()

    def plot_ad(self, pd_data):
        ax = self.canvas_ad.ax
        ax.clear()
        train_x = np.arange(len(pd_data['hat_train']))
        ax.scatter(train_x, pd_data['hat_train'], label='Training')
        if len(pd_data['hat_test']):
            test_x = np.arange(len(pd_data['hat_train']), len(pd_data['hat_train']) + len(pd_data['hat_test']))
            ax.scatter(test_x, pd_data['hat_test'], label='External test')
        ax.axhline(pd_data['warning_leverage'], color='r', linestyle='--', label='h*')
        ax.set(xlabel='Compound index', ylabel='Leverage', title='Descriptor-Space Applicability Domain')
        ax.legend()
        ax.grid(True)
        all_x = np.concatenate([
            train_x,
            np.arange(len(pd_data['hat_train']), len(pd_data['hat_train']) + len(pd_data['hat_test']))
        ]) if len(pd_data['hat_test']) else train_x
        all_y = np.concatenate([pd_data['hat_train'], pd_data['hat_test']]) if len(pd_data['hat_test']) else pd_data['hat_train']
        all_labels = np.concatenate([pd_data['sn_train'], pd_data['sn_test']]) if len(pd_data['sn_test']) else pd_data['sn_train']
        self.canvas_ad.set_interactive_data(all_x, all_y, all_labels)
        self.canvas_ad.draw()

    def plot_yrandomization(self, pd_data):
        ax = self.canvas_yrand.ax
        ax.clear()
        x = np.arange(1, len(pd_data['yrand_balacc']) + 1)
        ax.scatter(x, pd_data['yrand_balacc'], label='Random CV Balanced Accuracy', alpha=0.7)
        ax.scatter(x, pd_data['yrand_mcc'], label='Random CV MCC', alpha=0.7)
        ax.axhline(pd_data['orig_cv_balacc'], linestyle='--', label='Original CV Balanced Accuracy')
        ax.axhline(pd_data['orig_cv_mcc'], linestyle=':', label='Original CV MCC')
        ax.set(xlabel='Permutation run', ylabel='Metric value', title='Class-label Y-Randomization')
        ax.legend(fontsize=max(6, plt.rcParams['font.size'] - 1))
        ax.grid(True)
        self.canvas_yrand.clear_interactive_data()
        self.canvas_yrand.draw()

    def plot_heatmap(self, df_descriptors):
        self.canvas_heatmap.fig.clf()
        ax = self.canvas_heatmap.fig.add_subplot(111)
        self.canvas_heatmap.ax = ax
        corr = df_descriptors.corr()
        cax = ax.imshow(corr, cmap='coolwarm', aspect='auto', vmin=-1, vmax=1)
        self.canvas_heatmap.fig.colorbar(cax, ax=ax)
        ax.set_xticks(np.arange(len(corr.columns)))
        ax.set_yticks(np.arange(len(corr.columns)))
        ax.set_xticklabels(corr.columns, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(corr.columns, fontsize=8)
        ax.set_title('Selected Descriptor Correlation Heatmap')
        self.canvas_heatmap.fig.tight_layout()
        self.canvas_heatmap.clear_interactive_data()
        self.canvas_heatmap.draw()

    # ---------------------------
    # Export/window helpers
    # ---------------------------
    def export_results(self):
        if self.result_cache is None:
            return
        parent = self._get_dialog_parent()
        path, _ = QFileDialog.getSaveFileName(
            parent, 'Save Classification Results',
            'QSAR_Classification_Results.csv', 'CSV (*.csv)'
        )
        parent.raise_()
        parent.activateWindow()
        if path:
            try:
                self.result_cache.to_csv(path, index=False)
                QMessageBox.information(parent, 'Export Success', f'Results saved to {path}')
            except Exception as e:
                QMessageBox.critical(parent, 'Export Error', str(e))

    def export_plot(self):
        current_widget = self.tabs.currentWidget()
        if not isinstance(current_widget, FigureCanvas):
            return
        parent = self._get_dialog_parent()
        path, _ = QFileDialog.getSaveFileName(
            parent, 'Save Plot', 'plot.png', 'Images (*.png *.jpg *.jpeg)'
        )
        parent.raise_()
        parent.activateWindow()
        if path:
            try:
                current_widget.fig.savefig(
                    path, dpi=self.dpi_spin.value(), bbox_inches='tight'
                )
                QMessageBox.information(
                    parent, 'Success',
                    f'Plot saved successfully at {self.dpi_spin.value()} DPI.'
                )
            except Exception as e:
                QMessageBox.critical(parent, 'Export Error', f'Could not save image:\n{e}')

    def _open_popout(self, fullscreen=False):
        if self._active_popout is not None:
            self._active_popout.raise_()
            self._active_popout.activateWindow()
            return
        dialog = QDialog()
        dialog.setWindowTitle('RFA + MARS Classification QSAR')
        dialog.setWindowFlags(Qt.Window)
        vbox = QVBoxLayout(dialog)
        vbox.setContentsMargins(4, 4, 4, 4)
        self.splitter.setParent(dialog)
        vbox.addWidget(self.splitter)
        dialog.finished.connect(self._restore_from_popout)
        self._active_popout = dialog
        if fullscreen:
            dialog.showFullScreen()
        else:
            dialog.showMaximized()

    def _get_dialog_parent(self):
        """Return the currently active top-level window for use as a QFileDialog parent.
        When the plugin is maximized/popped out, returns the popout dialog so that
        native file dialogs are owned by the visible window and don't cause it to go
        behind the main application window on Windows."""
        if getattr(self, '_active_popout', None) is not None:
            return self._active_popout
        return self.widget

    def _restore_from_popout(self):
        self.splitter.setParent(self.widget)
        self.widget.layout().addWidget(self.splitter)
        self._active_popout = None
        self.widget.update()

    def toggle_maximize(self):
        self._open_popout(fullscreen=False)


# ===================================================================
# PLUGIN WRAPPER
# ===================================================================
class QsarRfaPlugin(BasePlugin):
    def __init__(self):
        super().__init__(PluginInfo(
            name='RFA + MARS Classification QSAR',
            version='2.0.0',
            description=(
                'RFA descriptor selection and optional MARS-like spline expansion '
                'for binary, ternary and multiclass QSAR classification with stratified validation'
            ),
            author='SMILES Team',
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
        self.logger.info('RFA+MARS Classification QSAR plugin initialized')
        return True

    def cleanup(self):
        if self.widget:
            self.widget.widget.deleteLater()
            self.widget = None
        self.logger.info('RFA+MARS Classification QSAR plugin cleaned up')

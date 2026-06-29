#!/usr/bin/env python
import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

META = {"sample_id", "label", "cancer_type", "batch", "sex", "age", "smoking_status"}
NEGATIVE_LABELS = {"0", "false", "control", "healthy", "normal", "negative", "no", "non_cancer", "non-cancer"}
POSITIVE_LABELS = {"1", "true", "cancer", "case", "positive", "yes"}


def parse_binary_labels(labels):
    normalized = labels.astype(str).str.strip().str.lower()
    y = []
    for label in normalized:
        if label in NEGATIVE_LABELS or label == "" or label == "nan":
            y.append(0)
        elif label in POSITIVE_LABELS:
            y.append(1)
        else:
            # For cohort labels such as lung_cancer, colorectal_cancer, etc.
            # treat non-control disease labels as positive in binary mode.
            y.append(1)
    y = pd.Series(y, index=labels.index, dtype=int)
    if y.nunique() != 2:
        counts = normalized.value_counts(dropna=False).to_dict()
        raise SystemExit(f"Training requires two binary classes after label parsing; observed labels: {counts}")
    return y


def threshold_at_specificity(y_true, scores, specificity):
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores, dtype=float)
    controls = scores[y_true == 0]
    if len(controls) == 0:
        return 0.5
    return float(np.quantile(controls, specificity))


def sensitivity_at_threshold(y_true, scores, threshold):
    y_true = np.asarray(y_true).astype(int)
    positives = y_true == 1
    if positives.sum() == 0:
        return None
    return float((scores[positives] >= threshold).mean())


def metric_block(y_true, scores, threshold=None):
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores, dtype=float)
    out = {}
    if len(np.unique(y_true)) == 2:
        out["roc_auc"] = float(roc_auc_score(y_true, scores))
        out["pr_auc"] = float(average_precision_score(y_true, scores))
    if threshold is not None:
        pred = scores > threshold
        out["tp"] = int(((y_true == 1) & pred).sum())
        out["tn"] = int(((y_true == 0) & ~pred).sum())
        out["fp"] = int(((y_true == 0) & pred).sum())
        out["fn"] = int(((y_true == 1) & ~pred).sum())
        out["sensitivity"] = out["tp"] / max(out["tp"] + out["fn"], 1)
        out["specificity"] = out["tn"] / max(out["tn"] + out["fp"], 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--out-model", required=True)
    ap.add_argument("--out-metadata", required=True)
    ap.add_argument("--out-metrics", required=True)
    ap.add_argument("--out-cv-predictions")
    ap.add_argument("--out-test-predictions")
    ap.add_argument("--max-features", type=int, default=5000)
    ap.add_argument("--test-size", type=float, default=0.25)
    ap.add_argument("--random-state", type=int, default=7)
    args = ap.parse_args()

    df = pd.read_csv(args.matrix, sep="\t")
    if "sample_id" in df.columns:
        df = df.sort_values("sample_id", kind="mergesort").reset_index(drop=True)
    if "label" not in df or df["label"].isna().all():
        raise SystemExit("label column is required for training")
    y = parse_binary_labels(df["label"])
    feature_cols = [c for c in df.columns if c not in META]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    usable = X.notna().any(axis=0)
    X = X.loc[:, usable]
    feature_cols = X.columns.tolist()
    if not feature_cols:
        raise SystemExit("No usable non-constant numeric features found for training")

    min_class = int(y.value_counts().min())
    can_holdout = min_class >= 2 and len(y) >= 8 and args.test_size > 0
    if can_holdout:
        train_idx, test_idx = train_test_split(
            np.arange(len(df)),
            test_size=args.test_size,
            stratify=y,
            random_state=args.random_state,
        )
    else:
        train_idx = np.arange(len(df))
        test_idx = np.array([], dtype=int)

    X_train, y_train = X.iloc[train_idx], y.iloc[train_idx]
    X_test = X.iloc[test_idx]
    y_test = y.iloc[test_idx] if len(test_idx) else pd.Series(dtype=int)
    selector_k = "all"
    if args.max_features and args.max_features > 0:
        selector_k = min(args.max_features, X_train.shape[1])

    candidates = {
        "logistic_l2": LogisticRegression(max_iter=1000, penalty="l2", solver="liblinear", class_weight="balanced"),
        "random_forest": RandomForestClassifier(n_estimators=200, max_features="sqrt", class_weight="balanced", random_state=7, n_jobs=-1),
        "gradient_boosting": GradientBoostingClassifier(random_state=7),
    }
    metrics = {}
    cv_probabilities = {}
    best_name, best_auc = None, -1.0
    train_min_class = int(y_train.value_counts().min())
    cv = StratifiedKFold(n_splits=min(5, train_min_class), shuffle=True, random_state=args.random_state) if train_min_class >= 2 else None
    for name, clf in candidates.items():
        steps = [
            ("impute", SimpleImputer(strategy="median")),
            ("select", SelectKBest(score_func=f_classif, k=selector_k)),
        ]
        if name == "logistic_l2":
            steps.append(("scale", StandardScaler(with_mean=False)))
        steps.append(("model", clf))
        pipe = Pipeline(steps)
        try:
            if cv is None:
                pipe.fit(X_train, y_train)
                probs = pipe.predict_proba(X_train)[:, 1]
            else:
                probs = cross_val_predict(pipe, X_train, y_train, cv=cv, method="predict_proba")[:, 1]
            cv_probabilities[name] = probs
            auc = roc_auc_score(y_train, probs) if y_train.nunique() == 2 else None
            pr = average_precision_score(y_train, probs) if y_train.nunique() == 2 else None
        except Exception as exc:
            metrics[name] = {"error": str(exc)}
            continue
        fixed = {}
        for spec in [0.95, 0.98, 0.99]:
            thr = threshold_at_specificity(y_train, probs, spec)
            fixed[str(spec)] = {"threshold": thr, "sensitivity": sensitivity_at_threshold(y_train, probs, thr)}
        metrics[name] = {"roc_auc": auc, "pr_auc": pr, "fixed_specificity": fixed, "evaluation_scope": "cross_validated_training_partition"}
        if auc is not None and auc > best_auc:
            best_name, best_auc = name, auc

    if best_name is None:
        raise SystemExit("No model trained successfully")
    best_steps = [
        ("impute", SimpleImputer(strategy="median")),
        ("select", SelectKBest(score_func=f_classif, k=selector_k)),
    ]
    if best_name == "logistic_l2":
        best_steps.append(("scale", StandardScaler(with_mean=False)))
    best_steps.append(("model", candidates[best_name]))
    best = Pipeline(best_steps)
    best.fit(X_train, y_train)
    threshold_95 = metrics[best_name]["fixed_specificity"]["0.95"]["threshold"]
    test_metrics = None
    test_probs = None
    if len(test_idx):
        test_probs = best.predict_proba(X_test)[:, 1]
        test_metrics = metric_block(y_test, test_probs, threshold_95)
        metrics[best_name]["heldout_test"] = test_metrics
    artifact = {
        "pipeline": best,
        "feature_cols": feature_cols,
        "threshold_95_specificity": threshold_95,
        "model_name": best_name,
        "threshold_source": "cross_validated_training_scores",
        "training_scope": "train_partition_only",
    }
    joblib.dump(artifact, args.out_model)
    Path(args.out_metadata).write_text(json.dumps({
        "model_name": best_name,
        "model_version": "0.1.0",
        "analysis_language": "research use only",
        "input_feature_count": len(feature_cols),
        "max_features": args.max_features,
        "class_counts": {str(k): int(v) for k, v in y.value_counts().sort_index().items()},
        "train_class_counts": {str(k): int(v) for k, v in y_train.value_counts().sort_index().items()},
        "test_class_counts": {str(k): int(v) for k, v in y_test.value_counts().sort_index().items()} if len(test_idx) else {},
        "test_size": args.test_size if len(test_idx) else 0,
        "random_state": args.random_state,
        "split_id_column": "sample_id" if "sample_id" in df.columns else None,
        "split_order": "sample_id_sorted" if "sample_id" in df.columns else "input_row_order",
        "train_sample_ids": df["sample_id"].iloc[train_idx].astype(str).tolist() if "sample_id" in df.columns else [],
        "test_sample_ids": df["sample_id"].iloc[test_idx].astype(str).tolist() if len(test_idx) and "sample_id" in df.columns else [],
        "threshold_95_specificity": threshold_95,
        "threshold_source": "cross_validated_training_scores",
        "training_scope": "train_partition_only",
    }, indent=2) + "\n")
    Path(args.out_metrics).write_text(json.dumps(metrics, indent=2) + "\n")
    if args.out_cv_predictions and best_name in cv_probabilities:
        cv_probs = cv_probabilities[best_name]
        cv_prob_all = pd.Series(cv_probs, index=train_idx).reindex(np.arange(len(df)))
        cv_class_all = pd.Series(
            np.where(cv_prob_all > threshold_95, "elevated_fragmentome_score", "lower_fragmentome_score"),
            index=np.arange(len(df)),
        )
        cv_class_all[cv_prob_all.isna()] = pd.NA
        pd.DataFrame({
            "sample_id": df["sample_id"],
            "label": df["label"],
            "partition": np.where(np.isin(np.arange(len(df)), train_idx), "train", "test"),
            "y_true": y,
            "cv_probability": cv_prob_all.to_numpy(),
            "cv_classification": cv_class_all.to_numpy(),
            "threshold_used": threshold_95,
            "model_name": best_name,
        }).to_csv(args.out_cv_predictions, sep="\t", index=False)
    if args.out_test_predictions and len(test_idx) and test_probs is not None:
        pd.DataFrame({
            "sample_id": df["sample_id"].iloc[test_idx].to_numpy(),
            "label": df["label"].iloc[test_idx].to_numpy(),
            "y_true": y_test.to_numpy(),
            "test_probability": test_probs,
            "test_classification": np.where(test_probs > threshold_95, "elevated_fragmentome_score", "lower_fragmentome_score"),
            "threshold_used": threshold_95,
            "model_name": best_name,
        }).to_csv(args.out_test_predictions, sep="\t", index=False)


if __name__ == "__main__":
    main()

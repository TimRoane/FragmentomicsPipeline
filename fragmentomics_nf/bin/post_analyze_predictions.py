#!/usr/bin/env python
import argparse
import html
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

NEGATIVE_LABELS = {"0", "false", "control", "healthy", "normal", "negative", "no", "non_cancer", "non-cancer"}


def parse_labels(labels):
    norm = labels.astype(str).str.strip().str.lower()
    return (~norm.isin(NEGATIVE_LABELS | {"", "nan"})).astype(int)


def write_empty_plot(path, title):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7, 4))
    plt.text(0.5, 0.5, "No labeled data available", ha="center", va="center")
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def fmt_number(value, digits=3):
    if value is None:
        return "NA"
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def fmt_int(value):
    if value is None:
        return "NA"
    try:
        if pd.isna(value):
            return "NA"
        return f"{int(value):,}"
    except Exception:
        return str(value)


def pct(value):
    if value is None:
        return "NA"
    try:
        if pd.isna(value):
            return "NA"
        return f"{100 * float(value):.1f}%"
    except Exception:
        return str(value)


def metric_card(label, value, note=""):
    note_html = f"<span>{html.escape(note)}</span>" if note else ""
    return (
        '<div class="metric-card">'
        f'<div class="metric-label">{html.escape(label)}</div>'
        f'<div class="metric-value">{html.escape(str(value))}</div>'
        f'{note_html}'
        '</div>'
    )


def plot_card(title, image, caption=""):
    caption_html = f"<p>{html.escape(caption)}</p>" if caption else ""
    return (
        '<section class="plot-card">'
        f'<h3>{html.escape(title)}</h3>'
        f'<img src="{html.escape(image)}" alt="{html.escape(title)}">'
        f'{caption_html}'
        '</section>'
    )


def records_table(records, columns, empty_text="No records available."):
    if not records:
        return f'<p class="muted">{html.escape(empty_text)}</p>'
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    rows = []
    for record in records:
        cells = []
        for key, _ in columns:
            value = record.get(key, "")
            if isinstance(value, float):
                value = fmt_number(value, 4)
            cells.append(f"<td>{html.escape(str(value))}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return '<div class="table-wrap"><table><thead><tr>' + head + "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"


def build_html_report(summary, df, cv_df, test_df):
    model_metadata = summary.get("model_metadata", {})
    model_metrics = summary.get("model_metrics", {})
    model_name = model_metadata.get("model_name", "model")
    selected_metrics = model_metrics.get(model_name, {}) if isinstance(model_metrics, dict) else {}
    heldout = selected_metrics.get("heldout_test", {})
    cv_summary = summary.get("cross_validated_predictions", {})
    full_summary = {
        "roc_auc": summary.get("roc_auc"),
        "pr_auc": summary.get("pr_auc"),
        "sensitivity": summary.get("sensitivity"),
        "specificity": summary.get("specificity"),
    }

    threshold = model_metadata.get("threshold_95_specificity")
    train_counts = model_metadata.get("train_class_counts", {})
    test_counts = model_metadata.get("test_class_counts", {})
    total_counts = model_metadata.get("class_counts", {})
    n_train = sum(int(v) for v in train_counts.values()) if train_counts else None
    n_test = sum(int(v) for v in test_counts.values()) if test_counts else None

    heldout_cards = [
        metric_card("ROC AUC", fmt_number(heldout.get("roc_auc")), "Held-out test"),
        metric_card("PR AUC", fmt_number(heldout.get("pr_auc")), "Held-out test"),
        metric_card("Sensitivity", pct(heldout.get("sensitivity")), "At selected threshold"),
        metric_card("Specificity", pct(heldout.get("specificity")), "At selected threshold"),
        metric_card("Test Samples", fmt_int(n_test), f"{fmt_int(test_counts.get('0'))} controls / {fmt_int(test_counts.get('1'))} cancers"),
        metric_card("Threshold", fmt_number(threshold, 4), "From CV training scores"),
    ]

    model_cards = [
        metric_card("Model", model_name),
        metric_card("Input Features", fmt_int(model_metadata.get("input_feature_count"))),
        metric_card("Selected Features", fmt_int(model_metadata.get("max_features"))),
        metric_card("Training Samples", fmt_int(n_train), f"{fmt_int(train_counts.get('0'))} controls / {fmt_int(train_counts.get('1'))} cancers"),
        metric_card("Total Cohort", fmt_int(summary.get("n_samples")), f"{fmt_int(total_counts.get('0'))} controls / {fmt_int(total_counts.get('1'))} cancers"),
        metric_card("Random State", model_metadata.get("random_state", "NA")),
    ]

    confusion = {
        "True Positive": heldout.get("tp"),
        "False Negative": heldout.get("fn"),
        "False Positive": heldout.get("fp"),
        "True Negative": heldout.get("tn"),
    }
    confusion_html = "".join(
        f'<div class="confusion-cell"><strong>{html.escape(k)}</strong><span>{fmt_int(v)}</span></div>'
        for k, v in confusion.items()
    )

    misclassified_records = []
    if test_df is not None and {"sample_id", "label", "y_true", "test_probability", "test_classification"}.issubset(test_df.columns):
        wrong = test_df[
            ((test_df["y_true"] == 1) & (test_df["test_classification"] != "elevated_fragmentome_score"))
            | ((test_df["y_true"] == 0) & (test_df["test_classification"] == "elevated_fragmentome_score"))
        ].copy()
        wrong["error_type"] = np.where(wrong["y_true"] == 1, "False negative", "False positive")
        wrong = wrong.sort_values("test_probability", ascending=False)
        misclassified_records = wrong[["sample_id", "label", "error_type", "test_probability", "test_classification"]].to_dict("records")

    top_records = []
    if {"sample_id", "label", "cancer_probability", "classification"}.issubset(df.columns):
        top_records = (
            df.sort_values("cancer_probability", ascending=False)
            [["sample_id", "label", "cancer_probability", "classification"]]
            .head(12)
            .to_dict("records")
        )

    warning_items = []
    for key in ["evaluation_warning", "score_warning", "classification_warning", "metric_warning", "cv_metric_warning", "test_metric_warning"]:
        if key in summary:
            warning_items.append(f"<li>{html.escape(summary[key])}</li>")
    warnings_html = (
        '<section class="notice"><h2>Interpretation Notes</h2><ul>' + "".join(warning_items) + "</ul></section>"
        if warning_items
        else ""
    )

    cv_cards = [
        metric_card("CV ROC AUC", fmt_number(cv_summary.get("roc_auc")), "Training partition only"),
        metric_card("CV PR AUC", fmt_number(cv_summary.get("pr_auc")), "Training partition only"),
        metric_card("CV Samples", fmt_int(cv_summary.get("n_samples"))),
        metric_card("Unique CV Scores", fmt_int(cv_summary.get("unique_score_count"))),
    ]
    full_cards = [
        metric_card("Full-Cohort ROC AUC", fmt_number(full_summary.get("roc_auc")), "Same cohort predictions"),
        metric_card("Full-Cohort PR AUC", fmt_number(full_summary.get("pr_auc")), "Same cohort predictions"),
        metric_card("Full-Cohort Sensitivity", pct(full_summary.get("sensitivity")), "Thresholded calls"),
        metric_card("Full-Cohort Specificity", pct(full_summary.get("specificity")), "Thresholded calls"),
    ]

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Fragmentomics Model Evaluation Report</title>
  <style>
    :root {{
      --ink: #18212f;
      --muted: #647184;
      --line: #d9e0e8;
      --surface: #ffffff;
      --surface-alt: #f5f7fa;
      --accent: #0f766e;
      --accent-2: #b42318;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background: #eef2f6;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.45;
    }}
    header {{
      background: #111827;
      color: white;
      padding: 34px 44px 30px;
      border-bottom: 5px solid var(--accent);
    }}
    header h1 {{ margin: 0 0 8px; font-size: 30px; letter-spacing: 0; }}
    header p {{ margin: 0; color: #cbd5e1; max-width: 980px; }}
    main {{ max-width: 1220px; margin: 0 auto; padding: 28px; }}
    section {{
      background: var(--surface);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 24px;
      margin-bottom: 22px;
      box-shadow: 0 1px 2px rgba(17, 24, 39, 0.04);
    }}
    h2 {{ margin: 0 0 14px; font-size: 21px; }}
    h3 {{ margin: 0 0 12px; font-size: 16px; }}
    p {{ margin-top: 0; }}
    .muted {{ color: var(--muted); }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
      gap: 12px;
      margin-top: 16px;
    }}
    .metric-card {{
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--surface-alt);
      padding: 15px;
      min-height: 104px;
    }}
    .metric-label {{
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      margin-bottom: 8px;
    }}
    .metric-value {{
      font-size: 26px;
      font-weight: 760;
      margin-bottom: 4px;
    }}
    .metric-card span {{ color: var(--muted); font-size: 13px; }}
    .plot-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
      gap: 16px;
    }}
    .plot-card {{ margin: 0; padding: 16px; box-shadow: none; }}
    img {{
      display: block;
      width: 100%;
      max-width: 100%;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: white;
    }}
    .confusion {{
      display: grid;
      grid-template-columns: repeat(2, minmax(140px, 1fr));
      gap: 10px;
      max-width: 520px;
      margin-top: 14px;
    }}
    .confusion-cell {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
      background: var(--surface-alt);
    }}
    .confusion-cell strong {{ display: block; color: var(--muted); font-size: 13px; }}
    .confusion-cell span {{ font-size: 28px; font-weight: 760; }}
    .table-wrap {{ overflow-x: auto; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
    th, td {{ padding: 10px 12px; border-bottom: 1px solid var(--line); text-align: left; }}
    th {{ background: var(--surface-alt); color: var(--muted); font-size: 12px; text-transform: uppercase; letter-spacing: 0.04em; }}
    .notice {{ border-left: 5px solid var(--accent-2); }}
    .notice ul {{ margin: 0; padding-left: 20px; color: var(--muted); }}
    .two-col {{
      display: grid;
      grid-template-columns: minmax(260px, 0.8fr) minmax(360px, 1.2fr);
      gap: 20px;
      align-items: start;
    }}
    footer {{ color: var(--muted); font-size: 13px; padding: 8px 2px 28px; }}
    code {{ background: #e8edf3; padding: 0.1rem 0.25rem; border-radius: 4px; }}
    @media (max-width: 760px) {{
      header {{ padding: 26px 22px; }}
      main {{ padding: 16px; }}
      .two-col {{ grid-template-columns: 1fr; }}
      .plot-grid {{ grid-template-columns: 1fr; }}
      .metric-value {{ font-size: 22px; }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>Fragmentomics Model Evaluation Report</h1>
    <p>Research-use fragmentomics summary. The held-out test set is the primary estimate of model performance; cross-validation and full-cohort predictions are supporting diagnostics.</p>
  </header>
  <main>
    <section>
      <h2>Held-Out Test Performance</h2>
      <p class="muted">The model was selected and thresholded using the training partition, then evaluated on samples withheld from model fitting.</p>
      <div class="metrics">{''.join(heldout_cards)}</div>
    </section>

    <section class="two-col">
      <div>
        <h2>Held-Out Confusion Matrix</h2>
        <p class="muted">Calls are based on the selected threshold, optimized for approximately 95% specificity on cross-validated training scores.</p>
        <div class="confusion">{confusion_html}</div>
      </div>
      <div>
        <h2>Held-Out Misclassifications</h2>
        {records_table(misclassified_records, [
            ("sample_id", "Sample"),
            ("label", "Label"),
            ("error_type", "Error"),
            ("test_probability", "Probability"),
            ("test_classification", "Call"),
        ], "No held-out test misclassifications at the selected threshold.")}
      </div>
    </section>

    <section>
      <h2>Held-Out Test Plots</h2>
      <div class="plot-grid">
        {plot_card("Held-Out ROC Curve", "test_roc_curve.png", "Ranking performance on samples not used for fitting.")}
        {plot_card("Held-Out Precision-Recall Curve", "test_precision_recall_curve.png", "Precision-recall is especially informative when the cancer class is less common.")}
        {plot_card("Held-Out Score Distribution", "test_score_histogram.png", "Distribution of held-out probabilities by label.")}
      </div>
    </section>

    <section>
      <h2>Model And Cohort</h2>
      <div class="metrics">{''.join(model_cards)}</div>
    </section>

    <section>
      <h2>Cross-Validated Training Diagnostics</h2>
      <p class="muted">These metrics are from cross-validated predictions inside the training partition. They are useful for model selection and thresholding, but are secondary to the held-out test set.</p>
      <div class="metrics">{''.join(cv_cards)}</div>
      <div class="plot-grid" style="margin-top: 16px;">
        {plot_card("Cross-Validated ROC Curve", "cv_roc_curve.png")}
        {plot_card("Cross-Validated Precision-Recall Curve", "cv_precision_recall_curve.png")}
      </div>
    </section>

    <section>
      <h2>Full-Cohort Prediction Diagnostics</h2>
      <p class="muted">These plots summarize the supplied prediction file. If the same cohort was used for training, treat these as ranking and QC views rather than unbiased performance estimates.</p>
      <div class="metrics">{''.join(full_cards)}</div>
      <div class="plot-grid" style="margin-top: 16px;">
        {plot_card("Score Histogram", "score_histogram.png")}
        {plot_card("Scores By Label", "scores_by_label.png")}
        {plot_card("Ranked Scores", "ranked_scores.png")}
        {plot_card("Full-Cohort ROC Curve", "roc_curve.png")}
        {plot_card("Full-Cohort Precision-Recall Curve", "precision_recall_curve.png")}
      </div>
    </section>

    <section>
      <h2>Highest Full-Cohort Scores</h2>
      {records_table(top_records, [
          ("sample_id", "Sample"),
          ("label", "Label"),
          ("cancer_probability", "Probability"),
          ("classification", "Call"),
      ])}
    </section>

    {warnings_html}

    <footer>
      Output files: <code>prediction_analysis_summary.json</code>, <code>heldout_test_predictions.tsv</code>, <code>cross_validated_predictions.tsv</code>, <code>predictions_with_labels.tsv</code>.
    </footer>
  </main>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--model-metrics")
    ap.add_argument("--model-metadata")
    ap.add_argument("--cv-predictions")
    ap.add_argument("--test-predictions")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pred = pd.read_csv(args.predictions, sep="\t")
    matrix_cols = ["sample_id"]
    matrix_head = pd.read_csv(args.matrix, sep="\t", nrows=0)
    for col in ["label", "cancer_type", "batch", "sex", "age", "smoking_status"]:
        if col in matrix_head.columns:
            matrix_cols.append(col)
    meta = pd.read_csv(args.matrix, sep="\t", usecols=matrix_cols, low_memory=False)
    df = meta.merge(pred, on="sample_id", how="right")
    df.to_csv(outdir / "predictions_with_labels.tsv", sep="\t", index=False)

    ranked = df.sort_values("cancer_probability", ascending=False)
    ranked.to_csv(outdir / "ranked_predictions.tsv", sep="\t", index=False)

    summary = {
        "n_samples": int(len(df)),
        "n_predictions": int(pred["cancer_probability"].notna().sum()),
        "score_min": float(pred["cancer_probability"].min()),
        "score_median": float(pred["cancer_probability"].median()),
        "score_max": float(pred["cancer_probability"].max()),
        "unique_score_count": int(pred["cancer_probability"].nunique(dropna=True)),
        "classification_counts": pred["classification"].value_counts(dropna=False).to_dict(),
    }
    if "label" in df.columns and df["label"].notna().any():
        summary["evaluation_warning"] = (
            "These metrics are computed on the supplied prediction file. If this is the same cohort used for training, "
            "ROC/PR curves are resubstitution-style diagnostics and can be strongly optimistic."
        )
    if summary["unique_score_count"] <= 2:
        summary["score_warning"] = (
            "Prediction scores have very few unique values. ROC/PR curves may look deceptively clean because ranking "
            "collapses into a small number of score levels."
        )
    if pred["classification"].nunique(dropna=False) == 1:
        summary["classification_warning"] = (
            "All samples received the same classification. Inspect threshold_used and score ranking; "
            "this can happen when predicting on the training cohort or when many scores tie at the threshold."
        )

    labeled = "label" in df.columns and df["label"].notna().any()
    if labeled:
        y_true = parse_labels(df["label"])
        y_pred = df["classification"].eq("elevated_fragmentome_score")
        tp = int((y_true.eq(1) & y_pred).sum())
        tn = int((y_true.eq(0) & ~y_pred).sum())
        fp = int((y_true.eq(0) & y_pred).sum())
        fn = int((y_true.eq(1) & ~y_pred).sum())
        summary.update({
            "class_counts": {str(k): int(v) for k, v in y_true.value_counts().sort_index().items()},
            "confusion_matrix": {"tp": tp, "tn": tn, "fp": fp, "fn": fn},
            "sensitivity": tp / max(tp + fn, 1),
            "specificity": tn / max(tn + fp, 1),
        })
        if y_true.nunique() == 2:
            try:
                from sklearn.metrics import average_precision_score, roc_auc_score
                summary["roc_auc"] = float(roc_auc_score(y_true, df["cancer_probability"]))
                summary["pr_auc"] = float(average_precision_score(y_true, df["cancer_probability"]))
            except Exception as exc:
                summary["metric_warning"] = str(exc)

    if args.model_metrics and Path(args.model_metrics).exists():
        summary["model_metrics"] = json.loads(Path(args.model_metrics).read_text())
    if args.model_metadata and Path(args.model_metadata).exists():
        summary["model_metadata"] = json.loads(Path(args.model_metadata).read_text())
    cv_df = None
    if args.cv_predictions and Path(args.cv_predictions).exists():
        cv_df = pd.read_csv(args.cv_predictions, sep="\t")
        cv_eval = cv_df.dropna(subset=["cv_probability"]) if "cv_probability" in cv_df.columns else cv_df
        if {"y_true", "cv_probability"}.issubset(cv_eval.columns) and cv_eval["y_true"].nunique() == 2:
            try:
                from sklearn.metrics import average_precision_score, roc_auc_score
                summary["cross_validated_predictions"] = {
                    "n_samples": int(len(cv_eval)),
                    "unique_score_count": int(cv_eval["cv_probability"].nunique(dropna=True)),
                    "roc_auc": float(roc_auc_score(cv_eval["y_true"], cv_eval["cv_probability"])),
                    "pr_auc": float(average_precision_score(cv_eval["y_true"], cv_eval["cv_probability"])),
                }
            except Exception as exc:
                summary["cv_metric_warning"] = str(exc)
        cv_df.to_csv(outdir / "cross_validated_predictions.tsv", sep="\t", index=False)
    test_df = None
    if args.test_predictions and Path(args.test_predictions).exists():
        test_df = pd.read_csv(args.test_predictions, sep="\t")
        if {"y_true", "test_probability"}.issubset(test_df.columns) and test_df["y_true"].nunique() == 2:
            try:
                from sklearn.metrics import average_precision_score, roc_auc_score
                summary["heldout_test_predictions"] = {
                    "n_samples": int(len(test_df)),
                    "unique_score_count": int(test_df["test_probability"].nunique(dropna=True)),
                    "roc_auc": float(roc_auc_score(test_df["y_true"], test_df["test_probability"])),
                    "pr_auc": float(average_precision_score(test_df["y_true"], test_df["test_probability"])),
                }
            except Exception as exc:
                summary["test_metric_warning"] = str(exc)
        test_df.to_csv(outdir / "heldout_test_predictions.tsv", sep="\t", index=False)

    Path(outdir / "prediction_analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    score_summary = df.groupby("label", dropna=False)["cancer_probability"].describe() if "label" in df.columns else pred["cancer_probability"].describe().to_frame().T
    score_summary.to_csv(outdir / "score_summary_by_label.tsv", sep="\t")

    if labeled:
        confusion = pd.DataFrame([summary["confusion_matrix"]])
        confusion.to_csv(outdir / "confusion_matrix.tsv", sep="\t", index=False)

    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 4.5))
    if "label" in df.columns:
        for label, sub in df.groupby("label", dropna=False):
            plt.hist(sub["cancer_probability"], bins=20, alpha=0.65, label=str(label))
        plt.legend(title="label")
    else:
        plt.hist(df["cancer_probability"], bins=20, alpha=0.8)
    plt.xlabel("Fragmentome model probability")
    plt.ylabel("Samples")
    plt.title("Prediction Score Distribution")
    plt.tight_layout()
    plt.savefig(outdir / "score_histogram.png", dpi=200)
    plt.close()

    if "label" in df.columns:
        plt.figure(figsize=(6, 4.5))
        df.boxplot(column="cancer_probability", by="label")
        plt.suptitle("")
        plt.title("Scores By Label")
        plt.ylabel("Fragmentome model probability")
        plt.tight_layout()
        plt.savefig(outdir / "scores_by_label.png", dpi=200)
        plt.close()
    else:
        write_empty_plot(outdir / "scores_by_label.png", "Scores By Label")

    plt.figure(figsize=(9, 4.5))
    plot_df = ranked.reset_index(drop=True)
    colors = None
    if "label" in plot_df.columns:
        y_true = parse_labels(plot_df["label"])
        colors = np.where(y_true == 1, "#c44e52", "#4c72b0")
    plt.bar(np.arange(len(plot_df)), plot_df["cancer_probability"], color=colors)
    plt.xlabel("Samples ranked by score")
    plt.ylabel("Fragmentome model probability")
    plt.title("Ranked Prediction Scores")
    plt.tight_layout()
    plt.savefig(outdir / "ranked_scores.png", dpi=200)
    plt.close()

    if labeled and y_true.nunique() == 2:
        from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay
        RocCurveDisplay.from_predictions(y_true, df["cancer_probability"])
        plt.tight_layout()
        plt.savefig(outdir / "roc_curve.png", dpi=200)
        plt.close()
        PrecisionRecallDisplay.from_predictions(y_true, df["cancer_probability"])
        plt.tight_layout()
        plt.savefig(outdir / "precision_recall_curve.png", dpi=200)
        plt.close()
    else:
        write_empty_plot(outdir / "roc_curve.png", "ROC Curve")
        write_empty_plot(outdir / "precision_recall_curve.png", "Precision-Recall Curve")

    if cv_df is not None and {"y_true", "cv_probability"}.issubset(cv_df.columns):
        from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay
        cv_eval = cv_df.dropna(subset=["cv_probability"])
        if cv_eval["y_true"].nunique() == 2:
            RocCurveDisplay.from_predictions(cv_eval["y_true"], cv_eval["cv_probability"])
            plt.title("Cross-Validated ROC Curve")
            plt.tight_layout()
            plt.savefig(outdir / "cv_roc_curve.png", dpi=200)
            plt.close()
            PrecisionRecallDisplay.from_predictions(cv_eval["y_true"], cv_eval["cv_probability"])
            plt.title("Cross-Validated Precision-Recall Curve")
            plt.tight_layout()
            plt.savefig(outdir / "cv_precision_recall_curve.png", dpi=200)
            plt.close()
        else:
            write_empty_plot(outdir / "cv_roc_curve.png", "Cross-Validated ROC Curve")
            write_empty_plot(outdir / "cv_precision_recall_curve.png", "Cross-Validated Precision-Recall Curve")
    else:
        write_empty_plot(outdir / "cv_roc_curve.png", "Cross-Validated ROC Curve")
        write_empty_plot(outdir / "cv_precision_recall_curve.png", "Cross-Validated Precision-Recall Curve")

    if test_df is not None and {"y_true", "test_probability"}.issubset(test_df.columns) and test_df["y_true"].nunique() == 2:
        from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay
        RocCurveDisplay.from_predictions(test_df["y_true"], test_df["test_probability"])
        plt.title("Held-Out Test ROC Curve")
        plt.tight_layout()
        plt.savefig(outdir / "test_roc_curve.png", dpi=200)
        plt.close()
        PrecisionRecallDisplay.from_predictions(test_df["y_true"], test_df["test_probability"])
        plt.title("Held-Out Test Precision-Recall Curve")
        plt.tight_layout()
        plt.savefig(outdir / "test_precision_recall_curve.png", dpi=200)
        plt.close()
        plt.figure(figsize=(7, 4.5))
        for label, sub in test_df.groupby("label", dropna=False):
            plt.hist(sub["test_probability"], bins=10, alpha=0.65, label=str(label))
        plt.xlabel("Held-out test probability")
        plt.ylabel("Samples")
        plt.title("Held-Out Test Score Distribution")
        plt.legend(title="label")
        plt.tight_layout()
        plt.savefig(outdir / "test_score_histogram.png", dpi=200)
        plt.close()
    else:
        write_empty_plot(outdir / "test_roc_curve.png", "Held-Out Test ROC Curve")
        write_empty_plot(outdir / "test_precision_recall_curve.png", "Held-Out Test Precision-Recall Curve")
        write_empty_plot(outdir / "test_score_histogram.png", "Held-Out Test Score Distribution")

    html_report = build_html_report(summary, df, cv_df, test_df)
    Path(outdir / "prediction_analysis_report.html").write_text(html_report)


if __name__ == "__main__":
    main()

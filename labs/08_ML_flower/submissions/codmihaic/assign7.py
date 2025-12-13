import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.feature_selection import VarianceThreshold


# =========================
# CONFIG
# =========================
HANDLE = "codmihaic"
INPUT_CSV = f"labs/08_ML_flower/submissions/{HANDLE}/expression_matrix_{HANDLE}.csv"  
OUT_DIR = f"labs/08_ML_flower/submissions/{HANDLE}/outputs"

RANDOM_STATE = 42
TEST_SIZE = 0.25
KMEANS_K = 3
SCALE_BEFORE_PCA = True


# =========================
# HELPERS
# =========================
def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)

def plot_confusion_matrix(cm, class_names, out_path):
    plt.figure(figsize=(6, 5))
    plt.imshow(cm, interpolation="nearest")
    plt.title("Random Forest – Confusion Matrix")
    plt.colorbar()
    tick_marks = np.arange(len(class_names))
    plt.xticks(tick_marks, class_names, rotation=45, ha="right")
    plt.yticks(tick_marks, class_names)

    thresh = cm.max() / 2.0 if cm.max() != 0 else 0.5
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(
                j, i, format(cm[i, j], "d"),
                ha="center", va="center",
                color="white" if cm[i, j] > thresh else "black"
            )

    plt.ylabel("True")
    plt.xlabel("Predicted")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def plot_pca_scatter(X_pca, clusters, out_path):
    plt.figure(figsize=(6, 5))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=clusters)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Unsupervised: PCA (2D) + KMeans clusters")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


# =========================
# MAIN
# =========================
def main():
    ensure_outdir(OUT_DIR)

    # 1) Load data
    df = pd.read_csv(INPUT_CSV)

    # last column is the Label
    X = df.iloc[:, :-1]
    y = df.iloc[:, -1]

    # 2) Encode labels
    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    # 3) Stratified train/test
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc
    )

    # =========================
    # TASK 2 — Random Forest (Supervised)
    # =========================
    rf = RandomForestClassifier(
        n_estimators=300,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    # (a) classification_report_<handle>.txt
    report_txt = classification_report(y_test, y_pred, target_names=le.classes_)
    report_path = os.path.join(OUT_DIR, f"classification_report_{HANDLE}.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_txt)

    # (b) confusion_rf_<handle>.png
    cm = confusion_matrix(y_test, y_pred)
    cm_path = os.path.join(OUT_DIR, f"confusion_rf_{HANDLE}.png")
    plot_confusion_matrix(cm, le.classes_, cm_path)

    # (c) feature_importance_<handle>.csv
    importances = pd.DataFrame({
        "gene": X.columns,
        "importance": rf.feature_importances_
    }).sort_values("importance", ascending=False)

    fi_path = os.path.join(OUT_DIR, f"feature_importance_{HANDLE}.csv")
    importances.to_csv(fi_path, index=False)

    # =========================
    # TASK 3 — Logistic Regression
    # =========================
    lr = Pipeline([
        ("scaler", StandardScaler()),
        ("logreg", LogisticRegression(max_iter=2000, random_state=RANDOM_STATE))
    ])

    lr.fit(X_train, y_train)
    y_pred_lr = lr.predict(X_test)

    # comparation
    rf_acc = accuracy_score(y_test, y_pred)
    lr_acc = accuracy_score(y_test, y_pred_lr)
    rf_f1 = f1_score(y_test, y_pred, average="macro")
    lr_f1 = f1_score(y_test, y_pred_lr, average="macro")

    print("\n=== Task 3: Logistic Regression vs Random Forest ===")
    print(f"RF  - accuracy: {rf_acc:.4f}, macro-F1: {rf_f1:.4f}")
    print(f"LR  - accuracy: {lr_acc:.4f}, macro-F1: {lr_f1:.4f}")

    # =========================
    # TASK 4 — Unsupervised: PCA + KMeans
    # =========================
    X_uns = X_train.values

    if SCALE_BEFORE_PCA:
        scaler = StandardScaler()
        X_uns = scaler.fit_transform(X_uns)

    pca = PCA(n_components=2, random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(X_uns)

    kmeans = KMeans(n_clusters=KMEANS_K, random_state=RANDOM_STATE, n_init="auto")
    clusters = kmeans.fit_predict(X_pca)

    # (d) sup_vs_unsup_scatter_<handle>.png
    scatter_path = os.path.join(OUT_DIR, f"sup_vs_unsup_scatter_{HANDLE}.png")
    plot_pca_scatter(X_pca, clusters, scatter_path)

    # (e) cluster_crosstab_<handle>.csv
    ct = pd.crosstab(
        pd.Series(y_train, name="LabelEncoded"),
        pd.Series(clusters, name="Cluster")
    )
    ct.index = [le.inverse_transform([i])[0] for i in ct.index]
    crosstab_path = os.path.join(OUT_DIR, f"cluster_crosstab_{HANDLE}.csv")
    ct.to_csv(crosstab_path, index=True)

    # =========================
    # TASK 5 — Semi-Supervised
    # =========================
    rng = np.random.RandomState(RANDOM_STATE)
    missing_ratio = 0.4
    mask_unknown = rng.rand(len(y_train)) < missing_ratio 
    X_labeled = X_train.loc[~mask_unknown]
    y_labeled = y_train[~mask_unknown]

    X_unlabeled = X_train.loc[mask_unknown]
    rf_semi_1 = RandomForestClassifier(
        n_estimators=300,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )

    rf_semi_1.fit(X_labeled, y_labeled)
    y_pred_semi_1 = rf_semi_1.predict(X_test)
    acc_semi_1 = accuracy_score(y_test, y_pred_semi_1)
    f1_semi_1 = f1_score(y_test, y_pred_semi_1, average="macro")
    pseudo_labels = rf_semi_1.predict(X_unlabeled)
    X_full = pd.concat([X_labeled, X_unlabeled], axis=0)
    y_full = np.concatenate([y_labeled, pseudo_labels], axis=0)

    rf_semi_2 = RandomForestClassifier(
        n_estimators=300,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )
    rf_semi_2.fit(X_full, y_full)
    y_pred_semi_2 = rf_semi_2.predict(X_test)
    acc_semi_2 = accuracy_score(y_test, y_pred_semi_2)
    f1_semi_2 = f1_score(y_test, y_pred_semi_2, average="macro")

    print("\n=== Task 5: Semi-supervised (Pseudo-labeling) ===")
    print(f"Unknown ratio: {missing_ratio:.0%} (train labels hidden)")
    print(f"RF trained only on labeled   -> acc: {acc_semi_1:.4f}, macro-F1: {f1_semi_1:.4f}")
    print(f"RF retrained w/ pseudo-label -> acc: {acc_semi_2:.4f}, macro-F1: {f1_semi_2:.4f}")


    print("✅ Done. Generated deliverables in:", OUT_DIR)
    print("-", os.path.basename(report_path))
    print("-", os.path.basename(cm_path))
    print("-", os.path.basename(fi_path))
    print("-", os.path.basename(crosstab_path))
    print("-", os.path.basename(scatter_path))

    # =========================
    # BONUS (+1p) — PCA before vs after low-variance gene removal
    # =========================
    scaler_bonus = StandardScaler()
    X_scaled = scaler_bonus.fit_transform(X_train)
    pca_before = PCA(n_components=2, random_state=RANDOM_STATE)
    X_pca_before = pca_before.fit_transform(X_scaled)

    variances = X_train.var(axis=0)
    low_var_genes = variances.sort_values().head(10).index.tolist()
    X_filtered = X_train.drop(columns=low_var_genes)
    X_filt_scaled = scaler_bonus.fit_transform(X_filtered)
    pca_after = PCA(n_components=2, random_state=RANDOM_STATE)
    X_pca_after = pca_after.fit_transform(X_filt_scaled)

    # --- Plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(X_pca_before[:, 0], X_pca_before[:, 1], c=clusters)
    axes[0].set_title("PCA before low-variance filtering")
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")

    axes[1].scatter(X_pca_after[:, 0], X_pca_after[:, 1], c=clusters)
    axes[1].set_title("PCA after removing 10 low-variance genes")
    axes[1].set_xlabel("PC1")
    axes[1].set_ylabel("PC2")

    plt.tight_layout()
    bonus_path = os.path.join(OUT_DIR, f"bonus_pca_comparison_{HANDLE}.png")
    plt.savefig(bonus_path, dpi=200)
    plt.close()

    print("🎁 Bonus PCA comparison saved to:", bonus_path)


if __name__ == "__main__":
    main()

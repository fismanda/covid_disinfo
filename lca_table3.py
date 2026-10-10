#!/usr/bin/env python3
"""
Latent class analysis of COVID-19 misinformation beliefs (4-class model)
with posterior-weighted baseline characteristics and Wald tests.

Fits a 4-class LCA via EM algorithm (100 random starts) to four common
misinformation items from the EKOS COVID-19 Monitor, Wave 27:
  A. COVID-19 deaths exaggerated
  B. Vaccine deaths hidden by government
  C. Vaccines cause infertility
  D. Vaccines alter DNA

Item coding (all are false statements):
  0 = Correct (disagree/strongly disagree)
  1 = Don't Know
  2 = Incorrect (agree/strongly agree)

Outputs a CSV of posterior-weighted proportions, means (SD), and
Wald test p-values (chi-square for categorical, F for continuous).

Usage:
    python lca_table3.py path/to/05021W27_client_data_file.sav [--output table3.csv]

Dependencies: numpy, pandas, pyreadstat, scipy
"""

import argparse
import numpy as np
import pandas as pd
import pyreadstat
from scipy import stats


# ---------------------------------------------------------------------------
# 1. Data loading and recoding
# ---------------------------------------------------------------------------

def load_and_recode(sav_path):
    """Load Wave 27 SPSS file and recode the four common items."""
    df, meta = pyreadstat.read_sav(sav_path)
    print(f"Wave 27 raw: {df.shape[0]} rows, {df.shape[1]} columns")

    def recode_false(x):
        """Recode false-statement items: 1,2->Correct; 98->DK; 3,4->Incorrect."""
        if x in (1, 2):
            return 0
        elif x == 98:
            return 1
        elif x in (3, 4):
            return 2
        else:
            return np.nan

    item_cols_raw = ["Q34A", "Q34B", "Q34C", "Q34D"]
    for col in item_cols_raw:
        df[col + "_r"] = df[col].apply(recode_false)

    item_cols = [c + "_r" for c in item_cols_raw]
    df = df.dropna(subset=item_cols).copy()
    df[item_cols] = df[item_cols].astype(int)
    print(f"Complete on 4 items: N = {len(df)}")
    return df, item_cols


# ---------------------------------------------------------------------------
# 2. EM algorithm for K-class LCA with 3-category items
# ---------------------------------------------------------------------------

def em_lca(X, K, n_starts=100, max_iter=500, tol=1e-8):
    """
    Fit a K-class latent class model via EM with multiple random starts.

    Parameters
    ----------
    X : ndarray, shape (N, J)
        Item response matrix (integer-coded, 0..C-1).
    K : int
        Number of latent classes.
    n_starts : int
        Number of random initialisations.
    max_iter : int
        Maximum EM iterations per start.
    tol : float
        Log-likelihood convergence tolerance.

    Returns
    -------
    pi : ndarray, shape (K,)
        Class prevalences.
    rho : ndarray, shape (K, J, C)
        Item-response probabilities.
    post : ndarray, shape (N, K)
        Posterior class-membership probabilities.
    ll : float
        Final log-likelihood.
    """
    N, J = X.shape
    C = int(X.max()) + 1
    best_ll = -np.inf
    best_pi = best_rho = best_post = None

    for start in range(n_starts):
        rng = np.random.RandomState(start * 42 + 7)
        pi = rng.dirichlet(np.ones(K))
        rho = np.zeros((K, J, C))
        for k in range(K):
            for j in range(J):
                rho[k, j, :] = rng.dirichlet(np.ones(C) * 2)

        ll_prev = -np.inf
        for _ in range(max_iter):
            # E-step
            log_lik = np.tile(np.log(pi), (N, 1))
            for j in range(J):
                for k in range(K):
                    log_lik[:, k] += np.log(rho[k, j, X[:, j]] + 1e-300)

            max_ll = log_lik.max(axis=1, keepdims=True)
            log_denom = max_ll + np.log(
                np.exp(log_lik - max_ll).sum(axis=1, keepdims=True)
            )
            post = np.exp(log_lik - log_denom)
            ll = np.sum(log_denom)

            if abs(ll - ll_prev) < tol:
                break
            ll_prev = ll

            # M-step
            Nk = post.sum(axis=0)
            pi = Nk / N
            for k in range(K):
                for j in range(J):
                    for c in range(C):
                        rho[k, j, c] = (
                            (post[:, k] * (X[:, j] == c)).sum() / (Nk[k] + 1e-300)
                        )

        if ll > best_ll:
            best_ll = ll
            best_pi = pi.copy()
            best_rho = rho.copy()

    # Final E-step with best parameters
    log_lik = np.tile(np.log(best_pi), (N, 1))
    for j in range(J):
        for k in range(K):
            log_lik[:, k] += np.log(best_rho[k, j, X[:, j]] + 1e-300)
    max_ll = log_lik.max(axis=1, keepdims=True)
    log_denom = max_ll + np.log(
        np.exp(log_lik - max_ll).sum(axis=1, keepdims=True)
    )
    best_post = np.exp(log_lik - log_denom)

    return best_pi, best_rho, best_post, best_ll


# ---------------------------------------------------------------------------
# 3. Auto-label classes
# ---------------------------------------------------------------------------

def label_classes(rho, K):
    """
    Assign interpretive labels based on average item-response probabilities.
    Falls back to rank-ordering by P(Correct) if labels collide.
    """
    labels = []
    for k in range(K):
        avg_correct = rho[k, :, 0].mean()
        avg_dk = rho[k, :, 1].mean()
        avg_incorrect = rho[k, :, 2].mean()
        if avg_correct > 0.7:
            labels.append("Well-Informed")
        elif avg_incorrect > 0.4:
            labels.append("Misinformed")
        elif avg_dk > 0.3:
            labels.append("Uncertain/Disengaged")
        else:
            labels.append("Selectively Skeptical")

    # Resolve duplicates by rank-ordering on P(Correct)
    if len(set(labels)) < K:
        order = np.argsort([-rho[k, :, 0].mean() for k in range(K)])
        manual = [
            "Well-Informed",
            "Selectively Skeptical",
            "Uncertain/Disengaged",
            "Misinformed",
        ]
        labels = [""] * K
        for rank, k in enumerate(order):
            labels[k] = manual[rank]

    return labels


# ---------------------------------------------------------------------------
# 4. Posterior-weighted Wald tests
# ---------------------------------------------------------------------------

def pw_chisq_test(binary_series, post_matrix, labels, class_order):
    """
    Posterior-weighted Wald chi-square test for a binary variable.

    W = sum_k N_k (p_k - p_bar)^2 / [p_bar (1 - p_bar)]  ~  chi2(K-1)
    """
    K = len(class_order)
    valid = binary_series.notna()
    x = binary_series[valid].values.astype(float)
    post_v = post_matrix[valid.values]

    Nk = post_v.sum(axis=0)
    N_total = Nk.sum()

    pk = np.array(
        [
            (post_v[:, labels.index(c)] * x).sum() / Nk[labels.index(c)]
            for c in class_order
        ]
    )
    p_bar = (post_v.sum(axis=1) * x).sum() / N_total

    if p_bar <= 0 or p_bar >= 1:
        return np.nan, np.nan

    Nk_ordered = np.array([Nk[labels.index(c)] for c in class_order])
    W = np.sum(Nk_ordered * (pk - p_bar) ** 2) / (p_bar * (1 - p_bar))
    p_val = 1 - stats.chi2.cdf(W, K - 1)
    return W, p_val


def pw_chisq_multi_test(cat_series, categories, post_matrix, labels, class_order):
    """
    Posterior-weighted chi-square for a multi-category variable.

    W = sum_c sum_k N_k (p_kc - p_c)^2 / p_c   ~  chi2((K-1)(C-1))
    """
    K = len(class_order)
    valid = cat_series.notna()
    x = cat_series[valid].values
    post_v = post_matrix[valid.values]

    Nk = post_v.sum(axis=0)
    N_total = Nk.sum()
    C_cat = len(categories)

    W = 0
    for cat_val in categories:
        indicator = (x == cat_val).astype(float)
        pk = np.array(
            [
                (post_v[:, labels.index(c)] * indicator).sum() / Nk[labels.index(c)]
                for c in class_order
            ]
        )
        p_bar = (post_v.sum(axis=1) * indicator).sum() / N_total
        if p_bar > 0:
            Nk_ordered = np.array([Nk[labels.index(c)] for c in class_order])
            W += np.sum(Nk_ordered * (pk - p_bar) ** 2) / p_bar

    df = (K - 1) * (C_cat - 1)
    p_val = 1 - stats.chi2.cdf(W, df)
    return W, p_val


def pw_F_test(cont_series, post_matrix, labels, class_order):
    """
    Posterior-weighted one-way ANOVA F-test for a continuous variable.

    F = MSB / MSW, where MSB and MSW use posterior weights in place of
    hard class assignments.
    """
    K = len(class_order)
    valid = cont_series.notna()
    x = cont_series[valid].values.astype(float)
    post_v = post_matrix[valid.values]

    Nk = post_v.sum(axis=0)
    N_total = Nk.sum()

    mu_k = np.array(
        [
            (post_v[:, labels.index(c)] * x).sum() / Nk[labels.index(c)]
            for c in class_order
        ]
    )
    mu_bar = (post_v.sum(axis=1) * x).sum() / N_total

    Nk_ordered = np.array([Nk[labels.index(c)] for c in class_order])
    SSB = np.sum(Nk_ordered * (mu_k - mu_bar) ** 2)
    dfB = K - 1

    SSW = 0
    for i_c, c in enumerate(class_order):
        k = labels.index(c)
        SSW += np.sum(post_v[:, k] * (x - mu_k[i_c]) ** 2)
    dfW = N_total - K

    if dfW <= 0 or SSW <= 0:
        return np.nan, np.nan

    F = (SSB / dfB) / (SSW / dfW)
    p_val = 1 - stats.f.cdf(F, dfB, dfW)
    return F, p_val


# ---------------------------------------------------------------------------
# 5. Posterior-weighted descriptives (helpers)
# ---------------------------------------------------------------------------

def pw_mean(var_series, weights):
    valid = ~var_series.isna()
    if valid.sum() == 0:
        return np.nan, 0
    return np.average(var_series[valid], weights=weights[valid]), int(valid.sum())


def pw_prop(condition, weights):
    valid = ~condition.isna()
    if valid.sum() == 0:
        return np.nan, 0
    return np.average(condition[valid].astype(float), weights=weights[valid]), int(
        valid.sum()
    )


def pw_sd(var_series, weights):
    valid = ~var_series.isna()
    if valid.sum() == 0:
        return np.nan
    v = var_series[valid].values
    w = weights[valid].values
    mean = np.average(v, weights=w)
    variance = np.average((v - mean) ** 2, weights=w)
    return np.sqrt(variance)


def fmt_p(p):
    if np.isnan(p):
        return ""
    if p < 0.001:
        return "<0.001"
    elif p < 0.01:
        return f"{p:.3f}"
    else:
        return f"{p:.2f}"


# ---------------------------------------------------------------------------
# 6. Compute class-specific estimates and tests
# ---------------------------------------------------------------------------

def compute_table(w27, pi, rho, post, labels, class_order, K):
    """Build the Table 3 DataFrame with p-values."""

    results = []
    p_values = {}
    N = len(w27)

    for label in class_order:
        k_idx = labels.index(label)
        weights = w27[f"post_{k_idx}"]
        row = {"Class": label, "N": f"{pi[k_idx]*N:.0f} ({pi[k_idx]*100:.1f}%)"}

        # Sex
        if "QGENDR" in w27.columns:
            male = (w27["QGENDR"] == 1).where(w27["QGENDR"].notna())
            m, _ = pw_prop(male, weights)
            row["Male (%)"] = f"{m*100:.1f}"
            female = (w27["QGENDR"] == 2).where(w27["QGENDR"].notna())
            m, _ = pw_prop(female, weights)
            row["Female (%)"] = f"{m*100:.1f}"

        # Age groups
        if "QAGEMIX" in w27.columns:
            for ages, lbl in [
                ([1, 2], "Age 18-34 (%)"),
                ([3, 4], "Age 35-54 (%)"),
                ([5], "Age 55-64 (%)"),
                ([6], "Age >= 65 (%)"),
            ]:
                cond = w27["QAGEMIX"].isin(ages).where(w27["QAGEMIX"].notna())
                m, _ = pw_prop(cond, weights)
                row[lbl] = f"{m*100:.1f}"

        # Education
        if "QEDUC" in w27.columns:
            univ = w27["QEDUC"].isin([7, 8]).where(w27["QEDUC"].notna())
            m, _ = pw_prop(univ, weights)
            row["Bachelor's or postgraduate (%)"] = f"{m*100:.1f}"

        # Rural
        if "QRURBAN" in w27.columns:
            rural = (w27["QRURBAN"] == 1).where(w27["QRURBAN"].notna())
            m, _ = pw_prop(rural, weights)
            row["Rural community (%)"] = f"{m*100:.1f}"

        # Social class
        if "SOCIALC1" in w27.columns:
            poor_wc = w27["SOCIALC1"].isin([1, 2]).where(w27["SOCIALC1"].notna())
            m, _ = pw_prop(poor_wc, weights)
            row["Poor or working social class (%)"] = f"{m*100:.1f}"

        # Vaccination
        if "Q15" in w27.columns:
            unvax = (w27["Q15"] == 4).where(w27["Q15"].notna())
            m, _ = pw_prop(unvax, weights)
            row["Unvaccinated (%)"] = f"{m*100:.1f}"
            vax3 = (w27["Q15"] == 3).where(w27["Q15"].notna())
            m, _ = pw_prop(vax3, weights)
            row["Three or more vaccine doses (%)"] = f"{m*100:.1f}"

        # Trust (1-7 scale)
        trust_map = {
            "Q3A": "Trust in doctors",
            "Q3E": "Trust in scientists",
            "Q3F": "Trust in public health experts",
            "Q3H": "Trust in journalists",
        }
        for var, name in trust_map.items():
            if var in w27.columns:
                vals = w27[var].replace(99, np.nan)
                m, _ = pw_mean(vals, weights)
                s = pw_sd(vals, weights)
                row[f"{name} (1-7)"] = f"{m:.2f} ({s:.2f})"

        # Approval of PHAC
        if "Q10D" in w27.columns:
            vals = w27["Q10D"].replace(99, np.nan)
            m, _ = pw_mean(vals, weights)
            s = pw_sd(vals, weights)
            row["Approval of PHAC (1-5)"] = f"{m:.2f} ({s:.2f})"

        # Media reliance
        media_map = {
            "Q35F": "Reliance on health professionals",
            "Q35C": "Reliance on public health agencies",
            "Q35G": "Reliance on traditional media",
            "Q35A": "Reliance on social media",
        }
        for var, name in media_map.items():
            if var in w27.columns:
                vals = w27[var].replace(99, np.nan)
                m, _ = pw_mean(vals, weights)
                s = pw_sd(vals, weights)
                row[f"{name} (1-5)"] = f"{m:.2f} ({s:.2f})"

        # Concern about disinformation
        if "Q9" in w27.columns:
            vals = w27["Q9"].replace(99, np.nan)
            m, _ = pw_mean(vals, weights)
            s = pw_sd(vals, weights)
            row["Concern about disinformation (1-5)"] = f"{m:.2f} ({s:.2f})"

        results.append(row)

    # --- Statistical tests ---
    print("\n=== Statistical tests (posterior-weighted Wald) ===\n")

    # Sex (binary)
    male_bin = (w27["QGENDR"] == 1).where(w27["QGENDR"].notna())
    W, p = pw_chisq_test(male_bin, post, labels, class_order)
    p_values["Sex"] = p
    print(f"Sex: W={W:.2f}, p={fmt_p(p)}")

    # Age (4-category)
    age_grp = w27["QAGEMIX"].copy().where(w27["QAGEMIX"].notna())
    age_grp = age_grp.map(
        lambda x: (
            1
            if x in (1, 2)
            else (2 if x in (3, 4) else (3 if x == 5 else (4 if x == 6 else np.nan)))
        )
        if pd.notna(x)
        else np.nan
    )
    W, p = pw_chisq_multi_test(age_grp, [1, 2, 3, 4], post, labels, class_order)
    p_values["Age"] = p
    print(f"Age group: W={W:.2f}, p={fmt_p(p)}")

    # Education (binary)
    univ = w27["QEDUC"].isin([7, 8]).where(w27["QEDUC"].notna())
    W, p = pw_chisq_test(univ, post, labels, class_order)
    p_values["Education"] = p
    print(f"Education: W={W:.2f}, p={fmt_p(p)}")

    # Rural (binary)
    rural = (w27["QRURBAN"] == 1).where(w27["QRURBAN"].notna())
    W, p = pw_chisq_test(rural, post, labels, class_order)
    p_values["Rural"] = p
    print(f"Rural: W={W:.2f}, p={fmt_p(p)}")

    # Social class (binary)
    poor_wc = w27["SOCIALC1"].isin([1, 2]).where(w27["SOCIALC1"].notna())
    W, p = pw_chisq_test(poor_wc, post, labels, class_order)
    p_values["Social class"] = p
    print(f"Social class: W={W:.2f}, p={fmt_p(p)}")

    # Vaccination (3-category)
    vax_valid = w27["Q15"].copy().where(w27["Q15"].isin([2, 3, 4]))
    W, p = pw_chisq_multi_test(vax_valid, [2, 3, 4], post, labels, class_order)
    p_values["Vaccination"] = p
    print(f"Vaccination: W={W:.2f}, p={fmt_p(p)}")

    # Trust variables (continuous, F-test)
    trust_vars = {
        "Q3A": "Trust: doctors",
        "Q3E": "Trust: scientists",
        "Q3F": "Trust: PH experts",
        "Q3H": "Trust: journalists",
    }
    for var, name in trust_vars.items():
        vals = w27[var].replace(99, np.nan)
        F, p = pw_F_test(vals, post, labels, class_order)
        p_values[name] = p
        print(f"{name}: F={F:.2f}, p={fmt_p(p)}")

    # Approval of PHAC
    vals = w27["Q10D"].replace(99, np.nan)
    F, p = pw_F_test(vals, post, labels, class_order)
    p_values["Approval PHAC"] = p
    print(f"Approval PHAC: F={F:.2f}, p={fmt_p(p)}")

    # Media reliance (continuous, F-test)
    media_vars = {
        "Q35F": "Reliance: health prof",
        "Q35C": "Reliance: PH agencies",
        "Q35G": "Reliance: trad media",
        "Q35A": "Reliance: social media",
    }
    for var, name in media_vars.items():
        vals = w27[var].replace(99, np.nan)
        F, p = pw_F_test(vals, post, labels, class_order)
        p_values[name] = p
        print(f"{name}: F={F:.2f}, p={fmt_p(p)}")

    # Concern about disinformation
    vals = w27["Q9"].replace(99, np.nan)
    F, p = pw_F_test(vals, post, labels, class_order)
    p_values["Concern disinfo"] = p
    print(f"Concern disinfo: F={F:.2f}, p={fmt_p(p)}")

    # --- Assemble DataFrame ---
    df = pd.DataFrame(results).set_index("Class").T

    pval_map = {
        "N": "",
        "Male (%)": p_values.get("Sex"),
        "Female (%)": "",
        "Age 18-34 (%)": p_values.get("Age"),
        "Age 35-54 (%)": "",
        "Age 55-64 (%)": "",
        "Age >= 65 (%)": "",
        "Bachelor's or postgraduate (%)": p_values.get("Education"),
        "Rural community (%)": p_values.get("Rural"),
        "Poor or working social class (%)": p_values.get("Social class"),
        "Unvaccinated (%)": p_values.get("Vaccination"),
        "Three or more vaccine doses (%)": "",
        "Trust in doctors (1-7)": p_values.get("Trust: doctors"),
        "Trust in scientists (1-7)": p_values.get("Trust: scientists"),
        "Trust in public health experts (1-7)": p_values.get("Trust: PH experts"),
        "Trust in journalists (1-7)": p_values.get("Trust: journalists"),
        "Approval of PHAC (1-5)": p_values.get("Approval PHAC"),
        "Reliance on health professionals (1-5)": p_values.get("Reliance: health prof"),
        "Reliance on public health agencies (1-5)": p_values.get(
            "Reliance: PH agencies"
        ),
        "Reliance on traditional media (1-5)": p_values.get("Reliance: trad media"),
        "Reliance on social media (1-5)": p_values.get("Reliance: social media"),
        "Concern about disinformation (1-5)": p_values.get("Concern disinfo"),
    }

    pval_col = []
    for idx in df.index:
        val = pval_map.get(idx, "")
        if isinstance(val, float) and not np.isnan(val):
            pval_col.append(fmt_p(val))
        else:
            pval_col.append("")
    df["P value"] = pval_col

    return df


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="4-class LCA of COVID-19 misinformation beliefs (Wave 27)"
    )
    parser.add_argument("sav_path", help="Path to Wave 27 SPSS (.sav) data file")
    parser.add_argument(
        "--output",
        "-o",
        default="table3_baseline_by_class.csv",
        help="Output CSV path (default: table3_baseline_by_class.csv)",
    )
    args = parser.parse_args()

    K = 4
    CLASS_ORDER = [
        "Well-Informed",
        "Selectively Skeptical",
        "Uncertain/Disengaged",
        "Misinformed",
    ]

    # Load and recode
    w27, item_cols = load_and_recode(args.sav_path)
    X = w27[item_cols].values

    # Fit LCA
    print("\nRunning EM with 100 random starts...")
    pi, rho, post, ll = em_lca(X, K)
    labels = label_classes(rho, K)
    print(f"Labels: {labels}")
    print(f"Class proportions: {np.round(pi, 3)}")
    print(f"Log-likelihood: {ll:.2f}")

    # Store posteriors in dataframe
    for k in range(K):
        w27[f"post_{k}"] = post[:, k]

    # Build table
    df = compute_table(w27, pi, rho, post, labels, CLASS_ORDER, K)

    print("\n" + "=" * 80)
    print("TABLE 3: Baseline Characteristics by Latent Class (Wave 27)")
    print("=" * 80)
    print(df.to_string())

    df.to_csv(args.output)
    print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()

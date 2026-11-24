"""DESeq2's default analysis for a design with one factor, vectorised with numpy.

The paper used DESeq2 1.44 in R. This module repeats its default steps so the
build runs in seconds and without R:

- size factors by median of ratios, or gene-by-sample normalisation factors from
  tximport's average transcript lengths (DESeq2's estimateNormFactors);
- gene-wise dispersions by Cox-Reid adjusted profile likelihood, with DESeq2's
  backtracking line search and its grid fallback;
- the parametric dispersion trend, the log-normal prior and MAP dispersions,
  keeping gene-wise values for dispersion outliers;
- IRLS with DESeq2's 1e-6 ridge, Wald tests for any pair of levels (the
  covariance uses the fitted means floored at 0.5, as DESeq2's fitBeta does);
- Cook's distances, flagged above the 0.99 quantile of F(p, m - p).

Independent filtering is a separate helper, because the paper adjusted the raw
p-values over all features instead (p.adjust in its own code).

Checked against the data by script (see the paper's README): on the authors'
ATAC count matrix it returns the published differentially accessible peaks
exactly (40, 50, 72 and 22), and on the RNA-seq matrix with the tximport
lengths it lands within 1% of the published gene counts.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MIN_DISP = 1e-8
MIN_MU = 0.5
KAPPA0 = 1.0
DISP_TOL = 1e-6
MAXIT = 100
ARMIJO = 1e-4
RIDGE = 1e-6  # DESeq2's betaPriorVar of 1e6 on the log2 scale, as a ridge


def size_factors(counts: np.ndarray) -> np.ndarray:
    """Median-of-ratios size factors (estimateSizeFactorsForMatrix)."""
    with np.errstate(divide="ignore"):
        logc = np.log(np.asarray(counts, dtype=float))
    loggeo = logc.mean(axis=1)
    ok = np.isfinite(loggeo)
    return np.exp(np.median(logc[ok] - loggeo[ok, None], axis=0))


def norm_factors_from_lengths(counts: np.ndarray, avg_tx_length: np.ndarray) -> np.ndarray:
    """DESeq2's normalisation factors when tximport's avgTxLength is present."""
    lengths = np.asarray(avg_tx_length, dtype=float)
    nm = lengths / np.exp(np.log(lengths).mean(1, keepdims=True))
    sf = size_factors(np.asarray(counts, dtype=float) / nm)
    nf = nm * sf[None, :]
    return nf / np.exp(np.log(nf).mean(1, keepdims=True))


def length_scaled_counts(counts: np.ndarray, avg_tx_length: np.ndarray) -> np.ndarray:
    """Counts with the per-sample length differences divided out, rounded.

    The alternative tximport offers when a tool cannot take offsets: each gene's
    lengths are scaled to a geometric mean of 1 across samples first.
    """
    lengths = np.asarray(avg_tx_length, dtype=float)
    nm = lengths / np.exp(np.log(lengths).mean(1, keepdims=True))
    return np.round(np.asarray(counts, dtype=float) / nm)


def _cell_sums(m: np.ndarray, gi: np.ndarray, k: int) -> np.ndarray:
    return np.stack([m[:, gi == g].sum(1) for g in range(k)], 1)


def _log_post(la, y, mu, gi, k, prior_mean=None, prior_var=None):
    """Cox-Reid adjusted NB log likelihood (plus log prior) at log dispersion la."""
    from scipy.special import gammaln

    a = np.exp(la)[:, None]
    an = 1.0 / a
    ll = (gammaln(y + an) - gammaln(an) - y * np.log(mu + an) - an * np.log1p(mu * a)).sum(1)
    w = 1.0 / (1.0 / mu + a)
    out = ll - 0.5 * np.log(_cell_sums(w, gi, k)).sum(1)
    if prior_mean is not None:
        out = out - 0.5 * (la - prior_mean) ** 2 / prior_var
    return out


def _dlog_post(la, y, mu, gi, k, prior_mean=None, prior_var=None):
    """Derivative of _log_post with respect to log dispersion."""
    from scipy.special import digamma

    a = np.exp(la)[:, None]
    an = 1.0 / a
    ll_part = (an ** 2 * (digamma(an) + np.log1p(mu * a) - mu * a / (1.0 + mu * a)
                          - digamma(y + an) + y / (mu + an))).sum(1)
    w = 1.0 / (1.0 / mu + a)
    dw = -1.0 / (1.0 / mu + a) ** 2
    cr = -0.5 * (_cell_sums(dw, gi, k) / _cell_sums(w, gi, k)).sum(1)
    out = (ll_part + cr) * a[:, 0]
    if prior_mean is not None:
        out = out - (la - prior_mean) / prior_var
    return out


def _line_search(la0, y, mu, gi, k, prior_mean=None, prior_var=None):
    """DESeq2's fitDisp: backtracking line search on log dispersion, all genes in lockstep."""
    min_log_alpha = np.log(MIN_DISP / 10)
    n_genes = la0.shape[0]
    a = la0.copy()

    def pm(ix):
        return None if prior_mean is None else prior_mean[ix]

    everyone = np.arange(n_genes)
    lp = _log_post(a, y, mu, gi, k, pm(everyone), prior_var)
    dlp = _dlog_post(a, y, mu, gi, k, pm(everyone), prior_var)
    initial_lp = lp.copy()
    kappa = np.full(n_genes, KAPPA0)
    iters = np.zeros(n_genes, dtype=int)
    accepts = np.zeros(n_genes, dtype=int)
    active = np.ones(n_genes, dtype=bool)
    for _ in range(MAXIT):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        iters[idx] += 1
        ai, dl, ka = a[idx], dlp[idx], kappa[idx]
        prop = ai + ka * dl
        ka = np.where(prop < -30.0, (-30.0 - ai) / dl, ka)
        ka = np.where(prop > 10.0, (10.0 - ai) / dl, ka)
        theta_k = -_log_post(ai + ka * dl, y[idx], mu[idx], gi, k, pm(idx), prior_var)
        theta_hat = -lp[idx] - ka * ARMIJO * dl ** 2
        acc = theta_k <= theta_hat
        kappa[idx[~acc]] = ka[~acc] / 2.0
        if not acc.any():
            continue
        ok = idx[acc]
        new_a = ai[acc] + ka[acc] * dl[acc]
        accepts[ok] += 1
        a[ok] = new_a
        lp_new = _log_post(new_a, y[ok], mu[ok], gi, k, pm(ok), prior_var)
        change = lp_new - lp[ok]
        lp[ok] = lp_new
        stop = (change < DISP_TOL) | (new_a < min_log_alpha)
        active[ok[stop]] = False
        cont = ok[~stop]
        if cont.size:
            dlp[cont] = _dlog_post(a[cont], y[cont], mu[cont], gi, k, pm(cont), prior_var)
            kn = np.minimum(ka[acc][~stop] * 1.1, KAPPA0)
            kappa[cont] = np.where(accepts[cont] % 5 == 0, kn / 2.0, kn)
    return a, iters, initial_lp, lp


def _grid(y, mu, gi, k, max_log_alpha, prior_mean=None, prior_var=None):
    """DESeq2's fitDispGrid: 20 points, then 20 finer points around the best."""
    n_genes = y.shape[0]
    grid = np.linspace(np.log(1e-8), max_log_alpha, 20)
    delta = grid[1] - grid[0]
    vals = np.stack([_log_post(np.full(n_genes, g), y, mu, gi, k, prior_mean, prior_var) for g in grid], 1)
    best = grid[np.argmax(vals, 1)]
    cand = best[:, None] + np.linspace(-1, 1, 20)[None, :] * delta
    vals = np.stack([_log_post(cand[:, j], y, mu, gi, k, prior_mean, prior_var) for j in range(20)], 1)
    return cand[np.arange(n_genes), np.argmax(vals, 1)]


def _parametric_trend(means, disps):
    """DESeq2's parametricDispersionFit: Gamma GLM (identity link) of disps on 1/means."""
    coefs = np.array([0.1, 1.0])
    for _outer in range(11):
        resid = disps / (coefs[0] + coefs[1] / means)
        good = (resid > 1e-4) & (resid < 15)
        x, yv = 1.0 / means[good], disps[good]
        design = np.stack([np.ones_like(x), x], 1)
        b = coefs.copy()
        dev_old = 2 * np.sum(-np.log(yv / (b[0] + b[1] * x)) + (yv - (b[0] + b[1] * x)) / (b[0] + b[1] * x))
        for _inner in range(100):
            fitted = b[0] + b[1] * x
            xtw = design.T * (1.0 / fitted ** 2)
            b = np.linalg.solve(xtw @ design, xtw @ yv)
            fitted = b[0] + b[1] * x
            dev = 2 * np.sum(-np.log(yv / fitted) + (yv - fitted) / fitted)
            if abs(dev - dev_old) / (abs(dev) + 0.1) < 1e-8:
                break
            dev_old = dev
        old, coefs = coefs, b
        if not np.all(coefs > 0):
            raise RuntimeError("parametric dispersion fit failed")
        if np.sum(np.log(coefs / old) ** 2) < 1e-6:
            break
    return coefs


def _mad(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def _trimmed_mean(a, trim):
    """R's mean(x, trim=) along axis 1."""
    n = a.shape[1]
    lo = int(np.floor(n * trim))
    return np.sort(a, axis=1)[:, lo:n - lo].mean(1)


@dataclass
class Fit:
    levels: list
    base_mean: np.ndarray
    disp: np.ndarray
    beta: np.ndarray  # features x levels, natural log, treatment coding against levels[0]
    sigma: np.ndarray  # features x levels x levels
    max_cooks: np.ndarray
    cooks_cutoff: float
    all_zero: np.ndarray

    def contrast(self, num: str, den: str) -> dict:
        """Wald test of level `num` against level `den`; log2 fold change is num over den."""
        from scipy.stats import norm

        c = np.zeros(len(self.levels))
        i, j = self.levels.index(num), self.levels.index(den)
        if i > 0:
            c[i] += 1
        if j > 0:
            c[j] -= 1
        lfc = (self.beta @ c) / np.log(2)
        se = np.sqrt(np.einsum("i,gij,j->g", c, self.sigma, c)) / np.log(2)
        stat = lfc / se
        pvalue = 2 * norm.sf(np.abs(stat))
        pvalue = np.where(self.max_cooks > self.cooks_cutoff, np.nan, pvalue)
        for arr in (lfc, se, stat, pvalue):
            arr[self.all_zero] = np.nan
        return {"baseMean": self.base_mean, "log2FoldChange": lfc, "lfcSE": se, "stat": stat,
                "pvalue": pvalue}


def fit(counts, groups, levels=None, norm_factors=None) -> Fit:
    """Run DESeq() with design ~group. counts: features x samples."""
    from scipy.special import polygamma
    from scipy.stats import f as fdist

    y_all = np.asarray(counts, dtype=float)
    groups = np.asarray(groups)
    levels = sorted(set(groups.tolist())) if levels is None else list(levels)
    k = len(levels)
    gi = np.array([levels.index(g) for g in groups])
    n = y_all.shape[1]
    X = np.zeros((n, k))
    X[:, 0] = 1
    for g in range(1, k):
        X[gi == g, g] = 1
    all_zero = y_all.sum(1) == 0
    y = y_all[~all_zero]
    if norm_factors is None:
        nf = np.broadcast_to(size_factors(y_all)[None, :], y.shape).copy()
    else:
        nf = np.asarray(norm_factors, dtype=float)[~all_zero]
    norm_c = y / nf
    base_mean = norm_c.mean(1)
    max_disp = max(10.0, n)

    # gene-wise dispersions, with means from the linear model on normalised counts
    gm = np.stack([norm_c[:, gi == g].mean(1) for g in range(k)], 1)
    mu = np.maximum(gm[:, gi] * nf, MIN_MU)
    mu_rough = np.maximum(1.0, gm[:, gi])
    rough = np.maximum((((norm_c - mu_rough) ** 2 - mu_rough) / mu_rough ** 2).sum(1) / (n - k), 0)
    xim = np.mean(1.0 / nf.mean(0))
    moments = (norm_c.var(1, ddof=1) - xim * base_mean) / base_mean ** 2
    alpha_init = np.clip(np.minimum(rough, moments), MIN_DISP, max_disp)
    la, iters, init_lp, last_lp = _line_search(np.log(alpha_init), y, mu, gi, k)
    disp_gene = np.minimum(np.exp(la), max_disp)
    no_increase = last_lp < init_lp + np.abs(init_lp) / 1e6
    disp_gene[no_increase] = alpha_init[no_increase]
    refit = ~((iters < MAXIT) & (iters != 1)) & (disp_gene > MIN_DISP * 10)
    if refit.any():
        disp_gene[refit] = np.exp(_grid(y[refit], mu[refit], gi, k, np.log(max_disp)))
    disp_gene = np.clip(disp_gene, MIN_DISP, max_disp)

    # trend, prior and MAP dispersions
    use = disp_gene > 100 * MIN_DISP
    coefs = _parametric_trend(base_mean[use], disp_gene[use])
    disp_fit = coefs[0] + coefs[1] / base_mean
    var_log = _mad((np.log(disp_gene) - np.log(disp_fit))[use]) ** 2
    prior_var = max(var_log - float(polygamma(1, (n - k) / 2)), 0.25)
    start = np.where(disp_gene > 0.1 * disp_fit, disp_gene, disp_fit)
    log_fit = np.log(disp_fit)
    la2, it2, _, _ = _line_search(np.log(start), y, mu, gi, k, log_fit, prior_var)
    disp_map = np.exp(la2)
    refit2 = ~(it2 < MAXIT)
    if refit2.any():
        disp_map[refit2] = np.exp(_grid(y[refit2], mu[refit2], gi, k, np.log(max_disp), log_fit[refit2], prior_var))
    disp_map = np.clip(disp_map, MIN_DISP, max_disp)
    outlier = np.log(disp_gene) > log_fit + 2 * np.sqrt(var_log)
    disp = np.where(outlier, disp_gene, disp_map)

    # GLM coefficients by IRLS with a ridge
    from scipy.special import gammaln

    ridge = np.diag(np.full(k, RIDGE / np.log(2) ** 2))
    beta = np.linalg.lstsq(X, np.log(norm_c + 0.1).T, rcond=None)[0].T
    a = disp[:, None]
    mu_hat = np.maximum(nf * np.exp(beta @ X.T), MIN_MU)
    dev_old = np.zeros(y.shape[0])
    active = np.ones(y.shape[0], dtype=bool)
    for t in range(MAXIT):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        m_ = mu_hat[idx]
        w = m_ / (1.0 + a[idx] * m_)
        z = np.log(m_ / nf[idx]) + (y[idx] - m_) / m_
        xtwx = np.einsum("ni,gn,nj->gij", X, w, X) + ridge
        xtwz = np.einsum("ni,gn->gi", X, w * z)
        b = np.linalg.solve(xtwx, xtwz[..., None])[..., 0]
        beta[idx] = b
        m_new = np.maximum(nf[idx] * np.exp(b @ X.T), MIN_MU)
        mu_hat[idx] = m_new
        r = 1.0 / a[idx]
        ll = (gammaln(y[idx] + r) - gammaln(r) - gammaln(y[idx] + 1) - r * np.log1p(m_new / r)
              + y[idx] * (np.log(m_new) - np.log(m_new + r)))
        dev = -2.0 * ll.sum(1)
        conv = np.abs(dev - dev_old[idx]) / (np.abs(dev) + 0.1)
        done = ((t > 0) & (conv < 1e-8)) | np.any(np.abs(b) > 30, axis=1) | np.isnan(conv)
        dev_old[idx] = dev
        active[idx[done]] = False
    mu_f = nf * np.exp(beta @ X.T)
    mu_w = np.maximum(mu_f, MIN_MU)
    w = mu_w / (1.0 + a * mu_w)
    xtwx = np.einsum("ni,gn,nj->gij", X, w, X)
    inv = np.linalg.inv(xtwx + ridge)
    sigma = inv @ xtwx @ inv
    hat = np.einsum("gn,ni,gij,nj->gn", w, X, inv, X)

    # Cook's distances (robust method-of-moments dispersions, cells of 3 or more)
    sizes = np.bincount(gi, minlength=k)
    cells = [g for g in range(k) if sizes[g] >= 3]
    variances = []
    for g in cells:
        cg = norm_c[:, gi == g]
        small = cg.shape[1] <= 3.5
        trim, scale = (1 / 3, 2.04) if small else (1 / 4, 1.86)
        centre = _trimmed_mean(cg, trim)
        variances.append(scale * _trimmed_mean((cg - centre[:, None]) ** 2, trim))
    v = np.max(np.stack(variances, 1), 1)
    alpha_cooks = np.maximum((v - base_mean) / base_mean ** 2, 0.04)
    pearson_sq = (y - mu_f) ** 2 / (mu_f + alpha_cooks[:, None] * mu_f ** 2)
    cooks = pearson_sq / k * hat / (1 - hat) ** 2
    max_cooks = cooks[:, np.isin(gi, cells)].max(1)

    def expand(arr, fill=np.nan):
        out = np.full((y_all.shape[0],) + arr.shape[1:], fill, dtype=float)
        out[~all_zero] = arr
        return out

    return Fit(levels=levels, base_mean=expand(base_mean, 0.0), disp=expand(disp), beta=expand(beta),
               sigma=expand(sigma), max_cooks=expand(max_cooks), cooks_cutoff=float(fdist.ppf(0.99, k, n - k)),
               all_zero=all_zero)


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg over the non-missing p-values (R's p.adjust, method 'fdr')."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    order = np.argsort(q, kind="stable")
    ranked = q[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def padj_independent_filtering(base_mean, pvalue, alpha=0.1) -> np.ndarray:
    """DESeq2 results() padj: BH after a baseMean filter chosen as DESeq2 chooses it."""
    from statsmodels.nonparametric.smoothers_lowess import lowess

    base_mean = np.asarray(base_mean, float)
    pvalue = np.asarray(pvalue, float)
    lower = np.mean(base_mean == 0)
    theta = np.linspace(lower, 0.95 if lower < 0.95 else 1.0, 50)
    padjs = []
    for cut in np.quantile(base_mean, theta):
        pa = np.full(pvalue.shape, np.nan)
        keep = base_mean >= cut
        pa[keep] = bh(pvalue[keep])
        padjs.append(pa)
    num_rej = np.array([np.nansum(pa < alpha) for pa in padjs])
    if num_rej.max() <= 10:
        return padjs[0]
    fitted = lowess(num_rej, theta, frac=1 / 5, it=3, delta=0.01 * (theta.max() - theta.min()),
                    return_sorted=False)
    residual = num_rej[num_rej > 0] - fitted[num_rej > 0]
    ok = num_rej > fitted.max() - np.sqrt(np.var(residual, ddof=1))
    return padjs[int(np.flatnonzero(ok)[0]) if ok.any() else 0]


def filter_by_expr(counts, groups, min_count=10, min_total_count=15, large_n=10, min_prop=0.7) -> np.ndarray:
    """edgeR's filterByExpr with a group design."""
    y = np.asarray(counts, dtype=float)
    lib = y.sum(0)
    _, sizes = np.unique(np.asarray(groups), return_counts=True)
    m = sizes.min()
    if m > large_n:
        m = large_n + (m - large_n) * min_prop
    cutoff = min_count / np.median(lib) * 1e6
    cpm = y / lib * 1e6
    return ((cpm >= cutoff).sum(1) >= m - 1e-14) & (y.sum(1) >= min_total_count - 1e-14)

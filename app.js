window.ACCT_APP='scales';
/* ===== Psychometrics engine (pure functions) ===== */
const PSY = (() => {
  const sum = a => a.reduce((s, v) => s + v, 0);
  const mean = a => a.length ? sum(a) / a.length : NaN;
  const variance = a => { const n = a.length; if (n < 2) return NaN; const m = mean(a); return sum(a.map(v => (v - m) ** 2)) / (n - 1); };
  const sd = a => Math.sqrt(variance(a));
  const cov = (x, y) => { const n = x.length, mx = mean(x), my = mean(y); let s = 0; for (let i = 0; i < n; i++) s += (x[i] - mx) * (y[i] - my); return s / (n - 1); };
  const pearson = (x, y) => { const c = cov(x, y), d = sd(x) * sd(y); return d ? c / d : NaN; };
  const col = (M, j) => M.map(r => r[j]);
  const rowSums = M => M.map(r => sum(r));
  const sorted = a => a.slice().sort((x, y) => x - y);
  const median = a => { const s = sorted(a), n = s.length; return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2; };

  /* ---- special functions ---- */
  function lgamma(z) { const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]; if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z); z -= 1; let x = c[0]; for (let i = 1; i < 9; i++) x += c[i] / (z + i); const t = z + 7.5; return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x); }
  function betacf(a, b, x) { const FP = 1e-300; let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap; if (Math.abs(d) < FP) d = FP; d = 1 / d; let h = d; for (let m = 1; m <= 300; m++) { const m2 = 2 * m; let aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d; if (Math.abs(d) < FP) d = FP; c = 1 + aa / c; if (Math.abs(c) < FP) c = FP; d = 1 / d; h *= d * c; aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d; if (Math.abs(d) < FP) d = FP; c = 1 + aa / c; if (Math.abs(c) < FP) c = FP; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 3e-14) break; } return h; }
  function ibeta(x, a, b) { if (x <= 0) return 0; if (x >= 1) return 1; const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x)); return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b; }
  const tP2 = (t, df) => ibeta(df / (df + t * t), df / 2, 0.5);
  function gammaP(a, x) { if (x <= 0) return 0; if (x < a + 1) { let ap = a, s = 1 / a, d = s; for (let i = 0; i < 500; i++) { ap++; d *= x / ap; s += d; if (Math.abs(d) < Math.abs(s) * 1e-14) break; } return s * Math.exp(-x + a * Math.log(x) - lgamma(a)); } let b = x + 1 - a, c = 1e300, d = 1 / b, h = d; for (let i = 1; i < 500; i++) { const an = -i * (i - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300; c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-14) break; } return 1 - Math.exp(-x + a * Math.log(x) - lgamma(a)) * h; }
  const chiP = (x, df) => 1 - gammaP(df / 2, x / 2);
  function erf(x) { const s = Math.sign(x); x = Math.abs(x); const t = 1 / (1 + 0.3275911 * x); const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return s * y; }
  const ncdf = z => 0.5 * (1 + erf(z / Math.SQRT2));
  function binomTail(n, k) { /* P(X >= k), p = .5 */ let s = 0; for (let i = k; i <= n; i++) s += Math.exp(lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1) - n * Math.LN2); return s; }
  const comb = (n, k) => Math.exp(lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1));

  /* ---- reliability ---- */
  function alphaFromVars(k, sumItemVar, totalVar) { return (k / (k - 1)) * (1 - sumItemVar / totalVar); }
  function itemStats(M) {
    const n = M.length, k = M[0].length, tot = rowSums(M), vt = variance(tot);
    const iv = [], im = [], isd = [], rit = [], aid = [];
    for (let j = 0; j < k; j++) {
      const c = col(M, j); im.push(mean(c)); iv.push(variance(c)); isd.push(sd(c));
      const rest = tot.map((t, i) => t - c[i]);
      rit.push(pearson(c, rest));
    }
    const sIV = sum(iv), alpha = alphaFromVars(k, sIV, vt);
    for (let j = 0; j < k; j++) {
      if (k < 3) { aid.push(NaN); continue; }
      const rest = tot.map((t, i) => t - M[i][j]);
      aid.push(alphaFromVars(k - 1, sIV - iv[j], variance(rest)));
    }
    return { n, k, total: tot, totalMean: mean(tot), totalVar: vt, totalSd: Math.sqrt(vt), itemMean: im, itemVar: iv, itemSd: isd, sumItemVar: sIV, rit, alphaIfDeleted: aid, alpha };
  }
  function kr20(M) {
    const k = M[0].length, tot = rowSums(M), vt = variance(tot);
    const p = []; for (let j = 0; j < k; j++) p.push(mean(col(M, j)));
    const spq = sum(p.map(x => x * (1 - x))), n = M.length, vp = vt * (n - 1) / n; /* p·q are population variances, so the total variance must be too */
    return { k, n, p, spq, totalVar: vp, value: (k / (k - 1)) * (1 - spq / vp), mean: mean(tot) };
  }
  const kr21 = (k, M, v) => (k / (k - 1)) * (1 - (M * (k - M)) / (k * v));
  const spearmanBrown = (r, factor = 2) => factor * r / (1 + (factor - 1) * r);
  const guttman = (v1, v2, vt) => 2 * (1 - (v1 + v2) / vt);
  function splitHalf(M, mode = 'oddeven') {
    const k = M[0].length;
    const A = [], Bi = [];
    for (let j = 0; j < k; j++) { const inA = mode === 'oddeven' ? j % 2 === 0 : j < Math.ceil(k / 2); (inA ? A : Bi).push(j); }
    const h1 = M.map(r => sum(A.map(j => r[j]))), h2 = M.map(r => sum(Bi.map(j => r[j]))), tot = h1.map((v, i) => v + h2[i]);
    const r = pearson(h1, h2), v1 = variance(h1), v2 = variance(h2), vt = variance(tot);
    return { A, B: Bi, r, sb: spearmanBrown(r), guttman: guttman(v1, v2, vt), v1, v2, vt };
  }
  function kappa(T) {
    const k = T.length, N = sum(T.map(r => sum(r)));
    let po = 0; for (let i = 0; i < k; i++) po += T[i][i]; po /= N;
    let pe = 0; for (let i = 0; i < k; i++) { const ri = sum(T[i]), ci = sum(T.map(r => r[i])); pe += (ri / N) * (ci / N); }
    return { N, po, pe, kappa: (po - pe) / (1 - pe) };
  }
  function omega(loadings) {
    const sl = sum(loadings), th = loadings.map(l => 1 - l * l), st = sum(th);
    const ave = mean(loadings.map(l => l * l));
    return { sl, sl2: sl * sl, theta: th, st, omega: sl * sl / (sl * sl + st), ave };
  }
  const sem = (s, r) => s * Math.sqrt(1 - r);

  /* ---- content validity ---- */
  const cvr = (ne, N) => (ne - N / 2) / (N / 2);
  const cvrModified = (big, mid, N) => (big + mid / 2 - N / 2) / (N / 2);
  const LAWSHE = { 5: .99, 6: .99, 7: .99, 8: .75, 9: .78, 10: .62, 11: .59, 12: .56, 13: .54, 14: .51, 15: .49, 20: .42, 25: .37, 30: .33, 35: .31, 40: .29 };
  function cvrCritical(N, alpha = 0.05) { /* Ayre & Scally (2014): exact one-tailed binomial */
    for (let ne = Math.ceil(N / 2); ne <= N; ne++) if (binomTail(N, ne) <= alpha) return { ne, cvr: cvr(ne, N) };
    return { ne: N, cvr: 1 };
  }
  function lawsheCritical(N) { if (LAWSHE[N]) return LAWSHE[N]; const ks = Object.keys(LAWSHE).map(Number); const lo = ks.filter(x => x < N).pop(), hi = ks.find(x => x > N); if (!lo) return .99; if (!hi) return LAWSHE[40]; return LAWSHE[lo] + (LAWSHE[hi] - LAWSHE[lo]) * (N - lo) / (hi - lo); }
  function cvi(R /* judges x items, 1..4 */) {
    const J = R.length, k = R[0].length, icvi = [], kstar = [];
    for (let j = 0; j < k; j++) {
      const A = R.filter(r => r[j] >= 3).length, ic = A / J; icvi.push(ic);
      const pc = comb(J, A) * Math.pow(0.5, J); kstar.push((ic - pc) / (1 - pc));
    }
    return { J, k, icvi, kstar, scviAve: mean(icvi), scviUA: icvi.filter(v => v === 1).length / k };
  }
  function aikenV(ratings, lo, c) { const n = ratings.length, S = sum(ratings.map(r => r - lo)); return { n, S, V: S / (n * (c - 1)) }; }

  /* ---- item analysis ---- */
  function extremeGroups(M, frac = 0.27) {
    const n = M.length, k = M[0].length, tot = rowSums(M);
    const idx = tot.map((t, i) => [t, i]).sort((a, b) => b[0] - a[0]);
    const g = Math.max(2, Math.round(n * frac));
    const hi = idx.slice(0, g).map(x => x[1]), lo = idx.slice(n - g).map(x => x[1]);
    const res = [];
    for (let j = 0; j < k; j++) {
      const a = hi.map(i => M[i][j]), b = lo.map(i => M[i][j]);
      const ma = mean(a), mb = mean(b), va = variance(a), vb = variance(b);
      const se = Math.sqrt(va / a.length + vb / b.length);
      const t = se ? (ma - mb) / se : (ma === mb ? 0 : Infinity);
      const df = se ? Math.pow(va / a.length + vb / b.length, 2) / (Math.pow(va / a.length, 2) / (a.length - 1) + Math.pow(vb / b.length, 2) / (b.length - 1)) : a.length + b.length - 2;
      res.push({ mHi: ma, mLo: mb, t, df, p: isFinite(t) ? tP2(t, df) : 0 });
    }
    return { g, items: res };
  }
  function difficultyDiscrimination(M, frac = 0.27) { /* binary */
    const n = M.length, k = M[0].length, tot = rowSums(M);
    const idx = tot.map((t, i) => [t, i]).sort((a, b) => b[0] - a[0]);
    const g = Math.max(1, Math.round(n * frac));
    const hi = idx.slice(0, g).map(x => x[1]), lo = idx.slice(n - g).map(x => x[1]);
    const out = [];
    for (let j = 0; j < k; j++) {
      const p = mean(col(M, j)), U = sum(hi.map(i => M[i][j])), L = sum(lo.map(i => M[i][j]));
      const rest = tot.map((t, i) => t - M[i][j]);
      out.push({ p, U, L, D: (U - L) / g, rpb: pearson(col(M, j), rest) });
    }
    return { g, items: out };
  }

  /* ---- scores ---- */
  const zScore = (x, m, s) => (x - m) / s;
  const tScore = z => 50 + 10 * z;
  function percentileRank(arr, x) { const n = arr.length; const below = arr.filter(v => v < x).length, eq = arr.filter(v => v === x).length; return 100 * (below + 0.5 * eq) / n; }
  function skewKurt(a) {
    const n = a.length, m = mean(a), s = sd(a);
    const sk = n > 2 && s ? n / ((n - 1) * (n - 2)) * sum(a.map(v => ((v - m) / s) ** 3)) : NaN;
    const ku = n > 3 && s ? n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * sum(a.map(v => ((v - m) / s) ** 4)) - 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)) : NaN;
    const ses = Math.sqrt(6 * n * (n - 1) / ((n - 2) * (n + 1) * (n + 3)));
    const sek = 2 * ses * Math.sqrt((n * n - 1) / ((n - 3) * (n + 5)));
    return { skew: sk, kurt: ku, seSkew: ses, seKurt: sek };
  }
  function levels(min, max, nLevels) { const w = (max - min) / nLevels; const out = []; for (let i = 0; i < nLevels; i++) out.push([min + i * w, min + (i + 1) * w]); return { width: w, ranges: out }; }

  /* ---- DIF ---- */
  function mantelHaenszel(strata /* [{a,b,c,d}] a=ref correct b=ref wrong c=focal correct d=focal wrong */) {
    let num = 0, den = 0, A = 0, EA = 0, VA = 0;
    strata.forEach(s => {
      const N = s.a + s.b + s.c + s.d; if (N < 2) return;
      num += s.a * s.d / N; den += s.b * s.c / N;
      const nR = s.a + s.b, nF = s.c + s.d, m1 = s.a + s.c, m0 = s.b + s.d;
      A += s.a; EA += nR * m1 / N; VA += nR * nF * m1 * m0 / (N * N * (N - 1));
    });
    const alphaMH = num / den, delta = -2.35 * Math.log(alphaMH);
    const chi = VA ? Math.pow(Math.abs(A - EA) - 0.5, 2) / VA : NaN;
    const p = chiP(chi, 1);
    const ad = Math.abs(delta);
    const cls = (ad < 1 || p >= 0.05) ? 'A' : ad >= 1.5 ? 'C' : 'B';
    return { alphaMH, delta, chi, p, cls };
  }

  /* ---- matrix / factor analysis ---- */
  function corrMatrix(M) { const k = M[0].length, C = []; const cols = []; for (let j = 0; j < k; j++) cols.push(col(M, j)); for (let i = 0; i < k; i++) { C.push([]); for (let j = 0; j < k; j++) C[i].push(i === j ? 1 : pearson(cols[i], cols[j])); } return C; }
  function jacobiEigen(A0) {
    const n = A0.length, A = A0.map(r => r.slice()), V = A0.map((r, i) => r.map((_, j) => i === j ? 1 : 0));
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
      if (off < 1e-18) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-15) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
      }
    }
    const vals = A.map((r, i) => r[i]);
    const order = vals.map((v, i) => [v, i]).sort((a, b) => b[0] - a[0]);
    return { values: order.map(o => o[0]), vectors: order.map(o => V.map(r => r[o[1]])) };
  }
  function inverse(A0) {
    const n = A0.length, A = A0.map((r, i) => r.concat(r.map((_, j) => i === j ? 1 : 0)));
    for (let c = 0; c < n; c++) {
      let piv = c; for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[piv][c])) piv = r;
      if (Math.abs(A[piv][c]) < 1e-12) return null;
      [A[c], A[piv]] = [A[piv], A[c]];
      const d = A[c][c]; for (let j = 0; j < 2 * n; j++) A[c][j] /= d;
      for (let r = 0; r < n; r++) if (r !== c) { const f = A[r][c]; if (f) for (let j = 0; j < 2 * n; j++) A[r][j] -= f * A[c][j]; }
    }
    return A.map(r => r.slice(n));
  }
  function logDet(A) { const e = jacobiEigen(A).values; return e.every(v => v > 0) ? sum(e.map(Math.log)) : NaN; }
  function kmo(R) {
    const Ri = inverse(R); if (!Ri) return null; const k = R.length;
    let r2 = 0, p2 = 0; const msa = [];
    for (let i = 0; i < k; i++) {
      let ri = 0, pi = 0;
      for (let j = 0; j < k; j++) if (i !== j) { const pij = -Ri[i][j] / Math.sqrt(Ri[i][i] * Ri[j][j]); ri += R[i][j] ** 2; pi += pij ** 2; }
      r2 += ri; p2 += pi; msa.push(ri / (ri + pi));
    }
    return { kmo: r2 / (r2 + p2), msa };
  }
  function bartlett(R, n) { const p = R.length, ld = logDet(R); const chi = -(n - 1 - (2 * p + 5) / 6) * ld; const df = p * (p - 1) / 2; return { chi, df, p: chiP(chi, df) }; }
  function rng(seed) { let s = seed >>> 0; return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }; }
  function parallelAnalysis(n, p, reps = 40, seed = 7) {
    const r = rng(seed), all = [];
    const gauss = () => { let u = 0, v = 0; while (u === 0) u = r(); v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
    for (let t = 0; t < reps; t++) { const X = []; for (let i = 0; i < n; i++) { const row = []; for (let j = 0; j < p; j++) row.push(gauss()); X.push(row); } all.push(jacobiEigen(corrMatrix(X)).values); }
    const meanEig = [], p95 = [];
    for (let j = 0; j < p; j++) { const v = sorted(all.map(e => e[j])); meanEig.push(mean(v)); p95.push(v[Math.min(v.length - 1, Math.floor(0.95 * v.length))]); }
    return { meanEig, p95 };
  }
  function varimax(L0, maxIt = 100) {
    const p = L0.length, m = L0[0].length; if (m < 2) return L0.map(r => r.slice());
    const h = L0.map(r => Math.sqrt(sum(r.map(x => x * x))) || 1);
    let L = L0.map((r, i) => r.map(x => x / h[i]));
    for (let it = 0; it < maxIt; it++) {
      let changed = 0;
      for (let a = 0; a < m - 1; a++) for (let b = a + 1; b < m; b++) {
        let A = 0, B = 0, C = 0, D = 0;
        for (let i = 0; i < p; i++) { const x = L[i][a], y = L[i][b], u = x * x - y * y, v = 2 * x * y; A += u; B += v; C += u * u - v * v; D += 2 * u * v; }
        const num = D - 2 * A * B / p, den = C - (A * A - B * B) / p;
        const phi = Math.atan2(num, den) / 4;
        if (Math.abs(phi) > 1e-7) { changed++; const c = Math.cos(phi), s = Math.sin(phi); for (let i = 0; i < p; i++) { const x = L[i][a], y = L[i][b]; L[i][a] = c * x + s * y; L[i][b] = -s * x + c * y; } }
      }
      if (!changed) break;
    }
    return L.map((r, i) => r.map(x => x * h[i]));
  }
  function efa(M, nf) {
    const R = corrMatrix(M), e = jacobiEigen(R), p = R.length;
    const m = nf || Math.max(1, e.values.filter(v => v > 1).length);
    const L = []; for (let i = 0; i < p; i++) { const row = []; for (let f = 0; f < m; f++) row.push(e.vectors[f][i] * Math.sqrt(Math.max(0, e.values[f]))); L.push(row); }
    for (let f = 0; f < m; f++) { const s = sum(L.map(r => r[f])); if (s < 0) L.forEach(r => r[f] = -r[f]); }
    const rot = varimax(L);
    const varExp = e.values.map(v => v / p * 100);
    return { R, eigen: e.values, varExp, nf: m, loadings: L, rotated: rot };
  }

  /* ---- descriptive helpers ---- */
  function describe(a) { const s = sorted(a); return { n: a.length, mean: mean(a), sd: sd(a), min: s[0], max: s[s.length - 1], median: median(a), ...skewKurt(a) }; }

  return { sum, mean, variance, sd, cov, pearson, median, sorted, col, rowSums, tP2, chiP, ncdf, ibeta, lgamma,
    alphaFromVars, itemStats, kr20, kr21, spearmanBrown, guttman, splitHalf, kappa, omega, sem,
    cvr, cvrModified, cvrCritical, lawsheCritical, cvi, aikenV, extremeGroups, difficultyDiscrimination,
    zScore, tScore, percentileRank, skewKurt, levels, mantelHaenszel, corrMatrix, jacobiEigen, inverse, kmo, bartlett, parallelAnalysis, varimax, efa, describe };
})();
if (typeof module !== 'undefined') module.exports = PSY;

/* ===== Core: environment, state, shell, router, helpers, gamification ===== */
"use strict";
const C = JSON.parse(document.getElementById('content').textContent);
const $ = (s, r = document) => r.querySelector(s), $$ = (s, r = document) => [...r.querySelectorAll(s)];
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const fmt = (x, d = 2) => { if (x == null || !isFinite(x)) return '—'; let s = (Math.abs(x) < 1e-12 ? 0 : x).toFixed(d); return s === '-0.00' ? '0.00' : s; };
const fmtp = p => p < 0.001 ? '< .001' : fmt(p, 3).replace(/^0/, '');
const shuffle = a => { a = a.slice(); for (let i = a.length - 1; i > 0; i--) { const j = Math.random() * (i + 1) | 0; [a[i], a[j]] = [a[j], a[i]]; } return a; };
const pick = (a, n) => shuffle(a).slice(0, n);
const uid = (p = '') => p + Math.random().toString(36).slice(2, 8) + Date.now().toString(36).slice(-4);
const today = () => { const d = new Date(); return d.getFullYear() + '-' + String(d.getMonth() + 1).padStart(2, '0') + '-' + String(d.getDate()).padStart(2, '0'); };
const dayDiff = (a, b) => Math.round((new Date(b + 'T12:00') - new Date(a + 'T12:00')) / 864e5);
const toNum = s => Number(String(s).replace(/[٠-٩]/g, d => '٠١٢٣٤٥٦٧٨٩'.indexOf(d)).replace(/٫/g, '.').replace(/,/g, '.'));
function parseNums(s) { return String(s).replace(/[٠-٩]/g, d => '٠١٢٣٤٥٦٧٨٩'.indexOf(d)).replace(/[،؛]/g, ',').replace(/٫/g, '.').split(/[\s,;]+/).filter(Boolean).map(Number).filter(v => isFinite(v)); }
function parseMatrix(s) {
  const lines = String(s).replace(/[٠-٩]/g, d => '٠١٢٣٤٥٦٧٨٩'.indexOf(d)).replace(/٫/g, '.').split(/\r?\n/).map(l => l.trim()).filter(Boolean);
  let head = null; const rows = [];
  lines.forEach((l, i) => {
    const cells = l.split(/[\t,;،]+|\s+/).filter(x => x !== '');
    const nums = cells.map(Number);
    if (nums.some(v => !isFinite(v))) { if (i === 0) head = cells; return; }
    rows.push(nums);
  });
  if (!rows.length) return { error: 'لم أجد صفوفًا رقمية. ضع كل مستجيب في سطر، وافصل بين الدرجات بمسافة أو فاصلة.' };
  const k = rows[0].length;
  if (rows.some(r => r.length !== k)) return { error: 'عدد القيم غير متساوٍ في كل الصفوف. تأكد أن كل سطر فيه درجات جميع العبارات.' };
  return { rows, k, n: rows.length, head };
}

const ICON = {
  heart: '<path d="M12 20s-7-4.4-7-10a4 4 0 0 1 7-2.6A4 4 0 0 1 19 10c0 5.6-7 10-7 10z"/>',
  home: '<path d="M3 11 12 4l9 7v9a1 1 0 0 1-1 1h-5v-6H9v6H4a1 1 0 0 1-1-1z"/>',
  book: '<path d="M4 4h6a3 3 0 0 1 3 3v13a2 2 0 0 0-2-2H4zM20 4h-6a3 3 0 0 0-3 3v13a2 2 0 0 1 2-2h7z"/>',
  flask: '<path d="M9 3h6M10 3v6l-5 9a2 2 0 0 0 1.8 3h10.4a2 2 0 0 0 1.8-3l-5-9V3"/><path d="M7.5 15h9"/>',
  build: '<path d="M4 20h16M6 16V9l6-5 6 5v7"/><path d="M10 16v-4h4v4"/>',
  sigma: '<path d="M18 4H6l6 8-6 8h12"/>',
  cards: '<rect x="3" y="6" width="14" height="14" rx="2"/><path d="M7 3h12a2 2 0 0 1 2 2v12"/>',
  trophy: '<path d="M8 21h8M12 17v4M7 4h10v5a5 5 0 0 1-10 0z"/><path d="M17 5h3v2a3 3 0 0 1-3 3M7 5H4v2a3 3 0 0 0 3 3"/>',
  user: '<circle cx="12" cy="8" r="4"/><path d="M4 21a8 8 0 0 1 16 0"/>',
  search: '<circle cx="11" cy="11" r="7"/><path d="m20 20-4-4"/>',
  moon: '<path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/>',
  more: '<circle cx="5" cy="12" r="1.6"/><circle cx="12" cy="12" r="1.6"/><circle cx="19" cy="12" r="1.6"/>',
  check: '<path d="m5 12 5 5 9-10"/>', x: '<path d="M6 6l12 12M18 6 6 18"/>', back: '<path d="m9 6 6 6-6 6"/>',
  fire: '<path d="M12 2c1 4 5 5.5 5 10.5A5 5 0 0 1 7 12.5c0-2 1-3.5 2-4.5 0 2 1 3 2 3 0-3-1-5 1-9z" fill="currentColor" stroke="none"/>',
  star: '<path d="m12 2 3 6.5 7 .8-5.2 4.8 1.5 7L12 17.6 5.7 21l1.5-7L2 9.3l7-.8z" fill="currentColor" stroke="none"/>',
  bolt: '<path d="M13 2 4 14h7l-1 8 9-12h-7z" fill="currentColor" stroke="none"/>',
  copy: '<rect x="8" y="8" width="12" height="12" rx="2"/><path d="M16 8V5a1 1 0 0 0-1-1H5a1 1 0 0 0-1 1v10a1 1 0 0 0 1 1h3"/>',
  download: '<path d="M12 3v12m0 0-5-5m5 5 5-5M4 21h16"/>', plus: '<path d="M12 5v14M5 12h14"/>', trash: '<path d="M4 7h16M10 11v6M14 11v6M6 7l1 13h10l1-13M9 7V4h6v3"/>',
  chart: '<path d="M4 20V10M10 20V4M16 20v-7M22 20H2"/>', form: '<rect x="4" y="3" width="16" height="18" rx="2"/><path d="M8 8h8M8 12h8M8 16h5"/>',
  link: '<path d="M10 14a4 4 0 0 0 6 0l3-3a4 4 0 0 0-6-6l-1 1M14 10a4 4 0 0 0-6 0l-3 3a4 4 0 0 0 6 6l1-1"/>',
  lock: '<rect x="5" y="11" width="14" height="10" rx="2"/><path d="M8 11V7a4 4 0 0 1 8 0v4"/>', play: '<path d="M7 4v16l13-8z"/>',
  timer: '<circle cx="12" cy="13" r="8"/><path d="M12 9v4l2 2M9 2h6"/>', redo: '<path d="M3 12a9 9 0 1 0 3-6.7L3 8"/><path d="M3 3v5h5"/>',
  match: '<rect x="3" y="4" width="7" height="6" rx="1.5"/><rect x="14" y="14" width="7" height="6" rx="1.5"/><path d="M10 7h4a3 3 0 0 1 3 3v4"/>',
  calc: '<rect x="5" y="3" width="14" height="18" rx="2"/><path d="M8 7h8M8 11h2M12 11h2M16 11v6M8 15h2M12 15h2"/>',
  shield: '<path d="M12 3 4 6v6c0 5 3.5 8 8 9 4.5-1 8-4 8-9V6z"/>', users: '<circle cx="9" cy="8" r="3.5"/><path d="M2 20a7 7 0 0 1 14 0M16 4.5a3.5 3.5 0 0 1 0 7M18 13.5a6 6 0 0 1 4 6.5"/>',
  upload: '<path d="M12 21V9m0 0-5 5m5-5 5 5M4 3h16"/>', eye: '<path d="M2 12s4-7 10-7 10 7 10 7-4 7-10 7S2 12 2 12z"/><circle cx="12" cy="12" r="3"/>',
  whats: '<path d="M4 20l1.5-4A8 8 0 1 1 8 18.5z"/>', gear: '<circle cx="12" cy="12" r="3"/><path d="M19.4 15a1.7 1.7 0 0 0 .3 1.8l.1.1a2 2 0 1 1-2.8 2.8l-.1-.1a1.7 1.7 0 0 0-1.8-.3 1.7 1.7 0 0 0-1 1.5V21a2 2 0 1 1-4 0v-.1a1.7 1.7 0 0 0-1.1-1.5 1.7 1.7 0 0 0-1.8.3l-.1.1a2 2 0 1 1-2.8-2.8l.1-.1a1.7 1.7 0 0 0 .3-1.8 1.7 1.7 0 0 0-1.5-1H3a2 2 0 1 1 0-4h.1a1.7 1.7 0 0 0 1.5-1.1 1.7 1.7 0 0 0-.3-1.8l-.1-.1a2 2 0 1 1 2.8-2.8l.1.1a1.7 1.7 0 0 0 1.8.3H9a1.7 1.7 0 0 0 1-1.5V3a2 2 0 1 1 4 0v.1a1.7 1.7 0 0 0 1 1.5 1.7 1.7 0 0 0 1.8-.3l.1-.1a2 2 0 1 1 2.8 2.8l-.1.1a1.7 1.7 0 0 0-.3 1.8V9a1.7 1.7 0 0 0 1.5 1H21a2 2 0 1 1 0 4h-.1a1.7 1.7 0 0 0-1.5 1z"/>'
};
const ic = (n, cls = '') => `<svg class="${cls}" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${ICON[n] || ''}</svg>`;

/* ---------- storage ---------- */
const LS = { get(k, d) { try { const v = localStorage.getItem(k); return v ? JSON.parse(v) : d; } catch (e) { return d; } }, set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { } } };
const KEY = 'scaleslab-v2';
const RANKS = [[0, 'مساعد باحث'], [150, 'باحث مبتدئ'], [400, 'باحث'], [800, 'باحث أول'], [1400, 'خبير قياس'], [2200, 'محكّم مقاييس'], [3200, 'عالِم قياس نفسي']];
function fresh() { return { v: 2, xp: 0, streak: 0, best: 0, last: null, goal: 40, hist: {}, read: {}, quiz: {}, cards: {}, miss: [], badges: [], quests: null, correct: 0, answered: 0, perfect: 0, labs: {}, probs: 0, matches: 0, fs: 1, theme: '', name: '', updated: 0 }; }
let S = Object.assign(fresh(), LS.get(KEY, {}));
const Cloud = { push: null }; /* set by account layer */
let saveT = null;
function save() { S.updated = Date.now(); LS.set(KEY, S); clearTimeout(saveT); if (Cloud.push) saveT = setTimeout(() => Cloud.push(S), 1500); }

/* ---------- gamification ---------- */
const rankOf = xp => { let r = RANKS[0], nx = null; RANKS.forEach((x, i) => { if (xp >= x[0]) { r = x; nx = RANKS[i + 1] || null; } }); return { name: r[1], min: r[0], next: nx }; };
const QPOOL = [
  { k: 'xp', n: 30, t: 'اجمع 30 نقطة' }, { k: 'read', n: 2, t: 'أنهِ درسين' }, { k: 'quiz', n: 1, t: 'أنهِ اختبار جزء' },
  { k: 'lab', n: 1, t: 'استخدم معملًا حسابيًا' }, { k: 'cards', n: 10, t: 'راجع 10 بطاقات' }, { k: 'items', n: 5, t: 'أضف 5 عبارات لمقياسك' }, { k: 'prob', n: 3, t: 'حل 3 مسائل حسابية' }];
function tick() {
  const t = today();
  if (S.last && S.last !== t && dayDiff(S.last, t) > 1) S.streak = 0;
  if (!S.quests || S.quests.day !== t) {
    let x = (t.split('-').join('') | 0) % 9973; const chosen = [];
    while (chosen.length < 3) { x = (x * 9301 + 49297) % 233280; const q = QPOOL[x % QPOOL.length]; if (!chosen.some(c => c.k === q.k)) chosen.push(Object.assign({ p: 0, done: false }, q)); }
    S.quests = { day: t, list: chosen };
  }
}
function quest(k, amt = 1) { if (!S.quests) return; S.quests.list.forEach(q => { if (q.k !== k || q.done) return; q.p += amt; if (q.p >= q.n) { q.done = true; S.xp += 10; toast('أنجزت مهمة اليوم: ' + q.t + ' (+10)'); } }); }
function addXP(n) {
  if (typeof MOTION !== 'undefined') MOTION.xp(n);
  const t = today(); S.xp += n; S.hist[t] = (S.hist[t] || 0) + n; quest('xp', n);
  if (S.last !== t) { S.streak = (S.last && dayDiff(S.last, t) === 1) ? S.streak + 1 : 1; S.last = t; S.best = Math.max(S.best, S.streak); if (S.streak > 1) toast('سلسلتك الآن ' + S.streak + ' أيام متتالية'); }
  checkBadges(); save(); renderStats();
}
const BADGES = [
  { id: 'first', t: 'أول درس', d: 'أنهِ أول درس', f: () => Object.keys(S.read).length >= 1 },
  { id: 'quiz1', t: 'أول اختبار', d: 'أنهِ اختبار جزء', f: () => Object.keys(S.quiz).length >= 1 },
  { id: 'lab1', t: 'أول معمل', d: 'احسب أول نتيجة في المعامل', f: () => Object.keys(S.labs).length >= 1 },
  { id: 'lab6', t: 'متمكن من المعادلات', d: 'استخدم 6 معامل مختلفة', f: () => Object.keys(S.labs).length >= 6 },
  { id: 'scale1', t: 'باني المقاييس', d: 'أنشئ أول مقياس', f: () => (typeof SCALES !== 'undefined') && SCALES.list().length >= 1 },
  { id: 'resp30', t: 'جامع البيانات', d: 'اجمع 30 استجابة لمقياسك', f: () => (typeof SCALES !== 'undefined') && SCALES.list().some(s => (s.responses || []).length >= 30) },
  { id: 's3', t: 'ثلاثة أيام', d: 'سلسلة 3 أيام', f: () => S.best >= 3 },
  { id: 's7', t: 'أسبوع بحثي', d: 'سلسلة 7 أيام', f: () => S.best >= 7 },
  { id: 'perfect', t: 'ثبات تام', d: 'اختبار كامل بلا أخطاء', f: () => S.perfect >= 1 },
  { id: 'c100', t: 'مئة إجابة صحيحة', d: '100 إجابة صحيحة', f: () => S.correct >= 100 },
  { id: 'prob10', t: 'حاسب ماهر', d: 'حل 10 مسائل حسابية', f: () => S.probs >= 10 },
  { id: 'half', t: 'نصف الطريق', d: 'أنهِ نصف الأجزاء', f: () => LEARN.doneCount() >= LEARN.total() / 2 },
  { id: 'all', t: 'عالِم قياس', d: 'أنهِ كل الأجزاء', f: () => LEARN.doneCount() >= LEARN.total() },
  { id: 'acct', t: 'عضو مسجل', d: 'أنشئ حسابًا لحفظ تقدمك', f: () => !!(window.ACCT && ACCT.signedIn()) }
];
function checkBadges() { BADGES.forEach(b => { if (!S.badges.includes(b.id) && b.f()) { S.badges.push(b.id); S.xp += 15; if (typeof MOTION !== 'undefined') setTimeout(() => MOTION.stamp('وسام: ' + b.t), 250); setTimeout(() => toast('وسام جديد: ' + b.t + ' (+15)'), 300); } }); }
function labUse(k) { const first = !S.labs[k]; S.labs[k] = (S.labs[k] || 0) + 1; quest('lab', 1); if (first) addXP(5); else save(); }

/* ---------- UI helpers ---------- */
let toastT; function toast(m) { const t = $('#toast'); t.textContent = m; t.hidden = false; clearTimeout(toastT); toastT = setTimeout(() => t.hidden = true, 2600); }
function modal(html, onMount) {
  const bg = document.createElement('div'); bg.className = 'modal-bg'; bg.setAttribute('data-lenis-prevent', ''); bg.innerHTML = `<div class="modal" role="dialog" aria-modal="true">${html}</div>`;
  document.body.appendChild(bg);
  const close = () => { bg.remove(); document.removeEventListener('keydown', onKey); };
  const onKey = e => { if (e.key === 'Escape') close(); };
  document.addEventListener('keydown', onKey);
  bg.addEventListener('click', e => { if (e.target === bg) close(); });
  $$('[data-close]', bg).forEach(b => b.onclick = close);
  onMount && onMount(bg, close);
  if (typeof MOTION !== 'undefined' && MOTION.on) gsap.from($('.modal', bg), { y: 24, autoAlpha: 0, scale: 0.97, duration: 0.32, ease: 'power3.out' });
  const f = $('input,textarea,button', bg); f && f.focus();
  return close;
}
function confirmBox(msg, yes, onYes) { modal(`<h3>تأكيد</h3><p>${esc(msg)}</p><div class="row"><button class="btn danger" id="cfY">${esc(yes)}</button><button class="btn ghost" data-close>إلغاء</button></div>`, (bg, close) => { $('#cfY', bg).onclick = () => { close(); onYes(); }; }); }
async function copyText(text) {
  try { await navigator.clipboard.writeText(text); toast('نُسخ'); return true; } catch (e) { }
  modal(`<h3>انسخ النص</h3><p class="muted small">حدّد النص ثم انسخه.</p><textarea id="cpTa" class="data" readonly></textarea><div class="row" style="margin-top:8px"><button class="btn" data-close>تم</button></div>`, bg => { const ta = $('#cpTa', bg); ta.value = text; ta.focus(); ta.select(); });
  return false;
}
async function downloadFile(name, text, mime = 'text/plain;charset=utf-8') {
  if (ENV.claude) {
    try { const dl = await window.claude.use('downloads'); if (dl) { await dl.save({ filename: name, data: text }); return; } } catch (e) { }
    return copyText(text);
  }
  try { const blob = new Blob([(/\.(csv|txt)$/i.test(name) ? '﻿' : '') + text], { type: mime }); const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); }
  catch (e) { copyText(text); }
}
function ring(pct, label, sub) { const R = 42, Cc = 2 * Math.PI * R; return `<div class="ring" role="img" aria-label="${esc(sub)}"><svg viewBox="0 0 96 96"><circle cx="48" cy="48" r="${R}" fill="none" stroke="var(--surface-2)" stroke-width="9"/><circle cx="48" cy="48" r="${R}" fill="none" stroke="var(--amber)" stroke-width="9" stroke-linecap="round" stroke-dasharray="${Cc}" stroke-dashoffset="${Cc * (1 - Math.min(1, pct))}"/></svg><div class="c"><div><b>${label}</b><small>${esc(sub)}</small></div></div></div>`; }
const frac = (a, b) => `<span class="fr"><span>${a}</span><span>${b}</span></span>`;
function tableHTML(head, rows, cls = '') { return `<div class="tw"><table class="${cls}"><thead><tr>${head.map(h => `<th>${h}</th>`).join('')}</tr></thead><tbody>${rows.map(r => `<tr${r.cls ? ` class="${r.cls}"` : ''}>${(r.cells || r).map(c => `<td${typeof c === 'object' && c && c.n ? ' class="n"' : ''}>${typeof c === 'object' && c ? c.v : c}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`; }
const N = (v, d = 2) => ({ n: 1, v: typeof v === 'number' ? fmt(v, d) : v });
function relVerdict(r) { if (!isFinite(r)) return ['bad', 'غير قابل للحساب']; if (r >= 0.9) return ['good', 'ممتاز (≥ 0.90)']; if (r >= 0.8) return ['good', 'جيد (0.80–0.89)']; if (r >= 0.7) return ['warn', 'مقبول للبحث (0.70–0.79)']; return ['bad', 'ضعيف (أقل من 0.70)']; }
function barSVG(vals, labels, opts = {}) {
  const W = 560, H = opts.h || 180, P = 26, n = vals.length, mx = Math.max(opts.max || 0, ...vals.map(v => Math.abs(v))) || 1, bw = (W - 2 * P) / n;
  let g = `<line x1="${P}" y1="${H - P}" x2="${W - P}" y2="${H - P}" stroke="var(--line)"/>`;
  if (opts.ref != null) { const y = H - P - (H - 2 * P) * opts.ref / mx; g += `<line x1="${P}" x2="${W - P}" y1="${y}" y2="${y}" stroke="var(--amber)" stroke-dasharray="4 4"/><text x="${W - P}" y="${y - 4}" text-anchor="end">${esc(opts.refLabel || '')}</text>`; }
  vals.forEach((v, i) => { const h = (H - 2 * P) * Math.max(0, v) / mx, x = P + i * bw; g += `<g><title>${esc(labels[i])}: ${fmt(v, opts.d ?? 2)}</title><rect x="${x + 2}" y="${H - P - h}" width="${Math.max(2, bw - 4)}" height="${h}" rx="3" fill="${opts.color ? opts.color(v, i) : 'var(--primary)'}"/></g>`; if (n <= 30) g += `<text x="${x + bw / 2}" y="${H - P + 14}" text-anchor="middle">${esc(labels[i])}</text>`; });
  if (opts.line) { const pts = opts.line.map((v, i) => `${P + i * bw + bw / 2},${H - P - (H - 2 * P) * Math.max(0, v) / mx}`).join(' '); g += `<polyline points="${pts}" fill="none" stroke="var(--amber)" stroke-width="2" stroke-dasharray="5 4"/>`; }
  return `<div class="chart"><svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(opts.label || 'رسم بياني')}">${g}</svg></div>`;
}
function histSVG(arr) {
  if (!arr.length) return '';
  const s = PSY.sorted(arr), mn = s[0], mx = s[s.length - 1], k = Math.max(5, Math.min(12, Math.ceil(1 + 3.322 * Math.log10(arr.length))));
  const w = (mx - mn) / k || 1, bins = Array(k).fill(0); arr.forEach(v => bins[Math.min(k - 1, Math.floor((v - mn) / w))]++);
  return barSVG(bins, bins.map((_, i) => fmt(mn + i * w, 0)), { d: 0, label: 'توزيع الدرجات الكلية' });
}

/* ---------- environment ---------- */
const ENV = { claude: !!(window.claude && typeof window.claude.use === 'function') };
ENV.web = !ENV.claude;

/* ---------- shell & router ---------- */
const NAV = [
  ['home', 'الرئيسية', 'home'], ['know', 'اعرف نفسك', 'heart'], ['learn', 'الدروس', 'book'], ['labs', 'المعامل', 'flask'], ['build', 'ابنِ مقياسك', 'build'], ['exam', 'أسئلة موضوعية', 'check'], ['bank', 'بنك الأسئلة', 'cards'],
  ['formulas', 'المعادلات', 'sigma'], ['practice', 'التدريب', 'cards'], ['neuro', 'المختبر العصبي', 'flask'], ['respond', 'أجب عن مقياس', 'form'], ['schedules', 'جداول التعزيز', 'trophy'], ['me', 'إنجازاتي', 'trophy'], ['account', 'حسابي', 'user'], ['admin', 'لوحة المسؤول', 'shield']
];
const LOGO_SVG = `<svg viewBox="0 0 40 40" aria-hidden="true"><defs><linearGradient id="lgS" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#7c5cff"/><stop offset="1" stop-color="#1f1452"/></linearGradient></defs><rect width="40" height="40" rx="11" fill="url(#lgS)"/><rect x="8.5" y="7" width="23" height="27" rx="3.4" fill="#fff"/><rect x="14.5" y="5" width="11" height="4.2" rx="1.6" fill="#bef264"/><circle cx="27" cy="13.5" r="2.3" fill="#65a30d" stroke="#65a30d" stroke-width="1.3"/><path d="M25.9 13.5 l.8 .8 1.5-1.6" stroke="#fff" stroke-width="1" fill="none" stroke-linecap="round" stroke-linejoin="round"/><rect x="12" y="12.4" width="11" height="2.2" rx="1.1" fill="#b9b0e8"/><circle cx="27" cy="20.1" r="2.3" fill="#65a30d" stroke="#65a30d" stroke-width="1.3"/><path d="M25.9 20.1 l.8 .8 1.5-1.6" stroke="#fff" stroke-width="1" fill="none" stroke-linecap="round" stroke-linejoin="round"/><rect x="15" y="19.0" width="8" height="2.2" rx="1.1" fill="#b9b0e8"/><circle cx="27" cy="26.7" r="2.3" fill="#fff" stroke="#7c5cff" stroke-width="1.3"/><rect x="13" y="25.599999999999998" width="10" height="2.2" rx="1.1" fill="#b9b0e8"/></svg>`;
const GUEST_SVG = `<svg viewBox="0 0 24 24" aria-hidden="true" fill="none"><circle cx="12" cy="9" r="3.6" fill="currentColor"/><path d="M5.2 19.2c1.3-3.1 3.9-4.7 6.8-4.7s5.5 1.6 6.8 4.7" stroke="currentColor" stroke-width="2.2" stroke-linecap="round"/></svg>`;
const avatarInner = nm => nm && nm.trim() ? esc(nm.trim()[0]) : GUEST_SVG;
const VIEWS = {};
let current = 'home';
function shell() {
  $('#app').innerHTML = `
  <aside class="rail" aria-label="التنقل">
    <a class="brand" href="#home"><span class="logo">${LOGO_SVG}</span><span>مختبر المقاييس النفسية<small>تعلّم · احسب · ابنِ مقياسك</small></span></a>
    <nav>${NAV.slice(0, 9).map(n => `<button data-go="${n[0]}">${ic(n[2])}${n[1]}</button>`).join('')}<div class="navsep">أدوات أخرى</div>${NAV.slice(9).map(n => `<button data-go="${n[0]}">${ic(n[2])}${n[1]}</button>`).join('')}</nav>
    <div class="me" id="railMe"></div>
  </aside>
  <div class="col">
    <header class="topbar">
      <a class="brand" href="#home"><span class="logo">${LOGO_SVG}</span><span>مختبر المقاييس النفسية</span></a>
      <span class="pill fire" id="pStreak" title="أيام متتالية"></span><span class="pill xp" id="pXp" title="نقاطك"></span>
      <button class="ibtn" id="bSearch" aria-label="بحث">${ic('search')}</button>
      <button class="ibtn" id="bTheme" aria-label="تبديل الوضع الليلي">${ic('moon')}</button>
      <button class="avatar" id="bAcct" aria-label="حسابي"></button>
    </header>
    <main class="main" id="view" tabindex="-1"></main>
  </div>
  <nav class="tabbar" aria-label="التنقل السفلي">${[['home', 'الرئيسية', 'home'], ['learn', 'الدروس', 'book'], ['labs', 'المعامل', 'flask'], ['build', 'ابنِ مقياسك', 'build'], ['more', 'المزيد', 'more']].map(n => `<button data-go="${n[0]}">${ic(n[2])}${n[1]}</button>`).join('')}</nav>`;
  $$('[data-go]').forEach(b => b.onclick = () => b.dataset.go === 'more' ? moreSheet() : go(b.dataset.go));
  $('#bSearch').onclick = searchSheet;
  $('#bTheme').onclick = () => { const cur = document.documentElement.dataset.theme || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'); S.theme = cur === 'dark' ? 'light' : 'dark'; applyPrefs(); save(); };
  $('#bAcct').onclick = () => go('account');
}
function moreSheet() {
  modal(`<h3>المزيد</h3><div class="stack">${[['bank', 'بنك الأسئلة', 'cards', 'كل أسئلة المادة: بحث وتصفية واختبار سريع'], ['schedules', 'جداول التعزيز وكروت العلماء', 'trophy', 'جرّب FR وVR وFI وVI وشوف السجل التراكمي، وكروتك'], ['know', 'اعرف نفسك', 'heart', 'رحلة بشكل لعبة: المزاج والذاكرة والتركيز والحالة النفسية وحياتك الرقمية'], ['neuro', 'المختبر العصبي (TRIBE v2)', 'flask', 'تماسك عبارات مقياسك وفق التنبؤ العصبي — استكشافي'], ['exam', 'أسئلة موضوعية', 'check', 'صح وخطأ، واختيار من متعدد، وامتحان تجريبي ونموذج محلول'], ['formulas', 'المعادلات', 'sigma', 'كل معادلات الصدق والثبات بأمثلة'], ['practice', 'التدريب', 'cards', 'بطاقات ومطابقة ومسائل حسابية'], ['respond', 'أجب عن مقياس', 'form', 'استخدم كود المقياس الذي أرسله الباحث'], ['me', 'إنجازاتي', 'trophy', 'النقاط والأوسمة والإعدادات'], ['account', 'حسابي', 'user', 'تسجيل الدخول وحفظ التقدم'], ['admin', 'لوحة المسؤول', 'shield', 'للمسؤول فقط: فتح وقفل المراحل ومتابعة النتائج']].map(n => `<button class="linkcard" data-m="${n[0]}"><span class="ic">${ic(n[2])}</span><span><b>${n[1]}</b><span>${n[3]}</span></span></button>`).join('')}</div><div class="row" style="margin-top:12px"><button class="btn ghost" data-close>إغلاق</button></div>`,
    (bg, close) => $$('[data-m]', bg).forEach(b => b.onclick = () => { close(); go(b.dataset.m); }));
}
function go(v) { if (v === 'labs' && typeof LABS !== 'undefined') LABS.clear(); if (v === 'learn' && typeof LEARN !== 'undefined') LEARN.clearOpen(); if (v === 'exam' && typeof EXAM !== 'undefined') EXAM.clear(); if (location.hash !== '#' + v) location.hash = v; else route(); }
function route() {
  tick();
  let h = (location.hash || '#home').slice(1);
  let arg = null;
  const m = h.match(/^form-([A-Za-z0-9]+)$/); if (m) { h = 'respond'; arg = m[1]; } const tk = h.match(/^take\/([A-Za-z0-9_-]+)$/); if (tk) { h = 'take'; arg = tk[1]; }
  if (!VIEWS[h]) h = 'home';
  current = h;
  $$('[data-go]').forEach(b => b.setAttribute('aria-current', b.dataset.go === h || (b.dataset.go === 'more' && ['qbuild', 'take', 'bank', 'know', 'schedules', 'admin', 'neuro', 'exam', 'formulas', 'practice', 'respond', 'me', 'account'].includes(h)) ? 'page' : 'false'));
  const v = $('#view'); v.innerHTML = '<div class="container" id="vc"></div>';
  if (typeof MOTION !== 'undefined') MOTION.scrollTop();
  VIEWS[h]($('#vc'), arg);
  if (typeof MOTION !== 'undefined') MOTION.page($('#vc'));
  renderStats();
  if (!routeKeepScroll) scrollTo(0, 0); routeKeepScroll = false;
}
let routeKeepScroll = false;
function rerender() { routeKeepScroll = true; route(); }
function renderStats() {
  const done = (S.hist[today()] || 0) > 0;
  $('#pStreak').innerHTML = ic('fire') + `<span class="num">${S.streak}</span>`; $('#pStreak').style.opacity = done ? 1 : .55;
  $('#pXp').innerHTML = ic('bolt') + `<span class="num">${S.xp}</span>`;
  const nm = (window.ACCT && ACCT.displayName()) || S.name || '';
  $('#bAcct').innerHTML = avatarInner(nm); $('#bAcct').classList.toggle('guest', !(nm && nm.trim()));
  const rk = rankOf(S.xp);
  const rm = $('#railMe'); if (rm) rm.innerHTML = `<div class="row" style="gap:8px"><span class="avatar${nm && nm.trim() ? '' : ' guest'}">${avatarInner(nm)}</span><div class="grow small"><b>${esc(nm || 'زائر')}</b><div class="muted">${esc(rk.name)}</div></div></div><div class="row small"><span class="pill xp">${ic('bolt')}<span class="num">${S.xp}</span></span><span class="pill fire">${ic('fire')}<span class="num">${S.streak}</span></span></div>${window.ACCT && !ACCT.signedIn() && ACCT.available() ? `<button class="btn sm soft" data-go2="account">احفظ تقدمك بحساب</button>` : ''}<div class="row"><button class="btn sm ghost" id="rTheme">${ic('moon')} الوضع</button><button class="btn sm ghost" id="rFs">أ± الخط</button></div>`;
  $$('[data-go2]').forEach(b => b.onclick = () => go(b.dataset.go2));
  const rt = $('#rTheme'); if (rt) rt.onclick = () => $('#bTheme').click();
  const rf = $('#rFs'); if (rf) rf.onclick = () => { S.fs = S.fs >= 1.2 ? 1 : +(S.fs + 0.1).toFixed(1); applyPrefs(); save(); };
}
function applyPrefs() { if (S.theme) document.documentElement.dataset.theme = S.theme; else delete document.documentElement.dataset.theme; document.documentElement.style.setProperty('--fs', S.fs || 1); }

/* ---------- search ---------- */
function searchSheet() {
  modal(`<h3>بحث في الموقع</h3><input id="sq" type="search" placeholder="مثال: ألفا، صدق المحتوى، راش" aria-label="كلمة البحث"><div id="sr" class="stack" style="margin-top:10px;max-height:60vh;overflow:auto"></div>`, (bg, close) => {
    const q = $('#sq', bg), out = $('#sr', bg);
    q.oninput = () => {
      const k = q.value.trim(); if (k.length < 2) { out.innerHTML = '<p class="muted small">اكتب حرفين على الأقل.</p>'; return; }
      const hits = [];
      LEARN.parts().forEach(p => p.sections.forEach((s, si) => { const t = LEARN.text(s); const i = t.indexOf(k); if (i >= 0 || s.title.includes(k)) hits.push({ p, si, s, snip: i >= 0 ? t.slice(Math.max(0, i - 50), i + k.length + 70) : '' }); }));
      const g = C.glossary.filter(x => x.ar.includes(k) || x.en.toLowerCase().includes(k.toLowerCase())).slice(0, 5);
      const labs = LABS.list.filter(l => l.t.includes(k) || l.d.includes(k)).slice(0, 4);
      out.innerHTML = labs.map(l => `<button class="linkcard" data-lab="${l.id}"><span class="ic">${ic('flask')}</span><span><b>${esc(l.t)}</b><span>معمل حسابي</span></span></button>`).join('') +
        g.map(x => `<div class="card"><b>${esc(x.ar)}</b> <span class="ltr muted small">${esc(x.en)}</span><p class="small" style="margin:4px 0 0">${esc(x.d)}</p></div>`).join('') +
        hits.slice(0, 30).map((h, i) => `<button class="linkcard" data-hit="${i}"><span class="ic">${ic('book')}</span><span><b>${esc(h.s.title || h.p.title)}</b><span>${esc(h.p.title)}${h.snip ? ' — …' + esc(h.snip).replace(esc(k), `<mark>${esc(k)}</mark>`) + '…' : ''}</span></span></button>`).join('') || '<p class="muted">لا نتائج. جرّب كلمة أقصر.</p>';
      $$('[data-hit]', out).forEach(b => b.onclick = () => { const h = hits[+b.dataset.hit]; close(); LEARN.open(h.p.key, h.si); });
      $$('[data-lab]', out).forEach(b => b.onclick = () => { close(); LABS.open(b.dataset.lab); });
    };
    q.oninput();
  });
}

/* ---------- home ---------- */
const RULES12 = ['الوضوح المفاهيمي', 'المراجعة الأدبية', 'تحديد الأبعاد', 'جدول المواصفات', 'فيض الفقرات', 'جودة الصياغة', 'سلم الاستجابة', 'صدق المحتوى', 'التجريب الاستطلاعي', 'تحليل الفقرات', 'الثبات', 'الصدق العاملي والتقنين'];
function heroX(rk, nx) {
  const lo = (RANKS.filter(r => r[0] <= S.xp).pop() || RANKS[0])[0], hi = rk.next ? rk.next[0] : lo + 1, pct = rk.next ? Math.min(100, (S.xp - lo) / (hi - lo) * 100) : 100;
  const lvl = RANKS.filter(r => r[0] <= S.xp).length, nScales = (typeof SCALES !== 'undefined') ? SCALES.list().length : 0;
  return `<section class="hx" aria-labelledby="hxT"><div class="hx-in">
    <p class="hx-kick"><span class="chip lime">${ic('star')} ${esc(rk.name)}</span><span>المستوى <b class="num">${lvl}</b> من ${RANKS.length}</span></p>
    <h1 class="hx-title" id="hxT">قِسْ ما لا يُرى، وابنِ مقياسك بيدك</h1>
    <p class="hx-sub">دروس قصيرة، ومعامل تحسب الصدق والثبات في ثوانٍ، ومعالج يبني مقياسك ويجمع بياناته ويحللها. كل خطوة تكسبك نقاطًا وتقربك من رتبة «عالِم قياس نفسي».</p>
    <div class="hx-cta"><button class="btn lime magnet" id="hCont">${ic('play')}${S.xp ? 'تابع' : 'ابدأ'}: ${esc(nx.short)}</button><button class="btn glass magnet" data-go3="build">${ic('build')}ابنِ مقياسك</button><button class="btn glass" data-go3="exam">${ic('timer')}تحدّي الامتحان</button></div>
    <div class="hx-stats"><div><b class="num" data-count="${S.xp}">${S.xp}</b><small>نقطة</small></div><div><b class="num" data-count="${S.streak}">${S.streak}</b><small>يوم متتالٍ</small></div><div><b class="num" data-count="${LEARN.doneCount()}">${LEARN.doneCount()}</b><small>من ${LEARN.total()} دروس</small></div><div><b class="num" data-count="${nScales}">${nScales}</b><small>مقياس بنيته</small></div></div>
    <div class="ruler" aria-hidden="true">${Array.from({ length: 61 }, (_, i) => `<i class="${i % 10 === 0 ? 'l' : i % 5 === 0 ? 'm' : ''}"></i>`).join('')}<span class="needle" style="left:${100 - pct}%"><b class="num">${S.xp}</b></span></div>
    <p class="ruler-cap">${rk.next ? `باقٍ <b class="num">${hi - S.xp}</b> نقطة لتصبح «${esc(rk.next[1])}»` : 'وصلت إلى أعلى رتبة. أنت عالِم قياس!'}</p>
  </div>
  <aside class="hx-live" aria-label="مثال حي: كلما زادت الاستجابات المتسقة ارتفع معامل ألفا">
    <div class="lq-head"><b>مقياس التسويف الأكاديمي</b><span>ن = <b class="lq-n num">0</b></span></div>
    ${['أؤجل المذاكرة حتى الليلة السابقة للامتحان', 'أنشغل بالهاتف عندما أجلس للمذاكرة', 'أترك المهام الصعبة إلى آخر لحظة', 'أفاجأ بقرب موعد التسليم'].map(t => `<div class="lq-row"><span>${t}</span><div class="lq-dots">${[1, 2, 3, 4, 5].map(() => '<i class="lq-dot"></i>').join('')}</div></div>`).join('')}
    <div class="lq-gauge"><svg viewBox="0 0 120 70" aria-hidden="true"><path d="M10 62a50 50 0 0 1 100 0" class="lq-bg"/><path d="M10 62a50 50 0 0 1 100 0" class="lq-arc"/></svg><div><small>ثبات ألفا كرونباخ</small><b class="lq-a num">0.00</b><small class="lq-ok">≥ 0.70 مقبول</small></div></div>
  </aside></section>`;
}
function journeyX() {
  const lit = (typeof BUILD !== 'undefined' && BUILD.bestSteps) ? BUILD.bestSteps() : 0;
  const pos = i => { const row = i < 6 ? 0 : 1, k = row ? i - 6 : i; return [row ? 90 + k * 164 : 910 - k * 164, row ? 235 : 75]; };
  const pts = RULES12.map((_, i) => pos(i));
  let d = `M ${pts[0][0]} ${pts[0][1]}`; for (let i = 1; i < 12; i++) { const [x, y] = pts[i], [px, py] = pts[i - 1]; d += i === 6 ? ` C ${px - 120} ${py}, ${px - 120} ${y}, ${x} ${y}` : ` L ${x} ${y}`; }
  return `<section class="card journey"><div class="row"><h3 class="grow" style="margin:0">${ic('trophy')} رحلة الباحث: 12 مستوى</h3><span class="chip ${lit === 12 ? 'good' : 'lime'}">${lit}/12</span></div>
    <p class="small muted" style="margin:0">كل مستوى قاعدة من قواعد بناء المقياس. أكمل خطواته في «ابنِ مقياسك» لتضيء الخريطة.</p>
    <svg class="jmap" viewBox="0 0 1000 310" role="img" aria-label="أكملت ${lit} من 12 مستوى"><path class="trail-bg" d="${d}"/><path class="trail" d="${d}" style="--lit:${lit / 11}"/>${pts.map(([x, y], i) => `<g class="jn ${i < lit ? 'lit' : i === lit ? 'now' : ''}" transform="translate(${x} ${y})"><circle r="30"/><text class="n" y="9">${i < lit ? '✓' : i + 1}</text><text class="l" y="${y > 100 ? 62 : -44}">${RULES12[i]}</text></g>`).join('')}<g class="jtoken" transform="translate(${pts[0][0]} ${pts[0][1]})"><circle r="17" class="t1"/><circle r="9" class="t2"/></g></svg>
    <ol class="jlist">${RULES12.map((r, i) => `<li class="${i < lit ? 'lit' : i === lit ? 'now' : ''}"><span>${i < lit ? '✓' : i + 1}</span>${r}</li>`).join('')}</ol>
    <div class="row"><button class="btn soft sm" data-go3="build">${ic('build')} ${lit ? 'تابع رحلتك' : 'ابدأ المستوى الأول'}</button><button class="btn ghost sm" id="hRules2">اقرأ القواعد</button></div></section>`;
}
VIEWS.home = el => {
  const t = today(), tx = S.hist[t] || 0, rk = rankOf(S.xp), nx = LEARN.next();
  const acctNudge = window.ACCT && ACCT.available() && !ACCT.signedIn();
  el.innerHTML = `
  ${heroX(rk, nx)}
  ${(() => { let k = null; try { k = JSON.parse(localStorage.getItem('know-v1') || 'null'); } catch (e) { } const n = k && k.done ? Object.keys(k.done).length : 0; return n >= 8 ? '' : `<button class="linkcard" data-go3="know" style="--c:var(--coral);--s:var(--bad-soft)"><span class="ic">${ic('heart')}</span><span><b>${n ? 'كمّل رحلة «اعرف نفسك» (' + n + ' من 8)' : 'ابدأ هنا: رحلة «اعرف نفسك»'}</b><span>اختيارية، حوالي 15 دقيقة بشكل لعبة: المزاج، والذاكرة، والتركيز، والحالة النفسية، وحياتك الرقمية (M6)، وعلاقتك بالذكاء الاصطناعي (M2).</span></span></button>`; })()}
  ${acctNudge ? `<div class="notice amber row"><div class="grow"><b>احفظ تقدمك على كل أجهزتك</b><span class="small">أنشئ حسابًا مجانيًا بالبريد الإلكتروني (اختياري): تحفظ نقاطك وسلسلتك ومقاييسك، وتظهر في لوحة المتصدرين.</span></div><button class="btn amber sm" data-go3="account">إنشاء حساب</button></div>` : ''}
  <div class="grid2">
    <section class="card"><h3>${ic('bolt')} مهام اليوم</h3>${S.quests.list.map(q => `<div class="quest"><span>${esc(q.t)} ${q.done ? '<span class="chip good">تم</span>' : ''}</span><span class="muted num small">+10</span><div class="pbar amber"><i style="width:${Math.min(100, q.p / q.n * 100)}%"></i></div></div>`).join('')}</section>
    <section class="card"><h3>${ic('book')} تقدمك في الدروس</h3><p class="muted small" style="margin:0 0 8px">${LEARN.doneCount()} من ${LEARN.total()} أجزاء مكتملة</p><div class="pbar good"><i style="width:${LEARN.doneCount() / LEARN.total() * 100}%"></i></div><div class="row" style="margin-top:12px"><button class="btn soft sm" data-go3="learn">كل الدروس</button><button class="btn ghost sm" data-go3="practice">تدرّب</button></div></section>
  </div>
  <div class="section-title"><h2>اختصارات</h2></div>
  <div class="grid3">
    ${[['labs', 'flask', 'المعامل الحسابية', 'ألفا، KR‑20، التجزئة النصفية، كابا، صدق المحتوى، تحليل الفقرات', '--primary', '--primary-soft'],
      ['build', 'build', 'ابنِ مقياسك', 'معالج من 12 خطوة، جمع البيانات مباشرة، وتحليل فوري', '--good', '--good-soft'],
      ['formulas', 'sigma', 'المعادلات', 'كل المعادلات برموزها ومتى تُستخدم وأمثلة محلولة', '--violet', '--violet-soft'],
      ['practice', 'cards', 'التدريب', 'بطاقات مراجعة، ومطابقة، ومسائل حسابية مولدة', '--amber', '--amber-soft'],
      ['bank', 'cards', 'بنك الأسئلة', 'أكثر من 200 سؤال بإجاباتها: صح وخطأ، واختيار، ومسائل، ومقالي', '--lime', '--lime-soft'],
      ['exam', 'check', 'أسئلة موضوعية والامتحان', 'صح وخطأ واختيار من متعدد بنمط امتحان القسم، ونموذج محلول', '--bad', '--bad-soft']]
      .map(c => `<button class="linkcard" data-go3="${c[0]}" style="--c:var(${c[4]});--s:var(${c[5]})"><span class="ic">${ic(c[1])}</span><span><b>${c[2]}</b><span>${c[3]}</span></span></button>`).join('')}
  </div>
  ${journeyX()}
  <p class="foot">إعداد الطالب: محمد عثمان — الدبلوم الخاص في علم النفس · إشراف: الدكتورة نجاة زكي</p>`;
  $('#hCont').onclick = () => LEARN.open(nx.key);
  $('#hRules2').onclick = () => LEARN.open('rules12');
  $$('[data-go3]', el).forEach(b => b.onclick = () => go(b.dataset.go3));
  if (typeof MOTION !== 'undefined') { MOTION.hero($('.hx', el)); MOTION.journey($('.journey', el)); MOTION.magnetic(el); }
};

/* ===== Learn: path, lessons, interactive stepper, quizzes, practice, profile ===== */
const LEARN = (() => {
  const ORDER = C.order;
  const PARTS = []; ORDER.forEach(([u, sub, keys], ui) => keys.forEach(k => { const p = C.parts[k]; if (p) { p.unit = ui; p.short = p.title.split(' — ')[0]; PARTS.push(p); } }));
  PARTS.forEach((p, i) => { p.idx = i; p.sections.forEach((s, j) => s.id = p.key + '.' + j); });
  const byKey = k => PARTS.find(p => p.key === k);
  const strip = h => h.replace(/<img[^>]*>/g, '').replace(/<[^>]+>/g, ' ').replace(/&amp;/g, '&').replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/\s+/g, ' ').trim();
  const text = s => s._t ?? (s._t = strip(s.html));
  const isLib = p => ['gloss', 'sch', 'refs', 'apx'].includes(p.key);

  /* ---- question pools ---- */
  C.bank.forEach((b, i) => b.id = 'b' + i);
  C.glossary.forEach((g, i) => g.id = 'g' + i);
  C.scholars.forEach((s, i) => s.id = 's' + i);
  C.mcq.forEach((m, i) => m.id = 'm' + i);
  const glossOK = C.glossary.filter(g => g.d && g.d.length > 8 && g.en);
  const qMcq = m => ({ id: m.id, kind: 'اختيار من متعدد', q: m.q, opts: shuffle(m.o), a: m.o[m.a], explain: m.w });
  const qBank = b => ({ id: b.id, kind: 'سؤال مراجعة', q: b.q, opts: shuffle([b.a, ...pick(C.bank.filter(x => x.a !== b.a), 3).map(x => x.a)]), a: b.a, explain: b.a });
  const qGlossEn = g => ({ id: g.id + 'e', kind: 'المقابل الإنجليزي', q: 'ما المقابل الإنجليزي لمصطلح «' + g.ar + '»؟', opts: shuffle([g.en, ...pick(glossOK.filter(x => x.en !== g.en), 3).map(x => x.en)]), a: g.en, ltr: true, explain: g.ar + ' = ' + g.en + '. ' + g.d });
  const qGlossDef = g => ({ id: g.id + 'd', kind: 'المصطلح من تعريفه', q: 'أي مصطلح يطابق هذا التعريف؟ «' + g.d + '»', opts: shuffle([g.ar, ...pick(glossOK.filter(x => x.ar !== g.ar), 3).map(x => x.ar)]), a: g.ar, explain: g.ar + ' (' + g.en + ')' });
  const qSch = s => ({ id: s.id, kind: 'من صاحب الإسهام؟', q: '«' + s.c + '» — من صاحب هذا الإسهام؟', opts: shuffle([s.ar, ...pick(C.scholars.filter(x => x.ar !== s.ar), 3).map(x => x.ar)]), a: s.ar, explain: s.ar + ' (' + s.en + (s.y ? '، ' + s.y : '') + ')' });
  const qQuick = k => ({ id: 'k' + k.k, kind: 'ورقة المراجعة', q: 'ما الجملة التي تلخص محور «' + k.k + '»؟', opts: shuffle([k.v, ...pick(C.quick.filter(x => x !== k), 3).map(x => x.v)]), a: k.v, explain: k.v });
  const pool = {};
  function poolFor(p) {
    if (pool[p.key]) return pool[p.key];
    let g = [];
    if (p.key === 'intro') g = C.quick.map(k => () => qQuick(k));
    else if (p.key === 'gloss') g = glossOK.flatMap(x => [() => qGlossEn(x), () => qGlossDef(x)]);
    else if (p.key === 'sch') g = C.scholars.filter(s => s.c).map(s => () => qSch(s));
    else if (p.key === 'refs' || p.key === 'apx') g = [];
    else {
      const t = p.sections.map(text).join(' ');
      g = C.mcq.filter(m => m.part === p.key).map(m => () => qMcq(m));
      g = g.concat(C.bank.filter(b => b.part === p.key).map(b => () => qBank(b)));
      glossOK.filter(x => x.ar.length > 3 && t.includes(x.ar)).forEach(x => { g.push(() => qGlossDef(x)); g.push(() => qGlossEn(x)); });
      C.scholars.filter(s => s.c && t.includes(s.ar)).forEach(s => g.push(() => qSch(s)));
    }
    return pool[p.key] = g;
  }
  const hasQuiz = p => poolFor(p).length >= 4;
  const readCount = p => p.sections.filter(s => S.read[s.id]).length;
  const partDone = p => hasQuiz(p) ? (S.quiz[p.key]?.best || 0) >= 0.6 && readCount(p) >= Math.ceil(p.sections.length / 2) : readCount(p) === p.sections.length;
  const progress = p => { const r = readCount(p) / p.sections.length; return hasQuiz(p) ? r * 0.5 + (S.quiz[p.key]?.best || 0) * 0.5 : r; };
  const learnParts = () => PARTS.filter(p => !isLib(p));
  const next = () => learnParts().find(p => !partDone(p)) || PARTS[0];

  /* ---- path view ---- */
  VIEWS.learn = el => {
    if (openKey) return renderLesson();
    const np = next();
    el.innerHTML = `<div class="section-title"><h2>الدروس</h2><span class="muted small">${doneCount()} من ${total()} أجزاء مكتملة</span></div>` +
      ORDER.map(([u, sub, keys], ui) => `<section class="stack"><div class="unit-h"><span class="n">${ui === 0 ? '٠' : ui === ORDER.length - 1 ? ic('book') : ui}</span><div><b>${esc(u)}</b><small>${esc(sub)}</small></div></div>` +
        keys.map(byKey).filter(Boolean).map(p => { const d = partDone(p), pr = progress(p), q = S.quiz[p.key]; return `<button class="partcard ${d ? 'done' : ''} ${p === np && !d ? 'current' : ''}" data-part="${p.key}"><span class="disc">${d ? ic('check') : p.idx + 1}</span><span><b>${esc(p.title)}</b><span class="meta">${p.sections.length} أقسام · ${readCount(p)} مقروءة${q ? ' · أفضل نتيجة ' + Math.round(q.best * 100) + '٪' : ''}</span><div class="pbar ${d ? 'good' : ''}"><i style="width:${pr * 100}%"></i></div></span><span class="muted">${ic('back')}</span></button>`; }).join('') + '</section>').join('');
    $$('[data-part]', el).forEach(b => b.onclick = () => open(b.dataset.part));
  };
  let openKey = null, openSec = null;
  function open(key, si) { openKey = key; openSec = si ?? null; if (location.hash === '#learn') renderLesson(); else location.hash = 'learn'; }
  function renderLesson() {
    const p = byKey(openKey); if (!p) return;
    const el = $('#view'); el.innerHTML = '<div class="container" id="vc"></div>'; const c = $('#vc');
    c.innerHTML = `<button class="back" id="lBack">${ic('back')} كل الدروس</button>
      <header class="card lesson-head"><span class="chip blue">${esc(ORDER[p.unit][0])}</span><h1>${esc(p.title)}</h1>${p.en ? `<div class="muted small ltr" style="text-align:right">${esc(p.en)}</div>` : ''}
      <div class="row small muted"><span>${readCount(p)}/${p.sections.length} أقسام مقروءة</span>${S.quiz[p.key] ? `<span class="chip good">أفضل نتيجة ${Math.round(S.quiz[p.key].best * 100)}٪</span>` : ''}</div><div class="pbar good"><i style="width:${progress(p) * 100}%"></i></div></header>
      <div class="stack" id="secs">${p.sections.map((s, i) => { const rd = !!S.read[s.id], op = openSec === i; const n = chunks(s).length; return `<article class="sec ${rd ? 'read' : ''}" id="sec-${i}"><div class="sec-h"><span class="ck">${rd ? ic('check') : ''}</span><span class="t">${esc(s.title || p.title)}<small>${n} بطاقات · ${rd ? 'مكتمل' : '≈ ' + Math.max(1, Math.round(text(s).length / 900)) + ' د قراءة'}</small></span></div><div class="row" style="padding:0 14px 12px"><button class="btn sm" data-play="${i}">${ic('play')}${rd ? 'أعد الدرس' : 'ابدأ الدرس'}</button><button class="btn sm ghost" data-sec="${i}">${op ? 'إخفاء النص' : 'قراءة كاملة'}</button></div>${op ? `<div class="body"><div class="prose">${decorate(s.html)}</div>${rd ? '' : `<button class="btn soft" data-read="${i}">${ic('check')} أنهيت القراءة (+5)</button>`}</div>` : ''}</article>`; }).join('')}</div>
      <div class="lfoot">${hasQuiz(p) ? `<button class="btn amber" id="lQuiz">${ic('bolt')} اختبر نفسك في هذا الجزء</button>` : ''}${p.idx < PARTS.length - 1 ? `<button class="btn ghost" id="lNext">الجزء التالي ${ic('back')}</button>` : ''}</div>`;
    $('#lBack').onclick = () => { openKey = null; route(); };
    $$('[data-sec]', c).forEach(b => b.onclick = () => { const i = +b.dataset.sec; openSec = openSec === i ? null : i; renderLesson(); const e = $('#sec-' + i); e && openSec != null && e.scrollIntoView({ block: 'start' }); });
    $$('[data-play]', c).forEach(b => b.onclick = () => startStepper(p, +b.dataset.play));
    $$('[data-read]', c).forEach(b => b.onclick = () => { const s = p.sections[+b.dataset.read]; markRead(s); openSec = +b.dataset.read + 1 < p.sections.length ? +b.dataset.read + 1 : null; renderLesson(); });
    const q = $('#lQuiz'); if (q) q.onclick = () => startPartQuiz(p);
    const nx = $('#lNext'); if (nx) nx.onclick = () => open(PARTS[p.idx + 1].key);
    if (openSec != null) { const e = $('#sec-' + openSec); e && e.scrollIntoView({ block: 'start' }); } else scrollTo(0, 0);
  }
  function markRead(s, xp = 5) { if (!S.read[s.id]) { S.read[s.id] = 1; quest('read', 1); addXP(xp); } }
  const GL = new Map(); C.glossary.forEach(g => { if (g.ar.length > 2 && !GL.has(g.ar)) GL.set(g.ar, g); });
  function decorate(h) { return h.replace(/<b>([^<]{3,60})<\/b>/g, (m, t) => { const g = GL.get(t.trim()); return g ? `<b class="gl" data-g="${g.id}" tabindex="0">${t}</b>` : m; }); }
  function showDef(id) { const g = C.glossary.find(x => x.id === id); if (!g) return; let d = $('#defsheet'); if (!d) { d = document.createElement('div'); d.id = 'defsheet'; d.className = 'defsheet'; document.body.appendChild(d); } d.hidden = false; d.innerHTML = `<div class="in"><span class="chip blue">من المسرد</span><b style="font-family:var(--f-display);font-size:18px">${esc(g.ar)}</b><div class="ltr muted small" style="text-align:right">${esc(g.en)}</div><div>${esc(g.d)}</div><button class="btn ghost sm" id="defx">إغلاق</button></div>`; $('#defx').onclick = () => d.hidden = true; }
  document.addEventListener('click', e => { const b = e.target.closest && e.target.closest('b.gl'); if (b) showDef(b.dataset.g); else { const d = $('#defsheet'); if (d && !d.hidden && !e.target.closest('#defsheet')) d.hidden = true; } });
  document.addEventListener('keydown', e => { if (e.key === 'Enter' && e.target.matches && e.target.matches('b.gl')) showDef(e.target.dataset.g); });

  /* ---- stepper (cards + retrieval checks) ---- */
  function chunks(sec) {
    if (sec._c) return sec._c;
    const t = document.createElement('template'); t.innerHTML = sec.html; const out = []; let cur = [], len = 0;
    const flush = () => { if (cur.length) { out.push(cur.map(e => e.outerHTML).join('')); cur = []; len = 0; } };
    [...t.content.children].forEach(e => {
      const tag = e.tagName, l = e.textContent.length;
      if (e.classList.contains('en')) { cur.push(e); return; }
      if (tag === 'H4') { flush(); cur.push(e); return; }
      const big = tag === 'ASIDE' || tag === 'FIGURE' || e.classList.contains('tw');
      if (big) { if (cur.some(x => x.tagName !== 'H4') && len > 200) flush(); cur.push(e); flush(); return; }
      if (len + l > 650 && len > 0) flush(); cur.push(e); len += l;
    });
    flush(); return sec._c = out.length ? out : [sec.html];
  }
  const terms = {};
  function checksFor(html, p) {
    if (!terms[p.key]) { const T = new Set(); p.sections.forEach(s => { const t = document.createElement('template'); t.innerHTML = s.html; t.content.querySelectorAll('p b,li b').forEach(b => { const x = b.textContent.trim(); if (x.length >= 3 && x.length <= 40 && !/[:：]$/.test(x)) T.add(x); }); }); terms[p.key] = [...T]; }
    const t = document.createElement('template'); t.innerHTML = html; t.content.querySelectorAll('.en').forEach(x => x.remove()); const out = [];
    t.content.querySelectorAll('p,li').forEach(el => {
      if (el.closest('table')) return;
      const txt = el.textContent.replace(/\s+/g, ' ').trim();
      el.querySelectorAll('b').forEach(b => {
        const term = b.textContent.trim();
        if (term.length < 3 || term.length > 40 || /[:：]$/.test(term) || txt.startsWith(term + ':') || /^[\d.\s=<>≥≤%()−-]+$/.test(term)) return;
        const s = txt.split(/(?<=[.؛!؟])\s+/).find(x => x.includes(term)); if (!s || s.length < term.length + 30 || s.length > 300) return;
        const others = terms[p.key].filter(x => x !== term && !x.includes(term) && !term.includes(x) && !s.includes(x)); if (others.length < 3) return;
        out.push({ kind: 'أكمل الفراغ', stem: s, term, opts: shuffle([term, ...pick(others.sort((a, b) => Math.abs(a.length - term.length) - Math.abs(b.length - term.length)).slice(0, 10), 3)]), a: term });
      });
    });
    return out;
  }
  let L = null;
  function startStepper(p, si) {
    const sec = p.sections[si], ch = chunks(sec), steps = []; let since = [], made = 0;
    ch.forEach((c, i) => { steps.push({ t: 'read', html: c }); since.push(c); if ((since.length >= 2 || i === ch.length - 1) && made < 4) { const cs = since.flatMap(x => checksFor(x, p)); if (cs.length) { steps.push({ t: 'check', q: cs[Math.random() * cs.length | 0] }); made++; } since = []; } });
    L = { p, si, sec, steps, i: 0, ok: 0, checks: made, xp: 0 }; stepShow();
  }
  const top = () => `<div class="qtop"><button class="x" id="lx" aria-label="إنهاء الدرس">${ic('x')}</button><div class="qprog"><i style="width:${L.i / L.steps.length * 100}%"></i></div><span class="pill xp num">+${L.xp}</span></div>`;
  function stepShow() {
    if (L.i >= L.steps.length) return stepDone();
    const st = L.steps[L.i];
    if (st.t === 'read') {
      const k = L.steps.slice(0, L.i + 1).filter(x => x.t === 'read').length, tot = L.steps.filter(x => x.t === 'read').length;
      openOv(top() + `<div class="lcard"><div class="row small" style="margin-bottom:6px"><span class="chip blue">بطاقة ${k} من ${tot}</span><span class="muted">${esc(L.sec.title)}</span></div><div class="prose">${decorate(st.html)}</div></div><div class="sheet"><button class="btn" id="lnext">${L.steps[L.i + 1] && L.steps[L.i + 1].t === 'check' ? 'فهمت، اختبرني' : 'تابع'}</button></div>`);
      $('#lnext').onclick = () => { L.i++; stepShow(); };
    } else {
      const q = st.q; let picked = null, done = false;
      const stem = esc(q.stem).split(esc(q.term)).join('<span class="blank">' + esc(q.term) + '</span>');
      openOv(top() + `<div class="qkind">اختبر نفسك · ${q.kind}</div><div class="qtext">${stem}</div><div class="opts">${q.opts.map((o, i) => `<button class="opt" data-o="${i}" aria-pressed="false"><span class="k">${i + 1}</span><span>${esc(o)}</span></button>`).join('')}</div><div class="sheet" id="lsheet"><button class="btn" id="lcheck" disabled>تحقّق</button></div>`);
      $$('.opt').forEach(b => b.onclick = () => { if (done) return; $$('.opt').forEach(x => x.setAttribute('aria-pressed', 'false')); b.setAttribute('aria-pressed', 'true'); picked = +b.dataset.o; $('#lcheck').disabled = false; });
      $('#lcheck').onclick = () => {
        if (done || picked == null) return; done = true; const ok = q.opts[picked] === q.a; S.answered++;
        $$('.opt').forEach((b, i) => { if (q.opts[i] === q.a) b.classList.add('ok'); else if (i === picked) b.classList.add('no'); });
        $$('.blank').forEach(b => { b.style.color = 'var(--good)'; b.style.borderBottomStyle = 'solid'; });
        if (ok) { L.ok++; S.correct++; L.xp += 3; } else { const it = { id: 'c' + L.sec.id + q.a, kind: q.kind, q: q.stem.split(q.term).join('＿＿＿'), opts: q.opts, a: q.a, explain: q.stem }; if (!S.miss.some(m => m.id === it.id)) S.miss.push({ id: it.id, q: it }); if (S.miss.length > 120) S.miss.shift(); }
        save();
        $('#lsheet').innerHTML = `<div class="fb ${ok ? 'ok' : 'no'}"><b>${ok ? 'إجابة صحيحة' : 'المعلومة الصحيحة:'}</b><div class="small">${esc(q.stem)}</div>${ok ? '' : '<div class="muted small">أُضيفت إلى «مراجعة الأخطاء» في التدريب.</div>'}</div><button class="btn" id="lnext">تابع</button>`;
        $('#lnext').onclick = () => { L.i++; stepShow(); }; $('#lnext').focus();
      };
    }
    $('#lx').onclick = closeOv;
  }
  function stepDone() {
    const first = !S.read[L.sec.id]; let g = L.xp + (first ? 5 : 2); if (L.checks >= 2 && L.ok === L.checks) g += 5;
    if (first) { S.read[L.sec.id] = 1; quest('read', 1); } addXP(g);
    const p = L.p, nx = p.sections[L.si + 1];
    openOv(`<div class="result"><div class="ring" style="width:84px;height:84px;display:grid;place-items:center;border-radius:50%;background:var(--good-soft);color:var(--good)">${ic('check')}</div><h2>اكتمل الدرس</h2><p class="muted" style="margin:0">${esc(L.sec.title)}</p><div class="tiles3"><div class="tile hl"><small>نقاط</small><b>+${g}</b></div><div class="tile good"><small>استرجاع صحيح</small><b>${L.ok}/${L.checks}</b></div><div class="tile"><small>بطاقات</small><b>${L.steps.filter(x => x.t === 'read').length}</b></div></div><div class="row" style="justify-content:center">${nx ? '<button class="btn" id="lnx">الدرس التالي</button>' : (hasQuiz(p) ? '<button class="btn amber" id="lqz">اختبر نفسك في الجزء</button>' : '')}<button class="btn ghost" id="lbk">العودة</button></div></div>`);
    if (nx) $('#lnx').onclick = () => startStepper(p, L.si + 1);
    const qz = $('#lqz'); if (qz) qz.onclick = () => startPartQuiz(p);
    $('#lbk').onclick = closeOv;
  }

  /* ---- overlay & quiz engine ---- */
  function openOv(html) { let o = $('#ov'); if (!o) { o = document.createElement('div'); o.id = 'ov'; o.className = 'ov'; o.innerHTML = '<div class="in" id="ovIn"></div>'; document.body.appendChild(o); } $('#ovIn').innerHTML = html; o.hidden = false; document.body.style.overflow = 'hidden'; o.scrollTop = 0; }
  function closeOv() { const o = $('#ov'); if (o) o.hidden = true; document.body.style.overflow = ''; clearInterval(Q && Q.tm); Q = null; L = null; if (openKey && current === 'learn') renderLesson(); else rerender(); }
  let Q = null;
  function startPartQuiz(p) { const g = poolFor(p); runQuiz({ title: p.short, qs: pick(g, Math.min(8, g.length)).map(f => f()), onEnd: r => { const prev = S.quiz[p.key]?.best || 0; S.quiz[p.key] = { best: Math.max(prev, r.acc), last: today() }; quest('quiz', 1); } }); }
  function runQuiz(cfg) { Q = Object.assign({ i: 0, ok: 0, combo: 0, best: 0, xp: 0, wrong: [] }, cfg); showQ(); }
  function showQ() {
    if (Q.timer && Date.now() >= Q.deadline) return endQuiz();
    if (Q.i >= Q.qs.length) { if (Q.endless) Q.qs.push(...Q.more()); else return endQuiz(); }
    const it = Q.qs[Q.i]; let picked = null, answered = false;
    openOv(`<div class="qtop"><button class="x" id="qx" aria-label="إنهاء">${ic('x')}</button><div class="qprog"><i style="width:${(Q.timer ? Math.max(0, (Q.deadline - Date.now()) / 60000) : Q.i / Q.qs.length) * 100}%"></i></div>${Q.timer ? '<span class="pill xp num" id="qtime"></span>' : `<span class="muted num">${Q.i + 1}/${Q.qs.length}</span>`}</div>
      <div class="qkind">${esc(it.kind)}${Q.combo >= 2 ? ` · <span style="color:var(--good)">سلسلة × ${Q.combo}</span>` : ''}</div><div class="qtext">${esc(it.q)}</div>
      <div class="opts">${it.opts.map((o, i) => `<button class="opt" data-o="${i}" aria-pressed="false" ${it.ltr ? 'dir="ltr"' : ''}><span class="k">${i + 1}</span><span>${esc(o)}</span></button>`).join('')}</div>
      <div class="sheet" id="sheet"><button class="btn" id="qcheck" disabled>تحقّق</button></div>`);
    $('#qx').onclick = () => Q.i > 0 && !Q.timer ? endQuiz(true) : closeOv();
    $$('.opt').forEach(b => b.onclick = () => { if (answered) return; $$('.opt').forEach(x => x.setAttribute('aria-pressed', 'false')); b.setAttribute('aria-pressed', 'true'); picked = +b.dataset.o; $('#qcheck').disabled = false; if (Q.timer) $('#qcheck').click(); });
    $('#qcheck').onclick = () => {
      if (answered || picked == null) return; answered = true; const ok = it.opts[picked] === it.a; S.answered++;
      $$('.opt').forEach((b, i) => { if (it.opts[i] === it.a) b.classList.add('ok'); else if (i === picked) b.classList.add('no'); });
      if (ok) { Q.ok++; Q.combo++; Q.best = Math.max(Q.best, Q.combo); S.correct++; Q.xp += Q.timer ? 5 : 10 + (Q.combo >= 3 ? 5 : 0); S.miss = S.miss.filter(m => m.id !== it.id); }
      else { Q.combo = 0; Q.wrong.push(it); if (!S.miss.some(m => m.id === it.id)) S.miss.push({ id: it.id, q: it }); if (S.miss.length > 120) S.miss.shift(); }
      save();
      if (Q.timer) { Q.i++; setTimeout(() => Q && showQ(), ok ? 250 : 700); return; }
      $('#sheet').innerHTML = `<div class="fb ${ok ? 'ok' : 'no'}"><b>${ok ? 'أحسنت' : 'الإجابة الصحيحة: ' + esc(it.a)}</b>${it.explain && it.explain !== it.a ? `<div class="small muted">${esc(it.explain)}</div>` : ''}</div><button class="btn" id="qnext">تابع</button>`;
      $('#qnext').onclick = () => { Q.i++; showQ(); }; $('#qnext').focus();
    };
    if (Q.timer) { const upd = () => { const s = Math.max(0, Math.ceil((Q.deadline - Date.now()) / 1000)); const e = $('#qtime'); if (e) e.textContent = s + 'ث'; if (s <= 0 && Q) endQuiz(); }; upd(); clearInterval(Q.tm); Q.tm = setInterval(upd, 500); }
  }
  function endQuiz(early) {
    clearInterval(Q.tm); const n = Q.timer ? Q.i : Q.qs.length, acc = n ? Q.ok / n : 0;
    let bonus = 0; if (!Q.timer && !early && Q.ok === Q.qs.length && n >= 4) { bonus = 15; S.perfect++; }
    Q.onEnd && Q.onEnd({ acc, ok: Q.ok, n });
    addXP(Q.xp + bonus);
    const stars = Math.round(acc * 3);
    openOv(`<div class="result"><h2>${acc >= 0.9 ? 'نتيجة ممتازة' : acc >= 0.6 ? 'اجتزت الاختبار' : 'جولة تدريبية جيدة'}</h2><div style="font-size:30px;color:var(--amber);letter-spacing:4px">${'★'.repeat(stars)}${'☆'.repeat(3 - stars)}</div><div class="tiles3"><div class="tile hl"><small>نقاط</small><b>+${Q.xp + bonus}</b></div><div class="tile good"><small>الدقة</small><b>${Math.round(acc * 100)}%</b></div><div class="tile"><small>أطول سلسلة</small><b>${Q.best}</b></div></div>${Q.wrong.length ? `<details class="steps" style="width:100%;text-align:start"><summary>راجع أخطاءك (${Q.wrong.length})</summary>${Q.wrong.map(w => `<p><b>${esc(w.q)}</b><br><span style="color:var(--good)">✓ ${esc(w.a)}</span></p>`).join('')}</details>` : ''}<button class="btn" id="qdone" style="min-width:220px">متابعة</button></div>`);
    $('#qdone').onclick = closeOv;
  }
  document.addEventListener('keydown', e => { const o = $('#ov'); if (!o || o.hidden) return; if (e.key === 'Escape') { closeOv(); return; } if (/^[1-4]$/.test(e.key)) { const b = $(`.opt[data-o="${+e.key - 1}"]`); b && b.click(); } if (e.key === 'Enter') { const n = $('#qnext') || $('#lnext') || $('#qcheck') || $('#lcheck'); if (n && !n.disabled) { e.preventDefault(); n.click(); } } });

  /* ---- practice ---- */
  const BOX = [0, 1, 2, 4, 8, 16, 32];
  const allCards = () => [...C.bank.map(b => ({ id: b.id, f: b.q, b: b.a })), ...C.quick.map(k => ({ id: 'k' + k.k, f: k.k, b: k.v })), ...glossOK.map(g => ({ id: g.id, f: g.ar, b: g.en + ' — ' + g.d }))];
  const dueCards = () => { const t = today(); return allCards().filter(c => !S.cards[c.id] || S.cards[c.id].due <= t); };
  VIEWS.practice = el => {
    const items = [['cards', 'cards', 'بطاقات المراجعة المتباعدة', 'تقلب البطاقة وتقيّم نفسك؛ ما تعرفه يعود بعد أيام.', dueCards().length + ' بطاقة مستحقة', '--primary', '--primary-soft'],
      ['probs', 'calc', 'مسائل حسابية مولّدة', 'ألفا، KR‑20، سبيرمان–براون، جتمان، كابا، CVR، الخطأ المعياري، الدرجة التائية.', 'أرقام جديدة كل مرة', '--good', '--good-soft'],
      ['match', 'match', 'سباق المطابقة', 'طابق المصطلح العربي بمقابله الإنجليزي خلال 60 ثانية.', '+2 لكل زوج', '--amber', '--amber-soft'],
      ['mistakes', 'redo', 'مراجعة الأخطاء', 'الأسئلة التي أخطأت فيها سابقًا حتى تتقنها.', S.miss.length + ' سؤال', '--bad', '--bad-soft'],
      ['rapid', 'timer', 'تحدّي الدقيقة', 'أكبر عدد من الإجابات الصحيحة في 60 ثانية.', 'من كل الأجزاء', '--violet', '--violet-soft'],
      ['scholars', 'users', 'من صاحب الإسهام؟', 'العلماء وإسهاماتهم في القياس النفسي.', '8 أسئلة', '--primary', '--primary-soft']];
    el.innerHTML = `<div class="section-title"><h2>التدريب</h2></div><div class="grid3">${items.map(([k, i, t, d, e, c, s]) => `<button class="linkcard" data-pr="${k}" style="--c:var(${c});--s:var(${s})"><span class="ic">${ic(i)}</span><span><b>${t}</b><span>${d}</span><span class="chip" style="margin-top:6px">${e}</span></span></button>`).join('')}</div>`;
    $$('[data-pr]', el).forEach(b => b.onclick = () => practice(b.dataset.pr));
  };
  function practice(k) {
    if (k === 'cards') return cardsMode();
    if (k === 'match') return matchMode();
    if (k === 'probs') return runQuiz({ qs: Array.from({ length: 6 }, genProblem), onEnd: r => { S.probs += r.ok; quest('prob', r.ok); } });
    if (k === 'scholars') return runQuiz({ qs: pick(C.scholars.filter(s => s.c), 8).map(qSch) });
    if (k === 'mistakes') { if (!S.miss.length) return toast('لا أخطاء محفوظة. أحسنت'); return runQuiz({ qs: shuffle(S.miss.map(m => m.q)).slice(0, 10).map(q => Object.assign({}, q, { opts: shuffle(q.opts) })) }); }
    if (k === 'rapid') { const all = PARTS.flatMap(p => poolFor(p)); const more = () => pick(all, 20).map(g => g()); return runQuiz({ qs: more(), more, endless: true, timer: true, deadline: Date.now() + 60000 }); }
  }
  function genProblem() {
    const R = (a, b, d = 0) => +(a + Math.random() * (b - a)).toFixed(d), types = ['alpha', 'kr20', 'sb', 'gut', 'kappa', 'cvr', 'sem', 'tscore', 'kr21', 'levels'];
    const t = types[Math.random() * types.length | 0];
    let q, ans, ex;
    if (t === 'alpha') { const k = R(5, 20), sv = R(k * 0.6, k * 1.0, 1), tv = R(sv * 2.5, sv * 5, 1); ans = PSY.alphaFromVars(k, sv, tv); q = `مقياس من ${k} عبارة، مجموع تباينات العبارات = ${sv}، وتباين الدرجة الكلية = ${tv}. احسب ألفا كرونباخ.`; ex = `α = (${k}/${k - 1}) × (1 − ${sv}/${tv})`; }
    else if (t === 'kr20') { const k = R(10, 30), spq = R(k * 0.12, k * 0.22, 2), tv = R(spq * 3, spq * 6, 1); ans = PSY.alphaFromVars(k, spq, tv); q = `اختبار صح/خطأ من ${k} فقرة، Σpq = ${spq}، وتباين الدرجات = ${tv}. احسب KR‑20.`; ex = `KR‑20 = (${k}/${k - 1}) × (1 − ${spq}/${tv})`; }
    else if (t === 'kr21') { const k = R(10, 30), M = R(k * 0.4, k * 0.7, 1), v = R(k * 0.6, k * 1.4, 1); ans = PSY.kr21(k, M, v); q = `اختبار ثنائي من ${k} فقرة، متوسط الدرجات ${M} وتباينها ${v}. احسب KR‑21.`; ex = `KR‑21 = (${k}/${k - 1}) × [1 − ${M}×(${k}−${M}) ÷ (${k}×${v})]`; }
    else if (t === 'sb') { const r = R(0.4, 0.85, 2); ans = PSY.spearmanBrown(r); q = `معامل الارتباط بين نصفي اختبار = ${r}. احسب معامل الثبات بمعادلة سبيرمان–براون.`; ex = `2 × ${r} ÷ (1 + ${r})`; }
    else if (t === 'gut') { const v1 = R(3, 9, 1), v2 = R(3, 9, 1), vt = +(v1 + v2 + 2 * R(0.6, 0.9, 2) * Math.sqrt(v1 * v2)).toFixed(1); ans = PSY.guttman(v1, v2, vt); q = `تباين النصف الأول ${v1}، والثاني ${v2}، والدرجة الكلية ${vt}. احسب معامل جتمان.`; ex = `2 × [1 − (${v1} + ${v2}) ÷ ${vt}]`; }
    else if (t === 'kappa') { const po = R(0.7, 0.95, 2), pe = R(0.4, 0.6, 2); ans = (po - pe) / (1 - pe); q = `نسبة الاتفاق الملاحظ بين مقدِّرين ${po}، والمتوقع بالصدفة ${pe}. احسب كابا كوهين.`; ex = `(${po} − ${pe}) ÷ (1 − ${pe})`; }
    else if (t === 'cvr') { const N = R(8, 20), ne = R(Math.ceil(N / 2), N); ans = PSY.cvr(ne, N); q = `من ${N} محكمين قال ${ne} إن العبارة «أساسية». احسب نسبة صدق المحتوى (CVR) للوشي.`; ex = `(${ne} − ${N}/2) ÷ (${N}/2)`; }
    else if (t === 'sem') { const s = R(5, 20, 1), r = R(0.7, 0.95, 2); ans = PSY.sem(s, r); q = `الانحراف المعياري للدرجات ${s} ومعامل الثبات ${r}. احسب الخطأ المعياري للقياس.`; ex = `${s} × √(1 − ${r})`; }
    else if (t === 'tscore') { const m = R(40, 120), s = R(5, 20), x = R(m - 2 * s, m + 2 * s); ans = PSY.tScore((x - m) / s); q = `متوسط ${m} وانحراف معياري ${s}. ما الدرجة التائية T لمن حصل على ${x}؟`; ex = `z = (${x} − ${m}) ÷ ${s}؛ T = 50 + 10z`; }
    else { const lo = 1, hi = R(4, 7), n = 3; ans = (hi - lo) / n; q = `تدريج من ${lo} إلى ${hi}. ما طول الفئة لتقسيم المتوسطات إلى ${n} مستويات؟`; ex = `(${hi} − ${lo}) ÷ ${n}`; }
    const d = t === 'tscore' ? 1 : 2, A = fmt(ans, d), deltas = t === 'tscore' ? [-3.7, 2.4, 5.1] : [-0.07, 0.05, 0.12];
    const opts = shuffle([A, ...deltas.map(x => fmt(ans + x * (t === 'tscore' ? 1 : (Math.abs(ans) > 1 ? ans : 1)), d))].filter((v, i, a) => a.indexOf(v) === i)).slice(0, 4);
    if (!opts.includes(A)) opts[0] = A;
    return { id: 'pr' + t, kind: 'مسألة حسابية', q, opts, a: A, ltr: true, explain: ex + ' = ' + A };
  }
  function cardsMode() {
    const deck = shuffle(dueCards()).slice(0, 20); if (!deck.length) return toast('لا بطاقات مستحقة اليوم. عُد غدًا');
    let i = 0, known = 0;
    const show = () => {
      if (i >= deck.length) { addXP(deck.length); quest('cards', 0); openOv(`<div class="result"><h2>أنهيت المراجعة</h2><div class="tiles3"><div class="tile"><small>بطاقات</small><b>${deck.length}</b></div><div class="tile good"><small>عرفتها</small><b>${known}</b></div><div class="tile hl"><small>نقاط</small><b>+${deck.length}</b></div></div><button class="btn" id="qdone" style="min-width:220px">متابعة</button></div>`); $('#qdone').onclick = closeOv; return; }
      const c = deck[i], st = S.cards[c.id];
      openOv(`<div class="qtop"><button class="x" id="qx" aria-label="إنهاء">${ic('x')}</button><div class="qprog"><i style="width:${i / deck.length * 100}%"></i></div><span class="muted num">${i + 1}/${deck.length}</span></div><div class="qkind">بطاقة مراجعة · الصندوق ${(st && st.b) || 0}</div><div class="card3d"><button class="flip" id="flip" aria-label="اقلب البطاقة"><div class="face"><b>${esc(c.f)}</b><small class="muted">اضغط لترى الإجابة</small></div><div class="face back"><span>${esc(c.b)}</span></div></button></div><div class="sheet row" id="cb" style="justify-content:center"><button class="btn ghost" id="flipBtn">اقلب البطاقة</button></div>`);
      const flip = () => { $('#flip').classList.add('on'); $('#cb').innerHTML = '<button class="btn danger" id="no" style="flex:1">لم أعرفها</button><button class="btn" id="yes" style="flex:1">عرفتها</button>'; $('#yes').onclick = () => rate(true); $('#no').onclick = () => rate(false); };
      $('#flip').onclick = flip; $('#flipBtn').onclick = flip; $('#qx').onclick = closeOv;
    };
    const rate = ok => { const c = deck[i], s = S.cards[c.id] || { b: 0 }; s.b = ok ? Math.min(6, s.b + 1) : 0; const d = new Date(); d.setDate(d.getDate() + BOX[s.b]); s.due = d.getFullYear() + '-' + String(d.getMonth() + 1).padStart(2, '0') + '-' + String(d.getDate()).padStart(2, '0'); S.cards[c.id] = s; if (ok) known++; else deck.push(c); quest('cards', 1); save(); i++; show(); };
    show();
  }
  let MT = null;
  function matchMode() {
    const pool = shuffle(glossOK.filter(g => g.en.length < 34 && g.ar.length < 30)); let k = 0, pairs = 0, sel = null; const end0 = Date.now() + 60000; let end = end0;
    const round = () => {
      const set = pool.slice(k, k + 5); k += 5; if (set.length < 5) { k = 0; return round(); }
      const Ls = shuffle(set.map(g => ({ id: g.id, t: g.ar }))), Rs = shuffle(set.map(g => ({ id: g.id, t: g.en }))); let left = 5;
      openOv(`<div class="qtop"><button class="x" id="qx" aria-label="إنهاء">${ic('x')}</button><div class="qprog"><i id="mbar" style="width:100%;background:var(--amber)"></i></div><span class="pill xp num" id="mt"></span></div><div class="qkind">طابق المصطلح بمقابله · <span class="num" id="mp">${pairs}</span> زوجًا</div><div class="mgrid" style="margin-top:12px">${Ls.map((x, i) => `<button class="opt" data-side="L" data-id="${x.id}" style="grid-column:1;grid-row:${i + 1}">${esc(x.t)}</button>`).join('')}${Rs.map((x, i) => `<button class="opt" dir="ltr" data-side="R" data-id="${x.id}" style="grid-column:2;grid-row:${i + 1}">${esc(x.t)}</button>`).join('')}</div>`);
      $('#qx').onclick = finish;
      $$('.mgrid .opt').forEach(b => b.onclick = () => {
        if (!sel || sel.dataset.side === b.dataset.side) { sel && sel.setAttribute('aria-pressed', 'false'); sel = b; b.setAttribute('aria-pressed', 'true'); return; }
        if (sel.dataset.id === b.dataset.id) { [sel, b].forEach(x => { x.classList.add('ok'); setTimeout(() => x.classList.add('gone'), 250); }); pairs++; $('#mp').textContent = pairs; if (--left === 0) setTimeout(round, 300); }
        else { [sel, b].forEach(x => { x.classList.add('no'); setTimeout(() => x.classList.remove('no'), 400); }); end -= 3000; }
        sel.setAttribute('aria-pressed', 'false'); sel = null;
      });
    };
    const finish = () => { clearInterval(MT); MT = null; if (pairs) { S.matches++; addXP(pairs * 2); } openOv(`<div class="result"><h2>انتهى السباق</h2><div class="tiles3"><div class="tile hl"><small>أزواج</small><b>${pairs}</b></div><div class="tile"><small>نقاط</small><b>+${pairs * 2}</b></div><div class="tile"><small>الزمن</small><b>60ث</b></div></div><p class="muted small">كل خطأ يخصم 3 ثوانٍ.</p><button class="btn" id="qdone" style="min-width:220px">متابعة</button></div>`); $('#qdone').onclick = closeOv; };
    round(); clearInterval(MT); MT = setInterval(() => { const s = Math.max(0, (end - Date.now()) / 1000), t = $('#mt'), b = $('#mbar'); if (t) t.textContent = Math.ceil(s) + 'ث'; if (b) b.style.width = (s / 60 * 100) + '%'; if (s <= 0) finish(); }, 250);
  }

  /* ---- profile ---- */
  VIEWS.me = el => {
    const rk = rankOf(S.xp), acc = S.answered ? Math.round(S.correct / S.answered * 100) : 0;
    const days = []; for (let i = 6; i >= 0; i--) { const d = new Date(); d.setDate(d.getDate() - i); const k = d.getFullYear() + '-' + String(d.getMonth() + 1).padStart(2, '0') + '-' + String(d.getDate()).padStart(2, '0'); days.push([S.hist[k] || 0, d.toLocaleDateString('ar-EG', { weekday: 'short' })]); }
    el.innerHTML = `<div class="section-title"><h2>إنجازاتي</h2><span class="chip amber">${ic('star')} ${esc(rk.name)}</span></div>
    <div class="kpis">${[['النقاط', S.xp], ['السلسلة الحالية', S.streak + ' يوم'], ['أطول سلسلة', S.best + ' يوم'], ['الدقة', acc + '%'], ['دروس مكتملة', Object.keys(S.read).length], ['أجزاء منجزة', doneCount() + '/' + total()]].map(([a, b]) => `<div class="tile"><small>${a}</small><b>${b}</b></div>`).join('')}</div>
    <section class="card"><h3>${ic('chart')} نقاطك في آخر سبعة أيام</h3>${barSVG(days.map(d => d[0]), days.map(d => d[1]), { d: 0, ref: S.goal, refLabel: 'الهدف ' + S.goal, max: S.goal, label: 'النقاط اليومية', color: v => v >= S.goal ? 'var(--amber)' : 'var(--primary)' })}</section>
    <section class="card"><h3>${ic('trophy')} الأوسمة (${S.badges.length}/${BADGES.length})</h3><div class="badges">${BADGES.map(b => `<div class="badge ${S.badges.includes(b.id) ? '' : 'off'}"><div class="m">${ic('star')}</div><b>${b.t}</b><small>${b.d}</small></div>`).join('')}</div></section>
    <section class="card"><h3>${ic('gear')} الإعدادات</h3><div class="stack">
      <div class="row"><span class="grow">هدف النقاط اليومي</span><div class="seg">${[20, 40, 80].map(g => `<button data-goal="${g}" aria-pressed="${S.goal === g}">${g}</button>`).join('')}</div></div>
      <div class="row"><span class="grow">حجم الخط</span><div class="seg">${[[1, 'عادي'], [1.1, 'كبير'], [1.2, 'أكبر']].map(([z, t]) => `<button data-fs="${z}" aria-pressed="${(S.fs || 1) === z}">${t}</button>`).join('')}</div></div>
      <div class="row"><span class="grow">مسح كل التقدم على هذا الجهاز</span><button class="btn sm danger" id="reset">مسح</button></div></div></section>`;
    $$('[data-goal]', el).forEach(b => b.onclick = () => { S.goal = +b.dataset.goal; save(); rerender(); });
    $$('[data-fs]', el).forEach(b => b.onclick = () => { S.fs = +b.dataset.fs; applyPrefs(); save(); rerender(); });
    $('#reset').onclick = () => confirmBox('سيُمسح تقدمك ونقاطك وسلسلتك على هذا الجهاز نهائيًا.', 'نعم، امسح', () => { S = fresh(); tick(); save(); toast('بدأت من جديد'); rerender(); });
  };

  const doneCount = () => learnParts().filter(partDone).length, total = () => learnParts().length;
  return { parts: () => PARTS, text, open, next: () => { const p = next(); return { key: p.key, short: p.short }; }, doneCount, total, runQuiz, openOv, closeOv, renderLesson, isOpen: () => !!openKey, clearOpen: () => { openKey = null; } };
})();

/* ===== Labs (calculators) and formulas reference ===== */
const LABS = (() => {
  const EX_LIK = `ع1 ع2 ع3 ع4 ع5
5 4 5 4 5
4 4 5 5 4
2 1 2 2 1
3 3 3 4 3
5 5 4 5 5
1 2 2 1 2
4 5 4 4 4
2 2 1 2 3
3 4 3 3 3
4 3 4 4 5`;
  const EX_BIN = `1 1 1 1 1 0
1 1 1 0 0 0
1 0 1 0 0 0
1 1 1 1 1 1
1 0 0 0 0 0
1 1 1 1 0 0
0 0 0 0 0 1
1 1 1 1 1 0`;
  const EX_BIG = (() => { /* simulated two-dimension Likert data for practice (seeded) */
    let s = 11; const r = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const g = () => Math.sqrt(-2 * Math.log(r() || 1e-9)) * Math.cos(2 * Math.PI * r());
    const rows = []; for (let i = 0; i < 60; i++) { const f1 = g(), f2 = 0.3 * f1 + 0.95 * g(); const row = []; for (let j = 0; j < 8; j++) { const f = j < 4 ? f1 : f2, x = 3 + 1.1 * (0.75 * f + 0.6 * g()); row.push(Math.max(1, Math.min(5, Math.round(x)))); } rows.push(row.join(' ')); }
    return rows.join('\n');
  })();

  const list = [
    { id: 'alpha', t: 'ألفا كرونباخ وتحليل العبارات', d: 'من البيانات الخام أو من التباينات: ألفا، والارتباط المصحح، وألفا عند حذف العبارة.', ic: 'calc' },
    { id: 'kr', t: 'كودر–ريتشاردسون 20 و21', d: 'ثبات الاختبارات ثنائية الدرجة (صح/خطأ).', ic: 'calc' },
    { id: 'split', t: 'التجزئة النصفية', d: 'سبيرمان–براون وجتمان من البيانات أو من المعاملات.', ic: 'calc' },
    { id: 'retest', t: 'ثبات إعادة الاختبار', d: 'ارتباط بيرسون بين تطبيقين مع رسم الانتشار.', ic: 'chart' },
    { id: 'kappa', t: 'كابا كوهين', d: 'اتفاق مقدِّرين على تصنيف اسمي (2×2 أو 3×3 أو 4×4).', ic: 'users' },
    { id: 'content', t: 'صدق المحتوى', d: 'CVR للوشي بالقيمة الحرجة، وCVI للين، ومعامل V لأيكن.', ic: 'shield' },
    { id: 'items', t: 'تحليل الفقرات', d: 'الصعوبة والتمييز للفقرات الثنائية، والمجموعتان الطرفيتان لفقرات ليكرت.', ic: 'chart' },
    { id: 'scores', t: 'الدرجات والمستويات', d: 'الدرجة الزائية والتائية والمئين، والخطأ المعياري، والمتوسط النظري وطول الفئة.', ic: 'sigma' },
    { id: 'omega', t: 'أوميغا والصدق التقاربي', d: 'أوميغا ماكدونالد ومتوسط التباين المستخلص من التشبعات.', ic: 'sigma' },
    { id: 'dif', t: 'الأداء التفاضلي (مانتل–هانزل)', d: 'نسبة الأرجحية ودلتا وتصنيف هيئة الاختبارات التربوية.', ic: 'users' },
    { id: 'efa', t: 'البنية العاملية الاستكشافية', d: 'KMO وبارتلت والجذور الكامنة والتحليل الموازي والتشبعات بعد التدوير.', ic: 'chart' },
    { id: 'sample', t: 'حجم العينة', d: 'أحجام العينات المتداولة لكل مرحلة من مراحل بناء المقياس.', ic: 'users' }
  ];
  let openId = null;
  VIEWS.labs = el => {
    if (openId) return renderLab(el, openId);
    el.innerHTML = `<div class="section-title"><h2>المعامل الحسابية</h2><span class="muted small">كل معمل يحسب النتيجة ويعرض خطوات الحل والتفسير وصياغة النتيجة</span></div>
      <div class="grid3">${list.map(l => `<button class="linkcard" data-lab="${l.id}"><span class="ic">${ic(l.ic)}</span><span><b>${l.t}</b><span>${l.d}</span>${S.labs[l.id] ? '<span class="chip good" style="margin-top:6px">استخدمته</span>' : ''}</span></button>`).join('')}</div>
      <div class="notice small"><b>طريقة إدخال البيانات</b>ضع كل مستجيب في سطر، وافصل بين درجات العبارات بمسافة أو فاصلة أو Tab (يمكن اللصق مباشرة من Excel أو SPSS). يُقبل سطر عناوين في البداية، وتُقبل الأرقام العربية.</div>`;
    $$('[data-lab]', el).forEach(b => b.onclick = () => open(b.dataset.lab));
  };
  function open(id) { openId = id; if (location.hash === '#labs') renderLab($('#vc'), id); else location.hash = 'labs'; }
  const head = (l) => `<button class="back" id="labBack">${ic('back')} كل المعامل</button><div class="section-title"><h2>${esc(l.t)}</h2></div>`;
  const used = {};
  function mark(id) { if (!used[id]) { used[id] = 1; labUse(id); } }
  function renderLab(el, id) {
    const l = list.find(x => x.id === id); if (!l) { openId = null; return VIEWS.labs(el); }
    el.innerHTML = head(l) + '<div id="labBody" class="stack"></div>';
    $('#labBack').onclick = () => { openId = null; rerender(); };
    R[id]($('#labBody'));
  }
  const matrixInput = (id, val, help) => `<div class="field"><label class="fl" for="${id}">${help || 'البيانات: كل مستجيب في سطر، والأعمدة عبارات'}</label><textarea id="${id}" class="data" spellcheck="false">${esc(val)}</textarea></div>`;
  const out = (el, html) => { el.innerHTML = html; };
  const apa = s => `<div class="apa"><b>صياغة النتيجة:</b> ${s} <button class="btn sm ghost" data-copy>${ic('copy')} نسخ</button></div>`;
  const bindCopy = el => $$('[data-copy]', el).forEach(b => b.onclick = () => copyText(b.parentElement.textContent.replace('صياغة النتيجة:', '').replace('نسخ', '').trim()));
  const num = (id, label, val, step = 'any') => `<div class="field"><label class="fl" for="${id}">${label}</label><input type="number" id="${id}" value="${val}" step="${step}"></div>`;
  const v = id => toNum($('#' + id).value);
  const seg = (name, opts, cur) => `<div class="seg" role="group">${opts.map(([k, t]) => `<button data-${name}="${k}" aria-pressed="${k === cur}">${t}</button>`).join('')}</div>`;
  const dec = x => fmt(x, 3).replace(/^0\./, '.').replace(/^-0\./, '-.');

  const R = {};
  /* ---------- alpha ---------- */
  R.alpha = el => {
    let mode = 'raw';
    const draw = () => {
      el.innerHTML = `<div class="card stack">${seg('m', [['raw', 'من البيانات الخام'], ['var', 'من التباينات']], mode)}<div id="in"></div></div><div id="res" class="stack"></div>`;
      $$('[data-m]', el).forEach(b => b.onclick = () => { mode = b.dataset.m; draw(); });
      const inp = $('#in', el), res = $('#res', el);
      if (mode === 'var') {
        inp.innerHTML = `<div class="grid3">${num('k', 'عدد العبارات k', 10)}${num('sv', 'مجموع تباينات العبارات Σσ²ᵢ', 8)}${num('tv', 'تباين الدرجة الكلية σ²ₓ', 30)}</div>`;
        const run = () => { const k = v('k'), sv = v('sv'), tv = v('tv'); if (!(k > 1 && tv > 0)) return out(res, '<div class="notice warn">أدخل k أكبر من 1 وتباينًا كليًا موجبًا.</div>'); const a = PSY.alphaFromVars(k, sv, tv), [c, t] = relVerdict(a); mark('alpha');
          out(res, `<div class="kpis"><div class="tile hl"><small>ألفا كرونباخ</small><b>${fmt(a, 3)}</b></div></div><div class="verdict ${c}">${t}</div><div class="f">α = (${k} ÷ ${k - 1}) × (1 − ${sv} ÷ ${tv}) = ${fmt(k / (k - 1), 3)} × ${fmt(1 - sv / tv, 3)} = ${fmt(a, 3)}</div>`); };
        $$('input', inp).forEach(i => i.oninput = run); run(); return;
      }
      inp.innerHTML = matrixInput('mx', EX_LIK) + `<div class="row"><button class="btn sm ghost" data-ex="lik">مثال ليكرت (10 × 5)</button><button class="btn sm ghost" data-ex="big">بيانات محاكاة (60 × 8)</button></div>`;
      const ta = $('#mx', el);
      $$('[data-ex]', el).forEach(b => b.onclick = () => { ta.value = b.dataset.ex === 'big' ? EX_BIG : EX_LIK; run(); });
      const run = () => {
        const m = parseMatrix(ta.value); if (m.error) return out(res, `<div class="notice warn">${m.error}</div>`);
        if (m.n < 3 || m.k < 2) return out(res, '<div class="notice warn">تحتاج 3 مستجيبين وعبارتين على الأقل.</div>');
        const s = PSY.itemStats(m.rows), [c, t] = relVerdict(s.alpha), names = m.head && m.head.length === m.k ? m.head : s.itemMean.map((_, i) => 'ع' + (i + 1));
        const sh = PSY.splitHalf(m.rows); mark('alpha');
        out(res, `<div class="kpis"><div class="tile hl"><small>ألفا كرونباخ</small><b>${fmt(s.alpha, 3)}</b></div><div class="tile"><small>المستجيبون n</small><b>${m.n}</b></div><div class="tile"><small>العبارات k</small><b>${m.k}</b></div><div class="tile"><small>متوسط الدرجة الكلية</small><b>${fmt(s.totalMean)}</b></div><div class="tile"><small>الانحراف المعياري</small><b>${fmt(s.totalSd)}</b></div><div class="tile"><small>سبيرمان–براون</small><b>${fmt(sh.sb, 3)}</b></div></div>
          <div class="verdict ${c}">الثبات ${t}${s.alpha > 0.95 ? ' — قيمة مرتفعة جدًا قد تعني تكرار العبارات' : ''}</div>
          <section class="card"><h3>${ic('chart')} إحصاءات العبارات</h3>${tableHTML(['العبارة', 'المتوسط', 'الانحراف', 'الارتباط المصحح', 'ألفا عند الحذف', 'الحكم'], names.map((nm, j) => { const bad = s.rit[j] < 0.3, up = s.alphaIfDeleted[j] > s.alpha + 0.005; return { cls: bad ? 'flag' : '', cells: [esc(nm), N(s.itemMean[j]), N(s.itemSd[j]), N(s.rit[j]), N(s.alphaIfDeleted[j], 3), bad ? '<span class="chip bad">راجع العبارة (r < 0.30)</span>' : up ? '<span class="chip warn">حذفها يرفع ألفا</span>' : '<span class="chip good">جيدة</span>'] }; }))}
          ${barSVG(s.rit, names, { ref: 0.3, refLabel: '0.30', max: 1, label: 'الارتباط المصحح', color: x => x < 0.3 ? 'var(--bad)' : 'var(--primary)' })}</section>
          <details class="steps"><summary>خطوات الحل</summary><div class="f">Σσ²ᵢ = ${s.itemVar.map(x => fmt(x)).join(' + ')} = ${fmt(s.sumItemVar, 3)}</div><div class="f">σ²ₓ = ${fmt(s.totalVar, 3)}</div><div class="f">α = (${m.k} ÷ ${m.k - 1}) × (1 − ${fmt(s.sumItemVar, 3)} ÷ ${fmt(s.totalVar, 3)}) = ${fmt(s.alpha, 3)}</div><p class="small muted">التباينات محسوبة بقسمة مجموع المربعات على (n − 1) كما في SPSS. الارتباط المصحح هو ارتباط العبارة بمجموع باقي العبارات بعد استبعادها.</p></details>
          ${apa(`بلغ معامل ألفا كرونباخ للمقياس (ن = ${m.n}، ${m.k} عبارات) α = ${dec(s.alpha)}، وهي قيمة تدل على ثبات ${t.split(' ')[0]}، وتراوحت معاملات الارتباط المصححة بين العبارات والدرجة الكلية بين ${dec(Math.min(...s.rit))} و${dec(Math.max(...s.rit))}.`)}`);
        bindCopy(res);
      };
      ta.oninput = run; run();
    };
    draw();
  };
  /* ---------- KR ---------- */
  R.kr = el => {
    el.innerHTML = `<div class="card stack"><h3>KR‑20 من البيانات الثنائية (0/1)</h3>${matrixInput('mx', EX_BIN, 'كل مستجيب في سطر: 1 = صحيح/موافق، 0 = خطأ')}</div><div id="r20" class="stack"></div>
      <div class="card stack"><h3>KR‑21 من المتوسط والتباين</h3><div class="grid3">${num('k', 'عدد الفقرات k', 20)}${num('M', 'متوسط الدرجات M', 12)}${num('vv', 'تباين الدرجات σ²', 16)}</div><div id="r21"></div></div>`;
    const run20 = () => {
      const m = parseMatrix($('#mx', el).value), r = $('#r20', el); if (m.error) return out(r, `<div class="notice warn">${m.error}</div>`);
      if (m.rows.some(x => x.some(y => y !== 0 && y !== 1))) return out(r, '<div class="notice warn">KR‑20 للفقرات الثنائية فقط: استخدم 0 و1. لفقرات ليكرت استخدم معمل ألفا.</div>');
      const k = PSY.kr20(m.rows), [c, t] = relVerdict(k.value), k21 = PSY.kr21(k.k, k.mean, k.totalVar); mark('kr');
      out(r, `<div class="kpis"><div class="tile hl"><small>KR‑20</small><b>${fmt(k.value, 3)}</b></div><div class="tile"><small>KR‑21 (للمقارنة)</small><b>${fmt(k21, 3)}</b></div><div class="tile"><small>Σpq</small><b>${fmt(k.spq, 3)}</b></div><div class="tile"><small>تباين الدرجات</small><b>${fmt(k.totalVar, 3)}</b></div></div><div class="verdict ${c}">الثبات ${t}</div>
        ${tableHTML(['الفقرة', 'p (نسبة الصواب)', 'q = 1 − p', 'p × q'], k.p.map((p, j) => ['ف' + (j + 1), N(p, 3), N(1 - p, 3), N(p * (1 - p), 3)]))}
        <details class="steps"><summary>خطوات الحل</summary><div class="f">KR‑20 = (${k.k} ÷ ${k.k - 1}) × (1 − ${fmt(k.spq, 3)} ÷ ${fmt(k.totalVar, 3)}) = ${fmt(k.value, 3)}</div><p class="small muted">يُحسب تباين الدرجات الكلية هنا بالقسمة على n ليتسق مع p × q (وهما تباينان على مستوى المجتمع)، وبهذا يساوي KR‑20 قيمة ألفا للبيانات الثنائية.</p></details>`);
    };
    const run21 = () => { const k = v('k'), M = v('M'), vv = v('vv'), r = $('#r21', el); if (!(k > 1 && vv > 0)) return out(r, '<div class="notice warn">أدخل قيمًا صحيحة.</div>'); const x = PSY.kr21(k, M, vv); mark('kr'); out(r, `<div class="kpis"><div class="tile hl"><small>KR‑21</small><b>${fmt(x, 3)}</b></div></div><div class="f">KR‑21 = (${k} ÷ ${k - 1}) × [1 − ${M} × (${k} − ${M}) ÷ (${k} × ${vv})] = ${fmt(x, 3)}</div>`); };
    $('#mx', el).oninput = run20; $$('#k,#M,#vv', el).forEach(i => i.oninput = run21); run20(); run21();
  };
  /* ---------- split ---------- */
  R.split = el => {
    el.innerHTML = `<div class="card stack"><h3>من البيانات</h3>${matrixInput('mx', EX_LIK)}<div class="row"><span class="small muted">طريقة التقسيم:</span><div class="seg"><button data-sm="oddeven" aria-pressed="true">فردي/زوجي</button><button data-sm="half" aria-pressed="false">النصف الأول/الثاني</button></div></div></div><div id="rs" class="stack"></div>
      <div class="grid2"><div class="card stack"><h3>سبيرمان–براون من معامل الارتباط</h3>${num('r', 'الارتباط بين النصفين r', 0.7)}${num('fac', 'معامل الإطالة (2 = اختبار كامل)', 2)}<div id="rsb"></div></div>
      <div class="card stack"><h3>جتمان من التباينات</h3><div class="grid3">${num('v1', 'σ² النصف الأول', 5)}${num('v2', 'σ² النصف الثاني', 4.5)}${num('vt', 'σ² الكلي', 18)}</div><div id="rg"></div></div></div>`;
    let sm = 'oddeven';
    const run = () => { const m = parseMatrix($('#mx', el).value), r = $('#rs', el); if (m.error) return out(r, `<div class="notice warn">${m.error}</div>`); if (m.k < 2) return out(r, '<div class="notice warn">تحتاج عبارتين على الأقل.</div>'); const s = PSY.splitHalf(m.rows, sm === 'half' ? 'half' : 'oddeven'); mark('split');
      out(r, `<div class="kpis"><div class="tile"><small>الارتباط بين النصفين</small><b>${fmt(s.r, 3)}</b></div><div class="tile hl"><small>سبيرمان–براون</small><b>${fmt(s.sb, 3)}</b></div><div class="tile hl"><small>جتمان</small><b>${fmt(s.guttman, 3)}</b></div></div>
        <p class="small muted">النصف الأول: العبارات ${s.A.map(x => x + 1).join('، ')}؛ النصف الثاني: ${s.B.map(x => x + 1).join('، ')}. ${Math.abs(s.v1 - s.v2) / Math.max(s.v1, s.v2) > 0.25 ? 'تباينا النصفين مختلفان بوضوح، فمعامل جتمان أنسب.' : 'تباينا النصفين متقاربان، فالمعاملان متقاربان.'}</p>
        <details class="steps"><summary>خطوات الحل</summary><div class="f">r_SB = 2 × ${fmt(s.r, 3)} ÷ (1 + ${fmt(s.r, 3)}) = ${fmt(s.sb, 3)}</div><div class="f">r_G = 2 × [1 − (${fmt(s.v1)} + ${fmt(s.v2)}) ÷ ${fmt(s.vt)}] = ${fmt(s.guttman, 3)}</div></details>`); };
    const runSB = () => { const r = v('r'), f = v('fac') || 2; const x = PSY.spearmanBrown(r, f); mark('split'); out($('#rsb', el), `<div class="tile hl"><small>الثبات المصحح</small><b>${fmt(x, 3)}</b></div><div class="f">${f} × ${r} ÷ (1 + ${fmt(f - 1, 2)} × ${r}) = ${fmt(x, 3)}</div>`); };
    const runG = () => { const a = v('v1'), b = v('v2'), t = v('vt'); const x = PSY.guttman(a, b, t); mark('split'); out($('#rg', el), `<div class="tile hl"><small>معامل جتمان</small><b>${fmt(x, 3)}</b></div><div class="f">2 × [1 − (${a} + ${b}) ÷ ${t}] = ${fmt(x, 3)}</div>`); };
    $$('[data-sm]', el).forEach(b => b.onclick = () => { sm = b.dataset.sm; $$('[data-sm]', el).forEach(x => x.setAttribute('aria-pressed', x === b)); run(); });
    $('#mx', el).oninput = run; $$('#r,#fac', el).forEach(i => i.oninput = runSB); $$('#v1,#v2,#vt', el).forEach(i => i.oninput = runG); run(); runSB(); runG();
  };
  /* ---------- retest ---------- */
  R.retest = el => {
    el.innerHTML = `<div class="card stack">${matrixInput('mx', '23 25\n18 20\n30 29\n25 27\n15 17\n28 26\n20 22\n33 32\n22 21\n27 29\n19 18\n24 26', 'كل مستجيب في سطر: درجة التطبيق الأول ثم الثاني')}</div><div id="rr" class="stack"></div>`;
    const run = () => {
      const m = parseMatrix($('#mx', el).value), r = $('#rr', el); if (m.error) return out(r, `<div class="notice warn">${m.error}</div>`); if (m.k !== 2) return out(r, '<div class="notice warn">أدخل عمودين فقط: التطبيق الأول والثاني.</div>');
      const x = m.rows.map(a => a[0]), y = m.rows.map(a => a[1]), rr = PSY.pearson(x, y), [c, t] = relVerdict(rr), n = x.length, tt = rr * Math.sqrt((n - 2) / (1 - rr * rr)), p = PSY.tP2(tt, n - 2); mark('retest');
      const W = 320, H = 240, P = 30, mn = Math.min(...x, ...y), mx = Math.max(...x, ...y), sc = vv => P + (W - 2 * P) * (vv - mn) / ((mx - mn) || 1);
      const pts = x.map((a, i) => `<circle cx="${sc(a)}" cy="${H - sc(y[i]) + 0}" r="4" fill="var(--primary)" opacity=".75"><title>${a} ، ${y[i]}</title></circle>`).join('');
      out(r, `<div class="kpis"><div class="tile hl"><small>معامل الثبات r</small><b>${fmt(rr, 3)}</b></div><div class="tile"><small>n</small><b>${n}</b></div><div class="tile"><small>متوسط التطبيق 1</small><b>${fmt(PSY.mean(x))}</b></div><div class="tile"><small>متوسط التطبيق 2</small><b>${fmt(PSY.mean(y))}</b></div><div class="tile"><small>الدلالة p</small><b>${fmtp(p)}</b></div></div><div class="verdict ${c}">الاستقرار ${t}</div>
        <div class="card chart" style="max-width:420px"><svg viewBox="0 0 ${W} ${H}" role="img" aria-label="رسم الانتشار"><line x1="${P}" y1="${H - P}" x2="${W - P}" y2="${H - P}" stroke="var(--line)"/><line x1="${P}" y1="${P}" x2="${P}" y2="${H - P}" stroke="var(--line)"/><line x1="${sc(mn)}" y1="${H - sc(mn)}" x2="${sc(mx)}" y2="${H - sc(mx)}" stroke="var(--amber)" stroke-dasharray="4 4"/>${pts}<text x="${W / 2}" y="${H - 6}" text-anchor="middle">التطبيق الأول</text><text x="10" y="${H / 2}" transform="rotate(-90 10 ${H / 2})" text-anchor="middle">التطبيق الثاني</text></svg></div>
        ${apa(`بلغ معامل ثبات إعادة التطبيق r(${n - 2}) = ${dec(rr)}، p ${fmtp(p).startsWith('<') ? fmtp(p) : '= ' + fmtp(p)}، مما يدل على استقرار ${t.split(' ')[0]} للدرجات عبر الزمن.`)}`); bindCopy(r);
    };
    $('#mx', el).oninput = run; run();
  };
  /* ---------- kappa ---------- */
  R.kappa = el => {
    let size = 2, T = [[20, 6], [4, 20]];
    const draw = () => {
      el.innerHTML = `<div class="card stack"><div class="row"><span class="small muted">عدد الفئات:</span>${seg('sz', [[2, '2'], [3, '3'], [4, '4']], size)}</div><p class="small muted" style="margin:0">الصفوف: تصنيفات المقدِّر الأول، والأعمدة: تصنيفات المقدِّر الثاني.</p><div class="tw"><table><thead><tr><th></th>${T.map((_, j) => `<th>م2: فئة ${j + 1}</th>`).join('')}</tr></thead><tbody>${T.map((r, i) => `<tr><th>م1: فئة ${i + 1}</th>${r.map((c, j) => `<td><input type="number" min="0" data-c="${i},${j}" value="${c}" aria-label="خلية ${i + 1}،${j + 1}" style="min-width:70px"></td>`).join('')}</tr>`).join('')}</tbody></table></div></div><div id="rk" class="stack"></div>`;
      $$('[data-sz]', el).forEach(b => b.onclick = () => { size = +b.dataset.sz; T = Array.from({ length: size }, (_, i) => Array.from({ length: size }, (_, j) => i === j ? 10 : 2)); draw(); });
      $$('[data-c]', el).forEach(i => i.oninput = () => { const [a, b] = i.dataset.c.split(',').map(Number); T[a][b] = Math.max(0, toNum(i.value) || 0); run(); });
      run();
    };
    const run = () => { const k = PSY.kappa(T), r = $('#rk', el); if (!k.N) return out(r, ''); mark('kappa');
      const lab = k.kappa < 0 ? 'أقل من الصدفة' : k.kappa <= 0.2 ? 'ضعيف جدًا' : k.kappa <= 0.4 ? 'ضعيف' : k.kappa <= 0.6 ? 'متوسط' : k.kappa <= 0.8 ? 'جيد' : 'شبه تام';
      out(r, `<div class="kpis"><div class="tile hl"><small>كابا κ</small><b>${fmt(k.kappa, 3)}</b></div><div class="tile"><small>الاتفاق الملاحظ Pₒ</small><b>${fmt(k.po, 3)}</b></div><div class="tile"><small>الاتفاق المتوقع Pₑ</small><b>${fmt(k.pe, 3)}</b></div><div class="tile"><small>N</small><b>${k.N}</b></div></div><div class="verdict ${k.kappa > 0.6 ? 'good' : k.kappa > 0.4 ? 'warn' : 'bad'}">اتفاق ${lab} (لانديس وكوخ، 1977)</div>
        <details class="steps"><summary>خطوات الحل</summary><div class="f">Pₒ = مجموع القطر ÷ N = ${fmt(k.po, 3)}</div><div class="f">Pₑ = Σ (مجموع الصف ÷ N) × (مجموع العمود ÷ N) = ${fmt(k.pe, 3)}</div><div class="f">κ = (${fmt(k.po, 3)} − ${fmt(k.pe, 3)}) ÷ (1 − ${fmt(k.pe, 3)}) = ${fmt(k.kappa, 3)}</div></details>
        ${apa(`بلغ معامل كابا كوهين لاتفاق المقدِّرين κ = ${dec(k.kappa)} (نسبة الاتفاق الملاحظ ${fmt(k.po * 100, 1)}%)، وهو اتفاق ${lab}.`)}`); bindCopy(r); };
    draw();
  };
  /* ---------- content validity ---------- */
  R.content = el => {
    let mode = 'cvr';
    const draw = () => {
      el.innerHTML = `<div class="card stack">${seg('cm', [['cvr', 'CVR (لوشي)'], ['cvi', 'CVI (لين)'], ['aiken', 'معامل V (أيكن)']], mode)}<div id="cin"></div></div><div id="cres" class="stack"></div>`;
      $$('[data-cm]', el).forEach(b => b.onclick = () => { mode = b.dataset.cm; draw(); });
      const ci = $('#cin', el), cr = $('#cres', el);
      if (mode === 'cvr') {
        ci.innerHTML = `<div class="grid3">${num('N', 'عدد المحكمين N', 10)}</div>${matrixInput('mx', '9 1\n10 0\n7 2\n6 3\n8 2\n5 4', 'كل عبارة في سطر: عدد من قالوا «أساسية» ثم (اختياريًا) عدد «مفيدة»')}<label class="row small"><input type="checkbox" id="mod"> استخدم الصيغة المعدلة (تحتسب «مفيدة» بنصف وزن)</label>`;
        const run = () => { const Nn = Math.round(v('N')), m = parseMatrix($('#mx', el).value); if (!(Nn >= 2)) return out(cr, '<div class="notice warn">أدخل عدد المحكمين.</div>'); if (m.error) return out(cr, `<div class="notice warn">${m.error}</div>`); const md = $('#mod', el).checked;
          const exact = PSY.cvrCritical(Nn), law = PSY.lawsheCritical(Nn); mark('content');
          const rows = m.rows.map((r, i) => { const x = md ? PSY.cvrModified(r[0], r[1] || 0, Nn) : PSY.cvr(r[0], Nn); const ok = x >= exact.cvr - 1e-9; return { cls: ok ? '' : 'flag', cells: ['ع' + (i + 1), N(r[0], 0), N(x), ok ? '<span class="chip good">تُقبل</span>' : '<span class="chip bad">تُراجع أو تُحذف</span>'] }; });
          const vals = m.rows.map(r => md ? PSY.cvrModified(r[0], r[1] || 0, Nn) : PSY.cvr(r[0], Nn)), keep = vals.filter(x => x >= exact.cvr - 1e-9);
          out(cr, `<div class="kpis"><div class="tile hl"><small>القيمة الحرجة (أير وسكالي)</small><b>${fmt(exact.cvr)}</b></div><div class="tile"><small>جدول لوشي الأصلي</small><b>${fmt(law)}</b></div><div class="tile"><small>CVI = متوسط CVR للمقبولة</small><b>${keep.length ? fmt(PSY.mean(keep), 3) : '—'}</b></div><div class="tile"><small>عبارات مقبولة</small><b>${keep.length}/${vals.length}</b></div></div>${tableHTML(['العبارة', 'أساسية', 'CVR', 'القرار'], rows)}
            <details class="steps"><summary>المعادلة والقيمة الحرجة</summary><div class="f">CVR = (nₑ − N/2) ÷ (N/2)${md ? '  (المعدلة: nₑ + nₘ/2)' : ''}</div><p class="small">القيمة الحرجة المحسوبة هنا من التوزيع ذي الحدين (اختبار من طرف واحد عند 0.05) كما أعاد حسابها أير وسكالي (2014): يلزم ${exact.ne} محكمين على الأقل من ${Nn} يقولون «أساسية». وجدول لوشي (1975) الأصلي يعطي ${fmt(law)} لهذا العدد.</p></details>`); };
        $$('#N,#mx,#mod', el).forEach(i => i.oninput = run); $('#mod', el).onchange = run; run();
      } else if (mode === 'cvi') {
        ci.innerHTML = matrixInput('mx', '4 4 3 4 4 3\n3 4 4 4 3 4\n2 3 2 3 4 2\n4 4 4 4 4 4\n3 3 4 2 3 3', 'كل عبارة في سطر، والأعمدة تقديرات المحكمين للصلة (1–4)');
        const run = () => { const m = parseMatrix($('#mx', el).value); if (m.error) return out(cr, `<div class="notice warn">${m.error}</div>`); const Rj = m.rows[0].map((_, j) => m.rows.map(r => r[j])); const c = PSY.cvi(Rj); mark('content');
          out(cr, `<div class="kpis"><div class="tile hl"><small>S‑CVI/Ave</small><b>${fmt(c.scviAve, 3)}</b></div><div class="tile"><small>S‑CVI/UA</small><b>${fmt(c.scviUA, 3)}</b></div><div class="tile"><small>المحكمون</small><b>${c.J}</b></div></div>
            ${tableHTML(['العبارة', 'I‑CVI', 'كابا المعدل k*', 'الحكم'], c.icvi.map((x, i) => { const ok = x >= (c.J <= 5 ? 1 : 0.78); return { cls: ok ? '' : 'flag', cells: ['ع' + (i + 1), N(x), N(c.kstar[i]), ok ? '<span class="chip good">ممتازة</span>' : '<span class="chip bad">تُراجع</span>'] }; }))}
            <div class="notice small">المحك: I‑CVI = 1 عند 5 محكمين فأقل، و≥ 0.78 عند 6 محكمين فأكثر (لين، 1986)، وS‑CVI/Ave ≥ 0.90 (بوليت وبيك، 2006). كابا المعدل يصحح احتمال الاتفاق بالصدفة.</div>`); };
        $('#mx', el).oninput = run; run();
      } else {
        ci.innerHTML = `<div class="grid3">${num('lo', 'أدنى تقدير', 1)}${num('c', 'عدد فئات التقدير c', 4)}</div>` + matrixInput('mx', '4 4 3 4 3\n3 3 4 4 4\n2 3 2 3 2', 'كل عبارة في سطر، والأعمدة تقديرات المحكمين');
        const run = () => { const lo = v('lo'), c = v('c'), m = parseMatrix($('#mx', el).value); if (m.error) return out(cr, `<div class="notice warn">${m.error}</div>`); mark('content');
          const rows = m.rows.map((r, i) => { const a = PSY.aikenV(r, lo, c); return { cls: a.V < 0.7 ? 'flag' : '', cells: ['ع' + (i + 1), N(a.S, 0), N(a.V, 3), a.V >= 0.8 ? '<span class="chip good">مرتفع</span>' : a.V >= 0.7 ? '<span class="chip warn">مقبول</span>' : '<span class="chip bad">منخفض</span>'] }; });
          out(cr, tableHTML(['العبارة', 'S = Σ(التقدير − الأدنى)', 'V', 'الحكم'], rows) + `<div class="f">V = S ÷ [n × (c − 1)]</div><p class="small muted">يُحكم على دلالة V بجداول أيكن (1985) بحسب عدد المحكمين والفئات؛ والمحك المتداول عمليًا 0.70–0.80.</p>`); };
        $$('#lo,#c,#mx', el).forEach(i => i.oninput = run); run();
      }
    };
    draw();
  };
  /* ---------- item analysis ---------- */
  R.items = el => {
    el.innerHTML = `<div class="card stack">${matrixInput('mx', EX_BIG)}<div class="row"><button class="btn sm ghost" data-ex="big">بيانات ليكرت محاكاة</button><button class="btn sm ghost" data-ex="bin">بيانات ثنائية</button><span class="small muted">نسبة المجموعتين الطرفيتين:</span><input type="number" id="frac" value="27" style="width:90px" aria-label="النسبة المئوية"></div></div><div id="ri" class="stack"></div>`;
    $$('[data-ex]', el).forEach(b => b.onclick = () => { $('#mx', el).value = b.dataset.ex === 'bin' ? EX_BIN : EX_BIG; run(); });
    const run = () => {
      const m = parseMatrix($('#mx', el).value), r = $('#ri', el); if (m.error) return out(r, `<div class="notice warn">${m.error}</div>`); if (m.n < 6) return out(r, '<div class="notice warn">تحتاج 6 مستجيبين على الأقل.</div>');
      const fr = Math.min(0.5, Math.max(0.1, v('frac') / 100)), bin = m.rows.every(x => x.every(y => y === 0 || y === 1)); mark('items');
      if (bin) {
        const d = PSY.difficultyDiscrimination(m.rows, fr);
        out(r, `<div class="notice small">بيانات ثنائية: الصعوبة p = نسبة الإجابة الصحيحة، والتمييز D = (ع − د) ÷ ن لكل مجموعة طرفية (ن = ${d.g}).</div>${tableHTML(['الفقرة', 'p (السهولة)', 'ع (العليا)', 'د (الدنيا)', 'D', 'الارتباط الثنائي المتسلسل النقطي', 'الحكم'], d.items.map((x, j) => { const dd = x.D; const lab = dd >= 0.4 ? ['good', 'تمييز جيد جدًا'] : dd >= 0.3 ? ['good', 'تمييز جيد'] : dd >= 0.2 ? ['warn', 'مقبول/يُحسَّن'] : ['bad', dd < 0 ? 'سالب: راجع المفتاح' : 'ضعيف']; return { cls: dd < 0.2 ? 'flag' : '', cells: ['ف' + (j + 1), N(x.p), N(x.U, 0), N(x.L, 0), N(x.D), N(x.rpb), `<span class="chip ${lab[0]}">${lab[1]}</span>`] }; }))}<p class="small muted">تصنيف إيبل للتمييز: ≥ 0.40 جيد جدًا، 0.30–0.39 جيد، 0.20–0.29 حدّي، أقل من 0.20 ضعيف.</p>`);
      } else {
        const s = PSY.itemStats(m.rows), e = PSY.extremeGroups(m.rows, fr);
        out(r, `<div class="notice small">بيانات متعددة الدرجات: مقارنة متوسطي المجموعتين الطرفيتين (ن = ${e.g} لكل مجموعة) باختبار ت لعينتين مستقلتين (ويلش)، مع الارتباط المصحح بالدرجة الكلية.</div>${tableHTML(['العبارة', 'المتوسط', 'متوسط العليا', 'متوسط الدنيا', 't', 'p', 'r المصحح', 'الحكم'], e.items.map((x, j) => { const ok = x.p < 0.05 && s.rit[j] >= 0.3; return { cls: ok ? '' : 'flag', cells: ['ع' + (j + 1), N(s.itemMean[j]), N(x.mHi), N(x.mLo), N(x.t), { n: 1, v: fmtp(x.p) }, N(s.rit[j]), ok ? '<span class="chip good">مميزة</span>' : '<span class="chip bad">تُراجع</span>'] }; }))}`);
      }
    };
    $$('#mx,#frac', el).forEach(i => i.oninput = run); run();
  };
  /* ---------- scores ---------- */
  R.scores = el => {
    el.innerHTML = `<div class="grid2"><div class="card stack"><h3>الدرجة المعيارية والخطأ المعياري</h3><div class="grid3">${num('m', 'المتوسط', 146.5)}${num('s', 'الانحراف المعياري', 15.95)}${num('rel', 'معامل الثبات', 0.86)}${num('x', 'درجة الفرد', 160)}</div><div id="r1"></div></div>
      <div class="card stack"><h3>المتوسط النظري ومستويات الاستجابة</h3><div class="grid3">${num('k', 'عدد العبارات', 42)}${num('lo', 'أدنى تقدير', 1)}${num('hi', 'أعلى تقدير', 5)}${num('nl', 'عدد المستويات', 3)}</div><div id="r2"></div></div></div>`;
    const r1 = () => { const m = v('m'), s = v('s'), rel = v('rel'), x = v('x'); if (!(s > 0)) return; const z = (x - m) / s, T = 50 + 10 * z, pr = PSY.ncdf(z) * 100, se = PSY.sem(s, rel); mark('scores');
      out($('#r1', el), `<div class="kpis"><div class="tile"><small>z</small><b>${fmt(z, 3)}</b></div><div class="tile hl"><small>T</small><b>${fmt(T, 1)}</b></div><div class="tile"><small>المئين (اعتدالي)</small><b>${fmt(pr, 1)}</b></div><div class="tile"><small>SEM</small><b>${fmt(se, 2)}</b></div></div><div class="f">فترة الثقة 95% = ${x} ± 1.96 × ${fmt(se, 2)} = [${fmt(x - 1.96 * se, 1)} ، ${fmt(x + 1.96 * se, 1)}]</div><p class="small muted">الدرجة الحقيقية للفرد تقع غالبًا في هذه الفترة؛ لذلك تُفسَّر الفروق الصغيرة حول درجات القطع بحذر.</p>`); };
    const r2 = () => { const k = v('k'), lo = v('lo'), hi = v('hi'), nl = Math.max(2, Math.round(v('nl'))); if (!(k > 0 && hi > lo)) return; const L = PSY.levels(lo, hi, nl), names = nl === 3 ? ['منخفض', 'متوسط', 'مرتفع'] : nl === 5 ? ['منخفض جدًا', 'منخفض', 'متوسط', 'مرتفع', 'مرتفع جدًا'] : L.ranges.map((_, i) => 'المستوى ' + (i + 1)); mark('scores');
      out($('#r2', el), `<div class="kpis"><div class="tile hl"><small>المتوسط النظري للمقياس</small><b>${fmt(k * (lo + hi) / 2, 1)}</b></div><div class="tile"><small>مدى الدرجة الكلية</small><b>${k * lo}–${k * hi}</b></div><div class="tile"><small>طول الفئة</small><b>${fmt(L.width, 2)}</b></div></div>${tableHTML(['المستوى', 'متوسط العبارة', 'الدرجة الكلية'], L.ranges.map((r, i) => [names[i], { n: 1, v: fmt(r[0], 2) + ' – ' + fmt(r[1], 2) }, { n: 1, v: fmt(r[0] * k, 1) + ' – ' + fmt(r[1] * k, 1) }]))}<p class="small muted">هذا تقسيم وصفي اصطلاحي، وليس معيارًا محكيًا ولا معياريًا؛ وللحكم المعياري استخدم الرتب المئينية لعينة التقنين.</p>`); };
    $$('#m,#s,#rel,#x', el).forEach(i => i.oninput = r1); $$('#k,#lo,#hi,#nl', el).forEach(i => i.oninput = r2); r1(); r2();
  };
  /* ---------- omega ---------- */
  R.omega = el => {
    el.innerHTML = `<div class="card stack"><div class="field"><label class="fl" for="ld">التشبعات المعيارية للعبارات (مفصولة بمسافة)</label><input type="text" id="ld" class="ltr" value="0.70 0.65 0.80 0.60 0.75"></div></div><div id="ro" class="stack"></div>`;
    const run = () => { const L = parseNums($('#ld', el).value).filter(x => Math.abs(x) <= 1), r = $('#ro', el); if (L.length < 2) return out(r, '<div class="notice warn">أدخل تشبعين على الأقل بين −1 و1.</div>'); const o = PSY.omega(L), [c, t] = relVerdict(o.omega); mark('omega');
      out(r, `<div class="kpis"><div class="tile hl"><small>أوميغا ω (= الثبات المركب CR)</small><b>${fmt(o.omega, 3)}</b></div><div class="tile ${o.ave >= 0.5 ? 'good' : 'bad'}"><small>متوسط التباين المستخلص AVE</small><b>${fmt(o.ave, 3)}</b></div><div class="tile"><small>√AVE</small><b>${fmt(Math.sqrt(o.ave), 3)}</b></div></div><div class="verdict ${c}">الثبات ${t}؛ ${o.ave >= 0.5 ? 'والصدق التقاربي متحقق (AVE ≥ 0.50)' : 'والصدق التقاربي غير متحقق (AVE < 0.50)'}</div>
        <details class="steps"><summary>خطوات الحل</summary><div class="f">Σλ = ${fmt(o.sl, 3)} ← (Σλ)² = ${fmt(o.sl2, 3)}</div><div class="f">Σθ = Σ(1 − λ²) = ${o.theta.map(x => fmt(x, 3)).join(' + ')} = ${fmt(o.st, 3)}</div><div class="f">ω = ${fmt(o.sl2, 3)} ÷ (${fmt(o.sl2, 3)} + ${fmt(o.st, 3)}) = ${fmt(o.omega, 3)}</div><div class="f">AVE = متوسط λ² = ${fmt(o.ave, 3)}</div><p class="small muted">للصدق التمييزي (فورنل ولاركر، 1981): يجب أن يزيد √AVE لكل بعد على ارتباطه بالأبعاد الأخرى.</p></details>`); };
    $('#ld', el).oninput = run; run();
  };
  /* ---------- DIF ---------- */
  R.dif = el => {
    el.innerHTML = `<div class="card stack"><p class="small muted" style="margin:0">قسّم المستجيبين إلى طبقات بحسب الدرجة الكلية، واكتب لكل طبقة في سطر: a = المرجعية صحيح/موافق، b = المرجعية خطأ، c = المستهدفة صحيح، d = المستهدفة خطأ.</p>${matrixInput('mx', '12 8 8 12\n20 10 14 16\n25 5 20 10\n18 2 16 4', 'a b c d لكل طبقة')}</div><div id="rd" class="stack"></div>`;
    const run = () => { const m = parseMatrix($('#mx', el).value), r = $('#rd', el); if (m.error) return out(r, `<div class="notice warn">${m.error}</div>`); if (m.k !== 4) return out(r, '<div class="notice warn">أدخل أربعة أعمدة a b c d لكل طبقة.</div>'); const x = PSY.mantelHaenszel(m.rows.map(q => ({ a: q[0], b: q[1], c: q[2], d: q[3] }))); mark('dif');
      const lab = { A: ['good', 'A: أداء تفاضلي ضئيل أو غير دال'], B: ['warn', 'B: أداء تفاضلي متوسط'], C: ['bad', 'C: أداء تفاضلي كبير'] }[x.cls];
      out(r, `<div class="kpis"><div class="tile"><small>α_MH</small><b>${fmt(x.alphaMH, 3)}</b></div><div class="tile hl"><small>Δ_MH</small><b>${fmt(x.delta, 2)}</b></div><div class="tile"><small>χ²_MH</small><b>${fmt(x.chi, 3)}</b></div><div class="tile"><small>p</small><b>${fmtp(x.p)}</b></div></div><div class="verdict ${lab[0]}">${lab[1]}</div>
        <details class="steps"><summary>المعادلات</summary><div class="f">α_MH = Σ(a·d/N) ÷ Σ(b·c/N)</div><div class="f">Δ_MH = −2.35 × ln(α_MH)</div><div class="f">χ²_MH = (|Σa − ΣE(a)| − 0.5)² ÷ ΣVar(a)</div><p class="small muted">Δ سالبة تعني أن المفردة أصعب على المجموعة المستهدفة. التصنيف: A إذا |Δ| < 1 أو غير دال، B إذا 1 ≤ |Δ| < 1.5، C إذا |Δ| ≥ 1.5 (هيئة الاختبارات التربوية؛ Zwick, 2012، بصورة مبسطة). والإشارة الإحصائية لا تكفي للحكم بالتحيز: تُراجع صياغة المفردة ومحتواها.</p></details>`); };
    $('#mx', el).oninput = run; run();
  };
  /* ---------- EFA ---------- */
  R.efa = el => {
    el.innerHTML = `<div class="card stack">${matrixInput('mx', EX_BIG)}<div class="row"><span class="small muted">عدد العوامل:</span><input type="number" id="nf" value="" placeholder="تلقائي" style="width:110px" aria-label="عدد العوامل"><button class="btn sm" id="go">${ic('play')} حلّل</button></div></div><div id="re" class="stack"></div>`;
    const run = () => {
      const m = parseMatrix($('#mx', el).value), r = $('#re', el); if (m.error) return out(r, `<div class="notice warn">${m.error}</div>`);
      if (m.k < 3 || m.n < m.k + 2) return out(r, '<div class="notice warn">تحتاج 3 عبارات على الأقل، وعدد مستجيبين أكبر من عدد العبارات.</div>');
      r.innerHTML = '<p class="muted">جارٍ الحساب…</p>';
      setTimeout(() => {
        const R0 = PSY.corrMatrix(m.rows), km = PSY.kmo(R0), bt = PSY.bartlett(R0, m.n), eig = PSY.jacobiEigen(R0).values, pa = PSY.parallelAnalysis(m.n, m.k, 30);
        const kaiser = eig.filter(x => x > 1).length, par = eig.filter((x, i) => x > pa.p95[i]).length || 1, nfIn = Math.round(v('nf')), nf = nfIn > 0 ? Math.min(nfIn, m.k - 1) : par;
        const fa = PSY.efa(m.rows, nf), L = fa.rotated; mark('efa');
        const kv = km ? km.kmo : NaN, kl = kv >= 0.9 ? 'رائع' : kv >= 0.8 ? 'جدير بالتقدير' : kv >= 0.7 ? 'متوسط' : kv >= 0.6 ? 'مقبول' : 'غير مقبول';
        out(r, `<div class="kpis"><div class="tile ${kv >= 0.6 ? 'good' : 'bad'}"><small>KMO</small><b>${fmt(kv, 3)}</b></div><div class="tile ${bt.p < 0.05 ? 'good' : 'bad'}"><small>بارتلت χ² (df = ${bt.df})</small><b>${fmt(bt.chi, 1)}</b></div><div class="tile"><small>دلالة بارتلت</small><b>${fmtp(bt.p)}</b></div><div class="tile"><small>عوامل كايزر (جذر > 1)</small><b>${kaiser}</b></div><div class="tile hl"><small>عوامل التحليل الموازي</small><b>${par}</b></div></div>
          <div class="verdict ${kv >= 0.6 && bt.p < 0.05 ? 'good' : 'bad'}">ملاءمة البيانات للتحليل العاملي: KMO ${kl}${bt.p < 0.05 ? '، واختبار بارتلت دال' : '، واختبار بارتلت غير دال'}</div>
          <section class="card"><h3>${ic('chart')} منحنى الجذور الكامنة مع التحليل الموازي</h3>${barSVG(eig, eig.map((_, i) => String(i + 1)), { line: pa.p95, ref: 1, refLabel: 'جذر = 1', label: 'الجذور الكامنة', color: (x, i) => x > pa.p95[i] ? 'var(--primary)' : 'var(--line)' })}<p class="small muted">الأعمدة: الجذور الكامنة للبيانات، والخط المتقطع: المئين 95 لجذور بيانات عشوائية بالحجم نفسه (هورن، 1965). يُحتفظ بالعوامل التي تتجاوز جذورها الخط.</p>
          ${tableHTML(['العامل', 'الجذر الكامن', '% التباين', 'المئين 95 العشوائي'], eig.slice(0, Math.min(8, eig.length)).map((x, i) => [String(i + 1), N(x, 3), N(fa.varExp[i], 1), N(pa.p95[i], 3)]))}</section>
          <section class="card"><h3>تشبعات المكونات الأساسية بعد تدوير فاريماكس (${nf} ${nf === 1 ? 'عامل' : 'عوامل'})</h3>${tableHTML(['العبارة', ...L[0].map((_, f) => 'ع' + (f + 1)), 'الشيوع h²', 'MSA'], L.map((row, i) => { const mx = Math.max(...row.map(Math.abs)), cross = row.filter(x => Math.abs(x) >= 0.3).length > 1; return { cls: mx < 0.4 ? 'flag' : '', cells: [(m.head && m.head[i]) ? esc(m.head[i]) : 'ع' + (i + 1), ...row.map(x => ({ n: 1, v: Math.abs(x) === mx && mx >= 0.4 ? `<b>${fmt(x)}</b>` : fmt(x) })), N(PSY.sum(row.map(x => x * x))), N(km ? km.msa[i] : NaN)].concat([]) }; }))}<p class="small muted">العبارة المظللة تشبعها الأعلى أقل من 0.40. والتشبعات المتقاطعة (≥ 0.30 على أكثر من عامل) تستدعي مراجعة الصياغة. هذا تحليل مكونات أساسية تقريبي للاستكشاف؛ وللنشر يُفضَّل تحليل المحاور الأساسية أو الأرجحية العظمى ثم التحليل العاملي التوكيدي على عينة مستقلة.</p></section>`);
      }, 20);
    };
    $('#go', el).onclick = run; run();
  };
  /* ---------- sample size ---------- */
  R.sample = el => {
    el.innerHTML = `<div class="card stack">${num('k', 'عدد عبارات المقياس', 30)}</div><div id="rsz"></div>`;
    const run = () => { const k = Math.max(1, Math.round(v('k'))); mark('sample');
      out($('#rsz', el), tableHTML(['المرحلة', 'الحجم المقترح', 'المرجع'], [
        ['المقابلات المعرفية', '5–15 فردًا من الفئة المستهدفة', 'Willis (2005)'],
        ['التجريب الاستطلاعي', '30 فردًا تقريبًا', 'متداول'],
        ['تحليل العبارات والثبات', `${Math.max(100, 5 * k)} فأكثر (5 أفراد لكل عبارة)`, 'Nunnally & Bernstein (1994)'],
        ['التحليل العاملي الاستكشافي', `${Math.max(200, 10 * k)} (10 لكل عبارة) وحد أدنى 200؛ 300 جيد، 500 جيد جدًا`, 'Comrey & Lee (1992)؛ Costello & Osborne (2005)'],
        ['التحليل العاملي التوكيدي', '200 فأكثر في عينة مستقلة عن الاستكشافية', 'Boateng et al. (2018)'],
        ['فحص الأداء التفاضلي', '200 فأكثر في كل مجموعة مقارنة', 'متداول في دراسات مانتل–هانزل'],
        ['إعادة التطبيق', '50 فأكثر بفاصل نحو أسبوعين', 'متداول']]) + '<p class="small muted">هذه قواعد إرشادية؛ والحجم المطلوب فعليًا يتوقف على قوة التشبعات وعدد العبارات لكل عامل ودرجة الاشتراكيات.</p>'); };
    $('#k', el).oninput = run; run();
  };

  /* ---------- formulas page ---------- */
  const FX = [
    ['الثبات', 'ألفا كرونباخ', "Cronbach's α", `α = ${frac('k', 'k − 1')} × [1 − ${frac('Σσ²ᵢ', 'σ²ₓ')}]`, [['k', 'عدد العبارات'], ['σ²ᵢ', 'تباين العبارة'], ['σ²ₓ', 'تباين الدرجة الكلية']], 'الاتساق الداخلي لفقرات متعددة الدرجات', 'alpha'],
    ['الثبات', 'كودر–ريتشاردسون 20', 'KR‑20', `KR‑20 = ${frac('k', 'k − 1')} × [1 − ${frac('Σpq', 'σ²ₓ')}]`, [['p', 'نسبة الصواب في الفقرة'], ['q', '1 − p']], 'الفقرات ثنائية الدرجة', 'kr'],
    ['الثبات', 'كودر–ريتشاردسون 21', 'KR‑21', `KR‑21 = ${frac('k', 'k − 1')} × [1 − ${frac('M(k − M)', 'k σ²ₓ')}]`, [['M', 'متوسط الدرجات']], 'فقرات ثنائية متقاربة الصعوبة', 'kr'],
    ['الثبات', 'سبيرمان–براون', 'Spearman–Brown', `r_SB = ${frac('2r', '1 + r')}`, [['r', 'الارتباط بين النصفين']], 'تصحيح التجزئة النصفية؛ وصيغته العامة n·r ÷ [1 + (n − 1)r] لإطالة المقياس n مرة', 'split'],
    ['الثبات', 'جتمان للتجزئة النصفية', 'Guttman split‑half', `r = 2 × [1 − ${frac('σ²₁ + σ²₂', 'σ²ₓ')}]`, [['σ²₁، σ²₂', 'تباينا النصفين']], 'عند عدم تساوي تبايني النصفين', 'split'],
    ['الثبات', 'إعادة الاختبار', 'Test–retest', `r = ${frac('Cov(X₁, X₂)', 'σ₁ σ₂')}`, [['X₁، X₂', 'درجات التطبيقين']], 'استقرار الدرجات عبر الزمن', 'retest'],
    ['الثبات', 'كابا كوهين', "Cohen's κ", `κ = ${frac('Pₒ − Pₑ', '1 − Pₑ')}`, [['Pₒ', 'الاتفاق الملاحظ'], ['Pₑ', 'الاتفاق المتوقع بالصدفة']], 'اتفاق مقدِّرين على تصنيف اسمي', 'kappa'],
    ['الثبات', 'أوميغا ماكدونالد', "McDonald's ω", `ω = ${frac('(Σλ)²', '(Σλ)² + Σθ')}`, [['λ', 'التشبع المعياري'], ['θ', '1 − λ²']], 'الاتساق الداخلي مع تشبعات غير متساوية', 'omega'],
    ['الثبات', 'الخطأ المعياري للقياس', 'SEM', 'SEM = σ √(1 − r)', [['σ', 'الانحراف المعياري'], ['r', 'معامل الثبات']], 'فترة الثقة للدرجة الفردية: X ± 1.96 SEM', 'scores'],
    ['الصدق', 'نسبة صدق المحتوى', 'Lawshe CVR', `CVR = ${frac('nₑ − N/2', 'N/2')}`, [['nₑ', 'عدد من قالوا «أساسية»'], ['N', 'عدد المحكمين']], 'قبول العبارة إذا تجاوزت القيمة الحرجة لعدد المحكمين', 'content'],
    ['الصدق', 'مؤشر صدق المحتوى للعبارة', 'I‑CVI', `I‑CVI = ${frac('عدد التقديرات 3 أو 4', 'عدد المحكمين')}`, [['S‑CVI/Ave', 'متوسط I‑CVI']], 'I‑CVI ≥ 0.78 مع 6 محكمين فأكثر', 'content'],
    ['الصدق', 'معامل V لأيكن', "Aiken's V", `V = ${frac('Σ(r − lo)', 'n(c − 1)')}`, [['r', 'تقدير المحكم'], ['lo', 'أدنى تقدير'], ['c', 'عدد الفئات']], 'اتفاق المحكمين على تدريج متعدد الدرجات', 'content'],
    ['الصدق', 'متوسط التباين المستخلص', 'AVE', `AVE = ${frac('Σλ²', 'k')}`, [['λ', 'التشبع المعياري']], 'الصدق التقاربي: AVE ≥ 0.50', 'omega'],
    ['تحليل الفقرات', 'معامل الصعوبة (السهولة)', 'Item difficulty', `p = ${frac('عدد الإجابات الصحيحة', 'عدد المجيبين')}`, [], 'كلما ارتفع p كانت الفقرة أسهل', 'items'],
    ['تحليل الفقرات', 'معامل التمييز', 'Discrimination D', `D = ${frac('ع − د', 'ن')}`, [['ع، د', 'الإجابات الصحيحة في المجموعتين العليا والدنيا'], ['ن', 'حجم كل مجموعة (27%)']], 'من −1 إلى +1؛ ≥ 0.30 جيد', 'items'],
    ['الدرجات', 'الدرجة الزائية والتائية', 'z, T', `z = ${frac('X − M', 'σ')}  ؛  T = 50 + 10z`, [], 'تحويل الدرجات الخام إلى معايير', 'scores'],
    ['الدرجات', 'طول الفئة للمستويات', 'Class width', `w = ${frac('أعلى تقدير − أدنى تقدير', 'عدد المستويات')}`, [], 'تقسيم وصفي للمتوسطات (مثل 1.33 لخماسي بثلاثة مستويات)', 'scores'],
    ['العدالة', 'دلتا مانتل–هانزل', 'MH D‑DIF', 'Δ_MH = −2.35 ln(α_MH)', [['α_MH', 'نسبة الأرجحية المشتركة']], 'تصنيف A/B/C للأداء التفاضلي', 'dif'],
    ['البنية العاملية', 'اختبار بارتلت للكروية', "Bartlett's test", 'χ² = −[n − 1 − (2p + 5)/6] ln|R|', [['p', 'عدد العبارات'], ['|R|', 'محدد مصفوفة الارتباط']], 'يجب أن يكون دالًا لملاءمة التحليل العاملي', 'efa']
  ];
  VIEWS.formulas = el => {
    const cats = [...new Set(FX.map(f => f[0]))];
    el.innerHTML = `<div class="section-title"><h2>المعادلات</h2><span class="muted small">كل معادلة مع رموزها ومتى تُستخدم، وزر لتجربتها في المعمل</span></div>` +
      cats.map(c => `<section class="stack"><h3 style="font-size:18px">${esc(c)}</h3><div class="grid2">${FX.filter(f => f[0] === c).map(f => `<article class="card formula-card"><div class="row"><b class="grow" style="font-family:var(--f-display)">${esc(f[1])}</b><span class="chip ltr">${esc(f[2])}</span></div><div class="fbox">${f[3]}</div>${f[4].length ? `<dl class="sym">${f[4].map(s => `<dt>${s[0]}</dt><dd>${esc(s[1])}</dd>`).join('')}</dl>` : ''}<p class="small muted" style="margin:0">${esc(f[5])}</p><div><button class="btn sm soft" data-try="${f[6]}">${ic('flask')} جرّب في المعمل</button></div></article>`).join('')}</div></section>`).join('') +
      `<section class="card"><h3>${ic('shield')} تفسير قيم الثبات</h3>${tableHTML(['القيمة', 'التفسير'], [['≥ 0.90', 'ممتاز، مناسب للقرارات الفردية'], ['0.80 – 0.89', 'جيد'], ['0.70 – 0.79', 'مقبول للأغراض البحثية'], ['< 0.70', 'ضعيف؛ يُراجع المقياس']])}<p class="small muted">محكات اصطلاحية؛ ولكابا: 0.41–0.60 متوسط، 0.61–0.80 جيد، > 0.80 شبه تام (لانديس وكوخ، 1977).</p></section>`;
    $$('[data-try]', el).forEach(b => b.onclick = () => open(b.dataset.try));
  };

  return { list, open, clear: () => { openId = null; } };
})();

/* ===== Audience adaptation: age group / generation × language register, without breaking measurement =====
   Shared by the scales lab and the statistics lab. Rules:
   - wording may change; the construct, the item's direction (positive/negative), the number of response points and the
     time frame may not — otherwise the versions measure different things and cannot be compared;
   - every adapted version is checked automatically here, then empirically after data collection (DIF / invariance).
   Generation labels are birth cohorts (approximate); what drives wording is reading level and life context. */
var AUD = (function () {
  var PROFILES = [
    { id: 'alpha', t: 'الأطفال 8–12 سنة', g: 'جيل ألفا (مواليد 2013 وما بعدها)', maxW: 9, neg: false, pts: [3, 4], ctx: 'المدرسة، والبيت، والأصدقاء، والألعاب', tips: ['جمل قصيرة جدًا (أقل من 9 كلمات) وكلمات مألوفة', 'تجنب النفي والعبارات المعكوسة؛ الأطفال يخطئون في فهمها', 'تدريج من 3 نقاط بكلمات بسيطة أو وجوه', 'زمن قريب ومحدد: «هذا الأسبوع» بدل «عادةً»', 'اقرأ العبارات بصوت مرتفع للصغار'] },
    { id: 'teen', t: 'المراهقون 13–17 سنة', g: 'جيل Z (مواليد 1997–2012)', maxW: 12, neg: true, pts: [4, 5], ctx: 'المدرسة، والأصدقاء، ووسائل التواصل، والأسرة', tips: ['جمل قصيرة (حتى 12 كلمة)', 'عبارات سالبة قليلة وواضحة، ولا نفي مزدوج', 'لغة غير وعظية ولا تبدو كأسئلة الكبار', 'أكّد السرية: المراهقون أكثر حساسية للمرغوبية الاجتماعية'] },
    { id: 'youth', t: 'الشباب 18–28 سنة', g: 'جيل Z (مواليد 1997–2012)', maxW: 15, neg: true, pts: [5, 7], ctx: 'الجامعة، والعمل الأول، والعلاقات، والحياة الرقمية', tips: ['جمل مباشرة (حتى 15 كلمة)', 'التدريج الخماسي مناسب', 'الأمثلة الرقمية مقبولة في التعليمات لا في العبارات'] },
    { id: 'millennial', t: 'البالغون 29–44 سنة', g: 'جيل الألفية Y (مواليد 1981–1996)', maxW: 18, neg: true, pts: [5, 7], ctx: 'العمل، والأسرة، والأبوة والأمومة', tips: ['الصياغة القياسية مناسبة', 'راعِ ضيق الوقت: مقياس أقصر يرفع جودة الاستجابة'] },
    { id: 'genx', t: 'منتصف العمر 45–60 سنة', g: 'الجيل X (مواليد 1965–1980)', maxW: 16, neg: true, pts: [4, 5], ctx: 'العمل، والأسرة، والصحة', tips: ['فصحى ميسرة غالبًا أوضح من العامية في السياق الرسمي', 'تجنب المصطلحات الرقمية الحديثة غير الشائعة', 'خط أكبر وتباعد مريح'] },
    { id: 'older', t: 'كبار السن 60+ سنة', g: 'جيل الطفرة وما قبله (مواليد قبل 1965)', maxW: 12, neg: false, pts: [3, 4], ctx: 'الصحة، والأسرة، والأنشطة اليومية', tips: ['جمل قصيرة وخط كبير', 'تجنب العبارات المعكوسة؛ تزيد أخطاء الاستجابة', 'تدريج من 3 أو 4 نقاط بكلمات واضحة', 'يمكن التطبيق بالمقابلة بدل الورقة'] }
  ];
  var REGISTERS = [
    { id: 'fusha', t: 'فصحى رسمية', d: 'للمقاييس المقننة والبحوث المنشورة' },
    { id: 'easy', t: 'فصحى ميسرة', d: 'فصحى بكلمات شائعة وتراكيب قصيرة؛ الأنسب لمعظم الفئات' },
    { id: 'ammiya', t: 'عامية مصرية', d: 'أقرب للغة الحياة اليومية؛ مناسبة للأطفال والمراهقين والاستبيانات غير الرسمية' }
  ];
  // phrase-level dictionary (longest first). Output always needs human review.
  var EASY = [['يصعب عليّ', 'من الصعب عليّ'], ['في كثير من الأحيان', 'كثيرًا'], ['لديّ', 'عندي'], ['لدي', 'عندي'], ['المحيطون بي', 'من حولي'], ['أدائي الدراسي', 'دراستي'], ['المهام', 'الأعمال المطلوبة'], ['أتجنب', 'أبتعد عن'], ['أستطيع', 'أقدر'], ['عندما', 'لما']];
  var AMM = [
    ['لا أستطيع', 'مش قادر'], ['ليست لديّ', 'معنديش'], ['ليست لدي', 'معنديش'], ['في كثير من', 'في كتير من'], ['المواقف اليومية', 'مواقف كل يوم'], ['لا أقدر', 'مش قادر'], ['ليس لديّ', 'معنديش'], ['ليس لدي', 'معنديش'], ['في كثير من الأحيان', 'كتير'], ['في المواقف اليومية', 'في يومي'],
    ['يصعب عليّ', 'صعب عليا'], ['أجد صعوبة في', 'بلاقي صعوبة في'], ['يلاحظ المحيطون بي', 'اللي حواليا بياخدوا بالهم من'], ['المحيطون بي', 'اللي حواليا'],
    ['أشعر بالضيق', 'بتضايق'], ['أشعر بأنني', 'بحس إني'], ['أشعر أنني', 'بحس إني'], ['أشعر بـ', 'بحس بـ'], ['أشعر', 'بحس'], ['أستطيع', 'أقدر'],
    ['أفكر', 'بفكر'], ['أحب', 'بحب'], ['أتصرف', 'بتصرف'], ['أتذكر', 'فاكر'], ['أؤجل', 'بأجّل'], ['أنشغل', 'بتلهي'], ['أتجنب', 'ببعد عن'], ['أقلق', 'بقلق'], ['أخاف', 'بخاف'],
    ['أنجز', 'بخلّص'], ['أنام', 'بنام'], ['أستيقظ', 'بصحى'], ['أترك', 'بسيب'], ['أبدأ', 'ببدأ'], ['أضع', 'بحط'], ['ألتزم', 'بلتزم'], ['أواجه', 'بواجه'], ['أهرب', 'بهرب'], ['أفاجأ', 'بتفاجئ'],
    ['أجد نفسي', 'بلاقي نفسي'], ['أتمنى', 'نفسي'], ['أريد', 'عايز'], ['يظهر', 'بيبان'], ['يؤثر', 'بيأثر'], ['يستمر', 'بيفضل'], ['تتكرر', 'بتتكرر'],
    ['لديّ', 'عندي'], ['لدي', 'عندي'], ['عندما', 'لما'], ['الذي', 'اللي'], ['التي', 'اللي'], ['الذين', 'اللي'], ['هذا', 'ده'], ['هذه', 'دي'], ['هؤلاء', 'دول'],
    ['كثيرًا', 'كتير'], ['كثيرا', 'كتير'], ['أحيانًا', 'ساعات'], ['أحيانا', 'ساعات'], ['دائمًا', 'دايمًا'], ['جدًا', 'أوي'], ['الآن', 'دلوقتي'], ['كيف', 'إزاي'], ['لماذا', 'ليه'], ['ماذا', 'إيه'],
    ['ليست', 'مش'], ['ليس', 'مش'], ['لست', 'مش'], ['أكثر من غيري', 'أكتر من غيري'], ['أكثر', 'أكتر'], ['أمام', 'قدام'], ['مع ذلك', 'ومع كده'], ['فقط', 'بس'], ['أيضًا', 'كمان']
  ];
  var LABELS = {
    fusha: { 3: ['تنطبق كثيرًا', 'تنطبق أحيانًا', 'لا تنطبق'], 4: ['تنطبق دائمًا', 'تنطبق غالبًا', 'تنطبق نادرًا', 'لا تنطبق أبدًا'], 5: ['تنطبق دائمًا', 'تنطبق غالبًا', 'تنطبق أحيانًا', 'تنطبق نادرًا', 'لا تنطبق أبدًا'] },
    easy: { 3: ['كثيرًا', 'أحيانًا', 'لا'], 4: ['دائمًا', 'غالبًا', 'قليلًا', 'أبدًا'], 5: ['دائمًا', 'غالبًا', 'أحيانًا', 'نادرًا', 'أبدًا'] },
    ammiya: { 3: ['كتير', 'ساعات', 'لأ'], 4: ['دايمًا', 'غالبًا', 'قليل', 'خالص'], 5: ['دايمًا', 'غالبًا', 'ساعات', 'نادرًا', 'خالص'] }
  };
  var AGREE = { fusha: ['أوافق بشدة', 'أوافق', 'محايد', 'أرفض', 'أرفض بشدة'], easy: ['موافق جدًا', 'موافق', 'مش متأكد', 'غير موافق', 'غير موافق أبدًا'], ammiya: ['موافق جدًا', 'موافق', 'مش متأكد', 'مش موافق', 'مش موافق خالص'] };
  var NEG = /(^|\s)(لا|لم|لن|ليس|ليست|لست|غير|عدم|بدون|دون|مش|ما)(\s|$)/;
  function prof(id) { return PROFILES.filter(function (p) { return p.id === id; })[0] || PROFILES[2]; }
  function esc(s) { return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); }
  // whole-word replacement; a first-person verb right after «لم/لن» is left for human review (its negated form differs)
  function swap(text, dict, notes) {
    var t = ' ' + text + ' ';
    dict.slice().sort(function (a, b) { return b[0].length - a[0].length; }).forEach(function (p) {
      var tail = /بـ$/.test(p[0]) ? '' : '(?=$|[\\s.،؟!])';
      t = t.replace(new RegExp('(^|\\s)(لم\\s|لن\\s)?' + esc(p[0]) + tail, 'g'), function (m, sp, neg) {
        if (neg && /^أ/.test(p[0])) { notes && notes.push('negverb'); return m; }
        return sp + (neg || '') + p[1];
      });
    });
    return t.trim();
  }
  function adapt(text, profileId, registerId) {
    var t = String(text || '').replace(/[«»"]/g, ''), notes = [];
    if (registerId === 'ammiya') t = swap(t, AMM, notes); else if (registerId === 'easy') t = swap(t, EASY, notes);
    t = t.replace(/بـ(?=[^\s])/g, 'ب').replace(/\s+/g, ' ').replace(/\s+([.،؟])/g, '$1').trim();
    var w = check(text, t, profileId, registerId);
    if (notes.length) w.push(['warn', 'فيها فعل منفي بـ«لم/لن»؛ صُغه بالعامية يدويًا مع الحفاظ على النفي']);
    if (registerId === 'ammiya' && /(لم|لن|سوف|إنني|لكنني)\s/.test(t)) w.push(['warn', 'ما زالت فيها تراكيب فصحى؛ راجعها']);
    return { text: t, warnings: w };
  }
  // measurement-safety checks between the original and an adapted item
  function check(orig, adapted, profileId, registerId) {
    var w = [], p = prof(profileId), wc = String(adapted).split(/\s+/).filter(Boolean).length, wo = String(orig).split(/\s+/).filter(Boolean).length;
    if (!String(adapted).trim()) return [['bad', 'عبارة فارغة']];
    if (wc > p.maxW) w.push(['warn', 'أطول من المناسب لهذه الفئة (' + wc + ' كلمة؛ الحد ' + p.maxW + ')']);
    if (wo && wc > wo * 1.6 + 2) w.push(['warn', 'أطول كثيرًا من الأصل؛ تأكد أنك لم تضف معنى جديدًا']);
    var negO = NEG.test(' ' + orig + ' '), negA = NEG.test(' ' + adapted + ' ');
    if (negO !== negA) w.push(['bad', 'تغيّر اتجاه العبارة (نفي/إثبات)؛ هذا يغيّر التصحيح ويكسر المقارنة']);
    if (negA && !p.neg) w.push(['warn', 'عبارة منفية لفئة يصعب عليها النفي؛ اكتب بديلًا مثبتًا وعدّل اتجاه التصحيح في كل النسخ']);
    if ((String(adapted).match(/ و/g) || []).length >= 2) w.push(['warn', 'قد تحمل أكثر من فكرة']);
    if (/[؟?]/.test(adapted)) w.push(['warn', 'صيغة سؤال؛ اجعلها عبارة تقريرية']);
    return w;
  }
  function labels(base, profileId, registerId) {
    var n = (base || []).length, isAgree = (base || []).some(function (x) { return /أوافق|موافق|أرفض/.test(x); });
    if (n === 5 && isAgree) return AGREE[registerId] ? AGREE[registerId].slice() : base;
    var L = (LABELS[registerId] || LABELS.fusha)[n];
    return L ? L.slice() : base;
  }
  function instructions(profileId, registerId) {
    var kid = profileId === 'alpha';
    if (registerId === 'ammiya') return kid ? 'هنقرا جمل قصيرة. اختار الإجابة اللي شبهك. مفيش إجابة صح وإجابة غلط.' : 'اقرا كل جملة كويس، واختار الإجابة اللي بتعبّر عنك فعلًا. مفيش إجابة صح وغلط، وإجاباتك سرية.';
    if (registerId === 'easy') return kid ? 'اقرأ كل جملة، واختر الإجابة التي تشبهك. لا توجد إجابة صحيحة وأخرى خاطئة.' : 'اقرأ كل عبارة، واختر الإجابة التي تعبّر عنك. لا توجد إجابات صحيحة أو خاطئة، وإجاباتك سرية.';
    return null;
  }
  function prompt(items, profileId, registerId, construct) {
    var p = prof(profileId), r = REGISTERS.filter(function (x) { return x.id === registerId; })[0] || REGISTERS[0];
    return 'أنت خبير في القياس النفسي وتكييف المقاييس. أعد صياغة عبارات مقياس «' + (construct || 'السمة') + '» لتناسب ' + p.t + ' (' + p.g + ') بلغة: ' + r.t + '.\n' +
      'قواعد إلزامية حتى لا يختل القياس:\n1) احتفظ بالمعنى والسمة المقاسة نفسها؛ لا تضف ولا تحذف فكرة.\n2) احتفظ باتجاه كل عبارة (المثبتة تبقى مثبتة والمنفية تبقى منفية).\n3) احتفظ بالإطار الزمني وشدة العبارة (لا تحوّل «أحيانًا» إلى «دائمًا»).\n4) فكرة واحدة لكل عبارة، وبحد أقصى ' + p.maxW + ' كلمة.\n5) سياقات مألوفة للفئة: ' + p.ctx + '، دون تغيير ما تقيسه العبارة.\n' +
      'أعد JSON فقط: {"items":[{"id":"...","text":"..."}]}\nالعبارات:\n' + items.map(function (x) { return JSON.stringify({ id: x.id, text: x.text }); }).join('\n');
  }
  /* empirical check after data: standardization DIF (Dorans & Kulick) of each item between two versions,
     conditioning on the rest score; returns per-item STD on the item's own scale. */
  function stdDIF(rowsRef, rowsFoc, nItems, lo, hi) {
    var first = stdPass(rowsRef, rowsFoc, nItems, lo, hi, {});
    var drop = {}; first.forEach(function (x, j) { if (x.level === 'flag' || x.level === 'review') drop[j] = 1; });
    return Object.keys(drop).length && Object.keys(drop).length < nItems - 1 ? stdPass(rowsRef, rowsFoc, nItems, lo, hi, drop) : first;
  }
  function stdPass(rowsRef, rowsFoc, nItems, lo, hi, drop) {
    var out = [], R = hi - lo || 1;
    for (var j = 0; j < nItems; j++) {
      var rest = function (r) { var s = 0; for (var k = 0; k < nItems; k++) if (k !== j && !drop[k]) s += r[k]; return s; };
      var all = rowsRef.map(rest).concat(rowsFoc.map(rest)), mn = Math.min.apply(null, all), mx = Math.max.apply(null, all), B = Math.max(1, Math.min(8, Math.round((mx - mn) / 2) + 1)), bw = (mx - mn) / B || 1;
      var bin = function (v) { return Math.min(B - 1, Math.floor((v - mn) / bw)); }, sR = {}, nR = {}, sF = {}, nF = {};
      rowsRef.forEach(function (r) { var b = bin(rest(r)); sR[b] = (sR[b] || 0) + r[j]; nR[b] = (nR[b] || 0) + 1; });
      rowsFoc.forEach(function (r) { var b = bin(rest(r)); sF[b] = (sF[b] || 0) + r[j]; nF[b] = (nF[b] || 0) + 1; });
      var num = 0, den = 0; Object.keys(nF).forEach(function (b) { if (!nR[b]) return; num += nF[b] * (sF[b] / nF[b] - sR[b] / nR[b]); den += nF[b]; });
      var std = den ? num / den : NaN, rel = Math.abs(std) / R;
      out.push({ std: std, level: !isFinite(std) ? 'na' : rel < 0.1 ? 'ok' : rel < 0.2 ? 'review' : 'flag' });
    }
    return out;
  }
  function selectorsHTML(sel, idp) {
    idp = idp || 'au';
    return '<div class="aud-sel"><label class="fl" for="' + idp + 'P">الفئة العمرية المستهدفة</label><select id="' + idp + 'P">' + PROFILES.map(function (p) { return '<option value="' + p.id + '"' + (sel && sel.p === p.id ? ' selected' : '') + '>' + p.t + ' — ' + p.g + '</option>'; }).join('') + '</select>' +
      '<label class="fl" for="' + idp + 'R">لغة الصياغة</label><select id="' + idp + 'R">' + REGISTERS.map(function (r) { return '<option value="' + r.id + '"' + (sel && sel.r === r.id ? ' selected' : '') + '>' + r.t + ' — ' + r.d + '</option>'; }).join('') + '</select></div>';
  }
  function tipsHTML(profileId) { var p = prof(profileId); return '<ul class="aud-tips">' + p.tips.map(function (t) { return '<li>' + t + '</li>'; }).join('') + '</ul>'; }
  return { PROFILES: PROFILES, REGISTERS: REGISTERS, prof: prof, adapt: adapt, check: check, labels: labels, instructions: instructions, prompt: prompt, stdDIF: stdDIF, selectorsHTML: selectorsHTML, tipsHTML: tipsHTML };
})();

/* ===== Construct library for real scale construction (shared by the Scales Lab and the Statistics Lab).
   Each construct is a blueprint: definition -> theory -> facets -> indicators (cognition c, emotion e, behaviour b,
   physiology p) -> direct items + indirect items (ambiguous situations, situational judgement, estimates about others,
   forced choice) + a cumulative ladder (Guttman) and, for attitudes, Thurstone statements with scale values.
   All item wordings are original drafts written for teaching and piloting; none is copied from a published instrument.
   Every draft still needs expert review (CVR/CVI), cognitive interviews and a pilot before use. ===== */
var CLIB = (function () {
  var K = { c: 'فكر', e: 'شعور', b: 'سلوك', p: 'جسد' };
  // d: [facet, kind, text, reversed]   amb/sjt: [stem, [[option, weight 0-2], ...]]   oth: [stem]   fc: [trait statement, neutral statement]
  var L = [
  { id: 'social_anx', n: 'القلق الاجتماعي', type: 'trait', leak: ['قلق اجتماعي', 'القلق الاجتماعي', 'رهاب'],
    def: 'خوف واضح ومستمر من مواقف التفاعل أو الأداء أمام الآخرين، مصدره توقع التقييم السلبي، ويظهر في أفكار عن الحكم السلبي، وأعراض جسدية، وتجنب أو تحمّل للموقف بضيق.',
    theory: 'النموذج المعرفي للقلق الاجتماعي (Clark & Wells, 1995) والنموذج المعرفي السلوكي (Rapee & Heimberg, 1997): الانتباه المتمركز على الذات، وتوقع التقييم السلبي، وسلوكيات الأمان، والمعالجة اللاحقة للموقف.',
    nomo: { pos: ['الخجل', 'الخوف من التقييم السلبي', 'الاكتئاب (ارتباط متوسط)'], neg: ['الانبساط', 'الكفاءة الاجتماعية المدركة'], none: ['القدرة العقلية العامة'] },
    facets: [['fne', 'الخوف من التقييم السلبي'], ['perf', 'قلق مواقف الأداء'], ['inter', 'قلق مواقف التفاعل'], ['post', 'المعالجة اللاحقة للموقف']],
    d: [['fne', 'c', 'أتوقع أن يلاحظ الناس أخطائي الصغيرة أثناء الكلام معهم.'], ['fne', 'e', 'يقلقني أن يكوّن الآخرون عني انطباعًا سيئًا من أول لقاء.'], ['fne', 'c', 'رأي الغرباء فيّ لا يشغلني كثيرًا.', 1],
       ['perf', 'p', 'عندما يُطلب مني الكلام أمام مجموعة أشعر بحرارة في وجهي.'], ['perf', 'b', 'أتجنب الأكل أو الشرب أمام أشخاص لا أعرفهم جيدًا.'], ['perf', 'b', 'أتطوع للإجابة في المحاضرة دون تردد.', 1],
       ['inter', 'b', 'أؤجل الاتصال بشخص لا أعرفه حتى لو كان الأمر مهمًا.'], ['inter', 'e', 'أرتبك عندما يبدأ معي شخص غريب حديثًا مفاجئًا.'], ['inter', 'b', 'أبدأ الحديث بسهولة مع شخص أقابله لأول مرة.', 1],
       ['post', 'c', 'بعد أي لقاء، أعيد في ذهني ما قلته وأبحث عن أخطائي.'], ['post', 'e', 'يبقى موقف محرج بسيط في ذهني أيامًا.'], ['post', 'c', 'أنسى سريعًا التفاصيل الصغيرة لما قلته في لقاء عادي.', 1]],
    amb: [['دخلتَ قاعة فيها مجموعة يتكلمون، فسكتوا فجأة. الأرجح أن:', [['كانوا يتكلمون عني', 2], ['انتهى حديثهم بالصدفة', 0], ['توقفوا ليرحبوا بي', 0]]],
          ['أرسلتَ رأيك في مجموعة الدفعة، ولم يعلّق أحد لساعة. الأرجح أن:', [['رأيي لم يعجبهم', 2], ['الكل مشغول الآن', 0], ['ربما لم ينتبهوا للرسالة', 0]]]],
    sjt: [['طلب منك الأستاذ أن تعرض بحثك أمام الدفعة الأسبوع القادم. أقرب تصرف لك:', [['أعتذر وأطلب تسليمه مكتوبًا', 2], ['أوافق لكن أتمرن كثيرًا وأقرأ من الورقة', 1], ['أوافق وأستعد بشكل عادي', 0]]],
          ['في مناسبة عائلية كبيرة، جلستَ بجوار شخص لا تعرفه. غالبًا:', [['أنشغل بالهاتف حتى لا أضطر للكلام', 2], ['أنتظر حتى يبدأ هو الكلام', 1], ['أعرّفه بنفسي وأبدأ حديثًا', 0]]]],
    oth: [['من كل 10 طلاب في سنك، كم واحدًا تتوقع أنه يتجنب الكلام في المحاضرة خوفًا من الإحراج؟']],
    fc: [['أفضّل أن أعمل في الخلفية بعيدًا عن الأنظار.', 'أفضّل أن أعمل في مكان هادئ ومنظم.'], ['أفكر كثيرًا فيما يظنه الآخرون بي.', 'أفكر كثيرًا في خططي للمستقبل.']],
    ladder: ['أتردد قليلًا قبل أن أكلم شخصًا جديدًا.', 'أتجنب أن أكون مركز الانتباه في الجلسات.', 'أعتذر أحيانًا عن مناسبات اجتماعية بسبب التوتر.', 'أتغيب عن محاضرات فيها مشاركة شفهية.', 'يمنعني التوتر من أنشطة مهمة في دراستي أو عملي.'] },

  { id: 'test_anx', n: 'قلق الاختبار', type: 'trait', leak: ['قلق الاختبار', 'قلق الامتحان'],
    def: 'استجابة انفعالية ومعرفية وجسدية لمواقف التقييم الأكاديمي، تشمل الانشغال بالفشل وتوتر الجسم وتشتت الانتباه، وتعطل الأداء بدرجات مختلفة.',
    theory: 'التمييز بين مكوني الانزعاج (القلق المعرفي) والانفعالية (الاستثارة الجسدية) عند Liebert & Morris (1967)، ونظرية التحكم الانتباهي (Eysenck et al., 2007).',
    nomo: { pos: ['القلق العام', 'الكمالية غير التكيفية'], neg: ['الكفاءة الذاتية الأكاديمية'], none: ['الانبساط'] },
    facets: [['worry', 'الانزعاج المعرفي'], ['emo', 'الانفعالية الجسدية'], ['interf', 'تداخل الأفكار مع الأداء'], ['avoid', 'التجنب والتأجيل']],
    d: [['worry', 'c', 'أثناء الامتحان أفكر فيما سيحدث لو رسبت.'], ['worry', 'c', 'أقارن نفسي بزملائي وأنا أحل.'], ['worry', 'c', 'أدخل الامتحان وأنا متوقع أن أحل بشكل جيد.', 1],
       ['emo', 'p', 'قبل الامتحان بساعات تتسارع ضربات قلبي.'], ['emo', 'p', 'لا أستطيع الأكل جيدًا صباح يوم الامتحان.'], ['emo', 'p', 'أنام بشكل طبيعي ليلة الامتحان.', 1],
       ['interf', 'c', 'أنسى معلومات كنت أعرفها جيدًا بمجرد أن أرى ورقة الأسئلة.'], ['interf', 'b', 'أقرأ السؤال مرات كثيرة قبل أن أفهمه.'], ['interf', 'c', 'أركز في السؤال الذي أمامي دون أن يشتتني الوقت.', 1],
       ['avoid', 'b', 'أؤجل مراجعة المادة الصعبة حتى آخر وقت.'], ['avoid', 'b', 'أتمنى لو أُلغي الامتحان حتى لو ذاكرت جيدًا.']],
    amb: [['وزّع المراقب ورقة الأسئلة، وأول سؤال بدا لك غير مألوف. الأرجح أن:', [['الامتحان كله سيكون صعبًا', 2], ['هذا سؤال واحد وسأعود إليه', 0], ['ربما صيغته فقط مختلفة', 0]]],
          ['خرج زميل من الامتحان مبكرًا. أول فكرة تخطر لك:', [['هو حل كل شيء وأنا متأخر', 2], ['ربما لم يذاكر', 0], ['لكل شخص سرعته', 0]]]],
    sjt: [['الليلة قبل امتحان مهم، وما زال جزء لم تراجعه. غالبًا:', [['أسهر حتى الفجر وأنا متوتر', 2], ['أراجع ملخصًا سريعًا ثم أنام', 0], ['أراجعه وأقلق من نسيان ما سبق', 1]]]],
    oth: [['من كل 10 طلاب، كم واحدًا تتوقع أنه ينسى معلومات يعرفها بسبب التوتر في الامتحان؟']],
    fc: [['أتوتر بشدة عندما أعرف أنني سأُقيَّم.', 'أتضايق عندما يتغير جدولي فجأة.']],
    ladder: ['أشعر ببعض التوتر قبل الامتحانات.', 'يقل نومي ليلة الامتحان.', 'أنسى معلومات أعرفها أثناء الامتحان.', 'أفكر في التغيب عن امتحانات بسبب التوتر.', 'تغيبت فعلًا عن امتحان بسبب التوتر.'] },

  { id: 'procr', n: 'التسويف الأكاديمي', type: 'trait', leak: ['تسويف', 'التسويف', 'أسوّف'],
    def: 'تأجيل طوعي وغير ضروري لبدء أو إنهاء مهام أكاديمية مقصودة، رغم توقع أن التأجيل سيضر، ويصاحبه غالبًا ضيق أو شعور بالذنب.',
    theory: 'تعريف Steel (2007) ونظرية الدافعية الزمنية، ونموذج Sirois & Pychyl (2013) للتسويف كتنظيم قصير المدى للمزاج على حساب الهدف.',
    nomo: { pos: ['الاندفاعية', 'ضعف الضبط الذاتي', 'قلق الاختبار (ضعيف إلى متوسط)'], neg: ['يقظة الضمير', 'الكفاءة الذاتية الأكاديمية'], none: ['الذكاء'] },
    facets: [['delay', 'تأجيل البدء'], ['gap', 'الفجوة بين النية والفعل'], ['mood', 'التأجيل لتحسين المزاج'], ['cons', 'الضيق من العواقب']],
    d: [['delay', 'b', 'أبدأ في الواجبات قبل موعد تسليمها بيوم أو يومين.'], ['delay', 'b', 'أبدأ المهام الكبيرة بمجرد أن تُطلب مني.', 1], ['delay', 'b', 'أترك قراءة المراجع حتى يقترب الامتحان.'],
       ['gap', 'c', 'أخطط لأذاكر في وقت معين ثم أجد نفسي أفعل شيئًا آخر.'], ['gap', 'b', 'أنفذ جدول المذاكرة الذي أضعه في الغالب.', 1], ['gap', 'c', 'أقول لنفسي «سأبدأ بعد قليل» مرات كثيرة في اليوم.'],
       ['mood', 'b', 'عندما تكون المهمة مملة أفتح الهاتف لأرتاح قليلًا فيطول الأمر.'], ['mood', 'e', 'أشعر براحة مؤقتة عندما أقرر تأجيل مهمة صعبة.'],
       ['cons', 'e', 'ألوم نفسي بعد أن يضيع وقت كنت قد خصصته للمذاكرة.'], ['cons', 'b', 'سلّمت عملًا أقل من مستواي لأنني بدأت متأخرًا.']],
    amb: [['أمامك تكليف مطلوب بعد أسبوعين. أول ما يخطر لك:', [['ما زال هناك وقت طويل', 2], ['سأقسمه على الأيام من الآن', 0], ['سأبدأ بالجزء الأسهل اليوم', 0]]]],
    sjt: [['جلستَ للمذاكرة، ووصلتك رسالة من صديق يقترح خروجًا الآن. غالبًا:', [['أخرج وأقول لنفسي سأذاكر ليلًا', 2], ['أؤجل الخروج لساعتين حتى أنهي جزءًا', 0], ['أخرج لكن أشعر بالذنب طول الوقت', 1]]],
          ['كُلّفت ببحث لا تعرف من أين تبدأ فيه. أقرب تصرف:', [['أنتظر حتى تتضح الفكرة وحدها', 2], ['أسأل زميلًا أو الأستاذ اليوم', 0], ['أبحث سريعًا عن خطوات البدء', 0]]]],
    oth: [['من كل 10 طلاب في دفعتك، كم واحدًا تتوقع أنه يؤجل المذاكرة لآخر أسبوع قبل الامتحان؟']],
    fc: [['أنجز أفضل تحت ضغط الوقت الأخير.', 'أنجز أفضل عندما يكون المكان هادئًا.'], ['أبدأ المهام عندما «أكون في المزاج».', 'أبدأ المهام بعد أن أفهم المطلوب بدقة.']],
    ladder: ['أؤجل بعض المهام البسيطة.', 'أبدأ معظم الواجبات قرب موعدها.', 'سلّمت بعض الأعمال متأخرًا بسبب التأجيل.', 'خسرت درجات واضحة بسبب التأجيل.', 'رسبت أو أعدت مقررًا بسبب التأجيل.'] },

  { id: 'aseff', n: 'الكفاءة الذاتية الأكاديمية', type: 'trait', leak: ['كفاءة ذاتية', 'الكفاءة الذاتية'],
    def: 'اعتقاد الطالب في قدرته على تنظيم وتنفيذ الأفعال اللازمة لتحقيق أداء أكاديمي محدد، وهو اعتقاد خاص بالمهمة وليس تقديرًا عامًا للذات.',
    theory: 'نظرية الكفاءة الذاتية عند Bandura (1997): مصادرها الخبرات الناجحة، والخبرة البديلة، والإقناع اللفظي، والحالة الانفعالية. ويُفضّل قياسها مرتبطة بمهام محددة.',
    nomo: { pos: ['الأداء الأكاديمي', 'الدافعية الداخلية'], neg: ['قلق الاختبار', 'التسويف'], none: ['الطول (ضابط تمييزي)'] },
    facets: [['task', 'إنجاز مهام صعبة'], ['self', 'التنظيم الذاتي للتعلم'], ['persist', 'المثابرة أمام العقبات']],
    d: [['task', 'c', 'أستطيع فهم أصعب موضوع في المقرر إذا أعطيته وقتًا كافيًا.'], ['task', 'c', 'أستطيع حل مسألة إحصائية جديدة بخطوات منطقية.'], ['task', 'c', 'أشك في قدرتي على اجتياز المقررات الصعبة.', 1],
       ['self', 'b', 'أستطيع وضع جدول مذاكرة والالتزام به معظم الأسبوع.'], ['self', 'c', 'أعرف كيف أذاكر بالطريقة التي تناسبني.'], ['self', 'c', 'أحتاج من يتابعني حتى أذاكر.', 1],
       ['persist', 'b', 'عندما لا أفهم درسًا أبحث عن شرح آخر حتى أفهمه.'], ['persist', 'c', 'إذا حصلت على درجة سيئة أعرف ما أعدّله في المرة القادمة.'], ['persist', 'b', 'أترك المسألة إذا لم أصل للحل في المحاولة الأولى.', 1]],
    amb: [['حصلت على درجة أقل من المتوقع في اختبار قصير. الأرجح أن السبب:', [['أنني لست جيدًا في هذه المادة', 2], ['طريقة مذاكرتي تحتاج تعديلًا', 0], ['الاختبار ركّز على جزء لم أراجعه', 0]]]],
    sjt: [['أعلن الأستاذ مشروعًا اختياريًا صعبًا يرفع الدرجة. غالبًا:', [['لا أشارك لأنني لن أنجح فيه', 2], ['أشارك إذا شاركني زميل قوي', 1], ['أشارك وأبحث عن المطلوب', 0]]]],
    oth: [['من كل 10 طلاب في دفعتك، كم واحدًا تتوقع أنه يستطيع النجاح في أصعب مقرر بمجهوده؟']],
    fc: [['أثق أنني أستطيع تعلم أي مهارة دراسية جديدة.', 'أستمتع بتعلم المهارات الدراسية الجديدة.']],
    ladder: ['أستطيع إنجاز المهام السهلة وحدي.', 'أستطيع فهم معظم الدروس بمجهودي.', 'أستطيع حل مسائل لم أرَ مثلها.', 'أستطيع التفوق في أصعب مقرر.', 'أستطيع شرح أصعب الموضوعات لغيري.'] },

  { id: 'lonely', n: 'الوحدة النفسية', type: 'trait', leak: ['وحدة', 'الوحدة', 'وحيد'],
    def: 'خبرة ذاتية مؤلمة تنتج عن إدراك فجوة بين العلاقات الاجتماعية المرغوبة والعلاقات الموجودة فعلًا، كمًّا أو كيفًا، ولا تساوي العزلة الموضوعية.',
    theory: 'نموذج التناقض المعرفي (Perlman & Peplau, 1981) والتمييز بين الوحدة الاجتماعية والوحدة الانفعالية (Weiss, 1973)، ونموذج Cacioppo لليقظة المفرطة للتهديد الاجتماعي.',
    nomo: { pos: ['الاكتئاب', 'القلق الاجتماعي'], neg: ['المساندة الاجتماعية المدركة', 'الرضا عن الحياة'], none: ['عدد ساعات الدراسة'] },
    facets: [['emo', 'الوحدة الانفعالية (غياب علاقة قريبة)'], ['soc', 'الوحدة الاجتماعية (غياب جماعة)'], ['vig', 'اليقظة للرفض']],
    d: [['emo', 'e', 'أفتقد شخصًا أحكي له تفاصيل يومي.'], ['emo', 'c', 'يوجد على الأقل شخص واحد يفهمني بلا شرح طويل.', 1], ['emo', 'e', 'أشعر بأنني وحدي حتى وأنا وسط أهلي.'],
       ['soc', 'b', 'أعرف إلى من أتصل إذا أردت الخروج في عطلة.', 1], ['soc', 'e', 'أشعر أنني لا أنتمي لأي مجموعة في الكلية.'], ['soc', 'b', 'أقضي أغلب أوقات الفراغ وحدي دون أن أختار ذلك.'],
       ['vig', 'c', 'إذا تأخر صديق في الرد أظن أنه لم يعد يهتم بي.'], ['vig', 'c', 'أتوقع أن الناس يقبلونني عندما أتعرف عليهم.', 1]],
    amb: [['رأيت صورًا لأصدقائك خرجوا معًا ولم يدعوك. الأرجح أن:', [['لم يعد لي مكان بينهم', 2], ['كانت خروجة مرتبة على عجل', 0], ['ظنوا أنني مشغول', 0]]]],
    sjt: [['يوم إجازة طويلة، وليس لديك خطط. غالبًا:', [['أبقى في البيت أتابع ما يفعله الآخرون على الإنترنت', 2], ['أراسل صديقًا وأقترح لقاء', 0], ['أنتظر أن يدعوني أحد', 1]]]],
    oth: [['من كل 10 أشخاص في سنك، كم واحدًا تتوقع أنه يشعر بأنه بلا علاقة قريبة فعلًا؟']],
    fc: [['أشعر أن علاقاتي أقل عمقًا مما أتمنى.', 'أشعر أن وقتي أقل مما أحتاج.']],
    ladder: ['أتمنى أحيانًا علاقات أكثر.', 'أقضي معظم عطلاتي وحدي.', 'لا أجد من أحكي له مشكلة شخصية.', 'أشعر بالوحدة معظم الأيام.', 'تمنعني الوحدة من الاستمتاع بأي نشاط.'] },

  { id: 'resil', n: 'الصمود النفسي', type: 'trait', leak: ['صمود', 'الصمود', 'المرونة النفسية'],
    def: 'القدرة على التعافي والعودة لمستوى أداء ورضا مقبول بعد الضغوط والانتكاسات، وتظهر في التعامل النشط، والمعنى، والمرونة المعرفية، واستخدام المساندة.',
    theory: 'الصمود كعملية تكيف ديناميكية (Masten, 2001؛ Bonanno, 2004) وليس سمة ثابتة فقط، ولذلك تُفضل المؤشرات المرتبطة بخبرات فعلية بعد الضغط.',
    nomo: { pos: ['الرضا عن الحياة', 'الكفاءة الذاتية'], neg: ['الاكتئاب', 'الاجترار'], none: ['الطول'] },
    facets: [['recov', 'سرعة التعافي'], ['cope', 'التعامل النشط'], ['meaning', 'إيجاد المعنى'], ['support', 'طلب المساندة']],
    d: [['recov', 'b', 'بعد خيبة أمل كبيرة أعود لروتيني خلال أيام قليلة.'], ['recov', 'e', 'تظل المواقف الصعبة تثقل يومي لفترة طويلة.', 1],
       ['cope', 'b', 'عندما تواجهني مشكلة أكتب خطوات صغيرة لحلها.'], ['cope', 'b', 'أنتظر أن تنحل المشكلات وحدها.', 1], ['cope', 'c', 'أبحث عن أكثر من طريقة للخروج من موقف صعب.'],
       ['meaning', 'c', 'تعلمت من أصعب تجربة مررت بها شيئًا مفيدًا.'], ['meaning', 'c', 'أرى الضغوط الحالية مرحلة مؤقتة.'],
       ['support', 'b', 'أطلب المساعدة من شخص أثق به عندما أحتاجها.'], ['support', 'b', 'أخفي مشكلاتي عن الجميع حتى تشتد.', 1]],
    amb: [['تعطّل حاسوبك قبل تسليم مشروع بيوم. أول ما يخطر لك:', [['ضاع كل شيء ولن ألحق', 2], ['سأبحث عن نسخة أو جهاز بديل', 0], ['سأشرح الموقف للأستاذ وأطلب مهلة', 0]]]],
    sjt: [['لم تُقبل في منحة كنت تنتظرها. خلال الأسبوع التالي غالبًا:', [['أتوقف عن التقديم لأي شيء', 2], ['أحزن يومين ثم أبحث عن فرصة أخرى', 0], ['أتابع حياتي لكن دون حماس', 1]]]],
    oth: [['من كل 10 أشخاص يمرون بفشل كبير، كم واحدًا تتوقع أنه يعود لحياته الطبيعية خلال شهر؟']],
    fc: [['أعرف كيف أنهض بعد الفشل.', 'أعرف كيف أنظم وقتي.']],
    ladder: ['أتعافى من المضايقات اليومية بسهولة.', 'أعود لروتيني بعد خلاف مع صديق.', 'أتعافى بعد فشل دراسي خلال أسابيع.', 'أجد معنى حتى في خسارة كبيرة.', 'أساعد غيري على التعافي من أزماتهم.'] },

  { id: 'perf', n: 'الكمالية غير التكيفية', type: 'trait', leak: ['كمالية', 'الكمالية', 'مثالي'],
    def: 'وضع معايير مرتفعة جدًا مع نقد ذاتي قاسٍ عند عدم بلوغها، وقلق مستمر من الخطأ، وتقييم الذات بناءً على الأداء.',
    theory: 'التمييز بين الجهود الكمالية والمخاوف الكمالية (Stoeber & Otto, 2006)، والنموذج متعدد الأبعاد (Hewitt & Flett, 1991). المخاوف الكمالية هي الجزء غير التكيفي المرتبط بالضيق.',
    nomo: { pos: ['القلق', 'الاكتئاب', 'الاحتراق'], neg: ['الشفقة بالذات'], none: ['الانبساط'] },
    facets: [['crit', 'النقد الذاتي'], ['mist', 'الانشغال بالأخطاء'], ['doubt', 'الشك في الأداء'], ['worth', 'ربط القيمة بالإنجاز']],
    d: [['crit', 'c', 'إذا أخطأت في جزء صغير أشعر أن العمل كله فاشل.'], ['crit', 'e', 'أتقبل أخطائي كجزء طبيعي من التعلم.', 1], ['mist', 'b', 'أراجع الرسالة قبل إرسالها مرات كثيرة.'],
       ['mist', 'e', 'يضايقني خطأ إملائي صغير لساعات.'], ['doubt', 'c', 'حتى بعد أن أنهي عملًا جيدًا أشعر أنه ناقص.'], ['doubt', 'b', 'أؤجل تسليم العمل لأنني غير راضٍ عنه بعد.'],
       ['worth', 'c', 'أشعر بقيمتي فقط عندما أحقق نتيجة ممتازة.'], ['worth', 'c', 'أحترم نفسي حتى عندما تكون نتيجتي متوسطة.', 1]],
    amb: [['قال الأستاذ إن بحثك «جيد». أول ما تفهمه:', [['يعني أنه ليس ممتازًا وفيه عيوب', 2], ['يعني أنه جيد فعلًا', 0], ['ربما يقصد أنه يحتاج لمسات بسيطة', 1]]]],
    sjt: [['أنهيت عرضًا تقديميًا قبل الموعد بيوم وهو جيد. غالبًا:', [['أعيد تصميمه كله في الليلة الأخيرة', 2], ['أراجعه مرة وأسلّمه', 0], ['أعدل التفاصيل حتى آخر دقيقة', 1]]]],
    oth: [['من كل 10 طلاب، كم واحدًا تتوقع أنه يعتبر أي درجة أقل من الامتياز فشلًا؟']],
    fc: [['لا أرتاح حتى يكون عملي بلا أي خطأ.', 'لا أرتاح حتى يكون مكاني مرتبًا.']],
    ladder: ['أحب أن يكون عملي جيدًا.', 'أراجع عملي أكثر من مرة.', 'أتضايق كثيرًا من أخطائي الصغيرة.', 'أؤجل التسليم لأنه «ليس كاملًا».', 'تمنعني المعايير العالية من إنهاء أعمال كثيرة.'] },

  { id: 'burn', n: 'الاحتراق الأكاديمي', type: 'trait', leak: ['احتراق', 'الاحتراق'],
    def: 'حالة مزمنة من الإنهاك المرتبط بالدراسة، مع تبلد أو بعد عن الدراسة، وانخفاض الشعور بالكفاءة الأكاديمية، نتيجة ضغوط دراسية مستمرة.',
    theory: 'نموذج Maslach للاحتراق بأبعاده الثلاثة كما طُبق على الطلاب (Schaufeli et al., 2002)، ونموذج متطلبات–موارد الدراسة.',
    nomo: { pos: ['الضغط المدرك', 'الاكتئاب'], neg: ['الانخراط الدراسي', 'الصمود'], none: ['الطول'] },
    facets: [['exh', 'الإنهاك'], ['cyn', 'التبلد والبعد عن الدراسة'], ['ineff', 'انخفاض الكفاءة']],
    d: [['exh', 'p', 'أستيقظ متعبًا عندما أفكر في يوم دراسي طويل.'], ['exh', 'e', 'أشعر أن طاقتي تنفد قبل منتصف الأسبوع.'], ['exh', 'p', 'أنهي اليوم الدراسي وعندي طاقة لأنشطة أخرى.', 1],
       ['cyn', 'c', 'لم أعد أرى فائدة لما أدرسه.'], ['cyn', 'e', 'فقدت حماسي للمقررات التي كنت أحبها.'], ['cyn', 'c', 'أشعر أن ما أدرسه له معنى لمستقبلي.', 1],
       ['ineff', 'c', 'أشعر أنني أبذل مجهودًا بلا نتيجة.'], ['ineff', 'c', 'أستطيع حل مشكلات دراستي بفاعلية.', 1]],
    amb: [['أعلن الأستاذ تكليفًا إضافيًا. أول رد فعل داخلي:', [['لا أستطيع تحمل أي شيء آخر', 2], ['سأنظم وقتي له', 0], ['ربما يكون مفيدًا', 0]]]],
    sjt: [['يوم إجازة بعد أسبوع امتحانات. غالبًا:', [['أنام معظم اليوم ولا أرغب في أي شيء', 2], ['أخرج أو أمارس هواية', 0], ['أرتاح نصف اليوم وأخطط للأسبوع القادم', 0]]]],
    oth: [['من كل 10 طلاب، كم واحدًا تتوقع أنه يشعر بإنهاك دائم من الدراسة؟']],
    fc: [['أشعر أن الدراسة تستنزف طاقتي.', 'أشعر أن المواصلات تستهلك وقتي.']],
    ladder: ['أتعب أحيانًا من الدراسة.', 'أشعر بالإنهاك معظم الأسابيع.', 'فقدت اهتمامي بأغلب المقررات.', 'أفكر في ترك الدراسة.', 'تغيبت لفترات طويلة بسبب الإنهاك.'] },

  { id: 'phone', n: 'الاستخدام الإشكالي للهاتف الذكي', type: 'trait', leak: ['إدمان', 'الهاتف الذكي', 'إدمان الهاتف'],
    def: 'نمط استخدام للهاتف يصعب التحكم فيه ويؤدي إلى تأثيرات سلبية في النوم أو الدراسة أو العلاقات، مع انشغال وانزعاج عند تعذر الاستخدام.',
    theory: 'نموذج التفاعل بين الشخص والانفعال والمعرفة والتنفيذ (I-PACE، Brand et al., 2016). ويحذّر Billieux et al. (2015) من نقل معايير الإدمان آليًا، ولذلك تركز المؤشرات على فقدان التحكم والتأثير الوظيفي.',
    nomo: { pos: ['الاندفاعية', 'الخوف من فوات الأشياء (FOMO)', 'ضعف جودة النوم'], neg: ['الضبط الذاتي'], none: ['المستوى الاقتصادي (ضابط)'] },
    facets: [['ctrl', 'فقدان التحكم'], ['impair', 'التأثير الوظيفي'], ['withd', 'الانزعاج عند الانقطاع']],
    d: [['ctrl', 'b', 'أفتح الهاتف لأمر بسيط فأجد نصف ساعة مرت.'], ['ctrl', 'b', 'أستطيع ترك الهاتف في غرفة أخرى أثناء المذاكرة.', 1], ['ctrl', 'b', 'أفحص الهاتف أول ما أستيقظ قبل أي شيء.'],
       ['impair', 'b', 'نمت متأخرًا أكثر من مرة هذا الأسبوع بسبب الهاتف.'], ['impair', 'b', 'فاتني جزء من شرح المحاضرة لأنني كنت على الهاتف.'], ['impair', 'b', 'أنام في موعدي حتى لو كان الهاتف بجانبي.', 1],
       ['withd', 'e', 'أتوتر إذا نسيت الهاتف في البيت.'], ['withd', 'p', 'أشعر أحيانًا أن الهاتف يهتز وهو لا يهتز.']],
    amb: [['انقطع الإنترنت في البيت لساعتين. أول شعور:', [['ضيق شديد وأبحث عن أي شبكة', 2], ['فرصة لأفعل شيئًا آخر', 0], ['انزعاج بسيط ثم أتأقلم', 1]]]],
    sjt: [['تجلس مع أسرتك على العشاء، وسمعت صوت إشعار. غالبًا:', [['أنظر فورًا', 2], ['أنظر بعد العشاء', 0], ['أنظر سريعًا دون أن يلاحظوا', 1]]]],
    oth: [['من كل 10 أشخاص في سنك، كم واحدًا تتوقع أنه لا يستطيع ترك الهاتف ساعة كاملة؟']],
    fc: [['أمسك الهاتف دون أن أقرر ذلك.', 'أشرب القهوة دون أن أقرر ذلك.']],
    ladder: ['أستخدم الهاتف كثيرًا في الفراغ.', 'أفحصه كل بضع دقائق.', 'يقلل نومي بانتظام.', 'يؤثر على درجاتي بوضوح.', 'سبّب لي مشكلات حقيقية مع أسرتي أو عملي.'] },

  { id: 'emreg', n: 'صعوبات تنظيم الانفعال', type: 'trait', leak: ['تنظيم الانفعال', 'انفعال'],
    def: 'صعوبة في الوعي بالانفعالات وفهمها وتقبلها، والتحكم في السلوك الاندفاعي والاستمرار في الهدف أثناء الضيق، واستخدام استراتيجيات مرنة لتعديل الانفعال.',
    theory: 'النموذج متعدد الأبعاد لـ Gratz & Roemer (2004)، ونموذج العملية لـ Gross (1998) الذي يفرق بين إعادة التقييم والكبت.',
    nomo: { pos: ['القلق', 'الاندفاعية'], neg: ['الوعي الانفعالي', 'الصمود'], none: ['الطول'] },
    facets: [['aware', 'ضعف الوعي الانفعالي'], ['impulse', 'صعوبة التحكم في الاندفاع'], ['goal', 'صعوبة الاستمرار في الهدف'], ['strat', 'محدودية الاستراتيجيات']],
    d: [['aware', 'c', 'أعرف بالضبط ما أشعر به عندما أتضايق.', 1], ['aware', 'c', 'أحتار في تسمية ما أشعر به عندما أكون منزعجًا.'],
       ['impulse', 'b', 'عندما أغضب أقول كلامًا أندم عليه.'], ['impulse', 'b', 'أستطيع التوقف قبل أن أرد وأنا منفعل.', 1],
       ['goal', 'b', 'إذا تضايقت لا أستطيع التركيز في المذاكرة لساعات.'], ['goal', 'c', 'أكمل مهامي حتى وأنا حزين.', 1],
       ['strat', 'c', 'عندما أتضايق أشعر أنه لا يوجد ما يخفف ذلك.'], ['strat', 'b', 'عندي أكثر من طريقة لتهدئة نفسي.', 1]],
    amb: [['ردّ عليك صديق برد مقتضب. أول تفسير:', [['هو غاضب مني', 2], ['هو مشغول', 0], ['ربما يمر بيوم سيئ', 0]]]],
    sjt: [['تلقيت نقدًا حادًا أمام زملائك. في الساعة التالية غالبًا:', [['أنشغل بالموقف ولا أستطيع فعل أي شيء', 2], ['أتضايق ثم أفكر في المفيد في النقد', 0], ['أرد بحدة ثم أندم', 2]]]],
    oth: [['من كل 10 أشخاص، كم واحدًا تتوقع أنه يفقد السيطرة على نفسه عندما يغضب؟']],
    fc: [['يصعب عليّ أن أهدأ بعد موقف مزعج.', 'يصعب عليّ أن أستيقظ مبكرًا.']],
    ladder: ['أتضايق أحيانًا من مواقف صغيرة.', 'يأخذ هدوئي وقتًا بعد الضيق.', 'يتعطل يومي بسبب انفعال.', 'أتصرف باندفاع أندم عليه كثيرًا.', 'تسبب انفعالاتي مشكلات كبيرة في حياتي.'] },

  { id: 'selfest', n: 'تقدير الذات', type: 'trait', leak: ['تقدير الذات', 'تقدير ذاتي'],
    def: 'تقييم الفرد العام لقيمته كشخص، أي درجة تقبله واحترامه لنفسه، وهو مختلف عن الكفاءة في مهام محددة.',
    theory: 'تعريف Rosenberg (1965) للتقدير العام للذات، ونظرية المقياس الاجتماعي (Leary, 1999): التقدير مؤشر للقبول الاجتماعي المدرك.',
    nomo: { pos: ['الرضا عن الحياة', 'التفاؤل'], neg: ['الاكتئاب', 'القلق الاجتماعي'], none: ['الذكاء'] },
    facets: [['worth', 'الإحساس بالقيمة'], ['accept', 'تقبل الذات']],
    d: [['worth', 'c', 'أرى أن لدي صفات جيدة يقدّرها الآخرون.'], ['worth', 'c', 'أشعر أحيانًا أنني أقل من غيري في كل شيء.', 1], ['worth', 'e', 'أنظر لنفسي باحترام.'],
       ['accept', 'c', 'أتمنى لو كنت شخصًا آخر.', 1], ['accept', 'e', 'أنا راضٍ عن نفسي بوجه عام.'], ['accept', 'c', 'أتقبل عيوبي دون أن تقلل من قيمتي.']],
    amb: [['مدحك شخص على عمل قمت به. أول ما تفكر فيه:', [['ربما يجاملني فقط', 2], ['أشكره وأفرح', 0], ['أتذكر الأجزاء الناقصة', 1]]]],
    sjt: [['اختلفتَ مع مجموعة في رأي وأنت مقتنع به. غالبًا:', [['أتراجع حتى لا أبدو مخطئًا', 2], ['أشرح رأيي بهدوء', 0], ['أسكت وأحتفظ برأيي', 1]]]],
    oth: [['من كل 10 أشخاص في سنك، كم واحدًا تتوقع أنه غير راضٍ عن نفسه؟']],
    fc: [['أشعر أنني أستحق الاحترام كما أنا.', 'أشعر أنني أستحق إجازة الآن.']],
    ladder: ['أرى في نفسي بعض الصفات الجيدة.', 'أحترم نفسي في معظم المواقف.', 'أتقبل عيوبي.', 'أثق في قيمتي حتى عند الفشل.', 'لا يهتز إحساسي بقيمتي من النقد.'] },

  { id: 'rumin', n: 'الاجترار الفكري', type: 'trait', leak: ['اجترار', 'الاجترار'],
    def: 'تفكير متكرر وسلبي في مشاعر الضيق وأسبابها ونتائجها، دون الانتقال إلى حل أو فعل، ويطيل الحالة المزاجية السلبية.',
    theory: 'نظرية أساليب الاستجابة (Nolen-Hoeksema, 1991) والتمييز بين التأمل والاجترار القسري (Treynor et al., 2003).',
    nomo: { pos: ['الاكتئاب', 'القلق'], neg: ['حل المشكلات', 'اليقظة الذهنية'], none: ['الطول'] },
    facets: [['brood', 'الاجترار القسري'], ['why', 'أسئلة «لماذا» المتكررة'], ['stuck', 'صعوبة الانتقال للفعل']],
    d: [['brood', 'c', 'أظل أفكر في موقف مزعج حتى بعد أن ينتهي بيوم.'], ['brood', 'c', 'أترك المواقف المزعجة خلفي بسرعة.', 1], ['brood', 'p', 'تبقيني الأفكار مستيقظًا عندما أريد النوم.'],
       ['why', 'c', 'أسأل نفسي كثيرًا «لماذا يحدث لي هذا؟».'], ['why', 'c', 'أعيد تحليل كلام الآخرين معي مرات عديدة.'],
       ['stuck', 'b', 'أفكر في المشكلة كثيرًا لكنني لا أبدأ في حلها.'], ['stuck', 'b', 'عندما أفكر في مشكلة أنتقل سريعًا لخطوة عملية.', 1]],
    amb: [['تذكرتَ فجأة موقفًا محرجًا قديمًا. غالبًا ما يحدث:', [['أعيشه في ذهني مرات كثيرة', 2], ['أبتسم وأكمل ما أفعله', 0], ['أتضايق قليلًا ثم أنساه', 0]]]],
    sjt: [['قبل النوم تذكرت خلافًا مع صديق. غالبًا:', [['أظل ساعة أو أكثر أفكر فيه', 2], ['أقرر أن أكلمه غدًا وأنام', 0], ['أكتب ما يضايقني ثم أنام', 0]]]],
    oth: [['من كل 10 أشخاص، كم واحدًا تتوقع أنه يقضي وقتًا طويلًا في التفكير في مشكلاته دون حلها؟']],
    fc: [['تدور الأفكار نفسها في رأسي كثيرًا.', 'تدور الأغاني نفسها في رأسي كثيرًا.']],
    ladder: ['أفكر في المواقف المزعجة أحيانًا.', 'أعيد التفكير فيها لساعات.', 'تؤثر على نومي.', 'تعطل يومي الدراسي.', 'لا أستطيع إيقافها معظم الأيام.'] },

  /* ===== M2: original construct under development (Mohamed Osman) ===== */
  { id: 'm2_aidep', n: 'الاعتماد المعرفي على الذكاء الاصطناعي', type: 'trait', dev: true, leak: ['اعتماد', 'الاعتماد المعرفي', 'تفريغ معرفي'],
    def: 'نمط من تفويض العمليات المعرفية (التفكير، والتخطيط، والكتابة، والحكم) إلى أدوات الذكاء الاصطناعي التوليدي يتجاوز ما تتطلبه المهمة، ويصاحبه ضعف التحقق من مخرجاتها، وتراجع الثقة في القدرة الذاتية بدونها، وقد يمتد إلى اعتماد انفعالي وفقدان للتحكم.',
    theory: 'التفريغ المعرفي (Risko & Gilbert, 2016) ونظرية العقل الممتد (Clark & Chalmers, 1998) يفسران متى يكون التفويض تكيفيًا؛ وتحيز الأتمتة (Parasuraman & Manzey, 2010) ومعايرة الثقة (Lee & See, 2004) يفسران متى يصبح ضارًا؛ والدراسات الحديثة للاعتماد على الذكاء الاصطناعي التوليدي (مثل AIDep-22 لدى طلاب صينيين، ومقياس Goh & Hartanto، 2025) تؤكد تعدد الأبعاد وارتباطه السلبي بالتفكير النقدي (ارتباطي حتى الآن).',
    nomo: { pos: ['الاستخدام الإشكالي للهاتف', 'التسويف', 'انخفاض الحاجة للمعرفة (Need for Cognition)'], neg: ['التفكير النقدي', 'الكفاءة الذاتية الأكاديمية'], none: ['تكرار الاستخدام وحده (يجب أن يتميز عنه)'] },
    facets: [['deleg', 'التفويض المعرفي الزائد'], ['verify', 'ضعف التحقق ومعايرة الثقة'], ['eff', 'تآكل الثقة في القدرة الذاتية'], ['emo', 'الاعتماد الانفعالي'], ['ctrl', 'فقدان التحكم']],
    d: [['deleg', 'b', 'أطلب من الذكاء الاصطناعي حل المسألة قبل أن أحاول فيها بنفسي.'], ['deleg', 'b', 'أكتب أفكاري الأولى بنفسي ثم أستعين بالأدوات للتحسين.', 1], ['deleg', 'b', 'أطلب خطة اليوم أو المذاكرة جاهزة بدل أن أفكر فيها.'],
       ['verify', 'b', 'أنقل الإجابة التي أحصل عليها كما هي دون مراجعة مصادرها.'], ['verify', 'b', 'أتحقق من المعلومات المهمة في مرجع آخر قبل استخدامها.', 1], ['verify', 'c', 'إذا كانت الإجابة مكتوبة بثقة أفترض أنها صحيحة.'],
       ['eff', 'c', 'أشك في قدرتي على كتابة بحث جيد دون مساعدة هذه الأدوات.'], ['eff', 'e', 'أشعر بالارتباك إذا طُلب مني حل مسألة في الامتحان بلا أي أداة.'], ['eff', 'c', 'أثق أنني أستطيع التفكير في المشكلات الصعبة وحدي.', 1],
       ['emo', 'e', 'أرتاح نفسيًا عندما أستشير المساعد الذكي في قرار شخصي.'], ['emo', 'b', 'أحكي للمساعد الذكي أمورًا لا أحكيها لأحد.'],
       ['ctrl', 'b', 'أستخدم هذه الأدوات أكثر مما كنت أنوي.'], ['ctrl', 'e', 'أتضايق إذا كانت الخدمة متوقفة وعندي مهمة.']],
    amb: [['أعطاك المساعد الذكي إجابة تختلف عما في كتاب المقرر. الأرجح أن:', [['الكتاب قديم والمساعد أدق', 2], ['أحدهما خطأ ولا بد أن أتحقق', 0], ['المساعد أخطأ غالبًا', 0]]],
          ['طلب الأستاذ مقالًا قصيرًا برأيك الشخصي. أول ما تفكر فيه:', [['سأطلب مسودة ثم أعدلها', 2], ['سأكتب نقاطي أولًا', 0], ['سأقرأ عن الموضوع ثم أكتب', 0]]]],
    sjt: [['أمامك مسألة إحصاء جديدة والوقت متاح. أقرب تصرف:', [['أصورها وأطلب الحل مباشرة', 2], ['أحاول عشر دقائق ثم أطلب تلميحًا', 0], ['أحاول ثم أقارن حلي بحل الأداة', 0]]],
          ['انقطعت الخدمة الذكية ليلة تسليم واجب. غالبًا:', [['أشعر بالعجز وأؤجل حتى تعود', 2], ['أكمل بالمراجع والكتاب', 0], ['أطلب المساعدة من زميل', 0]]]],
    oth: [['من كل 10 طلاب، كم واحدًا تتوقع أنه يسلّم واجبًا كتبه الذكاء الاصطناعي دون أن يراجعه؟']],
    fc: [['أفضّل أن أحصل على الإجابة جاهزة وأوفر وقتي.', 'أفضّل أن أعمل على المهمة في الصباح الباكر.'], ['أشعر أن عقلي يعمل أفضل مع المساعد الذكي.', 'أشعر أن عقلي يعمل أفضل بعد النوم الجيد.']],
    ladder: ['أستخدم الأدوات الذكية للبحث السريع.', 'أستخدمها لكتابة أجزاء من واجباتي.', 'نادرًا ما أبدأ مهمة دون أن أسألها أولًا.', 'أشعر بأنني لا أستطيع الدراسة بدونها.', 'تراجعت قدرتي على التفكير وحدي بوضوح.'] },

  /* ===== M6: original construct under development (Mohamed Osman) ===== */
  { id: 'm6_avail', n: 'ضغط التواجد الرقمي الدائم', type: 'trait', dev: true, leak: ['تواجد رقمي', 'ضغط رقمي', 'الإتاحة'],
    def: 'شعور داخلي بالإلزام بأن يكون الفرد متاحًا ومستجيبًا عبر تطبيقات التواصل في كل وقت، يظهر في إلحاح الرد الفوري، والمراقبة المستمرة، والانشغال الذهني بالمجال الرقمي، والقلق من مؤشرات التواجد (تمت القراءة، آخر ظهور، متصل الآن)، وتوقع أن الآخرين ينتظرون الرد.',
    theory: 'يمتد المفهوم من «اليقظة الرقمية» (Reinecke et al., 2018) بأبعادها: البروز والاستجابية والمراقبة، ومن «الضغط عن بعد» (Barber & Santuzzi, 2015)، ويضيف بعدين ثقافيين: قلق مؤشرات التواجد في الثقافات الجمعية حيث يُفسَّر التأخر في الرد على أنه رسالة اجتماعية، والتوقعات المدركة للإتاحة من الأسرة والأصدقاء.',
    nomo: { pos: ['الخوف من فوات الأشياء (FOMO)', 'الاستخدام الإشكالي للهاتف', 'القلق الاجتماعي'], neg: ['الحدود الرقمية الصحية', 'جودة النوم'], none: ['عدد ساعات الاستخدام وحده (يجب أن يتميز عنه)'] },
    facets: [['react', 'إلحاح الرد الفوري'], ['monitor', 'المراقبة المستمرة'], ['salience', 'البروز الذهني'], ['signals', 'قلق مؤشرات التواجد'], ['expect', 'التوقعات المدركة للإتاحة']],
    d: [['react', 'b', 'أقطع ما أفعله لأرد على رسالة بمجرد وصولها.'], ['react', 'e', 'أشعر بثقل إذا تركت رسالة دون رد أكثر من ساعة.'], ['react', 'b', 'أرد على الرسائل عندما يناسبني الوقت دون ضيق.', 1],
       ['monitor', 'b', 'أفتح المحادثات لأتأكد أنه لم يصلني شيء جديد.'], ['monitor', 'b', 'أتابع المجموعات وأنا أتكلم مع شخص أمامي.'],
       ['salience', 'c', 'أفكر فيما يحدث في المحادثات حتى وأنا بعيد عن الهاتف.'], ['salience', 'c', 'أنسى الهاتف تمامًا عندما أكون مشغولًا بنشاط ممتع.', 1],
       ['signals', 'e', 'يضايقني أن يرى شخص رسالتي ولا يرد.'], ['signals', 'b', 'أخفي «آخر ظهور» أو علامة القراءة حتى لا يُطلب مني الرد.'], ['signals', 'c', 'أفسر تأخر الرد كأن الشخص غاضب أو غير مهتم.'],
       ['expect', 'c', 'أظن أن أسرتي أو أصدقائي ينتظرون مني ردًا سريعًا.'], ['expect', 'e', 'أشعر بالذنب إذا اكتشف أحدهم أنني كنت متصلًا ولم أرد عليه.'], ['expect', 'c', 'من حقي ألا أكون متاحًا في كل وقت، والناس يتفهمون ذلك.', 1]],
    amb: [['رأى صديقك رسالتك منذ ساعتين ولم يرد. الأرجح أن:', [['هو متضايق مني', 2], ['هو مشغول وسيرد لاحقًا', 0], ['نسي الرد', 0]]],
          ['ظهر لك أن قريبًا «متصل الآن» لكنه لم يرد عليك. أول ما تفكر فيه:', [['يتجاهلني', 2], ['ربما يرد على شيء عاجل', 0], ['قد يكون فتح التطبيق بلا قصد', 0]]]],
    sjt: [['في منتصف محاضرة مهمة اهتز الهاتف برسالة من مجموعة الأسرة. غالبًا:', [['أرد فورًا ولو بكلمة', 2], ['أنظر لأعرف المرسل ثم أكمل', 1], ['أقرؤها بعد المحاضرة', 0]]],
          ['قررت أن تبتعد عن الهاتف يومًا للمذاكرة. في منتصف اليوم غالبًا:', [['أفتحه لأطمئن أنه لا يوجد شيء مهم', 2], ['أكمل اليوم كما خططت', 0], ['أخبر المقربين مسبقًا ثم أبتعد مرتاحًا', 0]]]],
    oth: [['من كل 10 أشخاص في سنك، كم واحدًا تتوقع أنه يعتبر عدم الرد خلال ساعة تجاهلًا؟']],
    fc: [['أشعر أنني يجب أن أكون متاحًا لمن يحتاجني في أي وقت.', 'أشعر أنني يجب أن أكون منظمًا في كل أموري.'], ['أنتبه لعلامة القراءة في رسائلي.', 'أنتبه لنظافة مكان عملي.']],
    ladder: ['أرد على الرسائل في نفس اليوم.', 'أرد على معظمها خلال دقائق.', 'أقطع مذاكرتي أو عملي لأرد.', 'أتضايق بوضوح إذا تأخرت في الرد.', 'لا أستطيع الاستمتاع بأي نشاط دون متابعة الرسائل.'] },

  /* ===== attitudes (Thurstone / Likert / Osgood) ===== */
  { id: 'att_help', n: 'الاتجاه نحو طلب المساعدة النفسية', type: 'attitude', obj: 'زيارة أخصائي نفسي', leak: [],
    def: 'تقييم الفرد العام، المعرفي والوجداني، لطلب المساعدة من متخصص نفسي عند المرور بمشكلة نفسية.',
    theory: 'نماذج سلوك طلب المساعدة (Fischer & Farina, 1995) ووصمة الذات (Vogel et al., 2006)، مع أبعاد ثقافية مقترحة: التفسير الروحي للمشكلة، وسمعة الأسرة.',
    nomo: { pos: ['النية لطلب المساعدة', 'المعرفة بالصحة النفسية'], neg: ['وصمة الذات'], none: ['التحصيل الدراسي'] },
    facets: [['cog', 'الجانب المعرفي'], ['aff', 'الجانب الوجداني'], ['soc', 'الجانب الاجتماعي']],
    d: [['cog', 'c', 'الأخصائي النفسي يستطيع مساعدة الناس على حل مشكلات حقيقية.'], ['cog', 'c', 'زيارة الأخصائي مضيعة للوقت والمال.', 1], ['aff', 'e', 'سأشعر بالارتياح إذا تكلمت مع متخصص عن مشكلة تضايقني.'], ['aff', 'e', 'يحرجني أن يعرف أحد أنني ذهبت لأخصائي نفسي.', 1], ['soc', 'c', 'من يذهب لأخصائي نفسي شخص شجاع.'], ['soc', 'c', 'سمعة الأسرة قد تتأثر إذا عُرف أن أحد أفرادها يُعالج نفسيًا.', 1]],
    th: [[10.8, 'العلاج النفسي حق أساسي يجب أن يحصل عليه كل من يحتاجه.'], [9.6, 'لو مررت بأزمة نفسية سأبحث عن متخصص قبل أي شيء.'], [8.7, 'الأخصائي النفسي يقدم مساعدة لا يقدمها الأصدقاء.'], [7.5, 'طلب المساعدة النفسية مفيد في أغلب الحالات.'], [6.4, 'قد ألجأ لأخصائي إذا لم تنفع الحلول الأخرى.'], [5.5, 'لا أعرف إن كان العلاج النفسي يفيد أم لا.'], [4.6, 'أفضّل حل مشكلاتي داخل الأسرة.'], [3.5, 'كثير من الناس يذهبون لأخصائي بلا داعٍ.'], [2.4, 'المشكلات النفسية تحل بقوة الإرادة وحدها.'], [1.3, 'العلاج النفسي للمجانين فقط.']] },
  { id: 'att_elearn', n: 'الاتجاه نحو التعليم الإلكتروني', type: 'attitude', obj: 'التعليم الإلكتروني', leak: [],
    def: 'تقييم الطالب العام لاستخدام المنصات والمحاضرات الإلكترونية بديلًا أو مكملًا للتعليم الحضوري.',
    theory: 'نموذج قبول التكنولوجيا (Davis, 1989): المنفعة المدركة وسهولة الاستخدام المدركة.',
    nomo: { pos: ['الكفاءة الذاتية التكنولوجية'], neg: ['قلق التكنولوجيا'], none: ['الطول'] },
    facets: [['use', 'المنفعة المدركة'], ['ease', 'السهولة المدركة'], ['aff', 'الارتياح']],
    d: [['use', 'c', 'المحاضرات المسجلة تساعدني على فهم ما فاتني.'], ['use', 'c', 'أتعلم أقل عندما تكون المحاضرة عن بعد.', 1], ['ease', 'c', 'أستخدم منصات التعلم بسهولة.'], ['ease', 'e', 'أرتبك عند رفع الواجبات على المنصات.', 1], ['aff', 'e', 'أستمتع بالدراسة عبر الإنترنت.'], ['aff', 'e', 'أشعر بالعزلة في المقررات الإلكترونية.', 1]],
    th: [[10.6, 'التعليم الإلكتروني أفضل طريقة للتعلم في هذا العصر.'], [9.2, 'أفضّل أن تكون كل المحاضرات متاحة إلكترونيًا.'], [8.1, 'التعليم الإلكتروني يوفر وقتًا ومجهودًا كبيرًا.'], [6.9, 'التعليم الإلكتروني مفيد كمكمل للمحاضرات.'], [5.6, 'التعليم الإلكتروني له مزايا وعيوب متساوية.'], [4.4, 'التعليم الإلكتروني يناسب بعض المقررات فقط.'], [3.2, 'أفقد تركيزي بسرعة في المحاضرات الإلكترونية.'], [2.1, 'التعليم الإلكتروني لا يصلح للتعليم الجامعي.'], [1.2, 'التعليم الإلكتروني إهدار كامل للعملية التعليمية.']] }
  ];
  var BY = {}; L.forEach(function (c) { BY[c.id] = c; });
  var FN = function (c, f) { var x = c.facets.find(function (a) { return a[0] === f; }); return x ? x[1] : f; };
  var OTH_OPTS = ['0–2', '3–4', '5–6', '7–8', '9–10'];
  var FIT = { // method suitability by construct type
    trait: { likert: 1, rasch: 1, cat: 1, aig: 1, guttman: 1, sjt: 1, forced: 1, thurstone: 0, osgood: 0 },
    attitude: { likert: 1, thurstone: 1, osgood: 1, guttman: 0.5, aig: 1, rasch: 0.5, cat: 0.5, sjt: 0, forced: 0 } };
  function fit(m, c) {
    var v = (FIT[c.type] || {})[m]; if (v === 1) return null;
    if (c.type === 'trait' && (m === 'thurstone' || m === 'osgood')) return { block: true, msg: (m === 'thurstone' ? 'طريقة ثرستون' : 'التمايز الدلالي لأوزجود') + ' صُممت لقياس الاتجاهات نحو موضوع، وليس سمات مثل «' + c.n + '». اختر ليكرت أو راش أو الحكم الموقفي أو الاختيار القسري.' };
    if (c.type === 'attitude' && (m === 'sjt' || m === 'forced')) return { block: true, msg: 'الاتجاه نحو «' + c.obj + '» تقييم لموضوع، ولا يناسبه ' + (m === 'sjt' ? 'الحكم الموقفي' : 'الاختيار القسري') + '. اختر ثرستون أو ليكرت أو أوزجود.' };
    return { block: false, msg: 'هذه الطريقة تصلح بشروط لهذا النوع من المفاهيم؛ راجع التعليق في الملاحظات.' };
  }
  // interleave items across facets so every k covers all facets
  function spread(arr, key) { var g = {}, order = []; arr.forEach(function (x) { var f = key(x); if (!g[f]) { g[f] = []; order.push(f); } g[f].push(x); }); var out = [], i = 0, more = true; while (more) { more = false; order.forEach(function (f) { if (g[f][i]) { out.push(g[f][i]); more = true; } }); i++; } return out; }
  function directItems(c) { return spread(c.d, function (x) { return x[0]; }).map(function (x) { return { text: x[2], dim: FN(c, x[0]), facet: x[0], ind: x[1], reverse: !!x[3], kind: 'd' }; }); }
  function indirectItems(c) {
    var out = [];
    (c.amb || []).forEach(function (a) { out.push({ text: a[0], dim: 'غير مباشر: تفسير المواقف', options: a[1].map(function (o) { return o[0]; }), w: a[1].map(function (o) { return o[1]; }), kind: 'amb' }); });
    (c.sjt || []).forEach(function (a) { out.push({ text: a[0], dim: 'غير مباشر: موقف سلوكي', options: a[1].map(function (o) { return o[0]; }), w: a[1].map(function (o) { return o[1]; }), kind: 'sjt' }); });
    (c.oth || []).forEach(function (a) { out.push({ text: a[0], dim: 'غير مباشر: تقدير عن الآخرين', options: OTH_OPTS, w: [0, 0.5, 1, 1.5, 2], kind: 'oth' }); });
    (c.fc || []).forEach(function (a, i) { var flip = i % 2 === 1, o = flip ? [a[1], a[0]] : [a[0], a[1]]; out.push({ text: 'أي العبارتين تنطبق عليك أكثر؟', dim: 'غير مباشر: اختيار قسري', options: o, w: flip ? [0, 2] : [2, 0], kind: 'fc' }); });
    return out;
  }
  var KIND = { d: 'مباشر', amb: 'تفسير موقف غامض', sjt: 'موقف سلوكي', oth: 'تقدير عن الآخرين', fc: 'اختيار قسري', lad: 'سلّم تراكمي', th: 'عبارة ثرستون' };
  function build(m, cid, path, k) {
    var c = BY[cid]; if (!c) return null;
    var f = fit(m, c); if (f && f.block) return { items: [], warn: f.msg, c: c };
    var items = [], notes = [], instr = '';
    if (m === 'guttman') { items = c.ladder.map(function (t, i) { return { text: t, dim: 'تراكمي', value: i + 1, kind: 'lad' }; }); notes.push('سلّم تراكمي مكتوب من الأخف للأشد. تأكد بعد التطبيق من معامل إعادة الإنتاج ≥ 0.90.'); }
    else if (m === 'thurstone') { items = (c.th || []).map(function (x) { return { text: x[1], dim: c.n, value: x[0], kind: 'th' }; }); notes.push('القيم المقترحة للعبارات تقديرية، ولا بد من حسابها من جديد بتحكيم 20 إلى 50 محكّمًا (الوسيط والمدى الربيعي).'); }
    else if (m === 'osgood') { var P = [['جيد', 'سيئ', 'التقييم'], ['مفيد', 'ضار', 'التقييم'], ['مريح', 'مزعج', 'التقييم'], ['مهم', 'تافه', 'التقييم'], ['قوي', 'ضعيف', 'القوة'], ['عميق', 'سطحي', 'القوة'], ['صعب', 'سهل', 'القوة'], ['نشط', 'خامل', 'النشاط'], ['سريع', 'بطيء', 'النشاط'], ['حديث', 'تقليدي', 'النشاط']]; items = P.slice(0, Math.max(6, Math.min(P.length, k))).map(function (p) { return { text: p[0] + ' ـــ ' + p[1], dim: p[2], kind: 'd' }; }); instr = 'قيّم «' + c.obj + '» على كل زوج من الصفات بوضع علامة في الدرجة الأقرب لرأيك.'; }
    else if (m === 'sjt') { items = indirectItems(c).filter(function (x) { return x.kind === 'sjt' || x.kind === 'amb'; }); notes.push('أوزان البدائل (0 إلى 2) مقترحة من النظرية، وتُراجع بتحكيم خبراء أو بمفتاح تجريبي من عينة.'); }
    else if (m === 'forced') { items = indirectItems(c).filter(function (x) { return x.kind === 'fc'; }); notes.push('العبارتان في كل زوج متقاربتان في المرغوبية الاجتماعية. افحص ذلك بتقدير المحكّمين للمرغوبية قبل التطبيق.'); }
    else {
      var dir = directItems(c), ind = indirectItems(c);
      if (path === 'indirect') items = ind;
      else if (path === 'both') items = dir.slice(0, Math.max(3, k)).concat(ind);
      else items = dir.slice(0, Math.max(3, k));
      if (path !== 'indirect' && k > dir.length) notes.push('المخطط الحالي فيه ' + dir.length + ' عبارة مباشرة. أكمل الباقي بالذكاء الاصطناعي أو بنفسك، مع الالتزام بالجوانب نفسها.');
      if (m === 'rasch' || m === 'cat') notes.push('صعوبة العبارات تُقدّر من بيانات التجريب (راش)، ولا تُفترض مسبقًا.');
    }
    if (!instr) instr = path === 'indirect' || m === 'sjt' || m === 'forced' ? 'هذا الاستبيان عن المواقف اليومية وطريقة تفكير الناس فيها. لا توجد إجابات صحيحة أو خاطئة؛ اختر الأقرب لك.' : 'اقرأ كل عبارة واختر مدى انطباقها عليك خلال الشهر الماضي. لا توجد إجابات صحيحة أو خاطئة.';
    if (path !== 'direct' && m !== 'guttman' && m !== 'thurstone' && m !== 'osgood') notes.push('للعبارات غير المباشرة غرض غير ظاهر، فيلزم شرح ذلك للمفحوص بعد التطبيق، وأن تذكر الموافقة أن «بعض الأسئلة غرضها لا يظهر من أول نظرة».');
    if (c.dev) notes.push('مفهوم أصيل قيد التطوير في المختبر (محمد عثمان). العبارات مسودة أولى للتحكيم.');
    if (f && !f.block) notes.push(f.msg);
    return { items: items, instructions: instr, notes: notes.join(' '), c: c };
  }
  /* ---------- item quality engine ---------- */
  var ABS = /(دائمًا|دائما|أبدًا|ابدا|إطلاقًا|مطلقًا|كل الناس|لا أحد|جميع الناس)/;
  var JARGON = /(اكتئاب|اضطراب|وسواس|ذهان|صدمة|رهاب|هستيري|عصاب)/;
  var NEG = /(^|\s)(لا|لم|لن|ليس|ليست)\s/g;
  function norm(s) { return String(s).replace(/[ًٌٍَُِّْـ]/g, '').replace(/[أإآ]/g, 'ا').replace(/ة/g, 'ه').replace(/ى/g, 'ي'); }
  function quality(it, c) {
    var t = it.text || '', f = [], nt = norm(t);
    var leak = (c ? [c.n].concat(c.leak || []) : []).filter(Boolean);
    if (it.kind !== 'th' && leak.some(function (w) { return w.length > 2 && nt.indexOf(norm(w)) >= 0; })) f.push(['bad', 'العبارة تذكر اسم المفهوم، فتقيس تسمية الذات لا السمة']);
    var negs = (t.match(NEG) || []).length; if (negs >= 2) f.push(['bad', 'نفي مزدوج']); else if (it.reverse && /^(لا|لم|لن|ليس|نادرًا|نادرا)\s/.test(t.trim())) f.push(['warn', 'عكسية بالنفي: الأفضل صياغة الطرف المقابل']);
    if (ABS.test(t)) f.push(['warn', 'كلمة مطلقة تقلل التباين']);
    var t2 = t.replace(/\s(وأنا|وأنت|وهو|وهي|وأنتم)\s/g, ' '); if (/\s(و|أو)(أ|ي|ت|ن)[^\s]{2,}/.test(t2) && /(أ|ي|ت|ن)[^\s]{2,}.*\s(و|أو)(أ|ي|ت|ن)/.test(t2) && it.kind === 'd') f.push(['warn', 'قد تكون مزدوجة: فكرتان في عبارة']);
    var w = t.trim().split(/\s+/).length; if (w > 22) f.push(['warn', 'طويلة (' + w + ' كلمة)']);
    if (c && c.type === 'trait' && it.kind === 'd' && JARGON.test(t) && !leak.some(function (x) { return JARGON.test(x); })) f.push(['warn', 'مصطلح إكلينيكي يكشف الغرض']);
    if (it.kind === 'd' && /؟\s*$/.test(t)) f.push(['warn', 'صيغة سؤال في عبارة تقرير ذاتي']);
    return f;
  }
  function coverage(items, c) {
    var dir = items.filter(function (x) { return x.kind === 'd' && x.facet; }), kinds = ['c', 'e', 'b', 'p'];
    var rows = c.facets.map(function (fa) { var cells = kinds.map(function (k) { return dir.filter(function (x) { return x.facet === fa[0] && x.ind === k; }).length; }); return [fa[1]].concat(cells, [cells.reduce(function (a, b) { return a + b; }, 0)]); });
    var rev = dir.filter(function (x) { return x.reverse; }).length, pr = dir.length ? rev / dir.length : 0;
    var empty = rows.filter(function (r) { return r[5] === 0; }).map(function (r) { return r[0]; });
    return { head: ['الجانب', 'فكر', 'شعور', 'سلوك', 'جسد', 'المجموع'], rows: rows, rev: rev, n: dir.length, pr: pr, empty: empty };
  }
  /* ---------- scoring: direct %, indirect %, combined, gap ---------- */
  function score(items, answers, nOpts) { // answers[i] = chosen option index (0-based) for options items, or Likert index 0..nOpts-1 (agreement order high->low handled by caller)
    var dS = 0, dM = 0, iS = 0, iM = 0;
    items.forEach(function (it, i) { var a = answers[i]; if (a == null || a < 0) return;
      if (it.w) { iS += it.w[a]; iM += Math.max.apply(null, it.w); }
      else if (it.kind === 'd') { var v = it.reverse ? a : (nOpts - 1 - a); dS += v; dM += nOpts - 1; } });
    var d = dM ? dS / dM * 100 : null, ind = iM ? iS / iM * 100 : null;
    return { direct: d, indirect: ind, combined: d != null && ind != null ? (d + ind) / 2 : (d != null ? d : ind), gap: d != null && ind != null ? d - ind : null };
  }
  function keyRows(items) { return items.map(function (it, i) { return [i + 1, KIND[it.kind] || '', it.dim || '', it.w ? it.options.map(function (o, j) { return o + ' = ' + it.w[j]; }).join(' · ') : (it.value != null ? 'القيمة ' + it.value : (it.reverse ? 'معكوسة' : 'موجبة'))]; }); }
  return { list: L, get: function (id) { return BY[id]; }, build: build, fit: fit, quality: quality, coverage: coverage, score: score, keyRows: keyRows, KIND: KIND, K: K, directItems: directItems, indirectItems: indirectItems, FN: FN };
})();

/* ===== Scale builder: 12-step wizard, data collection, live analysis ===== */
const SCALES = (() => {
  const K = 'scaleslab-scales-v1';
  let all = LS.get(K, []);
  const Cl = { push: null };
  const list = () => all;
  const get = id => all.find(s => s.id === id);
  let t = null;
  function save(s) { s.updated = Date.now(); const i = all.findIndex(x => x.id === s.id); if (i < 0) all.unshift(s); else all[i] = s; LS.set(K, all); clearTimeout(t); if (Cl.push) t = setTimeout(() => Cl.push(s), 1200); }
  function remove(id) { all = all.filter(s => s.id !== id); LS.set(K, all); Cl.remove && Cl.remove(id); }
  function merge(remote) { remote.forEach(r => { const l = get(r.id); if (!l || (r.updated || 0) > (l.updated || 0)) { const i = all.findIndex(x => x.id === r.id); if (i < 0) all.push(r); else all[i] = r; } }); LS.set(K, all); }
  return { list, get, save, remove, merge, cloud: Cl };
})();

const BUILD = (() => {
  const PRESETS = {
    agree5: ['أوافق بشدة', 'أوافق', 'محايد', 'أرفض', 'أرفض بشدة'],
    apply5: ['تنطبق دائمًا', 'تنطبق غالبًا', 'تنطبق أحيانًا', 'تنطبق نادرًا', 'لا تنطبق أبدًا'],
    freq5: ['دائمًا', 'غالبًا', 'أحيانًا', 'نادرًا', 'أبدًا'],
    agree4: ['أوافق بشدة', 'أوافق', 'أرفض', 'أرفض بشدة'],
    three: ['نعم', 'أحيانًا', 'لا'],
    binary: ['نعم', 'لا'],
    agreeBin: ['أوافق', 'لا أوافق'],
    osgood: ['3', '2', '1', '0', '1', '2', '3'],
    apply4: ['تنطبق دائمًا', 'تنطبق غالبًا', 'تنطبق نادرًا', 'لا تنطبق أبدًا'],
    sjt4: ['أ', 'ب', 'ج', 'د'],
    forced2: ['العبارة الأولى', 'العبارة الثانية']
  };
  const METHODS = {
    thurstone: { t: 'ثرستون', en: 'Thurstone Equal-Appearing Intervals', yr: 1928, era: 'كلاسيكية', d: 'يكتب الباحث عبارات كثيرة تتدرج من الرفض الشديد إلى القبول الشديد، ويقدّرها محكمون على متصل من 1 إلى 11، وتُختار العبارات ذات الاتفاق العالي والقيم الموزعة. يختار المفحوص «أوافق» لما ينطبق، ودرجته وسيط قيم العبارات التي وافق عليها.', preset: 'agreeBin', n: k => [100, '50–100 محكم لتقدير قيم العبارات + 100 مفحوص على الأقل'] },
    likert: { t: 'ليكرت', en: 'Likert Summated Rating', yr: 1932, era: 'كلاسيكية (الأكثر استخدامًا)', d: 'عبارات موجبة وسالبة على تدريج متدرج، وتُجمع الدرجات بعد عكس السالبة. سهلة البناء والتصحيح، وتُحلل بألفا كرونباخ والتحليل العاملي.', preset: 'apply5', n: k => [Math.max(200, 10 * k), '5–10 أفراد لكل عبارة وبحد أدنى 200 للتحليل العاملي، و30–50 للتجربة الاستطلاعية'] },
    guttman: { t: 'جتمان', en: 'Guttman Cumulative Scaling', yr: 1944, era: 'كلاسيكية', d: 'عبارات مرتبة تصاعديًا في الشدة ويُجاب عنها بنعم/لا؛ من يوافق على عبارة قوية يوافق على ما قبلها. يُتحقق منها بمعامل إعادة الإنتاج (≥ 0.90).', preset: 'binary', n: k => [100, '100 فرد على الأقل لحساب معامل إعادة الإنتاج'] },
    osgood: { t: 'التمايز الدلالي', en: 'Osgood Semantic Differential', yr: 1957, era: 'كلاسيكية', d: 'يُقيَّم مفهوم واحد على أزواج من الصفات المتضادة (اكتبها: الصفة الموجبة | الصفة السالبة) بتدريج سباعي، وتمثل ثلاثة أبعاد: التقييم، والقوة، والنشاط.', preset: 'osgood', n: k => [100, '100 فرد على الأقل، ويزيد للتحليل العاملي'] },
    rasch: { t: 'راش / IRT', en: 'Rasch Model / Item Response Theory', yr: 1960, era: 'حديثة', d: 'عبارات تغطي مستويات مختلفة من شدة السمة (سهلة وصعبة التأييد)، ويُقدَّر موقع الفرد والعبارة على متصل واحد باللوجيت مستقلًا عن العينة، مع فحص المطابقة والأداء التفاضلي DIF.', preset: 'apply4', n: k => [250, '150–250 فردًا لتقديرات مستقرة (Linacre)، و500 لنماذج IRT متعددة المعالم'] },
    cat: { t: 'الاختبار التكيفي', en: 'Computerized Adaptive Testing (CAT)', yr: 1980, era: 'حديثة', d: 'بنك عبارات معايَر بنظرية الاستجابة للمفردة؛ يختار الحاسوب العبارة التالية حسب إجابات الفرد فيقيس بدقة بعدد أقل، ويتوقف عند وصول الخطأ المعياري لحد معين (مثل 0.30). هنا تبني البنك وتعايره، ويُعرض ترتيب العبارات من الأسهل للأصعب.', preset: 'apply4', n: k => [500, '500 فرد على الأقل لمعايرة بنك العبارات'] },
    sjt: { t: 'الحكم الموقفي', en: 'Situational Judgment Tests (SJT)', yr: 1990, era: 'حديثة', d: 'يُعرض موقف واقعي وأربعة تصرفات، ويختار المفحوص أفضلها. تقيس الكفايات السلوكية وتقل فيها المرغوبية الاجتماعية، والتصحيح بمفتاح يحدده الخبراء (1 للتصرف الأفضل و0 لغيره). اكتب العبارة: الموقف | أ | ب | ج | د.', preset: 'sjt4', n: k => [200, '200 فرد على الأقل + 5–10 خبراء لبناء مفتاح التصحيح'] },
    forced: { t: 'الاختيار القسري', en: 'Forced-Choice / Thurstonian IRT', yr: 2011, era: 'حديثة', d: 'عبارتان متساويتان في المرغوبية الاجتماعية من بعدين مختلفين، ويختار المفحوص الأكثر انطباقًا. تقلل التزييف، لكن درجاتها «إبسَاتية» (نسبية داخل الفرد) فلا تصلح معها ألفا والتحليل العاملي التقليدي؛ وتُحلل حديثًا بنموذج Thurstonian IRT (Brown & Maydeu-Olivares, 2011). اكتب العبارة: العبارة الأولى | العبارة الثانية؛ والعبارة الأولى تُحسب لبعد العبارة.', preset: 'forced2', n: k => [300, '300 فرد على الأقل لنموذج Thurstonian IRT'] },
    aig: { t: 'توليد بالذكاء الاصطناعي', en: 'AI-assisted Automatic Item Generation', yr: 2020, era: 'الأحدث', d: 'تُستخدم النماذج اللغوية لتوليد عبارات كثيرة من تعريف السمة، ثم يراجعها الخبراء ويحكّمونها، وتُفحص إحصائيًا (التحليل العاملي، راش). يوفر الوقت لكنه لا يغني عن التحكيم والتجريب. انسخ «طلب التوليد» من خطوة فيض الفقرات.', preset: 'agree5', n: k => [Math.max(200, 10 * k), 'مثل ليكرت: 10 أفراد لكل عبارة وبحد أدنى 200، بعد التحكيم'] }
  };
  const MORDER = ['thurstone', 'likert', 'guttman', 'osgood', 'rasch', 'cat', 'sjt', 'forced', 'aig'];
  const parts = t => (t || '').split('|').map(x => x.trim());
  /* ready-made templates (adapted from the statistics lab builder) */
  function tplItems(m, c, dims, k) {
    const q = `«${c}»`;
    const pos = [`أشعر بـ${q} في كثير من المواقف اليومية`, `يلاحظ المحيطون بي ${q} في تصرفاتي`, `أفكر في الأمور المرتبطة بـ${q} أكثر من غيري`, `يظهر ${q} لدي بوضوح عند التعامل مع الآخرين`, `يصعب عليّ التحكم في ${q} عندما أكون تحت ضغط`, `يؤثر ${q} في أدائي الدراسي أو المهني`, `أشعر بـ${q} عند مواجهة مواقف جديدة`, `أتصرف بطريقة تعكس ${q} حتى في المواقف البسيطة`, `تتكرر لدي مواقف يظهر فيها ${q}`, `أجد نفسي مدفوعًا بـ${q} في قراراتي`, `يستمر ${q} معي لفترات طويلة`, `أتذكر مواقف كثيرة شعرت فيها بـ${q}`];
    const neg = [`نادرًا ما أشعر بـ${q}`, `يمر يومي دون أن أفكر في ${q}`, `أتعامل مع الآخرين بعيدًا عن ${q}`, `أتخذ قراراتي بهدوء بعيدًا عن ${q}`];
    const D = dims.length ? dims : [null], out = [];
    const add = (text, i, extra = {}) => out.push(Object.assign({ text, dim: D[i % D.length], dir: 1 }, extra));
    if (m === 'osgood') { const P = [['جيد', 'سيئ'], ['مفيد', 'ضار'], ['ممتع', 'ممل'], ['مهم', 'تافه'], ['مريح', 'مزعج'], ['قوي', 'ضعيف'], ['عميق', 'سطحي'], ['كبير', 'صغير'], ['نشط', 'خامل'], ['سريع', 'بطيء'], ['حاد', 'هادئ']]; for (let i = 0; i < k; i++) { const p = P[i % P.length]; add(`${p[0]} | ${p[1]}`, i); } return out; }
    if (m === 'guttman') { const L = ['أحيانًا', 'مرة في الأسبوع على الأقل', 'عدة مرات في الأسبوع', 'يوميًا تقريبًا', 'لدرجة أنه يؤثر في بعض مهامي', 'لدرجة أنه يمنعني من بعض الأنشطة', 'لدرجة أنه يمنعني من معظم الأنشطة', 'لدرجة أنني أحتاج مساعدة متخصصة']; for (let i = 0; i < Math.min(k, L.length); i++) add(`أشعر بـ${q} ${L[i]}`, 0); return out; }
    if (m === 'thurstone') { const V = [[11, `${q} أمر ضروري ولا يمكن الاستغناء عنه`], [10, `${q} من أهم ما يميز الإنسان الناجح`], [9, `أرى أن ${q} له فوائد كثيرة`], [8, `${q} مفيد في أغلب المواقف`], [7, `${q} مقبول بشكل عام`], [6, `لا أمانع وجود ${q} في حياتي`], [5, `ليس لدي رأي واضح حول ${q}`], [4, `أتحفظ أحيانًا على ${q}`], [3, `${q} يسبب مشكلات أكثر مما يفيد`], [2, `أرفض ${q} في معظم المواقف`], [1, `${q} أمر ضار تمامًا ولا فائدة منه`]]; const kk = Math.min(k, V.length); for (let i = 0; i < kk; i++) { const j = Math.round(i * (V.length - 1) / Math.max(1, kk - 1)); add(V[j][1], 0, { sv: V[j][0] }); } return out; }
    if (m === 'sjt') { const C = ['في قاعة المحاضرات', 'أثناء العمل في مجموعة', 'مع أحد أفراد الأسرة', 'في مقابلة مهمة', 'مع صديق مقرب', 'قبل امتحان مهم', 'في موقف أمام الناس', 'عند تلقي نقد من شخص آخر']; for (let i = 0; i < k; i++) add(`${C[i % C.length]} واجهت موقفًا يتطلب ${q}. ماذا تفعل؟ | أتعامل مع الموقف مباشرة وبهدوء | أطلب المشورة من شخص أثق به | أؤجل التعامل مع الموقف | أنسحب من الموقف تمامًا`, i, { key: 0 }); return out; }
    if (m === 'forced') { for (let i = 0; i < k; i++) add(`${pos[i % pos.length]} | ${pos[(i + 5) % pos.length]}`, i); return out; }
    for (let i = 0; i < k; i++) { const rev = ['likert', 'aig'].includes(m) && i % 4 === 3; add(rev ? neg[Math.floor(i / 4) % neg.length] : pos[(i - Math.floor(i / 4)) % pos.length], i, { dir: rev ? -1 : 1 }); }
    return out;
  }
  function aiPrompt(s) {
    const m = METHODS[meth(s)], k = Math.round(s.target * s.mult) || 20;
    return `أنت خبير في بناء المقاييس النفسية باللغة العربية. اكتب ${k} عبارة لقياس «${s.construct || s.name}» لدى ${s.group || 'الفئة المستهدفة'} بطريقة ${m.t} (${m.en}).
التعريف النظري: ${s.defC || '—'}
الأبعاد: ${s.dims.map(d => d.name + (d.def ? ' (' + d.def + ')' : '') + ' ' + d.w + '%').join('، ')}
قواعد: فصحى واضحة، فكرة واحدة لكل عبارة، بلا نفي مزدوج، عبارات تقريرية بضمير المتكلم، مناسبة للفئة.${['likert', 'aig'].includes(meth(s)) ? ' اجعل نحو ربع العبارات سالبة وضع قبلها (−).' : ''}${meth(s) === 'sjt' ? ' اكتب كل عبارة: الموقف | تصرف أ | تصرف ب | تصرف ج | تصرف د، واجعل التصرف أ هو الأفضل.' : ''}${meth(s) === 'forced' ? ' اكتب كل عبارة زوجًا متساويًا في المرغوبية: العبارة الأولى | العبارة الثانية.' : ''}${meth(s) === 'osgood' ? ' اكتب كل عبارة زوج صفات متضادة: الصفة الموجبة | الصفة السالبة.' : ''}
اكتب كل بعد بعنوانه، ثم عباراته كل عبارة في سطر، بلا ترقيم.`;
  }
  const meth = s => s.method || 'likert';
  const STEPS = ['الوضوح المفاهيمي', 'المراجعة الأدبية', 'تحديد الأبعاد', 'جدول المواصفات', 'فيض الفقرات', 'جودة الصياغة', 'سلم الاستجابة', 'صدق المحتوى', 'التجريب الاستطلاعي', 'جمع البيانات وتحليل الفقرات', 'الثبات', 'الصدق العاملي والتقنين'];
  const blank = () => ({ id: uid('s'), name: 'مقياس جديد', construct: '', group: '', purpose: 'بحث', defC: '', defO: '', refs: [], dims: [{ id: uid('d'), name: 'البعد الأول', def: '', w: 100 }], target: 20, mult: 1.5, items: [], scale: { preset: 'apply5', labels: PRESETS.apply5.slice() }, instr: 'بين يديك مجموعة من العبارات التي تصف بعض المواقف والمشاعر. اقرأ كل عبارة بعناية، ثم اختر البديل الذي ينطبق عليك فعلًا. لا توجد إجابات صحيحة وأخرى خاطئة، وإجاباتك سرية ولن تُستخدم إلا لأغراض البحث العلمي.', demo: { gender: true, age: true }, judges: 10, pilotNotes: '', responses: [], versions: {}, step: 0, created: Date.now(), updated: Date.now() });
  function example() {
    const s = blank(); s.name = 'مقياس التسويف الأكاديمي (مثال)'; s.construct = 'التسويف الأكاديمي'; s.group = 'طلاب الجامعة';
    s.defC = 'ميل الطالب إلى تأجيل البدء في المهام الأكاديمية أو إتمامها رغم إدراكه لعواقب التأجيل، مع ما يصاحب ذلك من ضعف في إدارة الوقت وانفعالات سلبية.';
    s.defO = 'الدرجة الكلية التي يحصل عليها الطالب على المقياس المكوّن من (12) عبارة موزعة على ثلاثة أبعاد (التأجيل السلوكي، وإدارة الوقت، والقلق والتجنب)، وتدل الدرجة المرتفعة على تسويف أعلى.';
    s.refs = [{ t: 'مقياس التسويف لدى الطلاب (PASS)', a: 'Solomon & Rothblum (1984)', n: 'يقيس تكرار التسويف وأسبابه؛ غير مقنن على البيئة المصرية' }];
    const D = [['التأجيل السلوكي', 'تأخير البدء في المهام والمذاكرة', 35], ['إدارة الوقت', 'سوء تقدير الوقت والتخطيط', 30], ['القلق والتجنب', 'الانفعالات المصاحبة والهروب من المهمة', 35]].map(([n, d, w]) => ({ id: uid('d'), name: n, def: d, w }));
    s.dims = D; s.target = 12; s.mult = 1;
    const I = [['أؤجل المذاكرة حتى الليلة السابقة للامتحان', 0, 1], ['أبدأ في تكليفاتي فور استلامها', 0, -1], ['أنشغل بالهاتف عندما أجلس للمذاكرة', 0, 1], ['أترك المهام الصعبة إلى آخر لحظة', 0, 1], ['أسيء تقدير الوقت الذي تحتاجه المهمة', 1, 1], ['أضع جدولًا زمنيًا للمذاكرة وألتزم به', 1, -1], ['أفاجأ بقرب موعد التسليم', 1, 1], ['أنجز أكثر من مهمة في وقت واحد فيتعطل الجميع', 1, 1], ['أشعر بالضيق عندما أفكر في تكليف لم أبدأه', 2, 1], ['أهرب من المذاكرة بأنشطة أخرى', 2, 1], ['أواجه المهام المزعجة دون تأخير', 2, -1], ['أتجنب التفكير في المقررات الصعبة', 2, 1]];
    s.items = I.map(([t, d, dir]) => ({ id: uid('i'), text: t, dim: D[d].id, dir, src: 'استنباطي', ess: null, use: null, keep: true }));
    return s;
  }
  const kept = s => s.items.filter(i => i.keep !== false);
  const pts = s => s.scale.labels.length;
  const isBin = s => pts(s) === 2;
  const scoreOf = (s, it, raw) => { if (raw == null) return null; if (meth(s) === 'sjt') return (pts(s) - raw) === (it.key || 0) ? 1 : 0; if (meth(s) === 'forced') return raw; const hi = isBin(s) ? 1 : pts(s), lo = isBin(s) ? 0 : 1; return it.dir < 0 ? hi + lo - raw : raw; };

  /* ---- wording checks ---- */
  const NEG = /(^|\s)(لا|لم|لن|ليس|ليست|غير|عدم|بدون|دون)(\s|$)/g;
  function wordingWarnings(s, it) {
    const t = (it.text || '').trim(), w = [];
    if (!t) return [['bad', 'عبارة فارغة']];
    if (meth(s) === 'osgood') return t.includes('|') ? [] : [['warn', 'اكتبها: الصفة الموجبة | الصفة السالبة']];
    if (meth(s) === 'sjt') return parts(t).filter(Boolean).length === 1 + pts(s) ? [] : [['bad', `اكتبها: الموقف | ${pts(s)} تصرفات مفصولة بـ |`]];
    if (meth(s) === 'forced') return parts(t).filter(Boolean).length === 2 ? [] : [['bad', 'اكتبها: العبارة الأولى | العبارة الثانية']];
    const negs = (t.match(NEG) || []).length;
    if (negs >= 2) w.push(['bad', 'نفي مزدوج محتمل']);
    if (/[؟?]/.test(t)) w.push(['warn', 'صيغة سؤال؛ اجعلها عبارة تقريرية']);
    if (t.length > 90) w.push(['warn', 'طويلة (أكثر من 90 حرفًا)']);
    if (t.length < 12) w.push(['warn', 'قصيرة وغامضة']);
    if (/(^|\s)(دائما|دائمًا|أبدا|أبدًا|مطلقا|مطلقًا|إطلاقا|إطلاقًا|جميع|كل)(\s|$)/.test(t)) w.push(['warn', 'كلمة مطلقة (دائمًا/أبدًا/كل)']);
    if (/\sأو\s/.test(t) || (/\sو(?!ال)[ء-ي]{3,}/.test(t) && t.split(/\sو/).length > 2)) w.push(['warn', 'قد تحمل فكرتين']);
    const tw = new Set(t.split(/\s+/)); kept(s).forEach(o => { if (o === it || !o.text) return; const ow = new Set(o.text.trim().split(/\s+/)); const inter = [...tw].filter(x => ow.has(x)).length, uni = new Set([...tw, ...ow]).size; if (uni && inter / uni >= 0.7) w.push(['bad', 'مكررة تقريبًا مع عبارة أخرى']); });
    return w;
  }

  /* ---- step status ---- */
  function status(s) {
    const r = s.responses.filter(x => !x.pilot), k = kept(s);
    const dimCount = d => k.filter(i => i.dim === d.id).length;
    const sumW = s.dims.reduce((a, d) => a + (+d.w || 0), 0);
    const words = k.some(i => wordingWarnings(s, i).some(x => x[0] === 'bad'));
    let alpha = NaN; try { if (r.length >= 5 && k.length >= 2) alpha = analyze(s).total.alpha; } catch (e) { }
    return [
      !!(s.name && s.construct && s.defC && s.defO),
      s.refs.length >= 1,
      s.dims.length >= 1 && s.dims.every(d => d.name && d.def),
      Math.abs(sumW - 100) < 0.5 && s.target > 0,
      s.items.length >= Math.round(s.target * s.mult) && s.dims.every(d => s.items.filter(i => i.dim === d.id).length >= 3),
      k.length > 0 && !words,
      pts(s) >= 2,
      s.items.length > 0 && s.items.every(i => i.ess != null),
      !!s.pilotNotes || s.responses.some(x => x.pilot),
      r.length >= 30,
      alpha >= 0.7,
      r.length >= Math.max(100, 5 * k.length)
    ];
  }

  /* ---- analysis ---- */
  function matrix(s, items) {
    const rows = [], meta = [];
    s.responses.filter(x => !x.pilot).forEach(x => { const row = items.map(it => scoreOf(s, it, x.vals[it.id])); if (row.every(v => v != null && isFinite(v))) { rows.push(row); meta.push(x); } });
    return { rows, meta };
  }
  function analyze(s) {
    const items = kept(s), M = matrix(s, items), out = { items, n: M.rows.length, meta: M.meta };
    if (M.rows.length < 3 || items.length < 2) return out;
    out.rows = M.rows;
    out.total = PSY.itemStats(M.rows);
    out.split = items.length >= 2 ? PSY.splitHalf(M.rows) : null;
    if (isBin(s) || meth(s) === 'sjt') out.kr20 = PSY.kr20(M.rows);
    out.ext = M.rows.length >= 6 ? PSY.extremeGroups(M.rows) : null;
    if (meth(s) === 'guttman' && isBin(s)) {
      const k = items.length, p = items.map((_, j) => PSY.mean(M.rows.map(r => r[j]))), ord = p.map((v, j) => [v, j]).sort((a, b) => b[0] - a[0]).map(x => x[1]);
      let err = 0; M.rows.forEach(r => { const t = PSY.sum(r); ord.forEach((j, pos) => { const ideal = pos < t ? 1 : 0; if (r[j] !== ideal) err++; }); });
      const n = M.rows.length; out.guttman = { cr: 1 - err / (n * k), mmr: PSY.mean(p.map(v => Math.max(v, 1 - v))), order: ord, p };
    }
    if (meth(s) === 'thurstone') { out.thurstone = M.rows.map(r => { const v = items.map((it, j) => r[j] === 1 && it.sv != null ? it.sv : null).filter(x => x != null); return v.length ? PSY.median(v) : null; }); out.thMissing = items.filter(it => it.sv == null).length; }
    if (['rasch', 'cat', 'sjt'].includes(meth(s))) { const lo = isBin(s) || meth(s) === 'sjt' ? 0 : 1, hi = isBin(s) || meth(s) === 'sjt' ? 1 : pts(s); const d = items.map((_, j) => { const pp = Math.min(0.99, Math.max(0.01, (PSY.mean(M.rows.map(r => r[j])) - lo) / (hi - lo))); return Math.log((1 - pp) / pp); }); const m = PSY.mean(d); out.rasch = d.map(x => x - m); }
    out.dims = s.dims.map(d => { const idx = items.map((it, j) => it.dim === d.id ? j : -1).filter(j => j >= 0); if (idx.length < 2) return { d, idx, n: idx.length }; const sub = M.rows.map(r => idx.map(j => r[j])); const st = PSY.itemStats(sub); let om = null; if (idx.length >= 3 && M.rows.length > idx.length) { const e = PSY.efa(sub, 1); om = PSY.omega(e.loadings.map(l => Math.abs(l[0]))); } return { d, idx, n: idx.length, st, omega: om, totals: sub.map(r => PSY.sum(r)) }; });
    return out;
  }

  /* ---- codes ---- */
  const b64e = str => btoa(unescape(encodeURIComponent(str))).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
  const b64d = str => decodeURIComponent(escape(atob(str.replace(/-/g, '+').replace(/_/g, '/'))));
  const hash = str => { let h = 7; for (let i = 0; i < str.length; i++) h = (h * 31 + str.charCodeAt(i)) % 1000003; return h.toString(36); };
  function version(s) { const ids = kept(s).map(i => i.id), f = activeForm(s); const v = hash(ids.join(',') + (f ? '|' + f.id : '')); s.versions = s.versions || {}; s.versions[v] = ids; s.vform = s.vform || {}; s.vform[v] = f ? f.id : 'base'; return v; }
  /* audience versions: same items, adapted wording per age group / register; responses are tagged with their version */
  const activeForm = s => (s.forms || []).find(f => f.id === s.activeForm) || null;
  const formText = (s, it, f) => { f = f === undefined ? activeForm(s) : f; return f && f.texts && f.texts[it.id] ? f.texts[it.id] : it.text; };
  function formDef(s) { const v = version(s), f = activeForm(s); return { i: s.id, v, n: s.name, in: (f && f.instr) || s.instr, l: (f && f.labels && f.labels.length === s.scale.labels.length) ? f.labels : s.scale.labels, b: isBin(s) ? 1 : 0, d: s.demo, m: meth(s), q: kept(s).map(i => formText(s, i, f)) }; }
  const scaleCode = s => 'SC1.' + b64e(JSON.stringify(formDef(s)));
  function parseScaleCode(c) { c = c.trim(); if (!c.startsWith('SC1.')) throw new Error('الكود لا يبدأ بـ SC1. تأكد أنك نسخت كود المقياس كاملًا.'); const d = JSON.parse(b64d(c.slice(4))); if (!d.q || !d.l || !d.i) throw new Error('كود المقياس غير مكتمل.'); if (d.q.length > 200 || d.l.length > 7) throw new Error('كود غير صالح.'); return d; }
  function responseCode(def, vals, demo) { const p = { i: def.i, v: def.v, r: vals.join(''), g: demo.g || '', a: demo.a || '', t: Date.now() }; p.c = hash(p.i + p.v + p.r); return 'RS1.' + b64e(JSON.stringify(p)); }
  function parseResponseCode(c) { c = c.trim(); if (!c.startsWith('RS1.')) return null; try { const p = JSON.parse(b64d(c.slice(4))); if (hash(p.i + p.v + p.r) !== p.c) return { error: 'تالف' }; return p; } catch (e) { return { error: 'تالف' }; } }
  function addResponse(s, ver, values, extra) { const ids = (s.versions || {})[ver]; if (!ids) return 'version'; if (values.length !== ids.length) return 'length'; const vals = {}; ids.forEach((id, j) => vals[id] = values[j]); s.responses.push(Object.assign({ id: uid('r'), ts: Date.now(), vals, form: (s.vform || {})[ver] || 'base' }, extra)); return 'ok'; }

  /* ---- form renderer (shared by kiosk, preview, respond) ---- */
  function formHTML(def, opts = {}) {
    const L = def.l, k = L.length;
    return `<form id="fForm" class="stack" novalidate>
      <div class="card"><h3>${esc(def.n)}</h3><p style="margin:0">${esc(def.in)}</p></div>
      ${(def.d && (def.d.gender || def.d.age)) ? `<div class="card grid2">${def.d.gender ? `<div class="field"><span class="fl">النوع</span><div class="seg" role="radiogroup" aria-label="النوع"><button type="button" data-g="ذكر" aria-pressed="false">ذكر</button><button type="button" data-g="أنثى" aria-pressed="false">أنثى</button></div></div>` : ''}${def.d.age ? `<div class="field"><label class="fl" for="fAge">العمر</label><input type="number" id="fAge" min="10" max="99" inputmode="numeric"></div>` : ''}</div>` : ''}
      <div class="form-preview card" style="padding:0">${def.q.map((q, i) => { const og = def.m === 'osgood' && q.includes('|'), pr = og ? q.split('|').map(x => x.trim()) : null, ch = (def.m === 'sjt' || def.m === 'forced') && q.includes('|') ? q.split('|').map(x => x.trim()) : null, lead = ch ? (def.m === 'sjt' ? ch[0] : 'أي العبارتين تنطبق عليك أكثر؟') : null, LL = ch ? (def.m === 'sjt' ? ch.slice(1) : ch).slice(0, k) : L; return `<fieldset class="fq" style="border:none;margin:0"><legend class="qt"><span class="num muted">${i + 1}.</span> ${og ? `<b>${esc(pr[0])}</b> <span class="muted">↔</span> <b>${esc(pr[1])}</b>` : esc(ch ? lead : q)}</legend><div class="lik${ch ? ' opts' : ''}" style="--k:${ch ? 1 : k}">${LL.map((lab, j) => `<label><input type="radio" name="q${i}" value="${def.b ? (j === 0 ? 1 : 0) : k - j}" ${opts.preview ? 'disabled' : ''}><span>${ch && def.m === 'sjt' ? esc(L[j] + ') ') : ''}${esc(lab)}</span></label>`).join('')}</div></fieldset>`; }).join('')}</div>
      ${opts.preview ? '' : `<div class="notice warn" id="fMiss" hidden></div><button class="btn" type="submit">${ic('check')} إرسال الإجابات</button>`}
    </form>`;
  }
  function bindForm(root, def, onSubmit) {
    let g = '';
    $$('[data-g]', root).forEach(b => b.onclick = () => { g = b.dataset.g; $$('[data-g]', root).forEach(x => x.setAttribute('aria-pressed', x === b)); });
    $('#fForm', root).onsubmit = e => {
      e.preventDefault();
      const vals = def.q.map((_, i) => { const c = $(`input[name="q${i}"]:checked`, root); return c ? +c.value : null; });
      const miss = vals.map((v, i) => v == null ? i + 1 : 0).filter(Boolean), m = $('#fMiss', root);
      if (miss.length) { m.hidden = false; m.textContent = 'لم تُجب عن العبارات: ' + miss.join('، ') + '. أجب عنها ثم أرسل.'; const f = $(`input[name="q${miss[0] - 1}"]`, root); f && f.closest('fieldset').scrollIntoView({ block: 'center' }); return; }
      const age = $('#fAge', root); const a = age && age.value ? String(Math.max(10, Math.min(99, Math.round(toNum(age.value))))) : '';
      onSubmit(vals, { g, a });
    };
  }

  /* ---- views ---- */
  let openId = null;
  VIEWS.build = el => { if (openId && SCALES.get(openId)) return editor(el, SCALES.get(openId)); openId = null; listView(el); };
  function listView(el) {
    const L = SCALES.list();
    el.innerHTML = `<div class="section-title"><h2>ابنِ مقياسك</h2><div class="row"><button class="btn" id="bNew">${ic('plus')} مقياس جديد</button><button class="btn ghost" id="bEx">ابدأ من مثال</button><button class="btn ghost" id="bImp">${ic('upload')} استيراد ملف مقياس</button><input type="file" id="fImp" accept=".json,application/json" hidden></div></div>
      <div class="notice small"><b>كيف يعمل؟</b>يسير المعالج على القواعد الاثنتي عشرة لبناء المقياس: من تعريف السمة إلى التحكيم والتجريب، ثم جمع الاستجابات (على جهازك مباشرة، أو بكود ترسله للمستجيبين، أو من Google Forms)، ويحسب التحليل فور وصول البيانات.</div>

      ${L.length ? `<div class="grid2">${L.map(s => { const st = status(s), done = st.filter(Boolean).length, r = s.responses.filter(x => !x.pilot).length; return `<article class="card stack" style="gap:8px"><div class="row"><b class="grow" style="font-family:var(--f-display);font-size:17px">${esc(s.name)}</b><span class="chip ${done === 12 ? 'good' : 'blue'}">${done}/12 خطوة</span></div><div class="row small muted"><span>${kept(s).length} عبارة</span><span>${s.dims.length} أبعاد</span><span>${r} استجابة</span></div><div class="pbar good"><i style="width:${done / 12 * 100}%"></i></div><div class="row"><button class="btn sm" data-open="${s.id}">فتح</button><button class="btn sm ghost" data-dup="${s.id}">نسخة</button><button class="btn sm ghost" data-del="${s.id}" aria-label="حذف">${ic('trash')}</button></div></article>`; }).join('')}</div>` : `<div class="card empty"><p>لم تنشئ مقياسًا بعد.</p><p class="small">ابدأ بمقياس جديد، أو افتح مثال «التسويف الأكاديمي» لترى الخطوات مكتملة.</p></div>`}`;
    $('#bNew').onclick = () => modal(`<h3>مقياس جديد: اختر طريقة العمل</h3><div class="stack">
      <button class="linkcard" data-nw="quick" style="--c:var(--amber);--s:var(--amber-soft)"><span class="ic">${ic('bolt')}</span><span><b>البناء السريع (مثل مختبر الإحصاء)</b><span>اختر طريقة من 9 طرق مع عالِمها ونموذج حقيقي لها، وابنِ العبارات، وأرسل المقياس بلينك أو Google Forms، واجمع الردود في Google Sheet.</span></span></button>
      <button class="linkcard" data-nw="full"><span class="ic">${ic('build')}</span><span><b>المعالج الكامل: 12 خطوة</b><span>من التعريف والأبعاد وجدول المواصفات إلى التحكيم والتجريب والتحليل والمعايير ودليل المقياس.</span></span></button>
      <button class="linkcard" data-nw="ex" style="--c:var(--good);--s:var(--good-soft)"><span class="ic">${ic('star')}</span><span><b>ابدأ من مثال جاهز</b><span>مقياس التسويف الأكاديمي (12 عبارة في 3 أبعاد) لتجرب كل الخطوات.</span></span></button>
      <button class="linkcard" data-nw="real" style="--c:var(--violet);--s:var(--violet-soft)"><span class="ic">${ic('book')}</span><span><b>نماذج حقيقية لكل طريقة ولكل عالم</b><span>ثرستون، وليكرت (روزنبرج)، وجتمان (بوجاردس)، وأوسجود، وراش، وPROMIS، والاختيار القسري.</span></span></button></div>
      <div class="row" style="margin-top:10px"><button class="btn ghost" data-close>إلغاء</button></div>`, (bg, close) => $$('[data-nw]', bg).forEach(x => x.onclick = () => { close(); const k = x.dataset.nw; if (k === 'quick') return go('qbuild'); if (k === 'real') return LEARN.open('real'); const s = k === 'ex' ? example() : blank(); SCALES.save(s); checkBadges(); openId = s.id; rerender(); }));
    $('#bEx').onclick = () => { const s = example(); SCALES.save(s); checkBadges(); openId = s.id; rerender(); };
    $('#bImp').onclick = () => $('#fImp').click();
    $('#fImp').onchange = e => { const f = e.target.files[0]; if (!f) return; const rd = new FileReader(); rd.onload = () => { try { const s = JSON.parse(rd.result); if (!s.items || !s.dims || !s.scale) throw 0; s.id = uid('s'); s.responses = Array.isArray(s.responses) ? s.responses : []; SCALES.save(s); toast('استُورد المقياس'); rerender(); } catch (err) { toast('الملف ليس ملف مقياس صالحًا'); } }; rd.readAsText(f); };
    $$('[data-open]', el).forEach(b => b.onclick = () => { openId = b.dataset.open; rerender(); });
    $$('[data-dup]', el).forEach(b => b.onclick = () => { const s = JSON.parse(JSON.stringify(SCALES.get(b.dataset.dup))); s.id = uid('s'); s.name += ' (نسخة)'; s.responses = []; s.versions = {}; SCALES.save(s); rerender(); });
    $$('[data-del]', el).forEach(b => b.onclick = () => confirmBox('سيُحذف المقياس وكل استجاباته من هذا الجهاز' + (window.ACCT && ACCT.signedIn() ? ' ومن حسابك' : '') + '.', 'حذف المقياس', () => { SCALES.remove(b.dataset.del); rerender(); }));
  }

  function editor(el, s) {
    const st = status(s), step = Math.min(11, Math.max(0, s.step || 0));
    el.innerHTML = `<button class="back" id="bBack">${ic('back')} كل المقاييس</button>
      <div class="card stack" style="gap:10px"><div class="row"><input type="text" id="sName" value="${esc(s.name)}" aria-label="اسم المقياس" style="font-family:var(--f-display);font-weight:600;font-size:18px" class="grow"><span class="chip ${st.filter(Boolean).length === 12 ? 'good' : 'blue'}">${st.filter(Boolean).length}/12</span><button class="btn sm ghost" id="bExp">${ic('download')} حفظ ملف المقياس</button></div>
      <nav class="steps-rail" aria-label="خطوات البناء">${STEPS.map((t, i) => `<button data-step="${i}" class="${st[i] ? 'done' : ''}" ${i === step ? 'aria-current="step"' : ''}><span class="n">${st[i] ? '✓' : i + 1}</span>${t}</button>`).join('')}</nav></div>
      <section class="card stack" id="stepBody"></section>
      <div class="row"><button class="btn ghost" id="bPrev" ${step === 0 ? 'disabled' : ''}>${ic('back')} السابق</button><span class="grow"></span><button class="btn" id="bNext" ${step === 11 ? 'disabled' : ''}>التالي</button></div>`;
    $('#bBack').onclick = () => { openId = null; rerender(); };
    $('#sName').onchange = e => { s.name = e.target.value.trim() || 'مقياس بلا اسم'; SCALES.save(s); };
    $('#bExp').onclick = () => downloadFile((s.name || 'scale').replace(/[^؀-ۿa-zA-Z0-9]+/g, '-') + '.json', JSON.stringify(s, null, 1), 'application/json');
    $$('[data-step]', el).forEach(b => b.onclick = () => { s.step = +b.dataset.step; SCALES.save(s); rerender(); });
    $('#bPrev').onclick = () => { s.step = step - 1; SCALES.save(s); rerender(); };
    $('#bNext').onclick = () => { s.step = step + 1; SCALES.save(s); rerender(); };
    const body = $('#stepBody');
    const intro = (i, txt) => `<div class="row"><span class="chip blue">القاعدة ${i + 1}</span><h3 class="grow" style="margin:0">${STEPS[i]}</h3>${st[i] ? '<span class="chip good">مكتملة</span>' : ''}</div><p class="muted small" style="margin:0">${txt}</p>`;
    const persist = () => { SCALES.save(s); };
    STEP[step](body, s, intro, persist, st);
  }

  const STEP = [];
  STEP[0] = (b, s, intro, persist) => {
    b.innerHTML = intro(0, 'حدد السمة التي تقيسها والفئة والغرض، ثم اكتب تعريفًا نظريًا يحدد حدودها وتعريفًا إجرائيًا يربطها بالمقياس.') + `
      <div class="grid2"><div class="field"><label class="fl" for="c1">السمة (البناء)</label><input type="text" id="c1" value="${esc(s.construct)}" placeholder="مثال: الاتزان الانفعالي"></div><div class="field"><label class="fl" for="c2">الفئة المستهدفة</label><input type="text" id="c2" value="${esc(s.group)}" placeholder="مثال: طلاب كلية التربية"></div></div>
      <div class="field"><label class="fl" for="c3">الغرض من المقياس</label><select id="c3">${['بحث', 'تشخيص وفرز', 'انتقاء', 'تقويم برنامج'].map(x => `<option ${s.purpose === x ? 'selected' : ''}>${x}</option>`).join('')}</select></div>
      <div class="field"><label class="fl" for="c4">التعريف النظري (المفاهيمي)</label><textarea id="c4" placeholder="السمة هي ...، وتتميز بـ ...، وتختلف عن ... في ...">${esc(s.defC)}</textarea></div>
      <div class="field"><label class="fl" for="c5">التعريف الإجرائي</label><textarea id="c5">${esc(s.defO)}</textarea><button class="btn sm ghost" id="tpl" type="button" style="justify-self:start">اقترح صياغة من الأبعاد</button></div>
      <div class="notice small"><b>تذكير</b>التعريف الإجرائي القوي يذكر عدد العبارات والأبعاد ودلالة الدرجة المرتفعة، لا «الدرجة على المقياس» فقط.</div>`;
    const bind = (id, f) => $('#' + id, b).onchange = e => { s[f] = e.target.value.trim(); persist(); };
    bind('c1', 'construct'); bind('c2', 'group'); bind('c3', 'purpose'); bind('c4', 'defC'); bind('c5', 'defO');
    $('#tpl', b).onclick = () => { const n = kept(s).length || s.target; s.defO = `الدرجة الكلية التي يحصل عليها ${s.group ? 'المستجيب من ' + s.group : 'المستجيب'} على مقياس ${s.construct || '...'} المكوّن من (${n}) عبارة موزعة على ${s.dims.length} ${s.dims.length > 2 ? 'أبعاد' : 'بعد'} (${s.dims.map(d => d.name).join('، ')})، وتدل الدرجة المرتفعة على مستوى أعلى من ${s.construct || 'السمة'}.`; $('#c5', b).value = s.defO; persist(); };
  };
  STEP[1] = (b, s, intro, persist) => {
    const draw = () => {
      b.innerHTML = intro(1, 'سجّل المقاييس والدراسات السابقة التي راجعتها، وما استخلصته منها من أبعاد، وسبب الحاجة إلى مقياس جديد.') +
        (s.refs.length ? tableHTML(['المقياس/الدراسة', 'المؤلف والسنة', 'ملاحظات وأبعاد', ''], s.refs.map((r, i) => [esc(r.t), esc(r.a), esc(r.n), `<button class="btn sm ghost" data-rm="${i}" aria-label="حذف">${ic('trash')}</button>`])) : '<p class="muted small">لم تضف مراجع بعد.</p>') +
        `<div class="grid3"><div class="field"><label class="fl" for="r1">المقياس أو الدراسة</label><input type="text" id="r1"></div><div class="field"><label class="fl" for="r2">المؤلف والسنة</label><input type="text" id="r2"></div><div class="field"><label class="fl" for="r3">الأبعاد والملاحظات</label><input type="text" id="r3"></div></div><button class="btn sm" id="radd" style="justify-self:start">${ic('plus')} إضافة</button>`;
      $('#radd', b).onclick = () => { const t = $('#r1', b).value.trim(); if (!t) return toast('اكتب اسم المقياس أو الدراسة'); s.refs.push({ t, a: $('#r2', b).value.trim(), n: $('#r3', b).value.trim() }); persist(); draw(); };
      $$('[data-rm]', b).forEach(x => x.onclick = () => { s.refs.splice(+x.dataset.rm, 1); persist(); draw(); });
    };
    draw();
  };
  STEP[2] = (b, s, intro, persist) => {
    const draw = () => {
      b.innerHTML = intro(2, 'قسّم السمة إلى أبعاد متجانسة داخليًا ومتمايزة فيما بينها، واكتب تعريفًا لكل بعد. يُنصح بثلاث عبارات على الأقل لكل بعد في الصورة النهائية.') +
        s.dims.map((d, i) => `<div class="grid2" style="border-top:1px solid var(--line);padding-top:10px"><div class="field"><label class="fl" for="dn${i}">اسم البعد ${i + 1}</label><input type="text" id="dn${i}" data-dn="${i}" value="${esc(d.name)}"></div><div class="field"><label class="fl" for="dd${i}">تعريفه</label><div class="row"><input type="text" id="dd${i}" data-dd="${i}" value="${esc(d.def)}" class="grow"><button class="btn sm ghost" data-drm="${i}" ${s.dims.length < 2 ? 'disabled' : ''} aria-label="حذف البعد">${ic('trash')}</button></div></div></div>`).join('') +
        `<button class="btn sm" id="dadd" style="justify-self:start">${ic('plus')} إضافة بعد</button>`;
      $$('[data-dn]', b).forEach(x => x.onchange = () => { s.dims[+x.dataset.dn].name = x.value.trim(); persist(); });
      $$('[data-dd]', b).forEach(x => x.onchange = () => { s.dims[+x.dataset.dd].def = x.value.trim(); persist(); });
      $$('[data-drm]', b).forEach(x => x.onclick = () => { const d = s.dims[+x.dataset.drm]; const moveTo = s.dims.find(o => o !== d); s.items.forEach(it => { if (it.dim === d.id) it.dim = moveTo.id; }); s.dims.splice(+x.dataset.drm, 1); persist(); draw(); });
      $('#dadd', b).onclick = () => { s.dims.push({ id: uid('d'), name: 'بعد جديد', def: '', w: 0 }); persist(); draw(); };
    };
    draw();
  };
  STEP[3] = (b, s, intro, persist) => {
    const draw = () => {
      const sumW = s.dims.reduce((a, d) => a + (+d.w || 0), 0), pool = Math.round(s.target * s.mult);
      b.innerHTML = intro(3, 'وزّع الوزن النسبي على الأبعاد (المجموع 100%)، وحدد عدد العبارات المستهدف في الصورة النهائية ومضاعف المجمّع الأولي.') +
        `<div class="grid3">${`<div class="field"><label class="fl" for="tg">العبارات في الصورة النهائية</label><input type="number" id="tg" value="${s.target}" min="2"></div><div class="field"><label class="fl" for="ml">مضاعف المجمّع الأولي (1–4)</label><input type="number" id="ml" value="${s.mult}" step="0.5" min="1" max="4"></div>`}<div class="tile ${Math.abs(sumW - 100) < 0.5 ? 'good' : 'bad'}"><small>مجموع الأوزان</small><b>${sumW}%</b></div></div>` +
        tableHTML(['البعد', 'الوزن %', 'العبارات النهائية', 'المجمّع الأولي', 'المكتوب حاليًا'], s.dims.map((d, i) => { const fin = Math.round(s.target * (+d.w || 0) / 100), have = s.items.filter(x => x.dim === d.id).length, need = Math.round(fin * s.mult); return [esc(d.name), `<input type="number" data-w="${i}" value="${d.w}" style="width:90px" aria-label="وزن ${esc(d.name)}">`, N(fin, 0), N(need, 0), { n: 1, v: `${have} ${have >= need ? '✓' : ''}` }]; })) +
        `<div class="row"><button class="btn sm ghost" id="eq">وزّع بالتساوي</button><span class="small muted">المجمّع الأولي الكلي: <b class="num">${pool}</b> عبارة</span></div>`;
      $('#tg', b).onchange = e => { s.target = Math.max(2, Math.round(toNum(e.target.value) || 2)); persist(); draw(); };
      $('#ml', b).onchange = e => { s.mult = Math.min(4, Math.max(1, toNum(e.target.value) || 1)); persist(); draw(); };
      $$('[data-w]', b).forEach(x => x.onchange = () => { s.dims[+x.dataset.w].w = Math.max(0, toNum(x.value) || 0); persist(); draw(); });
      $('#eq', b).onclick = () => { const n = s.dims.length, base = Math.floor(100 / n); s.dims.forEach((d, i) => d.w = base + (i < 100 - base * n ? 1 : 0)); persist(); draw(); };
    };
    draw();
  };
  STEP[4] = (b, s, intro, persist) => {
    const draw = () => {
      b.innerHTML = intro(4, 'اكتب عبارات أكثر من المطلوب (من مصادر استنباطية واستقرائية)، وحدد لكل عبارة بعدها واتجاهها: (+) موجبة، (−) سالبة تُعكس درجتها عند التصحيح.') +
        s.dims.map(d => { const its = s.items.filter(i => i.dim === d.id), need = Math.round(s.target * (+d.w || 0) / 100 * s.mult); return `<div class="stack" style="gap:4px"><div class="row"><b class="grow">${esc(d.name)}</b><span class="chip ${its.length >= need ? 'good' : 'warn'}">${its.length}/${need}</span></div>${its.map(it => itemRow(s, it)).join('')}<div class="row"><input type="text" placeholder="اكتب عبارة جديدة لهذا البعد ثم اضغط إضافة" data-new="${d.id}" class="grow" aria-label="عبارة جديدة في ${esc(d.name)}"><button class="btn sm" data-add="${d.id}">${ic('plus')} إضافة</button></div></div>`; }).join('<hr style="border:none;border-top:1px solid var(--line)">') +
        `<div class="row"><span class="chip">${METHODS[meth(s)].t}</span><button class="btn sm soft" id="tplGo">${ic('bolt')} مسودة من قوالب جاهزة</button><button class="btn sm lime" id="aiGen" hidden>${ic('bolt')} ولّد العبارات بحسابك في Claude</button><button class="btn sm ghost" id="aiCopy">${ic('copy')} انسخ طلب التوليد بالذكاء الاصطناعي</button></div><p class="small muted" style="margin:0">القوالب والطلب يتبعان طريقة البناء المختارة في الخطوة 7 (سلم الاستجابة). المسودة نقطة بداية فقط: عدّلها لتناسب سمتك، ثم حكّمها.</p>` +
        `<details class="steps"><summary>إضافة عبارات كثيرة دفعة واحدة</summary><p class="small muted">كل عبارة في سطر. اكتب (−) في بداية العبارة السالبة.</p><select id="bulkDim" aria-label="البعد">${s.dims.map(d => `<option value="${d.id}">${esc(d.name)}</option>`).join('')}</select><textarea id="bulk" style="margin-top:6px"></textarea><button class="btn sm" id="bulkAdd" style="margin-top:6px">إضافة الكل</button></details>`;
      bindItems(b, s, persist, draw);
      $$('[data-add]', b).forEach(x => x.onclick = () => { const inp = $(`[data-new="${x.dataset.add}"]`, b), t = inp.value.trim(); if (!t) return; s.items.push({ id: uid('i'), text: t, dim: x.dataset.add, dir: 1, src: 'استنباطي', ess: null, use: null, keep: true }); quest('items', 1); persist(); draw(); const n = $(`[data-new="${x.dataset.add}"]`, b); n && n.focus(); });
      $$('[data-new]', b).forEach(x => x.onkeydown = e => { if (e.key === 'Enter') { e.preventDefault(); $(`[data-add="${x.dataset.new}"]`, b).click(); } });
      $('#tplGo', b).onclick = () => { const k = Math.max(4, Math.round(s.target * s.mult) || 12); const its = tplItems(meth(s), s.construct || s.name || 'السمة', s.dims.map(d => d.id), k); its.forEach(x => s.items.push(Object.assign({ id: uid('i'), src: 'قالب جاهز', ess: null, use: null, keep: true }, x))); quest('items', its.length); persist(); draw(); toast(`أُضيفت ${its.length} عبارة من القوالب`); };
      if (typeof AI !== 'undefined') AI.available().then(ok => { const g = $('#aiGen', b); if (!ok || !g) return; g.hidden = false; g.onclick = async () => {
        g.disabled = true; const old = g.innerHTML; g.textContent = 'Claude يكتب العبارات…';
        try { const k = Math.max(6, Math.round(s.target * s.mult) || 12); const its = await AI.items(s, k); if (!its || !its.length) throw { code: 'invalid_json' };
          its.forEach(x => { const d = s.dims.find(dd => dd.name === String(x.dim || '').trim()) || s.dims[0]; s.items.push({ id: uid('i'), text: String(x.text || '').trim().slice(0, 300), dim: d.id, dir: x.neg ? -1 : 1, src: 'ذكاء اصطناعي', ess: null, use: null, keep: true }); });
          quest('items', its.length); persist(); draw(); toast(`أضاف Claude ${its.length} عبارة؛ راجعها وعدّلها قبل التحكيم`);
        } catch (e) { toast(AI.ERR[e && e.code] || 'تعذر التوليد الآن؛ استخدم القوالب أو انسخ الطلب'); g.disabled = false; g.innerHTML = old; } }; });
      $('#aiCopy', b).onclick = () => copyText(aiPrompt(s)).then(() => toast('نُسخ الطلب؛ الصقه في Claude ثم انسخ العبارات إلى «إضافة عبارات كثيرة»'));
      $('#bulkAdd', b).onclick = () => { const d = $('#bulkDim', b).value; const lines = $('#bulk', b).value.split(/\n/).map(x => x.trim()).filter(Boolean); lines.forEach(l => { const neg = /^\(?[−-]\)?\s*/.test(l); s.items.push({ id: uid('i'), text: l.replace(/^\(?[−-]\)?\s*/, ''), dim: d, dir: neg ? -1 : 1, src: 'استنباطي', ess: null, use: null, keep: true }); }); quest('items', lines.length); persist(); draw(); };
    };
    draw();
  };
  function itemRow(s, it, opts = {}) {
    const w = wordingWarnings(s, it), n = s.items.indexOf(it) + 1;
    return `<div class="itemrow" data-row="${it.id}"><span class="no">${n}</span><div><input type="text" value="${esc(it.text)}" data-txt="${it.id}" aria-label="نص العبارة ${n}" ${it.keep === false ? 'style="text-decoration:line-through;opacity:.6"' : ''}>${w.length ? `<div class="warns">${w.map(x => `<span class="chip ${x[0]}">${x[1]}</span>`).join('')}</div>` : (opts.showOk ? '<div class="warns"><span class="chip good">لا ملاحظات</span></div>' : '')}</div><div class="row" style="gap:4px;flex-wrap:nowrap">${meth(s) === 'sjt' ? `<select data-key="${it.id}" aria-label="مفتاح التصحيح للعبارة ${n}" title="التصرف الأفضل" style="width:70px">${s.scale.labels.map((l, j) => `<option value="${j}" ${(it.key || 0) === j ? 'selected' : ''}>✔ ${esc(l)}</option>`).join('')}</select>` : ''}${meth(s) === 'thurstone' ? `<input type="number" min="1" max="11" step="0.1" data-sv="${it.id}" value="${it.sv ?? ''}" placeholder="قيمة" style="width:76px" aria-label="القيمة المقياسية للعبارة ${n}">` : ''}<button class="dirbtn ${it.dir < 0 ? 'neg' : 'pos'}" data-dir="${it.id}" aria-label="الاتجاه: ${it.dir < 0 ? 'سالبة' : 'موجبة'}" title="${it.dir < 0 ? 'سالبة (تُعكس)' : 'موجبة'}">${it.dir < 0 ? '−' : '+'}</button><button class="ibtn" data-irm="${it.id}" aria-label="حذف العبارة">${ic('trash')}</button></div></div>`;
  }
  function bindItems(b, s, persist, draw) {
    $$('[data-txt]', b).forEach(x => x.onchange = () => { const it = s.items.find(i => i.id === x.dataset.txt); it.text = x.value.trim(); persist(); draw(); });
    $$('[data-dir]', b).forEach(x => x.onclick = () => { const it = s.items.find(i => i.id === x.dataset.dir); it.dir = it.dir < 0 ? 1 : -1; persist(); draw(); });
    $$('[data-irm]', b).forEach(x => x.onclick = () => { s.items = s.items.filter(i => i.id !== x.dataset.irm); persist(); draw(); });
    $$('[data-key]', b).forEach(x => x.onchange = () => { const it = s.items.find(i => i.id === x.dataset.key); it.key = +x.value; persist(); });
    $$('[data-sv]', b).forEach(x => x.onchange = () => { const it = s.items.find(i => i.id === x.dataset.sv); const v = toNum(x.value); it.sv = isFinite(v) && v >= 1 && v <= 11 ? v : null; persist(); });
  }
  STEP[5] = (b, s, intro, persist) => {
    const draw = () => {
      const k = kept(s), bad = k.filter(i => wordingWarnings(s, i).some(x => x[0] === 'bad')).length, warn = k.filter(i => wordingWarnings(s, i).length).length;
      const bal = s.dims.map(d => { const its = k.filter(i => i.dim === d.id), neg = its.filter(i => i.dir < 0).length; return [esc(d.name), N(its.length, 0), N(neg, 0), its.length ? (neg === 0 ? '<span class="chip warn">لا توجد عبارات سالبة</span>' : neg / its.length > 0.5 ? '<span class="chip warn">السالبة أكثر من النصف</span>' : '<span class="chip good">متوازن</span>') : '—']; });
      b.innerHTML = intro(5, 'يفحص المساعد العبارات آليًا: النفي المزدوج، والعبارات المركبة، والطويلة، والكلمات المطلقة، والتكرار. الفحص الآلي مساعد فقط ولا يغني عن حكم الخبراء.') +
        `<div class="kpis"><div class="tile ${bad ? 'bad' : 'good'}"><small>ملاحظات جوهرية</small><b>${bad}</b></div><div class="tile"><small>عبارات بها تنبيهات</small><b>${warn}</b></div><div class="tile"><small>إجمالي العبارات</small><b>${k.length}</b></div></div>` +
        tableHTML(['البعد', 'العبارات', 'السالبة', 'توازن الاتجاه'], bal) +
        `<div class="stack" style="gap:0">${k.map(it => itemRow(s, it, { showOk: true })).join('')}</div>` +
        `<details class="steps"><summary>قواعد الصياغة</summary><ul class="small"><li>فكرة واحدة في كل عبارة، ولغة مناسبة للفئة.</li><li>تجنب النفي المزدوج والعبارات العامة.</li><li>عبارة تقريرية بضمير المتكلم، لا سؤال.</li><li>تجنب الإيحاء بالإجابة المرغوبة اجتماعيًا.</li><li>عبارات سالبة محسوبة (نحو الربع إلى الثلث) لضبط الموافقة المنهجية، مع فحص أدائها لاحقًا.</li></ul></details><section class="card stack aud-box" id="audBox"></section>`;
      bindItems(b, s, persist, draw);
      audPanel($('#audBox', b), s, persist);
    };
    draw();
  };
  /* ---- audience versions panel ---- */
  function audPanel(el, s, persist) {
    if (!el || typeof AUD === 'undefined') return;
    s.forms = s.forms || []; let cur = s._editForm && s.forms.find(f => f.id === s._editForm) || null;
    const k = kept(s), lo = isBin(s) ? 0 : 1;
    const draw = () => {
      const pickList = `<div class="row" style="flex-wrap:wrap;gap:6px"><span class="chip ${!cur ? 'blue' : ''}">الأصلية</span>${s.forms.map(f => `<button class="chip ${cur && cur.id === f.id ? 'lime' : ''}" data-fe="${f.id}">${esc(f.label)}</button>`).join('')}</div>`;
      let body = '';
      if (!cur) {
        body = `<div class="grid2">${AUD.selectorsHTML({ p: 'teen', r: 'ammiya' }, 'au')}<div id="auTips"></div></div><div class="row"><button class="btn sm" id="auNew">${ic('plus')} أنشئ نسخة لهذه الفئة</button></div>`;
      } else {
        const p = AUD.prof(cur.p), sameN = (cur.labels || []).length === s.scale.labels.length;
        body = `<div class="row"><b class="grow">${esc(cur.label)}</b><button class="btn sm lime" id="auAI" hidden>${ic('bolt')} حسّن الصياغة بحساب Claude</button><button class="btn sm ghost" id="auCp">${ic('copy')} انسخ طلب التكييف</button><button class="btn sm danger" id="auDel">${ic('trash')} احذف النسخة</button><button class="btn sm ghost" id="auBack">رجوع</button></div>
          ${AUD.tipsHTML(cur.p)}
          ${p.neg === false && k.some(i => i.dir < 0) ? '<div class="notice warn small"><b>تنبيه لهذه الفئة</b>في المقياس عبارات سالبة، وهذه الفئة يصعب عليها النفي. الحل الآمن: استبدلها بعبارات مثبتة <b>في كل النسخ</b> (بما فيها الأصلية) مع تعديل اتجاه التصحيح، حتى تبقى النسخ متكافئة.</div>' : ''}
          ${p.pts.indexOf(s.scale.labels.length) < 0 ? `<div class="notice warn small"><b>عدد نقاط التدريج</b>المناسب لهذه الفئة ${p.pts.join(' أو ')} نقاط، والمقياس ${s.scale.labels.length}. لا يمكن تغيير عدد النقاط داخل نسخة (يكسر المقارنة)؛ إن احتجت ذلك فابنِ مقياسًا مستقلًا لهذه الفئة بمعايير خاصة.</div>` : ''}
          <div class="tw"><table><thead><tr><th>#</th><th>العبارة الأصلية</th><th>صياغة هذه النسخة</th><th>فحص السلامة</th></tr></thead><tbody>${k.map((it, n) => { const w = AUD.check(it.text, cur.texts[it.id] || '', cur.p, cur.r); return `<tr><td class="num">${n + 1}</td><td class="small">${esc(it.text)}</td><td><input type="text" data-at="${it.id}" value="${esc(cur.texts[it.id] || '')}" aria-label="صياغة العبارة ${n + 1}"></td><td>${w.length ? w.map(x => `<span class="chip ${x[0]}">${x[1]}</span>`).join(' ') : '<span class="chip good">متكافئة شكليًا</span>'}</td></tr>`; }).join('')}</tbody></table></div>
          <div class="field"><span class="fl">بدائل الاستجابة لهذه النسخة (نفس العدد والترتيب)</span><div class="grid3">${(sameN ? cur.labels : s.scale.labels).map((l, i) => `<input type="text" data-al="${i}" value="${esc(l)}" aria-label="البديل ${i + 1}">`).join('')}</div></div>
          <div class="field"><label class="fl" for="auIn">تعليمات هذه النسخة</label><textarea id="auIn">${esc(cur.instr || s.instr)}</textarea></div>`;
      }
      el.innerHTML = `<div class="row"><h4 class="grow" style="margin:0">${ic('user')} نسخ حسب الفئة العمرية ولغة الصياغة</h4><span class="chip">${s.forms.length} نسخة</span></div>
        <p class="small muted" style="margin:0">اكتب نسخة من العبارات تناسب جيلًا معينًا (ألفا، Z، الألفية، X، كبار السن) بالفصحى أو الفصحى الميسرة أو العامية. يتغير <b>اللفظ</b> فقط؛ ويبقى المعنى، واتجاه العبارة، والإطار الزمني، وعدد نقاط التدريج ثابتة. وبعد جمع البيانات يقارن المعالج النسخ بالأداء التفاضلي للعبارات ليتأكد أنها تقيس السمة نفسها.</p>${pickList}${body}`;
      $$('[data-fe]', el).forEach(x => x.onclick = () => { cur = s.forms.find(f => f.id === x.dataset.fe); s._editForm = cur.id; draw(); });
      const tipsBox = $('#auTips', el), pSel = $('#auP', el); if (tipsBox && pSel) { const upd = () => tipsBox.innerHTML = AUD.tipsHTML(pSel.value); pSel.onchange = upd; upd(); }
      const nw = $('#auNew', el); if (nw) nw.onclick = () => {
        const p = $('#auP', el).value, r = $('#auR', el).value, P = AUD.prof(p), R0 = AUD.REGISTERS.find(x => x.id === r);
        const f = { id: uid('f'), p, r, label: P.t.split(' ')[0] + ' · ' + R0.t, texts: {}, labels: AUD.labels(s.scale.labels, p, r), instr: AUD.instructions(p, r) || s.instr };
        k.forEach(it => f.texts[it.id] = AUD.adapt(it.text, p, r).text);
        s.forms.push(f); cur = f; s._editForm = f.id; persist(); draw(); toast('أُنشئت النسخة بتحويل آلي أولي؛ راجع كل عبارة');
      };
      if (!cur) return;
      $$('[data-at]', el).forEach(x => x.onchange = () => { cur.texts[x.dataset.at] = x.value.trim(); persist(); draw(); });
      $$('[data-al]', el).forEach(x => x.onchange = () => { if ((cur.labels || []).length !== s.scale.labels.length) cur.labels = s.scale.labels.slice(); cur.labels[+x.dataset.al] = x.value.trim() || s.scale.labels[+x.dataset.al]; persist(); });
      $('#auIn', el).onchange = e => { cur.instr = e.target.value.trim(); persist(); };
      $('#auBack', el).onclick = () => { cur = null; s._editForm = null; draw(); };
      $('#auDel', el).onclick = () => confirmBox('ستُحذف النسخة (الاستجابات المسجلة بها تبقى موسومة باسمها).', 'حذف النسخة', () => { s.forms = s.forms.filter(f => f !== cur); if (s.activeForm === cur.id) s.activeForm = null; cur = null; s._editForm = null; persist(); draw(); });
      const items = k.map(it => ({ id: it.id, text: it.text }));
      $('#auCp', el).onclick = () => copyText(AUD.prompt(items, cur.p, cur.r, s.construct || s.name)).then(() => toast('نُسخ الطلب؛ الصقه في Claude ثم انقل الصياغات إلى الجدول'));
      if (typeof AI !== 'undefined') AI.available().then(ok => { const g = $('#auAI', el); if (!ok || !g) return; g.hidden = false; g.onclick = async () => {
        g.disabled = true; g.textContent = 'Claude يكيّف الصياغة…';
        try { const sm = await window.claude.use('sample'); const r = await sm.json(AUD.prompt(items, cur.p, cur.r, s.construct || s.name)); (r.items || []).forEach(x => { if (cur.texts[x.id] !== undefined && x.text) cur.texts[x.id] = String(x.text).slice(0, 300); }); persist(); draw(); toast('راجع فحص السلامة لكل عبارة'); }
        catch (e) { toast((AI.ERR || {})[e && e.code] || 'تعذر التكييف الآن'); g.disabled = false; }
      }; });
    };
    draw();
  }
  function versionsCompare(s, A) {
    if (typeof AUD === 'undefined' || !A.rows || !A.meta) return '';
    const groups = {}; A.meta.forEach((m, i) => { const f = m.form || 'base'; (groups[f] = groups[f] || []).push(A.rows[i]); });
    const keys = Object.keys(groups).filter(k2 => groups[k2].length >= 10); if (keys.length < 2) return Object.keys(groups).length > 1 ? '<div class="notice small"><b>مقارنة النسخ</b>تحتاج 10 استجابات على الأقل لكل نسخة لبدء فحص تكافؤ النسخ.</div>' : '';
    const name = k2 => k2 === 'base' ? 'الأصلية' : ((s.forms || []).find(f => f.id === k2) || { label: 'نسخة محذوفة' }).label;
    const ref = keys.indexOf('base') >= 0 ? 'base' : keys.sort((a, c) => groups[c].length - groups[a].length)[0];
    const lo = isBin(s) ? 0 : 1, hi = isBin(s) ? 1 : pts(s), nI = A.items.length;
    const sum = keys.map(k2 => { const st = PSY.itemStats(groups[k2]); return [esc(name(k2)), N(groups[k2].length, 0), N(st.alpha, 3), N(PSY.mean(groups[k2].map(r => PSY.sum(r))))]; });
    let h = `<section class="stack"><h4 style="margin:0">${ic('shield')} تكافؤ النسخ: هل تقيس السمة نفسها؟</h4>${tableHTML(['النسخة', 'ن', 'ألفا', 'متوسط الدرجة الكلية'], sum)}`;
    keys.filter(k2 => k2 !== ref).forEach(k2 => {
      const d = AUD.stdDIF(groups[ref], groups[k2], nI, lo, hi);
      h += `<p class="small" style="margin:0"><b>${esc(name(k2))}</b> مقارنة بـ <b>${esc(name(ref))}</b> (مؤشر التوحيد لدورانز وكوليك، بعد مطابقة المستجيبين على الدرجة الكلية):</p>` + tableHTML(['#', 'العبارة', 'الفرق المعياري STD', 'الحكم'], A.items.map((it, j) => ({ cls: d[j].level === 'flag' ? 'flag' : '', cells: [N(s.items.indexOf(it) + 1, 0), esc(formText(s, it, (s.forms || []).find(f => f.id === k2) || null)), N(d[j].std), d[j].level === 'ok' ? '<span class="chip good">متكافئة</span>' : d[j].level === 'review' ? '<span class="chip warn">تُراجع</span>' : d[j].level === 'flag' ? '<span class="chip bad">تحيز في الصياغة</span>' : '—'] })));
    });
    return h + '<p class="small muted">الفرق المعياري أقل من 10% من مدى التدريج = متكافئة؛ 10–20% = تُراجع؛ أكثر من 20% = الصياغة غيّرت ما تقيسه العبارة فتُعدّل. هذا فحص أولي؛ للتأكيد النهائي استخدم التحليل العاملي التوكيدي متعدد المجموعات (التكافؤ القياسي) في برنامج AMOS بموقع الإحصاء.</p></section>';
  }
  function methodPicker(s) {
    const m = METHODS[meth(s)], k = kept(s).length || s.target, r = m.n(k);
    return `<div class="field"><span class="fl">طريقة البناء (مرتبة من الأقدم إلى الأحدث)</span><div class="seg" style="flex-wrap:wrap">${MORDER.map(x => `<button data-me="${x}" aria-pressed="${meth(s) === x}">${METHODS[x].t} <span class="small muted num">${METHODS[x].yr}</span></button>`).join('')}</div>
      <div class="card" style="margin-top:8px;padding:12px"><div class="row"><b class="grow">${m.t} <span class="small muted" dir="ltr">${m.en}</span></b><span class="chip ${m.era.startsWith('كلاس') ? 'blue' : 'good'}">${m.era} · ${m.yr}</span></div><p class="small" style="margin:6px 0">${m.d}</p><div class="notice small" style="margin:0"><b>العينة المقترحة</b>${r[1]} ← الحد الأدنى ≈ <span class="num">${r[0]}</span>.</div></div></div>`;
  }
  STEP[6] = (b, s, intro, persist) => {
    const names = { agree5: 'موافقة خماسي', apply5: 'انطباق خماسي', freq5: 'تكرار خماسي', agree4: 'موافقة رباعي (بلا محايد)', apply4: 'انطباق رباعي', three: 'ثلاثي', binary: 'ثنائي (نعم/لا)', agreeBin: 'ثنائي (أوافق/لا أوافق)', osgood: 'سباعي ثنائي القطب', sjt4: 'أربعة تصرفات (أ–د)', forced2: 'اختيار بين عبارتين' };
    const draw = () => {
      b.innerHTML = intro(6, 'اختر طريقة بناء المقياس، ثم تدريجًا موحدًا لكل العبارات. البديل الأول يأخذ أعلى درجة في العبارة الموجبة، وتُعكس الدرجات تلقائيًا في العبارات السالبة.') +
        methodPicker(s) +
        `<div class="seg">${Object.keys(PRESETS).map(k => `<button data-ps="${k}" aria-pressed="${s.scale.preset === k}">${names[k]}</button>`).join('')}</div>
        <div class="grid3">${s.scale.labels.map((l, i) => `<div class="field"><label class="fl" for="lb${i}">البديل ${i + 1} (درجة ${isBin(s) ? (i === 0 ? 1 : 0) : s.scale.labels.length - i})</label><input type="text" id="lb${i}" data-lb="${i}" value="${esc(l)}"></div>`).join('')}</div>
        <div class="field"><label class="fl" for="ins">تعليمات المقياس</label><textarea id="ins">${esc(s.instr)}</textarea></div>
        <div class="row"><label class="row small"><input type="checkbox" id="dg" ${s.demo.gender ? 'checked' : ''}> سؤال النوع</label><label class="row small"><input type="checkbox" id="da" ${s.demo.age ? 'checked' : ''}> سؤال العمر</label></div>
        ${s.responses.length ? '<div class="notice warn small">لديك استجابات مسجلة؛ تغيير عدد البدائل الآن يجعل الاستجابات السابقة غير قابلة للمقارنة.</div>' : ''}`;
      $$('[data-me]', b).forEach(x => x.onclick = () => { s.method = x.dataset.me; const p = METHODS[s.method].preset; s.scale.preset = p; s.scale.labels = PRESETS[p].slice(); persist(); draw(); });
      $$('[data-ps]', b).forEach(x => x.onclick = () => { s.scale.preset = x.dataset.ps; s.scale.labels = PRESETS[x.dataset.ps].slice(); persist(); draw(); });
      $$('[data-lb]', b).forEach(x => x.onchange = () => { s.scale.labels[+x.dataset.lb] = x.value.trim() || ('بديل ' + (+x.dataset.lb + 1)); s.scale.preset = 'custom'; persist(); });
      $('#ins', b).onchange = e => { s.instr = e.target.value.trim(); persist(); };
      $('#dg', b).onchange = e => { s.demo.gender = e.target.checked; persist(); };
      $('#da', b).onchange = e => { s.demo.age = e.target.checked; persist(); };
    };
    draw();
  };
  STEP[7] = (b, s, intro, persist) => {
    const draw = () => {
      const Nj = s.judges, crit = PSY.cvrCritical(Nj), law = PSY.lawsheCritical(Nj);
      b.innerHTML = intro(7, 'اعرض العبارات على المحكمين، ثم سجّل لكل عبارة عدد من قالوا إنها «أساسية». يحسب الموقع نسبة صدق المحتوى (لوشي) ويقارنها بالقيمة الحرجة لعدد المحكمين.') +
        `<div class="grid3"><div class="field"><label class="fl" for="nj">عدد المحكمين</label><input type="number" id="nj" value="${Nj}" min="3" max="60"></div><div class="tile hl"><small>القيمة الحرجة (أير وسكالي)</small><b>${fmt(crit.cvr)}</b></div><div class="tile"><small>جدول لوشي</small><b>${fmt(law)}</b></div></div>
        <div class="row"><button class="btn sm ghost" id="cpJ">${ic('copy')} نسخ استمارة التحكيم</button><button class="btn sm ghost" id="apply">طبّق القرار: استبعد المرفوضة</button><button class="btn sm ghost" id="restore">أعد كل العبارات</button></div>` +
        tableHTML(['#', 'العبارة', 'أساسية', 'CVR', 'القرار'], s.items.map((it, i) => { const c = it.ess == null ? null : PSY.cvr(it.ess, Nj), ok = c != null && c >= crit.cvr - 1e-9; return { cls: c != null && !ok ? 'flag' : '', cells: [N(i + 1, 0), esc(it.text) + (it.keep === false ? ' <span class="chip bad">مستبعدة</span>' : ''), `<input type="number" min="0" max="${Nj}" data-ess="${it.id}" value="${it.ess ?? ''}" style="width:80px" aria-label="عدد المحكمين الذين قالوا أساسية للعبارة ${i + 1}">`, N(c), c == null ? '<span class="chip">بانتظار التحكيم</span>' : ok ? '<span class="chip good">تُقبل</span>' : '<span class="chip bad">تُراجع أو تُحذف</span>'] }; })) +
        (s.items.every(i => i.ess != null) ? `<div class="notice good small">CVI (متوسط CVR للعبارات المقبولة) = <b class="num">${fmt(PSY.mean(s.items.filter(i => PSY.cvr(i.ess, Nj) >= crit.cvr - 1e-9).map(i => PSY.cvr(i.ess, Nj))), 3)}</b></div>` : '');
      $('#nj', b).onchange = e => { s.judges = Math.max(3, Math.min(60, Math.round(toNum(e.target.value) || 10))); persist(); draw(); };
      $$('[data-ess]', b).forEach(x => x.onchange = () => { const it = s.items.find(i => i.id === x.dataset.ess); const v = x.value === '' ? null : Math.max(0, Math.min(s.judges, Math.round(toNum(x.value)))); it.ess = v; persist(); draw(); });
      $('#apply', b).onclick = () => { let n = 0; s.items.forEach(it => { if (it.ess != null && PSY.cvr(it.ess, s.judges) < crit.cvr - 1e-9) { it.keep = false; n++; } }); persist(); toast(n ? 'استُبعدت ' + n + ' عبارة من النموذج' : 'لا عبارات دون القيمة الحرجة'); draw(); };
      $('#restore', b).onclick = () => { s.items.forEach(it => it.keep = true); persist(); draw(); };
      $('#cpJ', b).onclick = () => copyText(`استمارة تحكيم: ${s.name}\nالسمة: ${s.construct}\nالتعريف الإجرائي: ${s.defO}\n\nأمام كل عبارة: (أساسية / مفيدة لكنها غير أساسية / غير ضرورية) — وملاحظاتك على الصياغة.\n\n` + s.dims.map(d => `البعد: ${d.name} — ${d.def}\n` + s.items.filter(i => i.dim === d.id).map((it) => `${s.items.indexOf(it) + 1}. ${it.text}  [   ] أساسية  [   ] مفيدة  [   ] غير ضرورية`).join('\n')).join('\n\n'));
    };
    draw();
  };
  STEP[8] = (b, s, intro, persist) => {
    const def = formDef(s), pilot = s.responses.filter(x => x.pilot).length;
    b.innerHTML = intro(8, 'جرّب النموذج على عينة صغيرة (نحو 30) واطلب من بعضهم التفكير بصوت مرتفع لكشف العبارات الغامضة. استجابات التجريب تُحفظ منفصلة ولا تدخل التحليل.') +
      `<div class="row"><button class="btn" id="kp">${ic('play')} طبّق تجريبيًا على هذا الجهاز</button><span class="chip">${pilot} استجابة تجريبية</span>${pilot ? '<button class="btn sm ghost" id="clrP">حذف استجابات التجريب</button>' : ''}</div>
      <div class="field"><label class="fl" for="pn">ملاحظات التجريب (عبارات غامضة، زمن الإجابة، تعليقات المستجيبين)</label><textarea id="pn">${esc(s.pilotNotes)}</textarea></div>
      <h4 style="margin:6px 0 0">معاينة النموذج كما يراه المستجيب</h4>${kept(s).length ? formHTML(def, { preview: true }) : '<p class="muted">أضف عبارات أولًا.</p>'}`;
    $('#pn', b).onchange = e => { s.pilotNotes = e.target.value.trim(); persist(); };
    $('#kp', b).onclick = () => kiosk(s, true);
    const cp = $('#clrP', b); if (cp) cp.onclick = () => { s.responses = s.responses.filter(x => !x.pilot); persist(); rerender(); };
  };
  STEP[9] = (b, s, intro, persist) => {
    const A = analyze(s), r = s.responses.filter(x => !x.pilot);
    const live = window.ACCT && ACCT.forms && ACCT.forms.available();
    b.innerHTML = intro(9, 'اجمع الاستجابات بالطريقة المناسبة لك، وتُضاف كلها إلى المقياس نفسه ويُحدّث التحليل فورًا.') +
      ((s.forms || []).length ? `<div class="row small"><label for="actF"><b>النسخة المستخدمة في الجمع:</b></label><select id="actF"><option value="">الأصلية</option>${s.forms.map(f => `<option value="${f.id}" ${s.activeForm === f.id ? 'selected' : ''}>${esc(f.label)}</option>`).join('')}</select><span class="muted">كل استجابة تُوسم بنسختها، ثم تُقارن النسخ تلقائيًا.</span></div>` : '') +
      `<div class="kpis"><div class="tile hl"><small>استجابات صالحة للتحليل</small><b>${A.n}</b></div><div class="tile"><small>كل الاستجابات</small><b>${r.length}</b></div><div class="tile"><small>آخر استجابة</small><b style="font-size:15px">${r.length ? new Date(Math.max(...r.map(x => x.ts))).toLocaleString('ar-EG', { dateStyle: 'short', timeStyle: 'short' }) : '—'}</b></div></div>
      <div class="grid2">
        <button class="linkcard" id="cK"><span class="ic">${ic('play')}</span><span><b>جمع مباشر على هذا الجهاز</b><span>سلّم الجهاز للمستجيبين واحدًا تلو الآخر (في قاعة أو معمل)؛ يعود النموذج فارغًا بعد كل إرسال.</span></span></button>
        <button class="linkcard" id="cL" style="--c:var(--good);--s:var(--good-soft)"><span class="ic">${ic('link')}</span><span><b>نموذج أونلاين مباشر ${live ? '<span class="live-dot" aria-hidden="true"></span>' : ''}</b><span>${live ? 'رابط تشاركه، وتصل الإجابات إلى حسابك لحظيًا.' : 'يحتاج تفعيل الحسابات (Firebase) على نسخة الموقع المنشورة على GitHub.'}</span></span></button>
        <button class="linkcard" id="cC" style="--c:var(--violet);--s:var(--violet-soft)"><span class="ic">${ic('whats')}</span><span><b>كود المقياس + أكواد الاستجابة</b><span>أرسل كود المقياس (واتساب مثلًا)؛ يجيب المستجيب من صفحة «أجب عن مقياس» ويرسل لك كودًا قصيرًا تلصقه هنا.</span></span></button>
        <button class="linkcard" id="cG" style="--c:var(--amber);--s:var(--amber-soft)"><span class="ic">${ic('form')}</span><span><b>Google Forms أو Excel</b><span>انسخ العبارات إلى Google Forms، ثم استورد ملف CSV الناتج، أو صدّر بياناتك.</span></span></button>
      </div>
      <div id="colPanel"></div>
      <div class="row small"><button class="btn sm ghost" id="simGo">${ic('bolt')} ولّد عينة افتراضية للتدريب</button>${r.some(x => x.sim) ? `<button class="btn sm ghost" id="simDel">حذف البيانات الافتراضية (${r.filter(x => x.sim).length})</button>` : ''}<span class="muted">للتدرب على التحليل قبل جمع بيانات حقيقية؛ احذفها قبل التحليل النهائي.</span></div>
      ${A.total ? itemTable(s, A) + versionsCompare(s, A) + (typeof THESIS !== 'undefined' ? THESIS.panel(s, A) : '') : `<div class="notice">التحليل يبدأ عند وصول 3 استجابات مكتملة على الأقل (ويُنصح بـ 30 فأكثر).</div>`}
      ${r.length ? `<details class="steps"><summary>إدارة الاستجابات (${r.length})</summary><div class="row" style="margin-top:8px"><button class="btn sm ghost" id="csvOut">${ic('download')} تصدير CSV</button><button class="btn sm danger" id="delAll">${ic('trash')} حذف كل الاستجابات</button></div></details>` : ''}`;
    $('#cK', b).onclick = () => kiosk(s, false); if (typeof THESIS !== 'undefined' && A.total) THESIS.bind(b, s, A);
    const af = $('#actF', b); if (af) af.onchange = () => { s.activeForm = af.value || null; persist(); toast('ستُجمع الاستجابات القادمة بنسخة: ' + af.options[af.selectedIndex].text); };
    $('#simGo', b).onclick = () => { if (!kept(s).length) return toast('أضف عبارات أولًا'); const n = Math.min(300, Math.max(30, METHODS[meth(s)].n(kept(s).length)[0])); simulate(s, n); persist(); rerender(); toast(`وُلدت ${n} استجابة افتراضية`); };
    const sd0 = $('#simDel', b); if (sd0) sd0.onclick = () => { s.responses = s.responses.filter(x => !x.sim); persist(); rerender(); };
    $('#cC', b).onclick = () => codesPanel(s, $('#colPanel', b));
    $('#cG', b).onclick = () => gformsPanel(s, $('#colPanel', b));
    $('#cL', b).onclick = () => livePanel(s, $('#colPanel', b));
    const co = $('#csvOut', b); if (co) co.onclick = () => downloadFile((s.name || 'scale') + '.csv', toCSV(s), 'text/csv;charset=utf-8');
    const da = $('#delAll', b); if (da) da.onclick = () => confirmBox('ستُحذف كل الاستجابات غير التجريبية لهذا المقياس.', 'حذف الاستجابات', () => { s.responses = s.responses.filter(x => x.pilot); persist(); rerender(); });
  };
  function simulate(s, n) {
    const its = kept(s), k = pts(s), bin = isBin(s), ver = version(s), m = meth(s);
    const gauss = () => { let u = 0, v = 0; while (!u) u = Math.random(); while (!v) v = Math.random(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
    const lam = its.map((_, i) => 0.55 + ((i * 13) % 7) / 20), thr = its.map((it, i) => m === 'guttman' ? -1.2 + 2.4 * i / Math.max(1, its.length - 1) : m === 'thurstone' && it.sv != null ? 0 : ((i * 37) % 10 - 5) / 10);
    for (let p = 0; p < n; p++) {
      const th = gauss();
      const vals = its.map((it, j) => {
        let z = lam[j] * th + Math.sqrt(1 - lam[j] ** 2) * gauss() - thr[j];
        if (m === 'thurstone' && it.sv != null) z = 1.2 - Math.abs((th * 2 + 6) - it.sv) / 1.6 + gauss() * 0.5;
        if (m === 'sjt') { const good = z > -0.3; const key = it.key || 0; const j2 = good ? key : [0, 1, 2, 3].filter(x => x !== key)[Math.floor(Math.random() * 3)] % k; return k - j2; }
        if (bin) return z > 0 ? 1 : 0;
        const cut = [-1.3, -0.45, 0.35, 1.15, 1.9, 2.4].slice(0, k - 1); const sc = 1 + cut.filter(c => z > c).length; return it.dir < 0 ? k + 1 - sc : sc;
      });
      addResponse(s, ver, vals, { src: 'افتراضية', sim: true, g: Math.random() < 0.5 ? 'ذكر' : 'أنثى', a: String(18 + Math.floor(Math.random() * 8)) });
    }
  }
  function itemTable(s, A) {
    const t = A.total, e = A.ext;
    return `<section class="stack"><h4 style="margin:0">تحليل الفقرات (ن = ${A.n})</h4>${tableHTML(['#', 'العبارة', 'البعد', 'المتوسط', 'الانحراف', 'r المصحح', 'ألفا عند الحذف', e ? 't الطرفيتين' : '', 'الحكم'].filter(Boolean), A.items.map((it, j) => { const bad = t.rit[j] < 0.3 || (e && e.items[j].p >= 0.05); const d = s.dims.find(x => x.id === it.dim); return { cls: bad ? 'flag' : '', cells: [N(s.items.indexOf(it) + 1, 0), esc(it.text) + (it.dir < 0 ? ' <span class="chip">معكوسة</span>' : ''), esc(d ? d.name : ''), N(t.itemMean[j]), N(t.itemSd[j]), N(t.rit[j]), N(t.alphaIfDeleted[j], 3), ...(e ? [N(e.items[j].t)] : []), bad ? `<span class="chip bad">${t.rit[j] < 0.3 ? 'r < 0.30' : 'تمييز غير دال'}</span>` : '<span class="chip good">جيدة</span>'] }; }))}<p class="small muted">العبارة الضعيفة تُراجع صياغتها واتجاه تصحيحها قبل حذفها؛ وألفا عند الحذف وحدها لا تكفي للحكم (القاعدة 10).</p></section>`;
  }
  STEP[10] = (b, s, intro) => {
    const A = analyze(s);
    if (!A.total) { b.innerHTML = intro(10, 'يُحسب الثبات بعد جمع البيانات.') + '<div class="notice">اجمع 3 استجابات على الأقل في الخطوة 10.</div>'; return; }
    const [c, tx] = relVerdict(A.total.alpha);
    b.innerHTML = intro(10, 'ثبات المقياس الكلي والأبعاد: ألفا كرونباخ، وأوميغا (من تشبعات العامل الأول لكل بعد)، والتجزئة النصفية' + (isBin(s) ? '، وKR‑20 للفقرات الثنائية' : '') + '.') +
      `<div class="kpis"><div class="tile hl"><small>ألفا للمقياس الكلي</small><b>${fmt(A.total.alpha, 3)}</b></div>${A.split ? `<div class="tile"><small>سبيرمان–براون</small><b>${fmt(A.split.sb, 3)}</b></div><div class="tile"><small>جتمان</small><b>${fmt(A.split.guttman, 3)}</b></div>` : ''}${A.kr20 ? `<div class="tile"><small>KR‑20</small><b>${fmt(A.kr20.value, 3)}</b></div>` : ''}<div class="tile"><small>ن</small><b>${A.n}</b></div></div><div class="verdict ${c}">الثبات الكلي ${tx}</div>` +
      tableHTML(['البعد', 'العبارات', 'ألفا', 'أوميغا', 'AVE', 'الحكم'], A.dims.map(d => d.st ? [esc(d.d.name), N(d.n, 0), N(d.st.alpha, 3), N(d.omega ? d.omega.omega : NaN, 3), N(d.omega ? d.omega.ave : NaN, 3), `<span class="chip ${relVerdict(d.st.alpha)[0]}">${relVerdict(d.st.alpha)[1]}</span>`] : [esc(d.d.name), N(d.n, 0), '—', '—', '—', '<span class="chip warn">عبارات أقل من 2</span>'])) +
      (A.dims.filter(d => d.totals).length >= 2 ? `<h4 style="margin:0">ارتباطات الأبعاد ببعضها وبالدرجة الكلية</h4>${(() => { const D = A.dims.filter(d => d.totals), tot = A.total.total; return tableHTML(['', ...D.map(d => esc(d.d.name)), 'الكلي'], D.map(a => [esc(a.d.name), ...D.map(b2 => N(a === b2 ? 1 : PSY.pearson(a.totals, b2.totals))), N(PSY.pearson(a.totals, tot))])); })()}` : '') +
      (A.guttman ? `<div class="kpis"><div class="tile ${A.guttman.cr >= 0.9 ? 'good' : 'bad'}"><small>معامل إعادة الإنتاج CR</small><b>${fmt(A.guttman.cr, 3)}</b></div><div class="tile"><small>الحد الأدنى الهامشي MMR</small><b>${fmt(A.guttman.mmr, 3)}</b></div><div class="tile"><small>تحسن الإنتاج (CR − MMR)</small><b>${fmt(A.guttman.cr - A.guttman.mmr, 3)}</b></div></div><div class="verdict ${A.guttman.cr >= 0.9 ? 'good' : 'bad'}">${A.guttman.cr >= 0.9 ? 'السلم تراكمي مقبول (CR ≥ 0.90)' : 'السلم لا يحقق شرط التراكم (CR < 0.90)'}</div><p class="small muted">ترتيب العبارات من الأسهل تأييدًا إلى الأصعب: ${A.guttman.order.map(j => s.items.indexOf(A.items[j]) + 1).join(' ← ')}</p>` : '') +
      (A.thurstone ? `<div class="notice small"><b>درجات ثرستون</b>درجة كل مستجيب = وسيط القيم المقياسية للعبارات التي وافق عليها. <b>تنبيه:</b> ألفا وKR‑20 أعلاه محسوبتان من عدد مرات الموافقة، ولا تناسبان طريقة ثرستون؛ لأن عباراتها تمتد من القبول إلى الرفض، فيوافق الفرد على العبارات القريبة من موقفه فقط (نموذج النقطة المثالية)، وقد تخرج ألفا منخفضة أو سالبة. الأنسب هنا ثبات إعادة الاختبار لدرجات ثرستون. ${A.thMissing ? `<span style="color:var(--bad)">لم تُحدَّد القيمة المقياسية لـ ${A.thMissing} عبارة (الخطوة 5).</span>` : ''}</div><div class="kpis"><div class="tile"><small>متوسط درجات ثرستون</small><b>${fmt(PSY.mean(A.thurstone.filter(x => x != null)))}</b></div><div class="tile"><small>الانحراف المعياري</small><b>${fmt(PSY.sd(A.thurstone.filter(x => x != null)))}</b></div></div>` : '') +
      (A.rasch ? `<section class="stack"><h4 style="margin:0">صعوبة تأييد العبارات (تقدير تقريبي باللوجيت)</h4>${barSVG(A.rasch.map(x => x + 3), A.items.map(it => String(s.items.indexOf(it) + 1)), { label: 'صعوبة التأييد', d: 2, color: x => x - 3 > 0 ? 'var(--violet)' : 'var(--primary)' })}<p class="small muted">القيمة الموجبة = عبارة أصعب تأييدًا. هذا تقدير تقريبي من نسب الدرجات (مركّز عند الصفر)؛ للتقدير الدقيق بنموذج راش ومؤشرات المطابقة استخدم برنامج Winsteps أو حزمة TAM في R.</p></section>` : '') +
      `<div class="notice small">ثبات إعادة التطبيق يحتاج تطبيقًا ثانيًا على العينة نفسها بعد نحو أسبوعين؛ احسبه في «معمل ثبات إعادة الاختبار».</div>`;
  };
  STEP[11] = (b, s, intro) => {
    const A = analyze(s);
    if (!A.total || A.n < A.items.length + 2) { b.innerHTML = intro(11, 'يحتاج التحليل العاملي عددًا من الاستجابات أكبر من عدد العبارات (ويُنصح بـ 200 فأكثر).') + `<div class="notice">لديك ${A.n} استجابة صالحة و${A.items.length} عبارة.</div>`; return; }
    const R0 = PSY.corrMatrix(A.rows), km = PSY.kmo(R0), bt = PSY.bartlett(R0, A.n), nf = Math.max(1, Math.min(s.dims.length, A.items.length - 1)), fa = PSY.efa(A.rows, nf), L = fa.rotated;
    const tot = A.total.total, srt = PSY.sorted(tot), m = A.total.totalMean, sd = A.total.totalSd, uniq = [...new Set(srt)];
    const lvl = PSY.levels(isBin(s) ? 0 : 1, isBin(s) ? 1 : pts(s), 3), k = A.items.length;
    const desc = PSY.describe(tot);
    b.innerHTML = intro(11, 'الصدق العاملي (KMO وبارتلت والتشبعات بعد التدوير) ثم المعايير (الرتب المئينية والدرجات التائية) ودليل المقياس.') +
      `<div class="kpis"><div class="tile ${km && km.kmo >= 0.6 ? 'good' : 'bad'}"><small>KMO</small><b>${fmt(km ? km.kmo : NaN, 3)}</b></div><div class="tile ${bt.p < 0.05 ? 'good' : 'bad'}"><small>بارتلت</small><b>${fmtp(bt.p)}</b></div><div class="tile"><small>عوامل الجذر > 1</small><b>${fa.eigen.filter(x => x > 1).length}</b></div><div class="tile"><small>تباين مفسر (${nf} ${nf === 1 ? 'عامل' : 'عوامل'})</small><b>${fmt(PSY.sum(fa.varExp.slice(0, nf)), 1)}%</b></div></div>
      <section class="stack"><h4 style="margin:0">التشبعات بعد التدوير ومطابقتها للأبعاد المفترضة</h4>${tableHTML(['#', 'العبارة', 'البعد المفترض', ...L[0].map((_, f) => 'ع' + (f + 1))], A.items.map((it, i) => { const row = L[i], mx = Math.max(...row.map(Math.abs)); return { cls: mx < 0.4 ? 'flag' : '', cells: [N(s.items.indexOf(it) + 1, 0), esc(it.text), esc((s.dims.find(d => d.id === it.dim) || {}).name || ''), ...row.map(x => ({ n: 1, v: Math.abs(x) === mx && mx >= 0.4 ? `<b>${fmt(x)}</b>` : fmt(x) }))] }; }))}<p class="small muted">تحليل مكونات أساسية مع تدوير فاريماكس لعدد عوامل يساوي عدد الأبعاد. تحقق أن عبارات كل بعد تتشبع على العامل نفسه بقيمة ≥ 0.40، ثم أكّد البنية بالتحليل العاملي التوكيدي على عينة مستقلة.</p></section>
      <section class="stack"><h4 style="margin:0">الإحصاء الوصفي وتوزيع الدرجة الكلية</h4><div class="kpis"><div class="tile"><small>المتوسط</small><b>${fmt(desc.mean)}</b></div><div class="tile"><small>الانحراف</small><b>${fmt(desc.sd)}</b></div><div class="tile"><small>الوسيط</small><b>${fmt(desc.median)}</b></div><div class="tile"><small>الالتواء ÷ خطأه</small><b>${fmt(desc.skew / desc.seSkew)}</b></div><div class="tile"><small>التفلطح ÷ خطأه</small><b>${fmt(desc.kurt / desc.seKurt)}</b></div></div>${histSVG(tot)}</section>
      <section class="stack"><h4 style="margin:0">جدول المعايير</h4>${tableHTML(['الدرجة الخام', 'الرتبة المئينية', 'الدرجة التائية', 'المستوى (بمتوسط العبارة)'], uniq.map(x => { const im = x / k; const li = lvl.ranges.findIndex(r => im <= r[1] + 1e-9); return [N(x, 0), N(PSY.percentileRank(tot, x), 1), N(50 + 10 * (x - m) / sd, 1), ['منخفض', 'متوسط', 'مرتفع'][Math.max(0, li)]]; }))}<div class="row"><button class="btn sm ghost" id="nCsv">${ic('download')} تصدير المعايير CSV</button></div></section>
      <section class="card stack" style="background:var(--surface-2)"><h4 style="margin:0">${ic('book')} دليل المقياس وتقرير النتائج</h4><p class="small muted" style="margin:0">نص جاهز بأقسام دليل المقياس ونتائجك الفعلية، لتنقله إلى رسالتك وتراجعه.</p><div class="row"><button class="btn" id="man">${ic('copy')} نسخ دليل المقياس</button><button class="btn ghost" id="manDl">${ic('download')} تنزيل</button></div></section>`;
    $('#nCsv', b).onclick = () => downloadFile((s.name || 'scale') + '-norms.csv', 'raw,percentile,T\n' + uniq.map(x => [x, fmt(PSY.percentileRank(tot, x), 1), fmt(50 + 10 * (x - m) / sd, 1)].join(',')).join('\n'), 'text/csv;charset=utf-8');
    const text = () => manual(s, A, { km, bt, fa, nf, desc });
    $('#man', b).onclick = () => copyText(text());
    $('#manDl', b).onclick = () => downloadFile((s.name || 'scale') + '-manual.txt', text());
  };

  function manual(s, A, F) {
    const k = A.items.length, crit = PSY.cvrCritical(s.judges);
    const lines = [];
    lines.push(`دليل مقياس: ${s.name}`, '', `1) السمة: ${s.construct}`, `التعريف النظري: ${s.defC}`, `التعريف الإجرائي: ${s.defO}`, '', `2) الفئة المستهدفة: ${s.group}. الغرض: ${s.purpose}.`, '', '3) المراجعة الأدبية:', ...s.refs.map(r => `- ${r.t} (${r.a}): ${r.n}`), '', '4) الأبعاد وجدول المواصفات:');
    s.dims.forEach(d => lines.push(`- ${d.name} (${d.w}%): ${d.def} — ${A.items.filter(i => i.dim === d.id).length} عبارة`));
    lines.push('', `5) صيغة الاستجابة: ${s.scale.labels.join(' / ')}؛ وتُعكس درجات العبارات السالبة (${A.items.filter(i => i.dir < 0).map(i => s.items.indexOf(i) + 1).join('، ') || 'لا يوجد'}).`, '', `6) صدق المحتوى: عُرض المقياس على ${s.judges} محكمين، والقيمة الحرجة لنسبة صدق المحتوى ${fmt(crit.cvr)} (أير وسكالي، 2014)؛ واستُبعدت ${s.items.filter(i => i.keep === false).length} عبارة.`);
    lines.push('', `7) العينة: ${A.n} مستجيبًا.`, '', `8) الثبات: ألفا كرونباخ للمقياس الكلي = ${fmt(A.total.alpha, 3)}` + (A.split ? `، والتجزئة النصفية بمعادلة سبيرمان–براون = ${fmt(A.split.sb, 3)} وبمعادلة جتمان = ${fmt(A.split.guttman, 3)}` : '') + '.');
    A.dims.forEach(d => { if (d.st) lines.push(`- ${d.d.name}: α = ${fmt(d.st.alpha, 3)}` + (d.omega ? `، ω = ${fmt(d.omega.omega, 3)}` : '')); });
    lines.push('', `9) تحليل الفقرات: تراوحت معاملات الارتباط المصححة بين ${fmt(Math.min(...A.total.rit))} و${fmt(Math.max(...A.total.rit))}.`);
    if (F) lines.push('', `10) الصدق العاملي: KMO = ${fmt(F.km ? F.km.kmo : NaN, 3)}، واختبار بارتلت χ²(${F.bt.df}) = ${fmt(F.bt.chi, 2)}، p ${fmtp(F.bt.p).startsWith('<') ? fmtp(F.bt.p) : '= ' + fmtp(F.bt.p)}؛ وفسرت ${F.nf} ${F.nf === 1 ? 'عامل' : 'عوامل'} ${fmt(PSY.sum(F.fa.varExp.slice(0, F.nf)), 1)}% من التباين.`, '', `11) الإحصاء الوصفي للدرجة الكلية: المتوسط ${fmt(F.desc.mean)}، والانحراف المعياري ${fmt(F.desc.sd)}، والمدى ${F.desc.min}–${F.desc.max}.`);
    lines.push('', '12) التعليمات:', s.instr, '', '13) حدود الاستخدام: تُفسَّر الدرجات في ضوء خصائص عينة التقنين، ولا يُتخذ قرار فردي من درجة المقياس وحدها.');
    return lines.join('\n');
  }
  function toCSV(s) {
    const items = kept(s), head = ['id', 'time', 'source', 'gender', 'age', ...items.map((_, j) => 'item' + (j + 1)), 'total_scored', ...s.dims.map((_, i) => 'dim' + (i + 1))];
    const q = x => /[",\n]/.test(String(x)) ? '"' + String(x).replace(/"/g, '""') + '"' : x;
    const rows = s.responses.filter(x => !x.pilot).map(x => { const raw = items.map(it => x.vals[it.id] ?? ''); const sc = items.map(it => scoreOf(s, it, x.vals[it.id])); const ok = sc.every(v => v != null); return [x.id, new Date(x.ts).toISOString(), x.src || '', x.g || '', x.a || '', ...raw, ok ? PSY.sum(sc) : '', ...s.dims.map(d => ok ? PSY.sum(items.map((it, j) => it.dim === d.id ? sc[j] : 0)) : '')].map(q).join(','); });
    return head.join(',') + '\n' + rows.join('\n');
  }

  /* ---- collection panels ---- */
  function kiosk(s, pilot) {
    const def = formDef(s); if (!def.q.length) return toast('أضف عبارات أولًا');
    SCALES.save(s);
    let count = 0;
    const show = () => {
      LEARN.openOv(`<div class="qtop"><span class="chip ${pilot ? 'warn' : 'good'}">${pilot ? 'تطبيق تجريبي' : 'جمع مباشر'} · ${count} في هذه الجلسة</span><span class="grow"></span><button class="btn sm ghost" id="kEnd">إنهاء الجلسة</button></div><div id="kf" style="margin-top:12px"></div>`);
      $('#kf').innerHTML = formHTML(def); bindForm($('#kf'), def, (vals, demo) => { addResponse(s, def.v, vals, { src: pilot ? 'pilot' : 'kiosk', pilot: !!pilot, g: demo.g, a: demo.a }); count++; SCALES.save(s); checkBadges();
        $('#kf').innerHTML = `<div class="result"><div class="ring" style="width:84px;height:84px;display:grid;place-items:center;border-radius:50%;background:var(--good-soft);color:var(--good)">${ic('check')}</div><h2>شكرًا لك</h2><p class="muted">سُجلت إجاباتك. سلّم الجهاز للمستجيب التالي.</p><button class="btn" id="kNext">المستجيب التالي</button></div>`; $('#kNext').onclick = show; scrollTo(0, 0); $('#ov').scrollTop = 0; });
      $('#kEnd').onclick = () => confirmBox('إنهاء جلسة الجمع والعودة إلى المقياس؟', 'إنهاء', () => { LEARN.closeOv(); });
    };
    show();
  }
  function codesPanel(s, el) {
    const code = scaleCode(s); SCALES.save(s);
    const wa = 'https://wa.me/?text=' + encodeURIComponent(`أرجو الإجابة عن مقياس «${s.name}»:\n1) افتح موقع مختبر المقاييس النفسية واختر «أجب عن مقياس».\n2) الصق هذا الكود:\n${code}\n3) بعد الإجابة أرسل لي كود الاستجابة.`);
    el.innerHTML = `<section class="card stack"><h3>${ic('whats')} كود المقياس</h3><p class="small muted" style="margin:0">أرسل هذا الكود للمستجيبين. يفتحون صفحة «أجب عن مقياس» في الموقع، ويلصقونه، ويجيبون، ثم يرسلون لك «كود الاستجابة».</p><div class="codebox" id="scCode">${esc(code)}</div><div class="row"><button class="btn sm" id="cpSc">${ic('copy')} نسخ الكود</button><a class="btn sm ghost" href="${esc(wa)}" target="_blank" rel="noopener noreferrer">${ic('whats')} مشاركة عبر واتساب</a></div>
      <div class="field"><label class="fl" for="rc">الصق أكواد الاستجابة (كود في كل سطر، أو كلها معًا)</label><textarea id="rc" class="data" placeholder="RS1.…"></textarea></div><button class="btn" id="rcAdd" style="justify-self:start">${ic('plus')} إضافة الاستجابات</button><div id="rcMsg"></div></section>`;
    $('#cpSc', el).onclick = () => copyText(code);
    $('#rcAdd', el).onclick = () => {
      const codes = ($('#rc', el).value.match(/RS1\.[A-Za-z0-9_-]+/g) || []); let ok = 0, dup = 0, bad = 0, other = 0;
      const seen = new Set(s.responses.map(r => r.code).filter(Boolean));
      codes.forEach(c => { const p = parseResponseCode(c); if (!p || p.error) { bad++; return; } if (p.i !== s.id) { other++; return; } if (seen.has(c)) { dup++; return; } const vals = p.r.split('').map(Number); const st = addResponse(s, p.v, vals, { src: 'code', g: p.g, a: p.a, code: c }); if (st === 'ok') { ok++; seen.add(c); } else bad++; });
      SCALES.save(s); checkBadges();
      $('#rcMsg', el).innerHTML = `<div class="notice ${ok ? 'good' : 'warn'} small">أُضيفت ${ok} استجابة${dup ? `، وتُجوهلت ${dup} مكررة` : ''}${other ? `، و${other} لمقياس آخر` : ''}${bad ? `، و${bad} كود تالف أو لنسخة قديمة من العبارات` : ''}.</div>`;
      if (ok) setTimeout(rerender, 900);
    };
  }
  function gformsPanel(s, el) {
    const items = kept(s);
    const txt = `${s.name}\n\n${s.instr}\n\n` + items.map((it, i) => `${i + 1}. ${it.text}`).join('\n') + `\n\nالبدائل (لكل عبارة، سؤال «اختيار من متعدد» أو «شبكة اختيار متعدد»):\n` + s.scale.labels.join('\n');
    el.innerHTML = `<section class="card stack"><h3>${ic('form')} Google Forms</h3><ol class="small" style="margin:0"><li>أنشئ نموذجًا جديدًا في Google Forms، واستخدم سؤالًا من نوع «شبكة اختيار متعدد»: الصفوف = العبارات بالترتيب، والأعمدة = البدائل.</li><li>انسخ النص أدناه والصق العبارات والبدائل.</li><li>بعد جمع الإجابات: من تبويب «الردود» افتح جدول البيانات، ثم «ملف ← تنزيل ← CSV».</li><li>ارفع ملف CSV هنا وحدد العمود الذي تبدأ منه العبارة الأولى.</li></ol>
      <div class="row"><button class="btn sm" id="cpG">${ic('copy')} نسخ العبارات والبدائل</button></div>
      <hr style="border:none;border-top:1px solid var(--line);width:100%"><h3>${ic('upload')} استيراد CSV (Google Forms أو Excel)</h3><input type="file" id="csvIn" accept=".csv,text/csv" aria-label="ملف CSV"><div id="csvMap"></div></section>`;
    $('#cpG', el).onclick = () => copyText(txt);
    $('#csvIn', el).onchange = e => { const f = e.target.files[0]; if (!f) return; const rd = new FileReader(); rd.onload = () => csvMapper(s, parseCSV(rd.result), $('#csvMap', el)); rd.readAsText(f, 'utf-8'); };
  }
  function parseCSV(t) { const rows = []; let row = [], cur = '', q = false; t = t.replace(/^﻿/, ''); for (let i = 0; i < t.length; i++) { const ch = t[i]; if (q) { if (ch === '"') { if (t[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += ch; } else if (ch === '"') q = true; else if (ch === ',' || ch === ';' || ch === '\t') { row.push(cur); cur = ''; } else if (ch === '\n' || ch === '\r') { if (ch === '\r' && t[i + 1] === '\n') i++; row.push(cur); rows.push(row); row = []; cur = ''; } else cur += ch; } if (cur || row.length) { row.push(cur); rows.push(row); } return rows.filter(r => r.some(c => c.trim() !== '')); }
  function csvMapper(s, rows, el) {
    if (rows.length < 2) { el.innerHTML = '<div class="notice warn">الملف فارغ أو بلا صفوف بيانات.</div>'; return; }
    const head = rows[0], items = kept(s), L = s.scale.labels, kpt = L.length;
    const guess = Math.max(0, head.findIndex(h => items[0] && h.includes(items[0].text.slice(0, 12))));
    el.innerHTML = `<div class="grid3"><div class="field"><label class="fl" for="mStart">عمود العبارة الأولى</label><select id="mStart">${head.map((h, i) => `<option value="${i}" ${i === guess ? 'selected' : ''}>${i + 1}: ${esc(h.slice(0, 40))}</option>`).join('')}</select></div><div class="field"><label class="fl" for="mG">عمود النوع (اختياري)</label><select id="mG"><option value="-1">لا يوجد</option>${head.map((h, i) => `<option value="${i}">${i + 1}: ${esc(h.slice(0, 40))}</option>`).join('')}</select></div><div class="field"><label class="fl" for="mA">عمود العمر (اختياري)</label><select id="mA"><option value="-1">لا يوجد</option>${head.map((h, i) => `<option value="${i}">${i + 1}: ${esc(h.slice(0, 40))}</option>`).join('')}</select></div></div><div id="mPrev"></div><button class="btn" id="mGo">${ic('upload')} استيراد</button>`;
    const conv = cell => { const c = String(cell).trim(); if (c === '') return null; const n = toNum(c); if (isFinite(n)) return n; const j = L.findIndex(l => l.trim() === c); if (j < 0) return null; return kpt === 2 ? (j === 0 ? 1 : 0) : kpt - j; };
    const build = () => { const st = +$('#mStart', el).value, g = +$('#mG', el).value, a = +$('#mA', el).value; const out = []; rows.slice(1).forEach(r => { const vals = items.map((_, j) => conv(r[st + j])); if (vals.every(v => v != null && (kpt === 2 ? v === 0 || v === 1 : v >= 1 && v <= kpt))) out.push({ vals, g: g >= 0 ? String(r[g] || '').slice(0, 10) : '', a: a >= 0 ? String(r[a] || '').slice(0, 3) : '' }); }); return out; };
    const prev = () => { const o = build(); $('#mPrev', el).innerHTML = `<p class="small ${o.length ? '' : 'muted'}">صفوف صالحة: <b class="num">${o.length}</b> من ${rows.length - 1}. ${o.length < rows.length - 1 ? 'الصفوف الناقصة أو التي فيها قيم خارج التدريج تُتجاهل.' : ''}</p>`; };
    $$('select', el).forEach(x => x.onchange = prev); prev();
    $('#mGo', el).onclick = () => { const o = build(); if (!o.length) return toast('لا صفوف صالحة؛ تحقق من عمود البداية'); const v = version(s); o.forEach(x => addResponse(s, v, x.vals, { src: 'csv', g: x.g, a: x.a })); SCALES.save(s); checkBadges(); toast('استُوردت ' + o.length + ' استجابة'); rerender(); };
  }
  function livePanel(s, el) {
    if (!(window.ACCT && ACCT.forms && ACCT.forms.available())) {
      el.innerHTML = `<section class="card stack"><h3>${ic('link')} النموذج الأونلاين المباشر</h3>${ENV.claude ? '<p>داخل Claude لا يستطيع زوار الرابط العام حفظ إجاباتهم على الخادم؛ لذلك استخدم «كود المقياس + أكواد الاستجابة» أو «الجمع المباشر على هذا الجهاز». النموذج الأونلاين المباشر متاح في نسخة الموقع المنشورة على GitHub بعد تفعيل الحسابات.</p>' : (window.ACCT && ACCT.available() ? '<p>سجّل الدخول من صفحة «حسابي» لتنشر نموذجًا أونلاين تصل إجاباته إلى حسابك لحظيًا.</p><button class="btn" id="goAcct">تسجيل الدخول</button>' : '<p>هذه الميزة تحتاج ربط الموقع بمشروع Firebase مجاني (مرة واحدة). الخطوات في ملف <b>SETUP-FIREBASE.md</b> في المستودع. حتى ذلك الحين استخدم الكود أو الجمع المباشر.</p>')}</section>`;
      const g = $('#goAcct', el); if (g) g.onclick = () => go('account'); return;
    }
    ACCT.forms.panel(s, el, { formDef, addResponse, save: () => SCALES.save(s) });
  }

  /* ---- respond view (respondent side) ---- */
  VIEWS.respond = (el, formId) => {
    if (formId && window.ACCT && ACCT.forms && ACCT.forms.available()) return ACCT.forms.respond(el, formId, { formHTML, bindForm });
    el.innerHTML = `<div class="section-title"><h2>أجب عن مقياس</h2></div><div class="card stack"><p style="margin:0">الصق «كود المقياس» الذي أرسله لك الباحث، ثم أجب عن العبارات. ستحصل في النهاية على «كود استجابة» ترسله للباحث.</p><div class="field"><label class="fl" for="scIn">كود المقياس</label><textarea id="scIn" class="data" placeholder="SC1.…" style="min-height:90px"></textarea></div><button class="btn" id="scGo" style="justify-self:start">${ic('play')} افتح المقياس</button><div id="scErr"></div></div><div id="scForm"></div>`;
    $('#scGo', el).onclick = () => {
      let def; try { const m = $('#scIn', el).value.match(/SC1\.[A-Za-z0-9_-]+/); if (!m) throw new Error('لم أجد كود مقياس في النص الملصق.'); def = parseScaleCode(m[0]); } catch (e) { $('#scErr', el).innerHTML = `<div class="notice warn small">${esc(e.message || 'كود غير صالح')}</div>`; return; }
      $('#scErr', el).innerHTML = ''; const f = $('#scForm', el); f.innerHTML = formHTML(def);
      bindForm(f, def, (vals, demo) => {
        const code = responseCode(def, vals, demo);
        const wa = 'https://wa.me/?text=' + encodeURIComponent(code);
        f.innerHTML = `<div class="card stack"><h3>${ic('check')} شكرًا لك، اكتملت إجاباتك</h3><p style="margin:0">انسخ كود الاستجابة التالي وأرسله للباحث. إجاباتك لا تُحفظ إلا في هذا الكود.</p><div class="codebox">${esc(code)}</div><div class="row"><button class="btn" id="cpR">${ic('copy')} نسخ كود الاستجابة</button><a class="btn ghost" href="${esc(wa)}" target="_blank" rel="noopener noreferrer">${ic('whats')} إرسال عبر واتساب</a></div></div>`;
        $('#cpR', f).onclick = () => copyText(code); f.scrollIntoView({ block: 'start' });
      });
      f.scrollIntoView({ block: 'start' });
    };
  };

  /* import a scale (and optional responses) from the quick builder; rows use the builder raw convention */
  function importScale(o) {
    const s = blank(); s.name = o.name || 'مقياس مستورد'; s.construct = o.construct || ''; s.group = o.group || '';
    s.method = METHODS[o.method] ? o.method : 'likert'; s.scale = { preset: 'custom', labels: o.labels.slice(0, 7) }; if (o.instr) s.instr = o.instr;
    s.demo = { gender: false, age: false };
    const names = [...new Set(o.items.map(x => x.dim || 'عام'))]; s.dims = names.map(n => ({ id: uid('d'), name: n, def: '', w: Math.round(100 / names.length) }));
    s.items = o.items.map(x => ({ id: uid('i'), text: x.text, dim: s.dims[names.indexOf(x.dim || 'عام')].id, dir: x.dir || 1, key: x.key || 0, sv: x.sv ?? null, src: 'البناء السريع', ess: null, use: null, keep: true }));
    s.target = s.items.length; s.mult = 1;
    const v = version(s); let added = 0;
    (o.rows || []).forEach(r => { if (addResponse(s, v, r.vals, Object.assign({ src: o.src || 'مستورد' }, r.extra || {})) === 'ok') added++; });
    s.step = o.rows && o.rows.length ? 9 : 4; SCALES.save(s); checkBadges(); openId = s.id; go('build');
    return added;
  }
  const bestSteps = () => SCALES.list().reduce((m, s) => Math.max(m, status(s).filter(Boolean).length), 0);
  return { analyze, kept, formDef, importScale, bestSteps, open: id => { openId = id; go('build'); }, clear: () => { openId = null; } };
})();

/* ===== Quick scale builder — ported from the statistics lab (#lab/build): methods, templates, simulation,
   direct link, Google Forms script, Google Sheets collection. Collected responses open in the 12-step builder. ===== */
const QUICK = (() => {
  let view = null;
  const labHead = t => `<div class="section-title"><h2>${t}</h2><button class="btn sm ghost" id="qFull">${ic('build')} المعالج الكامل (12 خطوة)</button></div><div class="notice small"><b>مثل «بناء المقياس» في مختبر الإحصاء</b>اختر طريقة البناء، وولّد العبارات، وأرسل المقياس بلينك مباشر أو نموذج Google Forms، واجمع الردود في Google Sheet؛ ثم تفتح الردود في «ابنِ مقياسك» لتحليل الفقرات والثبات والصدق العاملي والمعايير.</div>`;
  const f2 = (x, d = 2) => fmt(x, d);
  const alphaOf = cols => { const rows = cols[0].map((_, i) => cols.map(c => c[i])); return PSY.itemStats(rows).alpha; };
  const SCI = { thurstone: ['لويس ثرستون (مع شيف)', 'الاتجاه نحو الكنيسة (1929)'], likert: ['رنسيس ليكرت', 'مقياس روزنبرج لتقدير الذات (1965)'], guttman: ['لويس جتمان', 'مقياس بوجاردس للبعد الاجتماعي (1925)'], osgood: ['تشارلز أوسجود', 'تمايز مفهوم «التعلم» (مجيد، 2014)'], rasch: ['جورج راش · ديفيد أندريش', 'استبيان «حب العلوم» (رايت وماسترز، 1982)'], cat: ['فريدريك لورد · ديفيد فايس', 'بنوك فقرات PROMIS'], sjt: ['موتوويدلو ودونيت وكارتر (1990)', 'اختبار الحكم الموقفي في UCAT للطب'], forced: ['آلن إدواردز · آنا براون', 'دراسة القيم؛ وجدول إدواردز (1959)'], aig: ['هومل وزملاؤه · جوتس وزملاؤه', 'مولدات الفقرات بالمحولات (2022–2024)'] };
  const relLbl = a => a >= 0.9 ? ['ممتاز', 'ok'] : a >= 0.8 ? ['جيد', 'ok'] : a >= 0.7 ? ['مقبول', 'mid'] : ['ضعيف', 'no'];
const AI={p:null};
function aiGet(){if(!AI.p)AI.p=(async()=>{try{if(window.claude&&typeof window.claude.use==="function")return await window.claude.use("sample")}catch(e){}return null})();return AI.p}
const AI_ERR={not_granted:"لم يتم السماح للموقع باستخدام Claude من حسابك.",sampling_disabled:"Claude غير متاح لهذا الحساب.",rate_limited:"وصلت للحد المسموح الآن، حاول بعد قليل.",session_expired:"سجّل الدخول في Claude مرة أخرى.",refused:"رفض Claude هذا الطلب؛ عدّل النص وحاول.",invalid_json:"الرد لم يكن بالشكل المتوقع؛ جرّب مرة أخرى."};
function aiExternal(prompt,el){
 el.innerHTML=`<p class="muted" style="margin:6px 0">استخدم حسابك في Claude أو ChatGPT: اضغط «انسخ الطلب» ثم افتح الموقع والصق (Ctrl+V).</p><div class="row"><button class="btn sm gold" data-cp>انسخ الطلب</button><a class="btn sm ghost" style="text-decoration:none" href="https://claude.ai/new" target="_blank" rel="noopener">افتح Claude</a><a class="btn sm ghost" style="text-decoration:none" href="https://chatgpt.com/" target="_blank" rel="noopener">افتح ChatGPT</a></div>`;
 el.querySelector("[data-cp]").onclick=()=>{const done=()=>toast("تم نسخ الطلب — الصقه في Claude");(navigator.clipboard?navigator.clipboard.writeText(prompt).then(done,()=>fallbackCopy(prompt,done)):fallbackCopy(prompt,done))};
}
function fallbackCopy(t,cb){const ta=document.createElement("textarea");ta.value=t;document.body.appendChild(ta);ta.select();try{document.execCommand("copy");cb()}catch(e){toast("تعذر النسخ")}ta.remove()}
const METHODS=[
 {id:"thurstone",t:"طريقة ثرستون — الفترات المتساوية ظاهريًا",en:"Thurstone Equal-Appearing Intervals",yr:1928,era:"كلاسيكية",d:"يكتب الباحث عبارات كثيرة تتدرج من الرفض الشديد للقبول الشديد، ويحكم عليها محكمون على متصل من 1 إلى 11، وتُختار العبارات ذات الاتفاق العالي بقيم موزعة. المفحوص يختار «أوافق» على ما ينطبق، ودرجته = وسيط قيم العبارات التي وافق عليها.",resp:["أوافق","لا أوافق"],bin:true,n:k=>[100,"50–100 محكم لتقدير قيم العبارات + 100 مفحوص على الأقل"]},
 {id:"likert",t:"طريقة ليكرت — التقديرات التجميعية",en:"Likert Summated Rating",yr:1932,era:"كلاسيكية (الأكثر استخدامًا)",d:"عبارات موجبة وسالبة يستجيب لها المفحوص على تدريج خماسي، وتجمع الدرجات (مع عكس العبارات السالبة) في درجة كلية. سهلة البناء والتصحيح، وتحلل بألفا كرونباخ والتحليل العاملي.",resp:["موافق بشدة","موافق","محايد","غير موافق","غير موافق بشدة"],n:k=>[Math.max(200,10*k),"5–10 أفراد لكل بند وبحد أدنى 200 للتحليل العاملي، و30–50 للتجربة الاستطلاعية"]},
 {id:"guttman",t:"طريقة جتمان — التدرج التراكمي",en:"Guttman Cumulative Scaling",yr:1944,era:"كلاسيكية",d:"عبارات مرتبة تصاعديًا في الشدة؛ من يوافق على عبارة قوية يوافق بالضرورة على ما قبلها. يُتحقق منها بمعامل إعادة الإنتاج (≥ 0.90).",resp:["نعم","لا"],bin:true,n:k=>[100,"100 فرد على الأقل لحساب معامل إعادة الإنتاج"]},
 {id:"osgood",t:"التمايز الدلالي — أوسجود",en:"Osgood Semantic Differential",yr:1957,era:"كلاسيكية",d:"يُقيَّم مفهوم واحد على أزواج من الصفات المتضادة (جيد–سيئ، قوي–ضعيف…) بتدريج سباعي، وتمثل ثلاثة أبعاد: التقييم، القوة، النشاط.",resp:["1","2","3","4","5","6","7"],n:k=>[100,"100 فرد على الأقل، ويزيد للتحليل العاملي"]},
 {id:"rasch",t:"نموذج راش / نظرية الاستجابة للمفردة",en:"Rasch Model / Item Response Theory",yr:1960,era:"حديثة",d:"تُبنى البنود بحيث تغطي مستويات مختلفة من شدة السمة (سهلة وصعبة التأييد)، ويقدَّر موقع الفرد والبند على متصل واحد بالـ logit مستقلًا عن العينة. يتيح الأداء التفاضلي للمفردة DIF وفحص المطابقة.",resp:["لا تنطبق أبدًا","تنطبق نادرًا","تنطبق أحيانًا","تنطبق دائمًا"],n:k=>[250,"150–250 فردًا لتقديرات مستقرة (Linacre)، و500 لنماذج IRT متعددة المعالم"]},
 {id:"sjt",t:"اختبارات الحكم الموقفي",en:"Situational Judgment Tests (SJT)",yr:1990,era:"حديثة",d:"يُعرض موقف واقعي وعدة تصرفات، ويختار المفحوص الأفضل (والأسوأ). تقيس المهارات والكفايات السلوكية وتقل فيها المرغوبية الاجتماعية، وتحتاج مفتاحًا من الخبراء.",resp:[],n:k=>[200,"200 فرد على الأقل + 5–10 خبراء لبناء مفتاح التصحيح"]},
 {id:"forced",t:"الاختيار القسري",en:"Forced-Choice / Thurstonian IRT",yr:2011,era:"حديثة",d:"عبارتان (أو أكثر) متساويتان في المرغوبية من أبعاد مختلفة، ويختار المفحوص الأكثر انطباقًا. تقلل التزييف والمرغوبية الاجتماعية، وتُحلل حديثًا بنموذج Thurstonian IRT (Brown & Maydeu-Olivares, 2011).",resp:[],n:k=>[300,"300 فرد على الأقل لنموذج Thurstonian IRT"]},
 {id:"cat",t:"الاختبار التكيفي المحوسب",en:"Computerized Adaptive Testing (CAT)",yr:1980,era:"حديثة",d:"بنك بنود معايَر بنظرية الاستجابة للمفردة؛ يختار الحاسوب البند التالي حسب إجابات الفرد، فيقيس بدقة بعدد بنود أقل. يتوقف عند وصول الخطأ المعياري لحد معين (مثلًا 0.30).",resp:["لا تنطبق أبدًا","تنطبق نادرًا","تنطبق أحيانًا","تنطبق دائمًا"],n:k=>[500,"500 فرد على الأقل لمعايرة بنك البنود"]},
 {id:"aig",t:"توليد البنود بالذكاء الاصطناعي",en:"AI-assisted Automatic Item Generation",yr:2020,era:"الأحدث",d:"تُستخدم النماذج اللغوية لتوليد بنود كثيرة من تعريف السمة، ثم يراجعها الخبراء ويحكمونها، وتُفحص إحصائيًا (التحليل العاملي، راش). يوفر الوقت لكنه لا يغني عن التحكيم والتجريب.",resp:["موافق بشدة","موافق","محايد","غير موافق","غير موافق بشدة"],n:k=>[Math.max(200,10*k),"مثل ليكرت: 10 أفراد لكل بند وبحد أدنى 200، بعد التحكيم"]}
];
const B={m:null,res:null,aud:null};
const respOf=m=>B.aud&&typeof AUD!=="undefined"&&m.resp&&m.resp.length?AUD.labels(m.resp,B.aud.p,B.aud.r):m.resp;
const q=c=>`«${c}»`;
function tplItems(m,c,dims,k){
 const pos=[`أشعر بـ${q(c)} في كثير من المواقف اليومية.`,`يلاحظ المحيطون بي ${q(c)} في تصرفاتي.`,`أفكر في الأمور المرتبطة بـ${q(c)} أكثر من غيري.`,`يظهر ${q(c)} لدي بوضوح عند التعامل مع الآخرين.`,`يصعب عليّ التحكم في ${q(c)} عندما أكون تحت ضغط.`,`يؤثر ${q(c)} في أدائي الدراسي أو المهني.`,`أشعر بـ${q(c)} عند مواجهة مواقف جديدة.`,`أتصرف بطريقة تعكس ${q(c)} حتى في المواقف البسيطة.`,`تتكرر لدي مواقف يظهر فيها ${q(c)}.`,`أجد نفسي مدفوعًا بـ${q(c)} في قراراتي.`,`يستمر ${q(c)} معي لفترات طويلة.`,`أتذكر مواقف كثيرة شعرت فيها بـ${q(c)}.`];
 const neg=[`نادرًا ما أشعر بـ${q(c)}.`,`لا يؤثر ${q(c)} في علاقاتي بالآخرين.`,`أستطيع أن أمضي يومي دون أن أفكر في ${q(c)}.`,`لا يلاحظ أحد ${q(c)} في سلوكي.`];
 const D=dims.length?dims:["عام"];const out=[];
 if(m==="osgood"){const P=[["جيد","سيئ","التقييم"],["مفيد","ضار","التقييم"],["ممتع","ممل","التقييم"],["مهم","تافه","التقييم"],["مريح","مزعج","التقييم"],["إيجابي","سلبي","التقييم"],["قوي","ضعيف","القوة"],["كبير","صغير","القوة"],["صعب","سهل","القوة"],["عميق","سطحي","القوة"],["نشط","خامل","النشاط"],["سريع","بطيء","النشاط"],["متغير","ثابت","النشاط"],["حاد","هادئ","النشاط"]];for(let i=0;i<k;i++){const p=P[i%P.length];out.push({text:`${p[0]} ـــ ـــ ـــ ـــ ـــ ـــ ـــ ${p[1]}`,dim:p[2],reverse:i%3===2})}return out}
 if(m==="guttman"){const L=["أحيانًا","مرة في الأسبوع على الأقل","عدة مرات في الأسبوع","يوميًا تقريبًا","لدرجة أنه يؤثر على بعض مهامي","لدرجة أنه يمنعني من بعض الأنشطة","لدرجة أنه يمنعني من معظم الأنشطة","لدرجة أنني أحتاج مساعدة متخصصة"];for(let i=0;i<k;i++)out.push({text:`أشعر بـ${q(c)} ${L[Math.min(i,L.length-1)]}${i>=L.length?` (مستوى ${i+1})`:""}.`,dim:"تراكمي",value:i+1});return out}
 if(m==="thurstone"){const V=[[11,`${q(c)} أمر ضروري ولا يمكن الاستغناء عنه.`],[10,`${q(c)} من أهم ما يميز الإنسان الناجح.`],[9,`أرى أن ${q(c)} له فوائد كثيرة.`],[8,`${q(c)} مفيد في أغلب المواقف.`],[7,`${q(c)} مقبول بشكل عام.`],[6,`لا أمانع وجود ${q(c)} في حياتي.`],[5,`ليس لدي رأي واضح حول ${q(c)}.`],[4,`أتحفظ أحيانًا على ${q(c)}.`],[3,`${q(c)} يسبب مشكلات أكثر مما يفيد.`],[2,`أرفض ${q(c)} في معظم المواقف.`],[1,`${q(c)} أمر ضار تمامًا ولا فائدة منه.`]];for(let i=0;i<k;i++){const j=Math.round(i*(V.length-1)/Math.max(1,k-1));out.push({text:V[j][1],dim:"الاتجاه",value:V[j][0]})}return out}
 if(m==="sjt"){const C=["في قاعة المحاضرات","أثناء العمل في مجموعة","مع أحد أفراد الأسرة","في مقابلة مهمة","مع صديق مقرب","أثناء امتحان","في موقف عام أمام الناس","عند تلقي نقد من شخص آخر"];for(let i=0;i<k;i++)out.push({text:`${C[i%C.length]}، واجهت موقفًا يتطلب ${q(c)}. ماذا تفعل؟`,dim:D[i%D.length],options:["أتعامل مع الموقف مباشرة وبهدوء","أطلب المساعدة أو المشورة من شخص أثق به","أؤجل التعامل مع الموقف لوقت لاحق","أنسحب من الموقف تمامًا"],key:"الأفضل: أ — الأسوأ: د (يُراجع المفتاح مع الخبراء)"});return out}
 if(m==="forced"){for(let i=0;i<k;i++){const a=pos[i%pos.length],b=pos[(i+5)%pos.length];out.push({text:"أي العبارتين تنطبق عليك أكثر؟",dim:`${D[i%D.length]} / ${D[(i+1)%D.length]}`,options:[a,b]})}return out}
 const lv=["سهل التأييد","متوسط","صعب التأييد"];
 for(let i=0;i<k;i++){const rev=(m==="likert"||m==="aig")&&(i%4===3);const t=rev?neg[Math.floor(i/4)%neg.length]:pos[(i-Math.floor(i/4))%pos.length];const it={text:t,dim:D[i%D.length],reverse:rev};if(m==="rasch"||m==="cat"){it.level=lv[i%3];it.value=+(-1.5+3*(i%3)/2+((i*37)%10-5)/20).toFixed(2)}out.push(it)}
 return out;
}
function simulate(items,n,binary,steps){
 const k=items.length,lam=items.map((_,i)=>.55+((i*13)%7)/20),thr=items.map(it=>it.value!=null&&Math.abs(it.value)<4?it.value/2:0);
 const gauss=()=>{let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v)};
 const rows=[];for(let p=0;p<n;p++){const th=gauss();rows.push(items.map((it,j)=>{const z=lam[j]*th+Math.sqrt(1-lam[j]**2)*gauss()-thr[j];if(binary)return z>0?1:0;const cut=[-1.3,-.45,.35,1.15,1.9].slice(0,steps-1);return 1+cut.filter(c=>z>c).length}))}
 return rows;
}
function labBuild(){
 const step=()=>{if(!B.m){view.innerHTML=labHead("البناء السريع للمقياس")+`<div class="panel"><h3>1) اختر طريقة البناء</h3><p class="muted" style="margin-top:0">مرتبة من الأقدم إلى الأحدث.</p><div class="meth">${METHODS.map(m=>`<button class="opt" data-m="${m.id}"><b>${m.t}</b><span class="yr">${m.yr}</span> <span class="chip ${m.era.startsWith("كلاس")?"":"gold"}">${m.era}</span><span class="en" style="display:block;text-align:start">${m.en}</span><small>${m.d}</small><span class="sci">${ic("user")} ${SCI[m.id][0]} · <b>نموذج حقيقي:</b> ${SCI[m.id][1]}</span></button>`).join("")}</div><div class="row"><button class="btn sm soft" id="qReal">${ic("book")} شاهد النماذج الحقيقية لكل طريقة ولكل عالم</button></div></div><div class="panel" id="collect0"></div>`;renderCollect($("#collect0"));$("#qReal").onclick=()=>LEARN.open("real");$$("[data-m]").forEach(b=>b.onclick=()=>{B.m=METHODS.find(x=>x.id===b.dataset.m);B.res=null;step();window.scrollTo(0,0)});return}
  const m=B.m;
  view.innerHTML=labHead("البناء السريع للمقياس")+`<div class="panel"><div class="row" style="justify-content:space-between"><span class="chip">${m.t} · ${m.yr}</span><button class="btn sm ghost" id="chg">تغيير الطريقة</button></div><p class="muted">${m.d}</p>
  <h3>2) بيانات المقياس</h3><div class="g2"><div><label class="fl" for="bLib">المفهوم المراد قياسه</label><select id="bLib">${typeof CLIB!=="undefined"?CLIB.list.map(c=>`<option value="${c.id}" ${(B.lib||"social_anx")===c.id?"selected":""}>${esc(c.n)}${c.dev?" (مقياس أصيل قيد التطوير)":""}${c.type==="attitude"?" (اتجاه)":""}</option>`).join(""):""}<option value="" ${B.lib===""?"selected":""}>مفهوم آخر (أكتبه بنفسي)</option></select><input type="text" id="bC" value="${esc(B.c||"")}" placeholder="اكتب اسم المفهوم" ${B.lib===""?"":"hidden"} style="margin-top:6px"></div><div><label class="fl" for="bP">الفئة المستهدفة</label><input type="text" id="bP" value="${esc(B.p||"طلاب الجامعة")}"></div></div>
  <div class="g2"><div><label class="fl" for="bK">عدد البنود (الأسئلة)</label><input type="number" id="bK" min="3" max="60" value="${B.k||12}"></div><div><label class="fl" for="bN">حجم العينة المخطط</label><input type="number" id="bN" min="10" max="2000" value="${B.n||200}"></div></div>
  <label class="fl" for="bPath">مسار القياس</label><select id="bPath"><option value="direct" ${B.path==="direct"?"selected":""}>مباشر: عبارات سلوكية صريحة لا تذكر اسم المفهوم</option><option value="indirect" ${B.path==="indirect"?"selected":""}>غير مباشر: مواقف غامضة، ومواقف سلوكية، وتقدير عن الآخرين، واختيار قسري</option><option value="both" ${!B.path||B.path==="both"?"selected":""}>الاثنان معًا (موصى به): 3 درجات ومؤشر الفجوة</option></select><div id="bLibInfo"></div><label class="fl" for="bD">الأبعاد (افصل بينها بفاصلة — اختياري)</label><input type="text" id="bD" value="${esc(B.d!=null?B.d:"المعرفي، الوجداني، السلوكي")}">${typeof AUD!=="undefined"?`<details class="aud-q" ${B.aud?"open":""}><summary><b>الفئة العمرية ولغة الصياغة</b> <span class="muted">(${B.aud?esc(AUD.prof(B.aud.p).t)+" · "+esc((AUD.REGISTERS.find(x=>x.id===B.aud.r)||{}).t||""):"اختياري: الافتراضي فصحى للشباب"})</span></summary>${AUD.selectorsHTML(B.aud||{p:"youth",r:"fusha"},"qa")}<div id="qaTips"></div><p class="muted" style="font-size:13.5px;margin:4px 0 0">تتغير ألفاظ العبارات والبدائل فقط؛ المعنى واتجاه العبارة وعدد البدائل ثابتة حتى لا يختل القياس. وللمقارنة بين نسخ لفئات مختلفة ابنِ نسخًا من المقياس نفسه في المعالج الكامل (خطوة جودة الصياغة) ليفحص تكافؤها تلقائيًا.</p></details>`:""}<div id="bRec"></div>
  <div class="row" style="margin-top:12px"><button class="btn" id="bTpl">ابنِ المقياس من مخطط المفهوم</button><button class="btn gold" id="bAI" hidden>ابنِ الأسئلة بالذكاء الاصطناعي</button></div><div id="bExt"></div></div><div id="bOut"></div>`;
  $("#chg").onclick=()=>{B.m=null;step()};
  const vals=()=>{B.lib=$("#bLib")?$("#bLib").value:"";B.path=$("#bPath")?$("#bPath").value:"direct";const LC=B.lib&&typeof CLIB!=="undefined"?CLIB.get(B.lib):null;if(LC){B.c=LC.n;$("#bD").value=LC.facets.map(f=>f[1]).join("، ")}else B.c=$("#bC").value.trim()||"السمة";B.p=$("#bP").value.trim();B.k=Math.max(3,Math.min(60,+$("#bK").value||12));B.n=Math.max(10,Math.min(2000,+$("#bN").value||100));B.d=$("#bD").value;return{c:B.c,p:B.p,k:B.k,n:B.n,dims:B.d.split(/[،,]/).map(s=>s.trim()).filter(Boolean)}};
  const aiPrompt=()=>{const v=vals();return `أنت خبير في بناء المقاييس النفسية باللغة العربية. ابنِ ${v.k} بندًا لقياس «${v.c}» لدى ${v.p||"الفئة المستهدفة"} باستخدام ${m.t} (${m.en}).${v.dims.length?` الأبعاد: ${v.dims.join("، ")} — وزّع البنود عليها بالتساوي.`:""}
قواعد إلزامية: ابدأ بتعريف إجرائي للمفهوم من نظرية منشورة واذكر جوانبه. لا تذكر اسم المفهوم أو مرادفًا صريحًا له داخل أي عبارة؛ اكتب مؤشرات محددة (فكرة أو شعور أو سلوك أو عرض جسدي) في موقف محدد. فكرة واحدة لكل عبارة، بلا نفي مزدوج ولا كلمات مطلقة ولا مصطلحات إكلينيكية، ومناسبة للفئة المستهدفة. العبارات العكسية تصف الطرف المقابل للسمة وليست نفيًا.${B.path!=="direct"?" أضف عبارات غير مباشرة: مواقف غامضة لها تفسير مهدد وتفسيران محايدان، ومواقف سلوكية ببدائل موزونة 0–2، وسؤال تقدير عن الآخرين؛ وضع أوزان البدائل في w.":""}${m.id==="likert"||m.id==="aig"?" اجعل نحو 25% من البنود سالبة (معكوسة).":""}${m.id==="thurstone"?" أعطِ كل عبارة قيمة مقترحة من 1 إلى 11 في value، موزعة على المتصل.":""}${m.id==="guttman"?" رتّب العبارات تصاعديًا في الشدة وضع ترتيبها في value.":""}${m.id==="rasch"||m.id==="cat"?" نوّع صعوبة التأييد وضع صعوبة مقترحة بالـ logit من -2 إلى 2 في value.":""}${m.id==="osgood"?" كل بند زوج صفات متضادة بالشكل: صفة ـــ صفة، وdim أحد: التقييم، القوة، النشاط.":""}${m.id==="sjt"?" كل بند موقف واقعي قصير مع 4 تصرفات في options ومفتاح الأفضل والأسوأ في key.":""}${m.id==="forced"?" كل بند: text = «أي العبارتين تنطبق عليك أكثر؟» وoptions عبارتان متقاربتان في المرغوبية الاجتماعية من بعدين مختلفين.":""}
أعد JSON فقط بهذا الشكل:
{"instructions":"تعليمات للمفحوص","items":[{"text":"...","dim":"...","reverse":false,"value":null,"options":[],"w":[],"key":""}],"notes":"ملاحظات للباحث عن التحكيم والتحليل"}`};
  const rec=()=>{const v=vals();const r=m.n(v.k);$("#bRec").innerHTML=`<div class="verdict ${v.n>=r[0]?"ok":"mid"}">العينة المقترحة لهذه الطريقة: ${r[1]} ← الحد الأدنى المقترح ≈ ${r[0]}. عينتك (${v.n}) ${v.n>=r[0]?"كافية ✔":"أقل من المقترح"}.</div>`;aiExternal(aiPrompt(),$("#bExt"))};
  ["bC","bP","bK","bN","bD"].forEach(i=>$("#"+i).oninput=rec);const libInfo=()=>{const sel=$("#bLib");if(!sel)return;const LC=sel.value&&typeof CLIB!=="undefined"?CLIB.get(sel.value):null;$("#bC").hidden=!!LC;$("#bD").disabled=!!LC;if(LC)$("#bD").value=LC.facets.map(f=>f[1]).join("، ");const F=LC?CLIB.fit(m.id,LC):null;$("#bLibInfo").innerHTML=LC?`<p class="muted" style="font-size:13.5px;margin:6px 0">${esc(LC.def)}</p>${F?`<div class="verdict ${F.block?"no":"mid"}">${esc(F.msg)}</div>`:""}`:`<p class="muted" style="font-size:13.5px;margin:6px 0">لمفهوم غير موجود في المكتبة استخدم زرار الذكاء الاصطناعي (أو انسخ الطلب لـ Claude)، وهو ملتزم بنفس قواعد البناء.</p>`;B.res=null;const o=$("#bOut");if(o)o.innerHTML="";rec()};if($("#bLib")){$("#bLib").onchange=libInfo;$("#bPath").onchange=()=>{B.res=null;rec()};libInfo()}if($("#qaP")){const au=()=>{B.aud={p:$("#qaP").value,r:$("#qaR").value};$("#qaTips").innerHTML=AUD.tipsHTML(B.aud.p);B.res=null;const o=$("#bOut");if(o)o.innerHTML=""};$("#qaP").onchange=au;$("#qaR").onchange=au;$("#qaTips").innerHTML=AUD.tipsHTML((B.aud||{p:"youth"}).p)}rec();
  $("#bTpl").onclick=()=>{const v=vals();const au=B.aud;if(!B.lib||typeof CLIB==="undefined"){toast("المفهوم ده مش في المكتبة: استخدم زرار الذكاء الاصطناعي أو انسخ الطلب لـ Claude");return}const L=CLIB.build(m.id,B.lib,B.path,v.k);if(!L.items.length){$("#bOut").innerHTML=`<div class="panel"><div class="verdict no">${esc(L.warn||"لا توجد عبارات لهذه الطريقة")}</div></div>`;return}let its=L.items;if(au&&typeof AUD!=="undefined"&&(au.r!=="fusha"))its=its.map(it=>Object.assign({},it,{text:AUD.adapt(it.text,au.p,au.r).text,options:(it.options||[]).map(o=>AUD.adapt(o,au.p,au.r).text)}));B.res={instructions:(au&&typeof AUD!=="undefined"&&AUD.instructions(au.p,au.r))||L.instructions,items:its,notes:L.notes+(au&&typeof AUD!=="undefined"?" مُكيّف لفئة: "+AUD.prof(au.p).t+"، راجع كل عبارة فالتحويل الآلي مسودة.":""),src:"lib",lib:B.lib};show()};
  aiGet().then(s=>{const btn=$("#bAI");if(!s||!btn)return;btn.hidden=false;btn.onclick=async()=>{btn.disabled=true;const old=btn.textContent;btn.textContent="جارٍ البناء… (قد يستغرق دقيقة)";
   try{const r=await s.json(aiPrompt(),{cache:false});if(!r||!Array.isArray(r.items)||!r.items.length)throw{code:"invalid_json"};B.res={instructions:String(r.instructions||""),items:r.items.map(x=>({text:String(x.text||""),dim:String(x.dim||""),reverse:!!x.reverse,value:x.value==null||x.value===""?null:+x.value,options:Array.isArray(x.options)?x.options.map(String):[],w:Array.isArray(x.w)&&x.w.length?x.w.map(Number):null,kind:Array.isArray(x.options)&&x.options.length?"sjt":"d",key:String(x.key||"")})).filter(x=>x.text),notes:String(r.notes||""),src:"ai"};show()}
   catch(e){toast(AI_ERR[e&&e.code]||"تعذر البناء الآن، استخدم القوالب الجاهزة أو انسخ الطلب لـ Claude.")}
   btn.disabled=false;btn.textContent=old}});
  if(B.res)show();
 };
 function show(){const m=B.m,r=B.res,o=$("#bOut");if(!o)return;const v={c:B.c,n:B.n};
  const resp=respOf(m);const its=r.items;
  let h=`<div class="panel"><div class="row" style="justify-content:space-between"><h3 style="margin:0">3) مقياس ${esc(v.c)}</h3><span class="chip ${r.src==="ai"?"gold":""}">${r.src==="ai"?"مولّد بالذكاء الاصطناعي":"من مخطط المفهوم: مسودة علمية للتحكيم"}</span></div>
  <p><b>التعليمات:</b> ${esc(r.instructions||defInstr(m))}</p>`;
  const LC=(r.lib||B.lib)&&typeof CLIB!=="undefined"?CLIB.get(r.lib||B.lib):null;its.forEach((it,i)=>{const QF=typeof CLIB!=="undefined"?CLIB.quality(it,LC||{n:v.c,leak:[],type:"trait"}):[];h+=`<div class="qitem"><b>${i+1}.</b> ${esc(it.text)} ${it.reverse?'<span class="rev">(معكوس)</span>':""}${it.kind&&it.kind!=="d"&&typeof CLIB!=="undefined"?` <span class="chip">${CLIB.KIND[it.kind]||""}</span>`:""}${QF.map(q=>` <span class="chip ${q[0]==="bad"?"bad":"warn"}">${esc(q[1])}</span>`).join("")}<div class="muted" style="font-size:13.5px">${it.dim?"البعد: "+esc(it.dim):""}${it.value!=null&&isFinite(it.value)?` · ${m.id==="thurstone"?"القيمة المقترحة":m.id==="guttman"?"الترتيب":"الصعوبة (logit)"}: ${it.value}`:""}${it.level?" · "+it.level:""}</div>${it.options&&it.options.length?`<div class="sc">${it.options.map((op,j)=>`<span>${"أبجدهوز"[j]||j+1}) ${esc(op)}</span>`).join("")}</div>${it.key?`<div class="muted" style="font-size:13.5px">المفتاح: ${esc(it.key)}</div>`:""}`:resp.length?`<div class="sc">${resp.map(x=>`<span>${x}</span>`).join("")}</div>`:""}</div>`});
  h+=`<p class="muted" style="font-size:14px">${esc(r.notes||"")}</p><div class="box warn"><b>قبل التطبيق</b>هذه مسودة: اعرضها على 5–10 محكمين (صدق المحتوى CVR/CVI)، وجرّبها على عينة استطلاعية 30–50 فردًا، وعدّل الصياغة.</div>
  ${LC?(()=>{const cv=CLIB.coverage(its,LC);return `<details class="box" open><summary><b>مخطط المفهوم (Blueprint)</b></summary><p><b>التعريف الإجرائي:</b> ${esc(LC.def)}</p><p><b>الإطار النظري:</b> ${esc(LC.theory)}</p><p><b>الجوانب:</b> ${LC.facets.map(f=>esc(f[1])).join("، ")}</p><p><b>خطة الصدق مقدمًا:</b> يُتوقع ارتباط موجب مع: ${LC.nomo.pos.map(esc).join("، ")}${LC.nomo.neg.length?"؛ وسالب مع: "+LC.nomo.neg.map(esc).join("، "):""}${LC.nomo.none&&LC.nomo.none.length?"؛ وارتباط ضعيف أو منعدم مع: "+LC.nomo.none.map(esc).join("، "):""}. سجّل هذه التوقعات قبل جمع البيانات.</p>${cv.n?`<p><b>مصفوفة التغطية (العبارات المباشرة):</b></p><table class="otbl"><thead><tr>${cv.head.map(x=>`<th>${x}</th>`).join("")}</tr></thead><tbody>${cv.rows.map(rw=>`<tr>${rw.map(x=>`<td>${typeof x==="number"?x:esc(x)}</td>`).join("")}</tr>`).join("")}</tbody></table><p class="muted" style="font-size:13.5px">العبارات العكسية ${cv.rev} من ${cv.n} (${Math.round(cv.pr*100)}%)${cv.pr<.15?": أقل من المقترح (15% إلى 40%)":cv.pr>.5?": أكثر من المقترح":""}${cv.empty.length?". جوانب بلا عبارات: "+cv.empty.map(esc).join("، "):""}</p>`:""}</details><details class="box"><summary><b>مفتاح التصحيح (لا يظهر للمفحوص)</b></summary><table class="otbl"><thead><tr><th>#</th><th>النوع</th><th>الجانب</th><th>التصحيح</th></tr></thead><tbody>${CLIB.keyRows(its).map(rw=>`<tr>${rw.map(x=>`<td>${esc(String(x))}</td>`).join("")}</tr>`).join("")}</tbody></table><p class="muted" style="font-size:13.5px">${its.some(x=>x.w)&&its.some(x=>x.kind==="d")?"تُحسب 3 درجات: المباشرة (نسبة من الحد الأقصى)، وغير المباشرة (مجموع الأوزان)، والمجمّعة. والفرق بين المباشرة وغير المباشرة مؤشر للتجميل أو ضعف الاستبصار، ويُحلل كمتغير بحثي.":""}</p></details>`})():""}<div class="row"><button class="btn sm" id="bCopy">نسخ المقياس</button><button class="btn sm ghost" onclick="window.print()">طباعة / PDF</button></div></div>`;
  const sim=["sjt","forced"].includes(m.id)?null:true;
  if(sim)h+=`<div class="panel"><h3>4) جرّب التحليل على عينة افتراضية (${v.n} فردًا)</h3><p class="muted" style="margin-top:0">نولّد استجابات افتراضية لـ ${v.n} فردًا على بنودك (بعد عكس البنود السالبة) لتتدرب على حساب الثبات.</p><button class="btn sm gold" id="bSim">ولّد البيانات واحسب الثبات</button><div id="bSimOut"></div></div>`;
  o.innerHTML=h;o.scrollIntoView({behavior:"smooth",block:"start"});
  $("#bCopy").onclick=()=>{const t=`مقياس ${v.c}\n${r.instructions||defInstr(m)}\n\n`+its.map((it,i)=>`${i+1}. ${it.text}${it.reverse?" (معكوس)":""}${it.options&&it.options.length?"\n   "+it.options.map((x,j)=>("أبجدهوز"[j]||j+1)+") "+x).join("\n   "):""}`).join("\n")+(resp.length&&!(its[0].options&&its[0].options.length)?`\n\nالبدائل: ${resp.join(" / ")}`:"");const done=()=>toast("تم نسخ المقياس");navigator.clipboard?navigator.clipboard.writeText(t).then(done,()=>fallbackCopy(t,done)):fallbackCopy(t,done)};
  sharePanel(o);
  const sb=$("#bSim");if(sb)sb.onclick=()=>{const bin=!!m.bin,steps=Math.max(2,resp.length||5);if(its.some(x=>x.options&&x.options.length)&&its.some(x=>!(x.options&&x.options.length))){toast("المحاكاة للعبارات المباشرة فقط؛ اختر المسار المباشر للتجربة");return}const rows=simulate(its,Math.min(v.n,2000),bin,steps);const cols=its.map((_,j)=>rows.map(rw=>rw[j]));const a=alphaOf(cols);
   
   $("#bSimOut").innerHTML=`<div class="kpis" style="margin-top:10px"><div class="kpi hl"><small>${bin?"KR-20 / ألفا":"ألفا كرونباخ"}</small><b>${f2(a,3)}</b></div><div class="kpi"><small>عدد الأفراد</small><b>${rows.length}</b></div><div class="kpi"><small>عدد البنود</small><b>${its.length}</b></div></div><div class="verdict ${relLbl(a)[1]}">الثبات ${relLbl(a)[0]} (بيانات افتراضية للتدريب فقط — النتيجة الحقيقية تأتي من تطبيقك الفعلي)</div><div class="row"><button class="btn sm" id="bSimGo">حلّل البيانات الافتراضية بالتفصيل في «ابنِ مقياسك»</button></div>`;$("#bSimGo").onclick=()=>{const Q=buildQ();respToScale(Q,rows.map(rw=>({n:"",d:[],a:rw.map((x,j)=>simIdx(Q,Q.items[j].r?(bin?1-x:steps+1-x):x,bin,steps))})),0,"افتراضية")}};
 }
 step();
}
const PUBLIC_URL="https://mo-psych-stats-lab.github.io/psych-scales-lab/";
const QKEY="scaleslab-lastQ";
const b64e=u8=>{let bin="";for(let i=0;i<u8.length;i+=8192)bin+=String.fromCharCode.apply(null,u8.subarray(i,i+8192));return btoa(bin).replace(/\+/g,"-").replace(/\//g,"_").replace(/=+$/,"")};
const b64d=s=>{s=s.replace(/-/g,"+").replace(/_/g,"/");while(s.length%4)s+="=";const bin=atob(s);return Uint8Array.from(bin,c=>c.charCodeAt(0))};
async function packQ(obj){const u=new TextEncoder().encode(JSON.stringify(obj));if(window.CompressionStream){try{const buf=new Uint8Array(await new Response(new Blob([u]).stream().pipeThrough(new CompressionStream("deflate-raw"))).arrayBuffer());return"z"+b64e(buf)}catch(e){}}return"j"+b64e(u)}
async function unpackQ(s){const k=s[0],u=b64d(s.slice(1));if(k==="j")return JSON.parse(new TextDecoder().decode(u));if(k==="z"){if(!window.DecompressionStream)throw new Error("old browser");return JSON.parse(await new Response(new Blob([u]).stream().pipeThrough(new DecompressionStream("deflate-raw"))).text())}throw new Error("bad")}
const hashStr=s=>{let h=5381;for(let i=0;i<s.length;i++)h=((h<<5)+h+s.charCodeAt(i))>>>0;return h.toString(36)};
function buildQ(){const m=B.m,r=B.res,dirs={likert:"d",aig:"d",thurstone:"b",guttman:"b",osgood:"a",rasch:"a",cat:"a",sjt:"o",forced:"o"};
 const nm=$("#shName")?$("#shName").checked:true,dm=$$("#shDemo input:checked").map(i=>i.value);
 const Q={v:1,t:"مقياس "+(B.c||""),ins:r.instructions||defInstr(m),sc:respOf(m),s:dirs[m.id]||"a",m:m.id,nm,dm,items:r.items.map(it=>{const o={x:it.text};if(it.dim)o.g=it.dim;if(m.id==="thurstone"&&it.value!=null)o.v=it.value;if(it.options&&it.options.length)o.o=it.options;if(it.w)o.w=it.w;if(it.reverse)o.r=1;return o})};
 Q.id=hashStr(JSON.stringify([Q.t,Q.sc,Q.items]));const g=gsCfg();if(g&&$("#shGS")&&$("#shGS").checked)Q.ep=g.url;return Q}
/* --- Google Sheets as the responses server (Apps Script web app) --- */
const GSKEY="scaleslab-gs";
function gsCfg(){try{const g=JSON.parse(localStorage.getItem(GSKEY)||"null");return g&&g.url&&g.key?g:null}catch(e){return null}}
function gsScript(key){return String.raw`// خادم ردود «مختبر المقاييس النفسية» — Google Apps Script
// 1) اضغط Run على الدالة setup ووافق على الصلاحيات
// 2) Deploy ← New deployment ← Web app ← Execute as: Me ← Who has access: Anyone ← Deploy
// 3) انسخ Web app URL والصقه في الموقع
var SECRET = "__KEY__";

function setup() {
  var ss = getSS_();
  Logger.log("تم إنشاء شيت الردود: " + ss.getUrl());
  Logger.log("الخطوة التالية: Deploy ← New deployment ← Web app ← Who has access: Anyone");
}
function getSS_() {
  var p = PropertiesService.getScriptProperties(), id = p.getProperty("SSID");
  if (id) { try { return SpreadsheetApp.openById(id); } catch (e) {} }
  var ss = SpreadsheetApp.create("ردود مقاييس مختبر المقاييس النفسية");
  p.setProperty("SSID", ss.getId());
  return ss;
}
function json_(o) {
  return ContentService.createTextOutput(JSON.stringify(o)).setMimeType(ContentService.MimeType.JSON);
}
function sheetFor_(ss, qid, title, create) {
  var tag = " - " + qid;
  var sh = ss.getSheets().filter(function (s) { var n = s.getName(); return n.slice(-tag.length) === tag; })[0];
  if (!sh && create) sh = ss.insertSheet(String(title || "مقياس").replace(/[:\\\/?*\[\]]/g, "").slice(0, 60) + tag);
  return sh;
}
function doPost(e) {
  var lock = LockService.getScriptLock();
  lock.waitLock(20000);
  try {
    var d = JSON.parse(e.postData.contents);
    var qid = String(d.qid || "").replace(/[^a-z0-9]/g, "").slice(0, 20);
    if (!qid) return json_({ ok: false });
    var sh = sheetFor_(getSS_(), qid, d.t, true);
    if (sh.getLastRow() === 0) sh.appendRow(["التاريخ"].concat(d.h || []).concat(["raw"]));
    var row = [new Date()];
    if (d.nm) row.push(String(d.n || "").slice(0, 100));
    (d.d || []).forEach(function (x) { row.push(String(x).slice(0, 100)); });
    (d.s || []).forEach(function (x) { row.push(x); });
    row.push(JSON.stringify({ n: d.n || "", d: d.d || [], a: d.a || [] }));
    sh.appendRow(row);
    return json_({ ok: true });
  } finally {
    lock.releaseLock();
  }
}
function doGet(e) {
  var p = (e && e.parameter) || {};
  if (p.key !== SECRET) return json_({ ok: false, error: "key" });
  if (p.ping) return json_({ ok: true, sheet: getSS_().getUrl() });
  var sh = sheetFor_(getSS_(), String(p.qid || ""), "", false);
  if (!sh || sh.getLastRow() < 2) return json_({ ok: true, rows: [] });
  var v = sh.getRange(2, sh.getLastColumn(), sh.getLastRow() - 1, 1).getValues();
  var rows = [];
  v.forEach(function (r) { try { rows.push(JSON.parse(r[0])); } catch (x) {} });
  return json_({ ok: true, rows: rows });
}`.replace("__KEY__",key)}
function renderGS(el,onChange){const g=gsCfg();
 if(g){el.innerHTML=`<div class="verdict ok" style="margin:8px 0">✔ مربوط بـ Google Sheets — الردود هتتحفظ تلقائيًا في الشيت بتاعك.</div><div class="row"><button class="btn sm ghost" id="gsTest">اختبر الربط</button>${g.sheet?`<a class="btn sm ghost" style="text-decoration:none" target="_blank" rel="noopener" href="${esc(g.sheet)}">افتح شيت الردود</a>`:""}<button class="btn sm ghost" id="gsOff">فصل الربط</button></div>`;
  $("#gsOff",el).onclick=()=>confirmBox("فصل الربط؟ الردود المحفوظة في الشيت لن تُمسح.","فصل الربط",()=>{try{localStorage.removeItem(GSKEY)}catch(e){}renderGS(el,onChange);onChange&&onChange()});
  $("#gsTest",el).onclick=async()=>{const r=await gsPing(g.url,g.key);toast(r?"✔ الربط شغال":"تعذر الاتصال؛ تأكد من النشر بـ Anyone")};return}
 let key="";try{const u=new Uint8Array(16);crypto.getRandomValues(u);key=b64e(u)}catch(e){key=Math.random().toString(36).slice(2)+Date.now().toString(36)}
 const code=gsScript(key);
 el.innerHTML=`<div class="box note" style="margin-top:8px"><b>الإعداد مرة واحدة فقط (دقيقتين) — بحسابك في جوجل</b><ol style="margin:4px 0"><li>اضغط «نسخ كود السيرفر».</li><li>افتح <a href="https://script.google.com/home/projects/create" target="_blank" rel="noopener">script.google.com</a> (مشروع جديد)، امسح الموجود والصق الكود، واحفظ (Ctrl+S).</li><li>اختر الدالة <b>setup</b> من الأعلى واضغط <b>Run ▶</b>، ووافق على الصلاحيات (لو ظهر «Google hasn't verified this app»: Advanced ← Go to project ← Allow).</li><li>اضغط <b>Deploy ← New deployment</b> ← اختر النوع <b>Web app</b> ← Execute as: <b>Me</b> ← Who has access: <b>Anyone</b> ← <b>Deploy</b>.</li><li>انسخ <b>Web app URL</b> (ينتهي بـ /exec) والصقه تحت واضغط «اربط».</li></ol>الكود فيه <b>مفتاح سري</b> خاص بك: بيه بس تقدر تسحب الردود، وهو محفوظ على جهازك ومش موجود في لينك المفحوصين.</div>
 <div class="codebox">${esc(code)}</div><div class="row" style="margin-top:8px"><button class="btn sm gold" id="gsCp">نسخ كود السيرفر</button><a class="btn sm ghost" style="text-decoration:none" target="_blank" rel="noopener" href="https://script.google.com/home/projects/create">افتح Google Apps Script</a></div>
 <label class="fl" for="gsUrl">Web app URL</label><input type="text" id="gsUrl" style="direction:ltr" placeholder="https://script.google.com/macros/s/.../exec"><div class="row" style="margin-top:8px"><button class="btn sm" id="gsSave">اربط</button></div>`;
 $("#gsCp",el).onclick=()=>{const d=()=>toast("تم نسخ كود السيرفر");navigator.clipboard?navigator.clipboard.writeText(code).then(d,()=>fallbackCopy(code,d)):fallbackCopy(code,d)};
 $("#gsSave",el).onclick=async()=>{const url=$("#gsUrl",el).value.trim();if(!/^https:\/\/script\.google\.com\/macros\/s\/[\w-]+\/exec$/.test(url)){toast("اللينك لازم يبدأ بـ https://script.google.com/macros/s/ وينتهي بـ /exec");return}
  const b=$("#gsSave",el);b.disabled=true;b.textContent="جارٍ الاختبار…";const r=await gsPing(url,key);b.disabled=false;b.textContent="اربط";
  if(!r){toast("الاتصال فشل: تأكد إنك نشرت بـ Who has access: Anyone وإنك لزقت نفس الكود");return}
  try{localStorage.setItem(GSKEY,JSON.stringify({url,key,sheet:r.sheet||""}))}catch(e){}toast("✔ تم الربط بـ Google Sheets");renderGS(el,onChange);onChange&&onChange()};
}
async function gsPing(url,key){try{const r=await fetch(url+"?ping=1&key="+encodeURIComponent(key));const j=await r.json();return j&&j.ok?j:null}catch(e){return null}}
async function gsPost(ep,payload){try{const r=await fetch(ep,{method:"POST",body:JSON.stringify(payload)});const j=await r.json();return j&&j.ok?"ok":"fail"}catch(e){try{await fetch(ep,{method:"POST",mode:"no-cors",body:JSON.stringify(payload)});return"maybe"}catch(e2){return"fail"}}}
function sharePanel(o){const p=document.createElement("div");p.className="panel";
 p.innerHTML=`<h3>5) أرسل المقياس للمفحوصين</h3>
 <label class="row" style="gap:8px"><input type="checkbox" id="shName" checked> <b>إلزام المفحوص بكتابة اسمه</b></label>
 <p class="muted" style="margin:6px 0 2px">بيانات إضافية تطلبها (اختياري):</p><div class="row" id="shDemo">${["النوع","العمر","الكلية / التخصص","الفرقة الدراسية"].map(d=>`<label class="row" style="gap:6px"><input type="checkbox" value="${d}"> ${d}</label>`).join("")}</div>
 <div class="panel" style="background:var(--good-soft);border-color:var(--good);margin:12px 0 0"><h3 style="margin-bottom:4px">📊 حفظ الردود تلقائيًا في Google Sheets (موصى به)</h3><p class="muted" style="margin:0">لما المفحوص يدوس «إرسال» إجاباته تتسجل مباشرة في شيت على حسابك، من غير أكواد.</p><label class="row" style="gap:8px;margin-top:8px" id="shGSw"><input type="checkbox" id="shGS" checked> احفظ ردود هذا المقياس في Google Sheet</label><div id="gsBox"></div></div>
 <div class="row" style="margin-top:12px"><button class="btn" id="shLink">🔗 أنشئ لينك مباشر للمقياس</button><button class="btn ghost" id="shGF">📝 أنشئ نموذج Google Forms</button></div><div id="shOut"></div>`;
 o.appendChild(p);
 const gsUpd=()=>{$("#shGSw").hidden=!gsCfg()};renderGS($("#gsBox"),gsUpd);gsUpd();
 const c=document.createElement("div");c.className="panel";c.id="collect";o.appendChild(c);renderCollect(c);
 $("#shLink").onclick=async()=>{const Q=buildQ();try{localStorage.setItem(QKEY,JSON.stringify(Q))}catch(e){}const code=await packQ(Q);const url=PUBLIC_URL+"#take/"+code;
  $("#shOut").innerHTML=`<div class="box tip" style="margin-top:12px"><b>اللينك جاهز ✔</b>ابعته لأي حد؛ بيفتح من غير أي تسجيل، ويجاوب على الأسئلة${Q.nm?" بعد كتابة اسمه (إلزامي)":""}. ${Q.ep?"إجاباته هتتسجل تلقائيًا في Google Sheet بتاعك، وتسحبها تحت بزرار واحد.":"بعد ما يخلّص هيظهر له «كود ردّ» يبعته لك على واتساب، وأنت تجمع الأكواد تحت. (اربط Google Sheets فوق عشان الردود توصل تلقائيًا.)"}</div><input type="text" readonly value="${esc(url)}" style="direction:ltr" onclick="this.select()"><div class="row" style="margin-top:8px"><button class="btn sm gold" id="cpL">نسخ اللينك</button><a class="btn sm ghost" style="text-decoration:none" target="_blank" rel="noopener" href="https://wa.me/?text=${encodeURIComponent("برجاء الإجابة على "+Q.t+":\n"+url)}">إرسال على واتساب</a><a class="btn sm ghost" style="text-decoration:none" target="_blank" rel="noopener" href="${esc(url)}">جرّبه بنفسك</a></div>${url.length>6000?`<p class="muted" style="font-size:14px">⚠ اللينك طويل (${url.length} حرف) لكثرة البنود؛ يعمل في المتصفحات وواتساب، وإذا قُطع فقلل عدد البنود أو استخدم Google Forms.</p>`:""}`;
  $("#cpL").onclick=()=>{const d=()=>toast("تم نسخ اللينك");navigator.clipboard?navigator.clipboard.writeText(url).then(d,()=>fallbackCopy(url,d)):fallbackCopy(url,d)};renderCollect(c)};
 $("#shGF").onclick=()=>{const Q=buildQ();const code=gfScript(Q);
  $("#shOut").innerHTML=`<div class="box note" style="margin-top:12px"><b>إنشاء النموذج على Google Forms تلقائيًا (دقيقة واحدة)</b><ol style="margin:4px 0"><li>اضغط «نسخ الكود».</li><li>افتح <a href="https://script.google.com/home/projects/create" target="_blank" rel="noopener">script.google.com</a> بحسابك في جوجل (مشروع جديد).</li><li>امسح الموجود والصق الكود (Ctrl+V) ثم اضغط <b>Run ▶</b>.</li><li>وافق على الصلاحيات (Review permissions ← اختر حسابك ← Allow). لو ظهر «Google hasn't verified this app» اضغط Advanced ثم Go to project.</li><li>في <b>Execution log</b> هيظهر لك <b>لينك النموذج للمشاركين</b> ولينك التعديل، ونموذجك هيكون في Google Forms ومعاه Google Sheet للردود.</li></ol>${Q.nm?"الاسم مضاف كسؤال إلزامي.":""} بعد جمع الردود: من Google Sheets اعمل File ← Download ← CSV، واستورده في «ابنِ مقياسك» (خطوة جمع البيانات ← Google Forms أو Excel).</div><div class="codebox">${esc(code)}</div><div class="row" style="margin-top:8px"><button class="btn sm gold" id="cpG">نسخ الكود</button><a class="btn sm ghost" style="text-decoration:none" target="_blank" rel="noopener" href="https://script.google.com/home/projects/create">افتح Google Apps Script</a></div>`;
  $("#cpG").onclick=()=>{const d=()=>toast("تم نسخ كود Google Forms");navigator.clipboard?navigator.clipboard.writeText(code).then(d,()=>fallbackCopy(code,d)):fallbackCopy(code,d)}};
}
function gfScript(Q){const opts=Q.sc&&Q.sc.length?Q.sc:null;const items=Q.items.map(it=>({text:it.x,options:it.o||opts||["نعم","لا"]}));
 return `// كود إنشاء نموذج Google Forms — مختبر المقاييس النفسية
function createScaleForm() {
  var TITLE = ${JSON.stringify(Q.t)};
  var INSTRUCTIONS = ${JSON.stringify(Q.ins)};
  var REQUIRE_NAME = ${Q.nm?"true":"false"};
  var DEMOGRAPHICS = ${JSON.stringify(Q.dm)};
  var ITEMS = ${JSON.stringify(items,null,1)};

  var form = FormApp.create(TITLE);
  form.setDescription(INSTRUCTIONS);
  form.setProgressBar(true);
  if (REQUIRE_NAME) form.addTextItem().setTitle('الاسم').setRequired(true);
  DEMOGRAPHICS.forEach(function (d) { form.addTextItem().setTitle(d).setRequired(true); });
  ITEMS.forEach(function (it, i) {
    form.addMultipleChoiceItem().setTitle((i + 1) + '. ' + it.text).setChoiceValues(it.options).setRequired(true);
  });
  var sheet = SpreadsheetApp.create(TITLE + ' - الردود');
  form.setDestination(FormApp.DestinationType.SPREADSHEET, sheet.getId());
  Logger.log('لينك النموذج للمشاركين: ' + form.getPublishedUrl());
  Logger.log('لينك تعديل النموذج: ' + form.getEditUrl());
  Logger.log('شيت الردود: ' + sheet.getUrl());
}`}
function scoreOf(Q,i,idx){if(Q.items[i]&&Q.items[i].w)return Q.items[i].w[idx];const it=Q.items[i],n=(it.o&&it.o.length)||(Q.sc&&Q.sc.length)||2;if(Q.s==="b")return idx===0?1:0;if(Q.s==="d")return n-idx;return idx+1}
function renderCollect(c){let Q=null;try{Q=JSON.parse(localStorage.getItem(QKEY)||"null")}catch(e){}
 if(!Q){c.innerHTML=`<h3>6) جمّع الردود وحللها</h3><p class="muted" style="margin:0">أنشئ اللينك المباشر أولًا، ثم الصق هنا أكواد الردود اللي توصلك.</p>`;return}
 const g=gsCfg();
 c.innerHTML=`<h3>6) جمّع ردود «${esc(Q.t)}» وحللها</h3>${Q.ep&&g&&g.url===Q.ep?`<div class="row" style="margin-bottom:10px"><button class="btn gold" id="rcPull">📥 اسحب الردود من Google Sheet</button>${g.sheet?`<a class="btn sm ghost" style="text-decoration:none" target="_blank" rel="noopener" href="${esc(g.sheet)}">افتح الشيت</a>`:""}</div><p class="muted" style="margin-top:0">أو لو حد بعتلك كود رد يدوي (لو النت كان ضعيف عنده) الصقه هنا:</p>`:`<p class="muted" style="margin-top:0">الصق رسائل الواتساب أو أكواد الردود كما هي (كل كود يبدأ بـ RSP1)، حتى لو كانت كلها في رسالة واحدة.</p>`}<textarea id="rcIn" style="min-height:110px" placeholder="RSP1:..."></textarea><div class="row" style="margin-top:8px"><button class="btn" id="rcGo">حوّل الأكواد إلى بيانات وحللها</button></div><div id="rcOut"></div>`;
 const pull=$("#rcPull",c);if(pull)pull.onclick=async()=>{pull.disabled=true;pull.textContent="جارٍ السحب…";let j=null;try{const r=await fetch(g.url+"?key="+encodeURIComponent(g.key)+"&qid="+encodeURIComponent(Q.id));j=await r.json()}catch(e){}pull.disabled=false;pull.textContent="📥 اسحب الردود من Google Sheet";
  if(!j||!j.ok){$("#rcOut",c).innerHTML=`<div class="verdict no">تعذر الاتصال بالشيت${j&&j.error==="key"?" (المفتاح غير مطابق — هل غيّرت الكود؟)":""}. جرّب من نسخة GitHub للموقع، وتأكد من النشر بـ Anyone.</div>`;return}
  if(!j.rows.length){$("#rcOut",c).innerHTML=`<div class="verdict mid">لا توجد ردود بعد على هذا المقياس.</div>`;return}respToScale(Q,j.rows,0,'Google Sheet')};
 $("#rcGo",c).onclick=()=>{const txt=$("#rcIn",c).value,re=/RSP1:([a-z0-9]+):([A-Za-z0-9_-]+)/g;let mm,rows=[],bad=0,other=0;while((mm=re.exec(txt))){if(mm[1]!==Q.id){other++;continue}try{rows.push(JSON.parse(new TextDecoder().decode(b64d(mm[2]))))}catch(e){bad++}}
  if(!rows.length){$("#rcOut",c).innerHTML=`<div class="verdict no">لم أجد أكواد صالحة لهذا المقياس${other?` (${other} كود تخص مقياسًا آخر)`:""}.</div>`;return}
  respToScale(Q,rows,bad,'كود رد')};
}
async function renderTake(code){let Q;try{Q=await unpackQ(code)}catch(e){view.innerHTML=`<div class="panel verdict no">الرابط غير صالح أو ناقص. اطلب من الباحث إرسال الرابط كاملًا.</div>`;return}
 const opts=it=>it.o||Q.sc||[];
 view.innerHTML=`<div class="lhead"><h2>${esc(Q.t)}</h2><p style="margin:6px 0 0">${esc(Q.ins)}</p></div>
 ${Q.nm||Q.dm.length?`<div class="panel">${Q.nm?`<label class="fl" for="tkN">الاسم <span style="color:var(--bad)">*</span></label><input type="text" id="tkN" autocomplete="name">`:""}${Q.dm.map((d,i)=>`<label class="fl" for="tkD${i}">${esc(d)} <span style="color:var(--bad)">*</span></label><input type="text" id="tkD${i}">`).join("")}</div>`:""}
 ${Q.items.map((it,i)=>`<div class="tk-item" id="tk${i}"><b>${i+1}.</b> ${esc(it.x)}<div class="tk-opts ${it.o?"col":""}">${opts(it).map((op,k)=>`<label><input type="radio" name="tq${i}" value="${k}"> ${esc(op)}</label>`).join("")}</div></div>`).join("")}
 <div class="row"><button class="btn" id="tkGo">إرسال الإجابات</button></div><div id="tkOut"></div><p class="muted" style="font-size:13.5px;margin-top:16px">لا يتم حفظ إجاباتك على أي خادم؛ تُرسل فقط للباحث عن طريق الكود.</p>`;
 $("#tkGo").onclick=async()=>{let miss=0;const a=Q.items.map((_,i)=>{const s=$(`input[name="tq${i}"]:checked`);$("#tk"+i).classList.toggle("miss",!s);if(!s)miss++;return s?+s.value:-1});
  const n=Q.nm?$("#tkN").value.trim():"";const d=Q.dm.map((_,i)=>$("#tkD"+i).value.trim());
  if(Q.nm&&!n){toast("اكتب اسمك أولًا");$("#tkN").focus();return}if(d.some(x=>!x)){toast("أكمل البيانات المطلوبة");return}if(miss){toast(`باقي ${miss} سؤال بدون إجابة`);const f=$(".tk-item.miss");f&&f.scrollIntoView({behavior:"smooth",block:"center"});return}
  const rc=`RSP1:${Q.id}:${b64e(new TextEncoder().encode(JSON.stringify({n,d,a})))}`;
  $("#tkGo").disabled=true;
  if(Q.ep&&/^https:\/\/script\.google\.com\/macros\/s\/[\w-]+\/exec$/.test(Q.ep)){$("#tkGo").textContent="جارٍ الإرسال…";
   const h=[...(Q.nm?["الاسم"]:[]),...Q.dm,...Q.items.map((it,i)=>"q"+(i+1)+": "+it.x.slice(0,40))];
   const st=await gsPost(Q.ep,{qid:Q.id,t:Q.t,h,nm:Q.nm,n,d,a,s:a.map((x,i)=>scoreOf(Q,i,x))});
   if(st!=="fail"){$("#tkGo").textContent="تم الإرسال";$("#tkOut").innerHTML=`<div class="panel" style="margin-top:12px;text-align:center"><div style="font-size:48px">✅</div><h3 style="justify-content:center">تم إرسال إجاباتك بنجاح${n?" يا "+esc(n):""}. شكرًا لمشاركتك!</h3><p class="muted">يمكنك إغلاق الصفحة الآن.</p><details class="steps" style="text-align:start"><summary>لو طلب منك الباحث كود الرد</summary><div class="codebox">${esc(rc)}</div></details></div>`;$("#tkOut").scrollIntoView({behavior:"smooth"});return}
   toast("تعذر الإرسال التلقائي — ابعت الكود للباحث")}
  $("#tkGo").textContent="إرسال الإجابات";$("#tkOut").innerHTML=`<div class="panel" style="margin-top:12px"><h3>✔ شكرًا${n?" يا "+esc(n):""}! خطوة أخيرة</h3><p>ابعت كود الرد ده للباحث (الأسهل بزر الواتساب):</p><div class="codebox">${esc(rc)}</div><div class="row" style="margin-top:10px"><a class="btn gold" style="text-decoration:none" target="_blank" rel="noopener" href="https://wa.me/?text=${encodeURIComponent("ردي على "+Q.t+":\n"+rc)}">إرسال على واتساب</a><button class="btn ghost" id="tkCp">نسخ الكود</button><a class="btn ghost" style="text-decoration:none" href="mailto:?subject=${encodeURIComponent("رد "+Q.t)}&body=${encodeURIComponent(rc)}">إرسال بالإيميل</a></div></div>`;
  $("#tkCp").onclick=()=>{const dn=()=>toast("تم نسخ الكود");navigator.clipboard?navigator.clipboard.writeText(rc).then(dn,()=>fallbackCopy(rc,dn)):fallbackCopy(rc,dn)};$("#tkOut").scrollIntoView({behavior:"smooth"})};
}
function defInstr(m){return m.id==="thurstone"?"ضع علامة أمام العبارات التي توافق عليها فقط.":m.id==="guttman"?"أجب بنعم أو لا على كل عبارة.":m.id==="osgood"?"قيّم المفهوم على كل زوج من الصفات بوضع علامة في الخانة الأقرب لرأيك (1 = الصفة اليمنى، 7 = الصفة اليسرى).":m.id==="sjt"?"اقرأ كل موقف واختر التصرف الذي ستقوم به فعلًا.":m.id==="forced"?"اختر من كل زوج العبارة التي تنطبق عليك أكثر، حتى لو انطبقت الاثنتان.":"اقرأ كل عبارة واختر البديل الذي يعبر عنك بصدق، لا توجد إجابات صحيحة أو خاطئة."}
  /* simulated score -> chosen option index, following the link's scoring direction */
  function simIdx(Q, x, bin, steps) { if (Q.s === 'b' || bin) return x ? 0 : 1; if (Q.s === 'd') return steps - x; return x - 1; }
  /* collected link/Sheet responses -> a scale in the 12-step builder (raw values follow the builder convention) */
  function respToScale(Q, rows, bad, src) {
    const opts = it => it.o || Q.sc || [];
    const kind = Q.m === 'sjt' ? 'sjt' : Q.m === 'forced' ? 'forced' : Q.s;
    let labels, raw;
    if (kind === 'sjt') { labels = ['أ', 'ب', 'ج', 'د']; raw = idx => 4 - idx; }
    else if (kind === 'forced') { labels = ['العبارة الأولى', 'العبارة الثانية']; raw = idx => idx === 0 ? 1 : 0; }
    else if (kind === 'b') { labels = Q.sc.slice(0, 2); raw = idx => idx === 0 ? 1 : 0; }
    else if (kind === 'd') { labels = Q.sc.slice(); raw = idx => labels.length - idx; }
    else { labels = Q.sc.slice().reverse(); raw = idx => idx + 1; }
    const items = Q.items.map(it => {
      let text = it.x;
      if (kind === 'sjt') text = [it.x, ...opts(it).slice(0, 4)].join(' | ');
      else if (kind === 'forced') text = opts(it).slice(0, 2).join(' | ');
      else if (Q.m === 'osgood') { const p = it.x.split(/\s*(?:ـــ\s*)+/).filter(Boolean); if (p.length === 2) text = p.join(' | '); }
      return { text, dim: it.g || '', dir: it.r ? -1 : 1, sv: it.v ?? null, key: 0 };
    });
    const seen = new Set(); let dup = 0; const out = [];
    rows.forEach(rw => { const k = JSON.stringify(rw); if (seen.has(k)) { dup++; return; } seen.add(k); out.push({ vals: Q.items.map((_, i) => { const a = rw.a && rw.a[i]; return a == null || a < 0 ? null : raw(a); }), extra: rw.n ? { name: String(rw.n).slice(0, 100) } : {} }); });
    const n = BUILD.importScale({ name: Q.t, method: Q.m, labels, instr: Q.ins, items, rows: out, src });
    toast(`تم تحويل ${n} رد إلى «ابنِ مقياسك»${dup ? ` (حُذف ${dup} مكرر)` : ''}${bad ? ` — ${bad} كود تالف` : ''}`);
  }
  VIEWS.qbuild = el => { view = el; el.classList.add('qb'); el.onclick = e => { if (e.target.closest('#qFull')) go('build'); }; labBuild(); };
  VIEWS.take = (el, code) => { view = el; el.classList.add('qb'); renderTake(code || ''); };
  return { reset: () => { B.m = null; B.res = null; } };
})();

/* ===== Objective questions & exam bank — بناء المقاييس النفسية =====
   Styled after the Minia University, Faculty of Education, Educational Psychology special diploma exams:
   70 marks / 2 hours; ✓/✗ statements; four-option MCQ (أ ب ج د); numerical problems with data; essays with "answer N only". */
const EXAM = (() => {
  // [statement, true?, explanation]
  const TF = [
    ['القياس هو التقدير الكمي للسمة وفق قواعد محددة، أما التقويم فيتضمن إصدار حكم واتخاذ قرار.', 1, 'القياس وصف كمي، والتقييم إعطاء قيمة، والتقويم حكم يتبعه قرار وعلاج؛ فالتقويم أشمل.'],
    ['المقياس الثابت صادق بالضرورة.', 0, 'الثبات شرط لازم غير كافٍ للصدق؛ قد يقيس المقياس بثبات شيئًا غير المقصود. أما الصادق فثابت بالضرورة.'],
    ['المقياس الصادق ثابت بالضرورة.', 1, 'معامل الصدق لا يتجاوز الجذر التربيعي لمعامل الثبات (√r_xx)، فلا صدق بلا ثبات.'],
    ['معامل ألفا كرونباخ يُستخدم فقط مع الفقرات الثنائية (صح/خطأ).', 0, 'ألفا تصلح للفقرات المتدرجة والثنائية؛ ومعادلة KR‑20 حالة خاصة منها للفقرات الثنائية.'],
    ['معادلة كودر–ريتشاردسون 21 تفترض تساوي صعوبة الفقرات.', 1, 'KR‑21 تستخدم المتوسط والتباين فقط وتفترض تساوي الصعوبة، لذا تعطي قيمة أقل من KR‑20 أو مساوية لها.'],
    ['تُستخدم معادلة سبيرمان–براون لتصحيح معامل الارتباط بين نصفي المقياس.', 1, 'الارتباط بين النصفين يمثل ثبات نصف المقياس، فيُصحَّح بـ r = 2r½ ÷ (1 + r½).'],
    ['زيادة عدد فقرات المقياس المتجانسة تؤدي غالبًا إلى زيادة ثباته.', 1, 'هذا ما تتنبأ به معادلة سبيرمان–براون العامة: r_new = n r ÷ (1 + (n − 1) r).'],
    ['يرتفع معامل الثبات كلما زاد تجانس العينة في السمة المقاسة.', 0, 'العكس: تجانس العينة يقلل التباين فينخفض معامل الثبات؛ والعينة غير المتجانسة ترفعه.'],
    ['ثبات إعادة الاختبار يقيس الاستقرار عبر الزمن.', 1, 'يُطبق المقياس مرتين على العينة نفسها بفاصل زمني (عادة أسبوعان) ويُحسب الارتباط بين التطبيقين.'],
    ['الصدق الظاهري هو أقوى أنواع الصدق.', 0, 'الصدق الظاهري أضعفها؛ يعتمد على مظهر الفقرات لا على أدلة تجريبية أو تحكيم منهجي.'],
    ['يعتمد صدق المحتوى على آراء المحكمين في مدى تمثيل الفقرات للسمة.', 1, 'ويُقدَّر كميًا بنسبة صدق المحتوى CVR لـ Lawshe أو مؤشر CVI أو معامل أيكن V.'],
    ['الصدق التنبؤي نوع من صدق المحك يكون فيه قياس المحك لاحقًا لتطبيق المقياس.', 1, 'في الصدق التلازمي يُقاس المحك في الوقت نفسه، وفي التنبؤي يُقاس بعد فترة.'],
    ['صدق المقارنة الطرفية يقوم على المقارنة بين أعلى 27% وأدنى 27% من الدرجات.', 1, 'وتُختبر دلالة الفرق بين المجموعتين باختبار (ت)؛ والفقرة المميزة تفرق بينهما بدلالة.'],
    ['التحليل العاملي الاستكشافي يُستخدم للتحقق من بنية محددة مسبقًا بنموذج نظري.', 0, 'الاستكشافي يكشف البنية دون فرض مسبق؛ أما التوكيدي CFA فيختبر بنية محددة سلفًا.'],
    ['قيمة KMO الأقل من 0.50 تعني أن البيانات غير مناسبة للتحليل العاملي.', 1, 'KMO ≥ .60 مقبول و≥ .80 جيد جدًا؛ وأقل من .50 غير مقبول (Kaiser).'],
    ['محك كايزر يُبقي العوامل التي جذرها الكامن أكبر من الواحد الصحيح.', 1, 'ويُفضَّل معه التحليل الموازي ومنحنى المنحدر لأنه يبالغ أحيانًا في عدد العوامل.'],
    ['يُقبل تشبع الفقرة على العامل إذا بلغ 0.30 فأكثر.', 1, 'هذا محك جيلفورد الشائع؛ ويرفعه بعضهم إلى .40 مع العينات الصغيرة.'],
    ['معامل الصعوبة هو نسبة من أجابوا عن الفقرة إجابة خاطئة.', 0, 'معامل الصعوبة p = نسبة من أجابوا إجابة صحيحة (ولذلك يسميه بعضهم معامل السهولة)؛ فكلما ارتفع كانت الفقرة أسهل.'],
    ['أفضل قيمة لمعامل الصعوبة لتحقيق أعلى تمييز هي 0.50.', 1, 'عند p = .50 يبلغ التباين pq أقصاه (.25)، فيتيح أعلى قدرة على التمييز.'],
    ['معامل التمييز السالب يدل على أن الفقرة جيدة.', 0, 'التمييز السالب يعني أن المجموعة الدنيا أجابت أفضل من العليا؛ فالفقرة معيبة وتُحذف أو تُراجع.'],
    ['يُعد معامل التمييز 0.40 فأكثر ممتازًا وفق تصنيف إيبل.', 1, 'إيبل: ≥ .40 ممتاز، .30–.39 جيد، .20–.29 مقبول يحتاج تحسينًا، < .20 ضعيف.'],
    ['في طريقة ليكرت تُعكس درجات العبارات السالبة عند التصحيح.', 1, 'في التدريج الخماسي تصبح الدرجة (6 − الدرجة الخام) للعبارة السالبة.'],
    ['في طريقة ثرستون يستجيب المفحوص على تدريج خماسي لكل عبارة.', 0, 'في ثرستون يختار المفحوص العبارات التي يوافق عليها فقط، ودرجته وسيط القيم المقياسية لهذه العبارات.'],
    ['معامل إعادة الإنتاج في طريقة جتمان يجب ألا يقل عن 0.90.', 1, 'CR = 1 − (عدد الأخطاء ÷ عدد الاستجابات الكلي)، والقيمة ≥ .90 تدل على سلم تراكمي.'],
    ['مقياس التمايز الدلالي من إعداد أوسجود ويقوم على أزواج من الصفات المتضادة.', 1, 'وتنتظم الصفات في ثلاثة أبعاد: التقييم والقوة والنشاط.'],
    ['تتميز نظرية الاستجابة للمفردة بأن تقدير خصائص المفردة يعتمد على العينة المستخدمة.', 0, 'العكس؛ من أهم مزاياها استقلال تقدير معالم المفردة عن العينة، وتقدير القدرة عن مجموعة المفردات.'],
    ['النموذج أحادي المعلم (راش) يقدّر صعوبة المفردة فقط.', 1, 'ويضيف النموذج الثنائي معلم التمييز، والثلاثي معلم التخمين.'],
    ['المستوى الاسمي من مستويات القياس يسمح بإجراء العمليات الحسابية كالجمع والطرح.', 0, 'الاسمي للتصنيف فقط؛ والفئوي (الفتري) يسمح بالجمع والطرح، والنسبي يسمح بكل العمليات لوجود صفر مطلق.'],
    ['معظم المقاييس النفسية تقع في المستوى الفئوي (الفتري) لعدم وجود صفر مطلق.', 1, 'الصفر في المقاييس النفسية اعتباري لا يعني انعدام السمة.'],
    ['الدرجة التائية T متوسطها 50 وانحرافها المعياري 10.', 1, 'T = 50 + 10z، وتتخلص من الكسور والإشارات السالبة في الدرجة المعيارية z.'],
    ['الرتبة المئينية 80 تعني أن المفحوص أجاب عن 80% من الفقرات إجابة صحيحة.', 0, 'تعني أن 80% من أفراد المجموعة المعيارية حصلوا على درجة أقل منه.'],
    ['المعايير المحلية تُشتق من عينة تقنين ممثلة للمجتمع الذي سيُطبق عليه المقياس.', 1, 'ولذلك لا يصح استخدام معايير أجنبية على بيئة مختلفة دون تقنين.'],
    ['الخطأ المعياري للقياس يزيد كلما ارتفع معامل الثبات.', 0, 'SEM = SD √(1 − r)، فيقل الخطأ كلما ارتفع الثبات.'],
    ['تُعد المرغوبية الاجتماعية من مصادر التحيز في استجابات مقاييس التقرير الذاتي.', 1, 'وتُضبط بمقاييس الكذب أو المرغوبية، والاختيار القسري، وضمان السرية.'],
    ['يُفضَّل أن تتضمن العبارة الواحدة أكثر من فكرة لاختصار عدد الفقرات.', 0, 'العبارة المركبة (ذات الفكرتين) تُربك المفحوص؛ القاعدة: فكرة واحدة لكل عبارة.'],
    ['يُستخدم جدول المواصفات لضمان تمثيل الفقرات لأبعاد السمة بأوزان مناسبة.', 1, 'وهو أساس صدق المحتوى؛ يحدد عدد الفقرات لكل بعد بحسب وزنه النسبي.'],
    ['يُكتب في مرحلة «فيض الفقرات» عدد من الفقرات مساوٍ تمامًا للعدد النهائي المطلوب.', 0, 'يُكتب عادة ضعف العدد أو أكثر (1.5–2 مرة) لأن بعض الفقرات ستُحذف بعد التحكيم والتحليل.'],
    ['التجريب الاستطلاعي يكشف غموض التعليمات والعبارات ويقدّر زمن التطبيق.', 1, 'ويُجرى على عينة صغيرة (30–50) من خارج عينة الدراسة الأساسية.'],
    ['معامل الاتفاق كابا لكوهين يصحح نسبة الاتفاق المشاهدة من أثر الصدفة.', 1, 'κ = (Po − Pe) ÷ (1 − Pe).'],
    ['الصدق التقاربي يعني ارتباط المقياس ارتباطًا منخفضًا بمقاييس سمات مختلفة.', 0, 'هذا هو الصدق التمييزي (الاختلافي)؛ أما التقاربي فارتباط مرتفع بمقاييس السمة نفسها.'],
    ['نسبة صدق المحتوى CVR تتراوح بين −1 و +1.', 1, 'CVR = (ne − N/2) ÷ (N/2)؛ تساوي +1 إذا اتفق كل المحكمين و−1 إذا لم يعتبرها أحد أساسية.'],
    ['معامل أوميغا ماكدونالد يفترض تساوي تشبعات الفقرات على العامل.', 0, 'ألفا هي التي تفترض تكافؤ التشبعات (tau‑equivalence)؛ أوميغا لا تفترض ذلك فتكون أدق غالبًا.']
  ];
  // [stem, options, correct index, explanation]
  const MCQ = [
    ['التقدير الكمي لسمة أو خاصية وفق قواعد محددة هو:', ['التقويم', 'القياس', 'التقييم', 'الاختبار'], 1, 'القياس = التكميم وفق قواعد.'],
    ['مستوى القياس الذي يتميز بوجود صفر مطلق هو:', ['الاسمي', 'الرتبي', 'الفئوي', 'النسبي'], 3, 'النسبي (كالطول والوزن) له صفر حقيقي.'],
    ['من أدوات تحديد مدى تمثيل الفقرات لأبعاد السمة:', ['جدول المواصفات', 'معامل الصعوبة', 'الدرجة التائية', 'منحنى المنحدر'], 0, 'جدول المواصفات يوزع الفقرات على الأبعاد بأوزانها.'],
    ['معادلة كودر–ريتشاردسون 20 تُستخدم مع الفقرات:', ['المتدرجة خماسيًا', 'الثنائية (صفر/واحد)', 'المقالية', 'ذات التدريج السباعي'], 1, 'KR‑20 للفقرات الثنائية، وألفا تعميمها.'],
    ['إذا كان الارتباط بين نصفي اختبار 0.60، فإن ثباته بمعادلة سبيرمان–براون يساوي:', ['0.60', '0.70', '0.75', '0.80'], 2, '2(.60) ÷ (1 + .60) = 1.20 ÷ 1.60 = .75.'],
    ['معادلة جتمان للتجزئة النصفية لا تشترط:', ['تساوي تبايني النصفين', 'وجود فقرات', 'تطبيق الاختبار', 'حساب التباين الكلي'], 0, 'سبيرمان–براون تفترض تساوي تبايني النصفين، وجتمان لا تفترضه.'],
    ['من العوامل التي تؤدي إلى انخفاض معامل الثبات:', ['زيادة عدد الفقرات', 'تجانس العينة في السمة', 'وضوح التعليمات', 'موضوعية التصحيح'], 1, 'التجانس يقلل التباين فيخفض معامل الثبات.'],
    ['ثبات إعادة الاختبار يُعرف أيضًا بمعامل:', ['التكافؤ', 'الاستقرار', 'الاتساق الداخلي', 'التمييز'], 1, 'إعادة الاختبار = الاستقرار؛ والصور المتكافئة = التكافؤ.'],
    ['النوع من الصدق الذي يعتمد على آراء المحكمين:', ['الصدق العاملي', 'صدق المحتوى', 'الصدق التنبؤي', 'الصدق التقاربي'], 1, 'يُقدر بـ CVR وCVI.'],
    ['إذا اعتبر 8 من 10 محكمين أن الفقرة «أساسية»، فإن CVR يساوي:', ['0.80', '0.60', '0.40', '0.30'], 1, '(8 − 5) ÷ 5 = 0.60.'],
    ['يُحسب صدق المقارنة الطرفية باختبار الفرق بين:', ['الذكور والإناث', 'أعلى 27% وأدنى 27% من الدرجات', 'التطبيق الأول والثاني', 'النصف الفردي والزوجي'], 1, 'ويُستخدم اختبار (ت) لعينتين مستقلتين.'],
    ['ارتباط المقياس ارتباطًا مرتفعًا بمقياس آخر يقيس السمة نفسها يمثل الصدق:', ['التمييزي', 'التقاربي', 'الظاهري', 'التنبؤي'], 1, 'التقاربي = ارتباط مرتفع مع المشابه.'],
    ['الغرض الرئيسي من التحليل العاملي الاستكشافي في بناء المقاييس:', ['حساب الثبات', 'الكشف عن البنية الكامنة للفقرات', 'تقدير الصعوبة', 'حساب المعايير'], 1, 'يكشف عدد العوامل وتشبع الفقرات بها.'],
    ['قيمة KMO تساوي 0.85 تعني أن كفاية العينة:', ['غير مقبولة', 'ضعيفة', 'جيدة جدًا', 'منعدمة'], 2, 'Kaiser: .80–.89 جديرة بالثناء (جيدة جدًا).'],
    ['اختبار بارتليت للكروية يختبر فرض أن مصفوفة الارتباطات:', ['مصفوفة وحدة', 'سالبة', 'كبيرة', 'متماثلة'], 0, 'الدلالة تعني أنها ليست مصفوفة وحدة فيصلح التحليل العاملي.'],
    ['التدوير المتعامد الأكثر استخدامًا في التحليل العاملي:', ['أوبليمين', 'بروماكس', 'فاريماكس', 'كوارتيمين المائل'], 2, 'فاريماكس متعامد؛ أوبليمين وبروماكس مائلان.'],
    ['معامل الصعوبة لفقرة أجاب عنها 30 من 60 إجابة صحيحة:', ['0.30', '0.50', '0.60', '2.00'], 1, '30 ÷ 60 = .50.'],
    ['إذا أجاب 18 من المجموعة العليا (20) و6 من الدنيا (20)، فمعامل التمييز:', ['0.30', '0.45', '0.60', '0.90'], 2, '(18 − 6) ÷ 20 = .60.'],
    ['تُحذف الفقرة إذا كان معامل تمييزها:', ['0.45', '0.35', '0.25', '0.10'], 3, 'أقل من .20 ضعيف وفق إيبل.'],
    ['في مقياس ليكرت خماسي، الدرجة المعكوسة لاستجابة «4» على عبارة سالبة:', ['1', '2', '4', '5'], 1, '6 − 4 = 2.'],
    ['الطريقة التي يقدر فيها المحكمون العبارات على متصل من 11 فئة:', ['ليكرت', 'ثرستون', 'جتمان', 'أوسجود'], 1, 'الفترات المتساوية ظاهريًا (1928).'],
    ['الطريقة التي تقوم على الترتيب التراكمي للعبارات:', ['ليكرت', 'ثرستون', 'جتمان', 'التمايز الدلالي'], 2, 'من يوافق على عبارة قوية يوافق على ما قبلها.'],
    ['أبعاد التمايز الدلالي عند أوسجود هي:', ['المعرفي والوجداني والسلوكي', 'التقييم والقوة والنشاط', 'الصدق والثبات والموضوعية', 'الشدة والتكرار والمدة'], 1, 'Evaluation, Potency, Activity.'],
    ['النموذج الذي يضيف معلم التخمين في نظرية الاستجابة للمفردة:', ['أحادي المعلم', 'ثنائي المعلم', 'ثلاثي المعلم', 'نموذج راش'], 2, 'معالم: الصعوبة b، والتمييز a، والتخمين c.'],
    ['من مزايا نظرية الاستجابة للمفردة مقارنة بالنظرية الكلاسيكية:', ['سهولة الحسابات', 'استقلال معالم المفردة عن العينة', 'لا تحتاج عينات كبيرة', 'لا تحتاج حاسبًا'], 1, 'وكذلك تقدير خطأ القياس لكل مستوى قدرة.'],
    ['الدرجة التائية المقابلة لـ z = 1.5:', ['51.5', '60', '65', '75'], 2, 'T = 50 + 10(1.5) = 65.'],
    ['الدرجة التساعية (Stanine) متوسطها وانحرافها:', ['50 و10', '5 و2', '100 و15', '0 و1'], 1, 'Stanine = 5 + 2z، من 1 إلى 9.'],
    ['إذا كان الانحراف المعياري 10 ومعامل الثبات 0.91، فالخطأ المعياري للقياس:', ['1', '3', '9', '0.9'], 1, '10 × √(1 − .91) = 10 × .30 = 3.'],
    ['أول قواعد بناء المقياس النفسي:', ['حساب الثبات', 'الوضوح المفاهيمي وتعريف السمة', 'التحليل العاملي', 'اشتقاق المعايير'], 1, 'التعريف النظري والإجرائي للسمة.'],
    ['العينة الاستطلاعية في التجريب الأولي للمقياس تتراوح عادة بين:', ['5–10', '30–50', '500–1000', '1000 فأكثر'], 1, 'لكشف الغموض وتقدير الزمن.'],
    ['الحد الأدنى المقبول لمعامل ألفا في البحوث النفسية غالبًا:', ['0.30', '0.50', '0.70', '0.99'], 2, 'Nunnally: .70 للبحوث، و.80–.90 لاتخاذ قرارات فردية.'],
    ['صعوبة فقرة ثنائية p = 0.70؛ تباين الفقرة pq يساوي:', ['0.21', '0.30', '0.49', '0.70'], 0, '.70 × .30 = .21.'],
    ['ميل المفحوص إلى الموافقة على العبارات بغض النظر عن محتواها:', ['المرغوبية الاجتماعية', 'أسلوب الموافقة (الإذعان)', 'التطرف', 'الهالة'], 1, 'يُضبط بتضمين عبارات سالبة.'],
    ['معامل الاتفاق بين مصححَين مع تصحيح أثر الصدفة:', ['ألفا', 'كابا كوهين', 'سبيرمان–براون', 'KR‑21'], 1, 'κ = (Po − Pe) ÷ (1 − Pe).'],
    ['معامل أيكن (Aiken V) يُستخدم لتقدير:', ['الثبات', 'صدق المحتوى من تقديرات المحكمين', 'التمييز', 'الصعوبة'], 1, 'V = Σs ÷ [n(c − 1)].'],
    ['تُستخدم الطريقة الطرفية 27% لأنها:', ['تعطي أكبر تمييز مع أكبر عدد ممكن من الأفراد', 'أسهل في الحساب فقط', 'تلغي الحاجة للثبات', 'تقدر الصدق التنبؤي'], 0, 'نسبة كيلي (Kelley, 1939).'],
    ['الأداء التفاضلي للمفردة (DIF) يعني:', ['اختلاف صعوبة الفقرة بين مجموعتين متساويتين في السمة', 'ارتفاع ثبات الفقرة', 'تشبع الفقرة على عاملين', 'اختلاف الزمن'], 0, 'ويُكشف بطريقة مانتل–هانزل أو الانحدار اللوجستي.'],
    ['أي مما يلي ليس من شروط صياغة العبارة الجيدة؟', ['فكرة واحدة', 'لغة مناسبة للفئة', 'النفي المزدوج', 'الوضوح والإيجاز'], 2, 'النفي المزدوج يربك المفحوص.'],
    ['معادلة KR‑21 بالنسبة لـ KR‑20 تعطي قيمة:', ['أكبر دائمًا', 'أقل أو مساوية', 'مساوية دائمًا', 'سالبة'], 1, 'لافتراضها تساوي صعوبة الفقرات.'],
    ['عدد الفقرات اللازمة لرفع ثبات اختبار من 0.60 إلى 0.75 يتطلب إطالته:', ['مرة واحدة', 'مرتين', '3 مرات', '4 مرات'], 1, 'n = .75(1 − .60) ÷ [.60(1 − .75)] = .30 ÷ .15 = 2.']
  ];
  // numerical problems: q (html), a (html)
  const PROB = [
    { t: 'ألفا كرونباخ', q: 'طُبق مقياس من 5 فقرات على عينة، فكانت تباينات الفقرات: 0.8، 1.1، 0.9، 1.2، 1.0، وتباين الدرجة الكلية 14. احسب معامل ألفا كرونباخ وفسّره.', a: 'مجموع تباينات الفقرات = 0.8 + 1.1 + 0.9 + 1.2 + 1.0 = 5.0<br>α = [k ÷ (k − 1)] × [1 − (Σσ²ᵢ ÷ σ²ₓ)] = (5 ÷ 4) × (1 − 5 ÷ 14) = 1.25 × 0.643 = <b>0.80</b><br>التفسير: ثبات جيد (≥ 0.70 مقبول للبحوث، ≥ 0.80 جيد).' },
    { t: 'KR‑20', q: 'اختبار من 10 فقرات ثنائية، مجموع (p × q) للفقرات = 1.8، وتباين الدرجات الكلية = 7.5. احسب معامل الثبات بمعادلة KR‑20.', a: 'KR‑20 = [k ÷ (k − 1)] × [1 − (Σpq ÷ σ²ₓ)] = (10 ÷ 9) × (1 − 1.8 ÷ 7.5) = 1.111 × 0.76 = <b>0.84</b> ← ثبات جيد.' },
    { t: 'KR‑21', q: 'اختبار من 20 فقرة ثنائية، متوسط درجاته 12 وتباينها 16. احسب الثبات بمعادلة KR‑21، ولماذا تكون قيمتها غالبًا أقل من KR‑20؟', a: 'KR‑21 = [k ÷ (k − 1)] × [1 − M(k − M) ÷ (k σ²)] = (20 ÷ 19) × [1 − 12 × 8 ÷ (20 × 16)] = 1.053 × (1 − 0.30) = <b>0.74</b><br>تكون أقل أو مساوية لـ KR‑20 لأنها تفترض تساوي صعوبة الفقرات؛ وكلما تباينت الصعوبات زاد الفرق.' },
    { t: 'التجزئة النصفية', q: 'قُسم اختبار إلى نصفين (فردي وزوجي)، فكان الارتباط بينهما 0.60، وتباين النصف الأول 4، والثاني 5، وتباين الاختبار الكلي 15. احسب الثبات بمعادلتي سبيرمان–براون وجتمان.', a: 'سبيرمان–براون: r = 2r½ ÷ (1 + r½) = 1.20 ÷ 1.60 = <b>0.75</b><br>جتمان: r = 2 × [1 − (σ²₁ + σ²₂) ÷ σ²ₓ] = 2 × (1 − 9 ÷ 15) = 2 × 0.40 = <b>0.80</b><br>جتمان لا تفترض تساوي تبايني النصفين، لذلك تُفضَّل عند اختلافهما.' },
    { t: 'معادلة التنبؤ بالطول', q: 'اختبار من 20 فقرة ثباته 0.70. (أ) ما ثباته إذا أصبح 40 فقرة متجانسة؟ (ب) كم فقرة نحتاج لرفع ثباته إلى 0.85؟', a: '(أ) n = 2: r = 2(0.70) ÷ (1 + 0.70) = 1.40 ÷ 1.70 = <b>0.82</b><br>(ب) n = r_new(1 − r) ÷ [r(1 − r_new)] = 0.85 × 0.30 ÷ (0.70 × 0.15) = 0.255 ÷ 0.105 = 2.43 ← عدد الفقرات ≈ 20 × 2.43 = <b>49 فقرة</b>.' },
    { t: 'نسبة صدق المحتوى CVR', q: 'عُرضت أربع فقرات على 10 محكمين، فكان عدد من اعتبرها «أساسية»: 9، 7، 10، 5. احسب CVR لكل فقرة، وحدد الفقرات المقبولة إذا كانت القيمة الحرجة لـ 10 محكمين 0.62 (Lawshe).', a: 'CVR = (ne − N/2) ÷ (N/2)، N/2 = 5<br>الفقرة 1: (9 − 5) ÷ 5 = <b>0.80</b> ✔ تُقبل<br>الفقرة 2: (7 − 5) ÷ 5 = <b>0.40</b> ✘ تُراجع<br>الفقرة 3: (10 − 5) ÷ 5 = <b>1.00</b> ✔ تُقبل<br>الفقرة 4: (5 − 5) ÷ 5 = <b>0.00</b> ✘ تُحذف<br>مؤشر صدق المحتوى CVI للفقرات المقبولة = (0.80 + 1.00) ÷ 2 = 0.90.<br><span class="muted">ملاحظة: الجدول المصحح (Ayre & Scally, 2014) يضع القيمة الحرجة لـ 10 محكمين عند 0.80.</span>' },
    { t: 'الصعوبة والتمييز', q: 'طُبق اختبار على 74 طالبًا، وأُخذت أعلى 27% وأدنى 27% (20 طالبًا لكل مجموعة). الفقرة (أ): أجاب عنها 16 من العليا و6 من الدنيا. الفقرة (ب): 10 من العليا و9 من الدنيا. احسب معاملي الصعوبة والتمييز لكل فقرة واحكم عليهما.', a: 'الفقرة (أ): p = (16 + 6) ÷ 40 = <b>0.55</b> (مناسبة)، D = (16 − 6) ÷ 20 = <b>0.50</b> (ممتاز) ← تُبقى.<br>الفقرة (ب): p = (10 + 9) ÷ 40 = <b>0.475</b> (مناسبة)، D = (10 − 9) ÷ 20 = <b>0.05</b> (ضعيف جدًا) ← تُحذف أو تُعاد صياغتها.' },
    { t: 'معامل كابا', q: 'صنّف مصححان 50 استجابة إلى (مقبول/مرفوض): اتفقا على «مقبول» في 22 وعلى «مرفوض» في 18، وقبل الأول ورفض الثاني 4، ورفض الأول وقبل الثاني 6. احسب معامل كابا كوهين وفسّره.', a: 'Po = (22 + 18) ÷ 50 = 0.80<br>الأول قَبِل 26، والثاني قَبِل 28 ← Pe = (26 × 28 + 24 × 22) ÷ 50² = (728 + 528) ÷ 2500 = 0.502<br>κ = (0.80 − 0.502) ÷ (1 − 0.502) = <b>0.60</b> ← اتفاق متوسط (Landis & Koch: .41–.60).' },
    { t: 'الدرجات المعيارية', q: 'متوسط درجات مقياس 28 وانحرافه المعياري 4. حصل أحمد على 34 وسارة على 22. حوّل الدرجتين إلى z وT والتساعية وفسّر.', a: 'أحمد: z = (34 − 28) ÷ 4 = <b>1.5</b>، T = 50 + 10(1.5) = <b>65</b>، التساعية = 5 + 2(1.5) = <b>8</b> ← مرتفع.<br>سارة: z = (22 − 28) ÷ 4 = <b>−1.5</b>، T = <b>35</b>، التساعية = <b>2</b> ← منخفض.' },
    { t: 'الخطأ المعياري للقياس', q: 'مقياس انحرافه المعياري 10 ومعامل ثباته 0.91. احسب الخطأ المعياري للقياس، ثم فترة الثقة 95% للدرجة الحقيقية لمفحوص درجته 60.', a: 'SEM = SD √(1 − r) = 10 × √0.09 = <b>3</b><br>فترة الثقة 95% = 60 ± 1.96 × 3 = 60 ± 5.88 ← من <b>54.1 إلى 65.9</b>.' }
  ];
  const ESSAY = [
    { q: 'تناول بالشرح خطوات بناء المقياس النفسي.', a: '<p>تمر عملية البناء باثنتي عشرة قاعدة في ثلاث مراحل:</p><p><b>أولًا: التخطيط</b></p><ol><li><b>الوضوح المفاهيمي:</b> تحديد السمة والفئة والغرض، وكتابة تعريف نظري يحدد حدودها، وتعريف إجرائي يربطها بدرجة المقياس.</li><li><b>المراجعة الأدبية:</b> مراجعة النظريات والمقاييس السابقة لاشتقاق الأبعاد وتجنب التكرار.</li><li><b>تحديد الأبعاد:</b> تقسيم السمة إلى مكونات مع تعريف كل بعد.</li><li><b>جدول المواصفات:</b> توزيع الفقرات على الأبعاد بأوزان نسبية تضمن تمثيل المحتوى.</li></ol><p><b>ثانيًا: البناء</b></p><ol start="5"><li><b>فيض الفقرات:</b> كتابة عدد أكبر من المطلوب (1.5–2 مرة) من مصادر استنباطية (النظرية) واستقرائية (أسئلة مفتوحة للعينة).</li><li><b>جودة الصياغة:</b> فكرة واحدة، لغة مناسبة، بلا نفي مزدوج، مع عبارات سالبة لضبط أسلوب الموافقة.</li><li><b>سلم الاستجابة:</b> اختيار التدريج (ليكرت خماسي مثلًا) وكتابة التعليمات.</li><li><b>صدق المحتوى:</b> عرض المقياس على محكمين وحساب CVR وCVI.</li></ol><p><b>ثالثًا: التجريب والتحقق</b></p><ol start="9"><li><b>التجريب الاستطلاعي:</b> على 30–50 فردًا لكشف الغموض وتقدير الزمن.</li><li><b>التحليل الإحصائي للفقرات:</b> الارتباط المصحح بالدرجة الكلية، والمقارنة الطرفية، والصعوبة والتمييز.</li><li><b>الثبات:</b> ألفا، والتجزئة النصفية، وإعادة الاختبار، وأوميغا.</li><li><b>الصدق العاملي والتقنين:</b> التحليل العاملي الاستكشافي ثم التوكيدي، واشتقاق المعايير (المئينية والتائية) من عينة تقنين ممثلة.</li></ol>' },
    { q: 'قارن بين طرق ليكرت وثرستون وجتمان في بناء المقاييس من حيث: طريقة البناء، والاستجابة، والتصحيح، والمزايا والعيوب.', a: '<div class="tw"><table><thead><tr><th>وجه المقارنة</th><th>ثرستون (1928)</th><th>ليكرت (1932)</th><th>جتمان (1944)</th></tr></thead><tbody><tr><td>البناء</td><td>عبارات كثيرة يقدرها محكمون (50–100) على متصل من 11 فئة، وتُختار ذات الاتفاق العالي والقيم الموزعة</td><td>عبارات موجبة وسالبة تُنتقى بتحليل الفقرات (الارتباط بالكلية، المقارنة الطرفية)</td><td>عبارات مرتبة تصاعديًا في الشدة بحيث تكون تراكمية</td></tr><tr><td>الاستجابة</td><td>يختار المفحوص العبارات التي يوافق عليها فقط</td><td>تدريج متعدد (خماسي غالبًا) لكل عبارة</td><td>نعم / لا</td></tr><tr><td>التصحيح</td><td>وسيط (أو متوسط) القيم المقياسية للعبارات المختارة</td><td>مجموع الدرجات بعد عكس السالبة</td><td>عدد العبارات المقبولة، ويُتحقق بمعامل إعادة الإنتاج ≥ .90</td></tr><tr><td>المزايا</td><td>فترات متساوية ظاهريًا؛ قيم العبارات محددة مسبقًا</td><td>سهلة البناء والتصحيح؛ ثبات مرتفع؛ تعطي شدة الاتجاه</td><td>أحادية البعد؛ الدرجة تعني نمط استجابة محددًا</td></tr><tr><td>العيوب</td><td>جهد ووقت كبيران؛ تأثر القيم باتجاهات المحكمين</td><td>الدرجة الواحدة قد تأتي من أنماط مختلفة؛ أسلوب الموافقة والحياد</td><td>صعوبة تحقيق التراكم في السمات المعقدة؛ عدد عبارات محدود</td></tr></tbody></table></div>' },
    { q: 'تناول بالشرح مفهوم الصدق وأنواعه وطرق تقدير كل نوع.', a: '<p><b>الصدق:</b> أن يقيس المقياس ما وُضع لقياسه، وهو أهم خصائص المقياس، ويتعلق بتفسير الدرجات واستخدامها.</p><ol><li><b>الصدق الظاهري:</b> مظهر الفقرات يوحي بقياس السمة؛ أضعف الأنواع، يفيد في تقبل المفحوصين.</li><li><b>صدق المحتوى:</b> تمثيل الفقرات لمجال السمة، ويُقدّر بجدول المواصفات وتحكيم الخبراء (CVR لـ Lawshe، وCVI، ومعامل أيكن V).</li><li><b>صدق المحك:</b> ارتباط المقياس بمحك خارجي؛ <i>تلازمي</i> (المحك في الوقت نفسه) أو <i>تنبؤي</i> (المحك لاحقًا).</li><li><b>صدق البناء (التكوين الفرضي):</b> مدى قياس المقياس للتكوين النظري، ويُقدّر بـ: التحليل العاملي (الاستكشافي والتوكيدي)، والصدق التقاربي والتمييزي (مصفوفة السمات المتعددة/الطرق المتعددة لكامبل وفيسك)، والمقارنة الطرفية (27%)، والاتساق الداخلي (ارتباط الفقرة بالبعد والكلية)، والفروق بين المجموعات المتضادة.</li></ol><p>العلاقة بالثبات: الثبات شرط لازم غير كافٍ للصدق، والحد الأقصى لمعامل الصدق = √(معامل الثبات).</p>' },
    { q: 'ما المقصود بثبات المقياس؟ تناول طرق تقديره والعوامل المؤثرة فيه.', a: '<p><b>الثبات:</b> اتساق الدرجات واستقرارها إذا أعيد القياس في ظروف مماثلة، أي خلو الدرجة من خطأ القياس العشوائي (نسبة التباين الحقيقي إلى التباين الكلي).</p><p><b>طرق تقديره:</b></p><ol><li><b>إعادة الاختبار:</b> تطبيقان بفاصل (أسبوعان) ← معامل الاستقرار.</li><li><b>الصور المتكافئة:</b> صورتان متكافئتان ← معامل التكافؤ.</li><li><b>التجزئة النصفية:</b> ارتباط النصفين مصححًا بسبيرمان–براون، أو معادلة جتمان/رولون.</li><li><b>الاتساق الداخلي:</b> ألفا كرونباخ، وKR‑20 وKR‑21 للفقرات الثنائية، وأوميغا ماكدونالد.</li><li><b>ثبات المصححين:</b> معامل الاتفاق وكابا كوهين.</li></ol><p><b>العوامل المؤثرة:</b> طول المقياس (يزيد الثبات بزيادته)، وتجانس العينة (يقل الثبات بتجانسها)، وصعوبة الفقرات (المتوسطة أفضل)، ووضوح التعليمات، وموضوعية التصحيح، والفاصل الزمني في إعادة الاختبار، وظروف التطبيق والتخمين.</p>' },
    { q: 'تناول بالشرح تحليل المفردات (الفقرات) ودوره في بناء المقاييس.', a: '<p>يهدف تحليل المفردات إلى انتقاء الفقرات الجيدة وحذف الضعيفة أو تعديلها، ويشمل:</p><ol><li><b>معامل الصعوبة p:</b> نسبة الإجابات الصحيحة؛ المقبول .20–.80 والأمثل نحو .50.</li><li><b>معامل التمييز D:</b> قدرة الفقرة على التفريق بين المرتفعين والمنخفضين: D = (صع − صد) ÷ ن، باستخدام 27% العليا والدنيا؛ ويحكم عليه بتصنيف إيبل (≥ .40 ممتاز، < .20 ضعيف).</li><li><b>الارتباط المصحح بالدرجة الكلية:</b> للفقرات المتدرجة؛ المقبول ≥ .30.</li><li><b>المقارنة الطرفية:</b> اختبار (ت) بين المجموعتين لكل فقرة.</li><li><b>فاعلية البدائل (المموهات):</b> في الاختيار من متعدد؛ البديل الجيد يجذب من الدنيا أكثر من العليا.</li><li><b>ألفا عند حذف الفقرة:</b> إذا ارتفع الثبات بحذفها فهي مرشحة للحذف.</li></ol>' },
    { q: 'وضّح الفرق بين القياس والتقييم والتقويم، واذكر مستويات القياس مع أمثلة.', a: '<p><b>القياس:</b> التقدير الكمي للسمة وفق قواعد (درجة الطالب 35). <b>التقييم:</b> إعطاء قيمة للدرجة بمقارنتها بمحك أو معيار (35 تعني «جيد»). <b>التقويم:</b> عملية أشمل تتضمن الحكم واتخاذ القرار والعلاج والتطوير. العلاقة هرمية: القياس ⊂ التقييم ⊂ التقويم.</p><p><b>مستويات القياس (ستيفنز):</b></p><ol><li><b>الاسمي:</b> تصنيف بلا ترتيب (النوع، التخصص).</li><li><b>الرتبي:</b> ترتيب دون تساوي المسافات (ترتيب الطلاب).</li><li><b>الفئوي (الفتري):</b> مسافات متساوية وصفر اعتباري (درجات الحرارة المئوية، ومعظم المقاييس النفسية).</li><li><b>النسبي:</b> صفر مطلق (الطول، الوزن، زمن الرجع).</li></ol>' },
    { q: 'قارن بين النظرية الكلاسيكية في القياس ونظرية الاستجابة للمفردة.', a: '<div class="tw"><table><thead><tr><th>وجه المقارنة</th><th>النظرية الكلاسيكية</th><th>نظرية الاستجابة للمفردة (IRT)</th></tr></thead><tbody><tr><td>النموذج</td><td>X = T + E (الدرجة الملاحظة = الحقيقية + الخطأ)</td><td>احتمال الإجابة دالة لوجستية في القدرة θ ومعالم المفردة</td></tr><tr><td>معالم المفردة</td><td>الصعوبة والتمييز تعتمد على العينة</td><td>الصعوبة b والتمييز a والتخمين c مستقلة عن العينة</td></tr><tr><td>تقدير الفرد</td><td>يعتمد على مجموعة الفقرات</td><td>مستقل عن مجموعة الفقرات المطبقة</td></tr><tr><td>خطأ القياس</td><td>واحد لكل الأفراد (SEM)</td><td>يختلف حسب مستوى القدرة (دالة المعلومات)</td></tr><tr><td>العينة</td><td>تكفي عينات صغيرة نسبيًا</td><td>تحتاج عينات كبيرة (250–500 فأكثر)</td></tr><tr><td>التطبيقات</td><td>بناء المقاييس التقليدية</td><td>بنوك الأسئلة، والاختبارات التكيفية، والمعادلة، وكشف DIF</td></tr></tbody></table></div><p>النماذج: أحادي المعلم (راش)، وثنائي المعلم، وثلاثي المعلم.</p>' },
    { q: 'تناول بالشرح المعايير وأنواعها وخطوات تقنين المقياس.', a: '<p><b>المعايير:</b> درجات مرجعية مشتقة من أداء عينة تقنين ممثلة، تُفسَّر في ضوئها الدرجة الخام لأنها لا معنى لها في ذاتها.</p><p><b>أنواعها:</b> المئينيات (الرتبة المئينية)، والدرجات المعيارية z، والتائية T = 50 + 10z، والتساعية (Stanine = 5 + 2z)، ونسب الذكاء الانحرافية (100 ± 15)، والمعايير العمرية والصفية.</p><p><b>خطوات التقنين:</b> تحديد المجتمع، واختيار عينة تقنين ممثلة كبيرة (عشوائية طبقية)، وتوحيد شروط التطبيق والتعليمات والزمن والتصحيح، وحساب الصدق والثبات على عينة التقنين، واشتقاق جداول المعايير لكل فئة (النوع، العمر)، وإعداد دليل المقياس.</p>' },
    { q: 'ما مصادر الخطأ في استجابات المقاييس النفسية (أساليب الاستجابة)؟ وكيف يمكن ضبطها؟', a: '<ol><li><b>المرغوبية الاجتماعية:</b> الإجابة بما يقبله المجتمع ← ضمان السرية، ومقياس للمرغوبية (مارلو–كراون)، والاختيار القسري، وصياغة محايدة.</li><li><b>أسلوب الموافقة (الإذعان):</b> الموافقة دون النظر للمحتوى ← تضمين عبارات سالبة وعكسها.</li><li><b>التطرف أو التوسط:</b> اختيار البدائل الطرفية أو المحايدة دائمًا ← تدريج زوجي بلا محايد عند الحاجة، وتوضيح البدائل.</li><li><b>التزييف والتمارض:</b> ← مقاييس الصدق والكذب، والاختيار القسري، واختبارات الحكم الموقفي.</li><li><b>الإهمال والاستجابة العشوائية:</b> ← فقرات تحقق من الانتباه، واختصار المقياس.</li></ol>' },
    { q: 'اذكر شروط صياغة العبارة الجيدة في المقاييس النفسية مع أمثلة.', a: '<ol><li>تتضمن فكرة واحدة (لا «أحب المذاكرة وأنجح بسهولة»).</li><li>واضحة ومناسبة لمستوى الفئة، وبلغة سليمة مختصرة.</li><li>تجنب النفي المزدوج (لا «لا أرفض عدم التأجيل»).</li><li>تقريرية بضمير المتكلم وليست سؤالًا.</li><li>تجنب الكلمات المطلقة (دائمًا، أبدًا، كل) والعامة الغامضة (كثيرًا، أحيانًا داخل العبارة).</li><li>لا توحي بالإجابة المرغوبة اجتماعيًا.</li><li>ترتبط بالبعد الذي وُضعت له وفق التعريف الإجرائي.</li><li>تضمين نسبة من العبارات السالبة (ربع إلى ثلث) لضبط أسلوب الموافقة.</li></ol>' }
  ];
  // model exam composed in the department's format (70 marks / 2 hours)
  const MODEL = {
    tf: [1, 4, 7, 9, 13, 17, 19, 22, 23, 26, 30, 32, 35, 39, 40],
    mcq: [3, 4, 9, 10, 15, 17, 19, 20, 22, 25, 27, 28, 30, 36, 39],
    prob: [0, 5, 6],
    essay: [0, 1, 2, 4]
  };
  const L4 = ['أ', 'ب', 'ج', 'د'];
  const shuf = a => { a = a.slice(); for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(Math.random() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; } return a; };
  const range = n => Array.from({ length: n }, (_, i) => i);
  let mode = null, run = null, timer = null;
  const clearTimer = () => { if (timer) clearInterval(timer); timer = null; };

  function home(el) {
    clearTimer(); run = null;
    el.innerHTML = `<div class="section-title"><h2>أسئلة موضوعية — بناء المقاييس</h2></div>
      <div class="notice small"><b>على نمط امتحانات الدبلوم الخاص — كلية التربية، جامعة المنيا</b>امتحان تحريري مدته ساعتان ودرجته 70: صح وخطأ (✓ / ✗)، واختيار من متعدد (أ، ب، ج، د)، ومسائل حسابية على بيانات معطاة، وأسئلة مقالية «تناول بالشرح» و«قارن» و«أجب عن عدد منها فقط».</div>
      <div class="kpis"><div class="tile"><small>صح وخطأ</small><b>${TF.length}</b></div><div class="tile"><small>اختيار من متعدد</small><b>${MCQ.length}</b></div><div class="tile"><small>مسائل محلولة</small><b>${PROB.length}</b></div><div class="tile"><small>أسئلة مقالية بإجاباتها</small><b>${ESSAY.length}</b></div>${S.exam && S.exam.best != null ? `<div class="tile hl"><small>أفضل درجة في الامتحان التجريبي</small><b>${S.exam.best}/40</b></div>` : ''}</div>
      <div class="grid2">
        <button class="linkcard" data-x="quick"><span class="ic">${ic('bolt')}</span><span><b>تدريب سريع</b><span>10 صح وخطأ + 10 اختيار من متعدد عشوائية، مع التصحيح والشرح فورًا.</span></span></button>
        <button class="linkcard" data-x="mock" style="--c:var(--bad);--s:var(--bad-soft)"><span class="ic">${ic('timer')}</span><span><b>امتحان تجريبي بالتوقيت</b><span>شكل الامتحان: 15 صح وخطأ + 15 اختيار من متعدد + مسألة، في 40 دقيقة، ثم النتيجة والتصحيح.</span></span></button>
        <button class="linkcard" data-x="model" style="--c:var(--good);--s:var(--good-soft)"><span class="ic">${ic('check')}</span><span><b>نموذج امتحان محلول</b><span>نموذج كامل من 70 درجة بنمط القسم، مع الإجابة النموذجية لكل سؤال، وقابل للطباعة.</span></span></button>
        <button class="linkcard" data-x="bank" style="--c:var(--violet);--s:var(--violet-soft)"><span class="ic">${ic('book')}</span><span><b>بنك الأسئلة كاملًا</b><span>كل الأسئلة مرتبة حسب النوع، والإجابة مخفية حتى تضغط.</span></span></button>
      </div>`;
    $$('[data-x]', el).forEach(b => b.onclick = () => { mode = b.dataset.x; rerender(); });
  }
  const back = () => `<button class="back" id="xBack">${ic('back')} الأسئلة الموضوعية</button>`;
  const bindBack = el => { $('#xBack', el).onclick = () => { mode = null; clearTimer(); rerender(); }; };

  function tfHTML(i, n, name, show) { const q = TF[i]; return `<fieldset class="fq" style="border:none;margin:0" data-tf="${i}"><legend class="qt"><span class="num muted">${n}.</span> ${esc(q[0])}</legend><div class="lik" style="--k:2">${[['✓ صح', 1], ['✗ خطأ', 0]].map(([l, v]) => `<label><input type="radio" name="${name}" value="${v}" ${show ? 'disabled' : ''}><span>${l}</span></label>`).join('')}</div><div class="xfb" hidden></div></fieldset>`; }
  function mcqHTML(i, n, name) { const q = MCQ[i]; return `<fieldset class="fq" style="border:none;margin:0" data-mc="${i}"><legend class="qt"><span class="num muted">${n}.</span> ${esc(q[0])}</legend><div class="lik opts" style="--k:1">${q[1].map((o, j) => `<label><input type="radio" name="${name}" value="${j}"><span>${L4[j]}) ${esc(o)}</span></label>`).join('')}</div><div class="xfb" hidden></div></fieldset>`; }
  function grade(root) {
    let right = 0, total = 0;
    $$('[data-tf]', root).forEach(f => { const i = +f.dataset.tf, q = TF[i], c = $('input:checked', f), ok = c && +c.value === q[1]; total++; if (ok) right++; mark(f, ok, `${q[1] ? '✓ العبارة صحيحة' : '✗ العبارة خاطئة'} — ${q[2]}`, c); });
    $$('[data-mc]', root).forEach(f => { const i = +f.dataset.mc, q = MCQ[i], c = $('input:checked', f), ok = c && +c.value === q[2]; total++; if (ok) right++; mark(f, ok, `الإجابة: (${L4[q[2]]}) ${q[1][q[2]]} — ${q[3]}`, c); });
    return { right, total };
  }
  function mark(f, ok, txt, c) { const fb = $('.xfb', f); fb.hidden = false; fb.className = 'xfb notice small ' + (ok ? 'good' : 'bad'); fb.innerHTML = `<b>${ok ? 'إجابة صحيحة' : c ? 'إجابة خاطئة' : 'لم تُجب'}</b>${esc(txt)}`; $$('input', f).forEach(x => x.disabled = true); }

  function quick(el) {
    if (!run) run = { tf: shuf(range(TF.length)).slice(0, 10), mc: shuf(range(MCQ.length)).slice(0, 10) };
    el.innerHTML = back() + `<div class="section-title"><h2>تدريب سريع</h2><button class="btn sm ghost" id="xNew">أسئلة جديدة</button></div>
      <form id="xF" class="stack" novalidate><div class="card stack"><h3>السؤال الأول: ضع علامة (✓) أمام العبارة الصحيحة وعلامة (✗) أمام العبارة الخاطئة</h3>${run.tf.map((i, n) => tfHTML(i, n + 1, 't' + n)).join('')}</div>
      <div class="card stack"><h3>السؤال الثاني: اختر الإجابة الصحيحة</h3>${run.mc.map((i, n) => mcqHTML(i, n + 1, 'm' + n)).join('')}</div>
      <div id="xRes"></div><button class="btn" type="submit">${ic('check')} صحّح إجاباتي</button></form>`;
    bindBack(el);
    $('#xNew', el).onclick = () => { run = null; rerender(); };
    $('#xF', el).onsubmit = e => { e.preventDefault(); const r = grade(el); finish(el, r, 20); e.target.querySelector('button[type=submit]').hidden = true; };
  }
  function finish(el, r, of, extra = '') {
    S.correct = (S.correct || 0) + r.right; if (r.right === r.total) S.perfect = (S.perfect || 0) + 1;
    addXP(Math.round(r.right / 2) + 5);
    const pct = r.right / r.total * 100;
    $('#xRes', el).innerHTML = `<div class="card row" style="gap:16px">${ring(pct, `${r.right}/${r.total}`, 'النتيجة')}<div class="grow"><h3 style="margin:0 0 4px">${pct >= 85 ? 'ممتاز! أنت مستعد للامتحان' : pct >= 65 ? 'جيد؛ راجع الأسئلة الخاطئة' : 'تحتاج مراجعة؛ اقرأ الشرح تحت كل سؤال'}</h3><p class="small muted" style="margin:0">التصحيح والشرح يظهران تحت كل سؤال. ${extra}</p></div></div>`;
    $('#xRes', el).scrollIntoView({ behavior: 'smooth', block: 'center' });
  }

  function mock(el) {
    if (!run) run = { tf: shuf(range(TF.length)).slice(0, 15), mc: shuf(range(MCQ.length)).slice(0, 15), p: Math.floor(Math.random() * PROB.length), end: Date.now() + 40 * 60000, done: false };
    el.innerHTML = back() + `<div class="card row" style="position:sticky;top:64px;z-index:5"><b class="grow">امتحان تجريبي — مقرر بناء المقاييس النفسية</b><span class="pill" id="xClock" aria-live="off">${ic('timer')} 40:00</span></div>
      <form id="xF" class="stack" novalidate><div class="card stack"><h3>السؤال الأول (15 درجة): ضع علامة (✓) أو (✗)</h3>${run.tf.map((i, n) => tfHTML(i, n + 1, 't' + n)).join('')}</div>
      <div class="card stack"><h3>السؤال الثاني (15 درجة): اختر الإجابة الصحيحة</h3>${run.mc.map((i, n) => mcqHTML(i, n + 1, 'm' + n)).join('')}</div>
      <div class="card stack"><h3>السؤال الثالث (10 درجات): ${esc(PROB[run.p].t)}</h3><p style="margin:0">${esc(PROB[run.p].q)}</p><textarea id="xP" placeholder="اكتب خطوات الحل هنا (للمراجعة الذاتية)"></textarea><div id="xPA" hidden></div></div>
      <div id="xRes"></div><button class="btn" type="submit">${ic('check')} إنهاء الامتحان وتصحيحه</button></form>`;
    bindBack(el);
    const sub = () => { if (run.done) return; run.done = true; clearTimer(); const r = grade(el); const pa = $('#xPA', el); pa.hidden = false; pa.innerHTML = `<div class="notice good small"><b>الحل النموذجي</b>${PROB[run.p].a}</div><p class="small muted">قيّم حلك بنفسك مقارنة بالحل النموذجي (من 10).</p>`; S.exam = S.exam || {}; S.exam.best = Math.max(S.exam.best || 0, r.right); S.exam.tries = (S.exam.tries || 0) + 1; finish(el, r, 30, `درجتك في الموضوعي ${r.right} من 30، ويضاف إليها تقديرك للمسألة من 10.`); $('#xF button[type=submit]', el).hidden = true; };
    $('#xF', el).onsubmit = e => { e.preventDefault(); sub(); };
    const tickC = () => { const c = $('#xClock'); if (!c) return clearTimer(); const ms = Math.max(0, run.end - Date.now()), m = Math.floor(ms / 60000), s2 = Math.floor(ms % 60000 / 1000); c.innerHTML = `${ic('timer')} ${String(m).padStart(2, '0')}:${String(s2).padStart(2, '0')}`; if (ms <= 0) { toast('انتهى الوقت'); sub(); } };
    if (!run.done) { tickC(); timer = setInterval(tickC, 1000); }
  }

  function model(el) {
    const M = MODEL;
    el.innerHTML = back() + `<div class="row"><span class="grow"></span><button class="btn sm soft" id="xShow">إظهار كل الإجابات</button><button class="btn sm ghost" id="xPrint">${ic('download')} طباعة / PDF</button></div>
      <article class="card stack exam-sheet" id="xSheet">
        <div class="tw"><table class="exam-head"><tbody><tr><td>جامعة المنيا<br>كلية التربية<br>قسم علم النفس التربوي</td><td>الدبلوم الخاص — علم النفس<br><b>امتحان مقرر: بناء المقاييس النفسية</b><br>(نموذج تدريبي)</td><td>الزمن: ساعتان<br>الدرجة الكلية: 70 درجة</td></tr></tbody></table></div>
        <p style="text-align:center;margin:0"><b>أجب عن الأسئلة الآتية:</b></p>
        <h3>السؤال الأول: ضع علامة (✓) أمام العبارة الصحيحة وعلامة (✗) أمام العبارة الخاطئة (15 درجة)</h3>
        <ol>${M.tf.map(i => `<li>${esc(TF[i][0])} <span class="ans" hidden><b>${TF[i][1] ? '(✓)' : '(✗)'}</b> ${esc(TF[i][2])}</span></li>`).join('')}</ol>
        <h3>السؤال الثاني: اختر الإجابة الصحيحة مما بين الأقواس (15 درجة)</h3>
        <ol>${M.mcq.map(i => `<li>${esc(MCQ[i][0])}<div class="small">${MCQ[i][1].map((o, j) => `${L4[j]}) ${esc(o)}`).join(' &nbsp; ')}</div><span class="ans" hidden><b>الإجابة: (${L4[MCQ[i][2]]})</b> ${esc(MCQ[i][3])}</span></li>`).join('')}</ol>
        <h3>السؤال الثالث: (10 درجات)</h3><p>${esc(PROB[M.prob[0]].q)}</p><div class="ans" hidden>${PROB[M.prob[0]].a}</div>
        <h3>السؤال الرابع: (15 درجة)</h3><p>(أ) ${esc(PROB[M.prob[1]].q)}</p><div class="ans" hidden>${PROB[M.prob[1]].a}</div><p>(ب) ${esc(PROB[M.prob[2]].q)}</p><div class="ans" hidden>${PROB[M.prob[2]].a}</div>
        <h3>السؤال الخامس: أجب عن سؤالين فقط مما يأتي (15 درجة)</h3>
        <ol>${M.essay.map(i => `<li>${esc(ESSAY[i].q)}<div class="ans" hidden>${ESSAY[i].a}</div></li>`).join('')}</ol>
        <p style="text-align:center" class="muted">انتهت الأسئلة، مع أطيب التمنيات بالتوفيق</p>
        <p class="small muted" style="text-align:center;margin:0">نموذج تدريبي أعده مختبر المقاييس النفسية على نمط امتحانات القسم للمذاكرة، وليس امتحانًا رسميًا صادرًا عن الكلية.</p>
      </article>`;
    bindBack(el);
    let on = false;
    $('#xShow', el).onclick = e => { on = !on; $$('.ans', el).forEach(x => x.hidden = !on); e.target.textContent = on ? 'إخفاء الإجابات' : 'إظهار كل الإجابات'; };
    $('#xPrint', el).onclick = () => window.print();
  }

  function bank(el) {
    const sec = (t, items) => `<details class="steps" open><summary>${t}</summary><div class="stack" style="gap:6px;margin-top:8px">${items}</div></details>`;
    const card = (q, a) => `<details class="card" style="padding:10px 14px"><summary style="cursor:pointer">${q}</summary><div class="small" style="margin-top:8px">${a}</div></details>`;
    el.innerHTML = back() + `<div class="section-title"><h2>بنك الأسئلة</h2></div>` +
      sec(`صح وخطأ (${TF.length})`, TF.map((q, i) => card(`${i + 1}. ${esc(q[0])}`, `<b>${q[1] ? '✓ صحيحة' : '✗ خاطئة'}</b> — ${esc(q[2])}`)).join('')) +
      sec(`اختيار من متعدد (${MCQ.length})`, MCQ.map((q, i) => card(`${i + 1}. ${esc(q[0])}<div class="small muted">${q[1].map((o, j) => `${L4[j]}) ${esc(o)}`).join(' · ')}</div>`, `<b>(${L4[q[2]]}) ${esc(q[1][q[2]])}</b> — ${esc(q[3])}`)).join('')) +
      sec(`مسائل حسابية (${PROB.length})`, PROB.map((p, i) => card(`${i + 1}. <b>${esc(p.t)}:</b> ${esc(p.q)}`, p.a)).join('')) +
      sec(`أسئلة مقالية (${ESSAY.length})`, ESSAY.map((p, i) => card(`${i + 1}. ${esc(p.q)}`, p.a)).join(''));
    bindBack(el);
  }

  VIEWS.exam = el => { ({ quick, mock, model, bank })[mode] ? ({ quick, mock, model, bank })[mode](el) : home(el); };
  return { clear: () => { mode = null; run = null; clearTimer(); }, TF, MCQ, PROB, ESSAY, MODEL };
})();

/* ===== Question bank: every question on the site in one place, filterable by type and topic, with search and a quiz
   on the current selection. Adds questions drawn from Majeed (2014), ch. 10, and the real-examples lesson. ===== */
const BANK = (() => {
  const B = 'مجيد، 2014';
  // book-based true/false: [statement, true?, explanation]
  const BTF = [
    ['في طريقة ثرستون تُحسب القيمة المقياسية للعبارة بوسيط تقديرات المحكمين.', 1, `القيمة المقياسية = الوسيط، والمدى الربيعي Q مؤشر غموض العبارة (${B}، ص 358–359).`],
    ['العبارة التي يكون مداها الربيعي كبيرًا في طريقة ثرستون تُختار لأنها واضحة.', 0, 'المدى الربيعي الكبير يعني اختلاف المحكمين، فالعبارة غامضة وتُستبعد.'],
    ['استخدم ثرستون وشيف 300 محكم لتصنيف 130 عبارة عن الاتجاه نحو الكنيسة في 11 فئة.', 1, `فئات من A إلى K، والفئة F محايدة (${B}، ص 356).`],
    ['يرى المفحوص الأوزان (القيم المقياسية) بجانب عبارات مقياس ثرستون عند التطبيق.', 0, 'الأوزان يعرفها الباحث فقط، ويختار المفحوص العبارات التي يوافق عليها.'],
    ['درجة المفحوص في مقياس ثرستون هي وسيط (أو متوسط) القيم المقياسية للعبارات التي وافق عليها.', 1, 'وليست عدد العبارات التي وافق عليها.'],
    ['أشارت الدراسات إلى أن عددًا قليلًا من المحكمين (نحو 15–50) يعطي قيمًا مقياسية قريبة من قيم العدد الكبير.', 1, `بلغ الارتباط بين قيم مجموعتين من المحكمين نحو 0.99 (${B}، ص 359–360).`],
    ['من مزايا طريقة ليكرت أنها أقل جهدًا ووقتًا من طريقة ثرستون.', 1, `ولا تحتاج مجموعة كبيرة من المحكمين (${B}، ص 366).`],
    ['في طريقة ليكرت تُختار العبارات بعد مقارنة استجابات المجموعتين العليا والدنيا لكل عبارة.', 1, 'باختبار (ت) بين أعلى 25% وأدنى 25%؛ والعبارة المميزة تبقى.'],
    ['يبلغ ثبات مقياس ليكرت المكون من 20 عبارة عادة 0.85 فأكثر.', 1, `(${B}، ص 365).`],
    ['مقياس بوجاردس للبعد الاجتماعي (1925) من أوائل المقاييس ذات التدرج التراكمي.', 1, 'سبع درجات من قبول الزواج حتى الاستبعاد من الوطن.'],
    ['وحدات مقياس بوجاردس متساوية المسافات تمامًا.', 0, `من عيوبه أن الوحدات ليست متساوية (${B}، ص 355).`],
    ['في سلم جتمان يمكن من معرفة الدرجة الكلية للمفحوص استنتاج نمط إجاباته على العبارات.', 1, 'وهذه خاصية التدرج التراكمي.'],
    ['يُشترط في سلم جتمان ألا يقل معامل القياسية (Scalability) عن 0.60.', 1, `مع معامل إعادة إنتاج لا يقل عن 0.90 (${B}، ص 370).`],
    ['نسبة الخطأ المقبولة في سلم جتمان تصل إلى 25%.', 0, 'لا تتجاوز 10%، أي معامل إعادة إنتاج ≥ 0.90.'],
    ['وضع أوسجود طريقة التمايز الدلالي عام 1957 لدراسة المعاني.', 1, 'في كتاب «قياس المعنى» مع سوسي وتاننباوم.'],
    ['تقيس طريقة التمايز الدلالي الجانب المعرفي (المعلومات) أكثر من الجانب الوجداني.', 0, 'تقيس الجانب الوجداني (المشاعر نحو المفهوم) أكثر.'],
    ['تراوح ارتباط درجات التمايز الدلالي بدرجات مقياسي ثرستون وليكرت بين 0.74 و0.82.', 1, `(${B}، ص 372).`],
    ['مقياس الاتجاه نحو المعوقين (ATDP) له صورتان متكافئتان A وB.', 1, 'في كل منهما 30 عبارة، ومعامل التكافؤ 0.83.'],
    ['الصور المتكافئة تُستخدم لتقدير الثبات بمعامل التكافؤ.', 1, 'وتسمح بإعادة القياس دون أثر التذكر.'],
    ['دراسة القيم لألبورت وفيرنون وليندزي تقيس ست قيم بطريقة الاختيار القسري والترتيب.', 1, 'النظرية والاقتصادية والجمالية والاجتماعية والسياسية والدينية.'],
    ['درجات الاختيار القسري التقليدية تصلح لمقارنة الأفراد بعضهم ببعض مثل درجات ليكرت.', 0, 'هي درجات إبساتية تقارن قيم الفرد ببعضها؛ والمقارنة بين الأفراد تحتاج نموذج Thurstonian IRT.'],
    ['مقياس روزنبرج لتقدير الذات مثال على مقياس ليكرت ذي العبارات الموجبة والسالبة.', 1, 'عشر عبارات، خمس منها سالبة تُعكس.'],
    ['في المقياس اللفظي للاتجاهات قد لا تتطابق الاستجابة اللفظية مع السلوك الفعلي.', 1, `ذكر الكتاب أمثلة على التعارض بين القول والسلوك (${B}، ص 349–350).`],
    ['طريقة المقارنة المزدوجة تعرض على المفحوص كل زوج من الموضوعات ليختار المفضل منهما.', 1, `(${B}، ص 352).`],
    ['الاختبار التكيفي المحوسب يعرض على كل المفحوصين الفقرات نفسها بالترتيب نفسه.', 0, 'يختار الفقرة التالية حسب إجابات الفرد، فيختلف ما يراه كل مفحوص.']
  ];
  // book-based MCQ: [stem, options, correct index, explanation]
  const BMCQ = [
    ['القيمة المقياسية للعبارة في طريقة ثرستون هي:', ['متوسط استجابات المفحوصين', 'وسيط تقديرات المحكمين', 'المدى الربيعي', 'معامل التمييز'], 1, 'الوسيط، والمدى الربيعي للغموض.'],
    ['عدد الفئات التي صنّف فيها المحكمون عبارات ثرستون:', ['5', '7', '9', '11'], 3, 'من A إلى K.'],
    ['مفحوص وافق على عبارات ثرستون ذات القيم 1.5 و2.2 و3.3؛ درجته:', ['3', '2.2', '7', '1.5'], 1, 'وسيط القيم = 2.2.'],
    ['في مقياس ليكرت خماسي تُختار العبارة إذا:', ['وافق عليها كل المحكمين', 'ميزت دالًا بين المجموعتين العليا والدنيا', 'كانت قيمتها المقياسية 6', 'كان مداها الربيعي كبيرًا'], 1, 'تحليل الفقرات بالمقارنة الطرفية.'],
    ['أول من قاس البعد الاجتماعي بمقياس تراكمي (1925):', ['ثرستون', 'ليكرت', 'بوجاردس', 'أوسجود'], 2, 'مقياس بوجاردس.'],
    ['في مقياس بوجاردس، الدرجة التي تمثل أقصى قرب اجتماعي:', ['قبول الزواج من فرد منهم', 'قبوله زائرًا', 'قبوله مواطنًا', 'استبعاده'], 0, 'الزواج = أقرب درجة.'],
    ['في تطبيق بوجاردس على 1725 أمريكيًا، أعلى قبول للزواج كان لـ:', ['الكوريين', 'البولنديين', 'السويديين', 'الإنجليز'], 3, 'الإنجليز 93.7%.'],
    ['الحد الأدنى لمعامل إعادة الإنتاج في سلم جتمان:', ['0.60', '0.75', '0.80', '0.90'], 3, 'أي أخطاء لا تتجاوز 10%.'],
    ['20 مفحوصًا × 5 عبارات جتمان، وعدد الأخطاء 7؛ معامل إعادة الإنتاج:', ['0.93', '0.65', '0.70', '0.97'], 0, '1 − 7/100.'],
    ['معامل القياسية في سلم جتمان يساوي:', ['نسبة التحسن ÷ نسبة التحسن المحتملة', 'عدد الأخطاء ÷ عدد الاستجابات', '1 − CR', 'CR × 2'], 0, 'ولا يقل عن 0.60.'],
    ['عدد نقاط تدريج التمايز الدلالي عادة:', ['3', '5', '7', '11'], 2, 'سباعي.'],
    ['زوج الصفات (قوي – ضعيف) في التمايز الدلالي ينتمي إلى بعد:', ['التقييم', 'القوة', 'النشاط', 'الفهم'], 1, 'Potency.'],
    ['زوج الصفات (سريع – بطيء) ينتمي إلى بعد:', ['التقييم', 'القوة', 'النشاط', 'الصدق'], 2, 'Activity.'],
    ['معامل التكافؤ بين صورتي مقياس الاتجاه نحو المعوقين:', ['0.41', '0.65', '0.83', '0.99'], 2, '(مجيد، 2014، ص 386).'],
    ['عدد القيم التي تقيسها دراسة القيم لألبورت وفيرنون وليندزي:', ['4', '5', '6', '8'], 2, 'تصنيف شبرنجر.'],
    ['نوع الدرجات الناتجة عن الاختيار القسري التقليدي:', ['معيارية', 'إبساتية (نسبية داخل الفرد)', 'نسبية بصفر مطلق', 'محكية'], 1, 'مجموعها ثابت لكل مفحوص.'],
    ['في مقياس روزنبرج (رباعي)، الدرجة المعكوسة لاستجابة «أوافق» (3) على عبارة سالبة:', ['1', '2', '3', '4'], 1, '5 − 3 = 2.'],
    ['مدى الدرجة الكلية في مقياس روزنبرج لتقدير الذات:', ['0–10', '10–40', '10–50', '1–4'], 1, '10 عبارات × (1 إلى 4).'],
    ['نموذج سلم التقدير لمقاييس ليكرت في عائلة راش وضعه:', ['أندريش', 'لوشي', 'كوهين', 'سبيرمان'], 0, 'Andrich (1978).'],
    ['بنوك فقرات PROMIS مثال حقيقي على:', ['طريقة ثرستون', 'الاختبار التكيفي المحوسب المعاير بنظرية الاستجابة للمفردة', 'التمايز الدلالي', 'سلم جتمان'], 1, 'يختار الحاسوب الفقرة التالية.'],
    ['يتوقف الاختبار التكيفي غالبًا عندما:', ['تنتهي كل الفقرات', 'يصبح الخطأ المعياري لتقدير السمة صغيرًا بما يكفي', 'يخطئ المفحوص مرة', 'تمر 10 دقائق'], 1, 'محك الدقة.'],
    ['مفتاح تصحيح اختبارات الحكم الموقفي يُحدَّد غالبًا بـ:', ['الصدفة', 'إجماع الخبراء', 'معامل ألفا', 'المتوسط'], 1, 'Expert consensus.'],
    ['الطريقة التي يرتب فيها المفحوص موضوعات أو ألوانًا حسب تفضيله:', ['التصويت', 'الترتيب', 'التقدير', 'جتمان'], 1, `(${B}، ص 351–352).`],
    ['أي مما يلي من أهداف قياس الاتجاهات التي ذكرها الكتاب؟', ['استطلاع الرأي العام', 'دراسة الإعلان والتسويق', 'تقويم برامج التربية', 'كل ما سبق'], 3, `(${B}، ص 348–349).`],
    ['الفقرات المولدة بالذكاء الاصطناعي في بحوث بناء المقاييس:', ['تُستخدم مباشرة دون مراجعة', 'تمر بمراجعة الخبراء والتحليل الإحصائي', 'تغني عن العينة', 'لا تحتاج صدقًا'], 1, 'Hommel et al. (2022)؛ Götz et al. (2024).']
  ];
  const TOPICS = [['concepts', 'المفاهيم والمستويات', /القياس|التقويم|التقييم|مستوى|الاسمي|الرتبي|النسبي|الفئوي/], ['steps', 'خطوات البناء والصياغة', /قاعد|خطو|صياغ|جدول المواصفات|فيض|العبار|التجريب|تعريف/], ['methods', 'طرق البناء', /ثرستون|ليكرت|جتمان|أوسجود|التمايز|بوجاردس|قسري|موقفي|ذكاء اصطناعي|المقارنة المزدوجة|الترتيب|ألبورت|روزنبرج/], ['validity', 'الصدق', /صدق|محك|CVR|CVI|أيكن|تقارب|تمييزي|بناء/], ['reliability', 'الثبات', /ثبات|ألفا|KR|كودر|سبيرمان|جتمان للتجزئة|تجزئة|إعادة|كابا|أوميغا|الخطأ المعياري/], ['items', 'تحليل الفقرات', /صعوب|تمييز|الفقرة|المفردة|طرفي|27%/], ['factor', 'التحليل العاملي', /عامل|KMO|بارتليت|تشبع|تدوير|فاريماكس|كايزر/], ['norms', 'المعايير والدرجات', /معيار|مئين|تائي|التساعي|z|T =|قطع|أنجوف/], ['irt', 'نظرية الاستجابة والحديثة', /راش|الاستجابة للمفردة|IRT|تكيفي|PROMIS|DIF|التفاضلي|لوجيت/]];
  const topicOf = t => { for (const [k, , re] of TOPICS) if (re.test(t)) return k; return 'steps'; };
  const TYPES = [['tf', 'صح وخطأ'], ['mcq', 'اختيار من متعدد'], ['prob', 'مسائل'], ['essay', 'مقالي'], ['review', 'مراجعة وفهم']];
  const L4 = ['أ', 'ب', 'ج', 'د'];
  let ALL = null;
  function build() {
    if (ALL) return ALL;
    const out = [], add = o => { o.id = 'q' + out.length; o.topic = o.topic || topicOf(o.q + ' ' + (o.w || '')); out.push(o); };
    (EXAM.TF || []).forEach(x => add({ type: 'tf', q: x[0], a: x[1] ? 1 : 0, w: x[2], src: 'نمط الامتحان' }));
    BTF.forEach(x => add({ type: 'tf', q: x[0], a: x[1] ? 1 : 0, w: x[2], src: 'كتاب مجيد' }));
    (EXAM.MCQ || []).forEach(x => add({ type: 'mcq', q: x[0], o: x[1], a: x[2], w: x[3], src: 'نمط الامتحان' }));
    BMCQ.forEach(x => add({ type: 'mcq', q: x[0], o: x[1], a: x[2], w: x[3], src: 'كتاب مجيد' }));
    (C.mcq || []).forEach(x => add({ type: 'mcq', q: x.q, o: x.o, a: x.a, w: x.w || '', src: 'الدروس' }));
    (EXAM.PROB || []).forEach(x => add({ type: 'prob', q: `${x.t}: ${x.q}`, ans: x.a, src: 'مسائل محلولة' }));
    (EXAM.ESSAY || []).forEach(x => add({ type: 'essay', q: x.q, ans: x.a, src: 'نمط الامتحان' }));
    (C.bank || []).forEach(x => add({ type: 'review', q: x.q, ans: esc(x.a), src: 'الدروس' }));
    return (ALL = out);
  }
  const st = { type: '', topic: '', q: '', page: 0 };
  const PER = 25;
  function filtered() { const s = st.q.trim(); return build().filter(x => (!st.type || x.type === st.type) && (!st.topic || x.topic === st.topic) && (!s || (x.q + ' ' + (x.o || []).join(' ')).includes(s))); }
  function card(x, n) {
    let body = '';
    if (x.type === 'tf') body = `<div class="small"><b>${x.a ? '✓ صحيحة' : '✗ خاطئة'}</b> — ${esc(x.w)}</div>`;
    else if (x.type === 'mcq') body = `<div class="small"><b>(${L4[x.a]}) ${esc(x.o[x.a])}</b>${x.w ? ' — ' + esc(x.w) : ''}</div>`;
    else body = `<div class="small">${x.ans}</div>`;
    const tname = (TYPES.find(t => t[0] === x.type) || [])[1], top = (TOPICS.find(t => t[0] === x.topic) || [])[1];
    return `<details class="card bq"><summary><span class="num muted">${n}.</span> ${esc(x.q)}${x.o ? `<div class="small muted">${x.o.map((o, j) => `${L4[j]}) ${esc(o)}`).join(' · ')}</div>` : ''}<div class="row" style="gap:6px;margin-top:6px"><span class="chip">${tname}</span><span class="chip violet">${top}</span><span class="chip">${esc(x.src)}</span></div></summary><div style="margin-top:8px">${body}</div></details>`;
  }
  VIEWS.bank = el => {
    const all = build(), list = filtered(), pages = Math.max(1, Math.ceil(list.length / PER)); st.page = Math.min(st.page, pages - 1);
    const cnt = k => all.filter(x => x.type === k).length, tcnt = k => all.filter(x => x.topic === k).length;
    el.innerHTML = `<div class="section-title"><h2>بنك الأسئلة</h2><span class="chip lime">${all.length} سؤال</span></div>
      <div class="notice small"><b>كل أسئلة المادة في مكان واحد</b>أسئلة الدروس، وأسئلة بنمط امتحان القسم، وأسئلة من كتاب «أسس بناء الاختبارات والمقاييس النفسية والتربوية» (مجيد، 2014)، ومسائل محلولة، وأسئلة مقالية بإجاباتها. اختر النوع والموضوع، ثم اضغط على السؤال لترى الإجابة، أو اختبر نفسك في الأسئلة المختارة.</div>
      <div class="seg" style="flex-wrap:wrap"><button data-ty="" aria-pressed="${!st.type}">الكل (${all.length})</button>${TYPES.map(t => `<button data-ty="${t[0]}" aria-pressed="${st.type === t[0]}">${t[1]} (${cnt(t[0])})</button>`).join('')}</div>
      <div class="seg" style="flex-wrap:wrap"><button data-tp="" aria-pressed="${!st.topic}">كل الموضوعات</button>${TOPICS.map(t => `<button data-tp="${t[0]}" aria-pressed="${st.topic === t[0]}">${t[1]} (${tcnt(t[0])})</button>`).join('')}</div>
      <div class="row"><input type="search" id="bqS" class="grow" placeholder="ابحث في الأسئلة… (مثلًا: ألفا، ثرستون، صدق المحتوى)" value="${esc(st.q)}" aria-label="بحث في بنك الأسئلة"><button class="btn lime" id="bqQuiz">${ic('bolt')} اختبرني في المختار</button>${typeof REWARD !== 'undefined' ? `<button class="btn ghost" id="bqRev">مراجعة اليوم (${REWARD.dueCount(all.map(x => 's:' + x.q.slice(0, 60)))})</button>` : ''}</div>
      <p class="small muted" style="margin:0">${list.length} سؤال مطابق${pages > 1 ? ` · صفحة ${st.page + 1} من ${pages}` : ''}</p>
      <div class="stack" style="gap:8px">${list.slice(st.page * PER, st.page * PER + PER).map((x, i) => card(x, st.page * PER + i + 1)).join('') || '<p class="muted">لا توجد أسئلة مطابقة.</p>'}</div>
      ${pages > 1 ? `<div class="row"><button class="btn ghost sm" id="bqPrev" ${st.page === 0 ? 'disabled' : ''}>السابق</button><span class="grow"></span><button class="btn ghost sm" id="bqNext" ${st.page >= pages - 1 ? 'disabled' : ''}>التالي</button></div>` : ''}`;
    $$('[data-ty]', el).forEach(b => b.onclick = () => { st.type = b.dataset.ty; st.page = 0; rerender(); });
    $$('[data-tp]', el).forEach(b => b.onclick = () => { st.topic = b.dataset.tp; st.page = 0; rerender(); });
    const s = $('#bqS', el); let tm; s.oninput = () => { clearTimeout(tm); tm = setTimeout(() => { st.q = s.value; st.page = 0; rerender(); setTimeout(() => { const n = $('#bqS'); if (n) { n.focus(); n.setSelectionRange(n.value.length, n.value.length); } }, 0); }, 350); };
    const pv = $('#bqPrev', el), nx = $('#bqNext', el); if (pv) pv.onclick = () => { st.page--; rerender(); }; if (nx) nx.onclick = () => { st.page++; rerender(); };
    $('#bqQuiz', el).onclick = () => quiz(list.filter(x => x.type === 'tf' || x.type === 'mcq')); const rv = $('#bqRev', el); if (rv) rv.onclick = () => { const due = all.filter(x => (x.type === 'tf' || x.type === 'mcq') && REWARD.isDue('s:' + x.q.slice(0, 60))); if (!due.length) return toast('مفيش أسئلة مستحقة المراجعة النهارده. الأسئلة اللي بتغلط فيها بترجع بعد يوم، ثم 3، ثم 7، ثم 14 يوم.'); quiz(due, 1); };
  };
  function quiz(pool, min = 3) {
    if (pool.length < min) return toast('اختر نوعًا أو موضوعًا فيه 3 أسئلة موضوعية على الأقل');
    const qs = pool.slice().sort(() => Math.random() - 0.5).slice(0, 10); let i = 0, right = 0;
    modal('<div id="bqM"></div>', (bg, close) => {
      const box = $('#bqM', bg);
      const show = () => {
        if (i >= qs.length) { const pct = right / qs.length * 100; S.correct = (S.correct || 0) + right; if (right === qs.length) { S.perfect = (S.perfect || 0) + 1; if (typeof MOTION !== 'undefined') MOTION.confetti(); } addXP(right + 3); box.innerHTML = `<h3>النتيجة</h3><div class="row" style="gap:16px">${ring(pct, right + '/' + qs.length, 'النتيجة')}<p class="grow">${pct >= 80 ? 'ممتاز! استمر.' : 'راجع الشرح تحت كل سؤال في البنك.'}</p></div><div class="row"><button class="btn" id="bqAgain">اختبار جديد</button><button class="btn ghost" data-close>إغلاق</button></div>`; $('#bqAgain', box).onclick = () => { close(); quiz(pool); }; $$('[data-close]', box).forEach(b => b.onclick = close); return; }
        const x = qs[i], opts = x.type === 'tf' ? ['✓ صح', '✗ خطأ'] : x.o;
        box.innerHTML = `<div class="row"><span class="chip">${i + 1} / ${qs.length}</span><span class="grow"></span><button class="ibtn" data-close aria-label="إغلاق">✕</button></div><h3 style="margin:8px 0">${esc(x.q)}</h3><div class="opts">${opts.map((o, j) => `<button class="opt" data-j="${j}"><span class="k">${x.type === 'tf' ? '' : L4[j]}</span>${esc(o)}</button>`).join('')}</div><div id="bqF"></div>`;
        $$('[data-close]', box).forEach(b => b.onclick = close);
        $$('[data-j]', box).forEach(b => b.onclick = () => {
          const j = +b.dataset.j, correct = x.type === 'tf' ? (x.a ? 0 : 1) : x.a, ok = j === correct; if (ok) right++; if (typeof REWARD !== 'undefined') REWARD.onAnswer('s:' + x.q.slice(0, 60), ok);
          $$('[data-j]', box).forEach(o => { o.disabled = true; if (+o.dataset.j === correct) o.style.borderColor = 'var(--good)'; });
          if (!ok) b.style.borderColor = 'var(--bad)';
          $('#bqF', box).innerHTML = `<div class="notice small ${ok ? 'good' : 'bad'}"><b>${ok ? 'إجابة صحيحة' : 'إجابة خاطئة'}</b>${esc(x.w || '')}</div><button class="btn" id="bqN">${i + 1 < qs.length ? 'التالي' : 'النتيجة'}</button>`;
          $('#bqN', box).onclick = () => { i++; show(); };
        });
      };
      show();
    });
  }
  return { count: () => build().length };
})();

/* Destrieux 2009 cortical atlas (the parcellation used on TRIBE v2's fsaverage5 output): Arabic name, lobe and a
   cautious description of commonly reported functions. Functions are tendencies from the neuroimaging literature,
   not one-to-one mappings ("reverse inference" is weak). Subcortical structures (amygdala, hippocampus, striatum,
   thalamus) are not on the cortical surface and are therefore never predicted by TRIBE v2. */
const DESTRIEUX = (() => {
  const F = 'الجبهي', P = 'الجداري', T = 'الصدغي', O = 'القذالي', I = 'الجزيرة', L = 'الحوفي (الحزامي)', C = 'المركزي (الحسي–الحركي)';
  const R = {
    'G_and_S_frontomargin': ['التلفيف والتلم الجبهي الهامشي', F, 'جزء من القشرة الجبهية الأمامية؛ التخطيط والتقييم'],
    'G_and_S_occipital_inf': ['التلفيف والتلم القذالي السفلي', O, 'معالجة بصرية'],
    'G_and_S_paracentral': ['الفصيص المجاور للمركزي', C, 'حركة وإحساس الطرف السفلي'],
    'G_and_S_subcentral': ['التلفيف والتلم تحت المركزي', C, 'حركة الوجه والفم؛ إنتاج الكلام'],
    'G_and_S_transv_frontopol': ['القطب الجبهي المستعرض', F, 'التفكير المجرد وتنسيق الأهداف المتعددة'],
    'G_and_S_cingul-Ant': ['الحزامي الأمامي', L, 'كشف التعارض، والانفعال، والألم الانفعالي، وتنظيم الانتباه'],
    'G_and_S_cingul-Mid-Ant': ['الحزامي الأوسط الأمامي', L, 'الضبط المعرفي والتحكم في الاستجابة'],
    'G_and_S_cingul-Mid-Post': ['الحزامي الأوسط الخلفي', L, 'التوجه الحركي والانتباه للجسم'],
    'G_cingul-Post-dorsal': ['الحزامي الخلفي الظهري', L, 'شبكة الوضع الافتراضي؛ التفكير الذاتي'],
    'G_cingul-Post-ventral': ['الحزامي الخلفي البطني', L, 'الذاكرة الذاتية والمعالجة الانفعالية'],
    'G_cuneus': ['الوتد (الكونيوس)', O, 'معالجة بصرية أولية'],
    'G_front_inf-Opercular': ['التلفيف الجبهي السفلي – الجزء الغطائي', F, 'إنتاج اللغة (منطقة بروكا) والتحكم في الكلام'],
    'G_front_inf-Orbital': ['التلفيف الجبهي السفلي – الجزء الحجاجي', F, 'معاني الكلمات والتقييم الانفعالي'],
    'G_front_inf-Triangul': ['التلفيف الجبهي السفلي – الجزء المثلثي', F, 'معالجة المعنى في اللغة (منطقة بروكا)'],
    'G_front_middle': ['التلفيف الجبهي الأوسط', F, 'الذاكرة العاملة والتخطيط والوظائف التنفيذية'],
    'G_front_sup': ['التلفيف الجبهي العلوي', F, 'الوعي بالذات، والتحكم التنفيذي، والتفكير في الآخرين'],
    'G_Ins_lg_and_S_cent_ins': ['الجزيرة الطويلة', I, 'الإحساس الداخلي بالجسم (الحشوي)'],
    'G_insular_short': ['الجزيرة القصيرة (الأمامية)', I, 'الوعي الانفعالي، والقلق، وإدراك حالة الجسم'],
    'G_occipital_middle': ['التلفيف القذالي الأوسط', O, 'معالجة بصرية'],
    'G_occipital_sup': ['التلفيف القذالي العلوي', O, 'معالجة بصرية مكانية'],
    'G_oc-temp_lat-fusifor': ['التلفيف المغزلي', T, 'التعرف على الوجوه والكلمات المكتوبة'],
    'G_oc-temp_med-Lingual': ['التلفيف اللساني', O, 'معالجة بصرية والذاكرة البصرية'],
    'G_oc-temp_med-Parahip': ['التلفيف المجاور للحصين', T, 'الذاكرة والسياق المكاني'],
    'G_orbital': ['التلفيف الحجاجي', F, 'تقييم المكافأة والقيمة الانفعالية واتخاذ القرار'],
    'G_pariet_inf-Angular': ['التلفيف الزاوي', P, 'المعنى، والقراءة، وشبكة الوضع الافتراضي'],
    'G_pariet_inf-Supramar': ['التلفيف فوق الهامشي', P, 'الذاكرة الصوتية، والتعاطف، وتمييز الذات عن الآخر'],
    'G_parietal_sup': ['الفصيص الجداري العلوي', P, 'الانتباه المكاني والذاكرة العاملة'],
    'G_postcentral': ['التلفيف خلف المركزي', C, 'الإحساس الجسدي'],
    'G_precentral': ['التلفيف أمام المركزي', C, 'الحركة الإرادية'],
    'G_precuneus': ['مقدم الوتد', P, 'التفكير الذاتي والذكريات الشخصية والتخيل'],
    'G_rectus': ['التلفيف المستقيم', F, 'التقييم الانفعالي والاجتماعي'],
    'G_subcallosal': ['التلفيف تحت الجسم الثفني', L, 'تنظيم المزاج؛ مرتبط في البحوث بالاكتئاب'],
    'G_temp_sup-G_T_transv': ['تلفيف هيشل (السمعي الأولي)', T, 'السمع الأولي'],
    'G_temp_sup-Lateral': ['التلفيف الصدغي العلوي – الوحشي', T, 'فهم الكلام المسموع والإدراك الاجتماعي'],
    'G_temp_sup-Plan_polar': ['المستوى القطبي الصدغي', T, 'معالجة سمعية وصوتية'],
    'G_temp_sup-Plan_tempo': ['المستوى الصدغي', T, 'معالجة أصوات الكلام (منطقة فيرنيكه)'],
    'G_temporal_inf': ['التلفيف الصدغي السفلي', T, 'التعرف على الأشياء والمعاني'],
    'G_temporal_middle': ['التلفيف الصدغي الأوسط', T, 'معاني الكلمات وفهم اللغة'],
    'Lat_Fis-ant-Horizont': ['الشق الجانبي – الفرع الأفقي الأمامي', F, 'جزء من منطقة اللغة الجبهية'],
    'Lat_Fis-ant-Vertical': ['الشق الجانبي – الفرع العمودي الأمامي', F, 'جزء من منطقة اللغة الجبهية'],
    'Lat_Fis-post': ['الشق الجانبي – الفرع الخلفي', T, 'معالجة سمعية ولغوية'],
    'Pole_occipital': ['القطب القذالي', O, 'الإبصار الأولي'],
    'Pole_temporal': ['القطب الصدغي', T, 'المعرفة الدلالية والاجتماعية والانفعالية'],
    'S_calcarine': ['التلم المهمازي', O, 'القشرة البصرية الأولية'],
    'S_central': ['التلم المركزي', C, 'الحد بين الحركة والإحساس'],
    'S_cingul-Marginalis': ['التلم الحزامي الهامشي', L, 'جزء من الشبكة الحزامية'],
    'S_circular_insula_ant': ['التلم الدائري للجزيرة – الأمامي', I, 'الوعي الانفعالي وحالة الجسم'],
    'S_circular_insula_inf': ['التلم الدائري للجزيرة – السفلي', I, 'معالجة سمعية وحشوية'],
    'S_circular_insula_sup': ['التلم الدائري للجزيرة – العلوي', I, 'الإحساس الداخلي بالجسم'],
    'S_collat_transv_ant': ['التلم الجانبي المستعرض الأمامي', T, 'الذاكرة والتعرف البصري'],
    'S_collat_transv_post': ['التلم الجانبي المستعرض الخلفي', O, 'معالجة بصرية'],
    'S_front_inf': ['التلم الجبهي السفلي', F, 'الضبط المعرفي واللغة'],
    'S_front_middle': ['التلم الجبهي الأوسط', F, 'الوظائف التنفيذية'],
    'S_front_sup': ['التلم الجبهي العلوي', F, 'الانتباه والذاكرة العاملة'],
    'S_interm_prim-Jensen': ['تلم جنسن', P, 'جزء من الفص الجداري السفلي'],
    'S_intrapariet_and_P_trans': ['التلم داخل الجداري', P, 'الانتباه، والأرقام، والتحكم في حركة العين'],
    'S_oc_middle_and_Lunatus': ['التلم القذالي الأوسط والهلالي', O, 'معالجة بصرية'],
    'S_oc_sup_and_transversal': ['التلم القذالي العلوي والمستعرض', O, 'معالجة بصرية مكانية'],
    'S_occipital_ant': ['التلم القذالي الأمامي', O, 'معالجة بصرية للحركة'],
    'S_oc-temp_lat': ['التلم القذالي الصدغي الوحشي', T, 'التعرف البصري'],
    'S_oc-temp_med_and_Lingual': ['التلم القذالي الصدغي الإنسي واللساني', T, 'الذاكرة البصرية'],
    'S_orbital_lateral': ['التلم الحجاجي الوحشي', F, 'التقييم واتخاذ القرار'],
    'S_orbital_med-olfact': ['التلم الحجاجي الإنسي (الشمي)', F, 'الشم والتقييم'],
    'S_orbital-H_Shaped': ['التلم الحجاجي على شكل H', F, 'تقييم المكافأة'],
    'S_parieto_occipital': ['التلم الجداري القذالي', P, 'الحد بين الجداري والقذالي'],
    'S_pericallosal': ['التلم حول الجسم الثفني', L, 'جزء من الشبكة الحزامية'],
    'S_postcentral': ['التلم خلف المركزي', P, 'معالجة حسية جسدية'],
    'S_precentral-inf-part': ['التلم أمام المركزي – السفلي', C, 'التحكم الحركي في الكلام'],
    'S_precentral-sup-part': ['التلم أمام المركزي – العلوي', C, 'التحكم الحركي وحركة العين'],
    'S_suborbital': ['التلم تحت الحجاجي', F, 'التقييم الانفعالي'],
    'S_subparietal': ['التلم تحت الجداري', P, 'شبكة الوضع الافتراضي'],
    'S_temporal_inf': ['التلم الصدغي السفلي', T, 'التعرف البصري والمعاني'],
    'S_temporal_sup': ['التلم الصدغي العلوي', T, 'فهم الكلام، وقراءة نوايا الآخرين، والإدراك الاجتماعي'],
    'S_temporal_transverse': ['التلم الصدغي المستعرض', T, 'السمع']
  };
  const LOBES = [F, C, P, T, O, I, L];
  function info(col) { // col like "L_G_front_sup" or "R_S_central" or "region_12"
    const m = String(col).match(/^([LR])_(.+)$/); const hemi = m ? (m[1] === 'L' ? 'الأيسر' : 'الأيمن') : '', key = m ? m[2] : col, r = R[key];
    return r ? { key, hemi, ar: r[0], lobe: r[1], fn: r[2], known: true } : { key, hemi, ar: col, lobe: '—', fn: '', known: false };
  }
  return { R, LOBES, info, keys: Object.keys(R) };
})();

/* ===== Neuro lab: TRIBE v2 (Meta FAIR, 2026) used where it is defensible for scale construction.
   TRIBE v2 predicts an average subject's fMRI response (fsaverage5 cortical mesh) to video, audio or text. It does not
   read anyone's mind and cannot run in a browser (GPU + gated LLaMA 3.2 text encoder; licence CC BY-NC 4.0).
   Use here: an exploratory "predicted neural coherence" check — do items written for the same dimension evoke more
   similar predicted cortical profiles than items from different dimensions? The site exports the items with a Python
   script; the researcher runs it on a GPU and imports the per-region CSV; all analysis below runs in the browser. ===== */
const NEURO = (() => {
  let data = null, srcScale = null;
  const pearson = (a, b) => { const n = a.length; let ma = 0, mb = 0; for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; } ma /= n; mb /= n; let sab = 0, sa = 0, sb = 0; for (let i = 0; i < n; i++) { const x = a[i] - ma, y = b[i] - mb; sab += x * y; sa += x * x; sb += y * y; } return sa && sb ? sab / Math.sqrt(sa * sb) : 0; };
  function script(items) {
    return `# -*- coding: utf-8 -*-
# Running TRIBE v2 on scale items — Psychological Scales Lab
# Requirements: a GPU with 24 GB or more (a Colab T4 is not enough), acceptance of the LLaMA 3.2 licence on Hugging Face,
# and tribev2 installed from https://github.com/facebookresearch/tribev2  (licence CC BY-NC 4.0: non-commercial use only).
# Check the function names against the repository README; this script follows the documented flow:
# TribeModel.from_pretrained -> get_events_dataframe -> predict  (predictions: timesteps x vertices on fsaverage5).
# Language: the text encoder does not officially support Arabic. For each item, fill in "text_en" (an accurate translation)
# and set USE_ENGLISH = True, or compare both runs.
import json, os, tempfile
import numpy as np, pandas as pd
from nilearn import datasets
from tribev2 import TribeModel  # adjust if the README uses a different import path

USE_ENGLISH = False
ITEMS = json.loads(r'''${JSON.stringify(items, null, 1)}''')

model = TribeModel.from_pretrained("facebook/tribev2")
atlas = datasets.fetch_atlas_surf_destrieux()          # fsaverage5: 10242 vertices per hemisphere
names = [n.decode() if isinstance(n, bytes) else str(n) for n in atlas["labels"]]
lab = np.concatenate([atlas["map_left"], np.asarray(atlas["map_right"]) + 1000])

rows = []
for it in ITEMS:
    text = (it.get("text_en") or it["text"]) if USE_ENGLISH else it["text"]
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as f:
        f.write(text); path = f.name
    events = model.get_events_dataframe(text_path=path)
    preds, segments = model.predict(events=events)
    v = np.asarray(preds).mean(axis=0)                  # average predicted response over time
    row = {"item_id": it["id"], "dim": it["dim"]}
    for code in np.unique(lab):
        nm = names[code % 1000]
        if nm in ("Unknown", "Medial_wall"):
            continue
        row[("R_" if code >= 1000 else "L_") + nm] = float(v[lab == code].mean())
    rows.append(row); os.unlink(path)
    print("done:", it["id"])

pd.DataFrame(rows).to_csv("tribe_items_regions.csv", index=False, encoding="utf-8")
print("Upload tribe_items_regions.csv to the Neuro lab on the site")
`;
  }
  /* ---- running TRIBE: Colab notebook and a Hugging Face Space (ZeroGPU) that the site can call ---- */
  const CORE_PY = `def item_regions(model, text):
    """Average predicted cortical response to one item, aggregated to the 148 Destrieux regions."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as f:
        f.write(text); path = f.name
    try:
        events = model.get_events_dataframe(text_path=path)
        preds, segments = model.predict(events=events)
    finally:
        os.unlink(path)
    v = np.asarray(preds).mean(axis=0)
    out = {}
    for code in np.unique(LAB):
        nm = NAMES[code % 1000]
        if nm in ("Unknown", "Medial_wall"):
            continue
        out[("R_" if code >= 1000 else "L_") + nm] = float(v[LAB == code].mean())
    return out`;
  const ATLAS_PY = `atlas = datasets.fetch_atlas_surf_destrieux()          # fsaverage5: 10242 vertices per hemisphere
NAMES = [n.decode() if isinstance(n, bytes) else str(n) for n in atlas["labels"]]
LAB = np.concatenate([np.asarray(atlas["map_left"]), np.asarray(atlas["map_right"]) + 1000])`;
  function notebook(items) {
    const cell = (type, src) => ({ cell_type: type, metadata: {}, source: src.split('\n').map((l, i, a) => i < a.length - 1 ? l + '\n' : l), ...(type === 'code' ? { execution_count: null, outputs: [] } : {}) });
    const nb = { nbformat: 4, nbformat_minor: 5, metadata: { accelerator: 'GPU', colab: { gpuType: 'A100' }, kernelspec: { name: 'python3', display_name: 'Python 3' } }, cells: [
      cell('markdown', '# TRIBE v2 on scale items\nRuntime > Change runtime type > **A100** or **L4** (a T4 runs out of memory). Then Runtime > Run all.\n\nBefore running: accept the LLaMA 3.2 licence at huggingface.co/meta-llama/Llama-3.2-3B with your Hugging Face account, and create a read token.\n\nLicence of TRIBE v2: CC BY-NC 4.0 (non-commercial research and teaching).'),
      cell('code', '!pip -q install git+https://github.com/facebookresearch/tribev2 nilearn pandas'),
      cell('code', 'from huggingface_hub import notebook_login\nnotebook_login()   # paste your read token'),
      cell('code', 'import json, os, tempfile\nimport numpy as np, pandas as pd\nfrom nilearn import datasets\nfrom tribev2 import TribeModel\n\nUSE_ENGLISH = False   # True: use the "text_en" field (fill it with accurate translations)\nITEMS = json.loads(r\'\'\'' + JSON.stringify(items, null, 1) + '\'\'\')\n\nmodel = TribeModel.from_pretrained("facebook/tribev2")\n' + ATLAS_PY + '\n\n' + CORE_PY),
      cell('code', 'rows = []\nfor it in ITEMS:\n    text = (it.get("text_en") or it["text"]) if USE_ENGLISH else it["text"]\n    rows.append({"item_id": it["id"], "dim": it["dim"], **item_regions(model, text)})\n    print("done:", it["id"])\npd.DataFrame(rows).to_csv("tribe_items_regions.csv", index=False, encoding="utf-8")'),
      cell('code', 'from google.colab import files\nfiles.download("tribe_items_regions.csv")   # then upload it in the site\'s Neuro lab')] };
    return JSON.stringify(nb, null, 1);
  }
  const SPACE = {
    'app.py': `# Hugging Face Space (Gradio + ZeroGPU) that runs TRIBE v2 on scale items for the Psychological Scales Lab.
# Setup: Space SDK = Gradio, hardware = ZeroGPU (needs a PRO account), add a secret HF_TOKEN whose account accepted the
# LLaMA 3.2 licence. Non-commercial use only (CC BY-NC 4.0). Check names against the tribev2 README if the API changes.
import os, json, tempfile
import numpy as np, pandas as pd, gradio as gr, spaces
from huggingface_hub import login
from nilearn import datasets
from tribev2 import TribeModel

if os.environ.get("HF_TOKEN"):
    login(token=os.environ["HF_TOKEN"])
model = TribeModel.from_pretrained("facebook/tribev2")
${ATLAS_PY}
MAX_ITEMS, MAX_CHARS = 60, 400

${CORE_PY}

@spaces.GPU(duration=180)
def run_batch(texts):
    return [item_regions(model, t) for t in texts]

def regions(items_json):
    try:
        items = json.loads(items_json)
    except Exception:
        raise gr.Error("items must be a JSON list")
    if not isinstance(items, list) or not 1 <= len(items) <= MAX_ITEMS:
        raise gr.Error(f"send 1-{MAX_ITEMS} items")
    texts = []
    for it in items:
        t = str(it.get("text_en") or it.get("text") or "")[:MAX_CHARS]
        if not t.strip():
            raise gr.Error("empty item text")
        texts.append(t)
    rows = [{"item_id": str(it.get("id")), "dim": str(it.get("dim", "")), **r} for it, r in zip(items, run_batch(texts))]
    return pd.DataFrame(rows).to_csv(index=False)

demo = gr.Interface(fn=regions, inputs=gr.Textbox(label="items JSON", lines=8), outputs=gr.Textbox(label="regions CSV"),
                    title="TRIBE v2 - scale items to Destrieux regions", api_name="regions")
demo.queue(max_size=8).launch()
`,
    'requirements.txt': 'git+https://github.com/facebookresearch/tribev2\nnilearn\npandas\nnumpy\nspaces\n',
    'README.md': '---\ntitle: TRIBE Scale Items\nsdk: gradio\napp_file: app.py\nlicense: cc-by-nc-4.0\n---\nRuns TRIBE v2 (Meta FAIR) on questionnaire items and returns mean predicted activation per Destrieux cortical region.\n'
  };
  async function callSpace(base, items, onMsg) { // Gradio HTTP API: POST /gradio_api/call/<api>, then read the event stream
    base = base.trim().replace(/\/+$/, '');
    if (!/^https:\/\/[a-z0-9-]+\.hf\.space$/i.test(base)) throw new Error('اكتب رابط الخدمة بالشكل https://username-spacename.hf.space');
    onMsg('إرسال العبارات…');
    const r = await fetch(base + '/gradio_api/call/regions', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ data: [JSON.stringify(items)] }) });
    if (!r.ok) throw new Error('لم تقبل الخدمة الطلب (' + r.status + '). تأكد أنها تعمل وأنها عامة.');
    const { event_id } = await r.json(); onMsg('النموذج يعمل على البطاقة الرسومية… (أول تشغيل قد يستغرق دقائق)');
    const txt = await (await fetch(base + '/gradio_api/call/regions/' + event_id)).text();
    const blocks = txt.split(/\n\n/).map(b => ({ ev: (b.match(/^event: *(.+)$/m) || [])[1], data: (b.match(/^data: *(.*)$/m) || [])[1] }));
    const done = blocks.find(b => b.ev === 'complete'), err = blocks.find(b => b.ev === 'error');
    if (!done) throw new Error('توقفت الخدمة' + (err && err.data && err.data !== 'null' ? ': ' + err.data : ' (ربما انتهت حصة الـ GPU المجانية؛ جرّب لاحقًا).'));
    return JSON.parse(done.data)[0];
  }
  function parseCSV(txt) {
    const lines = txt.replace(/\r/g, '').split('\n').filter(l => l.trim()); if (lines.length < 3) throw new Error('الملف قصير');
    const head = lines[0].split(',').map(x => x.trim()); const iI = head.indexOf('item_id'), iD = head.indexOf('dim');
    if (iI < 0 || iD < 0) throw new Error('يلزم عمودا item_id و dim');
    const regs = head.map((h, i) => i).filter(i => i !== iI && i !== iD); if (regs.length < 5) throw new Error('عدد المناطق أقل من 5');
    const items = lines.slice(1).map(l => { const c = l.split(','); return { id: c[iI], dim: c[iD], v: regs.map(i => +c[i]) }; }).filter(x => x.v.every(isFinite));
    if (items.length < 4) throw new Error('عدد العبارات أقل من 4');
    return { items, regions: regs.map(i => head[i]) };
  }
  function simulate() { // illustrative only: real atlas names, a shared language response, dimension-specific emphases
    const regions = ['L_', 'R_'].flatMap(h => DESTRIEUX.keys.map(k => h + k));
    const lang = ['G_temp_sup-Lateral', 'S_temporal_sup', 'G_temporal_middle', 'G_front_inf-Triangul', 'G_front_inf-Opercular', 'G_pariet_inf-Angular', 'G_temp_sup-Plan_tempo'];
    const prof = { 'التأجيل السلوكي': ['G_precuneus', 'G_cingul-Post-dorsal', 'G_front_sup', 'G_pariet_inf-Angular', 'S_subparietal'], 'إدارة الوقت': ['G_front_middle', 'S_front_inf', 'G_parietal_sup', 'S_intrapariet_and_P_trans', 'S_front_sup'], 'القلق والتجنب': ['G_insular_short', 'G_and_S_cingul-Ant', 'S_circular_insula_ant', 'G_orbital', 'Pole_temporal', 'G_subcallosal'] };
    const w = (keys, col, left) => keys.includes(col.slice(2)) ? (col[0] === 'L' ? left : 1) : 0;
    const shared = regions.map(c => 1.6 * w(lang, c, 1.4) + (Math.random() - 0.5) * 0.6);
    const items = []; Object.keys(prof).forEach((d, k) => { const base = regions.map(c => 1.3 * w(prof[d], c, 1) + (Math.random() - 0.5) * 0.5); for (let i = 0; i < 4; i++) items.push({ id: 'sim' + k + i, dim: d, v: base.map((x, r) => 0.8 * x + shared[r] + (Math.random() - 0.5) * (i === 3 && k === 1 ? 2.2 : 0.7)) }); });
    return { items, regions, sim: true };
  }
  /* What is distinctive about each dimension: z of the dimension mean against all items, region by region.
     Raw activation is dominated by the language network for every sentence; the contrast is what tells dimensions apart. */
  function regionMap(D) {
    const it = D.items, R = D.regions.length, mu = [], sd = [];
    for (let r = 0; r < R; r++) { const x = it.map(i => i.v[r]), m = x.reduce((a, b) => a + b, 0) / x.length; mu[r] = m; sd[r] = Math.sqrt(x.reduce((a, b) => a + (b - m) ** 2, 0) / Math.max(1, x.length - 1)) || 1; }
    const z = v => v.map((x, r) => (x - mu[r]) / sd[r]);
    const rawTop = mu.map((x, r) => [r, x]).sort((a, b) => b[1] - a[1]).slice(0, 6).map(x => x[0]);
    const dims = [...new Set(it.map(x => x.dim))];
    const byDim = dims.map(d => { const m = it.filter(x => x.dim === d), zz = z(D.regions.map((_, r) => m.reduce((s, x) => s + x.v[r], 0) / m.length));
      const top = zz.map((x, r) => [r, x]).sort((a, b) => b[1] - a[1]).slice(0, 6), lobes = {};
      zz.forEach((x, r) => { const L = DESTRIEUX.info(D.regions[r]).lobe; (lobes[L] = lobes[L] || []).push(x); });
      return { d, top, lobe: Object.fromEntries(Object.entries(lobes).map(([k, a]) => [k, a.reduce((s, x) => s + x, 0) / a.length])) }; });
    const byItem = it.map(x => z(x.v).map((v, r) => [r, v]).sort((a, b) => b[1] - a[1]).slice(0, 3).map(p => p[0]));
    return { rawTop, byDim, byItem, known: D.regions.filter(c => DESTRIEUX.info(c).known).length };
  }
  const regName = c => { const f = DESTRIEUX.info(c); return f.known ? `${esc(f.ar)} <span class="muted">(${f.hemi})</span>` : esc(c); };
  function brainSVG(lobe) { // schematic lateral view; colour = mean z of the lobe for this dimension
    const ok = v => typeof v === 'number' && isFinite(v);
    const col = v => !ok(v) ? 'transparent' : v >= 0 ? `rgba(91,63,209,${Math.min(0.9, 0.12 + v * 1.2).toFixed(2)})` : `rgba(242,100,60,${Math.min(0.9, 0.12 - v * 1.2).toFixed(2)})`;
    const P = [['الجبهي', 'M30 118 C18 70 52 26 112 20 L128 22 L120 96 L74 128 Z', 72, 72], ['المركزي (الحسي–الحركي)', 'M128 22 L150 24 L142 104 L120 96 Z', 135, 62], ['الجداري', 'M150 24 C200 26 236 46 246 82 L196 98 L142 104 Z', 194, 64], ['القذالي', 'M246 82 C262 104 258 132 236 146 L204 128 L196 98 Z', 230, 116], ['الصدغي', 'M74 128 L120 96 L142 104 L196 98 L204 128 L236 146 C200 166 124 168 84 150 Z', 150, 136]];
    const g = P.map(([n, d, x, y]) => `<path d="${d}" fill="${col(lobe[n])}" stroke="currentColor" stroke-opacity=".55" stroke-width="1.2"><title>${n}: ${ok(lobe[n]) ? fmt(lobe[n]) : '—'}</title></path><text x="${x}" y="${y}" font-size="10" text-anchor="middle" fill="currentColor">${n.split(' ')[0]}</text>`).join('');
    const extra = [['الجزيرة', 'الجزيرة (داخل الشق الجانبي)', 182], ['الحوفي (الحزامي)', 'الحزامي (السطح الداخلي)', 182]].map(([k, t], i) => `<circle cx="${i ? 124 : 262}" cy="178" r="6" fill="${col(lobe[k])}" stroke="currentColor" stroke-opacity=".55"><title>${t}: ${ok(lobe[k]) ? fmt(lobe[k]) : '—'}</title></circle><text x="${i ? 114 : 252}" y="182" font-size="9" fill="currentColor" direction="rtl" text-anchor="start">${t}</text>`).join('');
    return `<svg viewBox="0 0 270 192" style="width:100%;max-width:300px;height:auto" role="img" aria-label="خريطة تقريبية لفصوص المخ (منظر جانبي)">${g}${extra}</svg>`;
  }
  function analyse(D) {
    const it = D.items, n = it.length, S = it.map(a => it.map(b => pearson(a.v, b.v)));
    const coh = lab => { let w = 0, nw = 0, bt = 0, nb = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { if (lab[i] === lab[j]) { w += S[i][j]; nw++; } else { bt += S[i][j]; nb++; } } return { w: nw ? w / nw : NaN, b: nb ? bt / nb : NaN }; };
    const labs = it.map(x => x.dim), obs = coh(labs), nci = obs.w - obs.b;
    let ge = 0; const P = 2000; for (let p = 0; p < P; p++) { const sh = labs.slice().sort(() => Math.random() - 0.5), c = coh(sh); if (c.w - c.b >= nci) ge++; }
    const dims = [...new Set(labs)];
    const centroid = d => { const m = it.filter(x => x.dim === d); return m[0].v.map((_, r) => m.reduce((s, x) => s + x.v[r], 0) / m.length); };
    const C = Object.fromEntries(dims.map(d => [d, centroid(d)]));
    const fit = it.map(x => { const own = (() => { const m = it.filter(y => y.dim === x.dim && y !== x); if (!m.length) return NaN; const c = x.v.map((_, r) => m.reduce((s, y) => s + y.v[r], 0) / m.length); return pearson(x.v, c); })(); const other = dims.filter(d => d !== x.dim).map(d => [d, pearson(x.v, C[d])]).sort((a, b) => b[1] - a[1])[0] || ['', NaN]; return { own, other }; });
    const perDim = dims.map(d => { const idx = labs.map((l, i) => l === d ? i : -1).filter(i => i >= 0); let s = 0, c = 0; idx.forEach((i, a) => idx.slice(a + 1).forEach(j => { s += S[i][j]; c++; })); return [d, idx.length, c ? s / c : NaN]; });
    return { S, obs, nci, p: (ge + 1) / (P + 1), fit, perDim, dims };
  }
  function mantel(S, E) { // correlation between two similarity matrices, permutation p
    const n = S.length, iu = []; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) iu.push([i, j]);
    const r0 = pearson(iu.map(([i, j]) => S[i][j]), iu.map(([i, j]) => E[i][j])); let ge = 0; const P = 2000;
    for (let p = 0; p < P; p++) { const pi = Array.from({ length: n }, (_, i) => i).sort(() => Math.random() - 0.5); if (pearson(iu.map(([i, j]) => S[pi[i]][pi[j]]), iu.map(([i, j]) => E[i][j])) >= r0) ge++; }
    return { r: r0, p: (ge + 1) / (P + 1) };
  }
  function heat(S, labels) {
    const n = S.length, c = Math.max(10, Math.min(28, Math.floor(420 / n)));
    const col = v => { const t = Math.max(-1, Math.min(1, v)); return t >= 0 ? `rgba(91,63,209,${0.12 + 0.85 * t})` : `rgba(242,100,60,${0.12 + 0.85 * -t})`; };
    let g = ''; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) g += `<rect x="${j * c}" y="${i * c}" width="${c - 1}" height="${c - 1}" fill="${col(S[i][j])}"><title>${esc(labels[i])} × ${esc(labels[j])}: ${fmt(S[i][j])}</title></rect>`;
    return `<svg viewBox="0 0 ${n * c} ${n * c}" style="width:min(100%,${n * c}px);height:auto" role="img" aria-label="مصفوفة تشابه الأنماط العصبية المتوقعة بين العبارات">${g}</svg>`;
  }
  function results(el) {
    const A = analyse(data), M = regionMap(data), it = data.items, dimOf = id => (it.find(x => x.id === id) || {}).dim;
    const label = x => { if (srcScale) { const s = SCALES.get(srcScale), item = s && s.items.find(i => i.id === x.id); if (item) return (s.items.indexOf(item) + 1) + '. ' + item.text; } return x.id; };
    let emp = '';
    if (srcScale && !data.sim && typeof BUILD !== 'undefined') { const s = SCALES.get(srcScale), An = s && BUILD.analyze(s); if (An && An.rows && An.rows.length >= 20) { const idx = it.map(x => An.items.findIndex(i => i.id === x.id)); if (idx.every(i => i >= 0)) { const cols = idx.map(j => An.rows.map(r => r[j])), E = cols.map(a => cols.map(b => pearson(a, b))), m = mantel(A.S, E); emp = `<div class="notice small"><b>مقارنة بالبيانات الفعلية</b>الارتباط بين مصفوفة التشابه العصبي المتوقع ومصفوفة الارتباط الفعلية بين العبارات (اختبار مانتل): r = ${fmt(m.r)}، p = ${fmtp(m.p)}. ${m.p < 0.05 ? 'العبارات المتشابهة في نمطها العصبي المتوقع ترتبط أكثر في استجابات المفحوصين.' : 'لم يظهر تطابق دال بين التنبؤ العصبي والبيانات الفعلية.'}</div>`; } } }
    el.innerHTML = `${data.sim ? '<div class="notice warn small"><b>بيانات محاكاة للتدريب</b>هذه الأرقام مولّدة لتجربة التحليل، وليست مخرجات TRIBE v2.</div>' : ''}
      <div class="kpis"><div class="tile hl"><small>مؤشر التماسك العصبي المتوقع</small><b>${fmt(A.nci, 3)}</b></div><div class="tile"><small>التشابه داخل الأبعاد</small><b>${fmt(A.obs.w, 3)}</b></div><div class="tile"><small>التشابه بين الأبعاد</small><b>${fmt(A.obs.b, 3)}</b></div><div class="tile ${A.p < 0.05 ? 'good' : ''}"><small>اختبار التباديل (2000)</small><b>${fmtp(A.p)}</b></div></div>
      <div class="verdict ${A.p < 0.05 ? 'good' : 'warn'}">${A.p < 0.05 ? 'عبارات البعد الواحد تُحدث أنماطًا عصبية متوقعة أكثر تشابهًا من عبارات الأبعاد الأخرى (دال).' : 'لا يميز التنبؤ العصبي بين الأبعاد تمييزًا دالًا؛ هذا لا ينفي صدق المقياس، لكنه لا يدعمه من هذه الزاوية.'}</div>
      ${emp}
      <div class="grid2"><section class="stack"><h4 style="margin:0">مصفوفة التشابه بين العبارات</h4>${heat(A.S, it.map(label))}<p class="small muted">بنفسجي = تشابه موجب، برتقالي = سالب. العبارات مرتبة حسب البعد.</p></section>
      <section class="stack"><h4 style="margin:0">التماسك داخل كل بعد</h4>${tableHTML(['البعد', 'العبارات', 'متوسط التشابه'], A.perDim.map(d => [esc(d[0]), N(d[1], 0), N(d[2], 3)]))}</section></div>
      <h3 style="margin:8px 0 0">أماكن النشاط المتوقع في المخ</h3>
      ${M.known < data.regions.length / 2 ? '<div class="notice warn small">أسماء الأعمدة ليست من أطلس Destrieux، فتظهر المناطق بأسمائها كما هي من غير ترجمة.</div>' : ''}
      <div class="notice small"><b>كيف تُقرأ؟</b>أي جملة تُنشّط شبكة اللغة أكثر من غيرها، فأعلى المناطق نشاطًا خامًا متشابهة في كل الأبعاد: ${M.rawTop.map(r => regName(data.regions[r])).join('، ')}. لذلك نعرض لكل بعد <b>المناطق التي تميّزه عن باقي العبارات</b> (درجة z)، وليس أعلى نشاط خام.</div>
      <div class="stack">${M.byDim.map(o => `<section class="card stack"><h4 style="margin:0">${esc(o.d)}</h4><div style="display:flex;flex-wrap:wrap;gap:14px;align-items:flex-start"><div style="flex:0 1 260px">${brainSVG(o.lobe)}</div><div style="flex:1 1 360px;min-width:0">${tableHTML(['المنطقة', 'الفص', 'z', 'وظائف شائعة في البحوث'], o.top.map(([r, v]) => { const f = DESTRIEUX.info(data.regions[r]); return [regName(data.regions[r]), esc(f.lobe), N(v, 2), esc(f.fn || '—')]; }))}</div></div></section>`).join('')}</div>
      <p class="small muted">بنفسجي = نشاط متوقع أعلى من متوسط العبارات، برتقالي = أقل. الوظائف المذكورة «ميول» من الأدبيات وليست تشخيصًا؛ فالقول بأن «نشاط الجزيرة = قلق» استدلال عكسي ضعيف، لأن المنطقة الواحدة تشارك في وظائف كثيرة.</p>
      <h4 style="margin:0">ملاءمة كل عبارة لبعدها، ومناطقها المميزة</h4>${tableHTML(['العبارة', 'البعد المقصود', 'التشابه مع بعدها', 'أقرب بعد آخر', 'التشابه معه', 'أبرز 3 مناطق', 'الحكم'], it.map((x, i) => { const f = A.fit[i], bad = isFinite(f.own) && f.other[1] > f.own; return { cls: bad ? 'flag' : '', cells: [esc(label(x)), esc(x.dim), N(f.own, 3), esc(f.other[0]), N(f.other[1], 3), '<span class="small">' + M.byItem[i].map(r => regName(data.regions[r])).join('<br>') + '</span>', bad ? '<span class="chip warn">أقرب لبعد آخر — راجع صياغتها</span>' : '<span class="chip good">منسجمة</span>'] }; }))}
      <div class="notice warn small"><b>حدود الدقة المكانية</b>TRIBE v2 يتنبأ على <b>سطح القشرة فقط</b> (حوالي 20 ألف نقطة على fsaverage5) لمخ «متوسط». لذلك لا تظهر هنا البنى العميقة مثل اللوزة والحُصين والمخطط والمهاد، رغم أهميتها في القلق والانفعال والذاكرة. وكل منطقة هنا متوسط مئات النقاط، فالدقة «منطقة تشريحية» وليست «نقطة».</div>
      <p class="small muted">التفسير الآمن: هذا مؤشر <b>استكشافي</b> مبني على تنبؤ لمخ «متوسط»، وليس قياسًا لأدمغة المفحوصين ولا دليل صدق مستقلًا. استخدمه لاكتشاف عبارات تستحق المراجعة، ثم احسم بالتحليل العاملي والتحكيم.</p>`;
  }
  VIEWS.neuro = el => {
    const L = (typeof SCALES !== 'undefined' ? SCALES.list() : []).filter(s => (s.items || []).filter(i => i.keep !== false).length >= 4);
    el.innerHTML = `<div class="section-title"><h2>المختبر العصبي · TRIBE v2</h2><span class="chip">ناشئ · استكشافي</span></div>
      <section class="card stack"><h3 style="margin:0">ما هو TRIBE v2؟</h3>
        <p style="margin:0">نموذج أساس أطلقه فريق Meta FAIR في مارس 2026، يتنبأ باستجابة الرنين المغناطيسي الوظيفي (fMRI) على سطح القشرة المخية لمخ «متوسط» عند مشاهدة فيديو أو سماع صوت أو قراءة نص. دُرّب على أكثر من 1000 ساعة من تسجيلات fMRI لـ 720 شخصًا.</p>
        ${tableHTML(['يستطيع', 'لا يستطيع'], [['التنبؤ بالنمط المتوقع لنشاط القشرة لمثير معين', 'قراءة أفكار شخص أو مشاعره'], ['مقارنة مثيرات أو عبارات افتراضيًا قبل التجربة', 'الحلول محل تجربة على بشر'], ['اقتراح فرضيات عن معالجة اللغة والانفعال', 'تشخيص اضطراب أو قياس سمة لدى فرد']])}
        <p class="small muted" style="margin:0">الرخصة CC BY-NC 4.0 (استخدام بحثي وتعليمي غير تجاري). يتطلب بطاقة رسومية قوية (24 جيجابايت فأكثر) وقبول رخصة نموذج LLaMA 3.2، لذلك لا يعمل داخل المتصفح. ومُشفّر النص لا يدعم العربية رسميًا، فقارن النتائج بترجمة إنجليزية دقيقة للعبارات.</p></section>
      <section class="card stack"><h3 style="margin:0">${ic('download')} الخطوة 1: صدّر عبارات مقياسك وشغّل النموذج</h3>
        <p class="small muted" style="margin:0">الفكرة: لو كانت عبارات البعد الواحد تقيس شيئًا واحدًا، فالمتوقع أن تُحدث أنماطًا عصبية متوقعة أكثر تشابهًا من عبارات الأبعاد الأخرى.</p>
        ${L.length ? `<div class="row"><select id="nScale" aria-label="المقياس">${L.map(s => `<option value="${s.id}">${esc(s.name)} (${s.items.filter(i => i.keep !== false).length} عبارة)</option>`).join('')}</select><button class="btn sm" id="nPy">${ic('download')} سكربت Python + العبارات</button><button class="btn sm ghost" id="nJs">ملف العبارات JSON</button></div>` : '<p class="muted small" style="margin:0">ابنِ مقياسًا فيه 4 عبارات على الأقل في «ابنِ مقياسك» أولًا.</p>'}
        <p class="small" style="margin:0">الناتج ملف <code>tribe_items_regions.csv</code>: لكل عبارة متوسط النشاط المتوقع في 148 منطقة قشرية (أطلس Destrieux، 74 في كل نصف كرة)، ثم ترفعه في الخطوة 2. اختر طريقة التشغيل:</p>
        ${tableHTML(['الطريقة', 'التكلفة', 'السهولة'], [['<b>أ. Google Colab</b>: افتح الدفتر واختر A100 أو L4 ثم «تشغيل الكل»', 'Colab Pro تقريبًا 10 دولار/شهر', 'سهلة جدًا'], ['<b>ب. خدمة Hugging Face (ZeroGPU)</b>: ترفع الملفات الثلاثة مرة واحدة، والموقع يتصل بها مباشرة', 'اشتراك PRO تقريبًا 9 دولار/شهر', 'إعداد مرة واحدة، بعدها زرار'], ['<b>ج. خادم الجامعة أو Kaggle</b>: سكربت Python', 'مجاني غالبًا (Kaggle 16 جيجا قد لا تكفي)', 'متوسطة']])}
        ${L.length ? `<div class="row"><button class="btn sm" id="nNb">${ic('download')} دفتر Colab</button><button class="btn sm ghost" id="nSp">ملفات خدمة Hugging Face</button></div>` : ''}
        <details class="small"><summary>خطوات خدمة Hugging Face (مرة واحدة)</summary><ol style="margin:6px 0 0"><li>اشترك في Hugging Face PRO، واقبل رخصة LLaMA 3.2 من صفحة meta-llama/Llama-3.2-3B.</li><li>أنشئ Space جديدة: SDK = Gradio، العتاد = ZeroGPU، واجعلها عامة.</li><li>ارفع app.py وrequirements.txt وREADME.md، وأضف في الإعدادات «Secret» اسمه HF_TOKEN بتوكن قراءة.</li><li>انسخ رابط الخدمة (https://username-spacename.hf.space) والصقه بالأسفل.</li></ol></details>
        ${L.length ? `<div class="row"><input id="nSpace" dir="ltr" placeholder="https://username-tribe.hf.space" style="flex:1;min-width:220px" aria-label="رابط خدمة Hugging Face"><button class="btn sm" id="nRun">شغّل TRIBE على عبارات المقياس</button></div><div id="nRunMsg" class="small muted"></div>` : ''}
        <p class="small muted" style="margin:0">العرض التجريبي الرسمي من Meta (aidemos.atmeta.com/tribev2) يعرض خرائط للمخ للمشاهدة، لكنه لا يصدّر أرقامًا لكل منطقة، فلا يصلح للتحليل هنا.</p></section>
      <section class="card stack"><h3 style="margin:0">${ic('upload')} الخطوة 2: حلّل النتائج</h3>
        <div class="row"><input type="file" id="nCsv" accept=".csv,text/csv"><button class="btn sm ghost" id="nSim">جرّب ببيانات محاكاة</button></div><div id="nOut"></div></section>`;
    const sel = $('#nScale', el), items = () => { const s = SCALES.get(sel.value); return s.items.filter(i => i.keep !== false).map(i => ({ id: i.id, dim: (s.dims.find(d => d.id === i.dim) || {}).name || '', text: i.text, text_en: '' })); };
    if (sel) {
      $('#nNb', el).onclick = () => downloadFile('tribe_items.ipynb', notebook(items()), 'application/x-ipynb+json');
      $('#nSp', el).onclick = async () => { for (const [n, t] of Object.entries(SPACE)) await downloadFile(n, t); };
      const sp = $('#nSpace', el), msg = $('#nRunMsg', el); try { sp.value = localStorage.getItem('tribeSpace') || ''; } catch (e) {}
      $('#nRun', el).onclick = async () => { const b = $('#nRun', el); b.disabled = true; try { try { localStorage.setItem('tribeSpace', sp.value.trim()); } catch (e) {} const csv = await callSpace(sp.value, items(), t => { msg.textContent = t; }); data = parseCSV(csv); srcScale = sel.value; msg.textContent = 'تم. النتائج بالأسفل.'; results($('#nOut', el)); } catch (err) { msg.textContent = err.message === 'Failed to fetch' ? 'تعذّر الاتصال بالخدمة (قد تكون نائمة أو محجوبة من هذه الصفحة). جرّب من موقع GitHub أو انتظر دقيقة.' : err.message; } finally { b.disabled = false; } };
    }
    if (sel) { $('#nPy', el).onclick = () => downloadFile('tribe_items.py', script(items()), 'text/x-python;charset=utf-8'); $('#nJs', el).onclick = () => downloadFile('tribe_items.json', JSON.stringify(items(), null, 1), 'application/json'); }
    $('#nCsv', el).onchange = e => { const f = e.target.files[0]; if (!f) return; const rd = new FileReader(); rd.onload = () => { try { data = parseCSV(rd.result); srcScale = sel ? sel.value : null; results($('#nOut', el)); } catch (err) { $('#nOut', el).innerHTML = `<div class="notice bad small">${esc(err.message)}</div>`; } }; rd.readAsText(f); };
    $('#nSim', el).onclick = () => { data = simulate(); srcScale = null; results($('#nOut', el)); };
  };
  return {};
})();

/* ===== «اعرف نفسك» — an optional, game-like self-assessment journey (shared by both labs).
   Stages: consent with a comprehension check → mood (affect grid) → working memory (digit span) with a hidden
   notification probe → executive control (Stroop) → wellbeing screeners (WHO-5, PHQ-8, GAD-7) → M6 (always-on digital
   pressure) → M2 (cognitive dependence on AI) with a trust-calibration task → results, debrief and safety information.
   Screeners are presented as cards but keep their item content and scoring. Results are a screening profile for learning
   and referral, never a diagnosis. Data stay in this browser unless the person exports them. ===== */
var KNOW = (function () {
  var KEY = 'know-v1', st = load();
  function load() { try { return JSON.parse(localStorage.getItem(KEY) || 'null') || { v: 1, done: {}, r: {} }; } catch (e) { return { v: 1, done: {}, r: {} }; } }
  function save() { try { localStorage.setItem(KEY, JSON.stringify(st)); } catch (e) { } }
  var E = function (s) { return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]; }); };
  var xp = function (n) { try { if (typeof addXP === 'function') addXP(n); } catch (e) { } };
  var say = function (m) { try { if (typeof toast === 'function') return toast(m); } catch (e) { } };
  var root = null;
  var STAGES = [['consent', 'البداية والموافقة'], ['mood', 'المزاج'], ['span', 'الذاكرة'], ['stroop', 'التركيز'], ['well', 'الحالة النفسية'], ['m6', 'حياتك الرقمية'], ['m2', 'أنت والذكاء الاصطناعي'], ['res', 'النتائج']];
  var CSS = '.kn{max-width:880px;margin:0 auto;display:flex;flex-direction:column;gap:16px}.kn-map{display:flex;flex-wrap:wrap;gap:6px}.kn-node{border:1px solid var(--line,#ddd);border-radius:999px;padding:4px 12px;font-size:13px;background:var(--surface,#fff);color:inherit;cursor:pointer;font:inherit}.kn-node.on{background:var(--accent,var(--violet,#5b3fd1));color:#fff;border-color:transparent}.kn-node.ok{border-color:var(--good,#1e7a3c)}.kn-node[disabled]{opacity:.45;cursor:not-allowed}.kn-card{border:1px solid var(--line,#ddd);border-radius:14px;padding:18px;background:var(--surface,#fff);display:flex;flex-direction:column;gap:12px}.kn-card h2,.kn-card h3{margin:0}.kn-bar{height:8px;border-radius:5px;background:var(--line,#e5e5e5);overflow:hidden}.kn-bar i{display:block;height:100%;background:var(--accent,var(--violet,#5b3fd1));transition:width .3s}.kn-opts{display:flex;flex-direction:column;gap:8px}.kn-opt{text-align:start;border:1px solid var(--line,#ddd);border-radius:10px;padding:10px 12px;background:transparent;color:inherit;font:inherit;cursor:pointer}.kn-opt:hover,.kn-opt:focus-visible{border-color:var(--accent,var(--violet,#5b3fd1));outline:none}.kn-grid{display:grid;grid-template-columns:repeat(9,1fr);gap:3px;max-width:360px;aspect-ratio:1;direction:ltr}.kn-grid button{border:0;border-radius:4px;cursor:pointer}.kn-big{font-size:46px;font-weight:700;text-align:center;min-height:70px;letter-spacing:6px;direction:ltr}.kn-row{display:flex;flex-wrap:wrap;gap:8px;align-items:center}.kn-btn{border:1px solid var(--accent,var(--violet,#5b3fd1));background:var(--accent,var(--violet,#5b3fd1));color:#fff;border-radius:10px;padding:8px 16px;font:inherit;cursor:pointer}.kn-btn.ghost{background:transparent;color:var(--accent,var(--violet,#5b3fd1))}.kn-tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px}.kn-tile{border:1px solid var(--line,#ddd);border-radius:12px;padding:12px;display:flex;flex-direction:column;gap:4px}.kn-tile b{font-size:22px}.kn-tile small{opacity:.75}.kn-lv{display:inline-block;border-radius:999px;padding:0 10px;font-size:12px}.kn-lv.g{background:#e1f2e6;color:#1e7a3c}.kn-lv.y{background:#fbf0d5;color:#8a5a00}.kn-lv.r{background:#fbe3e1;color:#b3261e}.kn-safe{border:2px solid #b3261e;border-radius:12px;padding:14px;background:rgba(179,38,30,.06)}.kn-notif{position:fixed;top:16px;inset-inline-end:16px;z-index:9999;background:#1f2937;color:#fff;border-radius:12px;padding:10px 14px;box-shadow:0 8px 30px rgba(0,0,0,.3);display:flex;gap:10px;align-items:center;max-width:320px}.kn-notif button{border:0;border-radius:8px;padding:4px 10px;cursor:pointer;font:inherit}.kn-sci h3{margin-top:14px}.kn-sci p,.kn-sci li{line-height:1.9}.kn-sci table{border-collapse:collapse;width:100%}.kn-sci td,.kn-sci th{border:1px solid var(--line,#ddd);padding:6px 8px;text-align:start;vertical-align:top}.kn-stroop{font-size:44px;font-weight:800;text-align:center;min-height:64px}.kn-scale{display:grid;grid-template-columns:repeat(auto-fit,minmax(110px,1fr));gap:6px}';
  function css() { if (document.getElementById('kn-css')) return; var s = document.createElement('style'); s.id = 'kn-css'; s.textContent = CSS; document.head.appendChild(s); }
  /* ---------- admin gates (Firestore gates/global and gates/{uid}) and sharing with the admin ---------- */
  var gate = { mode: 'seq', upto: 7, hold: false, released: false, note: '' };
  var fsOK = function () { try { return !!(window.ACCT && ACCT.fs && ACCT.fs.ready()); } catch (e) { return false; } };
  var IDX = function (id) { return STAGES.findIndex(function (s) { return s[0] === id; }); };
  function flagged() { var w = st.r.well; return !!w && (w.phq >= 10 || w.gad >= 10 || w.who * 4 <= 28); }
  async function loadGates() { if (!fsOK()) return false; try { var g = await ACCT.fs.get('gates/global'), me = await ACCT.fs.get('gates/' + ACCT.fs.uid()); var e = { mode: 'seq', upto: 7, hold: false, released: false, note: '' }; if (g) Object.assign(e, { mode: g.mode || 'seq', upto: g.upto != null ? g.upto : 7, hold: !!g.hold, note: g.note || '' }); if (me) { if (me.mode && me.mode !== 'inherit') { e.mode = me.mode; e.upto = me.upto != null ? me.upto : e.upto; } e.released = !!me.released; if (me.note) e.note = me.note; } gate = e; return true; } catch (x) { return false; } }
  async function sync() { if (!st.share || !fsOK()) return; var data = { name: (ACCT.fs.name() || '').slice(0, 60), email: ACCT.fs.email(), consent: true, results: JSON.stringify(st.r).slice(0, 59000), done: Object.keys(st.done).length }; try { await ACCT.fs.set('know/' + ACCT.fs.uid(), data); } catch (e) { setTimeout(function () { ACCT.fs.set('know/' + ACCT.fs.uid(), data).catch(function () { }); }, 2600); } }
  function lockReason(id) { var i = IDX(id); if (i <= 0) return ''; if (gate.mode === 'locked') return 'الرحلة مقفولة حاليًا من مسؤول الموقع.' + (gate.note ? ' ' + E(gate.note) : ''); if (gate.hold && flagged() && !gate.released && i > IDX('well')) return 'نتايجك في مرحلة «الحالة النفسية» محتاجة مراجعة من مسؤول الموقع قبل ما تكمل، وهو هيفتح لك المرحلة الجاية. لو محتاج كلام مع حد دلوقتي: خط الصحة النفسية <b dir="ltr">16328</b> أو <b dir="ltr">08008880700</b>، والإسعاف <b dir="ltr">123</b>.'; if (gate.mode === 'open') return ''; if (gate.mode === 'upto' && i > gate.upto) return 'المسؤول فاتح الرحلة لحد مرحلة «' + STAGES[gate.upto][1] + '» دلوقتي. باقي المراحل هتتفتح قريب.' + (gate.note ? ' ' + E(gate.note) : ''); return (st.done[STAGES[i - 1][0]] || st.done[id]) ? '' : 'خلّص المرحلة اللي قبلها الأول.'; }
  function unlocked(id) { return !lockReason(id); }
  function map(cur) { return '<nav class="kn-map" aria-label="مراحل الرحلة">' + STAGES.map(function (s, i) { return '<button class="kn-node ' + (s[0] === cur ? 'on' : '') + (st.done[s[0]] ? ' ok' : '') + '" data-st="' + s[0] + '" ' + (unlocked(s[0]) ? '' : 'disabled') + '>' + (i + 1) + '. ' + s[1] + (st.done[s[0]] ? ' ✓' : '') + '</button>'; }).join('') + '<button class="kn-node" data-st="sci">العلم وراء M2 وM6</button></nav>'; }
  function frame(cur, body) { root.innerHTML = '<div class="kn"><div class="kn-row" style="justify-content:space-between"><h2 style="margin:0">اعرف نفسك</h2><span class="muted small">رحلة اختيارية، حوالي 15 دقيقة</span></div>' + map(cur) + body + '</div>'; root.querySelectorAll('[data-st]').forEach(function (b) { b.onclick = function () { go(b.dataset.st); }; }); }
  function finish(id, data, pts) { st.r[id] = data; st.done[id] = Date.now(); save(); xp(pts || 15); sync(); var i = STAGES.findIndex(function (s) { return s[0] === id; }); go(STAGES[Math.min(STAGES.length - 1, i + 1)][0]); }
  function go(id) { if (id !== 'sci' && !unlocked(id)) { var why = lockReason(id); if (!root.innerHTML || /خلّص/.test(why) === false) frame(id, '<section class="kn-card"><h3>المرحلة دي مقفولة</h3><p style="margin:0">' + why + '</p></section>'); else say('خلّص المرحلة اللي قبلها الأول'); return; } (V[id] || V.consent)(); try { window.scrollTo({ top: 0 }); } catch (e) { } }

  /* ---------- generic card-by-card questionnaire ---------- */
  function cards(cur, title, intro, items, onDone) {
    var i = 0, ans = [];
    function draw() {
      var it = items[i];
      frame(cur, '<section class="kn-card"><div class="kn-row" style="justify-content:space-between"><h3>' + E(title) + '</h3><span class="muted small">' + (i + 1) + ' / ' + items.length + '</span></div><div class="kn-bar"><i style="width:' + (i / items.length * 100) + '%"></i></div>' + (i === 0 && intro ? '<p class="muted" style="margin:0">' + intro + '</p>' : '') + '<p style="font-size:18px;margin:4px 0">' + E(it.text) + '</p><div class="' + (it.grid ? 'kn-scale' : 'kn-opts') + '">' + it.options.map(function (o, k) { return '<button class="kn-opt" data-k="' + k + '">' + E(o) + '</button>'; }).join('') + '</div>' + (i > 0 ? '<div><button class="kn-btn ghost" id="knBack">رجوع</button></div>' : '') + '</section>');
      root.querySelectorAll('.kn-opt').forEach(function (b) { b.onclick = function () { ans[i] = +b.dataset.k; i++; if (i >= items.length) onDone(ans); else draw(); }; });
      var bk = document.getElementById('knBack'); if (bk) bk.onclick = function () { i--; draw(); };
    }
    draw();
  }
  var LIK = ['تنطبق تمامًا', 'تنطبق', 'إلى حد ما', 'لا تنطبق', 'لا تنطبق إطلاقًا'];
  var FREQ = ['ولا مرة', 'عدة أيام', 'أكثر من نصف الأيام', 'تقريبًا كل يوم'];
  var WHO = ['طوال الوقت', 'معظم الوقت', 'أكثر من نصف الوقت', 'أقل من نصف الوقت', 'أحيانًا', 'ولا مرة'];
  // Working Arabic translations for teaching. For research, use the published validated Arabic versions.
  var WHO5 = ['شعرت بالبهجة وبمزاج جيد', 'شعرت بالهدوء والاسترخاء', 'شعرت بالنشاط والحيوية', 'استيقظت وأنا أشعر بالانتعاش والراحة', 'كانت حياتي اليومية مليئة بأشياء تهمني'];
  var PHQ8 = ['قلة الاهتمام أو المتعة في فعل الأشياء', 'الشعور بالإحباط أو الاكتئاب أو اليأس', 'صعوبة في النوم أو النوم المتقطع أو النوم أكثر من اللازم', 'الشعور بالتعب أو قلة الطاقة', 'ضعف الشهية أو الإفراط في الأكل', 'الشعور بعدم الرضا عن نفسك، أو أنك فاشل أو خذلت نفسك أو أسرتك', 'صعوبة التركيز في أشياء مثل القراءة أو مشاهدة التلفزيون', 'البطء في الحركة أو الكلام لدرجة لاحظها الآخرون، أو العكس: التململ وكثرة الحركة'];
  var GAD7 = ['الشعور بالتوتر أو القلق أو العصبية', 'عدم القدرة على إيقاف القلق أو التحكم فيه', 'القلق الزائد بشأن أمور مختلفة', 'صعوبة الاسترخاء', 'التململ لدرجة صعوبة الجلوس بهدوء', 'سرعة الانزعاج أو الضيق', 'الخوف كأن شيئًا فظيعًا قد يحدث'];

  var V = {};
  V.consent = function () {
    frame('consent', '<section class="kn-card"><h3>قبل ما نبدأ</h3><p style="margin:0">دي رحلة من 7 مراحل قصيرة بشكل لعبة: المزاج، والذاكرة، والتركيز، والحالة النفسية، وحياتك الرقمية، وعلاقتك بالذكاء الاصطناعي. في الآخر بيظهر لك بروفايل، ومعاه شرح علمي لكل نتيجة.</p><ul style="margin:0;line-height:1.9"><li><b>ده فرز تعليمي مش تشخيص.</b> النتيجة العالية معناها «يستحق الاهتمام والتحدث مع متخصص»، مش «عندك اضطراب».</li><li><b>إجاباتك محفوظة على جهازك بس.</b> مفيش حاجة بتتبعت لأي حد إلا لو إنت اخترت تصدّرها.</li><li><b>بعض المهام غرضها مش ظاهر من أول نظرة،</b> وده جزء من القياس العلمي. هنشرح لك كل حاجة في صفحة النتائج.</li><li>تقدر توقف في أي وقت وترجع تكمل.</li></ul><h3>سؤالين نتأكد بيهم إن الصورة واضحة</h3><fieldset class="kn-opts" style="border:0;padding:0;margin:0"><legend style="margin-bottom:6px">1) النتيجة اللي هتطلع:</legend><label class="kn-opt"><input type="radio" name="q1" value="0"> تشخيص نهائي لحالتي</label><label class="kn-opt"><input type="radio" name="q1" value="1"> فرز تعليمي، والقرار للمتخصص</label></fieldset><fieldset class="kn-opts" style="border:0;padding:0;margin:0"><legend style="margin-bottom:6px">2) مين يقدر يشوف إجاباتي؟</legend><label class="kn-opt"><input type="radio" name="q2" value="1"> أنا بس على جهازي، إلا لو صدّرتها بنفسي</label><label class="kn-opt"><input type="radio" name="q2" value="0"> أي حد بيفتح الموقع</label></fieldset><label class="kn-row"><input type="checkbox" id="knOk"> فهمت، وموافق أبدأ</label>' + (window.ACCT && ACCT.fs ? '<label class="kn-row" style="align-items:flex-start"><input type="checkbox" id="knShare" ' + (st.share ? 'checked' : '') + (fsOK() ? '' : ' disabled') + '> <span>اختياري: أوافق إن مسؤول الموقع يشوف نتايجي <b>باسمي وإيميلي</b> علشان يقدر يساعدني أو يوجهني لمتخصص. أقدر ألغي ده في أي وقت من صفحة النتائج.' + (fsOK() ? '' : ' <span class="muted">(سجّل دخولك من «حسابي» الأول علشان تقدر تفعّل ده)</span>') + '</span></label>' : '') + '<div id="knMsg" class="muted small"></div><div><button class="kn-btn" id="knGo">ابدأ الرحلة</button></div></section>');
    document.getElementById('knGo').onclick = function () { var a = root.querySelector('[name=q1]:checked'), b = root.querySelector('[name=q2]:checked'), m = document.getElementById('knMsg');
      if (!a || !b) { m.textContent = 'جاوب على السؤالين الأول.'; return; }
      if (a.value !== '1' || b.value !== '1') { m.textContent = 'راجع الإجابة: النتيجة فرز تعليمي مش تشخيص، والإجابات على جهازك بس.'; return; }
      if (!document.getElementById('knOk').checked) { m.textContent = 'علّم على «موافق» علشان نبدأ.'; return; }
      var sh = document.getElementById('knShare'); st.share = !!(sh && sh.checked); save(); finish('consent', { at: new Date().toISOString(), share: st.share }, 5); };
  };
  V.mood = function () {
    var cells = ''; for (var y = 0; y < 9; y++) for (var x = 0; x < 9; x++) { var v = x - 4, a = 4 - y; cells += '<button aria-label="المتعة ' + v + '، الطاقة ' + a + '" data-v="' + v + '" data-a="' + a + '" style="background:hsl(' + (v >= 0 ? 150 : 10) + ',' + (30 + Math.abs(v) * 12) + '%,' + (88 - Math.abs(a) * 7) + '%)"></button>'; }
    frame('mood', '<section class="kn-card"><h3>المرحلة 1: مزاجك دلوقتي</h3><p class="muted" style="margin:0">دي «شبكة المشاعر» (Russell, Weiss & Mendelsohn, 1989). اليمين واليسار = مريح أو مزعج، وفوق وتحت = طاقة عالية أو هادي. اضغط على المربع الأقرب لإحساسك الآن.</p><div class="kn-row" style="gap:16px;align-items:flex-start"><div><div class="muted small" style="text-align:center">طاقة عالية</div><div class="kn-grid">' + cells + '</div><div class="muted small" style="text-align:center">هادئ / طاقة منخفضة</div><div class="kn-row muted small" style="justify-content:space-between;direction:ltr;max-width:360px"><span>مزعج</span><span>مريح</span></div></div></div></section>');
    root.querySelectorAll('.kn-grid button').forEach(function (b) { b.onclick = function () { finish('mood', { v: +b.dataset.v, a: +b.dataset.a }, 10); }; });
  };
  V.span = function () {
    var len = 3, fails = 0, best = 0, trial = 0, probe = { shown: false, opened: false, lat: null, dismissed: false }, probeAt = 1;
    function notif() { probe.shown = true; var t0 = performance.now(), n = document.createElement('div'); n.className = 'kn-notif'; n.setAttribute('role', 'status'); n.innerHTML = '<span>💬 <b>مجموعة الدفعة</b><br><small>رسالة جديدة</small></span><button id="knNo" style="background:#374151;color:#fff">لاحقًا</button><button id="knOpen" style="background:#fff">افتح</button>'; document.body.appendChild(n);
      var end = function (opened) { if (!n.parentNode) return; probe.lat = Math.round(performance.now() - t0); probe.opened = opened; probe.dismissed = !opened; n.remove(); };
      n.querySelector('#knOpen').onclick = function () { end(true); say('مفيش رسالة فعلًا 🙂 هنشرح ليه في النتائج'); }; n.querySelector('#knNo').onclick = function () { end(false); }; setTimeout(function () { if (n.parentNode) { probe.lat = null; n.remove(); } }, 7000); }
    function run() {
      var d = Array.from({ length: len }, function () { return Math.floor(Math.random() * 10); }), k = 0;
      frame('span', '<section class="kn-card"><h3>المرحلة 2: الذاكرة (مدى الأرقام)</h3><p class="muted" style="margin:0">ركّز: هتظهر أرقام واحد ورا التاني. بعدها اكتبها بنفس الترتيب. الطول بيزيد كل ما تنجح.</p><div class="kn-big" id="knD">…</div><div class="muted small">الطول الحالي: ' + len + ' أرقام · أفضل نتيجة: ' + best + '</div></section>');
      if (trial === probeAt && !probe.shown) setTimeout(notif, 900);
      var el = document.getElementById('knD'), tick = setInterval(function () { if (k < d.length) { el.textContent = d[k]; k++; setTimeout(function () { if (el.isConnected) el.textContent = ''; }, 650); } else { clearInterval(tick); ask(d); } }, 1000);
    }
    function ask(d) {
      var c = root.querySelector('.kn-card'); c.insertAdjacentHTML('beforeend', '<label for="knIn">اكتب الأرقام:</label><input id="knIn" inputmode="numeric" autocomplete="off" dir="ltr" style="font-size:22px;letter-spacing:4px;padding:8px;border-radius:8px;border:1px solid var(--line,#ccc)"><div><button class="kn-btn" id="knChk">تأكيد</button></div>');
      var inp = document.getElementById('knIn'); inp.focus();
      var chk = function () { var ok = inp.value.replace(/\D/g, '') === d.join(''); trial++; if (ok) { best = Math.max(best, len); len++; fails = 0; } else { fails++; } if (fails >= 2 || len > 10) return finish('span', { best: best, trials: trial, probe: probe }, 20); say(ok ? 'صح! الطول هيزيد' : 'قريب، جرّب تاني'); run(); };
      document.getElementById('knChk').onclick = chk; inp.onkeydown = function (e) { if (e.key === 'Enter') chk(); };
    }
    run();
  };
  V.stroop = function () {
    var C = [['أحمر', '#d62828'], ['أزرق', '#1d4ed8'], ['أخضر', '#15803d'], ['أصفر', '#ca8a04']], N = 24, t = [], i = 0, t0 = 0;
    var seq = Array.from({ length: N }, function (_, k) { var w = Math.floor(Math.random() * 4), ink = k % 2 === 0 ? w : (w + 1 + Math.floor(Math.random() * 3)) % 4; return { w: w, ink: ink }; }).sort(function () { return Math.random() - 0.5; });
    frame('stroop', '<section class="kn-card"><h3>المرحلة 3: التركيز (ستروب)</h3><p class="muted" style="margin:0">اضغط على <b>لون الحبر</b> مش الكلمة نفسها. مثلًا لو كلمة «أحمر» مكتوبة بالأزرق، اضغط «أزرق». بأسرع ما تقدر من غير أخطاء.</p><div class="kn-stroop" id="knS">جاهز؟</div><div class="kn-row" id="knB" style="justify-content:center">' + C.map(function (c, k) { return '<button class="kn-btn" data-c="' + k + '" style="background:' + c[1] + ';border-color:' + c[1] + '">' + c[0] + '</button>'; }).join('') + '</div><div class="kn-bar"><i id="knP" style="width:0"></i></div></section>');
    var s = document.getElementById('knS');
    function next() { if (i >= N) { var ok = t.filter(function (x) { return x.ok; }), m = function (a) { return a.length ? a.reduce(function (p, x) { return p + x.rt; }, 0) / a.length : NaN; }; var cg = m(ok.filter(function (x) { return x.cong; })), ic = m(ok.filter(function (x) { return !x.cong; })); return finish('stroop', { cong: Math.round(cg), incong: Math.round(ic), interf: Math.round(ic - cg), err: t.filter(function (x) { return !x.ok; }).length }, 20); }
      s.textContent = '+'; s.style.color = 'inherit'; setTimeout(function () { var q = seq[i]; s.textContent = C[q.w][0]; s.style.color = C[q.ink][1]; t0 = performance.now(); }, 450); }
    root.querySelectorAll('#knB button').forEach(function (b) { b.onclick = function () { if (!t0 || i >= N) return; var q = seq[i]; t.push({ rt: performance.now() - t0, ok: +b.dataset.c === q.ink, cong: q.w === q.ink }); t0 = 0; i++; document.getElementById('knP').style.width = (i / N * 100) + '%'; next(); }; });
    setTimeout(next, 900);
  };
  V.well = function () {
    var items = WHO5.map(function (t) { return { text: 'خلال الأسبوعين الماضيين: ' + t, options: WHO, grid: true, s: 'who' }; })
      .concat(PHQ8.map(function (t) { return { text: 'خلال الأسبوعين الماضيين، كم مرة ضايقك: ' + t + '؟', options: FREQ, grid: true, s: 'phq' }; }))
      .concat(GAD7.map(function (t) { return { text: 'خلال الأسبوعين الماضيين، كم مرة ضايقك: ' + t + '؟', options: FREQ, grid: true, s: 'gad' }; }));
    cards('well', 'المرحلة 4: الحالة النفسية', 'ثلاث أدوات فرز عالمية مجانية: WHO-5 للرفاهية النفسية، وPHQ-8 للمزاج، وGAD-7 للقلق. الصياغة هنا ترجمة عمل للتدريب؛ للبحث استخدم النسخ العربية المقننة.', items, function (a) {
      var who = 0, phq = 0, gad = 0; items.forEach(function (it, k) { if (it.s === 'who') who += 5 - a[k]; else if (it.s === 'phq') phq += a[k]; else gad += a[k]; });
      finish('well', { who: who, phq: phq, gad: gad }, 25); });
  };
  function libStage(id, cur, title, intro, extra) {
    var c = CLIB.get(id), dir = CLIB.directItems(c), ind = CLIB.indirectItems(c), items = [];
    var mix = dir.map(function (x) { return Object.assign({}, x, { options: LIK, grid: true }); }); ind.forEach(function (x, k) { mix.splice(Math.min(mix.length, 3 + k * 3), 0, x); }); items = mix.concat(extra || []);
    cards(cur, title, intro, items, function (a) {
      var main = items.filter(function (x) { return x.kind !== 'cal'; }), am = a.slice(0, main.length), sc = CLIB.score(main.map(function (x) { return x.grid && x.kind === 'd' ? Object.assign({}, x, { w: null }) : x; }), am, 5);
      var fac = {}; main.forEach(function (x, k) { if (x.kind === 'd') { var v = x.reverse ? am[k] : 4 - am[k]; (fac[x.dim] = fac[x.dim] || []).push(v / 4 * 100); } });
      var facets = Object.keys(fac).map(function (k) { return [k, Math.round(fac[k].reduce(function (p, v) { return p + v; }, 0) / fac[k].length)]; });
      var out = { direct: sc.direct, indirect: sc.indirect, combined: sc.combined, gap: sc.gap, facets: facets };
      if (extra && extra.length) { var cal = items.map(function (x, k) { return x.kind === 'cal' ? [x, a[k]] : null; }).filter(Boolean); var wrongAcc = cal.filter(function (p) { return !p[0].ok && p[1] === 0; }).length, wrongN = cal.filter(function (p) { return !p[0].ok; }).length, rightRej = cal.filter(function (p) { return p[0].ok && p[1] === 2; }).length, verify = cal.filter(function (p) { return p[1] === 1; }).length; out.cal = { wrongAcc: wrongAcc, wrongN: wrongN, rightRej: rightRej, verify: verify, n: cal.length }; }
      finish(cur, out, 30); });
  }
  V.m6 = function () { libStage('m6_avail', 'm6', 'المرحلة 5: حياتك الرقمية', 'أسئلة ومواقف عن الرسائل وتطبيقات التواصل في يومك. اختار الأقرب ليك، ومفيش إجابة صح أو غلط.'); };
  V.m2 = function () {
    var claims = [['معامل ألفا كرونباخ يقيس صدق المقياس.', false], ['الوسيط أقل تأثرًا بالقيم المتطرفة من المتوسط.', true], ['زيادة حجم العينة لا تؤثر على الدلالة الإحصائية.', false], ['التوزيع الطبيعي متماثل حول متوسطه.', true], ['ارتباط قيمته 0.9 يثبت أن المتغير الأول سبب في الثاني.', false], ['سعة الذاكرة العاملة محدودة.', true]];
    var extra = claims.map(function (c) { return { text: 'سألتَ مساعدًا ذكيًا سؤالًا في المقرر، فرد بثقة: «' + c[0] + '». ماذا تفعل؟', options: ['أقبلها وأستخدمها كما هي', 'أتحقق منها في مرجع قبل استخدامها', 'أرفضها لأنها غلط'], kind: 'cal', ok: c[1] }; });
    libStage('m2_aidep', 'm2', 'المرحلة 6: أنت والذكاء الاصطناعي', 'أسئلة عن طريقة استخدامك لأدوات الذكاء الاصطناعي في الدراسة والحياة، وفي الآخر 6 ردود من «مساعد ذكي» تقرر تتعامل معاها إزاي.', extra);
  };
  function lv(x, cut) { return x >= cut[1] ? ['r', 'مرتفع'] : x >= cut[0] ? ['y', 'متوسط'] : ['g', 'منخفض']; }
  function pct(x) { return x == null || !isFinite(x) ? '—' : Math.round(x) + '%'; }
  V.res = function () {
    var r = st.r, w = r.well || {}, sp = r.span || {}, stp = r.stroop || {}, m6 = r.m6, m2 = r.m2, md = r.mood;
    var phqL = w.phq == null ? null : w.phq >= 15 ? ['r', 'مرتفع'] : w.phq >= 10 ? ['y', 'متوسط'] : w.phq >= 5 ? ['y', 'خفيف'] : ['g', 'منخفض'];
    var gadL = w.gad == null ? null : w.gad >= 15 ? ['r', 'مرتفع'] : w.gad >= 10 ? ['y', 'متوسط'] : w.gad >= 5 ? ['y', 'خفيف'] : ['g', 'منخفض'];
    var whoP = w.who == null ? null : w.who * 4, flag = (w.phq >= 10) || (w.gad >= 10) || (whoP != null && whoP <= 28);
    var moodTxt = md ? (md.v >= 0 ? (md.a >= 0 ? 'مريح ونشيط' : 'مريح وهادئ') : (md.a >= 0 ? 'متوتر أو منزعج' : 'مُحبط أو مرهق')) : '—';
    var tile = function (t, v, s, l) { return '<div class="kn-tile"><small>' + t + '</small><b>' + v + '</b>' + (l ? '<span class="kn-lv ' + l[0] + '">' + l[1] + '</span>' : '') + '<small>' + (s || '') + '</small></div>'; };
    var gapTxt = function (g) { return g == null ? '' : Math.abs(g) < 12 ? 'المباشر وغير المباشر متقاربين' : g > 0 ? 'وصفك المباشر أعلى من مؤشراتك غير المباشرة' : 'مؤشراتك غير المباشرة أعلى من وصفك المباشر'; };
    var fac = function (m) { return m && m.facets ? '<ul style="margin:4px 0">' + m.facets.map(function (f) { return '<li>' + E(f[0]) + ': ' + f[1] + '%</li>'; }).join('') + '</ul>' : ''; };
    frame('res', (flag ? '<section class="kn-safe"><h3 style="margin-top:0">نتايجك بتقول إنك ممكن تكون مرهق اليومين دول</h3><p style="margin:0">ده مش تشخيص، بس يستاهل تتكلم مع متخصص نفسي أو مع حد تثق فيه. في مصر: خط الصحة النفسية <b dir="ltr">16328</b> أو <b dir="ltr">08008880700</b> (اتأكد من الرقم الحالي من وزارة الصحة)، وفي الطوارئ الإسعاف <b dir="ltr">123</b>. خارج مصر: <a href="https://findahelpline.com" target="_blank" rel="noopener">findahelpline.com</a> فيه خطوط الدعم في كل دولة.</p></section>' : '') +
      '<section class="kn-card"><h3>بروفايلك</h3><div class="kn-tiles">' + tile('المزاج الآن', E(moodTxt), 'من شبكة المشاعر') + tile('مدى الأرقام', sp.best != null ? sp.best : '—', 'المتوسط عند الراشدين حوالي 7 ± 2 (Miller, 1956)') + tile('تداخل ستروب', stp.interf != null && isFinite(stp.interf) ? stp.interf + ' مللي ثانية' : '—', 'الفرق بين الكلمات المتعارضة والمتطابقة؛ أقل = تحكم انتباهي أعلى') + tile('الرفاهية (WHO-5)', whoP != null ? whoP + '%' : '—', '50% أو أقل تستاهل متابعة', whoP != null ? (whoP <= 28 ? ['r', 'منخفضة جدًا'] : whoP <= 50 ? ['y', 'منخفضة'] : ['g', 'جيدة']) : null) + tile('المزاج (PHQ-8)', w.phq != null ? w.phq + ' / 24' : '—', '10 فأكثر = يستاهل تقييم متخصص', phqL) + tile('القلق (GAD-7)', w.gad != null ? w.gad + ' / 21' : '—', '10 فأكثر = يستاهل تقييم متخصص', gadL) + '</div></section>' +
      ['m6', 'm2'].map(function (k) { var m = r[k], nm = k === 'm6' ? 'ضغط التواجد الرقمي الدائم (M6)' : 'الاعتماد المعرفي على الذكاء الاصطناعي (M2)'; if (!m) return ''; var probe = k === 'm6' && sp.probe && sp.probe.shown ? '<p style="margin:0"><b>مؤشر سلوكي مخفي:</b> أثناء مرحلة الذاكرة ظهر لك إشعار رسالة وهمي. ' + (sp.probe.opened ? 'فتحته بعد ' + (sp.probe.lat / 1000).toFixed(1) + ' ثانية، وده مؤشر لـ«إلحاح الرد» حتى وإنت مشغول بمهمة.' : sp.probe.dismissed ? 'أجّلته، وده مؤشر لقدرتك على تأجيل الرد وإنت مركز.' : 'تجاهلته لحد ما اختفى.') + '</p>' : ''; var cal = k === 'm2' && m.cal ? '<p style="margin:0"><b>مهمة معايرة الثقة:</b> من ' + m.cal.wrongN + ' ردود غلط من «المساعد الذكي» قبلت ' + m.cal.wrongAcc + ' كما هي، واخترت التحقق في ' + m.cal.verify + ' من ' + m.cal.n + ' ردود. قبول الرد الغلط بثقة هو «تحيز الأتمتة»، والتحقق هو السلوك الأمثل حتى لو الرد صح.</p>' : '';
        return '<section class="kn-card"><h3>' + nm + '</h3><div class="kn-tiles">' + tile('الدرجة المباشرة', pct(m.direct), 'من العبارات الصريحة', m.direct != null ? lv(m.direct, [45, 65]) : null) + tile('الدرجة غير المباشرة', pct(m.indirect), 'من المواقف والاختيارات', m.indirect != null ? lv(m.indirect, [45, 65]) : null) + tile('الدرجة المجمّعة', pct(m.combined), gapTxt(m.gap)) + '</div><div><b>الجوانب:</b>' + fac(m) + '</div>' + probe + cal + '<p class="muted small" style="margin:0">مقياس أصيل قيد التطوير؛ لسه مفيش معايير، فالنسب دي من الحد الأقصى مش مقارنة بغيرك.</p></section>'; }).join('') +
      '<section class="kn-card"><h3>إيه اللي كان مخفي في الرحلة؟ (Debriefing)</h3><ul style="margin:0;line-height:1.9"><li><b>إشعار الرسالة الوهمي</b> في مرحلة الذاكرة كان مقياس سلوكي لإلحاح الرد (مكوّن في M6). مفيش رسالة حقيقية.</li><li><b>ردود «المساعد الذكي»</b> في المرحلة 6 كان فيها 3 جمل غلط عمدًا، علشان نقيس معايرة الثقة (مكوّن في M2).</li><li><b>أسئلة المواقف والتقدير عن الآخرين</b> مؤشرات غير مباشرة: الناس بتسقط سلوكها على غيرها، وبتفسّر المواقف الغامضة حسب حالتها.</li></ul><p class="muted small" style="margin:0">الإخفاء ده مسموح بيه في البحث بشرط الشرح بعده، وإنك تقدر تمسح بياناتك. لو مش موافق على استخدام الجزء ده، امسح نتايجك من تحت.</p></section>' +
      '<section class="kn-card"><h3>اتعلم من إجاباتك</h3><ul style="margin:0;line-height:1.9"><li>ليه نتيجة PHQ-8 العالية <b>مش تشخيص</b>؟ لأن القيمة التنبؤية للفرز بتعتمد على انتشار الحالة في المجتمع (نظرية بايز). جرّب محاكي معدل الانتشار في الدروس.</li><li>ليه بنحسب درجة مباشرة وغير مباشرة؟ ده تطبيق مصغّر لمصفوفة «السمات والطرق المتعددة» (Campbell & Fiske, 1959): الارتباط بينهم دليل صدق تقاربي.</li><li>تأثير ستروب (Stroop, 1935) من أشهر المهام في قياس الوظائف التنفيذية، وثباته بين الأفراد أقل مما يتخيل الناس (Hedge et al., 2018).</li></ul></section>' +
      (window.ACCT && ACCT.fs ? '<section class="kn-card"><h3>المشاركة مع مسؤول الموقع</h3><p style="margin:0">' + (st.share ? 'نتايجك متشاركة باسمك مع مسؤول الموقع علشان يقدر يساعدك.' : 'نتايجك على جهازك بس.') + '</p><div><button class="kn-btn ghost" id="knSh">' + (st.share ? 'الغِ المشاركة وامسح نسختي عند المسؤول' : 'شارك نتايجي باسمي مع المسؤول') + '</button></div></section>' : '') + '<section class="kn-card"><div class="kn-row"><button class="kn-btn" id="knExp">انسخ نتايجي (JSON)</button><button class="kn-btn ghost" id="knSci">اقرأ العلم وراء M2 وM6</button><button class="kn-btn ghost" id="knReset">امسح نتايجي وابدأ من الأول</button></div></section>');
    document.getElementById('knExp').onclick = function () { var t = JSON.stringify(st, null, 1); (navigator.clipboard ? navigator.clipboard.writeText(t) : Promise.reject()).then(function () { say('اتنسخت'); }, function () { say('النسخ اتمنع في المتصفح ده'); }); };
    document.getElementById('knSci').onclick = function () { go('sci'); };
    var shb = document.getElementById('knSh'); if (shb) shb.onclick = async function () { if (!fsOK()) return say('سجّل دخولك من «حسابي» الأول'); if (st.share) { st.share = false; save(); try { await ACCT.fs.del('know/' + ACCT.fs.uid()); } catch (e) { } say('اتلغت المشاركة'); } else { st.share = true; save(); await sync(); say('اتشاركت نتايجك مع المسؤول'); } V.res(); };
    document.getElementById('knReset').onclick = function () { st = { v: 1, done: {}, r: {} }; save(); go('consent'); };
    if (!st.done.res) { st.done.res = Date.now(); save(); xp(40); }
  };
  V.sci = function () { frame('sci', '<article class="kn-card kn-sci">' + SCI + '</article>'); };

  var SCI = '<h2>مقياسان أصيلان قيد التطوير: الأساس العلمي</h2><p class="muted">إعداد: محمد عثمان، الدبلومة الخاصة في علم النفس، بمساعدة Claude. مسودة بحثية للمناقشة مع المشرفة قبل التقنين.</p>' +
  '<h3>أولًا: الاعتماد المعرفي على الذكاء الاصطناعي (M2)</h3>' +
  '<p><b>المشكلة.</b> انتقلت أدوات الذكاء الاصطناعي التوليدي في أقل من ثلاث سنوات من أداة بحث إلى شريك في التفكير والكتابة واتخاذ القرار. السؤال العلمي هنا ليس «هل يستخدم الطلاب هذه الأدوات؟» لأن الإجابة نعم، بل «متى يتحول التفويض المعرفي من امتداد تكيفي للعقل إلى اعتماد يُضعف القدرة المستقلة على التفكير والحكم؟». والقياس الجيد لهذا التحول شرط لأي دراسة سببية لاحقة عن أثره في التعلم.</p>' +
  '<p><b>التعريف الإجرائي.</b> نمط من تفويض العمليات المعرفية (التفكير، والتخطيط، والكتابة، والحكم) إلى أدوات الذكاء الاصطناعي يتجاوز ما تتطلبه المهمة، ويصاحبه ضعف التحقق من المخرجات، وتراجع الثقة في القدرة الذاتية بدون الأداة، وقد يمتد إلى اعتماد انفعالي وفقدان للتحكم.</p>' +
  '<p><b>الإطار النظري.</b> نبدأ من نظرية التفريغ المعرفي (Risko & Gilbert, 2016) ونظرية العقل الممتد (Clark & Chalmers, 1998)، وكلتاهما ترى أن الاستعانة بأدوات خارجية سلوك طبيعي وغالبًا مفيد. لذلك <b>الاستخدام نفسه ليس هو المفهوم</b>. المفهوم يبدأ عندما يجتمع التفويض مع خللين وثّقتهما أدبيات العوامل البشرية: <b>تحيز الأتمتة</b>، أي قبول مخرجات النظام الآلي دون تحقق (Parasuraman & Manzey, 2010)، و<b>سوء معايرة الثقة</b>، أي عدم تطابق الثقة مع الدقة الفعلية للنظام (Lee & See, 2004). وتؤكد مقاييس حديثة جدًا، مثل AIDep-22 لدى طلاب جامعيين صينيين ومقياس Goh وHartanto (2025)، أن الاعتماد متعدد الأبعاد، وأنه يرتبط سلبيًا بالتفكير النقدي والأداء. لكن هذه العلاقة ارتباطية حتى الآن، ولم يثبت أنها سببية.</p>' +
  '<p><b>ما الجديد في M2؟</b> (1) الفصل الصريح بين التفويض التكيفي والضار، عن طريق عبارات عكسية تصف التفويض الذكي (أكتب أفكاري ثم أستعين بالأداة). (2) <b>مهمة سلوكية لمعايرة الثقة</b> بدل الاكتفاء بالتقرير الذاتي: ردود من «مساعد ذكي» بعضها خاطئ عمدًا، ويُقاس قبول الخطأ والتحقق. (3) مسار غير مباشر بمواقف غامضة وتقدير عن الآخرين، يقلل أثر المرغوبية الاجتماعية في موضوع حساس أكاديميًا (الغش). (4) بناؤه بالعربية من الأساس للسياق العربي، وليس ترجمة.</p>' +
  '<table><thead><tr><th>الجانب</th><th>ما يقيسه</th><th>المبرر النظري</th></tr></thead><tbody><tr><td>التفويض المعرفي الزائد</td><td>تفويض قبل المحاولة الذاتية</td><td>التفريغ يصبح ضارًا عندما يحل محل الجهد المعرفي اللازم للتعلم (أثر «الصعوبة المرغوبة»، Bjork)</td></tr><tr><td>ضعف التحقق ومعايرة الثقة</td><td>قبول المخرجات دون مراجعة، والثقة في الأسلوب الواثق</td><td>تحيز الأتمتة ومعايرة الثقة</td></tr><tr><td>تآكل الكفاءة الذاتية</td><td>الشك في القدرة بدون الأداة</td><td>مصادر الكفاءة الذاتية عند Bandura: غياب خبرات النجاح المستقلة</td></tr><tr><td>الاعتماد الانفعالي</td><td>الراحة والبوح للمساعد الذكي</td><td>العلاقات شبه الاجتماعية مع الوكلاء المحادثين</td></tr><tr><td>فقدان التحكم</td><td>استخدام يتجاوز النية، وضيق عند الانقطاع</td><td>نموذج I-PACE (Brand et al., 2016) مع الحذر من نقل معايير الإدمان آليًا</td></tr></tbody></table>' +
  '<h3>ثانيًا: ضغط التواجد الرقمي الدائم (M6)</h3>' +
  '<p><b>المشكلة.</b> صممت تطبيقات التواصل مؤشرات اجتماعية جديدة لم تكن موجودة قبل 2014، مثل علامة «تمت القراءة» و«آخر ظهور» و«متصل الآن». هذه المؤشرات تحوّل التأخر في الرد من سلوك محايد إلى رسالة اجتماعية قابلة للتفسير. في الثقافات الجمعية، حيث الحفاظ على العلاقة والسمعة قيمة مركزية، يُتوقع أن يكون عبء هذه المؤشرات أعلى. ولا يوجد على حد علمنا مقياس يقيس هذا العبء النفسي بأبعاده الثقافية.</p>' +
  '<p><b>التعريف الإجرائي.</b> شعور داخلي بالإلزام بأن يكون الفرد متاحًا ومستجيبًا عبر تطبيقات التواصل في كل وقت. يظهر في إلحاح الرد الفوري، والمراقبة المستمرة، والانشغال الذهني بالمجال الرقمي، والقلق من مؤشرات التواجد، وتوقع أن الآخرين ينتظرون الرد.</p>' +
  '<p><b>الإطار النظري.</b> يمتد المفهوم من «اليقظة الرقمية» عند Reinecke وزملائه (2018)، وهي مقياس مقنن بثلاثة أبعاد: البروز، والاستجابية، والمراقبة. ويستفيد من «الضغط عن بعد» في بيئة العمل (Barber & Santuzzi, 2015). ويضيف بعدين لم يقيسهما أي منهما: <b>قلق مؤشرات التواجد</b> (تفسير الرؤية بلا رد على أنها تجاهل)، و<b>التوقعات المدركة للإتاحة</b> من الأسرة والأصدقاء.</p>' +
  '<p><b>التمييز عن المفاهيم القريبة (الصدق التمييزي المتوقع).</b> M6 ليس عدد ساعات الاستخدام: قد يستخدم شخص الهاتف كثيرًا للترفيه دون أي ضغط للتواجد. وليس «الخوف من فوات الأشياء» (Przybylski et al., 2013)، لأن مركز FOMO هو ما يفعله الآخرون، أما مركز M6 فهو ما ينتظره الآخرون مني. وليس الاستخدام الإشكالي للهاتف، لأن M6 يمكن أن يكون مرتفعًا لدى شخص يتحكم في استخدامه لكنه يتألم من الالتزام. هذه التوقعات تُسجل قبل جمع البيانات وتُختبر.</p>' +
  '<p><b>القياس متعدد الطرق.</b> (1) عبارات مباشرة سلوكية لا تذكر اسم المفهوم. (2) مواقف غامضة (صديق رأى رسالتك ولم يرد) تقيس تحيز التفسير. (3) تقدير عن الآخرين يستفيد من تأثير الإجماع الزائف (Ross, Greene & House, 1977). (4) <b>مؤشر سلوكي مخفي</b>: إشعار رسالة وهمي يظهر أثناء مهمة ذاكرة، ويُسجل فتحه وزمن الاستجابة. هذا المؤشر يقيس الاستجابية سلوكيًا وليس بالتقرير الذاتي، ويُشرح للمشارك بعد التطبيق.</p>' +
  '<h3>ثالثًا: الفرضيات القابلة للاختبار</h3><ol><li>تتشبع عبارات كل مقياس على عوامل متمايزة تطابق الجوانب المقترحة، مع عامل عام (نموذج ثنائي العامل).</li><li>ترتبط الدرجة المباشرة بغير المباشرة ارتباطًا موجبًا متوسطًا (0.30 إلى 0.50): صدق تقاربي عبر الطرق.</li><li>يرتبط M2 سلبيًا بأداء اختبار تفكير نقدي، وبدقة مهمة معايرة الثقة، بعد ضبط تكرار الاستخدام.</li><li>يرتبط M6 بجودة نوم أسوأ وتشتت أعلى في مهام الانتباه، بعد ضبط عدد ساعات الاستخدام وFOMO.</li><li>يكون M6 أعلى في الثقافات الأكثر جمعية، مع تكافؤ قياس جزئي على الأقل بين الدول (طريقة Alignment).</li><li>تختلف الأجيال في متوسطي المقياسين، ويُفصل أثر الجيل عن السن بقياسات متكررة لدفعات متتالية.</li></ol>' +
  '<h3>رابعًا: خطة التقنين</h3><table><thead><tr><th>الدراسة</th><th>الهدف</th><th>العينة والأدوات</th></tr></thead><tbody><tr><td>1. كيفية</td><td>اكتشاف الجوانب من الواقع</td><td>مقابلات ومجموعات نقاش (20 إلى 30 مشاركًا) وتحليل موضوعي</td></tr><tr><td>2. صدق المحتوى</td><td>تحكيم العبارات</td><td>8 إلى 10 خبراء، وCVR وCVI، ومقابلات معرفية مع 10 من الفئة</td></tr><tr><td>3. استكشافية</td><td>البنية العاملية الأولية</td><td>300 مشارك أو أكثر، وتحليل عاملي استكشافي وEGA، وكشف التكرار بـUVA</td></tr><tr><td>4. توكيدية</td><td>تأكيد البنية والتكافؤ</td><td>عينة مستقلة 400 أو أكثر، وCFA، وتكافؤ بين الذكور والإناث والأجيال والدول</td></tr><tr><td>5. الصدق</td><td>الشبكة المرتبطة والمحكات</td><td>مقاييس FOMO واليقظة الرقمية والتفكير النقدي وجودة النوم، والمعدل الدراسي</td></tr><tr><td>6. الثبات والتغير</td><td>الاستقرار والحالة اليومية</td><td>إعادة التطبيق بعد أسبوعين، ونسخة يومية لأسبوعين (EMA)</td></tr></tbody></table>' +
  '<h3>خامسًا: حدود يجب الاعتراف بها</h3><ul><li>العبارات الحالية مسودة أولى، ولا توجد معايير بعد. الدرجات في الرحلة نسب من الحد الأقصى وليست مقارنة بمجتمع.</li><li>المؤشرات غير المباشرة والسلوكية ثباتها عادة أقل من التقرير الذاتي، ولذلك تُستخدم مكمّلة وليست بديلة.</li><li>مهمة معايرة الثقة قصيرة (6 ردود) وتحتاج بنكًا أكبر ومتوازنًا في الصعوبة.</li><li>مجال الذكاء الاصطناعي يتغير بسرعة، فيلزم مراقبة انجراف خصائص العبارات من دفعة لأخرى.</li></ul>' +
  '<h3>المراجع الأساسية</h3><ul style="direction:ltr;text-align:left"><li>Barber, L. K., & Santuzzi, A. M. (2015). Please respond ASAP: Workplace telepressure and employee recovery. <i>Journal of Occupational Health Psychology, 20</i>(2), 172–189.</li><li>Brand, M., Young, K. S., Laier, C., Wölfling, K., & Potenza, M. N. (2016). Integrating psychological and neurobiological considerations regarding Internet-use disorders: The I-PACE model. <i>Neuroscience & Biobehavioral Reviews, 71</i>, 252–266.</li><li>Campbell, D. T., & Fiske, D. W. (1959). Convergent and discriminant validation by the multitrait-multimethod matrix. <i>Psychological Bulletin, 56</i>(2), 81–105.</li><li>Clark, A., & Chalmers, D. (1998). The extended mind. <i>Analysis, 58</i>(1), 7–19.</li><li>Clark, L. A., & Watson, D. (2019). Constructing validity: New developments in creating objective measuring instruments. <i>Psychological Assessment, 31</i>(12), 1412–1427.</li><li>Goh, A. Y. H., & Hartanto, A. (2025). Generative AI dependency scale. <i>Computers in Human Behavior Reports</i> (verify the full reference before citing).</li><li>Lee, J. D., & See, K. A. (2004). Trust in automation: Designing for appropriate reliance. <i>Human Factors, 46</i>(1), 50–80.</li><li>Parasuraman, R., & Manzey, D. H. (2010). Complacency and bias in human use of automation. <i>Human Factors, 52</i>(3), 381–410.</li><li>Przybylski, A. K., Murayama, K., DeHaan, C. R., & Gladwell, V. (2013). Motivational, emotional, and behavioral correlates of fear of missing out. <i>Computers in Human Behavior, 29</i>(4), 1841–1848.</li><li>Reinecke, L., et al. (2018). Permanently online and permanently connected: Development and validation of the Online Vigilance Scale. <i>PLOS ONE, 13</i>(10), e0205384.</li><li>Risko, E. F., & Gilbert, S. J. (2016). Cognitive offloading. <i>Trends in Cognitive Sciences, 20</i>(9), 676–688.</li><li>Ross, L., Greene, D., & House, P. (1977). The false consensus effect. <i>Journal of Experimental Social Psychology, 13</i>(3), 279–301.</li><li>Wu, et al. (2025/2026). Development and validation of the AI dependence scale for Chinese undergraduates (AIDep-22). PMC12862064 (verify year and authors before citing).</li></ul>';

  function mount(el) { root = el; css(); var start = function () { var first = STAGES.find(function (s) { return !st.done[s[0]]; }); go(first ? first[0] : 'res'); }; start(); var tries = 0, tryLoad = function () { loadGates().then(function (ok) { if (ok) { if (root === el && el.isConnected) start(); sync(); } else if (++tries < 4) setTimeout(tryLoad, 1500); }); }; tryLoad(); }
  // scales lab: register a view; statistics lab: route on #know
  if (typeof VIEWS !== 'undefined') VIEWS.know = mount;
  else { var route = function () { if ((location.hash || '').slice(1).split('/')[0] !== 'know') return; var v = document.getElementById('view'); if (v) mount(v); }; addEventListener('hashchange', function () { setTimeout(route, 0); }); if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', function () { setTimeout(route, 0); }); else setTimeout(route, 0); }
  return { mount: mount, state: function () { return st; } };
})();

/* ===== Reinforcement layer shared by both labs.
   1) Variable-ratio reinforcement (VR-5): after a correct answer, a scientist card appears with p = 0.2.
   2) Spaced repetition: a wrong answer comes back after 1 day, then 3, 7 and 14 days (Leitner-style boxes).
   3) Reinforcement-schedules lab: the learner responds under FR, VR, FI and VI and sees the cumulative record.
   Everything is stored in this browser only. ===== */
var REWARD = (function () {
  var KEY = 'reward-v1', D = 864e5, BOX = [1, 3, 7, 14];
  function load() { try { return JSON.parse(localStorage.getItem(KEY) || 'null') || { cards: {}, sr: {}, n: 0 }; } catch (e) { return { cards: {}, sr: {}, n: 0 }; } }
  var st = load();
  function save() { try { localStorage.setItem(KEY, JSON.stringify(st)); } catch (e) { } }
  var E = function (s) { return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); };
  var CARDS = [
    ['galton', 'فرانسيس جالتون', 1884, 'أسس القياس العقلي، وابتكر معامل الارتباط مفهومًا، ومعمل القياس الأنثروبومتري.', 'نادر'],
    ['pearson', 'كارل بيرسون', 1896, 'صاغ معامل الارتباط الذي يحمل اسمه، واختبار كاي تربيع.', 'شائع'],
    ['spearman', 'تشارلز سبيرمان', 1904, 'العامل العام g، وأساس التحليل العاملي، ومعادلة سبيرمان–براون.', 'نادر'],
    ['binet', 'ألفريد بينيه', 1905, 'أول اختبار ذكاء عملي مع تيودور سيمون.', 'شائع'],
    ['thurstone', 'لويس ثرستون', 1928, 'قياس الاتجاهات بالفترات المتساوية، والقدرات العقلية الأولية.', 'نادر'],
    ['likert', 'رنسيس ليكرت', 1932, 'مقياس التقديرات المجمّعة، الأكثر استخدامًا في العالم.', 'شائع'],
    ['guttman', 'لويس جتمان', 1944, 'التدرج التراكمي ومعامل إعادة الإنتاج.', 'نادر'],
    ['cronbach', 'لي كرونباخ', 1951, 'معامل ألفا، والصدق البنائي مع ميل (1955)، ونظرية التعميم.', 'شائع'],
    ['skinner', 'بي. إف. سكنر', 1957, 'جداول التعزيز مع فيرستر. وأنت الآن تحت جدول نسبة متغيرة!', 'أسطوري'],
    ['osgood', 'تشارلز أوزجود', 1957, 'التمايز الدلالي: التقييم والقوة والنشاط.', 'نادر'],
    ['campbell', 'دونالد كامبل', 1959, 'مصفوفة السمات والطرق المتعددة مع فيسك.', 'نادر'],
    ['rasch', 'جورج راش', 1960, 'نموذج راش، وأساس القياس الموضوعي.', 'نادر'],
    ['lord', 'فريدريك لورد', 1968, 'النظرية الحديثة للاختبارات ونظرية الاستجابة للمفردة.', 'نادر'],
    ['messick', 'صامويل ميسك', 1989, 'المفهوم الموحد للصدق وعواقب الاستخدام.', 'أسطوري'],
    ['kane', 'مايكل كين', 2013, 'حجة الصدق: التفسير والاستخدام كدعوى تحتاج أدلة.', 'أسطوري'],
    ['clarkwatson', 'كلارك وواتسون', 2019, '«بناء الصدق»: الدليل الحديث لبناء المقاييس.', 'أسطوري']];
  var BY = {}; CARDS.forEach(function (c) { BY[c[0]] = c; });
  var css = '.rw-card{position:fixed;inset-inline-start:50%;top:18%;transform:translateX(50%);z-index:9999;width:min(320px,90vw);border-radius:16px;padding:16px;color:#fff;background:linear-gradient(160deg,#3b2a8a,#1f6f6a);box-shadow:0 20px 60px rgba(0,0,0,.35);text-align:center}[dir=ltr] .rw-card{transform:translateX(-50%)}.rw-card b{display:block;font-size:20px;margin:6px 0}.rw-card small{opacity:.85}.rw-card .rw-r{display:inline-block;border-radius:999px;padding:0 10px;background:rgba(255,255,255,.2);font-size:12px}.rw-card button{margin-top:10px;border:0;border-radius:8px;padding:6px 14px;cursor:pointer;font:inherit}.rw-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:8px}.rw-c{border:1px solid var(--line,#ddd);border-radius:12px;padding:10px;font-size:13px}.rw-c.lock{opacity:.35}.rw-c b{display:block}.rw-sv{width:100%;height:auto;border:1px solid var(--line,#ddd);border-radius:10px}';
  function injectCss() { if (document.getElementById('rw-css')) return; var s = document.createElement('style'); s.id = 'rw-css'; s.textContent = css; document.head.appendChild(s); }
  function showCard(c, isNew) { injectCss(); var n = document.createElement('div'); n.className = 'rw-card'; n.setAttribute('role', 'dialog'); n.setAttribute('aria-label', 'كارت عالِم'); n.innerHTML = '<small>' + (isNew ? 'كارت جديد!' : 'كارت مكرر') + '</small><b>' + E(c[1]) + ' · ' + c[2] + '</b><span class="rw-r">' + E(c[4]) + '</span><p style="margin:8px 0 0">' + E(c[3]) + '</p><small>مجموعتك: ' + Object.keys(st.cards).length + ' / ' + CARDS.length + '</small><br><button>تمام</button>'; document.body.appendChild(n); var close = function () { n.remove(); }; n.querySelector('button').onclick = close; setTimeout(close, 6000); }
  function draw() { var r = Math.random(), pool = CARDS.filter(function (c) { return r < 0.1 ? c[4] === 'أسطوري' : r < 0.4 ? c[4] === 'نادر' : c[4] === 'شائع'; }); return pool[Math.floor(Math.random() * pool.length)]; }
  function onAnswer(id, ok) {
    st.n++;
    var s = st.sr[id];
    if (!ok) st.sr[id] = { b: 0, due: Date.now() + BOX[0] * D };
    else if (s) { var b = s.b + 1; if (b >= BOX.length) delete st.sr[id]; else st.sr[id] = { b: b, due: Date.now() + BOX[b] * D }; }
    var got = null;
    if (ok && Math.random() < 0.2) { var c = draw(), isNew = !st.cards[c[0]]; st.cards[c[0]] = (st.cards[c[0]] || 0) + 1; got = c; setTimeout(function () { showCard(c, isNew); }, 250); try { if (typeof addXP === 'function') addXP(isNew ? 5 : 2); } catch (e) { } }
    save(); return got;
  }
  function isDue(id) { var s = st.sr[id]; return !!s && s.due <= Date.now(); }
  function dueCount(ids) { st = load(); return ids.filter(isDue).length; }
  function inReview(id) { return !!st.sr[id]; }
  function collectionHTML() { injectCss(); return '<div class="rw-grid">' + CARDS.map(function (c) { var have = st.cards[c[0]]; return '<div class="rw-c ' + (have ? '' : 'lock') + '"><b>' + (have ? E(c[1]) : '؟؟؟') + '</b><span>' + c[2] + ' · ' + E(c[4]) + '</span>' + (have ? '<p style="margin:4px 0 0">' + E(c[3]) + '</p>' : '') + '</div>'; }).join('') + '</div><p class="muted small">جمعت ' + Object.keys(st.cards).length + ' من ' + CARDS.length + '. الكروت بتظهر بعد بعض الإجابات الصحيحة، من غير ما تعرف امتى: ده جدول «نسبة متغيرة» (VR-5). جرّب معمل جداول التعزيز علشان تعرف ليه الجدول ده بالذات.</p>'; }
  /* ---------- schedules lab ---------- */
  var SCH = [['FR', 'نسبة ثابتة FR-5', 'مكافأة بعد كل 5 استجابات بالظبط.'], ['VR', 'نسبة متغيرة VR-5', 'مكافأة بعد 5 استجابات في المتوسط، ومش معروف امتى بالظبط.'], ['FI', 'فترة ثابتة FI-4ث', 'أول استجابة بعد مرور 4 ثواني من آخر مكافأة.'], ['VI', 'فترة متغيرة VI-4ث', 'أول استجابة بعد فترة متوسطها 4 ثواني، ومتغيرة.']];
  function labMount(el) {
    injectCss(); var rec = {}, cur = null;
    el.innerHTML = '<div class="stack"><div class="section-title"><h2>معمل جداول التعزيز</h2><span class="chip">سكنر وفيرستر 1957</span></div><section class="card stack"><p style="margin:0">اختار جدول ودوس «استجب» زي ما إنت عايز لمدة 20 ثانية. البرنامج بيكافئك حسب الجدول، وبيرسم <b>السجل التراكمي</b>: كل ما الخط يطلع بسرعة يبقى معدل استجابتك أعلى، والعلامات الصغيرة هي المكافآت.</p><div class="row" id="rwS">' + SCH.map(function (s) { return '<button class="btn sm ghost" data-s="' + s[0] + '">' + s[1] + '</button>'; }).join('') + '</div><p class="small muted" id="rwD"></p><div class="row"><button class="btn" id="rwB" disabled style="min-width:160px;min-height:56px">استجب</button><span id="rwT" class="chip"></span><span id="rwR" class="chip"></span></div><svg id="rwG" class="rw-sv" viewBox="0 0 600 260" role="img" aria-label="السجل التراكمي"></svg><div id="rwX"></div></section><section class="card stack"><h3 style="margin:0">مجموعة كروت العلماء</h3>' + collectionHTML() + '</section></div>';
    var COL = { FR: '#d97706', VR: '#7c3aed', FI: '#0f766e', VI: '#2563eb' };
    function plot() { var g = '', W = 600, H = 260, P = 30, maxY = Math.max(20, Math.max.apply(null, Object.keys(rec).map(function (k) { return rec[k].r.length; }).concat([0]))); g += '<line x1="' + P + '" y1="' + (H - P) + '" x2="' + (W - 10) + '" y2="' + (H - P) + '" stroke="currentColor" stroke-opacity=".3"/><line x1="' + P + '" y1="10" x2="' + P + '" y2="' + (H - P) + '" stroke="currentColor" stroke-opacity=".3"/><text x="' + (W / 2) + '" y="' + (H - 6) + '" font-size="11" text-anchor="middle" fill="currentColor">الزمن (20 ثانية)</text><text x="12" y="' + (H / 2) + '" font-size="11" fill="currentColor" transform="rotate(-90 12 ' + (H / 2) + ')" text-anchor="middle">الاستجابات التراكمية</text>';
      Object.keys(rec).forEach(function (k) { var R = rec[k], pts = ['' + P + ',' + (H - P)]; R.r.forEach(function (t, i) { pts.push((P + t / 20000 * (W - P - 10)).toFixed(1) + ',' + (H - P - (i + 1) / maxY * (H - P - 10)).toFixed(1)); }); g += '<polyline fill="none" stroke="' + COL[k] + '" stroke-width="2" points="' + pts.join(' ') + '"/>'; R.w.forEach(function (i) { var t = R.r[i]; var x = P + t / 20000 * (W - P - 10), y = H - P - (i + 1) / maxY * (H - P - 10); g += '<line x1="' + x + '" y1="' + y + '" x2="' + (x + 6) + '" y2="' + (y + 6) + '" stroke="' + COL[k] + '" stroke-width="2"/>'; }); g += '<text x="' + (W - 14) + '" y="' + (18 + 14 * Object.keys(rec).indexOf(k)) + '" font-size="12" text-anchor="end" fill="' + COL[k] + '">' + k + ': ' + R.r.length + ' استجابة، ' + R.w.length + ' مكافأة</text>'; });
      el.querySelector('#rwG').innerHTML = g; }
    function start(k) {
      var s = SCH.find(function (x) { return x[0] === k; }); el.querySelector('#rwD').textContent = s[2]; var b = el.querySelector('#rwB'); b.disabled = false; rec[k] = { r: [], w: [] }; cur = k;
      var t0 = performance.now(), cnt = 0, need = k === 'VR' ? 1 + Math.floor(Math.random() * 9) : 5, last = 0, wait = k === 'VI' ? 1000 + Math.random() * 6000 : 4000;
      b.onclick = function () { var t = performance.now() - t0; if (t > 20000) return; var R = rec[k]; R.r.push(t); cnt++; var win = false;
        if (k === 'FR' || k === 'VR') { if (cnt >= need) { win = true; cnt = 0; need = k === 'VR' ? 1 + Math.floor(Math.random() * 9) : 5; } }
        else if (t - last >= wait) { win = true; last = t; wait = k === 'VI' ? 1000 + Math.random() * 6000 : 4000; }
        if (win) { R.w.push(R.r.length - 1); el.querySelector('#rwR').textContent = '⭐ مكافأة! (' + R.w.length + ')'; }
        plot(); };
      var tick = setInterval(function () { var t = performance.now() - t0, tt = el.querySelector('#rwT'); if (!tt) return clearInterval(tick); tt.textContent = '⏱ ' + Math.max(0, 20 - t / 1000).toFixed(0) + ' ث'; if (t > 20000) { clearInterval(tick); b.disabled = true; done(); } }, 200);
      plot();
    }
    function done() { var n = Object.keys(rec).length; el.querySelector('#rwX').innerHTML = '<div class="notice small"><b>إيه اللي المفروض يظهر؟</b>في تجارب سكنر: <b>VR</b> بيدّي أعلى معدل استجابة وأكثره ثباتًا، وأقوى مقاومة للانطفاء، وده سر ماكينات القمار وإشعارات السوشيال ميديا. <b>FR</b> بيدّي معدل عالي مع وقفة بعد كل مكافأة. <b>FI</b> بيدّي شكل «الصدفة» (scallop): بطء بعد المكافأة وتسارع قبل الموعد، زي المذاكرة قبل الامتحان. <b>VI</b> معدل متوسط ثابت، زي متابعة الرسايل.' + (n < 4 ? ' جرّب الأربع جداول وقارن.' : ' قارن خطوطك: هل ظهر النمط ده عندك؟ ولاحظ إن كروت العلماء في الموقع شغالة بجدول VR.') + '</div>'; }
    el.querySelectorAll('[data-s]').forEach(function (b) { b.onclick = function () { start(b.dataset.s); }; });
    plot();
  }
  if (typeof VIEWS !== 'undefined') VIEWS.schedules = labMount;
  else { var route = function () { if ((location.hash || '').slice(1) !== 'schedules') return; var v = document.getElementById('view'); if (v) labMount(v); }; addEventListener('hashchange', function () { setTimeout(route, 0); }); if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', function () { setTimeout(route, 0); }); else setTimeout(route, 0); }
  return { onAnswer: onAnswer, isDue: isDue, dueCount: dueCount, inReview: inReview, collectionHTML: collectionHTML, labMount: labMount };
})();

/* ===== Thesis tools for a collected scale (Scales Lab, builder step 9).
   1) Careless responding: long-string, intra-individual response variability (IRV), Mahalanobis D² (Meade & Craig, 2012;
      Curran, 2016), and the effect of excluding flagged cases on reliability.
   2) Total omega from a one-factor solution (McDonald, 1999) next to alpha.
   3) Norms: percentile ranks and normalised T scores for every observed raw score.
   4) Redundant items: high zero-order and partial correlations together (a simplified screen inspired by Unique Variable
      Analysis, Christensen, Garrido & Golino, 2023).
   5) Short form by ant colony optimisation (Leite, Huang & Marcoulides, 2008): maximises reliability, the correlation with
      the full score and facet coverage at the same time.
   6) An Arabic, APA-style Word report of the psychometric results. ===== */
var THESIS = (function () {
  var f = function (x, d) { return x == null || !isFinite(x) ? '—' : (+x).toFixed(d == null ? 2 : d); };
  var H = function (s) { return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); };
  function alphaOf(rows) { try { return PSY.itemStats(rows).alpha; } catch (e) { return NaN; } }
  function sub(rows, idx) { return rows.map(function (r) { return idx.map(function (j) { return r[j]; }); }); }
  /* ---------- 1) careless responding ---------- */
  function careless(rows) {
    var n = rows.length, k = rows[0].length, R = PSY.corrMatrix(rows), means = [], sds = [];
    for (var j = 0; j < k; j++) { var c = PSY.col(rows, j); means.push(PSY.mean(c)); sds.push(PSY.sd(c) || 1); }
    var Ri = null; try { Ri = PSY.inverse(R); } catch (e) { }
    var out = rows.map(function (r, i) {
      var ls = 1, run = 1; for (var j = 1; j < k; j++) { if (r[j] === r[j - 1]) { run++; ls = Math.max(ls, run); } else run = 1; }
      var irv = PSY.sd(r), d2 = NaN;
      if (Ri && n > k + 2) { var z = r.map(function (v, j) { return (v - means[j]) / sds[j]; }); d2 = 0; for (var a = 0; a < k; a++) for (var b = 0; b < k; b++) d2 += z[a] * Ri[a][b] * z[b]; }
      var p = isFinite(d2) ? PSY.chiP(d2, k) : NaN;
      var why = []; if (ls >= Math.max(6, Math.ceil(k * 0.8))) why.push('سلسلة إجابات متطابقة (' + ls + ')'); if (irv === 0) why.push('تباين صفري'); if (p < 0.001) why.push('مسافة ماهالانوبس متطرفة');
      return { i: i, ls: ls, irv: irv, d2: d2, p: p, why: why };
    });
    var bad = out.filter(function (x) { return x.why.length; }), keep = rows.filter(function (_, i) { return !out[i].why.length; });
    return { all: out, bad: bad, aBefore: alphaOf(rows), aAfter: keep.length > k ? alphaOf(keep) : NaN, nAfter: keep.length };
  }
  /* ---------- 2) omega total ---------- */
  function omegaTotal(rows) { try { if (rows[0].length < 3 || rows.length <= rows[0].length) return null; var e = PSY.efa(rows, 1); return PSY.omega(e.loadings.map(function (l) { return Math.abs(l[0]); })); } catch (e) { return null; } }
  /* ---------- 3) norms ---------- */
  function norms(tot) {
    var n = tot.length, s = tot.slice().sort(function (a, b) { return a - b; }), vals = Array.from(new Set(s));
    var inv = function (p) { var a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239], b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572], c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783], d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416], q, r; if (p < 0.02425) { q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); } if (p > 1 - 0.02425) { q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); } q = p - 0.5; r = q * q; return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1); };
    return vals.map(function (v) { var below = s.filter(function (x) { return x < v; }).length, eq = s.filter(function (x) { return x === v; }).length, pr = (below + eq / 2) / n; return { raw: v, pr: pr * 100, T: 50 + 10 * inv(Math.min(0.999, Math.max(0.001, pr))), n: eq }; });
  }
  /* ---------- 4) redundancy ---------- */
  function redundancy(rows, items) {
    var k = rows[0].length, R = PSY.corrMatrix(rows), P = null, out = []; try { P = PSY.inverse(R); } catch (e) { }
    for (var a = 0; a < k; a++) for (var b = a + 1; b < k; b++) { var r = R[a][b], pr = P ? -P[a][b] / Math.sqrt(P[a][a] * P[b][b]) : NaN; if (r >= 0.6 && (pr >= 0.4 || !isFinite(pr)) || r >= 0.8) out.push({ a: a, b: b, r: r, pr: pr }); }
    return out.sort(function (x, y) { return y.r - x.r; });
  }
  /* ---------- 5) short form by ant colony optimisation ---------- */
  function shortForm(rows, items, K) {
    var k = rows[0].length; K = Math.max(3, Math.min(k - 1, K)); var full = PSY.rowSums(rows), tau = items.map(function () { return 1; }), dims = items.map(function (it) { return it.dim || '_'; }), D = Array.from(new Set(dims));
    function score(idx) { var S = sub(rows, idx), a = alphaOf(S), r = PSY.pearson(PSY.rowSums(S), full), cov = D.filter(function (d) { return idx.some(function (j) { return dims[j] === d; }); }).length / D.length; if (!isFinite(a)) return -1; return 0.45 * a + 0.4 * r + 0.15 * cov; }
    var best = null, bestS = -Infinity, ANTS = 24, IT = 50;
    for (var it = 0; it < IT; it++) {
      var gen = [];
      for (var an = 0; an < ANTS; an++) { var pool = tau.map(function (t, j) { return [j, t]; }), pick = []; while (pick.length < K) { var tot = pool.reduce(function (s, x) { return s + x[1]; }, 0), r = Math.random() * tot, acc = 0, q = 0; for (; q < pool.length; q++) { acc += pool[q][1]; if (acc >= r) break; } q = Math.min(q, pool.length - 1); pick.push(pool[q][0]); pool.splice(q, 1); } pick.sort(function (a, b) { return a - b; }); var sc = score(pick); gen.push([pick, sc]); if (sc > bestS) { bestS = sc; best = pick; } }
      gen.sort(function (a, b) { return b[1] - a[1]; }); tau = tau.map(function (t) { return Math.max(0.05, t * 0.85); }); gen.slice(0, 4).forEach(function (g) { g[0].forEach(function (j) { tau[j] += g[1]; }); });
    }
    var S = sub(rows, best); return { idx: best, alpha: alphaOf(S), r: PSY.pearson(PSY.rowSums(S), full), aFull: alphaOf(rows), K: K };
  }
  /* ---------- UI ---------- */
  function panel(s, A) {
    if (!A || !A.rows || A.rows.length < 10 || A.items.length < 3) return '<details class="steps"><summary>أدوات الرسالة</summary><p class="muted small">تحتاج 10 استجابات مكتملة و3 عبارات على الأقل.</p></details>';
    return '<details class="steps" id="thBox"><summary>أدوات الرسالة: الاستجابات المتهاونة، وأوميجا، والمعايير، والتكرار، والصورة المختصرة، وتقرير Word</summary><div id="thOut" class="stack" style="margin-top:10px"><p class="muted small">افتح القسم لحساب النتائج.</p></div></details>';
  }
  function compute(s, A) {
    var rows = A.rows, items = A.items, k = items.length, C = careless(rows), om = omegaTotal(rows), tot = PSY.rowSums(rows), NR = norms(tot), RD = redundancy(rows, items), SF = shortForm(rows, items, Math.ceil(k / 2));
    return { C: C, om: om, NR: NR, RD: RD, SF: SF, n: rows.length, k: k };
  }
  function html(s, A, X) {
    var items = A.items, t = function (i) { return (i + 1) + '. ' + H(items[i].text || items[i].t || ''); };
    var h = '<h4 style="margin:0">1) الاستجابات العشوائية والمتهاونة</h4>' +
      (X.C.bad.length ? tableHTML(['المستجيب', 'أطول سلسلة متطابقة', 'التباين الفردي', 'D² (p)', 'السبب'], X.C.bad.slice(0, 40).map(function (b) { return [String(b.i + 1), String(b.ls), f(b.irv), f(b.d2, 1) + ' (' + (b.p < 0.001 ? '<.001' : f(b.p, 3)) + ')', H(b.why.join('، '))]; })) : '<div class="notice good small">لا توجد استجابات مشتبه بها وفق المؤشرات الثلاثة.</div>') +
      '<p class="small" style="margin:0">ألفا قبل الاستبعاد = ' + f(X.C.aBefore, 3) + ' (ن = ' + X.n + ')، وبعد استبعاد ' + X.C.bad.length + ' حالة = ' + f(X.C.aAfter, 3) + ' (ن = ' + X.C.nAfter + '). لا تحذف آليًا: راجع كل حالة، واذكر المعيار في الرسالة.</p>' +
      '<h4 style="margin:0">2) أوميجا الكلي</h4><p class="small" style="margin:0">ω = ' + f(X.om ? X.om.omega : NaN, 3) + ' مقابل α = ' + f(X.C.aBefore, 3) + '. أوميجا لا تفترض تساوي التشبعات (التكافؤ التاوي) كما تفترضها ألفا، ولذلك يوصى بتقريرهما معًا (McDonald, 1999؛ Dunn et al., 2014).</p>' +
      '<h4 style="margin:0">3) جدول المعايير (الدرجة الكلية)</h4>' + tableHTML(['الدرجة الخام', 'عدد الأفراد', 'الرتبة المئينية', 'الدرجة التائية المعدلة'], X.NR.map(function (r) { return [f(r.raw, 0), String(r.n), f(r.pr, 1), f(r.T, 0)]; })) + '<p class="muted small" style="margin:0">معايير مؤقتة من عينتك فقط. المعايير الحقيقية تحتاج عينة ممثلة وكبيرة (عادة 500 فأكثر) ومقسمة حسب السن والنوع.</p>' +
      '<h4 style="margin:0">4) عبارات قد تكون مكررة</h4>' + (X.RD.length ? tableHTML(['العبارة الأولى', 'العبارة الثانية', 'r', 'r الجزئي'], X.RD.slice(0, 15).map(function (d) { return [t(d.a), t(d.b), f(d.r, 2), f(d.pr, 2)]; })) + '<p class="muted small" style="margin:0">ارتباط مرتفع مع ارتباط جزئي مرتفع يعني غالبًا أن العبارتين تقولان الشيء نفسه بصياغتين، فيرتفع الثبات بشكل زائف. احتفظ بالأوضح منهما.</p>' : '<div class="notice good small">لا توجد أزواج عبارات مكررة واضحة.</div>') +
      '<h4 style="margin:0">5) الصورة المختصرة (خوارزمية مستعمرة النمل)</h4><p class="small" style="margin:0">أفضل ' + X.SF.K + ' عبارات: ' + X.SF.idx.map(function (j) { return j + 1; }).join('، ') + '. ألفا المختصرة = ' + f(X.SF.alpha, 3) + ' (الكاملة ' + f(X.SF.aFull, 3) + ')، والارتباط بالدرجة الكاملة r = ' + f(X.SF.r, 3) + '.</p><p class="muted small" style="margin:0">الخوارزمية بتوازن الثبات والارتباط بالصورة الكاملة وتغطية الأبعاد. لازم تتأكد من الصورة المختصرة على عينة مستقلة، لأن اختيارها من نفس البيانات بيرفع النتايج بشكل متفائل.</p>' +
      '<div class="row"><button class="btn sm" id="thDoc">' + ic('download') + ' تقرير الخصائص السيكومترية (Word)</button><button class="btn sm ghost" id="thRe">أعد الحساب</button></div>';
    return h;
  }
  function report(s, A, X) {
    var items = A.items, st = PSY.itemStats(A.rows);
    var rows = items.map(function (it, j) { return '<tr><td>' + (j + 1) + '</td><td>' + H(it.text || '') + '</td><td>' + f(st.im ? st.im[j] : NaN) + '</td><td>' + f(st.isd ? st.isd[j] : NaN) + '</td><td>' + f(st.rit ? st.rit[j] : NaN) + '</td><td>' + f(st.aid ? st.aid[j] : NaN, 3) + '</td></tr>'; }).join('');
    var dims = (A.dims || []).filter(function (d) { return d.st; }).map(function (d) { return '<tr><td>' + H(d.d.name) + '</td><td>' + d.n + '</td><td>' + f(d.st.alpha, 3) + '</td><td>' + f(d.omega ? d.omega.omega : NaN, 3) + '</td></tr>'; }).join('');
    var doc = '<html dir="rtl" lang="ar"><head><meta charset="utf-8"><title>الخصائص السيكومترية</title><style>body{font-family:"Simplified Arabic","Traditional Arabic",Arial;font-size:14pt;line-height:1.8}h1{font-size:18pt;text-align:center}h2{font-size:15pt}table{border-collapse:collapse;width:100%;margin:8pt 0}td,th{border:1px solid #000;padding:4pt;font-size:12pt}caption{font-weight:bold;text-align:right;padding:4pt 0}</style></head><body>' +
      '<h1>الخصائص السيكومترية لمقياس «' + H(s.name || 'المقياس') + '»</h1>' +
      '<h2>أولًا: العينة</h2><p>طُبق المقياس على عينة استطلاعية قوامها (' + A.n + ') مستجيبًا. وفُحصت الاستجابات العشوائية والمتهاونة بثلاثة مؤشرات: أطول سلسلة إجابات متطابقة، والتباين الفردي في الاستجابة، ومسافة ماهالانوبس (Meade & Craig, 2012)، فتبيّن وجود (' + X.C.bad.length + ') حالة مشتبه بها.</p>' +
      '<h2>ثانيًا: وصف المقياس</h2><p>يتكون المقياس في صورته الحالية من (' + items.length + ') عبارة' + ((s.dims || []).length ? ' موزعة على (' + s.dims.length + ') أبعاد هي: ' + s.dims.map(function (d) { return H(d.name); }).join('، ') : '') + '.</p>' +
      '<h2>ثالثًا: الاتساق الداخلي وتحليل العبارات</h2><table><caption>جدول (1): المتوسطات والانحرافات المعيارية ومعاملات ارتباط العبارة بالدرجة الكلية المصححة وألفا عند حذف العبارة</caption><tr><th>م</th><th>العبارة</th><th>م</th><th>ع</th><th>ر المصحح</th><th>α عند الحذف</th></tr>' + rows + '</table>' +
      '<h2>رابعًا: الثبات</h2><p>بلغ معامل ألفا كرونباخ للمقياس الكلي (' + f(X.C.aBefore, 3) + ')، ومعامل أوميجا ماكدونالد (' + f(X.om ? X.om.omega : NaN, 3) + ')' + (A.split ? '، ومعامل التجزئة النصفية بمعادلة سبيرمان–براون (' + f(A.split.sb, 3) + ')' : '') + '. ويُعد أوميجا أدق من ألفا لأنه لا يفترض تساوي تشبعات العبارات (Dunn, Baguley & Brunsden, 2014).</p>' + (dims ? '<table><caption>جدول (2): ثبات الأبعاد</caption><tr><th>البعد</th><th>العبارات</th><th>α</th><th>ω</th></tr>' + dims + '</table>' : '') +
      '<h2>خامسًا: الصورة المختصرة</h2><p>باستخدام خوارزمية مستعمرة النمل (Leite et al., 2008) اختيرت (' + X.SF.K + ') عبارات هي: ' + X.SF.idx.map(function (j) { return j + 1; }).join('، ') + '، وبلغ ثباتها (' + f(X.SF.alpha, 3) + ') وارتباطها بالصورة الكاملة (' + f(X.SF.r, 3) + '). وتحتاج إلى تأكيد على عينة مستقلة.</p>' +
      '<h2>سادسًا: المعايير المؤقتة</h2><table><caption>جدول (3): الدرجات الخام والرتب المئينية والدرجات التائية المعدلة</caption><tr><th>الدرجة الخام</th><th>الرتبة المئينية</th><th>الدرجة التائية</th></tr>' + X.NR.map(function (r) { return '<tr><td>' + f(r.raw, 0) + '</td><td>' + f(r.pr, 1) + '</td><td>' + f(r.T, 0) + '</td></tr>'; }).join('') + '</table>' +
      '<p style="font-size:11pt">ملاحظة: أُعدّ هذا التقرير آليًا في مختبر المقاييس النفسية كمسودة، ويجب مراجعة الأرقام والصياغة، واستكمال أدلة الصدق (التحليل العاملي التوكيدي، والصدق التقاربي والتمييزي) قبل اعتماده في الرسالة.</p>' +
      '<h2>المراجع</h2><p dir="ltr" style="text-align:left;font-size:11pt">Dunn, T. J., Baguley, T., & Brunsden, V. (2014). From alpha to omega. British Journal of Psychology, 105(3), 399–412.<br>Leite, W. L., Huang, I.-C., & Marcoulides, G. A. (2008). Item selection for the development of short forms of scales using an ant colony optimization algorithm. Multivariate Behavioral Research, 43(3), 411–431.<br>McDonald, R. P. (1999). Test theory: A unified treatment. Erlbaum.<br>Meade, A. W., & Craig, S. B. (2012). Identifying careless responses in survey data. Psychological Methods, 17(3), 437–455.</p></body></html>';
    downloadFile('psychometric-report.doc', doc, 'application/msword');
  }
  function bind(el, s, A) {
    var box = el.querySelector('#thBox'); if (!box) return; var X = null;
    var run = function () { var o = el.querySelector('#thOut'); o.innerHTML = '<p class="muted small">جارٍ الحساب…</p>'; setTimeout(function () { try { X = compute(s, A); o.innerHTML = html(s, A, X); o.querySelector('#thDoc').onclick = function () { report(s, A, X); }; o.querySelector('#thRe').onclick = run; } catch (e) { o.innerHTML = '<div class="notice bad small">تعذر الحساب: ' + H(e.message) + '</div>'; } }, 30); };
    box.addEventListener('toggle', function () { if (box.open && !X) run(); });
  }
  return { panel: panel, bind: bind, careless: careless, norms: norms, redundancy: redundancy, shortForm: shortForm, omegaTotal: omegaTotal };
})();

/* ===== Admin page (#admin), shared by both labs (firebase accounts mode).
   The admin is whoever has a document admins/{uid} (created by hand in the Firebase console; see the setup card).
   Features: global stage gates for the know-yourself journey (step by step, fully open, open up to a stage, closed,
   automatic hold after high wellbeing scores), per-student overrides (open the next stage, release after review),
   the list of students who shared their results by name, flags, details and CSV export. ===== */
var ADMIN = (function () {
  var E = function (s) { return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); };
  var say = function (m) { try { toast(m); } catch (e) { } };
  var dl = function (name, text, mime) { if (typeof downloadFile === 'function') return downloadFile(name, text, mime); var a = document.createElement('a'); a.href = URL.createObjectURL(new Blob(['\ufeff' + text], { type: mime })); a.download = name; document.body.appendChild(a); a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 500); };
  var ST = ['البداية والموافقة', 'المزاج', 'الذاكرة', 'التركيز', 'الحالة النفسية', 'حياتك الرقمية (M6)', 'أنت والذكاء الاصطناعي (M2)', 'النتائج'];
  var MODES = [['seq', 'مرحلة مرحلة تلقائيًا (كل مرحلة تفتح بعد اللي قبلها)'], ['open', 'كل المراحل مفتوحة'], ['upto', 'مفتوحة لحد مرحلة معينة'], ['locked', 'مقفولة']];
  var CSS = '.ad{display:flex;flex-direction:column;gap:14px;max-width:1100px;margin:0 auto}.ad-card{border:1px solid var(--line,#ddd);border-radius:14px;padding:16px;background:var(--surface,#fff);display:flex;flex-direction:column;gap:10px}.ad-card h3{margin:0}.ad-row{display:flex;flex-wrap:wrap;gap:8px;align-items:center}.ad select,.ad input[type=text]{border:1px solid var(--line,#ccc);border-radius:8px;padding:6px 8px;font:inherit;background:var(--surface,#fff);color:inherit}.ad-btn{border:1px solid var(--primary,var(--accent,#5b3fd1));background:var(--primary,var(--accent,#5b3fd1));color:#fff;border-radius:9px;padding:6px 12px;font:inherit;cursor:pointer}.ad-btn.g{background:transparent;color:var(--primary,var(--accent,#5b3fd1))}.ad-btn.r{border-color:#b3261e;background:transparent;color:#b3261e}.ad-tw{overflow-x:auto}.ad table{border-collapse:collapse;width:100%;min-width:900px;font-size:14px}.ad th,.ad td{border-bottom:1px solid var(--line,#e5e5e5);padding:7px 8px;text-align:start;vertical-align:top}.ad tr.flag td{background:rgba(179,38,30,.07)}.ad-chip{display:inline-block;border-radius:999px;padding:0 9px;font-size:12px}.ad-chip.r{background:#fbe3e1;color:#b3261e}.ad-chip.g{background:#e1f2e6;color:#1e7a3c}.ad-chip.y{background:#fbf0d5;color:#8a5a00}.ad code{direction:ltr;unicode-bidi:isolate;background:rgba(0,0,0,.06);padding:1px 6px;border-radius:6px;word-break:break-all}.ad ol li{line-height:1.9}';
  function css() { if (document.getElementById('ad-css')) return; var s = document.createElement('style'); s.id = 'ad-css'; s.textContent = CSS; document.head.appendChild(s); }
  var fsOK = function () { try { return !!(window.ACCT && ACCT.fs && ACCT.fs.ready()); } catch (e) { return false; } };
  function flags(r) { var w = r.well || {}, out = []; if (w.phq >= 10) out.push('PHQ-8 = ' + w.phq); if (w.gad >= 10) out.push('GAD-7 = ' + w.gad); if (w.who != null && w.who * 4 <= 28) out.push('WHO-5 = ' + (w.who * 4) + '%'); return out; }
  var pct = function (x) { return x == null || !isFinite(x) ? '—' : Math.round(x) + '%'; };
  function setupCard(uid, err) {
    return '<section class="ad-card"><h3>تفعيل صلاحية المسؤول (مرة واحدة)</h3>' + (err ? '<p style="margin:0;color:#b3261e">' + E(err) + '</p>' : '') + '<p style="margin:0">حسابك لسه مش مسؤول. الصلاحية بتتعمل يدويًا من Firebase علشان محدش غيرك يقدر يعملها، ومن غير ما إيميلك يظهر في الكود العام:</p><ol style="margin:0"><li>افتح Firebase Console ← مشروع الموقع ده ← <b>Firestore Database</b> ← تبويب <b>Rules</b>، والصق ملف <code>firestore.rules</code> الموجود في GitHub، ودوس <b>Publish</b>.</li><li>في تبويب <b>Data</b> دوس <b>Start collection</b> واكتب الاسم <code>admins</code>.</li><li>في <b>Document ID</b> الصق المعرّف ده: <code id="adUid">' + E(uid) + '</code> <button class="ad-btn g" id="adCp">نسخ</button></li><li>ضيف حقل اسمه <code>role</code> نوعه string وقيمته <code>admin</code>، ودوس <b>Save</b>، وبعدين ارجع حدّث الصفحة دي.</li></ol><p class="muted" style="margin:0;font-size:13px">كل موقع ليه مشروع Firebase منفصل، فاعمل الخطوات دي في المشروعين لو عايز اللوحة في الاتنين.</p></section>';
  }
  function gateForm(g) {
    g = g || { mode: 'seq', upto: 7, hold: false, note: '' };
    return '<section class="ad-card"><h3>إعدادات الرحلة لكل الطلاب</h3><div class="ad-row"><label>طريقة الفتح <select id="gM">' + MODES.map(function (m) { return '<option value="' + m[0] + '"' + (g.mode === m[0] ? ' selected' : '') + '>' + m[1] + '</option>'; }).join('') + '</select></label><label id="gUL">لحد مرحلة <select id="gU">' + ST.map(function (t, i) { return '<option value="' + i + '"' + ((g.upto != null ? g.upto : 7) === i ? ' selected' : '') + '>' + (i + 1) + '. ' + t + '</option>'; }).join('') + '</select></label></div><label class="ad-row"><input type="checkbox" id="gH" ' + (g.hold ? 'checked' : '') + '> وقّف الطالب بعد «الحالة النفسية» لو نتايجه مرتفعة (PHQ-8 أو GAD-7 من 10 فأكثر، أو WHO-5 من 28% فأقل) لحد ما أراجعها وأسمح له</label><label>رسالة تظهر للطلاب (اختياري) <input type="text" id="gN" maxlength="300" style="width:100%" value="' + E(g.note || '') + '"></label><div><button class="ad-btn" id="gS">احفظ للكل</button></div></section>';
  }
  function studentRow(k, gt) {
    var r = {}; try { r = JSON.parse(k.results || '{}'); } catch (e) { }
    var w = r.well || {}, F = flags(r), g = gt || {}, mode = g.mode || 'inherit', when = k.updated && k.updated.toDate ? k.updated.toDate().toLocaleString('ar-EG') : '';
    var sel = '<select data-m="' + k.id + '"><option value="inherit"' + (mode === 'inherit' ? ' selected' : '') + '>زي الإعداد العام</option>' + MODES.map(function (m) { return '<option value="' + m[0] + '"' + (mode === m[0] ? ' selected' : '') + '>' + m[1].split(' (')[0] + '</option>'; }).join('') + '</select> <select data-u="' + k.id + '">' + ST.map(function (t, i) { return '<option value="' + i + '"' + ((g.upto != null ? g.upto : 7) === i ? ' selected' : '') + '>' + (i + 1) + '</option>'; }).join('') + '</select>';
    return '<tr class="' + (F.length ? 'flag' : '') + '"><td><b>' + E(k.name || '—') + '</b><br><span style="direction:ltr;unicode-bidi:isolate;font-size:12px">' + E(k.email || '') + '</span><br><span class="muted" style="font-size:12px">' + E(when) + '</span></td><td>' + (k.done || 0) + ' / 8</td><td>' + (w.who != null ? w.who * 4 + '%' : '—') + '</td><td>' + (w.phq != null ? w.phq : '—') + '</td><td>' + (w.gad != null ? w.gad : '—') + '</td><td>' + pct(r.m6 && r.m6.combined) + '</td><td>' + pct(r.m2 && r.m2.combined) + '</td><td>' + (F.length ? '<span class="ad-chip r">تحتاج مراجعة</span><br><span style="font-size:12px">' + E(F.join('، ')) + '</span>' : '<span class="ad-chip g">لا تنبيه</span>') + (g.released ? '<br><span class="ad-chip y">مسموح بعد المراجعة</span>' : '') + '</td><td><div class="ad-row">' + sel + '</div><div class="ad-row" style="margin-top:4px"><button class="ad-btn" data-next="' + k.id + '" data-done="' + (k.done || 0) + '">افتح المرحلة التالية</button><button class="ad-btn g" data-rel="' + k.id + '">' + (g.released ? 'الغِ السماح' : 'اسمح بعد المراجعة') + '</button><button class="ad-btn g" data-save="' + k.id + '">احفظ</button></div><details style="margin-top:4px"><summary>التفاصيل</summary><pre style="white-space:pre-wrap;direction:ltr;font-size:11px;max-height:220px;overflow:auto">' + E(JSON.stringify(r, null, 1)) + '</pre></details><button class="ad-btn r" data-del="' + k.id + '" style="margin-top:4px">امسح نتايجه</button></td></tr>';
  }
  var VIS = 'body:not(.is-admin) [data-go="admin"],body:not(.is-admin) [data-m="admin"],body:not(.is-admin) a[href="#admin"],body:not(.is-admin) button[onclick="go(\'admin\')"]{display:none!important}';
  (function () { var s = document.createElement('style'); s.textContent = VIS; document.head.appendChild(s); })();
  var known = null, lastWhy = '';
  async function ensureAdmin(tryCreate) {
    if (!fsOK()) return null; var uid = ACCT.fs.uid();
    try { var d = await ACCT.fs.get('admins/' + uid); if (d) { mark(true); return d; } } catch (e) { }
    if (!tryCreate) { mark(false); return null; }
    try { if (ACCT.fs.refresh) await ACCT.fs.refresh(); } catch (e) { }
    lastWhy = ''; if (!ACCT.fs.verified()) { lastWhy = 'unverified'; mark(false); return null; }
    try { var em = (ACCT.fs.email() || '').toLowerCase(), doc = { role: 'admin', email: em, name: (ACCT.fs.name() || '').slice(0, 60), added: '__ts__' }; await ACCT.fs.put('admins/' + uid, doc); try { await ACCT.fs.del('invites/' + em); } catch (e) { } mark(true); say('بقيت مسؤول ✓'); return doc; } catch (e) { lastWhy = /permission/i.test((e && (e.code || e.message)) || '') ? 'perm' : 'other:' + ((e && (e.code || e.message)) || ''); mark(false); return null; }
  }
  function mark(v) { known = v; document.body.classList.toggle('is-admin', !!v); }
  (function poll(n) { if (fsOK()) { var u = ACCT.fs.uid(); if (u !== poll.u) { poll.u = u; ensureAdmin(false); } } else if (known) mark(false); setTimeout(function () { poll(n + 1); }, n < 20 ? 1500 : 5000); })(0);
  function adminsCard(admins, invites, me) {
    return '<section class="ad-card"><h3>المسؤولين والدعوات</h3>' + '<div class="ad-tw"><table style="min-width:0"><thead><tr><th>المسؤول</th><th>الإيميل</th><th></th></tr></thead><tbody>' + admins.map(function (a) { return '<tr><td>' + E(a.name || '—') + (a.id === me ? ' <span class="ad-chip g">أنت</span>' : '') + '</td><td style="direction:ltr;unicode-bidi:isolate">' + E(a.email || '') + '</td><td>' + (a.id === me ? '' : '<button class="ad-btn r" data-rmadm="' + E(a.id) + '">شيل الصلاحية</button>') + '</td></tr>'; }).join('') + '</tbody></table></div>' +
      (invites.length ? '<p style="margin:0"><b>دعوات لسه متقبلتش:</b> ' + invites.map(function (v) { return '<span style="direction:ltr;unicode-bidi:isolate">' + E(v.id) + '</span> <button class="ad-btn g" data-rminv="' + E(v.id) + '">الغِ</button>'; }).join(' · ') + '</p>' : '') +
      '<div class="ad-row"><input type="text" id="invE" dir="ltr" placeholder="email@example.com" style="flex:1;min-width:200px" aria-label="إيميل المسؤول الجديد"><button class="ad-btn" id="invGo">ابعت دعوة مسؤول</button></div><p class="muted" style="margin:0;font-size:13px">الدعوة بتتسجل في قاعدة البيانات. ابعت للشخص لينك اللوحة ده: <code>' + E(location.origin + location.pathname + '#admin') + '</code>، ولما يسجّل دخول بنفس الإيميل (بعد تأكيده) ويفتح اللينك، بيبقى مسؤول تلقائيًا.</p></section>';
  }
  async function render(el) {
    css(); el.innerHTML = '<div class="ad"><h2 style="margin:0">لوحة المسؤول</h2><div id="adB"><p class="muted">جارٍ التحميل…</p></div></div>';
    var B = el.querySelector('#adB');
    if (!(window.ACCT && ACCT.fs)) { B.innerHTML = '<section class="ad-card"><p style="margin:0">لوحة المسؤول بتشتغل على نسخة GitHub من الموقع (بحسابات Firebase)، مش جوه Claude.</p></section>'; return; }
    if (!fsOK()) { B.innerHTML = '<section class="ad-card"><p style="margin:0">سجّل دخولك الأول من «حسابي» بنفس الإيميل اللي هتعمله مسؤول، وبعدين ارجع للصفحة دي. اللوحة بتشتغل على نسخة GitHub من الموقع: ' + (/github\.io$/.test(location.hostname) ? '' : '<a href="https://mo-psych-stats-lab.github.io/' + (document.title.indexOf('الإحصاء') >= 0 ? 'psych-stats-lab' : 'psych-scales-lab') + '/#admin" target="_blank" rel="noopener">افتحها من هنا</a>') + '</p></section>'; return; }
    var uid = ACCT.fs.uid(), adm = await ensureAdmin(true), err = '';
    if (!adm) { var why = lastWhy === 'unverified' ? 'إيميلك لسه مش مؤكد. أكّده من «حسابي» (ولو اللينك بيفتح صفحة خطأ استخدم خانة «اللينك بيفتح صفحة خطأ؟»)، وبعدين ارجع هنا ودوس «حاول تاني».' : lastWhy === 'perm' ? 'قاعدة البيانات رفضت: إما إن ملف القواعد الجديد (اللي فيه إيميلك) لسه متنشرش في Firestore ← Rules ← Publish في مشروع الموقع ده، أو إنك داخل بإيميل غير الإيميل المسؤول، أو مفيش دعوة لإيميلك.' : lastWhy ? 'خطأ: ' + E(lastWhy.replace(/^other:/, '')) : 'حسابك مش مسؤول.';
      B.innerHTML = '<section class="ad-card"><h3>الصفحة دي للمسؤول بس</h3><p style="margin:0">' + why + '</p><div class="ad-tw"><table style="min-width:0"><tr><th>داخل بإيميل</th><td style="direction:ltr;unicode-bidi:isolate">' + E(ACCT.fs.email()) + '</td></tr><tr><th>الإيميل مؤكد؟</th><td>' + (ACCT.fs.verified() ? 'أيوه ✓' : 'لأ ✗') + '</td></tr><tr><th>المشروع</th><td style="direction:ltr;unicode-bidi:isolate">' + E((window.FIREBASE_CONFIG || {}).projectId || '') + '</td></tr></table></div><div><button class="ad-btn" id="adTry">حاول تاني</button></div></section>'; B.querySelector('#adTry').onclick = function () { render(el); }; return; }
    var admins = [], invites = []; try { admins = await ACCT.fs.list('admins'); } catch (e) { } try { invites = await ACCT.fs.list('invites'); } catch (e) { }
    var know = [], gates = [], global = null;
    try { know = await ACCT.fs.list('know'); } catch (e) { err = 'تعذر قراءة النتائج: ' + (e.message || ''); }
    try { gates = await ACCT.fs.list('gates'); } catch (e) { }
    global = gates.find(function (g) { return g.id === 'global'; }) || null;
    var GT = {}; gates.forEach(function (g) { GT[g.id] = g; });
    var onlyF = false;
    var draw = function () {
      var list = know.slice().sort(function (a, b) { return (b.updated && b.updated.seconds || 0) - (a.updated && a.updated.seconds || 0); }); var nf = list.filter(function (k) { try { return flags(JSON.parse(k.results || '{}')).length; } catch (e) { return false; } }).length;
      if (onlyF) list = list.filter(function (k) { try { return flags(JSON.parse(k.results || '{}')).length; } catch (e) { return false; } });
      B.innerHTML = (err ? '<p style="color:#b3261e">' + E(err) + '</p>' : '') + adminsCard(admins, invites, uid) + gateForm(global) +
        '<section class="ad-card"><div class="ad-row" style="justify-content:space-between"><h3>الطلاب اللي شاركوا نتايجهم (' + know.length + ')</h3><div class="ad-row"><span class="ad-chip r">تحتاج مراجعة: ' + nf + '</span><label class="ad-row"><input type="checkbox" id="adF" ' + (onlyF ? 'checked' : '') + '> اللي محتاجين مراجعة بس</label><button class="ad-btn g" id="adCsv">تصدير CSV</button><button class="ad-btn g" id="adRe">تحديث</button></div></div>' +
        (list.length ? '<div class="ad-tw"><table><thead><tr><th>الطالب</th><th>المراحل</th><th>WHO-5</th><th>PHQ-8</th><th>GAD-7</th><th>M6</th><th>M2</th><th>التنبيه</th><th>التحكم في مراحله</th></tr></thead><tbody>' + list.map(function (k) { return studentRow(k, GT[k.id]); }).join('') + '</tbody></table></div>' : '<p class="muted" style="margin:0">لسه مفيش طالب شارك نتايجه. المشاركة اختيارية، والطالب بيوافق عليها في أول الرحلة أو في صفحة النتائج.</p>') +
        '<p class="muted" style="margin:0;font-size:13px">النتايج دي فرز مش تشخيص. التنبيه معناه «يستاهل تتواصل معاه»، ولو فيه خطر واضح اتبع بروتوكول الأزمة مع المشرف المرخّص.</p></section>';
      var gM = B.querySelector('#gM'), uL = B.querySelector('#gUL'), sh = function () { uL.style.display = gM.value === 'upto' ? '' : 'none'; }; gM.onchange = sh; sh();
      B.querySelector('#gS').onclick = async function () { var d = { mode: gM.value, upto: +B.querySelector('#gU').value, hold: B.querySelector('#gH').checked, released: false, note: B.querySelector('#gN').value.slice(0, 300) }; try { await ACCT.fs.set('gates/global', d); global = Object.assign({ id: 'global' }, d); say('اتحفظ للكل'); } catch (e) { say('تعذر الحفظ: ' + (e.message || '')); } };
      B.querySelector('#invGo').onclick = async function () { var em = B.querySelector('#invE').value.trim().toLowerCase(); if (!/^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(em)) return say('اكتب إيميل صحيح'); try { await ACCT.fs.put('invites/' + em, { by: uid, created: '__ts__' }); invites.push({ id: em }); say('اتسجلت الدعوة. ابعت له لينك اللوحة'); draw(); } catch (e) { say('تعذر إرسال الدعوة: ' + (e.message || '')); } };
      B.querySelectorAll('[data-rminv]').forEach(function (b) { b.onclick = async function () { try { await ACCT.fs.del('invites/' + b.dataset.rminv); invites = invites.filter(function (v) { return v.id !== b.dataset.rminv; }); draw(); } catch (e) { say('تعذر الإلغاء'); } }; });
      B.querySelectorAll('[data-rmadm]').forEach(function (b) { b.onclick = async function () { if (b.dataset.sure !== '1') { b.dataset.sure = '1'; b.textContent = 'متأكد؟ دوس تاني'; return; } try { await ACCT.fs.del('admins/' + b.dataset.rmadm); admins = admins.filter(function (v) { return v.id !== b.dataset.rmadm; }); draw(); } catch (e) { say('تعذر'); } }; });
      B.querySelector('#adF').onchange = function (e) { onlyF = e.target.checked; draw(); };
      B.querySelector('#adRe').onclick = function () { render(el); };
      B.querySelector('#adCsv').onclick = function () { var rows = [['name', 'email', 'done', 'who5_pct', 'phq8', 'gad7', 'm6_direct', 'm6_indirect', 'm6_combined', 'm2_direct', 'm2_indirect', 'm2_combined', 'digit_span', 'stroop_ms', 'probe_opened', 'cal_wrong_accepted']]; know.forEach(function (k) { var r = {}; try { r = JSON.parse(k.results || '{}'); } catch (e) { } var w = r.well || {}, sp = r.span || {}; rows.push([k.name, k.email, k.done, w.who != null ? w.who * 4 : '', w.phq, w.gad, r.m6 && r.m6.direct, r.m6 && r.m6.indirect, r.m6 && r.m6.combined, r.m2 && r.m2.direct, r.m2 && r.m2.indirect, r.m2 && r.m2.combined, sp.best, r.stroop && r.stroop.interf, sp.probe ? (sp.probe.opened ? 1 : 0) : '', r.m2 && r.m2.cal ? r.m2.cal.wrongAcc : ''].map(function (v) { return v == null ? '' : typeof v === 'number' ? Math.round(v * 100) / 100 : '"' + String(v).replace(/"/g, '""') + '"'; })); }); dl('know-results.csv', rows.map(function (r) { return r.join(','); }).join('\n'), 'text/csv;charset=utf-8'); };
      var saveG = async function (id, patch) { var cur = GT[id] || { mode: 'inherit', upto: 7, hold: false, released: false, note: '' }, d = Object.assign({ mode: cur.mode || 'inherit', upto: cur.upto != null ? cur.upto : 7, hold: false, released: !!cur.released, note: cur.note || '' }, patch); try { await ACCT.fs.set('gates/' + id, d); GT[id] = Object.assign({ id: id }, d); say('اتحفظ'); draw(); } catch (e) { say('تعذر الحفظ: ' + (e.message || '')); } };
      B.querySelectorAll('[data-save]').forEach(function (b) { b.onclick = function () { var id = b.dataset.save; saveG(id, { mode: B.querySelector('[data-m="' + id + '"]').value, upto: +B.querySelector('[data-u="' + id + '"]').value }); }; });
      B.querySelectorAll('[data-next]').forEach(function (b) { b.onclick = function () { var id = b.dataset.next, cur = GT[id], base = cur && cur.mode === 'upto' ? cur.upto + 1 : +b.dataset.done; saveG(id, { mode: 'upto', upto: Math.min(7, Math.max(1, base)), released: true }); }; });
      B.querySelectorAll('[data-rel]').forEach(function (b) { b.onclick = function () { var id = b.dataset.rel; saveG(id, { released: !(GT[id] && GT[id].released) }); }; });
      B.querySelectorAll('[data-del]').forEach(function (b) { b.onclick = async function () { var id = b.dataset.del; if (b.dataset.sure !== '1') { b.dataset.sure = '1'; b.textContent = 'متأكد؟ دوس تاني للمسح'; return; } try { await ACCT.fs.del('know/' + id); know = know.filter(function (k) { return k.id !== id; }); say('اتمسح'); draw(); } catch (e) { say('تعذر المسح'); } }; });
    };
    draw();
  }
  if (typeof VIEWS !== 'undefined') VIEWS.admin = function (el) { render(el); };
  else { var route = function () { if ((location.hash || '').slice(1).split('/')[0] !== 'admin') return; var v = document.getElementById('view'); if (v) render(v); }; addEventListener('hashchange', function () { setTimeout(route, 0); }); if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', function () { setTimeout(route, 0); }); else setTimeout(route, 0); }
  return { render: render };
})();

/* ===== Motion layer «رحلة الباحث»: motion that explains the content rather than decorating it.
   - Hero: a live questionnaire — Likert answers get picked one by one while n grows and the α gauge fills.
   - Headline: masked word reveal (SplitText mask). Page changes: violet wipe + clip-path reveals.
   - Journey: a researcher token travels the 12-rule path with scroll (MotionPath + ScrollTrigger scrub).
   - Cards: 3D tilt with a cursor-following light. Badges: a stamp, not confetti.
   Lenis smooth scroll synced with ScrollTrigger. Reduced-motion viewers get the static, complete UI. ===== */
const MOTION = (() => {
  const G = window.gsap, ST = window.ScrollTrigger, SPLIT = window.SplitText, MP = window.MotionPathPlugin;
  const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
  const on = !!G && !reduce;
  const hover = matchMedia('(hover: hover)').matches;
  let lenis = null;
  if (G) [ST, SPLIT, MP].forEach(p => p && G.registerPlugin(p));

  function initLenis() {
    if (!on || !window.Lenis || !ST) return;
    lenis = new window.Lenis({ duration: 1.1, smoothWheel: true, syncTouch: false });
    lenis.on('scroll', ST.update);
    G.ticker.add(t => lenis.raf(t * 1000));
    G.ticker.lagSmoothing(0);
    document.documentElement.classList.add('lenis');
  }
  const scrollTop = () => { if (lenis) lenis.scrollTo(0, { immediate: true }); else window.scrollTo(0, 0); };

  /* ---- page change: a violet wipe passes over, then blocks are revealed with a clip-path cut ---- */
  let wipe = null, first = true;
  function page(root) {
    if (!on || !root) return;
    ST && ST.getAll().forEach(t => t.vars && t.vars.id === 'rv' && t.kill());
    if (!first) {
      if (!wipe) { wipe = document.createElement('div'); wipe.className = 'wipe'; document.body.appendChild(wipe); }
      G.fromTo(wipe, { scaleX: 1, transformOrigin: '0% 50%' }, { scaleX: 0, duration: 0.55, ease: 'expo.inOut' });
    }
    first = false;
    const blocks = [...root.children].filter(n => n.nodeType === 1);
    const top = blocks.slice(0, 2), rest = blocks.slice(2);
    G.fromTo(top, { clipPath: 'inset(0 0 100% 0)', y: 24 }, { clipPath: 'inset(0 0 0% 0)', y: 0, duration: 0.8, ease: 'expo.out', stagger: 0.08, delay: 0.12, clearProps: 'clipPath,transform' });
    if (ST && rest.length) {
      G.set(rest, { clipPath: 'inset(0 0 100% 0)', y: 30 });
      ST.batch(rest, { id: 'rv', start: 'top 92%', once: true, onEnter: b => G.to(b, { clipPath: 'inset(0 0 0% 0)', y: 0, duration: 0.85, ease: 'expo.out', stagger: 0.08, clearProps: 'clipPath,transform' }) });
      ST.refresh();
    }
    root.querySelectorAll('.ring svg circle:last-child').forEach(c => { const d = c.getAttribute('stroke-dasharray'); if (!d) return; G.from(c, { strokeDashoffset: parseFloat(d), duration: 1.2, ease: 'power3.out', delay: 0.3 }); });
    tilt(root);
  }

  /* ---- +XP chip flies into the points pill ---- */
  function xp(n) {
    if (!n) return;
    const pill = document.querySelector('#pXp') || document.querySelector('.pill.xp');
    const f = document.createElement('div'); f.className = 'xpfloat'; f.textContent = '+' + n; document.body.appendChild(f);
    const r = pill ? pill.getBoundingClientRect() : { left: innerWidth / 2, top: 60, width: 0, height: 0 };
    if (!on) { f.style.left = (r.left + r.width / 2) + 'px'; f.style.top = (r.top + 6) + 'px'; setTimeout(() => f.remove(), 900); return; }
    G.set(f, { left: innerWidth / 2, top: innerHeight * 0.55, xPercent: -50, scale: 0.4, autoAlpha: 0 });
    G.timeline({ onComplete: () => f.remove() })
      .to(f, { scale: 1.25, autoAlpha: 1, duration: 0.3, ease: 'back.out(2.5)' })
      .to(f, { left: r.left + r.width / 2, top: r.top + r.height / 2, scale: 0.5, duration: 0.65, ease: 'power3.in' }, '+=0.15')
      .to(f, { autoAlpha: 0, duration: 0.12 })
      .add(() => pill && G.fromTo(pill, { scale: 1.3, backgroundColor: 'rgba(190,242,100,.6)' }, { scale: 1, backgroundColor: '', duration: 0.5, ease: 'elastic.out(1,.4)', clearProps: 'backgroundColor' }), '-=0.12');
  }

  /* ---- badge stamp: the medal slams in with a shock ring (used for badges and perfect scores) ---- */
  function stamp(label) {
    if (!on) return;
    const w = document.createElement('div'); w.className = 'stamp'; w.innerHTML = `<div class="ring"></div><div class="medal">${LOGO_SVG}<b>${esc(label || 'وسام جديد')}</b></div>`;
    document.body.appendChild(w);
    G.timeline({ onComplete: () => w.remove() })
      .fromTo(w.querySelector('.medal'), { scale: 2.6, rotate: -18, autoAlpha: 0 }, { scale: 1, rotate: 0, autoAlpha: 1, duration: 0.45, ease: 'power4.in' })
      .fromTo(w.querySelector('.ring'), { scale: 0.6, autoAlpha: 0.9 }, { scale: 2.4, autoAlpha: 0, duration: 0.7, ease: 'expo.out' }, '-=0.02')
      .to(w.querySelector('.medal'), { y: -6, duration: 0.1, yoyo: true, repeat: 1 }, '<')
      .to(w, { autoAlpha: 0, y: -20, duration: 0.4, delay: 1.1 });
  }
  const confetti = () => stamp('أداء كامل!');

  /* ---- hero: masked headline + live questionnaire filling the α gauge ---- */
  function hero(root) {
    if (!root) return;
    const live = root.querySelector('.hx-live');
    if (!on) return;
    const h = root.querySelector('.hx-title');
    if (h && SPLIT) { const sp = new SPLIT(h, { type: 'words', mask: 'words' }); G.from(sp.words, { yPercent: 120, duration: 0.9, ease: 'expo.out', stagger: 0.07 }); }
    G.from(root.querySelectorAll('.hx-kick, .hx-sub, .hx-cta > *, .hx-stats > div'), { y: 18, autoAlpha: 0, duration: 0.7, ease: 'expo.out', stagger: 0.05, delay: 0.3, clearProps: 'all' });
    root.querySelectorAll('[data-count]').forEach(n => { const v = +n.dataset.count || 0, o = { v: 0 }; G.to(o, { v, duration: 1.4, ease: 'expo.out', delay: 0.4, onUpdate: () => n.textContent = Math.round(o.v).toLocaleString('en-US') }); });
    const ticks = root.querySelectorAll('.ruler i'); if (ticks.length) G.from(ticks, { scaleY: 0, transformOrigin: '50% 100%', duration: 0.6, ease: 'expo.out', stagger: { each: 0.01, from: 'end' }, delay: 0.4 });
    const needle = root.querySelector('.ruler .needle'); if (needle) G.from(needle, { left: '100%', duration: 1.6, ease: 'expo.out', delay: 0.6 });
    if (live) questionnaire(live, root);
  }
  function questionnaire(live, root) {
    const rows = [...live.querySelectorAll('.lq-row')], arc = live.querySelector('.lq-arc'), val = live.querySelector('.lq-a'), nEl = live.querySelector('.lq-n'), L = arc ? arc.getTotalLength() : 0;
    const alphaAt = n => Math.min(0.91, 0.42 + 0.49 * (1 - Math.exp(-n / 9)));
    let n = 0; if (arc) G.set(arc, { strokeDasharray: L, strokeDashoffset: L });
    const round = () => {
      const tl = G.timeline({ onComplete: () => { if (document.body.contains(live)) round(); } });
      tl.to(live.querySelectorAll('.lq-dot.on'), { scale: 0.6, autoAlpha: 0.25, duration: 0.25, ease: 'power2.in', onComplete: () => live.querySelectorAll('.lq-dot.on').forEach(d => d.classList.remove('on')) });
      rows.forEach((r, i) => {
        const dots = r.querySelectorAll('.lq-dot'), pick = dots[Math.floor(Math.random() * 2) + (i % 2 ? 2 : 3)] || dots[3];
        tl.add(() => pick.classList.add('on')).fromTo(pick, { scale: 0.2 }, { scale: 1, autoAlpha: 1, duration: 0.45, ease: 'back.out(3)' }, '+=0.18');
      });
      tl.add(() => { n += 1; nEl.textContent = n; const a = alphaAt(n); G.to(arc, { strokeDashoffset: L * (1 - a), duration: 0.9, ease: 'expo.out' }); const o = { v: parseFloat(val.textContent) || 0 }; G.to(o, { v: a, duration: 0.9, ease: 'expo.out', onUpdate: () => val.textContent = o.v.toFixed(2) }); live.classList.toggle('good', a >= 0.7); }, '+=0.15')
        .to({}, { duration: 0.9 });
    };
    round();
  }

  /* ---- journey: path drawn to the current level, a token rides the path with scroll ---- */
  function journey(root) {
    if (!root) return;
    const p = root.querySelector('path.trail'), tok = root.querySelector('.jtoken');
    if (p) {
      const L = p.getTotalLength(), fr = Math.max(0, Math.min(1, +getComputedStyle(p).getPropertyValue('--lit') || 0));
      p.style.strokeDasharray = L; p.style.strokeDashoffset = L * (1 - fr);
      if (on && ST && fr > 0) G.fromTo(p, { strokeDashoffset: L }, { strokeDashoffset: L * (1 - fr), ease: 'none', scrollTrigger: { trigger: root, start: 'top 80%', end: 'center 50%', scrub: 0.8 } });
      if (tok && MP && on && ST) G.to(tok, { motionPath: { path: p, align: p, alignOrigin: [0.5, 0.5], start: 0, end: Math.max(0.0001, fr) }, ease: 'none', scrollTrigger: { trigger: root, start: 'top 80%', end: 'center 50%', scrub: 0.8 } });
      else if (tok) { const pt = p.getPointAtLength(L * fr); tok.setAttribute('transform', `translate(${pt.x} ${pt.y})`); }
    }
    if (!on || !ST) return;
    G.from(root.querySelectorAll('.jn'), { scale: 0, transformOrigin: '50% 50%', duration: 0.6, ease: 'back.out(2.4)', stagger: { each: 0.05, from: 'start' }, scrollTrigger: { trigger: root, start: 'top 80%', once: true } });
  }

  /* ---- 3D tilt with a cursor-following light ---- */
  function tilt(root) {
    if (!on || !hover) return;
    root.querySelectorAll('.linkcard, .hx-stats > div, .tile').forEach(el => {
      if (el.dataset.tilt) return; el.dataset.tilt = '1';
      const rx = G.quickTo(el, 'rotationX', { duration: 0.5, ease: 'power3' }), ry = G.quickTo(el, 'rotationY', { duration: 0.5, ease: 'power3' });
      G.set(el, { transformPerspective: 700 });
      el.addEventListener('pointermove', e => { const r = el.getBoundingClientRect(), px = (e.clientX - r.left) / r.width, py = (e.clientY - r.top) / r.height; rx((0.5 - py) * 8); ry((px - 0.5) * 10); el.style.setProperty('--mx', px * 100 + '%'); el.style.setProperty('--my', py * 100 + '%'); });
      el.addEventListener('pointerleave', () => { rx(0); ry(0); });
    });
  }
  function magnetic(root) {
    if (!on || !hover) return;
    root.querySelectorAll('.magnet').forEach(b => {
      const qx = G.quickTo(b, 'x', { duration: 0.4, ease: 'power3' }), qy = G.quickTo(b, 'y', { duration: 0.4, ease: 'power3' });
      b.addEventListener('pointermove', e => { const r = b.getBoundingClientRect(); qx((e.clientX - r.left - r.width / 2) * 0.22); qy((e.clientY - r.top - r.height / 2) * 0.3); });
      b.addEventListener('pointerleave', () => { qx(0); qy(0); });
    });
  }
  return { on, initLenis, scrollTop, page, xp, confetti, stamp, hero, journey, magnetic, tilt, get lenis() { return lenis; } };
})();

/* ===== AI assistant on the student's own account =====
   Inside claude.ai the page asks Claude through the viewer's own account (claude.use("sample")): nothing is billed to the site.
   On the public GitHub copy there is no in-page AI, so the panel builds the same request for the student to paste
   into their own Claude or ChatGPT account. */
const AI = (() => {
  let samp, ready = null;
  const get = () => {
    if (!ready) ready = (async () => { try { if (window.claude && typeof window.claude.use === 'function') return await window.claude.use('sample'); } catch (e) { } return null; })();
    return ready.then(s => (samp = s));
  };
  const ERR = { not_granted: 'لم تسمح للموقع باستخدام Claude من حسابك.', sampling_disabled: 'Claude غير متاح لحسابك الآن.', rate_limited: 'وصلت للحد المسموح؛ حاول بعد قليل.', session_expired: 'سجّل الدخول في Claude مرة أخرى.', refused: 'رفض Claude هذا الطلب؛ عدّل السؤال وحاول.', invalid_json: 'الرد لم يكن بالشكل المتوقع؛ حاول مرة أخرى.' };
  const SYS = 'أنت مساعد تعليمي في مادة «بناء المقاييس النفسية» لطالب في الدبلوم الخاص في علم النفس (جامعة المنيا). أجب بالعربية الفصحى المبسطة وباختصار منظم (نقاط قصيرة)، واذكر المعادلة عند الحاجة، واربط الإجابة بما يظهر في صفحة الطالب. لا تخترع مراجع أو أرقامًا غير موجودة. إذا كان السؤال خارج القياس النفسي والإحصاء فاعتذر بلطف.';
  const TASKS = [
    ['explain', 'اشرح لي هذه الصفحة', 'اشرح لي محتوى هذه الصفحة خطوة بخطوة كأنك تشرح لزميل، وركز على المفاهيم الأساسية.'],
    ['result', 'فسّر نتائجي', 'فسّر الأرقام والنتائج الظاهرة في هذه الصفحة (الثبات، الصدق، تحليل الفقرات...) وماذا أفعل بعدها.'],
    ['quiz', 'اختبرني', 'اكتب لي 3 أسئلة اختيار من متعدد (أ، ب، ج، د) على محتوى هذه الصفحة بنمط امتحانات الدبلوم، ثم الإجابات في النهاية.'],
    ['items', 'راجع عبارات مقياسي', 'راجع صياغة عبارات المقياس الظاهرة في الصفحة وفق قواعد الصياغة الجيدة (فكرة واحدة، بلا نفي مزدوج، لغة مناسبة)، واقترح بديلًا لكل عبارة ضعيفة.'],
    ['simple', 'بسّط لي مفهومًا', 'بسّط أهم مفهوم في هذه الصفحة بمثال من الحياة اليومية لطالب جامعي.']
  ];
  const ctxText = () => {
    const vc = document.querySelector('#vc'); let t = vc ? vc.innerText : '';
    t = t.replace(/\n{2,}/g, '\n').slice(0, 4500);
    return `الصفحة الحالية: ${(location.hash || '#home').slice(1).split('/')[0]}\n----- محتوى الصفحة -----\n${t}`;
  };
  const prompt = q => `${SYS}\n\n${ctxText()}\n\n----- سؤال الطالب -----\n${q}`;
  const fmtAns = s => esc(s).replace(/\*\*(.+?)\*\*/g, '<b>$1</b>').replace(/^#{1,4}\s*(.+)$/gm, '<b>$1</b>').replace(/\n/g, '<br>');

  function open(preset) {
    modal(`<div class="ai-head"><span class="ai-orb" aria-hidden="true"></span><div class="grow"><h3 style="margin:0">مساعد Claude</h3><p class="small muted" style="margin:0" id="aiMode">يتحقق من حسابك…</p></div><button class="ibtn" data-close aria-label="إغلاق">✕</button></div>
      <div class="seg ai-tasks" style="flex-wrap:wrap">${TASKS.map(t => `<button data-t="${t[0]}">${t[1]}</button>`).join('')}</div>
      <label class="fl" for="aiQ">أو اكتب سؤالك</label><textarea id="aiQ" placeholder="مثلًا: متى أستخدم KR‑20 بدل ألفا؟"></textarea>
      <div class="row"><button class="btn" id="aiGo">${ic('bolt')} اسأل</button><span class="small muted grow">الإجابة تعتمد على محتوى الصفحة المفتوحة.</span></div>
      <div id="aiOut" class="ai-out" aria-live="polite"></div>`, (bg, close) => {
      bg.classList.add('ai-sheet');
      const out = $('#aiOut', bg), q = $('#aiQ', bg), mode = $('#aiMode', bg);
      $$('[data-t]', bg).forEach(b => b.onclick = () => { q.value = TASKS.find(t => t[0] === b.dataset.t)[2]; ask(); });
      if (preset) q.value = preset;
      get().then(s => { mode.textContent = s ? 'يعمل بحسابك في Claude (أول مرة سيطلب إذنك)' : 'انسخ الطلب وافتحه في حسابك في Claude أو ChatGPT'; $('#aiGo', bg).innerHTML = s ? `${ic('bolt')} اسأل` : `${ic('copy')} انسخ الطلب`; });
      let ctl = null;
      async function ask() {
        const text = q.value.trim(); if (!text) { q.focus(); return; }
        const p = prompt(text), s = await get();
        if (!s) {
          await copyText(p); out.innerHTML = `<div class="notice small"><b>نُسخ الطلب مع محتوى الصفحة</b>الصقه في حسابك الخاص، وسيجيبك مباشرة.</div><div class="row"><a class="btn sm" href="https://claude.ai/new" target="_blank" rel="noopener">افتح Claude</a><a class="btn sm ghost" href="https://chatgpt.com/" target="_blank" rel="noopener">افتح ChatGPT</a></div>`;
          return;
        }
        if (ctl) ctl.abort(); ctl = new AbortController();
        const btn = $('#aiGo', bg); btn.disabled = true; out.innerHTML = '<div class="ai-ans"><span class="muted">يفكر…</span></div>';
        try {
          const r = await s(p, { signal: ctl.signal, cache: false, onText: ({ text }) => { $('.ai-ans', out).innerHTML = fmtAns(text); } });
          $('.ai-ans', out).innerHTML = fmtAns(r.text) + (r.truncated ? '<p class="small muted">(الإجابة طويلة وتوقفت؛ اطلب «أكمل».)</p>' : '');
          addXP(2); quest && quest('ai', 1);
        } catch (e) { if (e && e.name === 'AbortError') return; out.innerHTML = `<div class="notice bad small">${esc(ERR[e && e.code] || 'تعذر الحصول على إجابة الآن.')}</div>`; }
        btn.disabled = false;
      }
      $('#aiGo', bg).onclick = ask;
      q.addEventListener('keydown', e => { if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) ask(); });
    });
  }

  /* structured generation for the builder: returns [{text, dim, neg}] or null */
  async function items(s, k) {
    const sm = await get(); if (!sm || !sm.json) return null;
    const dims = s.dims.map(d => d.name + (d.def ? ` (${d.def})` : ''));
    const r = await sm.json(`أنت خبير في بناء المقاييس النفسية باللغة العربية. اكتب ${k} عبارة لقياس «${s.construct || s.name}» لدى ${s.group || 'طلاب الجامعة'}.\nالتعريف النظري: ${s.defC || '—'}\nالأبعاد: ${dims.join('، ')} — وزّع العبارات عليها.\nقواعد: فصحى واضحة، فكرة واحدة لكل عبارة، عبارة تقريرية بضمير المتكلم، بلا نفي مزدوج، ونحو ربع العبارات سالبة.\nأعد JSON فقط: {"items":[{"text":"...","dim":"اسم البعد كما هو","neg":false}]}`, { modelTier: 'default' });
    return r && Array.isArray(r.items) ? r.items : null;
  }

  function fab() {
    if (document.getElementById('aiFab')) return;
    const b = document.createElement('button'); b.id = 'aiFab'; b.className = 'ai-fab'; b.setAttribute('aria-label', 'مساعد Claude');
    b.innerHTML = `<span class="ai-orb" aria-hidden="true"></span><span>اسأل Claude</span>`; b.onclick = () => open();
    document.body.appendChild(b);
  }
  return { open, items, fab, available: () => get().then(Boolean), ERR };
})();

/* ===== Optional accounts =====
   Web (GitHub Pages): Firebase Authentication (email + password) + Firestore, only when
   window.FIREBASE_CONFIG is set (see SETUP-FIREBASE.md). Passwords never touch this code's
   storage: Firebase Auth stores them hashed (scrypt) on Google's servers.
   Inside claude.ai: the viewer's Claude account; progress syncs through the artifact db when allowed. */
(function () {
  const APP = window.ACCT_APP || 'scales';             /* 'scales' | 'stats' */
  const cfg = window.FIREBASE_CONFIG || null;
  const OTHER = APP === 'scales' ? 'stats' : 'scales';
  const st = { mode: 'none', user: null, ready: false, fb: null, db: null, lbOptIn: true, verified: false, name: '', failures: 0, lockUntil: 0, claudeDb: null, claudeUid: null, claudeName: '' };
  const h = s => (typeof esc === 'function' ? esc(s) : String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])));
  const getS = () => window.ACCT_HOOKS.get(), setS = x => window.ACCT_HOOKS.set(x), changed = () => window.ACCT_HOOKS.changed();
  const week7 = s => { let t = 0; for (let i = 0; i < 7; i++) { const d = new Date(); d.setDate(d.getDate() - i); const k = d.getFullYear() + '-' + String(d.getMonth() + 1).padStart(2, '0') + '-' + String(d.getDate()).padStart(2, '0'); t += (s.hist && s.hist[k]) || 0; } return t; };

  /* ---------- loading the vendored SDK (same-origin files, pinned 10.12.2) ---------- */
  function loadScript(src) { return new Promise((res, rej) => { const s = document.createElement('script'); s.src = src; s.async = false; s.onload = res; s.onerror = () => rej(new Error('load ' + src)); document.head.appendChild(s); }); }

  /* ---------- email action links (verify email / reset password) handled inside the site ----------
     Works when the Firebase email templates use this site as the custom action URL, and also as the continue URL. */
  const acs = () => ({ url: location.origin + location.pathname + '#account' });
  async function sendVer(user) { try { await user.sendEmailVerification(acs()); } catch (e) { if (/continue-uri|unauthorized-domain/.test(e && e.code || '')) await user.sendEmailVerification(); else throw e; } }
  async function sendReset(auth, em) { try { await auth.sendPasswordResetEmail(em, acs()); } catch (e) { if (/continue-uri|unauthorized-domain/.test(e && e.code || '')) await auth.sendPasswordResetEmail(em); else throw e; } }
  function actionBox(html) { const d = document.createElement('div'); d.setAttribute('role', 'dialog'); d.style.cssText = 'position:fixed;inset:0;z-index:9999;display:grid;place-items:center;background:rgba(0,0,0,.45);padding:16px'; d.innerHTML = '<div style="background:var(--surface,#fff);color:var(--ink,#111);border-radius:14px;padding:18px;max-width:420px;width:100%;display:flex;flex-direction:column;gap:10px;font:inherit">' + html + '</div>'; document.body.appendChild(d); return d; }
  async function handleAction(auth) {
    const q = new URLSearchParams(location.search), mode = q.get('mode'), code = q.get('oobCode'); if (!mode || !code) return;
    const clean = () => { try { history.replaceState(null, '', location.pathname + '#account'); } catch (e) { } };
    const bad = e => { const c = (e && e.code) || ''; return /expired|invalid-action-code/.test(c) ? 'اللينك ده انتهى أو اتستخدم قبل كده. لو طلبت أكتر من رسالة، استخدم آخر رسالة بس، أو اطلب رسالة جديدة من «حسابي».' : 'حصل خطأ: ' + (c || (e && e.message) || ''); };
    try {
      if (mode === 'verifyEmail') { await auth.applyActionCode(code); clean(); await new Promise(r => { const un = auth.onAuthStateChanged(u => { un(); r(u); }); }); if (auth.currentUser) { await auth.currentUser.reload(); await auth.currentUser.getIdToken(true); st.verified = !!auth.currentUser.emailVerified; changed(); } const d = actionBox('<h3 style="margin:0">تم تأكيد بريدك ✓</h3><p style="margin:0">حسابك مفعّل دلوقتي. لو مش داخل، سجّل دخولك من «حسابي».</p><button style="padding:8px 14px;border-radius:9px;border:0;cursor:pointer">تمام</button>'); d.querySelector('button').onclick = () => d.remove(); }
      else if (mode === 'resetPassword') { const em = await auth.verifyPasswordResetCode(code); const d = actionBox('<h3 style="margin:0">كلمة مرور جديدة</h3><p style="margin:0" dir="ltr">' + String(em).replace(/[<>&"]/g, '') + '</p><label for="rpw">كلمة المرور الجديدة (8 أحرف على الأقل)</label><input id="rpw" type="password" autocomplete="new-password" style="padding:8px;border-radius:8px;border:1px solid #ccc"><div id="rpm" style="font-size:14px"></div><button style="padding:8px 14px;border-radius:9px;border:0;cursor:pointer">احفظ</button>'); d.querySelector('button').onclick = async () => { const pw = d.querySelector('#rpw').value; if (pw.length < 8) { d.querySelector('#rpm').textContent = 'كلمة المرور لازم تكون 8 أحرف على الأقل'; return; } try { await auth.confirmPasswordReset(code, pw); clean(); d.querySelector('div').innerHTML = '<h3 style="margin:0">اتغيرت كلمة المرور ✓</h3><p style="margin:0">سجّل دخولك دلوقتي بكلمة المرور الجديدة من «حسابي».</p>'; setTimeout(() => d.remove(), 4000); } catch (e) { d.querySelector('#rpm').textContent = bad(e); } }; }
      else if (mode === 'recoverEmail') { await auth.checkActionCode(code); await auth.applyActionCode(code); clean(); const d = actionBox('<h3 style="margin:0">رجع بريدك القديم ✓</h3><p style="margin:0">يُفضّل تغيّر كلمة المرور دلوقتي من «حسابي».</p><button style="padding:8px 14px;border-radius:9px;border:0;cursor:pointer">تمام</button>'); d.querySelector('button').onclick = () => d.remove(); }
    } catch (e) { clean(); const d = actionBox('<h3 style="margin:0">مقدرناش نكمّل</h3><p style="margin:0">' + bad(e) + '</p><button style="padding:8px 14px;border-radius:9px;border:0;cursor:pointer">تمام</button>'); d.querySelector('button').onclick = () => d.remove(); }
  }
  async function initFirebase() {
    const base = (window.FIREBASE_SDK_BASE || 'vendor/firebase-10.12.2/');
    await loadScript(base + 'firebase-app-compat.js'); await loadScript(base + 'firebase-auth-compat.js'); await loadScript(base + 'firebase-firestore-compat.js');
    const app = firebase.initializeApp(cfg);
    if (window.FIREBASE_APPCHECK_SITE_KEY) { try { await loadScript(base + 'firebase-app-check-compat.js'); firebase.appCheck().activate(new firebase.appCheck.ReCaptchaV3Provider(window.FIREBASE_APPCHECK_SITE_KEY), true); } catch (e) { } }
    st.fb = firebase; st.db = firebase.firestore();
    firebase.auth().useDeviceLanguage();
    handleAction(firebase.auth());
    firebase.auth().onAuthStateChanged(async u => {
      st.user = u; st.verified = !!(u && u.emailVerified); st.name = u ? (u.displayName || '') : '';
      st.ready = true;
      if (u) { await pullAndMerge(); hookCloud(); } else unhookCloud();
      changed();
    });
  }
  const userRef = () => st.db.collection('users').doc(st.user.uid);
  async function pullAndMerge() {
    try {
      const snap = await userRef().get(); const d = snap.exists ? snap.data() : {};
      if (d.prefs && typeof d.prefs.lb === 'boolean') st.lbOptIn = d.prefs.lb;
      const remote = d[APP] && d[APP].state ? safeParse(d[APP].state) : null, local = getS();
      if (remote && (remote.xp || 0) > (local.xp || 0)) { setS(Object.assign({}, local, remote, { fs: local.fs, theme: local.theme })); }
      await pushState(getS(), true);
      if (APP === 'scales' && typeof SCALES !== 'undefined') {
        const qs = await userRef().collection('scales').get();
        const remoteScales = qs.docs.map(x => safeParse(x.data().json)).filter(Boolean);
        SCALES.merge(remoteScales);
        SCALES.list().forEach(s => pushScale(s));
      }
    } catch (e) {  }
  }
  const safeParse = s => { try { return JSON.parse(s); } catch (e) { return null; } };
  let lastPush = 0;
  async function pushState(s, force) {
    if (!st.user) return; const now = Date.now(); if (!force && now - lastPush < 4000) return; lastPush = now;
    const json = JSON.stringify(s); if (json.length > 400000) return;
    const ts = st.fb.firestore.FieldValue.serverTimestamp();
    try {
      await userRef().set({ name: (st.name || '').slice(0, 40), [APP]: { state: json, xp: s.xp | 0 }, updated: ts }, { merge: true });
      if (st.lbOptIn && st.verified) {
        const lbRef = st.db.collection('leaderboard').doc(st.user.uid), cur = await lbRef.get(), c = cur.exists ? cur.data() : {};
        const xpMine = Math.max(0, Math.min(1e6, s.xp | 0)), wMine = Math.max(0, Math.min(1e6, week7(s) | 0));
        const xpOther = (c.xp && c.xp[OTHER]) | 0, wOther = (c.w && c.w[OTHER]) | 0;
        await lbRef.set({ name: (st.name || 'باحث').slice(0, 40), xp: { [APP]: xpMine, [OTHER]: xpOther }, w: { [APP]: wMine, [OTHER]: wOther }, total: xpMine + xpOther, week: wMine + wOther, updated: ts });
      }
    } catch (e) {  }
  }
  async function pushScale(s) { if (!st.user || APP !== 'scales') return; const json = JSON.stringify(s); if (json.length > 900000) return; try { await userRef().collection('scales').doc(s.id).set({ json, updated: st.fb.firestore.FieldValue.serverTimestamp() }); } catch (e) {  } }
  function hookCloud() { if (typeof Cloud !== 'undefined') Cloud.push = s => pushState(s); if (typeof SCALES !== 'undefined') { SCALES.cloud.push = s => pushScale(s); SCALES.cloud.remove = id => st.user && userRef().collection('scales').doc(id).delete().catch(() => { }); } }
  function unhookCloud() { if (typeof Cloud !== 'undefined') Cloud.push = null; if (typeof SCALES !== 'undefined') { SCALES.cloud.push = null; SCALES.cloud.remove = null; } }

  /* ---------- Claude (artifact) mode ---------- */
  async function initClaude() {
    try {
      const [db, user] = await Promise.all([window.claude.use('db'), window.claude.use('user')]);
      if (!user) return; const id = await user.id(); st.claudeUid = id; st.claudeName = (await user.name()) || '';
      st.mode = 'claude'; st.ready = true;
      if (db && id) {
        const ref = db.collection('data/users/' + id).doc('progress-v2');
        const snap = await ref.get();
        if (snap.exists) { const r = safeParse(snap.data().json); if (r && (r.xp || 0) > (getS().xp || 0)) setS(Object.assign({}, getS(), r)); }
        let writable = true;
        const push = async s => { if (!writable) return; try { await ref.set({ json: JSON.stringify(s), updated: Date.now() }); } catch (e) { writable = false; } };
        await push(getS());
        if (writable) { st.claudeDb = db; if (typeof Cloud !== 'undefined') Cloud.push = push; }
        if (writable && typeof SCALES !== 'undefined') {
          const col = db.collection('data/users/' + id);
          try { const qs = await col.get(); SCALES.merge(qs.docs.filter(d => d.id.startsWith('scale_')).map(d => safeParse(d.data().json)).filter(Boolean)); } catch (e) { }
          SCALES.cloud.push = s => col.doc('scale_' + s.id).set({ json: JSON.stringify(s), updated: Date.now() }).catch(() => { });
          SCALES.cloud.remove = sid => col.doc('scale_' + sid).delete().catch(() => { });
        }
      }
      changed();
    } catch (e) { }
  }

  /* ---------- live forms (Firestore) ---------- */
  const forms = {
    available: () => st.mode === 'firebase' && !!st.user && st.verified,
    async panel(s, el, api) {
      const meta = s.live || null;
      const link = fid => location.origin + location.pathname + '#form-' + fid;
      const draw = () => {
        el.innerHTML = `<section class="card stack"><h3>${ic('link')} النموذج الأونلاين المباشر</h3>${s.live ? `<div class="row"><span class="live-dot" aria-hidden="true"></span><b class="grow">${s.live.open ? 'النموذج مفتوح ويستقبل الإجابات' : 'النموذج مغلق'}</b><span class="chip blue" id="liveN">${(s.liveSeen || []).length} استجابة وصلت</span></div><div class="codebox">${h(link(s.live.fid))}</div><div class="row"><button class="btn sm" id="lCp">${ic('copy')} نسخ الرابط</button><a class="btn sm ghost" target="_blank" rel="noopener noreferrer" href="https://wa.me/?text=${encodeURIComponent('أرجو الإجابة عن مقياس «' + s.name + '»: ' + link(s.live.fid))}">${ic('whats')} واتساب</a><button class="btn sm ghost" id="lTog">${s.live.open ? 'إغلاق النموذج' : 'إعادة فتحه'}</button><button class="btn sm ghost" id="lUpd">تحديث العبارات في النموذج</button></div><p class="small muted" style="margin:0">لا يحتاج المستجيب إلى حساب. الإجابات تُحفظ في حسابك فقط، ولا يقرؤها غيرك.</p>` : `<p style="margin:0">انشر نموذجًا له رابط تشاركه مع المستجيبين؛ تصل كل إجابة إلى هذا المقياس فور إرسالها.</p><button class="btn" id="lPub" style="justify-self:start">${ic('link')} نشر النموذج</button>`}</section>`;
        const pub = $('#lPub', el); if (pub) pub.onclick = async () => { try { const def = api.formDef(s); const fid = st.db.collection('forms').doc().id; await st.db.collection('forms').doc(fid).set({ owner: st.user.uid, title: s.name.slice(0, 120), def: JSON.stringify(def), open: true, created: st.fb.firestore.FieldValue.serverTimestamp() }); s.live = { fid, open: true }; s.liveSeen = []; api.save(); draw(); listen(); toast('نُشر النموذج'); } catch (e) { toast('تعذر النشر: ' + errMsg(e)); } };
        const cp = $('#lCp', el); if (cp) cp.onclick = () => copyText(link(s.live.fid));
        const tg = $('#lTog', el); if (tg) tg.onclick = async () => { try { await st.db.collection('forms').doc(s.live.fid).update({ open: !s.live.open }); s.live.open = !s.live.open; api.save(); draw(); } catch (e) { toast(errMsg(e)); } };
        const up = $('#lUpd', el); if (up) up.onclick = async () => { try { await st.db.collection('forms').doc(s.live.fid).update({ def: JSON.stringify(api.formDef(s)), title: s.name.slice(0, 120) }); api.save(); toast('حُدّثت عبارات النموذج'); } catch (e) { toast(errMsg(e)); } };
      };
      const listen = () => {
        if (!s.live) return; if (forms._un) forms._un();
        forms._un = st.db.collection('forms').doc(s.live.fid).collection('responses').onSnapshot(snap => {
          let added = 0; s.liveSeen = s.liveSeen || [];
          snap.docChanges().forEach(ch => { if (ch.type !== 'added') return; const id = ch.doc.id; if (s.liveSeen.includes(id)) return; const d = ch.doc.data(); if (!/^[0-7]{1,200}$/.test(d.r || '')) return; const res = api.addResponse(s, d.v, d.r.split('').map(Number), { src: 'live', g: String(d.g || '').slice(0, 10), a: String(d.a || '').slice(0, 3) }); if (res === 'ok') { s.liveSeen.push(id); added++; } });
          if (added) { api.save(); const n = $('#liveN'); if (n) n.textContent = s.liveSeen.length + ' استجابة وصلت'; toast('وصلت ' + added + ' استجابة جديدة'); }
        }, e => toast('انقطع الاستماع للنموذج: ' + errMsg(e)));
      };
      draw(); listen();
    },
    async respond(el, fid, api) {
      el.innerHTML = '<p class="muted">جارٍ تحميل النموذج…</p>';
      try {
        const doc = await st.db.collection('forms').doc(fid).get();
        if (!doc.exists || !doc.data().open) { el.innerHTML = '<div class="card"><h3>النموذج غير متاح</h3><p class="muted">ربما أغلقه الباحث أو أن الرابط غير صحيح.</p></div>'; return; }
        const def = safeParse(doc.data().def); if (!def || !Array.isArray(def.q) || !Array.isArray(def.l) || def.q.length > 200 || def.l.length > 7) throw new Error('bad');
        let sent = false; try { sent = localStorage.getItem('sent-' + fid) === '1'; } catch (e) { }
        if (sent) { el.innerHTML = '<div class="card"><h3>سبق أن أرسلت إجاباتك</h3><p class="muted">شكرًا لك؛ وصلت إجاباتك إلى الباحث.</p></div>'; return; }
        el.innerHTML = api.formHTML(def);
        api.bindForm(el, def, async (vals, demo) => {
          const btn = $('button[type=submit]', el); if (btn) btn.disabled = true;
          try { await st.db.collection('forms').doc(fid).collection('responses').add({ r: vals.join(''), v: String(def.v).slice(0, 20), g: String(demo.g || '').slice(0, 10), a: String(demo.a || '').slice(0, 3), ts: st.fb.firestore.FieldValue.serverTimestamp() }); try { localStorage.setItem('sent-' + fid, '1'); } catch (e) { }
            el.innerHTML = `<div class="card stack"><h3>${ic('check')} شكرًا لك</h3><p style="margin:0">وصلت إجاباتك إلى الباحث.</p></div>`; }
          catch (e) { if (btn) btn.disabled = false; toast('تعذر الإرسال: ' + errMsg(e)); }
        });
      } catch (e) { el.innerHTML = '<div class="card"><h3>تعذر تحميل النموذج</h3><p class="muted">تحقق من الرابط أو اتصال الإنترنت.</p></div>'; }
    }
  };

  /* ---------- errors (generic on sign-in to avoid account enumeration) ---------- */
  function errMsg(e) {
    const c = (e && e.code) || '';
    return ({
      'auth/email-already-in-use': 'هذا البريد مسجل بالفعل؛ جرّب تسجيل الدخول أو استعادة كلمة المرور.',
      'auth/invalid-email': 'صيغة البريد الإلكتروني غير صحيحة.',
      'auth/weak-password': 'كلمة المرور ضعيفة؛ استخدم 8 أحرف على الأقل مع حروف وأرقام.',
      'auth/invalid-credential': 'البريد أو كلمة المرور غير صحيحين.', 'auth/wrong-password': 'البريد أو كلمة المرور غير صحيحين.', 'auth/user-not-found': 'البريد أو كلمة المرور غير صحيحين.',
      'auth/too-many-requests': 'محاولات كثيرة؛ انتظر قليلًا ثم حاول مرة أخرى.',
      'auth/network-request-failed': 'تعذر الاتصال بالإنترنت.',
      'auth/requires-recent-login': 'لأسباب أمنية سجّل الخروج ثم الدخول مرة أخرى قبل هذا الإجراء.',
      'permission-denied': 'ليست لديك صلاحية لهذا الإجراء.'
    })[c] || 'حدث خطأ غير متوقع؛ حاول مرة أخرى.';
  }
  function pwStrength(p) { let s = 0; if (p.length >= 8) s++; if (p.length >= 12) s++; if (/[a-z]/i.test(p) && /\d/.test(p)) s++; if (/[^A-Za-z0-9]/.test(p)) s++; return s; }

  /* ---------- account page ---------- */

  /* ---------- real-name rule, Google sign-in, and verification by pasting the email link ---------- */
  const realName = nm => { const w = String(nm || '').trim().split(/\s+/); return nm.length <= 40 && w.length >= 2 && w.every(x => /^[A-Za-zء-ي٠-٩کگی'-]{2,}$/.test(x)) && !/\d/.test(nm); };
  const NAME_ERR = 'اكتب اسمك الحقيقي بالكامل (الاسم الأول واسم الأب على الأقل) بالحروف بس، من غير أرقام أو رموز.';
  function loadGsi() { return new Promise((res, rej) => { if (window.google && google.accounts && google.accounts.id) return res(); const s = document.createElement('script'); s.src = 'https://accounts.google.com/gsi/client'; s.async = true; s.onload = res; s.onerror = () => rej(new Error('gsi')); document.head.appendChild(s); }); }
  async function googleButton(box, onErr) {
    if (!box || st.mode !== 'firebase') return;
    const cid = window.GOOGLE_CLIENT_ID;
    if (cid) { try { await loadGsi(); google.accounts.id.initialize({ client_id: cid, callback: async r => { try { await st.fb.auth().signInWithCredential(firebase.auth.GoogleAuthProvider.credential(r.credential)); toast('أهلًا بيك'); } catch (e) { onErr(errMsg(e)); } } }); box.innerHTML = ''; google.accounts.id.renderButton(box, { theme: 'outline', size: 'large', text: 'continue_with', locale: 'ar', width: 280 }); return; } catch (e) { } }
    box.innerHTML = '<button class="btn ghost" type="button" style="width:100%"><svg viewBox="0 0 24 24" aria-hidden="true"><path fill="#4285F4" d="M22.5 12.3c0-.8-.1-1.5-.2-2.3H12v4.3h5.9a5 5 0 0 1-2.2 3.3v2.7h3.5c2.1-1.9 3.3-4.7 3.3-8z"/><path fill="#34A853" d="M12 23c3 0 5.5-1 7.3-2.7l-3.5-2.7c-1 .7-2.3 1.1-3.8 1.1-2.9 0-5.4-2-6.3-4.6H2.1v2.8A11 11 0 0 0 12 23z"/><path fill="#FBBC05" d="M5.7 14.1a6.6 6.6 0 0 1 0-4.2V7.1H2.1a11 11 0 0 0 0 9.8z"/><path fill="#EA4335" d="M12 5.4c1.6 0 3.1.6 4.2 1.7l3.1-3.1A11 11 0 0 0 2.1 7.1l3.6 2.8C6.6 7.3 9.1 5.4 12 5.4z"/></svg> متابعة بحساب Google</button>';
    box.querySelector('button').onclick = async () => { try { await st.fb.auth().signInWithPopup(new firebase.auth.GoogleAuthProvider()); toast('أهلًا بيك'); } catch (e) { const c = (e && e.code) || ''; onErr(/operation-not-allowed/.test(c) ? 'الدخول بحساب Google لسه مش مفعّل في إعدادات Firebase.' : /popup-closed|cancelled-popup/.test(c) ? 'اتقفلت نافذة Google قبل ما تكمّل.' : /network|internal/.test(c) ? 'نافذة Google مقدرتش تفتح على الشبكة دي. جرّب شبكة تانية، أو سجّل بالإيميل.' : errMsg(e)); } };
  }
  async function verifyFromLink(link) { const m = String(link).match(/oobCode(?:=|%3D)([A-Za-z0-9_-]+)/); if (!m) throw new Error('nolink'); await st.fb.auth().applyActionCode(m[1]); const u = st.fb.auth().currentUser; if (u) { await u.reload(); await u.getIdToken(true); st.verified = !!u.emailVerified; } }
  function renderPanel(el) {
    const s = getS(), rk = typeof rankOf === 'function' ? rankOf(s.xp).name : '';
    const perks = `<ul class="stack small" style="list-style:none;padding:0;margin:0;gap:6px">
      <li>${ic('check')} يُحفظ تقدمك (النقاط، والسلسلة، والدروس، و${APP === 'scales' ? 'مقاييسك واستجاباتها' : 'نتائج اختباراتك'}) على كل أجهزتك.</li>
      <li>${ic('check')} حساب واحد للموقعين: مختبر المقاييس ومختبر الإحصاء.</li>
      <li>${ic('check')} تظهر في لوحة المتصدرين الأسبوعية (يمكنك إخفاء اسمك).</li>
      ${APP === 'scales' ? `<li>${ic('check')} تنشر نموذج مقياسك أونلاين وتصلك الإجابات لحظيًا.</li>` : ''}
      <li>${ic('check')} وسام «عضو مسجل» و15 نقطة.</li></ul>`;
    if (st.mode === 'claude') {
      el.innerHTML = `<div class="section-title"><h2>حسابي</h2></div><section class="card stack"><div class="row"><span class="avatar${st.claudeName ? '' : ' guest'}">${typeof avatarInner === 'function' ? avatarInner(st.claudeName) : h((st.claudeName || '؟')[0])}</span><div class="grow"><b>${h(st.claudeName || 'حساب Claude')}</b><div class="small muted">${h(rk)}</div></div></div><p class="small" style="margin:0">داخل Claude يُستخدم حسابك في Claude مباشرة؛ لا تحتاج إلى تسجيل منفصل.</p><div class="notice small ${st.claudeDb ? 'good' : ''}">${st.claudeDb ? 'تقدمك ومقاييسك تُحفظ في حسابك وتتزامن بين أجهزتك.' : 'تقدمك يُحفظ في هذا المتصفح. لحفظه على كل أجهزتك استخدم نسخة الموقع المنشورة على GitHub وأنشئ حسابًا هناك.'}</div></section>`;
      return;
    }
    if (st.mode !== 'firebase') {
      el.innerHTML = `<div class="section-title"><h2>حسابي</h2></div><section class="card stack"><h3>${ic('user')} التسجيل غير مفعّل في هذه النسخة بعد</h3><p style="margin:0">تقدمك محفوظ الآن في هذا المتصفح (${s.xp} نقطة، سلسلة ${s.streak} يوم). يُفعَّل التسجيل بالبريد الإلكتروني بعد ربط الموقع بمشروع Firebase مجاني (ملف SETUP-FIREBASE.md).</p><div class="field"><label class="fl" for="nmLocal">اسمك في هذا المتصفح (اختياري)</label><input type="text" id="nmLocal" maxlength="40" value="${h(s.name || '')}"></div></section><section class="card">${perks}</section>`;
      $('#nmLocal', el).onchange = e => { const x = getS(); x.name = e.target.value.trim().slice(0, 40); setS(x); changed(); };
      return;
    }
    if (!st.ready) { el.innerHTML = '<p class="muted">جارٍ التحميل…</p>'; return; }
    if (!st.user) {
      let mode = 'up';
      const draw = () => {
        el.innerHTML = `<div class="section-title"><h2>حسابي</h2><span class="chip">اختياري</span></div>
        <div class="grid2"><section class="card stack"><div class="seg" role="tablist"><button data-am="up" aria-pressed="${mode === 'up'}">إنشاء حساب</button><button data-am="in" aria-pressed="${mode === 'in'}">تسجيل الدخول</button></div>
          <div id="aG"></div><div class="row small muted" style="justify-content:center">أو بالإيميل وكلمة مرور</div><form id="aForm" class="stack" novalidate>
            ${mode === 'up' ? `<div class="field"><label class="fl" for="aName">اسمك الحقيقي بالكامل (ثنائي أو ثلاثي)</label><input type="text" id="aName" maxlength="40" autocomplete="name" required placeholder="مثال: محمد أحمد عثمان"></div>` : ''}
            <div class="field"><label class="fl" for="aEmail">البريد الإلكتروني</label><input type="email" id="aEmail" autocomplete="email" required inputmode="email"></div>
            <div class="field"><label class="fl" for="aPw">كلمة المرور</label><input type="password" id="aPw" autocomplete="${mode === 'up' ? 'new-password' : 'current-password'}" minlength="8" required>${mode === 'up' ? '<div class="pbar" style="margin-top:6px"><i id="pwBar" style="width:0"></i></div><small class="muted" id="pwTxt">8 أحرف على الأقل، وحروف وأرقام.</small>' : ''}</div>
            ${mode === 'up' ? '<div class="field"><label class="fl" for="aPw2">تأكيد كلمة المرور</label><input type="password" id="aPw2" autocomplete="new-password" required></div>' : ''}
            <div id="aErr" class="notice bad small" hidden role="alert"></div>
            <button class="btn" type="submit">${mode === 'up' ? 'إنشاء الحساب' : 'دخول'}</button>
            ${mode === 'in' ? '<button class="btn ghost sm" type="button" id="aReset">نسيت كلمة المرور</button>' : ''}
          </form></section>
          <section class="card stack"><h3>${ic('star')} لماذا تنشئ حسابًا؟</h3><div class="tile hl"><small>تقدمك الحالي في هذا المتصفح</small><b style="direction:rtl;text-align:start;font-family:inherit">${s.xp} نقطة · ${s.streak} يوم</b></div>${perks}<p class="small muted" style="margin:0">الحساب اختياري تمامًا؛ كل أدوات الموقع تعمل بدونه.</p></section></div>`;
        $$('[data-am]', el).forEach(b => b.onclick = () => { mode = b.dataset.am; draw(); });
        const pw = $('#aPw', el), bar = $('#pwBar', el); if (bar) pw.oninput = () => { const k = pwStrength(pw.value); bar.style.width = k * 25 + '%'; bar.style.background = k >= 3 ? 'var(--good)' : k === 2 ? 'var(--amber)' : 'var(--bad)'; $('#pwTxt', el).textContent = ['ضعيفة جدًا', 'ضعيفة', 'متوسطة', 'جيدة', 'قوية'][k]; };
        const err = m => { const e = $('#aErr', el); e.hidden = false; e.textContent = m; };
        googleButton($('#aG', el), err);
        const rs = $('#aReset', el); if (rs) rs.onclick = async () => { const em = $('#aEmail', el).value.trim(); if (!/^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(em)) return err('اكتب بريدك الإلكتروني أولًا.'); try { await sendReset(st.fb.auth(), em); } catch (e) { } toast('إن كان البريد مسجلًا فستصلك رسالة لاستعادة كلمة المرور'); };
        $('#aForm', el).onsubmit = async e => {
          e.preventDefault(); if (Date.now() < st.lockUntil) return err('محاولات كثيرة؛ انتظر ' + Math.ceil((st.lockUntil - Date.now()) / 1000) + ' ثانية.');
          const em = $('#aEmail', el).value.trim(), p = $('#aPw', el).value;
          if (!/^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(em) || em.length > 254) return err('صيغة البريد الإلكتروني غير صحيحة.');
          const btn = $('button[type=submit]', el); btn.disabled = true;
          try {
            if (mode === 'up') {
              const nm = $('#aName', el).value.replace(/[<>]/g, '').trim().slice(0, 40);
              if (!realName(nm)) { btn.disabled = false; return err(NAME_ERR); }
              if (p.length < 8 || !(/[a-z؀-ۿ]/i.test(p) && /\d/.test(p))) { btn.disabled = false; return err('كلمة المرور: 8 أحرف على الأقل وتحتوي حروفًا وأرقامًا.'); }
              if (p !== $('#aPw2', el).value) { btn.disabled = false; return err('كلمتا المرور غير متطابقتين.'); }
              const cr = await st.fb.auth().createUserWithEmailAndPassword(em, p);
              await cr.user.updateProfile({ displayName: nm }); st.name = nm;
              try { await sendVer(cr.user); } catch (x) { }
              toast('أُنشئ حسابك. أرسلنا رسالة تأكيد لبريدك، ولو مش في الوارد بص في الـ Spam');
            } else { await st.fb.auth().signInWithEmailAndPassword(em, p); st.failures = 0; toast('مرحبًا بعودتك'); }
          } catch (x) { st.failures++; if (st.failures >= 5) { st.lockUntil = Date.now() + 30000 * Math.min(4, st.failures - 4); } btn.disabled = false; err(errMsg(x)); }
        };
      };
      draw(); return;
    }
    const u = st.user;
    el.innerHTML = `<div class="section-title"><h2>حسابي</h2></div>
      <div class="grid2"><section class="card stack"><div class="row"><span class="avatar">${typeof avatarInner === 'function' ? avatarInner(st.name || u.email) : h((st.name || u.email || '؟')[0])}</span><div class="grow"><b>${h(st.name || 'باحث')}</b><div class="small muted ltr" style="text-align:right">${h(u.email || '')}</div></div><span class="chip amber">${h(rk)}</span></div>
        ${st.verified ? '<div class="notice good small">بريدك مؤكد. تقدمك يُحفظ في حسابك تلقائيًا.</div>' : '<div class="notice warn small"><b>أكّد بريدك الإلكتروني</b>افتح الرسالة التي أرسلناها واضغط اللينك. <b>لو مش لاقيها، بص في مجلد الرسائل غير المرغوب فيها (Spam)</b> واضغط «ليست رسالة غير مرغوب فيها»، علشان الرسائل الجاية توصل الوارد. استخدم آخر رسالة بس. التأكيد مطلوب للظهور في لوحة المتصدرين ونشر النماذج.<div class="row" style="margin-top:6px"><button class="btn sm ghost" id="aVer">إعادة إرسال الرسالة</button><button class="btn sm ghost" id="aRel">تحققت، حدّث الحالة</button></div><details style="margin-top:8px"><summary><b>اللينك بيفتح صفحة خطأ؟</b></summary><p style="margin:6px 0">انسخ اللينك من الإيميل (اضغط عليه ضغطة طويلة ثم «نسخ الرابط»)، والصقه هنا، والموقع هيأكّد بريدك من غير ما يفتح صفحة Firebase.</p><div class="row"><input type="text" id="aLink" dir="ltr" class="grow" placeholder="https://...oobCode=..."><button class="btn sm" id="aLinkGo">أكّد</button></div></details></div>'}<div class="row" style="margin-top:8px"><a class="btn sm ghost" href="#admin" style="text-decoration:none">لوحة المسؤول</a><span class="muted small">للمسؤول فقط</span></div>
        ${realName(st.name || '') ? '' : '<div class="notice warn small"><b>اكتب اسمك الحقيقي بالكامل</b>مطلوب علشان مسؤول الموقع يقدر يتابعك ويساعدك.</div>'}<div class="field"><label class="fl" for="aNm">اسمك الحقيقي بالكامل</label><div class="row"><input type="text" id="aNm" maxlength="40" value="${h(st.name)}" class="grow"><button class="btn sm ghost" id="aNmS">حفظ</button></div></div>
        <label class="row small"><input type="checkbox" id="aLb" ${st.lbOptIn ? 'checked' : ''}> اعرض اسمي في لوحة المتصدرين</label>
        <div class="row"><button class="btn ghost sm" id="aOut">تسجيل الخروج</button><button class="btn sm danger" id="aDel">حذف الحساب</button></div></section>
        <section class="card stack"><h3>${ic('trophy')} لوحة المتصدرين</h3><div class="seg"><button data-lb="week" aria-pressed="true">آخر 7 أيام</button><button data-lb="total" aria-pressed="false">الكل</button></div><div id="lbList"><p class="muted small">جارٍ التحميل…</p></div></section></div>`;
    const v = $('#aVer', el); if (v) v.onclick = async () => { try { await sendVer(u); toast('أُرسلت رسالة تأكيد جديدة. استخدم آخر رسالة بس، لأن الرسائل القديمة بتبطل'); } catch (e) { toast(errMsg(e)); } };
    const lg = $('#aLinkGo', el); if (lg) lg.onclick = async () => { try { await verifyFromLink($('#aLink', el).value); if (st.verified) await pushState(getS(), true); toast(st.verified ? 'تم تأكيد بريدك ✓' : 'اتأكد، حدّث الصفحة'); changed(); } catch (e) { toast(e && e.message === 'nolink' ? 'ده مش لينك التأكيد. انسخ اللينك اللي في الإيميل كله.' : /expired|invalid-action-code/.test((e && e.code) || '') ? 'اللينك ده قديم أو اتستخدم قبل كده. اطلب رسالة جديدة واستخدم آخر واحدة.' : errMsg(e)); } };
    const r = $('#aRel', el); if (r) r.onclick = async () => { await u.reload(); st.verified = st.fb.auth().currentUser.emailVerified; if (st.verified) { await st.fb.auth().currentUser.getIdToken(true); await pushState(getS(), true); } changed(); };
    $('#aNmS', el).onclick = async () => { const nm = $('#aNm', el).value.replace(/[<>]/g, '').trim().slice(0, 40); if (!realName(nm)) return toast(NAME_ERR); try { await u.updateProfile({ displayName: nm }); st.name = nm; await pushState(getS(), true); toast('حُفظ الاسم'); changed(); } catch (e) { toast(errMsg(e)); } };
    $('#aLb', el).onchange = async e => { st.lbOptIn = e.target.checked; try { await userRef().set({ prefs: { lb: st.lbOptIn }, updated: st.fb.firestore.FieldValue.serverTimestamp() }, { merge: true }); if (!st.lbOptIn) await st.db.collection('leaderboard').doc(u.uid).delete(); else await pushState(getS(), true); toast(st.lbOptIn ? 'اسمك ظاهر في اللوحة' : 'أُخفي اسمك من اللوحة'); loadLb('week'); } catch (x) { toast(errMsg(x)); } };
    $('#aOut', el).onclick = () => st.fb.auth().signOut();
    $('#aDel', el).onclick = () => (typeof confirmBox === 'function' ? confirmBox : (m, y, f) => f())('سيُحذف حسابك وبياناتك المحفوظة على الخادم نهائيًا (تبقى نسخة على هذا الجهاز).', 'حذف الحساب', async () => {
      try { const uidv = u.uid; const sc = await userRef().collection('scales').get(); for (const d of sc.docs) await d.ref.delete(); const fq = await st.db.collection('forms').where('owner', '==', uidv).get(); for (const f of fq.docs) { const rs = await f.ref.collection('responses').get(); for (const r of rs.docs) await r.ref.delete(); await f.ref.delete(); } await st.db.collection('leaderboard').doc(uidv).delete().catch(() => { }); await userRef().delete(); await u.delete(); toast('حُذف الحساب'); } catch (e) { toast(errMsg(e)); }
    });
    const loadLb = async by => {
      const box = $('#lbList', el); if (!box) return;
      try {
        const qs = await st.db.collection('leaderboard').orderBy(by, 'desc').limit(20).get();
        const rows = qs.docs.map((d, i) => { const x = d.data(); const me = d.id === u.uid; return `<div class="lb-row ${me ? 'me' : ''}"><b class="num">${i + 1}</b><span>${h(x.name || 'باحث')}${me ? ' <span class="chip amber">أنت</span>' : ''}</span><span class="num">${(x[by] | 0)}</span></div>`; }).join('');
        box.innerHTML = rows || '<p class="muted small">لا أحد في اللوحة بعد. كن الأول!</p>';
        if (!st.verified) box.innerHTML += '<p class="small muted">أكّد بريدك لتظهر هنا.</p>';
      } catch (e) { box.innerHTML = `<p class="muted small">${h(errMsg(e))}</p>`; }
    };
    $$('[data-lb]', el).forEach(b => b.onclick = () => { $$('[data-lb]', el).forEach(x => x.setAttribute('aria-pressed', x === b)); loadLb(b.dataset.lb); });
    loadLb('week');
  }

  window.ACCT = {
    available: () => st.mode === 'firebase' || st.mode === 'claude',
    signedIn: () => (st.mode === 'firebase' && !!st.user) || (st.mode === 'claude' && !!st.claudeDb),
    displayName: () => st.mode === 'firebase' ? (st.user ? st.name : '') : st.mode === 'claude' ? st.claudeName : '',
    renderPanel, forms,
    // minimal Firestore access for the know-yourself journey and the admin page (firebase mode only)
    fs: {
      ready: () => st.mode === 'firebase' && !!st.user && !!st.db,
      uid: () => st.user ? st.user.uid : null, email: () => st.user ? (st.user.email || '') : '', name: () => st.name || '', verified: () => !!st.verified,
      get: async p => { const d = await st.db.doc(p).get(); return d.exists ? d.data() : null; },
      set: async (p, data) => st.db.doc(p).set(Object.assign({}, data, { updated: firebase.firestore.FieldValue.serverTimestamp() })),
      list: async c => (await st.db.collection(c).get()).docs.map(d => Object.assign({ id: d.id }, d.data())),
      del: async p => st.db.doc(p).delete(),
      refresh: async () => { if (!st.user) return false; try { await st.user.reload(); } catch (e) { } const u = st.fb.auth().currentUser; if (u) { await u.getIdToken(true); st.verified = !!u.emailVerified; } return st.verified; },
      put: async (p, data) => { const d = {}; Object.keys(data).forEach(k => { d[k] = data[k] === '__ts__' ? firebase.firestore.FieldValue.serverTimestamp() : data[k]; }); return st.db.doc(p).set(d); }
    },
    init() {
      if (window.claude && typeof window.claude.use === 'function') { initClaude(); return; }
      if (cfg && cfg.apiKey && cfg.projectId) { st.mode = 'firebase'; initFirebase().catch(e => { st.mode = 'none';  changed(); }); }
    }
  };
})();

/* ===== Boot ===== */
VIEWS.account = el => ACCT.renderPanel(el);
window.ACCT_HOOKS = {
  get: () => S,
  set: x => { S = Object.assign(fresh(), x); save(); },
  changed: () => { checkBadges(); renderStats(); if (current === 'account' || current === 'home') rerender(); }
};
applyPrefs(); tick(); shell();
if (typeof MOTION !== 'undefined') MOTION.initLenis();
if (typeof AI !== 'undefined') AI.fab();
window.addEventListener('hashchange', route);
ACCT.init();
route();
/* GoatCounter visit counter — only on the public GitHub Pages copy (dashboard: https://mo-psych-stats-lab.goatcounter.com) */
(function () {
  if (!/\.github\.io$/.test(location.hostname)) return;
  const path = () => { const h = location.hash || '#home'; return location.pathname + (h.startsWith('#take/') ? '#take' : h.startsWith('#form-') ? '#form' : h); };
  const hit = () => { try { window.goatcounter && window.goatcounter.count && window.goatcounter.count({ path: path(), title: document.title }); } catch (e) { } };
  window.goatcounter = { no_onload: true };
  const gs = document.createElement('script'); gs.async = true; gs.dataset.goatcounter = 'https://mo-psych-stats-lab.goatcounter.com/count'; gs.src = 'https://gc.zgo.at/count.js'; gs.referrerPolicy = 'no-referrer'; gs.onload = hit;
  document.head.appendChild(gs); window.addEventListener('hashchange', hit);
})();

#!/usr/bin/env python3
"""
Independent Mathematical Verification Script
Audit of H. P. Okolo's July 15, 2026 manuscript:
  "The Harmonic Divisor Fan and the Critical Exponent"
  Zenodo DOI: 10.5281/zenodo.21381824 (HDF_CRITICAL_EXPONENT.pdf)

Zero external dependencies (Python 3 standard library only).
Run: python3 verify.py
"""

import cmath
import math


# ---------------------------------------------------------------------------
# Helper: High-accuracy Hurwitz zeta zeta(s, a) via Euler-Maclaurin
# ---------------------------------------------------------------------------
# Even Bernoulli numbers B_{2k} for k = 1..8
_BERNOULLI_EVEN = [
    1.0 / 6.0,
    -1.0 / 30.0,
    1.0 / 42.0,
    -1.0 / 30.0,
    5.0 / 66.0,
    -691.0 / 2730.0,
    7.0 / 6.0,
    -3617.0 / 510.0,
]


def hurwitz_zeta(s: complex, a: float, n_terms: int | None = None) -> complex:
    """Compute Hurwitz zeta(s, a) for a > 0, s != 1 via Euler-Maclaurin."""
    if n_terms is None:
        n_terms = max(12, int(abs(s.imag) * 1.3) + 16)
    total = sum((n + a) ** (-s) for n in range(n_terms))
    x = n_terms + a
    # Integral + endpoint term
    total += (x ** (1.0 - s)) / (s - 1.0) + 0.5 * (x ** (-s))
    # Asymptotic Bernoulli tail
    poch = s
    x_pow = x ** (-s - 1.0)
    inv_x2 = 1.0 / (x * x)
    fact = 2.0
    for k, b2k in enumerate(_BERNOULLI_EVEN, start=1):
        total += (b2k / fact) * poch * x_pow
        poch *= (s + 2 * k - 1) * (s + 2 * k)
        fact *= (2 * k + 1) * (2 * k + 2)
        x_pow *= inv_x2
    return total


def cgamma(z: complex) -> complex:
    """Complex Gamma(z) via Lanczos approximation (g=7, n=9)."""
    p = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ]
    if z.real < 0.5:
        return math.pi / (cmath.sin(math.pi * z) * cgamma(1.0 - z))
    z -= 1.0
    x = p[0]
    for i in range(1, len(p)):
        x += p[i] / (z + i)
    t = z + 7.5
    return cmath.sqrt(2.0 * math.pi) * (t ** (z + 0.5)) * cmath.exp(-t) * x


def davenport_heilbronn(s: complex) -> complex:
    """
    Davenport-Heilbronn (1936) Dirichlet series f(s).
    Satisfies the EXACT functional equation:
        Lambda_DH(s) = (5/pi)^(s/2) Gamma((s+1)/2) f(s) = Lambda_DH(1 - s)
    yet has infinitely many zeros OFF the critical line Re(s) = 1/2.
    """
    kappa = (math.sqrt(10.0 - 2.0 * math.sqrt(5.0)) - 2.0) / (math.sqrt(5.0) - 1.0)
    z1 = hurwitz_zeta(s, 0.2)
    z2 = hurwitz_zeta(s, 0.4)
    z3 = hurwitz_zeta(s, 0.6)
    z4 = hurwitz_zeta(s, 0.8)
    l_chi = (5.0 ** (-s)) * (z1 + 1j * z2 - 1j * z3 - z4)
    l_chib = (5.0 ** (-s)) * (z1 - 1j * z2 + 1j * z3 - z4)
    return 0.5 * (1.0 - 1j * kappa) * l_chi + 0.5 * (1.0 + 1j * kappa) * l_chib


def lambda_dh(s: complex) -> complex:
    return ((5.0 / math.pi) ** (0.5 * s)) * cgamma(0.5 * (s + 1.0)) * davenport_heilbronn(s)


def secant_root(fn, z0: complex, z1: complex, iters: int = 30) -> complex:
    f0, f1 = fn(z0), fn(z1)
    for _ in range(iters):
        if abs(f1 - f0) < 1e-15:
            break
        z2 = z1 - f1 * (z1 - z0) / (f1 - f0)
        z0, f0 = z1, f1
        z1, f1 = z2, fn(z2)
        if abs(f1) < 1e-13:
            break
    return z1


def main() -> None:
    print("=" * 72)
    print("INDEPENDENT MATHEMATICAL CHECKS ON OKOLO (JULY 15, 2026) MANUSCRIPT")
    print("=" * 72)

    # 1. Part I, §1.4 (p. 5): |S|/k != log(d)/log(n)
    n, d = 30, 5
    frac_primes = 1.0 / 3.0
    log_ratio = math.log(d) / math.log(n)
    print("\n[1] Part I, Section 1.4 (p. 5): |S|/k vs log(d)/log(n)")
    print(f"    n = {n} = 2*3*5 (k=3), divisor d = {d}:")
    print(f"    |S|/k        = {frac_primes:.6f}")
    print(f"    log(d)/log(n)= {log_ratio:.6f}  (NOT equal; delta = {abs(log_ratio - frac_primes):.6f})")
    assert abs(log_ratio - frac_primes) > 0.1

    # 2. Part III, Lemma 3.4 (pp. 8-9): dropped log(x) factor
    x0 = math.e ** 2  # u = log(x0) = 2
    # For f(x) = x = e^u:
    #   H_def33 = -i(u d/du + 1/2) f = -i(2*e^2 + 0.5*e^2) = -2.5i * e^2
    #   H_lem34 = -i(x d/dx + 1/2) f = -i(e^2 + 0.5*e^2)   = -1.5i * e^2
    h_def33 = -1j * (math.log(x0) * x0 + 0.5 * x0)
    h_lem34 = -1j * (x0 + 0.5 * x0)
    print("\n[2] Part III, Lemma 3.4 (pp. 8-9): dropped u = log(x) factor")
    print(f"    Test function f(x) = x at x = e^2 (u = 2):")
    print(f"    -i(u d/du + 1/2) f(x) = {h_def33.imag:.6f}i")
    print(f"    -i(x d/dx + 1/2) f(x) = {h_lem34.imag:.6f}i  (DIFFERENT operators!)")
    assert abs(h_def33 - h_lem34) > 7.0

    # 3. Part III, Theorem 3 (p. 9): H = -i(x d/dx + 1/2) fails self-adjointness on L^2(R+, dx/x)
    # Using u = log(x), du = dx/x, real unit-norm Schwartz test function g(u) = pi^(-1/4) exp(-u^2/2)
    steps = 40000
    a, b = -10.0, 10.0
    du = (b - a) / steps
    norm2 = 0.0
    exp_h_stated = 0j
    for j in range(steps):
        u = a + (j + 0.5) * du
        g = (math.pi ** -0.25) * math.exp(-0.5 * u * u)
        gp = -u * g
        h_g = -1j * (gp + 0.5 * g)  # -i(d/du + 1/2) g
        norm2 += (g * g) * du
        exp_h_stated += g * h_g * du
    print("\n[3] Part III, Theorem 3 (p. 9): self-adjointness failure on L^2(R+, dx/x)")
    print(f"    Real normalized g(u) = pi^(-1/4) exp(-u^2/2): ||g||^2 = {norm2:.12f}")
    print(f"    <g, H g>_{{dx/x}} = {exp_h_stated.real:+.12f} {exp_h_stated.imag:+.12f}i  (exact: -0.5i, NOT real!)")
    assert abs(exp_h_stated - (-0.5j)) < 1e-10

    # 4. Part III, §3.4 (p. 11): V(t) sign error and omitted prime powers
    t_val = 1.0
    s_conv = 2.0 + 1j * t_val

    def von_mangoldt(m: int) -> float:
        for p in range(2, int(m ** 0.5) + 1):
            if m % p == 0:
                q = m
                while q % p == 0:
                    q //= p
                return math.log(p) if q == 1 else 0.0
        return math.log(m)

    neg_zlog = sum(von_mangoldt(m) * (m ** (-s_conv)) for m in range(2, 10000))
    okolo_form = sum(
        math.log(p) * (p ** -2.0) * math.sin(t_val * math.log(p))
        for p in range(2, 10000)
        if all(p % d != 0 for d in range(2, int(p ** 0.5) + 1))
    )
    print("\n[4] Part III, Section 3.4 (p. 11): V(t) = Im[-zeta'/zeta] sign & prime-power check (at s = 2 + i)")
    print(f"    True Im[-zeta'/zeta(2 + i)]           = {neg_zlog.imag:+.6f}")
    print(f"    Okolo formula +sum_p (log p)p^-2 sin  = {okolo_form:+.6f}  (OPPOSITE sign + missing prime powers)")
    assert neg_zlog.imag < 0 < okolo_form

    # 5. Part IV, §4.2-4.3 & Theorem 5 (pp. 14-15): xi(-2) != 0 and simple pole at s = 0
    zeta_3 = hurwitz_zeta(3.0 + 0j, 1.0).real
    xi_3 = (math.pi ** -1.5) * math.gamma(1.5) * zeta_3
    exact_ratio = zeta_3 / (2.0 * math.pi)
    eps = 1e-5
    xi_near_m2 = (math.pi ** (-(-2.0 + eps) / 2.0)) * cgamma((-2.0 + eps) / 2.0) * hurwitz_zeta(-2.0 + eps, 1.0)
    xi_near_0 = (math.pi ** (-eps / 2.0)) * cgamma(eps / 2.0) * hurwitz_zeta(eps, 1.0)
    print("\n[5] Part IV, Theorem 5 (pp. 14-15): Completed zeta xi(s) = pi^(-s/2) Gamma(s/2) zeta(s)")
    print(f"    xi(-2 + 1e-5)        = {xi_near_m2.real:.10f}  (NONZERO — trivial zero cancels pole!)")
    print(f"    xi(3) = zeta(3)/(2pi)= {xi_3:.10f}  (exact: {exact_ratio:.10f})")
    print(f"    eps * xi(eps) at 0   = {(eps * xi_near_0).real:.10f}  (exact: -1.0 -> SIMPLE POLE at s=0, not cancelled!)")
    assert abs(xi_near_m2.real - exact_ratio) < 1e-4
    assert abs((eps * xi_near_0).real - (-1.0)) < 1e-4

    # 6. Part V, Theorem 5.3 (p. 19): negative epsilon relative error is bounded by 1
    eps_neg = -0.01
    err_1000 = abs(math.exp(eps_neg * 1000.0) - 1.0)
    print("\n[6] Part V, Theorem 5.3 (p. 19): negative-epsilon relative amplitude error |exp(eps*u) - 1|")
    print(f"    eps = {eps_neg}, u = 1000: |exp(eps*u) - 1| = {err_1000:.10f}  (approaches 1, BOUNDED, not unbounded)")
    assert err_1000 < 1.0

    # 7. Hard Analytic Counterexample: Davenport-Heilbronn functional equation + off-line zero
    s_test = 0.73 + 4.2j
    fe_diff = abs(lambda_dh(s_test) - lambda_dh(1.0 - s_test))
    rho = secant_root(davenport_heilbronn, 0.8085 + 85.6993j, 0.8086 + 85.6994j)
    print("\n[7] Hard Counterexample: Davenport-Heilbronn function f(s)")
    print(f"    Functional equation |Lambda_DH(s) - Lambda_DH(1-s)| at s={s_test}: {fe_diff:.3e}")
    print(f"    Off-line zero found:   rho     = {rho.real:.10f} + {rho.imag:.10f}i")
    print(f"    Symmetric partner:     1 - rho = {(1-rho).real:.10f} {(1-rho).imag:+.10f}i  (conjugate: {(1-rho).real:.10f} + {rho.imag:.10f}i)")
    print(f"    |f(rho)| = {abs(davenport_heilbronn(rho)):.3e},  |f(1 - rho)| = {abs(davenport_heilbronn(1.0 - rho)):.3e}")
    assert fe_diff < 1e-12
    assert abs(davenport_heilbronn(rho)) < 1e-11
    assert abs(davenport_heilbronn(1.0 - rho)) < 1e-11
    assert abs(rho.real - 0.5) > 0.3

    print("\nALL 7 INDEPENDENT MATHEMATICAL CHECKS PASSED.")


if __name__ == "__main__":
    main()

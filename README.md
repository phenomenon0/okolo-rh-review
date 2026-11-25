# Independent Mathematical Audit of Okolo’s Riemann Hypothesis Proof

Live report: **https://phenomenon0.github.io/okolo-rh-review/**

This repository hosts a standalone, primary-source mathematical audit of Hanyelichukwu Paul Okolo’s July 15, 2026 Zenodo manuscript (*The Harmonic Divisor Fan and the Critical Exponent*, Zenodo DOI [10.5281/zenodo.21381824](https://doi.org/10.5281/zenodo.21381824), decrypted using the author’s published key `99hUeDd6ZNxkHMS`), together with a comparison against:

- Levent Alpöge’s November 20, 2026 Lean-verified proof of the Riemann Hypothesis ([alpo.ge/remark.pdf](https://alpo.ge/remark.pdf); [Lean FRO / Mathlib certification](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/)).
- A widely circulated ChatGPT review PDF (`Okolo_Riemann_Proof_Review.pdf`).

## Reproduce the Numerical Checks

`verify.py` uses only the Python 3 standard library (no third-party packages required):

```bash
python3 verify.py
```

It independently verifies:
1. The $|S|/k \neq \log d / \log n$ mismatch in Part I, §1.4.
2. The dropped $\log x$ factor in Part III, Lemma 3.4.
3. The non-real expectation $\langle g, H g\rangle_{dx/x} = -0.5i$ refuting self-adjointness in Part III, Theorem 3.
4. The sign and prime-power errors in $V(t)$ in Part III, §3.4.
5. Non-vanishing completed zeta $\xi(-2) = \xi(3) = \zeta(3)/(2\pi) \approx 0.19131 > 0$ and the uncancelled simple pole at $s = 0$ (residue $-1$) in Part IV, Theorem 5.
6. Boundedness of the negative-$\varepsilon$ relative error $\lvert e^{\varepsilon u}-1\rvert \to 1$ in Part V, Theorem 5.3.
7. An explicit off-line zero ($\rho \approx 0.808517 + 85.699348i$) of the Davenport–Heilbronn function, which satisfies the exact functional equation $\Lambda_{\text{DH}}(s) = \Lambda_{\text{DH}}(1-s)$ yet violates the Riemann Hypothesis.

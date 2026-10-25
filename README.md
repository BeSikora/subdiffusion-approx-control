<!--
TODO before publishing the repository:
  1. Adjust the script names in "Repository layout" and "Reproducing the
     results" to the actual files (the names below are a proposed layout).
  2. Fill in the exact package versions in requirements.txt.
  3. Fill in the measured run times in the table below.
  4. Add the DOI after acceptance.
-->

# Approximate controllability of semilinear subdiffusion equations — numerical experiments

This repository contains the code for the numerical experiments (Section 6) of

> B. Sikora, *Approximate controllability of semilinear subdiffusion equations
> with the Caputo time derivative*, submitted to Fractional Calculus and
> Applied Analysis (2026).

The code computes the approximate controls
$u_\varepsilon = L_T^\ast R(\varepsilon)\,p(y_\varepsilon)$
constructed in the proofs of Theorems 4.2 and 5.6 for the one-dimensional
problem

$$
\partial_t^\alpha y - \partial_{xx} y = f(y) + \mathbf 1_\omega u
\quad\text{in } (0,1)\times(0,T),\qquad y(0,t)=y(1,t)=0,\qquad y(\cdot,0)=y_0,
$$

and reproduces all tables and figures of Section 6.

## Method

* **Space:** spectral Galerkin method with the first $M$ eigenfunctions
  $e_n(x)=\sqrt2\sin(n\pi x)$, $\lambda_n=n^2\pi^2$; the nonlinearity is
  evaluated on a uniform grid of 1000 points and projected back onto the
  modes by the midpoint rule.
* **Time:** L1 scheme for the Caputo derivative on a uniform grid with $N_t$
  steps, implicit in the linear part and explicit in the nonlinearity. The
  discrete solution kernels are computed by the scheme itself, so no
  Mittag-Leffler functions are evaluated.
* **Control:** discretise-then-control. The operators $K$, $L_T$, $L_T^\ast$
  and the Gramian $\Gamma_0^T$ are replaced by their exact discrete
  counterparts; $R(\varepsilon)$ is applied through the eigendecomposition
  of the discrete Gramian.
* **Fixed point:** damped iteration
  $u^{(k+1)}=(1-\vartheta)u^{(k)}+\vartheta L_T^\ast R(\varepsilon)p(y^{(k)})$,
  stopped when the maximal change of the modal coefficients is below
  $10^{-10}$.
* **Unbounded nonlinearity** $f(s)=as+\sin s$: the linear part is moved into
  the generator ($\lambda_n\mapsto\lambda_n-a$), and the fixed point
  iteration involves only $h(s)=\sin s$ (Lemma 5.4, Theorem 5.6).

## Data of the experiments

| Quantity | Value |
|---|---|
| Domain, control region, final time | $\Omega=(0,1)$, $\omega=(0.2,0.4)$, $T=1$ |
| Initial state, target | $y_0(x)=\sin(\pi x)$, $y_1(x)=e^{-30(x-0.7)^2}$ |
| Orders | $\alpha\in\{0.6,\,0.75,\,0.9\}$ |
| Regularisation | $\varepsilon=10^{-1},10^{-2},\dots,10^{-10}$ |
| Discretisation (default) | $M=60$ modes, $N_t=400$ time steps |
| Damping | $\vartheta=0.5$ for $f(s)=\sin s$; $\vartheta=0.2$ for $f(s)=as+\sin s$ |
| Potentials | $a\in\{0,\,3,\,12\}$; for $a=12$ the first mode is unstable, $\lambda_1^a=\pi^2-12\approx-2.13$ |

Since all orders are larger than $\tfrac12$, the control space is
$L^2([0,T],L^2(\omega))$ throughout.

## Installation

```bash
git clone https://github.com/BeSikora/subdiffusion-approx-control.git
cd subdiffusion-approx-control
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Repository layout

```
.
├── subdiffusion_control.py   # discretisation, Gramian, fixed point iteration
├── experiments_bounded.py    # f(s) = sin s: Table 1, Figures 1–3
├── experiments_unbounded.py  # f(s) = a s + sin s: Table 2, Figure 4
├── run_all.py                # runs everything
├── results/                  # generated tables (.csv, .tex)
├── figures/                  # generated figures (.pdf)
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

## Reproducing the results

All results are produced by

```bash
python run_all.py
```

which writes the tables to `results/` and the figures to `figures/`.
Individual items can be reproduced as follows (numbering as in the
journal version of the paper).

| Item in the paper | Content | Command | Run time* |
|---|---|---|---|
| Table 1 | error and cost for $\alpha=0.75$, $N_t\in\{400,800,1600\}$, $M\in\{60,120\}$ | `python experiments_bounded.py --table1` | TODO |
| Figure 1 | error and cost vs. $\varepsilon$, $\alpha\in\{0.6,0.75,0.9\}$ | `python experiments_bounded.py --fig1` | TODO |
| Figure 2 | terminal states, $\alpha=0.75$, $\varepsilon\in\{10^{-2},10^{-4},10^{-8}\}$ | `python experiments_bounded.py --fig2` | TODO |
| Figure 3 | $\lvert u_\varepsilon(0.3,t)\rvert$ vs. $T-t$, $\varepsilon=10^{-6}$ | `python experiments_bounded.py --fig3` | TODO |
| Table 2 | error and cost for $f(s)=as+\sin s$, $a\in\{0,3,12\}$ | `python experiments_unbounded.py --table2` | TODO |
| Figure 4 | error, cost and terminal states for $a\in\{0,3,12\}$ | `python experiments_unbounded.py --fig4` | TODO |
| Section 6.3, text | generator $A$ vs. $A_a$ for $u=0$, $a=12$, $N_t=100,\dots,1600$ | `python experiments_unbounded.py --generator-check` | TODO |

\* on TODO (processor, RAM).

### Checking the output

The following values from the paper ($\alpha=0.75$, $f(s)=\sin s$,
$N_t=400$, $M=60$) allow a quick check:

| $\varepsilon$ | $\lVert y_\varepsilon(T)-y_1\rVert$ | $\lVert u_\varepsilon\rVert_{L^2(\omega\times(0,T))}$ |
|---|---|---|
| $10^{-2}$ | 3.79e-01 | 1.54e+00 |
| $10^{-6}$ | 5.00e-02 | 2.83e+01 |
| $10^{-10}$ | 1.82e-02 | 9.88e+02 |

The fixed point iteration should converge in 27–41 iterations for
$f(s)=\sin s$ ($\vartheta=0.5$) and in 82–122 iterations for $f(s)=as+\sin s$
($\vartheta=0.2$), and the discrete identity
$y_\varepsilon(T)=y_1-\varepsilon R(\varepsilon)p(y_\varepsilon)$ should hold
up to $10^{-9}$ and $2\cdot10^{-7}$, respectively.

For $\varepsilon=10^{-10}$ the cost is not converged with respect to the time
grid (see Table 1), and small differences in the last digits may occur with
other versions of the linear algebra libraries.

## How to cite

If you use this code, please cite the article. GitHub shows the citation in
the sidebar ("Cite this repository") based on `CITATION.cff`.

```bibtex
@article{Sikora2026subdiffusion,
  author  = {Sikora, Beata},
  title   = {Approximate controllability of semilinear subdiffusion
             equations with the {C}aputo time derivative},
  journal = {Fractional Calculus and Applied Analysis},
  year    = {2026},
  note    = {Submitted}
}
```

## License

The code is released under the MIT License; see [`LICENSE`](LICENSE).

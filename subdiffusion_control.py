"""
Approximate controllability of the semilinear subdiffusion equation

    d_t^alpha y - y_xx = f(y) + 1_omega u   in (0,1) x (0,T),
    y(0,t) = y(1,t) = 0,  y(x,0) = y0(x),

(Caputo derivative, 1/2 < alpha < 1), discretised by
  * spectral Galerkin in space: e_n(x) = sqrt(2) sin(n pi x), lambda_n = (n pi)^2,
    or lambda_n = (n pi)^2 - a for f(s) = a s + h(s) with the linear part moved
    into the generator (Section 6.3 of the paper),
  * the L1 scheme in time on a uniform grid (implicit in the linear part,
    explicit in the nonlinearity).

The control is computed as in the proof of Theorem 4.2 (discrete version):
    u_eps = L_T^* R(eps) p(y_eps),   R(eps) = (eps I + Gamma)^{-1},
    p(y)  = y1 - S(T) y0 - K f(y),
with a damped fixed-point iteration in y. At the fixed point the discrete
identity  y_eps(T) = y1 - eps R(eps) p(y_eps)  holds exactly.
"""
import numpy as np
from math import gamma
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Discretisation
# ----------------------------------------------------------------------
class Problem:
    def __init__(self, alpha, T=1.0, Nt=400, M=60, Nx=1000,
                 omega=(0.2, 0.4), f=np.sin, a=0.0):
        self.alpha, self.T, self.Nt, self.M = alpha, T, Nt, M
        self.dt = T / Nt
        self.f = f
        # spatial quadrature (midpoint rule) and basis
        self.x = (np.arange(Nx) + 0.5) / Nx
        self.wx = 1.0 / Nx
        n = np.arange(1, M + 1)
        self.lam = (n * np.pi) ** 2 - a      # a != 0: generator A_a = -d_xx - a
        self.E = np.sqrt(2.0) * np.sin(np.outer(n, np.pi * self.x))  # M x Nx
        self.chi = ((self.x > omega[0]) & (self.x < omega[1])).astype(float)
        # Gram matrix of the basis restricted to omega: W_nm = int_omega e_n e_m
        self.W = (self.E * self.chi) @ self.E.T * self.wx
        # L1 weights
        k = np.arange(Nt + 1)
        self.b = (k + 1) ** (1 - alpha) - k ** (1 - alpha)
        self.c = self.dt ** (-alpha) / gamma(2 - alpha)
        # discrete kernels: s_j (free response) and g_m (impulse response)
        self.s = self._march(np.ones(M), np.zeros((Nt, M)))           # (Nt+1) x M
        imp = np.zeros((Nt, M)); imp[0] = 1.0
        self.g = self._march(np.zeros(M), imp)[1:]                     # g[m] = y_{m+1}
        # discrete Gramian  Gamma = L L^*  w.r.t. <u,v> = dt * sum_j int_omega u_j v_j
        G = self.g[::-1]                       # G[j-1] = g_{Nt-j}, j = 1..Nt
        self.Gram = (G.T @ G) * self.W / self.dt
        self.Gram = 0.5 * (self.Gram + self.Gram.T)
        self.mu, self.V = np.linalg.eigh(self.Gram)
        self.mu = np.clip(self.mu, 0.0, None)

    # projections between physical and modal coordinates
    def to_modes(self, v):          # v: (..., Nx)
        return v @ self.E.T * self.wx

    def to_phys(self, a):           # a: (..., M)
        return a @ self.E

    def _march(self, a0, F):
        """L1 scheme for d^alpha a + lam a = F (modal, linear). F: Nt x M."""
        Nt, c, b = self.Nt, self.c, self.b
        a = np.zeros((Nt + 1, self.M)); a[0] = a0
        d = np.zeros((Nt, self.M))     # increments a_j - a_{j-1}
        for j in range(1, Nt + 1):
            hist = b[1:j][:, None] * d[j - 2::-1][:j - 1] if j > 1 else 0.0
            hist = hist.sum(axis=0) if j > 1 else 0.0
            a[j] = (F[j - 1] + c * a[j - 1] - c * hist) / (c + self.lam)
            d[j - 1] = a[j] - a[j - 1]
        return a

    def forward(self, a0, U):
        """Semilinear forward solve. U: Nt x Nx control on the grid (t_1..t_Nt).
        Returns modal trajectory (Nt+1) x M and the modal nonlinear forcing."""
        Nt, c, b = self.Nt, self.c, self.b
        BU = self.to_modes(U * self.chi)                    # Nt x M
        a = np.zeros((Nt + 1, self.M)); a[0] = a0
        d = np.zeros((Nt, self.M))
        Fnl = np.zeros((Nt, self.M))
        for j in range(1, Nt + 1):
            Fnl[j - 1] = self.to_modes(self.f(self.to_phys(a[j - 1])))
            if j > 1:
                hist = (b[1:j][:, None] * d[j - 2::-1][:j - 1]).sum(axis=0)
            else:
                hist = 0.0
            a[j] = (Fnl[j - 1] + BU[j - 1] + c * a[j - 1] - c * hist) / (c + self.lam)
            d[j - 1] = a[j] - a[j - 1]
        return a, Fnl

    def K(self, F):
        """K F = sum_j g_{Nt-j} F_j  (modal)."""
        return (self.g[::-1] * F).sum(axis=0)

    def Lstar(self, phi):
        """(L^* phi)_j = (1/dt) 1_omega sum_n g^n_{Nt-j} phi_n e_n(x)."""
        G = self.g[::-1]                                     # Nt x M
        return self.to_phys(G * phi[None, :]) * self.chi / self.dt

    def R(self, eps, v):
        return self.V @ ((self.V.T @ v) / (eps + self.mu))

    def control(self, a0, a1, eps, theta=0.5, tol=1e-10, maxit=500):
        """Damped fixed point iteration for y = Phi_eps(y)."""
        U = np.zeros((self.Nt, self.x.size))
        a, Fnl = self.forward(a0, U)
        for it in range(maxit):
            p = a1 - self.s[-1] * a0 - self.K(Fnl)
            U_new = self.Lstar(self.R(eps, p))
            U = (1 - theta) * U + theta * U_new
            a_new, Fnl = self.forward(a0, U)
            delta = np.max(np.abs(a_new - a))
            a = a_new
            if delta < tol:
                break
        cost = np.sqrt(self.dt * (U ** 2).sum() * self.wx)
        err = np.linalg.norm(a[-1] - a1)
        # check of the discrete identity y(T) = y1 - eps R(eps) p
        p = a1 - self.s[-1] * a0 - self.K(Fnl)
        ident = np.linalg.norm(a[-1] - (a1 - eps * self.R(eps, p)))
        return dict(a=a, U=U, err=err, cost=cost, it=it + 1, ident=ident)


# ----------------------------------------------------------------------
# Experiments
# ----------------------------------------------------------------------
def main():
    T, Nt, M = 1.0, 400, 60
    omega = (0.2, 0.4)
    y0 = lambda x: np.sin(np.pi * x)
    y1 = lambda x: np.exp(-30 * (x - 0.7) ** 2)         # target concentrated outside omega
    alphas = [0.6, 0.75, 0.9]
    exps = np.arange(1, 11)
    epss = 10.0 ** -exps

    results = {}
    for alpha in alphas:
        P = Problem(alpha, T, Nt, M, omega=omega)
        a0, a1 = P.to_modes(y0(P.x)), P.to_modes(y1(P.x))
        res = []
        for eps in epss:
            r = P.control(a0, a1, eps)
            res.append(r)
            print(f"alpha={alpha:4.2f} eps={eps:7.1e}  |y(T)-y1|={r['err']:9.3e}"
                  f"  |u|={r['cost']:9.3e}  it={r['it']:3d}  identity={r['ident']:.1e}")
        results[alpha] = (P, a0, a1, res)

    # Figure 1: error and cost vs eps
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for alpha, (P, a0, a1, res) in results.items():
        ax[0].loglog(epss, [r["err"] for r in res], "o-", label=rf"$\alpha={alpha}$")
        ax[1].loglog(epss, [r["cost"] for r in res], "o-", label=rf"$\alpha={alpha}$")
    ax[0].set_xlabel(r"$\varepsilon$"); ax[0].set_ylabel(r"$\|y_\varepsilon(T)-y_1\|$")
    ax[1].set_xlabel(r"$\varepsilon$"); ax[1].set_ylabel(r"$\|u_\varepsilon\|_{L^2(\omega\times(0,T))}$")
    for a_ in ax:
        a_.grid(True, which="both", alpha=0.3); a_.legend(); a_.invert_xaxis()
    fig.tight_layout(); fig.savefig("fig_error_cost.pdf"); fig.savefig("fig_error_cost.png", dpi=150)

    # Figure 2: final states vs target (alpha = 0.75)
    P, a0, a1, res = results[0.75]
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.plot(P.x, P.to_phys(a1), "k--", lw=2, label=r"target $y_1$")
    for k in [1, 3, 7]:
        ax.plot(P.x, P.to_phys(res[k]["a"][-1]), label=rf"$y_\varepsilon(T)$, $\varepsilon=10^{{-{exps[k]}}}$")
    ax.plot(P.x, P.to_phys(P.s[-1] * a0 + P.K(P.forward(a0, 0 * res[0]['U'])[1])),
            ":", label=r"uncontrolled $y(T)$")
    ax.axvspan(*omega, color="grey", alpha=0.2, label=r"$\omega$")
    ax.set_xlabel("$x$"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("fig_final_states.pdf"); fig.savefig("fig_final_states.png", dpi=150)

    # Figure 3: time profile of the control at x ~ 0.3 (eps = 1e-6), log-log in T - t;
    # the control is concentrated near t = T (dips = sign changes)
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for alpha, (P, a0, a1, res) in results.items():
        r = res[5]
        ix = np.argmin(np.abs(P.x - 0.3))
        tau = T - np.linspace(P.dt, T, Nt)            # T - t_j
        sel = tau > 0
        ax.loglog(tau[sel], np.abs(r["U"][sel, ix]), label=rf"$\alpha={alpha}$")
    ax.set_xlabel("$T-t$"); ax.set_ylabel(r"$|u_\varepsilon(0.3,t)|$")
    ax.grid(True, which="both", alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig("fig_control.pdf"); fig.savefig("fig_control.png", dpi=150)

    # Table for LaTeX
    with open("table_results.tex", "w") as fh:
        fh.write("% generated by subdiffusion_control.py\n")
        for k, eps in enumerate(epss):
            if exps[k] % 2:
                continue
            row = [f"$10^{{-{exps[k]}}}$"]
            for alpha in alphas:
                r = results[alpha][3][k]
                row += [f"{r['err']:.2e}", f"{r['cost']:.2e}"]
            fh.write(" & ".join(row) + r" \\" + "\n")


if __name__ == "__main__":
    main()

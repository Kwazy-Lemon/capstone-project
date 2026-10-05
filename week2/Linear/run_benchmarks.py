from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix
from linear_solvers import dc_solve, ac_sweep, backward_euler
ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results_python"
OUT.mkdir(exist_ok=True)
def export_csv(name, fields, rows):
    with (OUT / name).open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerows(rows)


def rc_matrices(R, C):
    g = 1.0 / R
    G = csc_matrix([[g, -g, 1], [-g, g, 0], [1, 0, 0]], dtype=float)
    M = csc_matrix(([C], ([1], [1])), shape=(3, 3))
    return G, M


def run_dc():
    # Benchmark A: 5-V source, R1=1 kohm in->out, R2=2 kohm out->ground.
    R1, R2, Vs = 1e3, 2e3, 5.0
    g1, g2 = 1/R1, 1/R2
    G = csc_matrix([[g1, -g1, 1], [-g1, g1+g2, 0], [1, 0, 0]])
    b = np.array([0, 0, Vs])
    x = dc_solve(G, b)
    actual = [x[0], x[1], x[2], (x[0]-x[1])/R1, x[1]/R2]
    exact = [5, 10/3, -1/600, 1/600, 1/600]
    labels = ["V_in", "V_out", "I_V1_in_to_ground", "I_R1_in_to_out", "I_R2_out_to_ground"]
    rows = [["divider", lab, "V" if i < 2 else "A", a, e, abs(a-e)]
            for i, (lab, a, e) in enumerate(zip(labels, actual, exact))]
    np.testing.assert_allclose(actual, exact, rtol=1e-12, atol=1e-14)
    # Report KCL and voltage-constraint residuals separately (different units).
    residual_rows = [["divider", "KCL_in", "A", (G@x-b)[0]],
                     ["divider", "KCL_out", "A", (G@x-b)[1]],
                     ["divider", "V_source_constraint", "V", (G@x-b)[2]]]
    G2 = csc_matrix([[1e-3, 0, 1], [0, 0.5e-3, -1], [1, -1, 0]])
    b2 = np.array([0, 1e-3, 2])
    x2 = dc_solve(G2, b2)
    actual2 = [*x2, x2[0]/1e3, x2[1]/2e3]
    exact2 = [4/3, -2/3, -4/3000, 4/3000, -1/3000]
    labels2 = ["V_n1", "V_n2", "I_V1_n1_to_n2", "I_R1_n1_to_ground", "I_R2_n2_to_ground"]
    rows += [["floating_source", lab, "V" if i < 2 else "A", a, e, abs(a-e)]
             for i, (lab, a, e) in enumerate(zip(labels2, actual2, exact2))]
    np.testing.assert_allclose(actual2, exact2, rtol=1e-12, atol=1e-14)
    residual_rows += [["floating_source", "KCL_n1", "A", (G2@x2-b2)[0]],
                      ["floating_source", "KCL_n2", "A", (G2@x2-b2)[1]],
                      ["floating_source", "V_source_constraint", "V", (G2@x2-b2)[2]]]
    export_csv("dc_results.csv", ["circuit", "quantity", "unit", "python", "analytic", "abs_error"], rows)
    export_csv("dc_residuals.csv", ["circuit", "equation", "unit", "residual"], residual_rows)
    return {"dc_benchmark_cases": 2, "dc_analytic_checks": "passed"}

def run_ac():
    R, C = 1e3, 1e-6
    G, M = rc_matrices(R, C)
    fc = 1 / (2*np.pi*R*C)
    f = np.unique(np.r_[np.logspace(0, 5, 301), fc])
    X = ac_sweep(G, M, np.array([0, 0, 1.0]), f)
    H = X[:, 1] / X[:, 0]
    exact = 1 / (1 + 1j*2*np.pi*f*R*C)
    error = np.abs(H-exact)
    np.testing.assert_allclose(H, exact, rtol=1e-12, atol=1e-14)
    k = np.argmin(abs(f-fc))
    np.testing.assert_allclose(abs(H[k]), 1/np.sqrt(2), rtol=1e-12)
    export_csv("ac_results.csv", ["frequency_Hz", "H_real", "H_imag", "gain_dB", "phase_deg", "abs_complex_error"],
               zip(f, H.real, H.imag, 20*np.log10(abs(H)), np.angle(H, deg=True), error))
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for axis, y, yexact, ylabel in [
        (ax[0], 20*np.log10(abs(H)), 20*np.log10(abs(exact)), "Gain (dB)"),
        (ax[1], np.angle(H, deg=True), np.angle(exact, deg=True), "Phase (degrees)")]:
        axis.semilogx(f, yexact, color="#152c40", label="Analytical RC response")
        axis.semilogx(f[::12], y[::12], "o", color="#078584", markersize=4, label="Python MNA")
        axis.axvline(fc, color="#c67b2b", linestyle="--", linewidth=1)
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=.22)
    ax[0].legend()
    ax[0].set_title(f"RC low-pass | R = 1 kOhm, C = 1 uF | fc = {fc:.3f} Hz")
    ax[1].set_xlabel("Frequency (Hz)")
    fig.tight_layout()
    fig.savefig(OUT / "ac_frequency_response.png", dpi=160)
    plt.close(fig)
    return {"ac_max_complex_error": float(error.max()), "cutoff_Hz": fc,
            "cutoff_gain_dB": float(20*np.log10(abs(H[k]))), "cutoff_phase_deg": float(np.angle(H[k], deg=True))}

def run_transient():
    R, C, Vs = 1e3, 1e-6, 5.0
    tau, tstop = R*C, 5*R*C
    G, M = rc_matrices(R, C)
    # t=0+ after a 0->5-V step: capacitor voltage is still zero.
    # Algebraic rows require V_in=5 and I_V1=(V_out-V_in)/R=-5 mA.
    x0 = np.array([Vs, 0, -Vs/R])
    source = lambda t: np.array([0, 0, Vs])
    ratios = [0.2, 0.1, 0.05, 0.025, 0.0125]
    rows, previous_error = [], None
    fig, ax = plt.subplots(figsize=(8, 4.7))
    tfine = np.linspace(0, tstop, 1001)
    ax.plot(tfine*1e3, Vs*(1-np.exp(-tfine/tau)), color="#152c40", lw=2, label="Analytical")
    for ratio in ratios:
        dt = ratio*tau
        t, X = backward_euler(G, M, source, x0, dt, tstop)
        exact = Vs*(1-np.exp(-t/tau))
        err = X[:, 1] - exact
        emax, rmse = float(np.max(abs(err))), float(np.sqrt(np.mean(err**2)))
        # Independent exact recurrence for the backward-Euler RC discretization.
        discrete = Vs*(1 - (1+dt/tau)**(-np.arange(t.size)))
        np.testing.assert_allclose(X[:, 1], discrete, rtol=1e-11, atol=2e-13)
        np.testing.assert_allclose(X[:, 0], Vs, rtol=1e-12)
        np.testing.assert_allclose(X[:, 2], -(Vs-X[:, 1])/R, rtol=1e-11, atol=1e-14)
        order = np.nan if previous_error is None else np.log2(previous_error/emax)
        rows.append([dt, ratio, t.size-1, emax, rmse, order])
        previous_error = emax
        tag = str(round(dt*1e6, 6)).replace(".", "p")
        export_csv(f"tran_dt_{tag}us.csv", ["time_s", "V_in", "V_out", "I_V1_A", "analytic_V_out", "error_V"],
                   zip(t, X[:, 0], X[:, 1], X[:, 2], exact, err))
        if ratio in (0.2, 0.05, 0.0125):
            ax.plot(t*1e3, X[:, 1], "--", label=f"BE: dt = {dt*1e6:g} us")
    assert all(rows[i+1][3] < rows[i][3] for i in range(len(rows)-1)), "Error must decrease as dt decreases."
    assert 0.9 < rows[-1][-1] < 1.1, "The observed asymptotic order should be approximately one."
    export_csv("timestep_error.csv", ["dt_s", "dt_over_tau", "steps", "max_error_V", "rmse_V", "observed_order"], rows)
    ax.set(xlabel="Time (ms)", ylabel="Output voltage (V)", title="RC step response | Backward Euler convergence")
    ax.grid(alpha=.22)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "transient_comparison.png", dpi=160)
    plt.close(fig)
    return {"transient_checks": "passed", "finest_dt_s": rows[-1][0],
            "finest_max_error_V": rows[-1][3], "observed_order": rows[-1][-1]}

def check_inductor_sign():
    # One additional physical check: series RL, V_out across L.
    # Unknowns [V_in, V_out, I_V1, I_L]; V_out - L*dI_L/dt = 0.
    R, L, Vs = 100.0, 10e-3, 5.0
    g = 1/R
    G = csc_matrix([[g,-g,1,0],[-g,g,0,1],[1,0,0,0],[0,1,0,0]])
    M = csc_matrix(([-L], ([3], [3])), shape=(4,4))
    xdc = dc_solve(G, [0,0,Vs,0])
    np.testing.assert_allclose(xdc, [Vs,0,-Vs/R,Vs/R], atol=1e-13)
    f = np.logspace(1,5,21)
    X = ac_sweep(G,M,[0,0,1,0],f)
    zL = 1j*2*np.pi*f*L
    np.testing.assert_allclose(X[:,1],zL/(R+zL),rtol=1e-12,atol=1e-13)
    dt, tstop = (L/R)/20, 5*(L/R)
    t,X = backward_euler(G,M,lambda t: np.array([0,0,Vs,0]),[Vs,Vs,0,0],dt,tstop)
    exact_discrete = Vs/R*(1-(1+dt*R/L)**(-np.arange(t.size)))
    np.testing.assert_allclose(X[:,3],exact_discrete,rtol=1e-11,atol=1e-13)
    return {"inductor_dc_ac_be_checks":"passed"}

def main():
    metrics = {}
    for run in (run_dc, run_ac, run_transient, check_inductor_sign):
        metrics.update(run())
    metrics["matlab_comparison_status"] = "PENDING: run matlab_reference.m in MATLAB after this script."
    (OUT / "validation_summary.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    print(f"Evidence saved in {OUT}")

if __name__ == "__main__":
    main()

import warnings
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import MatrixRankWarning, splu, spsolve


def solve_checked(A, b):
    """Solve a square sparse system. Do not form an explicit inverse."""
    A = csc_matrix(A)
    b = np.asarray(b).reshape(-1)
    if A.shape[0] != A.shape[1] or A.shape[0] != b.size:
        raise ValueError("A must be square; b must match the number of unknowns.")
    if not np.all(np.isfinite(A.data)) or not np.all(np.isfinite(b)):
        raise ValueError("A and b must contain finite numbers.")
    with warnings.catch_warnings():
        warnings.simplefilter("error", MatrixRankWarning)
        x = spsolve(A, b)
    if not np.all(np.isfinite(x)):
        raise RuntimeError("Linear solve returned non-finite values.")
    return x


def dc_solve(G, b_dc):
    """DC: dx/dt = 0, so G*x_dc = b_dc."""
    return solve_checked(G, b_dc)


def ac_sweep(G, M, b_ac, frequencies_hz):
    """Linear phasor AC analysis with the exp(+j*omega*t) convention.

    b_ac is a complex phasor vector, not the DC source vector.
    Output shape: (number of frequencies, number of unknowns).
    """
    G, M = csc_matrix(G), csc_matrix(M)
    f = np.asarray(frequencies_hz, dtype=float)
    if f.ndim != 1 or f.size == 0 or np.any(f <= 0) or not np.all(np.isfinite(f)):
        raise ValueError("Frequencies must be a nonempty vector of positive Hz values.")
    if M.shape != G.shape:
        raise ValueError("G and M must have the same shape.")
    return np.vstack([
        solve_checked(G + 1j * 2 * np.pi * frequency * M, b_ac)
        for frequency in f
    ])


def backward_euler(G, M, source_at, x0, dt, tstop):
    """Fixed-step linear transient analysis, including the supplied t=0 state.

    (G + M/dt)*x[n] = b(t[n]) + (M/dt)*x[n-1].
    source_at(t) returns a length-N vector. x0 must be a consistent MNA state.
    The caller defines switching at t=0 and all initial capacitor/inductor states.
    """
    if not np.isfinite(dt) or not np.isfinite(tstop) or dt <= 0 or tstop <= 0:
        raise ValueError("dt and tstop must be finite and positive.")
    steps = int(round(tstop / dt))
    if steps < 1 or not np.isclose(steps * dt, tstop, rtol=1e-10, atol=0):
        raise ValueError("For this fixed-step solver, tstop must be an integer multiple of dt.")
    G, M = csc_matrix(G), csc_matrix(M)
    x0 = np.asarray(x0, dtype=float).reshape(-1)
    if G.shape != M.shape or G.shape != (x0.size, x0.size):
        raise ValueError("G, M and x0 dimensions do not match.")
    if not np.all(np.isfinite(x0)):
        raise ValueError("x0 must be finite.")
    t = np.arange(steps + 1) * dt
    X = np.empty((steps + 1, x0.size))
    X[0] = x0                 # t=0 is an initial condition, not a computed step.
    history_matrix = M / dt
    A = (G + history_matrix).tocsc()
    lu = splu(A)              # Constant linear circuit + fixed dt: factor only once.
    for n in range(1, steps + 1):
        b = np.asarray(source_at(t[n]), dtype=float).reshape(-1)
        if b.size != x0.size or not np.all(np.isfinite(b)):
            raise ValueError("source_at must return a finite length-N vector.")
        X[n] = lu.solve(b + history_matrix @ X[n - 1])
    if not np.all(np.isfinite(X)):
        raise RuntimeError("Transient solve returned non-finite values.")
    return t, X

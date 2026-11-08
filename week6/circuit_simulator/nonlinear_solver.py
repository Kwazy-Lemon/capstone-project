import numpy as np

from scipy.sparse.linalg import spsolve


class NonlinearSolver:

    def __init__(self, assembler):
        self.assembler = assembler



def get_residual(self, x, G, b):

    # Linear contributions
    F = np.asarray(G @ x).ravel() - b

    # Nonlinear contributions
    for element in self.assembler.elements:
        if hasattr(element, "nonlinear_contribution"):
            f_entries, _ = element.nonlinear_contribution(self.assembler.node_map, x)

            for row, value in f_entries:
                F[row] += value

    return F


def get_jacobian(self, x, G):


    #unfinished for work 5

    # Add linear and nonlinear Jacobian contributions
    J = G + J_nonlinear

    return J


def solve(self, max_iter=100, tol=1e-9,
          residual_tol=1e-9, x0=None):

    # Assemble the linear DC system
    G, M, b = self.assembler.assemble_system()

    # Initial guess
    if x0 is None:  # default of 0
        x = np.zeros(G.shape[0])
    else:
        x = np.array(x0, dtype=float).copy()

    # Newton iteration

    for iteration in range(max_iter):

        # Calculate F(x)
        F = self.calculate_residual(x, G, b)

        # Calculate J(x)
        J = self.calculate_jacobian(x, G)

        # Solve J * dx = -F
        dx = spsolve(J, -F)

        # Update solution
        x_new = x + dx

        # Check convergence
        change = np.linalg.norm(dx, np.inf)
        residual = np.linalg.norm(self.calculate_residual(x_new, G, b),np.inf)

        if change < tol and residual < residual_tol:
            #Newton converged
            return x_new

        x = x_new

    raise RuntimeError("Newton-Raphson failed to converge.")
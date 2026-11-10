import numpy as np

from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix


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

        j_entries = []

        # collect differential conductance, their row and columns as 3 tuples into j_entries
        for element in self.assembler.elements:

            #only take nonlinear components
            if hasattr(element, "nonlinear_contribution"):

                #only extract j array, ignore f array
                _, entries = element.nonlinear_contribution(self.assembler.node_map, x)

                #extend each stamping into j_entries
                j_entries.extend(entries)

        # if no diode/or diode connected to 2 grounds condition
        if not j_entries:
            return G

        rows, cols, values = zip(*j_entries)
        J = coo_matrix((values, (rows, cols)), shape=G.shape)

        return G + J




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
            F = self.get_residual(x, G, b)

            # Calculate J(x)
            J = self.get_jacobian(x, G)

            # Solve J * dx = -F
            dx = spsolve(J, -F)

            # Update solution
            x_new = x + dx

            # Check convergence
            change = np.linalg.norm(dx, np.inf)
            residual = np.linalg.norm(self.get_residual(x_new, G, b), np.inf)

            if change < tol and residual < residual_tol:
                #Newton converged
                return x_new

            x = x_new

        raise RuntimeError("Newton-Raphson failed to converge.")

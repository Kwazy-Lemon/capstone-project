import numpy as np

from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix

from linear_solvers import ac_sweep as linear_ac_sweep


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




    def solve_dc(self, max_iter=100, tol=1e-9,
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

    def ac_sweep(self, b_ac, frequencies_hz, x0=None, max_iter=100,tol=1e-9, residual_tol=1e-9):

        # Non-linear solve the DC operating points of the circuit
        x_dc = self.solve_dc(max_iter=max_iter, tol=tol, residual_tol=residual_tol, x0=x0)

        # Jacobian at dc operating point
        G, M, _ = self.assembler.assemble_system()
        J_dc = self.get_jacobian(x_dc, G)

        # using ac-sweep from linear_solver
        x_ac = linear_ac_sweep(J_dc, M, b_ac, frequencies_hz)
        return x_dc, x_ac

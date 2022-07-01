import os
import sys
import numpy as np

try:
    from utils import to_vtk
except:
    to_vtk = None


sys.path.append("../build")
import example2d as example


class Problem:
    def __init__(self):
        self.prefix = "results"
        if not os.path.isdir(self.prefix):
            os.mkdir(self.prefix)

        self.X, self.conn, self.bcs, self.forces = self.X_conn_bcs_forces_2d()
        return

    def set_up_model(self):
        """
        Set up finite element model
        """
        # Create model
        model = example.Elasticity_Model(self.X, self.bcs)

        # Create element and add to model
        elem = example.Elasticity_CPS4(self.conn)
        model.add_element(elem)

        # Create constitutive and add to model
        q = 5.0
        E = 70e3
        nu = 0.3
        density = 1.0
        design_stress = 1e3
        con = example.TopoIsoConstitutive_CPS4(elem, q, E, nu, density, design_stress)
        model.add_constitutive(con)

        # Initialize model
        model.init()

        return model

    def set_up_filter(self):
        """
        Set up finite element model
        """
        # Create model
        fltr = example.Helmholtz_Model(self.X)

        # Create element and add to model
        r0 = 0.05
        elem = example.Helmholtz_CPS4(self.conn, r0)
        fltr.add_element(elem)

        # Create constitutive and add to model
        con = example.HelmholtzConstitutive_CPS4(elem)
        fltr.add_constitutive(con)

        # Initialize model
        fltr.init()

        return fltr

    def solve_linear(self, model, fltr, forces):
        """
        Solve the linear problem Ku = f
        """
        # Set density to 1
        rho_ad = fltr.new_solution()
        rho_np = np.array(rho_ad, copy=False)
        rho_np[:] = 1.0

        # Allocate and compute the jacobian matrix
        K = model.new_matrix()
        model.set_design_vars(rho_ad)
        model.jacobian(K)

        # Set up rhs
        rhs_ad = model.new_solution()
        rhs_np = np.array(rhs_ad, copy=False)
        rhs_np[:] = forces[:]

        # Create the AMG solver
        amg = model.new_amg(3, 1.333, K, True)

        # Solve for u
        u_ad = model.new_solution()
        u_np = np.array(u_ad, copy=False)
        amg.cg(rhs_ad, u_ad, 5, 100)

        # res == Ku ?
        res_ad = model.new_solution()
        res_np = np.array(res_ad, copy=False)

        res_ad.zero()
        model.set_solution(u_ad)
        model.residual(res_ad)

        print(res_np)
        print(rhs_np)

        return u_np

    def solve_nonlinear(self, model, fltr, forces, maxit=20, rtol=1e-8, atol=1e-30):
        """
        Solve a series of linearized problem K(u)u = f and update the nonlinear
        solution by Newton's scheme
        """
        # Set density to 1.0
        rho_ad = fltr.new_solution()
        rho_np = np.array(rho_ad, copy=False)
        rho_np[:] = 1.0

        # Allocate and compute the jacobian matrix
        J = model.new_matrix()
        model.set_design_vars(rho_ad)

        # Set up residual
        res_ad = model.new_solution()
        res_np = np.array(res_ad, copy=False)
        # rhs_np[:] = forces[:]

        # Set up initial solution
        u_ad = model.new_solution()
        u_np = np.array(u_ad, copy=False)
        u_np[:] = 0.0

        # Set up update direction
        p_ad = model.new_solution()
        p_np = np.array(p_ad, copy=False)
        p_np[:] = 0.0

        # Solve for u
        sol = []
        for i in range(maxit):
            # set solution, update K and res
            model.set_solution(u_ad)
            model.jacobian(J)
            model.residual(res_ad)
            res_np[:] -= forces[:]
            res_np[:] *= -1.0

            # Set up solver and solve: K(u) p = res
            num_levels = 3
            omega = 4.0 / 3.0
            epsilon = 0.0
            amg = model.new_amg(num_levels, omega, epsilon, J, True)
            amg.cg(res_ad, p_ad, 5, 100)

            # update the solution: u = u + p
            u_np[:] += 1.0 * p_np
            sol.append(u_np)
            print("|u| = ", np.linalg.norm(u_np))
            print("|p| = ", np.linalg.norm(p_np))

            if to_vtk:
                to_vtk(
                    self.conn,
                    self.X,
                    {"ux": u_np[:, 0], "uy": u_np[:, 1]},
                    vtk_name=os.path.join(self.prefix, f"nonlinear_{i}.vtk"),
                )

            res_norm = np.linalg.norm(res_np)
            u_norm = np.linalg.norm(u_np)

            if res_norm < rtol * u_norm or res_norm < atol:
                break

        return sol

    def X_conn_bcs_forces_2d(self):
        nx = 96
        ny = 32

        lx = 3.0
        ly = 1.0

        nnodes = (nx + 1) * (ny + 1)
        nelems = nx * ny
        nbcs = ny + 1

        conn = np.zeros((nelems, 4), dtype=np.int32)
        X = np.zeros((nnodes, 2), dtype=np.double)
        bcs = np.zeros((nbcs, 2), dtype=np.int32)
        forces = np.zeros((nnodes, 2), dtype=np.double)
        nodes = np.zeros((nx + 1, ny + 1), dtype=int)

        # Set the node locations
        for j in range(ny + 1):
            for i in range(nx + 1):
                nodes[i, j] = i + (nx + 1) * j
                X[nodes[i, j], 0] = lx * i / nx
                X[nodes[i, j], 1] = ly * j / ny

        # Set the connectivity
        for j in range(ny):
            for i in range(nx):
                elem = i + nx * j
                conn[elem, 0] = nodes[i, j]
                conn[elem, 1] = nodes[i + 1, j]
                conn[elem, 2] = nodes[i + 1, j + 1]
                conn[elem, 3] = nodes[i, j + 1]

        # Set the boundary conditions, using the bits of the value to indicate
        # which values to fix
        bcs_val = 1 | 2
        index = 0
        for j in range(ny + 1):
            bcs[index, 0] = nodes[0, j]
            bcs[index, 1] = bcs_val
            index += 1

        for j in range(ny + 1):
            if j < 0.2 * ny:
                forces[nodes[-1, j], 1] = -100.0

        return X, conn, bcs, forces


if __name__ == "__main__":
    prob = Problem()
    X, conn, bcs, forces = prob.X_conn_bcs_forces_2d()
    model = prob.set_up_model()
    fltr = prob.set_up_filter()
    u = prob.solve_nonlinear(model, fltr, forces)

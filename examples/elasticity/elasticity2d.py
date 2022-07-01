import numpy as np
from mpi4py import MPI
from paropt import ParOpt

try:
    from utils import to_vtk
except:
    to_vtk = None

import os
import sys

sys.path.append("build")
import example2d as example


class TopOpt(ParOpt.Problem):
    def __init__(self, model, fltr, vol, vol_target, num_levels=3, nvars=1, ncon=1):
        super(TopOpt, self).__init__(MPI.COMM_SELF, nvars, ncon)

        self.model = model
        self.fltr = fltr
        self.vol = vol
        self.vol_target = vol_target
        self.num_levels = num_levels

        # Create the solution objects
        self.K = self.model.new_matrix()
        self.u_ref = self.model.new_solution()
        self.u = np.array(self.u_ref, copy=False)
        self.f_ref = self.model.new_solution()
        self.f = np.array(self.f_ref, copy=False)

        # Create the buffer objects for the solution
        self.x_ref = self.fltr.new_solution()
        self.x = np.array(self.x_ref, copy=False)
        self.rho_ref = self.fltr.new_solution()
        self.rho = np.array(self.rho_ref, copy=False)
        self.fltr_res_ref = self.fltr.new_solution()

        # Temporary vectors for derivatives
        self.dfdrho_ref = self.fltr.new_solution()
        self.dfdrho = np.array(self.dfdrho_ref, copy=False)
        self.dfdx_ref = self.fltr.new_solution()
        self.dfdx = np.array(self.dfdx_ref, copy=False)

        # Set up the Jacobian matrix for the filter problem
        self.Kf = self.fltr.new_matrix()
        self.fltr.jacobian(self.Kf)

        # Set up the AMG for the filter problem
        omega = 4.0 / 3.0
        print_info = True
        self.fltr_amg = self.fltr.new_amg(self.num_levels, omega, self.Kf, print_info)

        # Set the scaling for the compliance
        self.compliance_scale = None

        # Set the iteration count
        self.prefix = "results"
        if not os.path.isdir(self.prefix):
            os.mkdir(self.prefix)
        self.vtk_iter = 0

        return

    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        lb[:] = 1e-3
        ub[:] = 1.0
        x[:] = 0.5
        return

    def evalObjCon(self, x):
        """
        Return the objective, constraint and fail flag
        """

        # Set the design variable values
        self.x[:, 0] = x[:]
        self.fltr.set_design_vars(self.x_ref)

        # Compute the filter residual
        self.fltr.residual(self.fltr_res_ref)

        # Scale the residual by -1.0 and solve for rho
        monitor = 5
        max_iters = 100
        self.fltr_res_ref.scale(-1.0)
        self.fltr_amg.cg(self.fltr_res_ref, self.rho_ref, monitor, max_iters)

        # Set the design variable values
        self.model.set_design_vars(self.rho_ref)

        # Set up AMG for the structural problem
        self.model.jacobian(self.K)
        omega = 4.0 / 3.0
        print_info = True
        amg = self.model.new_amg(self.num_levels, omega, self.K, print_info)
        amg.cg(self.f_ref, self.u_ref, monitor, max_iters)

        # Set the new solution
        self.model.set_solution(self.u_ref)

        # Compute the objective function
        obj = np.sum(self.u * self.f)
        if self.compliance_scale is None:
            self.compliance_scale = 1.0 / obj
        obj *= self.compliance_scale

        # Compute the volume
        volume = self.vol.eval_functional()
        print("Volume = ", volume)
        con = np.array([1.0 - volume / self.vol_target])

        # Export to vtk every 10 iters
        if self.vtk_iter % 10 == 0 and to_vtk is not None:
            nodal_sol = {
                "ux": self.u[:, 0],
                "uy": self.u[:, 1],
                "x": self.x[:, 0],
                "rho": self.rho[:, 0],
            }

            to_vtk(
                conn,
                X,
                nodal_sol=nodal_sol,
                vtk_name=os.path.join(self.prefix, f"result_{self.vtk_iter}.vtk"),
            )

        self.vtk_iter += 1

        fail = 0
        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """

        monitor = 5
        max_iters = 100

        # Compliance gradient
        # Compute df/d(rho)
        self.dfdrho[:] = 0.0
        self.model.add_adjoint_dfdx(self.u_ref, self.dfdrho_ref)

        # Filter adjoint
        self.fltr_amg.cg(self.dfdrho_ref, self.fltr_res_ref, monitor, max_iters)

        self.dfdx[:] = 0.0
        self.fltr.add_adjoint_dfdx(self.fltr_res_ref, self.dfdx_ref)
        g[:] = self.compliance_scale * self.dfdx[:, 0]

        # Volume gradient
        self.dfdrho[:] = 0.0
        self.vol.eval_dfdx(self.dfdrho_ref)

        # Filter adjoint
        self.fltr_amg.cg(self.dfdrho_ref, self.fltr_res_ref, monitor, max_iters)

        self.dfdx[:] = 0.0
        self.fltr.add_adjoint_dfdx(self.fltr_res_ref, self.dfdx_ref)
        A[0][:] = (1.0 / self.vol_target) * self.dfdx[:, 0]

        fail = 0
        return fail


def X_conn_bcs_forces_2d():
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
            forces[nodes[-1, j], 1] = -1e3

    return X, conn, bcs, forces, nnodes


def X_conn_bcs_forces_3d():
    nx = 16
    ny = 16
    nz = 16

    nnodes = (nx + 1) * (ny + 1) * (nz + 1)
    nelems = nx * ny * nz
    nbcs = (ny + 1) * (nz + 1)

    conn = np.zeros((nelems, 8), dtype=np.int32)
    X = np.zeros((nnodes, 3), dtype=np.double)
    bcs = np.zeros((nbcs, 2), dtype=np.int32)
    forces = np.zeros((nnodes, 3), dtype=np.double)
    nodes = np.zeros((nx + 1, ny + 1, nz + 1), dtype=int)

    # Set the node locations
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nodes[i, j, k] = i + (nx + 1) * (j + (ny + 1) * k)
                X[nodes[i, j, k], 0] = i / nx
                X[nodes[i, j, k], 1] = j / ny
                X[nodes[i, j, k], 2] = k / nz

    # Set the connectivity
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elem = i + nx * (j + ny * k)
                conn[elem, 0] = nodes[i, j, k]
                conn[elem, 1] = nodes[i + 1, j, k]
                conn[elem, 2] = nodes[i + 1, j + 1, k]
                conn[elem, 3] = nodes[i, j + 1, k]

                conn[elem, 4] = nodes[i, j, k + 1]
                conn[elem, 5] = nodes[i + 1, j, k + 1]
                conn[elem, 6] = nodes[i + 1, j + 1, k + 1]
                conn[elem, 7] = nodes[i, j + 1, k + 1]

    # Set the boundary conditions, using the bits of the value to indicate
    # which values to fix
    bcs_val = 1 | 2 | 4  # 2^0 | 2^1 | 2^2
    index = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            bcs[index, 0] = nodes[0, j, k]
            bcs[index, 1] = bcs_val
            index += 1

    forces[nodes[-1, -1, 0], 2] = 1e3
    forces[nodes[-1, 0, -1], 1] = -1e3

    return X, conn, bcs, forces, nnodes


if __name__ == "__main__":
    X, conn, bcs, forces, nnodes = X_conn_bcs_forces_2d()

    # Set the constitutive data
    q = 5.0
    E = 70e3
    nu = 0.3
    density = 1.0
    design_stress = 1e3

    # Get the model data
    model = example.Elasticity_Model(X, bcs)

    # Add the element to the model
    elem = example.Elasticity_CPS4(conn)
    model.add_element(elem)

    # Add the constitutive object to the model
    con = example.TopoIsoConstitutive_CPS4(elem, q, E, nu, density, design_stress)
    model.add_constitutive(con)

    # Initialize the model
    model.init()

    # Set up the volume functional
    volume = example.Elasticity_Functional()
    volume.add_functional(example.TopoVolume_CPS4(con))

    # Get the model data
    fltr = example.Helmholtz_Model(X)

    # Set up the element
    length_scale = 0.05
    r0 = length_scale / (2.0 * np.sqrt(3.0))
    helmholtz_elem = example.Helmholtz_CPS4(conn, r0)

    # Add the model to the element
    fltr.add_element(helmholtz_elem)

    # Set the constitutive model
    helmholtz_con = example.HelmholtzConstitutive_CPS4(helmholtz_elem)
    fltr.add_constitutive(helmholtz_con)

    # Initialize the model
    fltr.init()

    # Create the optimization problem
    vol_init = 3.0
    vol_target = 0.4 * vol_init
    problem = TopOpt(model, fltr, volume, vol_target, nvars=nnodes)

    problem.f[:] = forces[:]

    problem.checkGradients()

    options = {"algorithm": "mma", "mma_max_iterations": 1000}

    # Set up the optimizer
    opt = ParOpt.Optimizer(problem, options)

    # Set a new starting point
    opt.optimize()
    x, z, zw, zl, zu = opt.getOptimizedPoint()

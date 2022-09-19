import numpy as np
from mpi4py import MPI
import argparse
from paropt import ParOpt
import os
from a2d import pya2d, to_vtk, InpParser


class TopOpt(ParOpt.Problem):
    def __init__(
        self,
        X,
        conn,
        bcs,
        vol_target,
        filter_length,
        q=5.0,
        E=70e3,
        nu=0.3,
        density=1.0,
        design_stress=0.27e3,
        num_levels=3,
        max_iters=200,
        monitor=10,
        omega=4.0 / 3.0,
        epsilon=0.0,
        vtk_prefix="",
    ):

        self.X = X
        self.conn = conn

        self.num_levels = num_levels
        self.monitor = monitor
        self.max_iters = max_iters
        self.omega = omega
        self.epsilon = epsilon
        self.vtk_prefix = vtk_prefix

        # Initialize the super class
        ncon = 1
        nvars = X.shape[0]
        super(TopOpt, self).__init__(MPI.COMM_SELF, nvars, ncon)

        # Create the model for solving the elasticity problem
        self.model = pya2d.Elasticity_Model(X, bcs)

        # Set up the volume functional
        self.vol = pya2d.Elasticity_Functional()

        # Add the element to the model
        if "C3D8" in conn:
            element = pya2d.Elasticity_C3D8(conn["C3D8"].astype(np.int32))
            self.model.add_element(element)
            con = pya2d.TopoIsoConstitutive_C3D8(
                element, q, E, nu, density, design_stress
            )
            self.model.add_constitutive(con)
            self.vol.add_functional(pya2d.TopoVolume_C3D8(con))
        elif "C3D8R" in conn:
            element = pya2d.Elasticity_C3D8(conn["C3D8R"].astype(np.int32))
            self.model.add_element(element)
            con = pya2d.TopoIsoConstitutive_C3D8(
                element, q, E, nu, density, design_stress
            )
            self.model.add_constitutive(con)
            self.vol.add_functional(pya2d.TopoVolume_C3D8(con))
        elif "C3D10" in conn:
            element = pya2d.Elasticity_C3D10(conn["C3D10"].astype(np.int32))
            self.model.add_element(element)
            con = pya2d.TopoIsoConstitutive_C3D10(
                element, q, E, nu, density, design_stress
            )
            self.model.add_constitutive(con)
            self.vol.add_functional(pya2d.TopoVolume_C3D10(con))

        # Initialize the model
        self.model.init()

        # Create the filter model
        self.fltr = pya2d.Helmholtz_Model(X)

        # Set up the element
        r0 = filter_length / (2.0 * np.sqrt(3.0))

        if "C3D8" in conn:
            element = pya2d.Helmholtz_C3D8(conn["C3D8"].astype(np.int32), r0)
            self.fltr.add_element(element)
            con = pya2d.HelmholtzConstitutive_C3D8(element)
            self.fltr.add_constitutive(con)
        elif "C3D8R" in conn:
            element = pya2d.Helmholtz_C3D8(conn["C3D8R"].astype(np.int32), r0)
            self.fltr.add_element(element)
            con = pya2d.HelmholtzConstitutive_C3D8(element)
            self.fltr.add_constitutive(con)
        elif "C3D10" in conn:
            element = pya2d.Helmholtz_C3D10(conn["C3D10"].astype(np.int32), r0)
            self.fltr.add_element(element)
            con = pya2d.HelmholtzConstitutive_C310(element)
            self.fltr.add_constitutive(con)

        # Initialize the model
        self.fltr.init()

        self.vol_target = vol_target

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
        print_info = True
        self.fltr_amg = self.fltr.new_amg(
            self.num_levels, self.omega, self.epsilon, self.Kf, print_info
        )

        # Set the scaling for the compliance
        self.compliance_scale = None

        # Set the iteration count
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
        self.fltr_res_ref.scale(-1.0)
        self.fltr_amg.cg(self.fltr_res_ref, self.rho_ref, self.monitor, self.max_iters)

        # Set the design variable values
        self.model.set_design_vars(self.rho_ref)

        # Set up AMG for the structural problem
        self.model.jacobian(self.K)
        print_info = True
        amg = self.model.new_amg(
            self.num_levels, self.omega, self.epsilon, self.K, print_info
        )

        amg.cg(self.f_ref, self.u_ref, self.monitor, self.max_iters)

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
            for key in self.conn:
                # nodal_sol = {
                #     "ux": self.u[:, 0],
                #     "uy": self.u[:, 1],
                #     "uz": self.u[:, 2],
                #     "dv": self.x[:, 0],
                #     "rho": self.rho[:, 0],
                # }

                soln = np.array(amg.get_vec(), copy=False)
                nodal_sol = {
                    "ux": soln[:, 0],
                    "uy": soln[:, 1],
                    "uz": soln[:, 2],
                    "dv": self.x[:, 0],
                    "rho": self.rho[:, 0],
                }

                vtk_name = os.path.join(
                    self.vtk_prefix, f"result_{key}_{self.vtk_iter}.vtk"
                )
                to_vtk(self.conn[key], self.X, nodal_sol=nodal_sol, vtk_name=vtk_name)

        self.vtk_iter += 1

        exit(0)

        fail = 0
        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """

        # Compliance gradient
        # Compute df/d(rho)
        self.dfdrho[:] = 0.0
        self.model.add_adjoint_dfdx(self.u_ref, self.dfdrho_ref)

        # Filter adjoint
        self.fltr_amg.cg(
            self.dfdrho_ref, self.fltr_res_ref, self.monitor, self.max_iters
        )

        self.dfdx[:] = 0.0
        self.fltr.add_adjoint_dfdx(self.fltr_res_ref, self.dfdx_ref)
        g[:] = self.compliance_scale * self.dfdx[:, 0]

        # Volume gradient
        self.dfdrho[:] = 0.0
        self.vol.eval_dfdx(self.dfdrho_ref)

        # Filter adjoint
        self.fltr_amg.cg(
            self.dfdrho_ref, self.fltr_res_ref, self.monitor, self.max_iters
        )

        self.dfdx[:] = 0.0
        self.fltr.add_adjoint_dfdx(self.fltr_res_ref, self.dfdx_ref)
        A[0][:] = (1.0 / self.vol_target) * self.dfdx[:, 0]

        fail = 0
        return fail


parser = argparse.ArgumentParser(description="Perform compliance minimization")
parser.add_argument("--inp", type=str, default="", help="Input file")
parser.add_argument("--opt", type=str, default="mma", help="Optimizer type")
parser.add_argument("--vol_target", type=float, default=0.4, help="Volume target")
parser.add_argument("--filter_length", type=float, default=0.025, help="Filter length")
parser.add_argument(
    "--num_levels", type=int, default=3, help="Number of multigrid levels"
)
parser.add_argument("--vtk_prefix", type=str, default="", help="File output prefix")
args = parser.parse_args()

filename = args.inp
vol_target = args.vol_target
filter_length = args.filter_length
num_levels = args.num_levels
vtk_prefix = args.vtk_prefix

if os.path.isfile(filename):
    ndof = 3  # Assumed to be a 3d problem with u, v, w unknowns
    inp_parser = InpParser(filename)
    conn, X, groups, loads_list, bcs_list = inp_parser.parse()

    # Set the boundary conditions
    bcs = []
    for bc in bcs_list:
        # Figure out the set of nodes on the boundary
        nset = []
        if isinstance(bc[0], str):
            nset = groups[bc[0]]
        else:
            nset = [bc[0]]

        # Figure out which dofs to constrain
        bcs_val = 0
        if isinstance(bc[1], str):
            # All degrees of freedom constrained
            for i in range(ndof):
                bcs_val = bcs_val | 1 << i
        else:
            for dof in bc[1:]:
                bcs_val = bcs_val | 1 << dof

        for node in nset:
            bcs.append([node, bcs_val])

    problem = TopOpt(
        X,
        conn,
        bcs,
        vol_target,
        filter_length,
        q=5.0,
        E=70e3,
        nu=0.3,
        density=1.0,
        design_stress=0.27e3,
        num_levels=num_levels,
        vtk_prefix=vtk_prefix,
        omega=4.0 / 3.0,
    )

    # Set the forces
    for load in loads_list:
        # Find the node set for the load
        nset = []
        if isinstance(load[0], str):
            nset = groups[load[0]]
        else:
            nset = [load[0]]

        # Find the dof for the load
        dof = load[1]

        # Find the load amplitude
        amp = load[2]

        for node in nset:
            problem.f[node, dof] += amp

else:
    nx = 96
    ny = 48
    nz = 48

    lx = 2.0
    ly = 1.0
    lz = 1.0

    nnodes = (nx + 1) * (ny + 1) * (nz + 1)
    nelems = nx * ny * nz
    nbcs = (ny + 1) * (nz + 1)

    conn_c3d8 = np.zeros((nelems, 8), dtype=np.int32)
    X = np.zeros((nnodes, 3), dtype=np.double)
    bcs = np.zeros((nbcs, 2), dtype=np.int32)

    # Set the node locations
    nodes = np.zeros((nx + 1, ny + 1, nz + 1), dtype=int)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nodes[i, j, k] = i + (nx + 1) * (j + (ny + 1) * k)
                X[nodes[i, j, k], 0] = lx * i / nx
                X[nodes[i, j, k], 1] = ly * j / ny
                X[nodes[i, j, k], 2] = lz * k / nz

    # Set the connectivity
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elem = i + nx * (j + ny * k)
                conn_c3d8[elem, 0] = nodes[i, j, k]
                conn_c3d8[elem, 1] = nodes[i + 1, j, k]
                conn_c3d8[elem, 2] = nodes[i + 1, j + 1, k]
                conn_c3d8[elem, 3] = nodes[i, j + 1, k]

                conn_c3d8[elem, 4] = nodes[i, j, k + 1]
                conn_c3d8[elem, 5] = nodes[i + 1, j, k + 1]
                conn_c3d8[elem, 6] = nodes[i + 1, j + 1, k + 1]
                conn_c3d8[elem, 7] = nodes[i, j + 1, k + 1]

    # Set the boundary conditions, using the bits of the value to indicate
    # which values to fix
    bcs_val = 1 | 2 | 4  # 2^0 | 2^1 | 2^2
    index = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            bcs[index, 0] = nodes[0, j, k]
            bcs[index, 1] = bcs_val
            index += 1

    conn = {"C3D8": conn_c3d8}

    vol_target = 0.4 * lx * ly * lz
    filter_length = 0.025 * ly

    problem = TopOpt(
        X,
        conn,
        bcs,
        vol_target,
        filter_length,
        q=5.0,
        E=70e3,
        nu=0.3,
        density=1.0,
        design_stress=0.27e3,
        vtk_prefix=vtk_prefix,
    )

    # Set the forces
    problem.f[nodes[-1, -1, 0], 2] = 1e3
    problem.f[nodes[-1, 0, -1], 1] = -1e3

problem.checkGradients()

if args.opt == "tr":
    options = {
        "algorithm": "tr",
        "tr_init_size": 0.05,
        "tr_min_size": 1e-6,
        "tr_max_size": 10.0,
        "tr_eta": 0.25,
        "tr_infeas_tol": 1e-6,
        "tr_l1_tol": 1e-6,
        "tr_linfty_tol": 0.0,
        "tr_adaptive_gamma_update": True,
        "tr_max_iterations": 2,
        "max_major_iters": 100,
        "penalty_gamma": 1e3,
        "qn_subspace_size": 2,
        "hessian_reset_freq": 2,
        "qn_type": "bfgs",
        "abs_res_tol": 1e-8,
        "starting_point_strategy": "affine_step",
        "barrier_strategy": "mehrotra_predictor_corrector",
        "use_line_search": False,
    }
else:
    options = {"algorithm": "mma", "mma_max_iterations": 5}

# Set up the optimizer
opt = ParOpt.Optimizer(problem, options)

# Set a new starting point
opt.optimize()
x, z, zw, zl, zu = opt.getOptimizedPoint()

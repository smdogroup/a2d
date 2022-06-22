import numpy as np
from mpi4py import MPI
import example
from paropt import ParOpt


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
        self.amg_f = self.fltr.new_amg(self.num_levels, omega, self.Kf, print_info)

        # Set the scaling for the compliance
        self.compliance_scale = None

        return

    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        lb[:] = 1e-3
        ub[:] = 1.0
        x[:] = 0.95
        return

    def evalObjCon(self, x):
        """
        Return the objective, constraint and fail flag
        """

        # Set the design variables
        self.x[:, 0] = x[:]
        monitor = 5
        max_iters = 100
        self.amg_f.cg(self.x_ref, self.rho_ref, monitor, max_iters)

        # Set the design variable values
        self.model.set_design_vars(self.rho_ref)

        # Set up the AMG for the structural problem
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
        con = np.array([1.0 - volume / self.vol_target])

        fail = 0
        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """

        self.dfdrho[:] = 0.0
        self.model.add_adjoint_dfdx(self.u_ref, self.dfdrho_ref)

        monitor = 5
        max_iters = 100
        self.amg_f.cg(self.dfdrho_ref, self.dfdx_ref, monitor, max_iters)
        g[:] = -self.compliance_scale * self.dfdx[:, 0]

        self.dfdrho[:] = 0.0
        self.vol.eval_dfdx(self.dfdrho_ref)
        self.amg_f.cg(self.dfdrho_ref, self.dfdx_ref, monitor, max_iters)
        A[0][:] = -(1.0 / self.vol_target) * self.dfdx[:, 0]

        fail = 0
        return fail


nx = 64
ny = 64
nz = 64
nnodes = (nx + 1) * (ny + 1) * (nz + 1)
nelems = nx * ny * nz
nbcs = (ny + 1) * (nz + 1)

conn = np.zeros((nelems, 8), dtype=np.int32)
X = np.zeros((nnodes, 3), dtype=np.double)
bcs = np.zeros((nbcs, 2), dtype=np.int32)

# Set the node locations
nodes = np.zeros((nx + 1, ny + 1, nz + 1), dtype=int)
for k in range(nz + 1):
    for j in range(ny + 1):
        for i in range(nx + 1):
            nodes[i, j, k] = i + (nx + 1) * (j + (ny + 1) * k)
            X[nodes[i, j, k], 0] = 1.0 * i / nx
            X[nodes[i, j, k], 1] = 1.0 * j / ny
            X[nodes[i, j, k], 2] = 1.0 * k / nz

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

# Set the constitutive data
q = 5.0
E = 70e3
nu = 0.3
density = 1.0
design_stress = 1e3

# Get the model data
model = example.Elasticity_Model(X, bcs)

# Add the element to the model
hex = example.Elasticity_C3D8(conn)
model.add_element(hex)

# Add the constitutive object to the model
con = example.TopoIsoConstitutive_C3D8(hex, q, E, nu, density, design_stress)
model.add_constitutive(con)

# Initialize the model
model.init()

# Get the model data
fltr = example.Helmholtz_Model(X)

# Set up the element
length_scale = 0.05
r0 = length_scale / (2.0 * np.sqrt(3.0))
helmholtz_hex = example.Helmholtz_C3D8(conn, r0)

# Add the model to the element
fltr.add_element(helmholtz_hex)

# Initialize the model
fltr.init()

# Set up the volume functional
volume = example.Elasticity_Functional()
volume.add_functional(example.TopoVolume_C3D8(con))

vol_init = 1.0
vol_target = 0.4 * vol_init
problem = TopOpt(model, fltr, volume, vol_target, nvars=nnodes)

problem.f[nodes[-1, 0, 0], 2] = 1e3
problem.f[nodes[-1, -1, -1], 1] = 1e3

problem.checkGradients()

options = {
    "algorithm": "tr",
    "tr_init_size": 0.05,
    "tr_min_size": 1e-6,
    "tr_max_size": 10.0,
    "tr_eta": 0.25,
    "tr_infeas_tol": 1e-6,
    "tr_l1_tol": 1e-3,
    "tr_linfty_tol": 0.0,
    "tr_adaptive_gamma_update": True,
    "tr_max_iterations": 1000,
    "max_major_iters": 100,
    "penalty_gamma": 1e3,
    "qn_subspace_size": 10,
    "qn_type": "bfgs",
    "abs_res_tol": 1e-8,
    "starting_point_strategy": "affine_step",
    "barrier_strategy": "mehrotra_predictor_corrector",
    "use_line_search": False,
}

# Set up the optimizer
opt = ParOpt.Optimizer(problem, options)

# Set a new starting point
opt.optimize()
x, z, zw, zl, zu = opt.getOptimizedPoint()

# stress = example.Elasticity_Functional()
# stress.add_functional(example.TopoVonMisesAggregation_C3D8(con, 100.0))

# xref = helmholtz_model.new_solution()
# x = np.array(xref, copy=False)
# x[:] = 0.5
# model.set_design_vars(xref)

# # Create the new matrix
# J = model.new_matrix()
# model.jacobian(J)

# # Set the AMG information
# num_levels = 3
# omega = 4.0 / 3.0
# print_info = True
# amg = model.new_amg(num_levels, omega, J, print_info)

# # Create the res and solution arrays
# res = model.new_solution()
# res_array = np.array(res, copy=False)

# ans = model.new_solution()
# ans_array = np.array(ans, copy=False)

# res_array[:] = 1.0
# res_array[bcs[:, 0], :] = 0.0

# monitor = 5
# max_iters = 100
# amg.cg(res, ans, monitor, max_iters)
# model.set_solution(ans)

# print("Volume = ", volume.eval_functional())
# print("Von Mises = ", stress.eval_functional())

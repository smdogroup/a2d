import numpy as np
import sys
import example

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

# Set up the element
hex = example.Elasticity_C3D8(conn)
con = example.TopoIsoConstitutive_C3D8(hex, q, E, nu, density, design_stress)

# Add the model to the element
model.add_element(hex)

# Add the constitutive class
model.add_constitutive(con)

# Initialize the model
model.init()

# Get the model data
helmholtz_model = example.Helmholtz_Model(X, bcs)

# Set up the element
length_scale = 0.25
r0 = length_scale / (2.0 * np.sqrt(3.0))
helmholtz_hex = example.Helmholtz_C3D8(conn, r0)

# Add the model to the element
helmholtz_model.add_element(helmholtz_hex)

# Initialize the model
helmholtz_model.init()

# Set up the volume functional
volume = example.Elasticity_Functional()
volume.add_functional(example.TopoVolume_C3D8(con))

stress = example.Elasticity_Functional()
stress.add_functional(example.TopoVonMisesAggregation_C3D8(con, 100.0))

xref = helmholtz_model.new_solution()
x = np.array(xref, copy=False)
x[:] = 1.0
model.set_design_vars(xref)

# Create the new matrix
J = model.new_matrix()
model.jacobian(J)

# Set the AMG information
num_levels = 3
omega = 4.0 / 3.0
print_info = True
amg = model.new_amg(num_levels, omega, J, print_info)

# Create the res and solution arrays
res = model.new_solution()
res_array = np.array(res, copy=False)

ans = model.new_solution()
ans_array = np.array(ans, copy=False)

res_array[:] = 1.0
res_array[bcs[:, 0], :] = 0.0

print('Volume = ', volume.eval_functional())

monitor = 5
max_iters = 100
amg.cg(res, ans, monitor, max_iters)
amg.mg(res, ans, monitor, max_iters)

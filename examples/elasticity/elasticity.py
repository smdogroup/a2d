import numpy as np
import sys
import example

nx = 64
ny = 64
nz = 64
nnodes = (nx + 1) * (ny + 1) * (nz + 1)
nelems = nx * ny * nz
nbcs = (ny + 1) * (nz + 1)

model = example.ElasticityModel(nnodes, nbcs)

# Get the model data
bcs = np.array(model.get_bcs(), copy=False)
X_ = model.get_nodes()
X = np.array(X_, copy=False)

# Set up the element
hex = example.ElasticityHexElement(nelems)

# Add the model to the element
model.add_element(hex)

# Get the data from the element
conn = np.array(hex.get_conn(), copy=False)
data = np.array(hex.get_quad_data(), copy=False)

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

if data.shape[2] == 2:
    data[:, :, 0] = 1.23
    data[:, :, 1] = 2.45
else:
    data[:] = 1.0

# Set the nodes
model.set_nodes(X_)

J = model.new_matrix()
model.jacobian(J)

# Set the AMG information
num_levels = 4
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

monitor = 5
max_iters = 100
amg.cg(res, ans, monitor, max_iters)
amg.mg(res, ans, monitor, max_iters)

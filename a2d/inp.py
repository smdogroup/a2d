import numpy as np
from parse_inp import InpParser, to_vtk
import example

# Parse the inp file
inp_file = "sampleFile.inp"
inp_parser = InpParser(inp_file)
conn, X, groups = inp_parser.parse()

# Set the number of nodes
nnodes = X.shape[0]

# Select the node set
bcs = groups["back_face_nodes"]

hex_conn = conn["C3D8R"]
tet_conn = conn["C3D10"]

# Set the material properties
E = 70e3
nu = 0.3
mu = 0.5 * E / (1.0 + nu)
lambda_ = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))

# Allocate the type of BCs to use
nbcs = len(bcs)
mbcs = np.zeros((nbcs, 2), dtype=np.int32)

# Set the boundary conditions
bcs_val = 1 | 2 | 4  # 2^0 | 2^1 | 2^2 - each indicates a fixed dof
for index, k in enumerate(bcs):
    mbcs[index, 0] = k
    mbcs[index, 1] = bcs_val

# Create the model class
model = example.ElasticityModel(X, mbcs)

hex = example.ElasticityHexElement(hex_conn.astype(np.int32))
hex_data = np.array(hex.get_quad_data(), copy=False)
hex_data[:, :, 0] = mu
hex_data[:, :, 1] = lambda_

tet = example.ElasticityTetElement(tet_conn.astype(np.int32))
tet_data = np.array(tet.get_quad_data(), copy=False)
tet_data[:, :, 0] = mu
tet_data[:, :, 1] = lambda_

# Add the elements to the model
model.add_element(hex)
model.add_element(tet)
model.init()

J = model.new_matrix()
model.jacobian(J)

# Set the AMG information
num_levels = 2
omega = 4.0 / 3.0
print_info = True
amg = model.new_amg(num_levels, omega, J, print_info)

# Create the res and solution arrays
res = model.new_solution()
res_array = np.array(res, copy=False)

ans = model.new_solution()
ans_array = np.array(ans, copy=False)

res_array[:] = 1.0
res_array[bcs, :] = 0.0

monitor = 5
max_iters = 500
amg.cg(res, ans, monitor, max_iters)

nodal_sol = {"u": ans_array[:, 0], "v": ans_array[:, 1], "w": ans_array[:, 2]}
to_vtk(conn, X, nodal_sol, vtk_name="output.vtk")

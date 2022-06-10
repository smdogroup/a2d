import numpy as np

def print_flatten(mat):
  for v in mat.flatten():
    print(f'{v:.16f}, ', end='')
  print()

A = np.array([0.96155402, 0.02855176, 0.95787560, 0.45439794]).reshape(2,2)
B = np.array([0.80766462, 0.60212270, 0.86418474, 0.65304149]).reshape(2,2)

print_flatten(A.dot(B))
print_flatten(A.dot(B.T))
print_flatten(A.T.dot(B))
print_flatten(A.T.dot(B.T))

from scipy.io import mmread
import argparse
import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("mtx", nargs="*")
args = p.parse_args()

mat_list = []

for mtx in args.mtx:
    fig, ax = plt.subplots()
    mat = mmread(mtx)
    mat_list.append(mat)
    ax.spy(mat)

# Compute difference
if len(mat_list) == 2:
    mat1 = mat_list[0]
    mat2 = mat_list[1]

    diff = np.max(np.abs(mat1 - mat2))

    print("max entry of mat 1: %20.10f" % np.max(np.abs(mat1)))
    print("max entry of mat 2: %20.10f" % np.max(np.abs(mat2)))
    print("max entry diff: %20.10f" % diff)

    fig, ax = plt.subplots()
    ax.matshow((mat1 - mat2).todense())

plt.show()

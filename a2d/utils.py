import numpy as np
import os


def to_vtk(conn, X, nodal_sol={}, vtk_name="problem.vtk"):
    """
    Generate a vtk given conn, X, and optionally nodal_sol

    Inputs:
        nnodes: ndarray
        conn: ndarray or dictionary if a mixed mesh is used
        X: ndarray
        nodal_sol: nodal solution values, dictionary with the following
                    structure:

                    nodal_sol = {
                        "scalar1": [...],
                        "scalar2": [...],
                        ...
                    }

        vtk_name: name of the vtk
    """
    ELEMENT_INFO = {
        "CPS3": {
            "nnode": 3,
            "vtk_type": 5,
            "note": "Three-node plane stress element",
        },
        "C3D4": {
            "nnode": 4,
            "vtk_type": 10,
            "note": "Four-node tetrahedral element",
        },
        "C3D8R": {
            "nnode": 8,
            "vtk_type": 12,
            "note": "general purpose linear brick element",
        },
        "C3D10": {
            "nnode": 10,
            "vtk_type": 24,
            "note": "Ten-node tetrahedral element",
        },
        "tri": {
            "nnode": 3,
            "vtk_type": 5,
            "note": "triangle element",
        },
        "quad": {
            "nnode": 4,
            "vtk_type": 9,
            "note": "2d quadrilateral element",
        },
        "block": {
            "nnode": 8,
            "vtk_type": 12,
            "note": "3d block element",
        },
    }

    if isinstance(conn, np.ndarray):
        if conn.shape[1] == 3:
            conn = {"tri": conn}
        elif conn.shape[1] == 4:
            conn = {"quad": conn}
        elif conn.shape[1] == 8:
            conn = {"block": conn}

    # vtk requires a 3-dimensional data point
    if X.shape[1] == 2:
        X = np.append(X, np.zeros((X.shape[0], 1)), axis=1)

    nnodes = X.shape[0]
    nelems = np.sum([len(c) for c in conn.values()])

    # Create a empty vtk file and write headers
    with open(vtk_name, "w") as fh:
        fh.write("# vtk DataFile Version 3.0\n")
        fh.write("my example\n")
        fh.write("ASCII\n")
        fh.write("DATASET UNSTRUCTURED_GRID\n")

        # Write nodal points
        fh.write("POINTS {:d} double\n".format(nnodes))
        for x in X:
            row = f"{x}"[1:-1]  # Remove square brackets in the string
            fh.write(f"{row}\n")

        # Write connectivity
        size = np.sum(
            [
                len(econn) * (1 + ELEMENT_INFO[etype]["nnode"])
                for etype, econn in conn.items()
            ]
        )
        fh.write(f"CELLS {nelems} {size}\n")
        for etype, econn in conn.items():
            for c in econn:
                node_idx = f"{c}"[1:-1]  # remove square bracket [ and ]
                npts = ELEMENT_INFO[etype]["nnode"]
                fh.write(f"{npts} {node_idx}\n")

        # Write cell type
        fh.write(f"CELL_TYPES {nelems}\n")
        for etype, econn in conn.items():
            for c in econn:
                vtk_type = ELEMENT_INFO[etype]["vtk_type"]
                fh.write(f"{vtk_type}\n")

        # Write solution
        if nodal_sol:
            fh.write(f"POINT_DATA {nnodes}\n")
            for name, data in nodal_sol.items():
                fh.write(f"SCALARS {name} float 1\n")
                fh.write("LOOKUP_TABLE default\n")
                for val in data:
                    fh.write(f"{val}\n")
    print(f"[Info] Done generating {vtk_name}")
    return

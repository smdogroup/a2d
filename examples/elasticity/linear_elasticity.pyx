import numpy as np
cimport numpy as np

np.import_array()

cdef extern: # linear_elasticity.cpp":
    void compute_residual( int nelems,
                           int nnodes,
                           int *conn_data,
                           double *X_data,
                           double *mat_data,
                           double *U_data,
                           double *res_data )

    void compute_jacobian( int nelems,
                           int nnodes,
                           int *conn_data,
                           double *X_data,
                           double *mat_data,
                           double *U_data,
                           double *jac_data )

cdef compute_res(
        np.ndarray[int, ndim=2] conn,
        np.ndarray[float, ndim=2] X,
        np.ndarray[float, ndim=2] data,
        np.ndarray[float, ndim=2] U,
        np.ndarray[float, ndim=3] res):

    cdef int nelems = conn.shape[0]
    cdef int nnodes = X.shape[0]

    compute_residual(nelems, nnodes,
        <int*>conn.data, <double*>X.data, <double*>data.data,
        <double*>U.data, <double*>res.data)


cdef compute_jac(
        np.ndarray[int, ndim=2] conn,
        np.ndarray[float, ndim=2] X,
        np.ndarray[float, ndim=2] data,
        np.ndarray[float, ndim=2] U,
        np.ndarray[float, ndim=5] jac):

    cdef int nelems = conn.shape[0]
    cdef int nnodes = X.shape[0]

    compute_jacobian(nelems, nnodes,
        <int*>conn.data, <double*>X.data, <double*>data.data,
        <double*>U.data, <double*>jac.data)
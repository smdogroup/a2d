cimport numpy as np

np.import_array()

cdef extern from "elasticity_impl.h":
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

cpdef compute_res(
        np.ndarray[int, ndim=2] conn,
        np.ndarray[double, ndim=2] X,
        np.ndarray[double, ndim=3] data,
        np.ndarray[double, ndim=2] U,
        np.ndarray[double, ndim=3] res):

    cdef int nelems = conn.shape[0]
    cdef int nnodes = X.shape[0]

    compute_residual(nelems, nnodes,
        <int*>conn.data, <double*>X.data, <double*>data.data,
        <double*>U.data, <double*>res.data)


cpdef compute_jac(
        np.ndarray[int, ndim=2] conn,
        np.ndarray[double, ndim=2] X,
        np.ndarray[double, ndim=3] data,
        np.ndarray[double, ndim=2] U,
        np.ndarray[double, ndim=5] jac):

    cdef int nelems = conn.shape[0]
    cdef int nnodes = X.shape[0]

    compute_jacobian(nelems, nnodes,
        <int*>conn.data, <double*>X.data, <double*>data.data,
        <double*>U.data, <double*>jac.data)
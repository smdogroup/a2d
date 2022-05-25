#include "elasticity3d.h"
#include "multiarray.h"
using namespace A2D;

void compute_residual(int nelems,
                      int nnodes,
                      int *conn_data,
                      double *X_data,
                      double *mat_data,
                      double *U_data,
                      double *res_data);

void compute_jacobian(int nelems,
                      int nnodes,
                      int *conn_data,
                      double *X_data,
                      double *mat_data,
                      double *U_data,
                      double *jac_data);
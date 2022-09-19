#ifndef A2D_H
#define A2D_H

#include "a2dmatcore2d.h"
#include "a2dmatcore3d.h"
#include "a2dmemory.h"
#include "a2dobjs.h"
#include "a2dtmp2d.h"
#include "a2dtmp3d.h"
#include "a2dtypes.h"
#include "block_numeric.h"
#include "fem/basis.h"
#include "fem/constitutive.h"
#include "fem/elasticity.h"
#include "fem/element.h"
#include "fem/functional.h"
#include "fem/helmholtz.h"
#include "fem/model.h"
#include "fem/quadrature.h"
#include "multiarray.h"
#include "parallel.h"
#include "sparse/sparse_amd.h"
#include "sparse/sparse_amg.h"
#include "sparse/sparse_matrix.h"
#include "sparse/sparse_numeric.h"
#include "sparse/sparse_symbolic.h"
#include "utils/a2dmesh.h"
#include "utils/a2dprofiler.h"
#include "utils/a2dvtk.h"

#ifdef A2D_USE_KOKKOS
#include "a2dkokkos.h"
#endif

#endif  // A2D_H
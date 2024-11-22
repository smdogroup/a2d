#ifndef A2D_CORE_H
#define A2D_CORE_H

// Key objects

#include "ad/a2dmat.h"
#include "ad/a2dobj.h"
#include "ad/a2dstack.h"
#include "ad/a2dvec.h"

// Operations

#include "ad/a2dgemm.h"
#include "ad/a2dgreenstrain.h"
#include "ad/a2dhadamard.h"
#include "ad/a2disotropic.h"
#include "ad/a2dmatdet.h"
#include "ad/a2dmatinv.h"
#include "ad/a2dmatsum.h"
#include "ad/a2dmattovec.h"
#include "ad/a2dmattrace.h"
#include "ad/a2dmatvecmult.h"
#include "ad/a2dquaternion.h"
#include "ad/a2dscalarops.h"
#include "ad/a2dsymeigs.h"
#include "ad/a2dsymmatmulttrace.h"
#include "ad/a2dsymrk.h"
#include "ad/a2dsymsum.h"
#include "ad/a2dveccross.h"
#include "ad/a2dvecnorm.h"
#include "ad/a2dvecouter.h"
#include "ad/a2dvecsum.h"

// shell routines
#include "ad/shell/a2dshellassembleframe.h"
// #include "ad/shell/a2dmatconcat.h"
#include "ad/shell/a2dmatrotateframe.h"
#include "ad/shell/a2dshellstrain.h"
#include "ad/shell/a2dsymmatrotateframe.h"

#endif  //  A2D_CORE_H

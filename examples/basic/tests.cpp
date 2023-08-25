#include "ad/a2dgemm.h"
#include "ad/a2dmatdet.h"
#include "ad/a2dmatinv.h"
#include "ad/a2dmattrace.h"

int main() {
  A2D::Test::MatMatMultTestAll();
  A2D::Test::MatDetTestAll();
  A2D::Test::MatInvTestAll();
  A2D::Test::MatTraceTestAll();

  return (0);
}
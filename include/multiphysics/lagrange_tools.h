#ifndef A2D_LAGRANGE_TOOLS
#define A2D_LAGRANGE_TOOLS

namespace A2D {

static const double GaussLobattoPoints2[] = {-1.0, 1.0};

static const double GaussLobattoPoints3[] = {-1.0, 0.0, 1.0};

static const double GaussLobattoPoints4[] = {-1.0, -0.5, 0.5, 1.0};

static const double GaussLobattoPoints5[] = {-1.0, -0.7071067811865475, 0.0,
                                             0.7071067811865475, 1.0};

static const double GaussLobattoPoints6[] = {-1.0,
                                             -0.8090169943749475,
                                             -0.30901699437494745,
                                             0.30901699437494745,
                                             0.8090169943749475,
                                             1.0};

template <index_t order>
void lagrange_basis(const double pt, double N[], double Nx[]) {}

template <>
void lagrange_basis<2u>(const double pt, double N[], double Nx[]) {
  N[0] = 0.5 * (1.0 - pt);
  N[1] = 0.5 * (pt - 1.0);

  Nx[0] = -0.5;
  Nx[1] = 0.5;
}

template <>
void lagrange_basis<3u>(const double pt, double N[], double Nx[]) {
  N[0] = -0.5 * pt * (1.0 - pt);
  N[1] = (1.0 - pt) * (1.0 + pt);
  N[2] = 0.5 * (1.0 + pt) * pt;

  Nx[0] = -0.5 + pt;
  Nx[1] = -2.0 * pt;
  Nx[2] = 0.5 + pt;
}

template <>
void lagrange_basis<4u>(const double pt, double N[], double Nx[]) {
  N[0] = -(2.0 / 3.0) * (0.5 + pt) * (0.5 - pt) * (1.0 - pt);
  N[1] = (4.0 / 3.0) * (1.0 + pt) * (0.5 - pt) * (1.0 - pt);
  N[2] = (4.0 / 3.0) * (1.0 + pt) * (0.5 + pt) * (1.0 - pt);
  N[3] = -(2.0 / 3.0) * (1.0 + pt) * (0.5 + pt) * (0.5 - pt);

  Nx[0] = -2.0 * pt * pt + (4.0 / 3.0) * pt + 1.0 / 6.0;
  Nx[1] = 4.0 * pt * pt - (4.0 / 3.0) * pt - 4.0 / 3.0;
  Nx[2] = -4.0 * pt * pt - (4.0 / 3.0) * pt + 4.0 / 3.0;
  Nx[3] = 2.0 * pt * pt + (4.0 / 3.0) * pt - 1.0 / 6.0;
}

}  // namespace A2D

#endif  //  A2D_LAGRANGE_TOOLS
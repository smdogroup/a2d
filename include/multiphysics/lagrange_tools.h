#ifndef A2D_LAGRANGE_TOOLS
#define A2D_LAGRANGE_TOOLS

namespace A2D {

const double GaussQuadPts1[] = {0.0};
const double GaussQuadWts1[] = {2.0};

const double GaussQuadPts2[] = {-0.577350269189626, 0.577350269189626};
const double GaussQuadWts2[] = {1.0, 1.0};

const double GaussQuadPts3[] = {-0.774596669241483, 0.0, 0.774596669241483};
const double GaussQuadWts3[] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

const double GaussQuadPts4[] = {-0.861136311594053, -0.339981043584856,
                                0.339981043584856, 0.861136311594053};
const double GaussQuadWts4[] = {0.347854845137454, 0.652145154862546,
                                0.652145154862546, 0.347854845137454};

const double GaussQuadPts5[] = {-0.906179845938664, -0.538469310105683, 0.0,
                                0.538469310105683, 0.906179845938664};
const double GaussQuadWts5[] = {0.236926885056189, 0.478628670499366,
                                0.568888888888889, 0.478628670499366,
                                0.236926885056189};

const double GaussQuadPts6[] = {
    -0.9324695142031520278123016, -0.6612093864662645136613996,
    -0.2386191860831969086305017, 0.2386191860831969086305017,
    0.6612093864662645136613996,  0.9324695142031520278123016};
const double GaussQuadWts6[] = {
    0.1713244923791703450402961, 0.3607615730481386075698335,
    0.4679139345726910473898703, 0.4679139345726910473898703,
    0.3607615730481386075698335, 0.1713244923791703450402961};

const double GaussQuadPts7[] = {
    -0.9491079123427585245261897, -0.7415311855993944398638648,
    -0.4058451513773971669066064, 0.0,
    0.4058451513773971669066064,  0.7415311855993944398638648,
    0.9491079123427585245261897};
const double GaussQuadWts7[] = {
    0.1294849661688696932706114, 0.2797053914892766679014678,
    0.3818300505051189449503698, 0.4179591836734693877551020,
    0.3818300505051189449503698, 0.2797053914892766679014678,
    0.1294849661688696932706114};

const double GaussQuadPts8[] = {
    -0.9602898564975362316835609, -0.7966664774136267395915539,
    -0.5255324099163289858177390, -0.1834346424956498049394761,
    0.1834346424956498049394761,  0.5255324099163289858177390,
    0.7966664774136267395915539,  0.9602898564975362316835609};
const double GaussQuadWts8[] = {
    0.1012285362903762591525314, 0.2223810344533744705443560,
    0.3137066458778872873379622, 0.3626837833783619829651504,
    0.3626837833783619829651504, 0.3137066458778872873379622,
    0.2223810344533744705443560, 0.1012285362903762591525314};

template <index_t order>
constexpr const double* get_gauss_quadrature_wts() {
  if (order == 1) {
    return GaussQuadWts1;
  } else if (order == 2) {
    return GaussQuadWts2;
  } else if (order == 3) {
    return GaussQuadWts3;
  } else if (order == 4) {
    return GaussQuadWts4;
  } else if (order == 5) {
    return GaussQuadWts5;
  } else if (order == 6) {
    return GaussQuadWts6;
  } else if (order == 7) {
    return GaussQuadWts7;
  } else if (order == 8) {
    return GaussQuadWts8;
  }
  return NULL;
}

template <index_t order>
constexpr const double* get_gauss_quadrature_pts() {
  if (order == 1) {
    return GaussQuadPts1;
  } else if (order == 2) {
    return GaussQuadPts2;
  } else if (order == 3) {
    return GaussQuadPts3;
  } else if (order == 4) {
    return GaussQuadPts4;
  } else if (order == 5) {
    return GaussQuadPts5;
  } else if (order == 6) {
    return GaussQuadPts6;
  } else if (order == 7) {
    return GaussQuadPts7;
  } else if (order == 8) {
    return GaussQuadPts8;
  }
  return NULL;
}

const double GaussLobattoPoints2[] = {-1.0, 1.0};

const double GaussLobattoPoints3[] = {-1.0, 0.0, 1.0};

const double GaussLobattoPoints4[] = {-1.0, -0.5, 0.5, 1.0};

const double GaussLobattoPoints5[] = {-1.0, -0.7071067811865475, 0.0,
                                      0.7071067811865475, 1.0};

const double GaussLobattoPoints6[] = {-1.0,
                                      -0.8090169943749475,
                                      -0.30901699437494745,
                                      0.30901699437494745,
                                      0.8090169943749475,
                                      1.0};

template <index_t order>
constexpr const double* get_gauss_lobatto_pts() {
  if (order == 2) {
    return GaussLobattoPoints2;
  } else if (order == 3) {
    return GaussLobattoPoints3;
  } else if (order == 4) {
    return GaussLobattoPoints4;
  } else if (order == 5) {
    return GaussLobattoPoints5;
  } else if (order == 6) {
    return GaussLobattoPoints6;
  }
  return NULL;
}

template <index_t order>
void lagrange_basis(const double knots[], const double pt, double N[]) {
  // Loop over the shape functions
  for (int i = 0; i < order; i++) {
    N[i] = 1.0;
    for (int j = 0; j < order; j++) {
      if (i != j) {
        double d = 1.0 / (knots[i] - knots[j]);
        N[i] *= (pt - knots[j]) * d;
      }
    }
  }
}

template <index_t order>
void lagrange_basis(const double knots[], const double pt, double N[],
                    double Nx[]) {
  // Loop over the shape function knot locations
  for (int i = 0; i < order; i++) {
    N[i] = 1.0;
    Nx[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for (int j = 0; j < order; j++) {
      if (i != j) {
        double d = 1.0 / (knots[i] - knots[j]);
        N[i] *= (pt - knots[j]) * d;

        // Now add up the contribution to the derivative
        for (int k = 0; k < order; k++) {
          if (k != i && k != j) {
            d *= (pt - knots[k]) / (knots[i] - knots[k]);
          }
        }

        // Add the derivative contribution
        Nx[i] += d;
      }
    }
  }
}

template <index_t order>
void lagrange_basis(const double pt, double N[]) {
  constexpr const double* knots = get_gauss_lobatto_pts<order>();
  lagrange_basis<order>(knots, pt, N);
}

template <index_t order>
void lagrange_basis(const double pt, double N[], double Nx[]) {
  constexpr const double* knots = get_gauss_lobatto_pts<order>();
  lagrange_basis<order>(knots, pt, N, Nx);
}

template <>
void lagrange_basis<1u>(const double pt, double N[]) {
  N[0] = 1.0;
}

template <>
void lagrange_basis<2u>(const double pt, double N[]) {
  N[0] = 0.5 * (1.0 - pt);
  N[1] = 0.5 * (pt - 1.0);
}

template <>
void lagrange_basis<3u>(const double pt, double N[]) {
  N[0] = -0.5 * pt * (1.0 - pt);
  N[1] = (1.0 - pt) * (1.0 + pt);
  N[2] = 0.5 * (1.0 + pt) * pt;
}

template <>
void lagrange_basis<4u>(const double pt, double N[]) {
  N[0] = -(2.0 / 3.0) * (0.5 + pt) * (0.5 - pt) * (1.0 - pt);
  N[1] = (4.0 / 3.0) * (1.0 + pt) * (0.5 - pt) * (1.0 - pt);
  N[2] = (4.0 / 3.0) * (1.0 + pt) * (0.5 + pt) * (1.0 - pt);
  N[3] = -(2.0 / 3.0) * (1.0 + pt) * (0.5 + pt) * (0.5 - pt);
}

template <>
void lagrange_basis<1u>(const double pt, double N[], double Nx[]) {
  N[0] = 1.0;

  Nx[0] = 0.0;
}

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
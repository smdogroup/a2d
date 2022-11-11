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
void lagrange_basis(const double pt, double N[]) {}

template <index_t order>
void lagrange_basis(const double pt, double N[], double Nx[]) {}

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
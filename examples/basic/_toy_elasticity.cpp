#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "a2d.h"

using namespace A2D;

typedef index_t I;
typedef double T;

struct ArgParser {
  ArgParser(int argc, char* argv[], const char exe_name[]) {
    Timer t("ArgParser::ArgParser()");
    // save args
    std::vector<std::string> args(argv + 1, argv + argc);

    // Loop over args
    for (auto it = args.begin(); it != args.end(); it++) {
      if ((*it).find("-h") != std::string::npos ||
          (*it).find("-help") != std::string::npos) {
        std::printf(
            "Usage: ./%s [-nx nx] [-ny ny] [-nz nz] [-lx lx] [-ly ly] "
            "[-lz "
            "lz]\n",
            exe_name);
        std::printf("       -nx (int)    number of elements in x direction\n");
        std::printf("       -ny (int)    number of elements in y direction\n");
        std::printf("       -nz (int)    number of elements in z direction\n");
        std::printf("       -lx (double) domain length in x direction\n");
        std::printf("       -ly (double) domain length in y direction\n");
        std::printf("       -lz (double) domain length in z direction\n");
        exit(0);
      } else if ((*it).find("-nx") != std::string::npos) {
        it++;
        nx = std::stoi(*it);
      } else if ((*it).find("-ny") != std::string::npos) {
        it++;
        ny = std::stoi(*it);
      } else if ((*it).find("-nz") != std::string::npos) {
        it++;
        nz = std::stoi(*it);
      } else if ((*it).find("-lx") != std::string::npos) {
        it++;
        lx = std::stod(*it);
      } else if ((*it).find("-ly") != std::string::npos) {
        it++;
        ly = std::stod(*it);
      } else if ((*it).find("-lz") != std::string::npos) {
        it++;
        lz = std::stod(*it);
      } else {
        std::printf("Unknown argument: %s\n", (*it).c_str());
        exit(-1);
      }
    }

    // Print summary
    index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
    index_t nelems = nx * ny * nz;
    printf(
        "execution configuration\nnx = %d, ny = %d, nz = %d, lx = %.1f, ly = "
        "%.1f, lz = %.1f\n",
        nx, ny, nz, lx, ly, lz);
    printf("nnodes: %d\n", nnodes);
    printf("nelems: %d\n", nelems);
  }

  // Default values
  index_t nx = 64;
  index_t ny = 64;
  index_t nz = 64;
  double lx = 1.0;
  double ly = 1.0;
  double lz = 1.0;
};

void main_body(int argc, char* argv[]) {
  ArgParser p(argc, argv, "toy_elasticity");

  // Define problem dimension
  static const int SPATIAL_DIM = 3;
  using Basis = BasisOps<SPATIAL_DIM, HexTriLinearBasisFunc, Hex8ptQuadrature>;
  using ElasticityPDE = ElasticityPDEInfo<SPATIAL_DIM, I, T>;

  // Set mesher
  index_t nx = p.nx;
  index_t ny = p.ny;
  index_t nz = p.nz;
  double lx = p.lx;
  double ly = p.ly;
  double lz = p.lz;
  MesherBrick3D mesher(nx, ny, nz, lx, ly, lz);

  // Set PDE model and element
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const index_t nelems = nx * ny * nz;
  const index_t nbcs = (ny + 1) * (nz + 1);
  auto model = std::make_shared<FEModel<I, T, ElasticityPDE>>(nnodes, nbcs);
  auto element = std::make_shared<LinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  // Set the boundary conditions
  auto bcs = model->get_bcs();
  mesher.set_bcs(bcs);

  // Set the connectivity
  auto conn = element->get_conn();

  // Set the node locations
  auto X = model->get_nodes();
  mesher.set_X_conn(X, conn);

  // Set the node locations - Note: This must be done after setting the
  // connectivity!
  model->init();

  // Set the element
  T q = 5.0, E = 70e3, nu = 0.3;
  T density = 1.0, design_stress = 1e3;
  auto constitutive = std::make_shared<TopoIsoConstitutive<I, T, Basis>>(
      element, q, E, nu, density, design_stress);
  model->add_constitutive(constitutive);

  // Create the design vector
  auto x = std::make_shared<A2D::MultiArrayNew<T* [1]>>("x", model->nnodes);

  // Set the design variable values
  mesher.set_dv(*x);
  model->set_design_vars(x);

  // Set up the stress functional
  auto functional = std::make_shared<Functional<I, T, ElasticityPDE>>();
  auto agg_functional =
      std::make_shared<TopoVonMisesAggregation<I, T, Basis>>(constitutive);
  functional->add_functional(agg_functional);

  // Compute the Jacobian matrix
  auto J = model->new_matrix();
  model->jacobian(J);

  int num_levels = 3;
  double omega = 0.6667;
  double epsilon = 0.01;
  bool print_info = true;
  auto amg = model->new_amg(num_levels, omega, epsilon, J, print_info);

  // Set the residuals and apply the boundary conditions
  auto solution = model->new_solution();
  auto residual = model->new_solution();

  mesher.set_force(model, residual);

  // Compute the solution
  index_t monitor = 10;
  index_t max_iters = 100;
  A2D::BLAS::zero(*solution);
  amg->cg(*residual, *solution, monitor, max_iters);

  Timer t_vtk("Generate vtk");
  {
    ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())> vtk(conn,
                                                                           X);
    vtk.write_mesh();
    vtk.write_sol("x", *x, 0);
    vtk.write_sol("ux", *solution, 0);
    vtk.write_sol("uy", *solution, 1);
    vtk.write_sol("uz", *solution, 2);
  }
  return;
}

int main(int argc, char* argv[]) {
  Kokkos::initialize();

  Timer main_timer("main");
  { main_body(argc, argv); }

  Kokkos::finalize();
  return (0);
}

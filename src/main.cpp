#include "util.hpp"
#include "mesh.hpp"
#include <cstddef>
#include <cstdint>

std::int32_t main (std::int32_t argc, char *argv[])
{
    using namespace kanba;
    Pair<std::size_t> mesh_size = Pair<std::size_t>(10, 10);
    Pair<double> mesh_length = Pair<double>(1.0, 1.0);
    double delta_t = 0.01;
    Mesh mesh = Mesh(mesh_size, mesh_length, delta_t);
    SquareCylinder sc = SquareCylinder(Pair<std::size_t>(4, 4), Pair<std::size_t>(3, 3));
    mesh.add_square_cylinder(sc);
    mesh.set_uniform_flow(Pair<double>(1.0, 0.0));
    mesh.process_outer_boudary_condition();
    mesh.process_wall_boudary_condition();
    mesh.process_HSMAC(100, 0.5);
    return 0;
}

#include "config.hpp"
#include <cstdio>
#include <print>

//      ##### Header do gerenowania siatki elementów skończonych #####

// Funkcja do tworzenia prostej siatki regularnej
auto makeMesh(const lstr::MpiComm& comm)
{
    const auto node_dist_x =
        lstr::util::linspaceArray< number_of_elements_x >(-X_max, X_max + X_d); // Ilosc podzialow po kierunku X
    const auto node_dist_y =
        lstr::util::linspaceArray< number_of_elements_y >(-L / 2.0, L / 2.0); // Liczba podizalow po kierunku Y

    const auto mesh_generator = [&] {
        return lstr::mesh::makeSquareMesh(node_dist_x, node_dist_y);
    };

    return generateAndDistributeMesh(comm, mesh_generator, L3STER_WRAP_CTVAL(mesh_order));
}
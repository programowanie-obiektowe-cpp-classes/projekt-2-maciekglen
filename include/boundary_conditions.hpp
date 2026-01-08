#include "config.hpp"
#include <array>

//      ##### Warunki brzegowe Neumanna #####

// Neuman na scianach pionowych
auto makeKernellVwall()
{
    constexpr auto kernel_params_vwall =
        lstr::KernelParams{.dimension = n_dimensions, .n_equations = 2, .n_unknowns = 2};
    return lstr::wrapBoundaryEquationKernel< kernel_params_vwall >([&](const auto& in, auto& out) {
        const double nx        = in.normal[0];
        const double ny        = in.normal[1];
        auto& [operators, rhs] = out;
        auto& [A0, A1, A2]     = operators;

        // Pochodna po fi, v, eta
        A1(0, 0) = nx; // dfi_dx = 0
        A1(1, 1) = nx; // deta_dx = 0
    });
}

// Neuman na scianach poziomych
auto makeKernellHwall()
{
    constexpr auto kernel_params_hwall =
        lstr::KernelParams{.dimension = n_dimensions, .n_equations = 2, .n_unknowns = 2};
    return lstr::wrapBoundaryEquationKernel< kernel_params_hwall >([&](const auto& in, auto& out) {
        const double nx        = in.normal[0];
        const double ny        = in.normal[1];
        auto& [operators, rhs] = out;
        auto& [A0, A1, A2]     = operators;

        // Pochodna po fi, u, eta
        A2(0, 0) = ny; // dfi_dy = 0
        A2(1, 1) = ny; // deta_dy = 0
    });
}

// Neuman na scianie cylindrycznej
auto makeKernellCwall()
{
    constexpr auto kernel_params_cwall =
        lstr::KernelParams{.dimension = n_dimensions, .n_equations = 4, .n_unknowns = 4};
    return lstr::wrapBoundaryEquationKernel< kernel_params_cwall >([&](const auto& in, auto& out) {
        const double nx        = in.normal[0];
        const double ny        = in.normal[1];
        auto& [operators, rhs] = out;
        auto& [A0, A1, A2]     = operators;

        // Pochodna po fi
        A1(0, 0) = nx;
        A2(0, 0) = ny;

        // Pochodna po eta
        A1(3, 3) = nx;
        A2(3, 3) = ny;
    });
}

//      ##### Warunki brzegowe Dirichleta #####

lstr::BCDefinition< n_unknowns > definitionDirichlet()
{
    auto bc_def = lstr::BCDefinition< n_unknowns >{};
    bc_def.defineDirichlet({right}, {U, KSIX});       // Warunek sciany na prawym brzegu
    bc_def.defineDirichlet({top, bottom}, {V, KSIY}); // Warunek sciany na gornym i dolnym brzegu
    bc_def.defineDirichlet({left}, {ETA, V});         // Warunek wlotu
    return bc_def;
}

void makeDirichlet(auto&        algebraic_system,
                   const int    Parameter1,
                   const int    Parameter2,
                   const int    place1,
                   const int    place2,
                   const double value1,
                   const double value2)

{
    const auto bc_inds_vwall  = std::array{Parameter1, Parameter2}; // Jakie zmienne chcemy ustalic
    const auto bc_value_vwall = std::array{value1, value2};         // Jaka wartosc chcemy przypisac

    algebraic_system.setDirichletBCValues(bc_value_vwall, {place1, place2}, bc_inds_vwall);
}

void makeDirichlet(auto&        algebraic_system,
                   const int    Parameter1,
                   const int    Parameter2,
                   const int    place1,
                   const int    place2,
                   const int    place3,
                   const double value1,
                   const double value2)

{
    const auto bc_inds_vwall  = std::array{Parameter1, Parameter2}; // Jakie zmienne chcemy ustalic
    const auto bc_value_vwall = std::array{value1, value2};         // Jaka wartosc chcemy przypisac

    algebraic_system.setDirichletBCValues(bc_value_vwall, {place1, place2, place3}, bc_inds_vwall);
}

void makeDirichlet(auto&        algebraic_system,
                   const int    Parameter1,
                   const int    Parameter2,
                   const int    place,
                   const double value1,
                   const double value2)

{
    const auto bc_inds_vwall  = std::array{Parameter1, Parameter2}; // Jakie zmienne chcemy ustalic
    const auto bc_value_vwall = std::array{value1, value2};         // Jaka wartosc chcemy przypisac

    algebraic_system.setDirichletBCValues(bc_value_vwall, {place}, bc_inds_vwall);
}

void makeDirichlet(auto& algebraic_system, const int Parameter, const int place, const double value)

{
    const auto bc_inds_vwall  = std::array{Parameter}; // Jakie zmienne chcemy ustalic
    const auto bc_value_vwall = std::array{value};     // Jaka wartosc chcemy przypisac

    algebraic_system.setDirichletBCValues(bc_value_vwall, {place}, bc_inds_vwall);
}

void makeDirichlet(auto& algebraic_system, const int Parameter, const int place1, const int place2, const double value)

{
    const auto bc_inds_vwall  = std::array{Parameter}; // Jakie zmienne chcemy ustalic
    const auto bc_value_vwall = std::array{value};     // Jaka wartosc chcemy przypisac

    algebraic_system.setDirichletBCValues(bc_value_vwall, {place1, place2}, bc_inds_vwall);
}
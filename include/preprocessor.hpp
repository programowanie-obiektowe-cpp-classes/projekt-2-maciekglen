#include "boundary_conditions.hpp"
#include "config.hpp"
#include "initial_conditions.hpp"
#include "kernel.hpp"
#include "mesh.hpp"
#include "solution.hpp"
#include <array>
#include <string>

// Klasa przechowujaca indeksy obliczanyhc wartosci
class Indices
{
public:
    Indices(int n_)
        : fi{FI + n_ * n_unknowns},
          u{U + n_ * n_unknowns},
          v{V + n_ * n_unknowns},
          eta{ETA + n_ * n_unknowns},
          ksix{KSIX + n_ * n_unknowns},
          ksiy{KSIY + n_ * n_unknowns}
    {}
    const std::array< unsigned, 1 > fi;
    const std::array< unsigned, 1 > u;
    const std::array< unsigned, 1 > v;
    const std::array< unsigned, 1 > eta;
    const std::array< unsigned, 1 > ksix;
    const std::array< unsigned, 1 > ksiy;

    decltype(util::concatArrays(fi, u, v, eta, ksix, ksiy)) results_inds =
        util::concatArrays(fi, u, v, eta, ksix, ksiy);
};

class Preprocessor
{
public:
    // Utworzenie prostej siatki regularnej
    Preprocessor(const std::shared_ptr< lstr::MpiComm > comm_)
        : comm(comm_),
          mesh(makeMesh(*comm)),
          bc_def(definitionDirichlet()),
          kernel(makeKernel()),
          kernel_error(makeErrorKernel()),
          kernel_vwall(makeKernellVwall()),
          kernel_hwall(makeKernellHwall())
    {
        const auto     problem_def = ProblemDefinition< n_unknowns >{{domain}};
        constexpr auto algsys_opts = AlgebraicSystemParams{.cond_policy = lstr::CondensationPolicy::ElementBoundary};

        auto algebraic_system = makeAlgebraicSystem(comm, mesh, problem_def, bc_def, L3STER_WRAP_CTVAL(algsys_opts));

        // Rzad kwadratury
        constexpr auto asm_opts     = AssemblyOptions{.value_order = 1, .derivative_order = 3};
        const auto     asmopt_ctval = L3STER_WRAP_CTVAL(asm_opts);
        const auto     solver       = Klu2{};
        auto           exporter     = PvtuExporter{comm, *mesh};

        Indices present(0);
        Indices old(1);
        Indices old_2(2);

        const auto initial_cond = makeInitialCond(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        makeDirichlet(algebraic_system, U, KSIX, right, 0.0, 0.0);
        makeDirichlet(algebraic_system, V, KSIY, top, bottom, left, 0.0, 0.0);

        auto solution_manager = SolutionManager{*mesh, (time_order + 1) * n_unknowns};

        solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, present.results_inds, {});
        solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, old.results_inds, {});
        solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, old_2.results_inds, {});

        if (time_export == 0)
        {
            // Eksport wartosci w chwili poczatkowej
            auto export_def = ExportDefinition{"results/bous_2D_0000.pvtu"};
            export_def.defineField("Potencjal", old.fi);
            export_def.defineField("PredkoscX", old.u);
            export_def.defineField("WysokoscFali", old.eta);
            exporter.exportSolution(export_def, solution_manager);
        }

        // Naglowek
        std::println("{:^12}|{:^8}|{:^10}|{:^10}|{:^10}|{:^10}|{:^10}|{:^10}|{:^10}",
                     "Krok czasowy",
                     "Iteracja",
                     "Zbieznosc",
                     "Blad 1",
                     "Blad 2",
                     "Blad 3",
                     "Blad 4",
                     "Blad 5",
                     "Blad 6");
#if (time_order == 2)

        Solution solution(present,
                          old,
                          old_2,
                          mesh,
                          comm,
                          solution_manager,
                          kernel_error,
                          algebraic_system,
                          kernel,
                          kernel_vwall,
                          kernel_hwall,
                          asmopt_ctval,
                          solver,
                          exporter);

#endif

#if (time_order == 3)
        Indices old_3(3);
        solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, old_3.results_inds, {});
        Solution solution(present,
                          old,
                          old_2,
                          old_3,
                          mesh,
                          comm,
                          solution_manager,
                          kernel_error,
                          algebraic_system,
                          kernel,
                          kernel_vwall,
                          kernel_hwall,
                          asmopt_ctval,
                          solver,
                          exporter);

#endif
    }

    // Odczytanie siatki z pliku
    Preprocessor(const std::shared_ptr< lstr::MpiComm > comm_, std::string file)
        : comm(comm_),
          bc_def(definitionDirichlet()),
          kernel(makeKernel()),
          kernel_error(makeErrorKernel()),
          kernel_vwall(makeKernellVwall()),
          kernel_hwall(makeKernellHwall())
    {

        if (!std::filesystem::exists(file))
        {
            throw std::runtime_error("Podany plik z siatka nie istnieje: " + file);
        }

        const auto problem_def = ProblemDefinition< n_unknowns >{{domain}};
        mesh                   = readAndDistributeMesh< mesh_order >(
            *comm, file, mesh::gmsh_tag, {left, bottom, top, right}, {}, problem_def);

        constexpr auto algsys_opts = AlgebraicSystemParams{.cond_policy = lstr::CondensationPolicy::ElementBoundary};

        auto algebraic_system = makeAlgebraicSystem(comm, mesh, problem_def, bc_def, L3STER_WRAP_CTVAL(algsys_opts));

        // Rzad kwadratury
        constexpr auto asm_opts     = AssemblyOptions{.value_order = 1, .derivative_order = 3};
        const auto     asmopt_ctval = L3STER_WRAP_CTVAL(asm_opts);
        const auto     solver       = Klu2{};
        auto           exporter     = PvtuExporter{comm, *mesh};

        Indices present(0);
        Indices old(1);
        Indices old_2(2);

        const auto initial_cond = makeInitialCond(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        makeDirichlet(algebraic_system, U, KSIX, right, 0.0, 0.0);
        makeDirichlet(algebraic_system, V, KSIY, top, bottom, left, 0.0, 0.0);

        auto solution_manager = SolutionManager{*mesh, (time_order + 1) * n_unknowns};

        solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, present.results_inds, {});
        solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, old.results_inds, {});
        solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, old_2.results_inds, {});

        if (time_export == 0)
        {
            // Eksport wartosci w chwili poczatkowej
            auto export_def = ExportDefinition{"results/bous_2D_0000.pvtu"};
            export_def.defineField("Potencjal", old.fi);
            export_def.defineField("PredkoscX", old.u);
            export_def.defineField("WysokoscFali", old.eta);
            exporter.exportSolution(export_def, solution_manager);
        }

        // Naglowek
        std::println("{:^12}|{:^8}|{:^10}|{:^10}|{:^10}|{:^10}|{:^10}|{:^10}|{:^10}",
                     "Krok czasowy",
                     "Iteracja",
                     "Zbieznosc",
                     "Blad 1",
                     "Blad 2",
                     "Blad 3",
                     "Blad 4",
                     "Blad 5",
                     "Blad 6");
#if (time_order == 2)

        Solution solution(present,
                          old,
                          old_2,
                          mesh,
                          comm,
                          solution_manager,
                          kernel_error,
                          algebraic_system,
                          kernel,
                          kernel_vwall,
                          kernel_hwall,
                          asmopt_ctval,
                          solver,
                          exporter);

#endif

#if (time_order == 3)

        if (time_order == 3)
        {
            Indices old_3(3);
            solution_manager.setFields(*comm, *mesh, initial_cond, {domain}, old_3.results_inds, {});
            Solution solution(present,
                              old,
                              old_2,
                              old_3,
                              mesh,
                              comm,
                              solution_manager,
                              kernel_error,
                              algebraic_system,
                              kernel,
                              kernel_vwall,
                              kernel_hwall,
                              asmopt_ctval,
                              solver,
                              exporter);
        }
#endif
    }

    const std::shared_ptr< lstr::MpiComm > comm;

private:
    double                     a;
    decltype(makeMesh(*comm))  mesh;
    BCDefinition< n_unknowns > bc_def;

    decltype(makeKernel())       kernel;
    decltype(makeErrorKernel())  kernel_error;
    decltype(makeKernellVwall()) kernel_vwall;
    decltype(makeKernellHwall()) kernel_hwall;
};
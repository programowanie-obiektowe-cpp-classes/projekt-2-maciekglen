#include "config.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <utility>

using namespace lstr;

class Solution
{
public:
#if (time_order == 2)
    // BDF2
    Solution(auto& present,
             auto& old,
             auto& old_2,
             auto& mesh,
             auto& comm,
             auto& solution_manager,
             auto& error_kernel,
             auto& algebraic_system,
             auto& kernel,
             auto& kernel_vwall,
             auto& kernel_hwall,
             auto  asmopt_ctval,
             auto  solver,
             auto& exporter)
    {

        constexpr auto vwall_dofs  = std::array< size_t, 2 >{FI, ETA};
        constexpr auto vwall_ctval = L3STER_WRAP_CTVAL(vwall_dofs);
        constexpr auto hwall_dofs  = std::array< size_t, 2 >{FI, ETA};
        constexpr auto hwall_ctval = L3STER_WRAP_CTVAL(hwall_dofs);

        auto inds_refs = std::array< std::reference_wrapper< decltype(present.results_inds) >, 2 >{old_2.results_inds,
                                                                                                   old.results_inds};

        const auto inds_for_access = util::concatArrays(old_2.results_inds, old.results_inds, present.results_inds);

        auto val_getter = solution_manager.getFieldAccess(inds_for_access);
        auto error      = computeNormL2(*comm, error_kernel, *mesh, {domain}, val_getter);

        auto start_time_solution = std::chrono::high_resolution_clock::now(); // Początek pomiaru czasu calej symulacji

        for (int time_step = 1; time_step <= time_steps; ++time_step)
        {

            prev_err = 0.1;
            const auto inds_for_access =
                util::concatArrays(inds_refs[0].get(), inds_refs[1].get(), present.results_inds);
            for (int n_iter = 1; n_iter <= max_iter; ++n_iter)
            {
                // Dostep do starych wartosci
                val_getter = solution_manager.getFieldAccess(inds_for_access);

                // Poczatek robienia ukladu
                algebraic_system.beginAssembly();

                t = time_step * dt; // Aktualny czas

                double height = H * std::cos(omega * t - pi / 2);

                const auto bc_value_inlet = std::array{height}; // Okreslenie wysokosci fali na wlocie
                algebraic_system.setDirichletBCValues(bc_value_inlet, {left}, std::array{ETA});

                algebraic_system.assembleProblem(kernel_vwall, {right}, {}, vwall_ctval, asmopt_ctval);
                algebraic_system.assembleProblem(kernel_hwall, {top, bottom}, {}, hwall_ctval, asmopt_ctval);
                algebraic_system.assembleProblem(kernel, {domain}, val_getter, {}, asmopt_ctval);

                // Koniec robienia ukladu
                algebraic_system.endAssembly();

                // Rozwiazanie
                algebraic_system.solve(solver);

                // Aktualizacja starych rozwiazan
                algebraic_system.updateSolution(present.results_inds, solution_manager, present.results_inds);

                // Obliczenie zbieznosci
                error      = computeNormL2(*comm, error_kernel, *mesh, {domain}, val_getter);
                error_norm = error.norm();

                if (time_step % verbosity == 0)
                {
                    er0 = error[0];
                    er1 = error[1];
                    er2 = error[2];
                    er3 = error[3];
                    er4 = error[4];
                    er5 = error[5];

                    std::println("{:^12}|{:^8}|{:^10.2e}|{:^10.2e}|{:^10.2e}|{:^10.2e}|{:^"
                                 "10.2e}|{:^10.2e}|{:^10.2e}",
                                 time_step,
                                 n_iter,
                                 error_norm,
                                 er0,
                                 er1,
                                 er2,
                                 er3,
                                 er4,
                                 er5); // Wypisanie informacji o symulacji}
                }

                if (time_step + 1 % verbosity * 10 == 0)
                {
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
                }

                // Ustalenie poziomu zbieznosci
                if (error_norm < error_max or std::fabs(prev_err - error_norm) / prev_err < .001)

                    break;

                prev_err = error_norm;
            }

            // Aktualizacja wartosci z poprzedniego kroku
            algebraic_system.updateSolution(present.results_inds, solution_manager, inds_refs[0].get());

            std::ranges::rotate(inds_refs, inds_refs.begin() + 1);

            if (time_step >= time_export)
            {

                // Eksport obecnej wartosci do ParaView
                auto export_def = ExportDefinition{std::format("results/bous_2D_{:04}.pvtu", time_step)};

                // Eksport poszczegolnych wartosci
                export_def.defineField("Potencjal", present.fi);
                export_def.defineField("PredkoscX", present.u);
                export_def.defineField("WysokoscFali", present.eta);
                exporter.exportSolution(export_def, solution_manager);
            }

            // Ustalenie poziomu zbieznosci
            if (prev_err > 1)
            {
                throw std::runtime_error("BLAD - norma bledu przekroczyla 1 - brak zbieznosci");
            }
        }

        auto end_time_solution = std::chrono::high_resolution_clock::now(); // Koniec pomiaru czasu dla kroku czasowego
        std::chrono::duration< double > duration_time_step = end_time_solution - start_time_solution;
        if (t >= time_steps * dt)
        {
            std::println("Symulacja zakononcza sukcesem!\nCzas trwania symulacji: {}", duration_time_step);
        }
    }
#endif

#if (time_order == 3)
    // BDF3
    Solution(auto& present,
             auto& old,
             auto& old_2,
             auto& old_3,
             auto& mesh,
             auto& comm,
             auto& solution_manager,
             auto& error_kernel,
             auto& algebraic_system,
             auto& kernel,
             auto& kernel_vwall,
             auto& kernel_hwall,
             auto  asmopt_ctval,
             auto  solver,
             auto& exporter)
    {

        constexpr auto vwall_dofs  = std::array< size_t, 2 >{FI, ETA};
        constexpr auto vwall_ctval = L3STER_WRAP_CTVAL(vwall_dofs);
        constexpr auto hwall_dofs  = std::array< size_t, 2 >{FI, ETA};
        constexpr auto hwall_ctval = L3STER_WRAP_CTVAL(hwall_dofs);

        auto inds_refs = std::array< std::reference_wrapper< decltype(present.results_inds) >, 3 >{
            old_3.results_inds, old_2.results_inds, old.results_inds};

        const auto inds_for_access =
            util::concatArrays(old_3.results_inds, old_2.results_inds, old.results_inds, present.results_inds);

        auto val_getter = solution_manager.getFieldAccess(inds_for_access);
        auto error      = computeNormL2(*comm, error_kernel, *mesh, {domain}, val_getter);

        auto start_time_solution = std::chrono::high_resolution_clock::now(); // Początek pomiaru czasu calej symulacji

        for (int time_step = 1; time_step <= time_steps; ++time_step)
        {
            prev_err = 0.1;
            const auto inds_for_access =
                util::concatArrays(inds_refs[0].get(), inds_refs[1].get(), inds_refs[2].get(), present.results_inds);
            for (int n_iter = 1; n_iter <= max_iter; ++n_iter)
            {

                // Dostep do starych wartosci
                val_getter = solution_manager.getFieldAccess(inds_for_access);

                // Poczatek robienia ukladu
                algebraic_system.beginAssembly();

                t = time_step * dt; // Aktualny czas

                double height    = H * std::cos(omega * t - pi / 2);
                double potencjal = (H / k) * omega * std::sin(omega * t - pi / 2);

                const auto bc_value_inlet = std::array{height}; // Okreslenie wysokosci fali na wlocie
                algebraic_system.setDirichletBCValues(bc_value_inlet, {left}, std::array{ETA});

                algebraic_system.assembleProblem(kernel_vwall, {right}, {}, vwall_ctval, asmopt_ctval);
                algebraic_system.assembleProblem(kernel_hwall, {top, bottom}, {}, hwall_ctval, asmopt_ctval);
                algebraic_system.assembleProblem(kernel, {domain}, val_getter, {}, asmopt_ctval);

                // Koniec robienia ukladu
                algebraic_system.endAssembly();

                // Rozwiazanie
                algebraic_system.solve(solver);

                // Aktualizacja starych rozwiazan
                algebraic_system.updateSolution(present.results_inds, solution_manager, present.results_inds);

                // Obliczenie zbieznosci
                error      = computeNormL2(*comm, error_kernel, *mesh, {domain}, val_getter);
                error_norm = error.norm();

                if (time_step % verbosity == 0)
                {
                    er0 = error[0];
                    er1 = error[1];
                    er2 = error[2];
                    er3 = error[3];
                    er4 = error[4];
                    er5 = error[5];

                    std::println("{:^12}|{:^8}|{:^10.2e}|{:^10.2e}|{:^10.2e}|{:^10.2e}|{:^"
                                 "10.2e}|{:^10.2e}|{:^10.2e}",
                                 time_step,
                                 n_iter,
                                 error_norm,
                                 er0,
                                 er1,
                                 er2,
                                 er3,
                                 er4,
                                 er5); // Wypisanie informacji o symulacji}
                }

                if (time_step + 1 % verbosity * 10 == 0)
                {
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
                }

                // Ustalenie poziomu zbieznosci
                if (error_norm < error_max or std::fabs(prev_err - error_norm) / prev_err < .001)
                    break;
                prev_err = error_norm;
            }

            // Aktualizacja wartosci z poprzedniego kroku
            algebraic_system.updateSolution(present.results_inds, solution_manager, inds_refs[0].get());

            std::ranges::rotate(inds_refs, inds_refs.begin() + 1);

            if (time_step >= time_export)
            {

                // Eksport obecnej wartosci do ParaView
                auto export_def = ExportDefinition{std::format("results/bous_2D_{:04}.pvtu", time_step)};

                // Eksport poszczegolnych wartosci
                export_def.defineField("Potencjal", present.fi);
                export_def.defineField("PredkoscX", present.u);
                export_def.defineField("WysokoscFali", present.eta);
                exporter.exportSolution(export_def, solution_manager);
            }

            // Ustalenie poziomu zbieznosci
            if (prev_err > 10)
            {
                throw std::runtime_error("BLAD - norma bledu przekroczyla 1 - brak zbieznosci");
            }
        }

        auto end_time_solution = std::chrono::high_resolution_clock::now(); // Koniec pomiaru czasu dla kroku czasowego
        std::chrono::duration< double > duration_time_step = end_time_solution - start_time_solution;
        if (t >= time_steps * dt)
        {
            std::println("Symulacja zakononcza sukcesem!\nCzas trwania symulacji: {}", duration_time_step);
        }
    }
#endif

private:
    double prev_err = 0.1;               // Zmienna bledu
    double t        = 0.0;               // Zmienna czasu
    double error_norm;                   // Norma bledow
    double er0, er1, er2, er3, er4, er5; // Bledy rownan
};

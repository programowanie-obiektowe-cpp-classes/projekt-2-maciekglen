#ifndef config_hpp
#define config_hpp

#include "l3ster/l3ster.hpp"
#include <cmath>

//      ##### Plik do ustawiania różnych wartości używanych przez program #####

// Parametry programu
constexpr double dt = 0.005; // Krok czasowy

constexpr unsigned time_steps = 8000; // Liczba krokow czasowych
                                      //
constexpr unsigned max_iter = 5;      // Maksymalna liczba iteracji w kroku czasowym

constexpr double element_size = 0.02; // Wielkosc elementow skonczonych

constexpr lstr::el_o_t mesh_order = 3; // Rzad elementow

constexpr unsigned verbosity = 10; // Co ile krokow czasowych ma sie wyswietlac informacja o zbieznosci

constexpr unsigned time_export = 4000; // Od ktorego kroku czasowego nalezy zapisywac wyniki

constexpr double error_max = 1e-7;

#define time_order 3 // Rzad dyskretyzacji czasowej

constexpr bool linear_damping_function = true; // Jesli funkcja tlumiaca ma byc liniowa to true, jak kwadratowa to false

// Stałe

constexpr double L = 0.16; // Szerokosc kanalu

constexpr double X_max = 10.0 * L; // Dlugosc kanalu

constexpr double X_d = 25.0 * L; // Dlugosc warstwy tlumiacej

constexpr double nu_0 = 1e-6; // Lepkosc wody

constexpr double A = 0.01; // Stala zwiekszania lepkosci

constexpr double H = 1e-05; // Wysokosc fali

constexpr double h = 0.01; // Srednia glebokosc

// Pomocnicze

constexpr double g = 9.81; // Przyspieszenie ziemskie

constexpr double pi = M_PI; // Stala pi

constexpr double k = pi / L; // Liczba falowa

const double omega = std::sqrt(g * k * std::tanh(k * h)); // Czestosc kolowa fali

constexpr unsigned number_of_elements_x = (unsigned)((X_max * 2.0 + X_d) / element_size + 1.0);

constexpr unsigned number_of_elements_y = 3;

#if (time_order == 2)
// BDF2 wspolczynniki
constexpr double coef0 = 1.0, coef1 = -4.0 / 3.0, coef2 = 1.0 / 3.0, coeft = 2.0 / 3.0;

#elif (time_order == 3)
// BDF3 wspolczynniki
constexpr double coef0 = 1.0, coef1 = -18.0 / 11.0, coef2 = 9.0 / 11.0, coef3 = -2.0 / 11.0, coeft = 6.0 / 11.0;

#else
#error "Nieobsługiwany rząd dyskretyzacji czasowej"
#endif

// Przypisanie odpowiednich nazw
constexpr unsigned domain = 0, bottom = 1, top = 2, left = 3, right = 4;

// Ustalenie niewiadomych
// FI = Potencjal predkosci
// U = Skladowa x predkosci
// V = Skladowa y predkosci
// ETA = Wysokosc fali
// KSI = Pochodna wysokosci, ma dwie skladowe
// D = Pochodna potencjalu fi
// NU = Lepkosc zalezna od x
constexpr unsigned FI = 0, U = 1, V = 2, ETA = 3, KSIX = 4, KSIY = 5;

constexpr unsigned n_unknowns   = 6; // liczba niewiadomych
constexpr unsigned n_dimensions = 2; // Liczba wymiarow

#endif

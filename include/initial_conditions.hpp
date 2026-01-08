#include "config.hpp"

//      ##### Header do definicji warunkow poczatkowych #####

const auto makeInitialCond(double fi, double u, double v, double eta, double ksix, double ksiy)
{
    constexpr auto params = lstr::KernelParams{.dimension = n_dimensions, .n_equations = n_unknowns};
    return lstr::wrapDomainResidualKernel< params >([&](const auto& in, auto& out) {
        out[FI]   = fi;
        out[U]    = u;
        out[V]    = v;
        out[ETA]  = eta;
        out[KSIX] = ksix;
        out[KSIY] = ksiy;
    });
}
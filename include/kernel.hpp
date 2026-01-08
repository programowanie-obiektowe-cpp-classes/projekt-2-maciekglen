
//      ##### Header do definicji rozwiazywanych rownan #####

//           #### Kernel obliczanych rownan ####
#include "config.hpp"

#if (time_order == 2)
constexpr auto makeKernel()
{

    constexpr auto kernel_params = lstr::KernelParams{.dimension   = n_dimensions,
                                                      .n_equations = n_unknowns,
                                                      .n_unknowns  = n_unknowns,
                                                      .n_fields    = n_unknowns * (time_order + 1)};

    return lstr::wrapDomainEquationKernel< kernel_params >([](const auto& in, auto& out) {
        const auto& [field_vals, field_ders, point] = in;
        const auto& [x_ders, y_ders]                = field_ders; // Pochodne czastkowe

        const auto& [fi_old_2,
                     u_old_2,
                     v_old_2,
                     eta_old_2,
                     ksix_old_2,
                     ksiy_old_2, // 2 wstecz
                     fi_old,
                     u_old,
                     v_old,
                     eta_old,
                     ksix_old,
                     ksiy_old, // 1 wstecz
                     fi,
                     u,
                     v,
                     eta,
                     ksix,
                     ksiy] = field_vals; // Aktualne

        const auto& [dfi_old_2_dx,
                     du_old_2_dx,
                     dv_old_2_dx,
                     deta_old_2_dx,
                     dksix_old_2_dx,
                     dksiy_old_2_dx, // 2 wstecz
                     dfi_old_dx,
                     du_old_dx,
                     dv_old_dx,
                     deta_old_dx,
                     dksix_old_dx,
                     dksiy_old_dx, // 1 wstecz
                     dfi_dx,
                     du_dx,
                     dv_dx,
                     deta_dx,
                     dksix_dx,
                     dksiy_dx] = x_ders; // Aktualne

        const auto& [dfi_old_2_dy,
                     du_old_2_dy,
                     dv_old_2_dy,
                     deta_old_2_dy,
                     dksix_old_2_dy,
                     dksiy_old_2_dy, // 2 wstecz
                     dfi_old_dy,
                     du_old_dy,
                     dv_old_dy,
                     deta_old_dy,
                     dksix_old_dy,
                     dksiy_old_dy, // 1 wstecz
                     dfi_dy,
                     du_dy,
                     dv_dy,
                     deta_dy,
                     dksix_dy,
                     dksiy_dy] = y_ders; // Aktualne

        auto& [operators, rhs] = out;
        auto& [A0, A1, A2]     = operators;

        const auto x = in.point.space.x(); // Polozenie X punktu

        double nu     = nu_0;
        double dnu_dx = 0.0;
        double dnu_dy = 0.0;

        if (x >= X_max)
        {

            if (linear_damping_function == true)
            {
                nu     = A * (x - X_max) + nu_0; // Liniowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A;                      // Analitycznie obliczona pochodna lepkości po x}
            }
            else
            {
                nu     = A * (x - X_max) * (x - X_max) + nu_0; // Kwadratowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A * 2 * (x - X_max);                  // Analitycznie obliczona pochodna lepkości po x
            }
        }

        // Rownanie 0 - ciaglosci

        A0(0, U)   = deta_dx;
        A0(0, V)   = deta_dy;
        A0(0, ETA) = coef0 / (coeft * dt) + du_dx + dv_dy;

        A1(0, U)    = h + eta;
        A1(0, ETA)  = u;
        A1(0, KSIX) = (-2.0 * h * h / 5.0 * coef0) / (coeft * dt);

        A2(0, V)    = h + eta;
        A2(0, ETA)  = v;
        A2(0, KSIY) = (-2.0 * h * h / 5.0 * coef0) / (coeft * dt);

        // Rownanie 1 - pedu

        A0(1, FI)  = coef0 / (coeft * dt);
        A0(1, U)   = u - dnu_dx;
        A0(1, V)   = v - dnu_dy;
        A0(1, ETA) = g;

        A1(1, KSIX) = -g * h * h / 15.0;

        A2(1, KSIY) = -g * h * h / 15.0;

        A1(1, U) = -nu; // Czlon tlumiacy
        A2(1, V) = -nu; // Czlon tlumiacy

        // Rownanie 2 i 3 - definicja u i v
        A1(2, FI) = 1.0;
        A0(2, U)  = -1.0;

        A2(3, FI) = 1.0;
        A0(3, V)  = -1.0;

        // Rownanie 4 i 5 - definicja ksix i ksiy
        A1(4, ETA)  = 1.0;
        A0(4, KSIX) = -1.0;

        A2(5, ETA)  = 1.0;
        A0(5, KSIY) = -1.0;

        // Prawe strony rownan
        rhs[0] = -(coef1 * eta_old + coef2 * eta_old_2) / (coeft * dt) + eta * du_dx + u * deta_dx + eta * dv_dy +
                 v * deta_dy - (-2.0 * h * h / 5.0 * (coef1 * dksix_old_dx + coef2 * dksix_old_2_dx)) / (coeft * dt) -
                 (-2.0 * h * h / 5.0 * (coef1 * dksiy_old_dy + coef2 * dksiy_old_2_dy)) / (coeft * dt);
        rhs[1] = -(coef1 * fi_old + coef2 * fi_old_2) / (coeft * dt) + 0.5 * u * u + 0.5 * v * v;
    });
}
#endif

#if (time_order == 3)
constexpr auto makeKernel()
{

    constexpr auto kernel_params = lstr::KernelParams{.dimension   = n_dimensions,
                                                      .n_equations = n_unknowns,
                                                      .n_unknowns  = n_unknowns,
                                                      .n_fields    = n_unknowns * (time_order + 1)};

    return lstr::wrapDomainEquationKernel< kernel_params >([](const auto& in, auto& out) {
        const auto& [field_vals, field_ders, point] = in;
        const auto& [x_ders, y_ders]                = field_ders; // Pochodne czastkowe

        const auto& [fi_old_3,
                     u_old_3,
                     v_old_3,
                     eta_old_3,
                     ksix_old_3,
                     ksiy_old_3, // 3 wstecz
                     fi_old_2,
                     u_old_2,
                     v_old_2,
                     eta_old_2,
                     ksix_old_2,
                     ksiy_old_2, // 2 wstecz
                     fi_old,
                     u_old,
                     v_old,
                     eta_old,
                     ksix_old,
                     ksiy_old, // 1 wstecz
                     fi,
                     u,
                     v,
                     eta,
                     ksix,
                     ksiy] = field_vals; // Aktualne

        const auto& [dfi_old_3_dx,
                     du_old_3_dx,
                     dv_old_3_dx,
                     deta_old_3_dx,
                     dksix_old_3_dx,
                     dksiy_old_3_dx, // 3 wstecz
                     dfi_old_2_dx,
                     du_old_2_dx,
                     dv_old_2_dx,
                     deta_old_2_dx,
                     dksix_old_2_dx,
                     dksiy_old_2_dx, // 2 wstecz
                     dfi_old_dx,
                     du_old_dx,
                     dv_old_dx,
                     deta_old_dx,
                     dksix_old_dx,
                     dksiy_old_dx, // 1 wstecz
                     dfi_dx,
                     du_dx,
                     dv_dx,
                     deta_dx,
                     dksix_dx,
                     dksiy_dx] = x_ders; // Aktualne

        const auto& [dfi_old_3_dy,
                     du_old_3_dy,
                     dv_old_3_dy,
                     deta_old_3_dy,
                     dksix_old_3_dy,
                     dksiy_old_3_dy, // 3 wstecz
                     dfi_old_2_dy,
                     du_old_2_dy,
                     dv_old_2_dy,
                     deta_old_2_dy,
                     dksix_old_2_dy,
                     dksiy_old_2_dy, // 2 wstecz
                     dfi_old_dy,
                     du_old_dy,
                     dv_old_dy,
                     deta_old_dy,
                     dksix_old_dy,
                     dksiy_old_dy, // 1 wstecz
                     dfi_dy,
                     du_dy,
                     dv_dy,
                     deta_dy,
                     dksix_dy,
                     dksiy_dy] = y_ders; // Aktualne

        auto& [operators, rhs] = out;
        auto& [A0, A1, A2]     = operators;

        const auto x = in.point.space.x(); // Polozenie X punktu

        double nu     = nu_0;
        double dnu_dx = 0.0;
        double dnu_dy = 0.0;

        if (x >= X_max)
        {

            if (linear_damping_function == true)
            {
                nu     = A * (x - X_max) + nu_0; // Liniowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A;                      // Analitycznie obliczona pochodna lepkości po x}
            }
            else
            {
                nu     = A * (x - X_max) * (x - X_max) + nu_0; // Kwadratowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A * 2 * (x - X_max);                  // Analitycznie obliczona pochodna lepkości po x
            }
        }

        // Rownanie 0 - ciaglosci

        A0(0, U)   = deta_dx;
        A0(0, V)   = deta_dy;
        A0(0, ETA) = coef0 / (coeft * dt) + du_dx + dv_dy;

        A1(0, U)    = h + eta;
        A1(0, ETA)  = u;
        A1(0, KSIX) = (-2.0 * h * h / 5.0 * coef0) / (coeft * dt);

        A2(0, V)    = h + eta;
        A2(0, ETA)  = v;
        A2(0, KSIY) = (-2.0 * h * h / 5.0 * coef0) / (coeft * dt);

        // Rownanie 1 - pedu

        A0(1, FI)  = coef0 / (coeft * dt);
        A0(1, U)   = u - dnu_dx;
        A0(1, V)   = v - dnu_dy;
        A0(1, ETA) = g;

        A1(1, KSIX) = -g * h * h / 15.0;

        A2(1, KSIY) = -g * h * h / 15.0;

        A1(1, U) = -nu; // Czlon tlumiacy
        A2(1, V) = -nu; // Czlon tlumiacy

        // Rownanie 2 i 3 - definicja u i v
        A1(2, FI) = 1.0;
        A0(2, U)  = -1.0;

        A2(3, FI) = 1.0;
        A0(3, V)  = -1.0;

        // Rownanie 4 i 5 - definicja ksix i ksiy
        A1(4, ETA)  = 1.0;
        A0(4, KSIX) = -1.0;

        A2(5, ETA)  = 1.0;
        A0(5, KSIY) = -1.0;

        // Prawe strony rownan
        rhs[0] = -(coef1 * eta_old + coef2 * eta_old_2 + coef3 * eta_old_3) / (coeft * dt) + eta * du_dx + u * deta_dx +
                 eta * dv_dy + v * deta_dy -
                 (-2.0 * h * h / 5.0 * (coef1 * dksix_old_dx + coef2 * dksix_old_2_dx + coef3 * dksix_old_3_dx)) /
                     (coeft * dt) -
                 (-2.0 * h * h / 5.0 * (coef1 * dksiy_old_dy + coef2 * dksiy_old_2_dy + coef3 * dksiy_old_3_dy)) /
                     (coeft * dt);
        rhs[1] = -(coef1 * fi_old + coef2 * fi_old_2 + coef3 * fi_old_3) / (coeft * dt) + 0.5 * u * u + 0.5 * v * v;
    });
}
#endif

#if (time_order == 2)
//           #### Kernel obliczania bledu ####
constexpr auto makeErrorKernel()
{
    constexpr auto kernel_params = lstr::KernelParams{.dimension   = n_dimensions,
                                                      .n_equations = n_unknowns,
                                                      .n_unknowns  = n_unknowns,
                                                      .n_fields    = n_unknowns * (time_order + 1)};
    //           ####Kernel obliczania zbieznosci####
    constexpr auto error = [](const auto& in, auto& error) {
        const auto& point            = in.point;
        const auto& vals             = in.field_vals;
        const auto& [x_ders, y_ders] = in.field_ders;

        const auto& [fi_old_2,
                     u_old_2,
                     v_old_2,
                     eta_old_2,
                     ksix_old_2,
                     ksiy_old_2, // 2 wstecz
                     fi_old,
                     u_old,
                     v_old,
                     eta_old,
                     ksix_old,
                     ksiy_old, // 1 wstecz
                     fi,
                     u,
                     v,
                     eta,
                     ksix,
                     ksiy] = vals; // Aktualne

        const auto& [dfi_old_2_dx,
                     du_old_2_dx,
                     dv_old_2_dx,
                     deta_old_2_dx,
                     dksix_old_2_dx,
                     dksiy_old_2_dx, // 2 wstecz
                     dfi_old_dx,
                     du_old_dx,
                     dv_old_dx,
                     deta_old_dx,
                     dksix_old_dx,
                     dksiy_old_dx, // 1 wstecz
                     dfi_dx,
                     du_dx,
                     dv_dx,
                     deta_dx,
                     dksix_dx,
                     dksiy_dx] = x_ders; // Aktualne

        const auto& [dfi_old_2_dy,
                     du_old_2_dy,
                     dv_old_2_dy,
                     deta_old_2_dy,
                     dksix_old_2_dy,
                     dksiy_old_2_dy, // 2 wstecz
                     dfi_old_dy,
                     du_old_dy,
                     dv_old_dy,
                     deta_old_dy,
                     dksix_old_dy,
                     dksiy_old_dy, // 1 wstecz
                     dfi_dy,
                     du_dy,
                     dv_dy,
                     deta_dy,
                     dksix_dy,
                     dksiy_dy] = y_ders; // Aktualne

        const auto x = in.point.space.x(); // Polozenie X punktu

        double           nu     = nu_0;
        double           dnu_dx = 0.0;
        constexpr double dnu_dy = 0.0;

        if (x >= X_max)
        {

            if (linear_damping_function == true)
            {
                nu     = A * (x - X_max) + nu_0; // Liniowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A;                      // Analitycznie obliczona pochodna lepkości po x}
            }
            else
            {
                nu     = A * (x - X_max) * (x - X_max) + nu_0; // Kwadratowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A * 2 * (x - X_max);                  // Analitycznie obliczona pochodna lepkości po x
            }
        }

        error[0] =
            (coef0 * eta + coef1 * eta_old + coef2 * eta_old_2) / (coeft * dt) + (h + eta) * du_dx + u * deta_dx +
            (h + eta) * dv_dy + v * deta_dy +
            (-2.0 * h * h / 5.0 * (coef0 * dksix_dx + coef1 * dksix_old_dx + coef2 * dksix_old_2_dx)) / (coeft * dt) +
            (-2.0 * h * h / 5.0 * (coef0 * dksiy_dy + coef1 * dksiy_old_dy + coef2 * dksiy_old_2_dy)) / (coeft * dt);

        error[1] = (coef0 * fi + coef1 * fi_old + coef2 * fi_old_2) / (coeft * dt) + 0.5 * u * u + 0.5 * v * v +
                   g * eta - nu * du_dx - nu * dv_dy - dnu_dx * u - dnu_dy * v - g * h * h / 15 * dksix_dx -
                   g * h * h / 15 * dksiy_dy;

        error[2] = u - dfi_dx;
        error[3] = v - dfi_dy;

        error[4] = ksix - deta_dx;
        error[5] = ksiy - deta_dy;
    };

    return lstr::wrapDomainResidualKernel< kernel_params >(error);
}
#endif

#if (time_order == 3)
constexpr auto makeErrorKernel()
{
    constexpr auto kernel_params = lstr::KernelParams{.dimension   = n_dimensions,
                                                      .n_equations = n_unknowns,
                                                      .n_unknowns  = n_unknowns,
                                                      .n_fields    = n_unknowns * (time_order + 1)};
    //           ####Kernel obliczania zbieznosci####
    constexpr auto error = [](const auto& in, auto& error) {
        const auto& point            = in.point;
        const auto& vals             = in.field_vals;
        const auto& [x_ders, y_ders] = in.field_ders;

        const auto& [fi_old_3,
                     u_old_3,
                     v_old_3,
                     eta_old_3,
                     ksix_old_3,
                     ksiy_old_3, // 3 wstecz
                     fi_old_2,
                     u_old_2,
                     v_old_2,
                     eta_old_2,
                     ksix_old_2,
                     ksiy_old_2, // 2 wstecz
                     fi_old,
                     u_old,
                     v_old,
                     eta_old,
                     ksix_old,
                     ksiy_old, // 1 wstecz
                     fi,
                     u,
                     v,
                     eta,
                     ksix,
                     ksiy] = vals; // Aktualne

        const auto& [dfi_old_3_dx,
                     du_old_3_dx,
                     dv_old_3_dx,
                     deta_old_3_dx,
                     dksix_old_3_dx,
                     dksiy_old_3_dx, // 3 wstecz
                     dfi_old_2_dx,
                     du_old_2_dx,
                     dv_old_2_dx,
                     deta_old_2_dx,
                     dksix_old_2_dx,
                     dksiy_old_2_dx, // 2 wstecz
                     dfi_old_dx,
                     du_old_dx,
                     dv_old_dx,
                     deta_old_dx,
                     dksix_old_dx,
                     dksiy_old_dx, // 1 wstecz
                     dfi_dx,
                     du_dx,
                     dv_dx,
                     deta_dx,
                     dksix_dx,
                     dksiy_dx] = x_ders; // Aktualne

        const auto& [dfi_old_3_dy,
                     du_old_3_dy,
                     dv_old_3_dy,
                     deta_old_3_dy,
                     dksix_old_3_dy,
                     dksiy_old_3_dy, // 3 wstecz
                     dfi_old_2_dy,
                     du_old_2_dy,
                     dv_old_2_dy,
                     deta_old_2_dy,
                     dksix_old_2_dy,
                     dksiy_old_2_dy, // 2 wstecz
                     dfi_old_dy,
                     du_old_dy,
                     dv_old_dy,
                     deta_old_dy,
                     dksix_old_dy,
                     dksiy_old_dy, // 1 wstecz
                     dfi_dy,
                     du_dy,
                     dv_dy,
                     deta_dy,
                     dksix_dy,
                     dksiy_dy] = y_ders; // Aktualne

        const auto x = in.point.space.x(); // Polozenie X punktu

        double nu     = nu_0;
        double dnu_dx = 0.0;
        double dnu_dy = 0.0;

        if (x >= X_max)
        {
            if (linear_damping_function == true)
            {
                nu     = A * (x - X_max) + nu_0; // Liniowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A;                      // Analitycznie obliczona pochodna lepkości po x}
            }
            else
            {
                nu     = A * (x - X_max) * (x - X_max) + nu_0; // Kwadratowe ustalenie lepkosci w warstwie tlumiacej
                dnu_dx = A * 2 * (x - X_max);                  // Analitycznie obliczona pochodna lepkości po x
            }
        }

        error[0] = (coef0 * eta + coef1 * eta_old + coef2 * eta_old_2 + coef3 * eta_old_3) / (coeft * dt) +
                   (h + eta) * du_dx + u * deta_dx + (h + eta) * dv_dy + v * deta_dy +
                   (-2.0 * h * h / 5.0 *
                    (coef0 * dksix_dx + coef1 * dksix_old_dx + coef2 * dksix_old_2_dx + coef3 * dksix_old_3_dx)) /
                       (coeft * dt) +
                   (-2.0 * h * h / 5.0 *
                    (coef0 * dksiy_dy + coef1 * dksiy_old_dy + coef2 * dksiy_old_2_dy + coef3 * dksiy_old_3_dy)) /
                       (coeft * dt);

        error[1] = (coef0 * fi + coef1 * fi_old + coef2 * fi_old_2 + coef3 * fi_old_3) / (coeft * dt) + 0.5 * u * u +
                   0.5 * v * v + g * eta - nu * du_dx - nu * dv_dy - dnu_dx * u - dnu_dy * v -
                   g * h * h / 15 * dksix_dx - g * h * h / 15 * dksiy_dy;

        error[2] = u - dfi_dx;
        error[3] = v - dfi_dy;

        error[4] = ksix - deta_dx;
        error[5] = ksiy - deta_dy;
    };

    return lstr::wrapDomainResidualKernel< kernel_params >(error);
}
#endif
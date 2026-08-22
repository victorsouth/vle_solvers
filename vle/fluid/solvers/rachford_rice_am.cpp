/// @file rachford_rice_am.cpp
/// @brief Реализация численной схемы RR (Ньютон, бисекция); объявления -
///     в `rachford_rice_am.h`.

#include "../../vle_solvers.h"

namespace vlelib {
;

namespace {
;

/// @brief Предел итераций Ньютона при решении RR.
constexpr int k_newton_max_iters = 20;
/// @brief Предел итераций бисекции при решении RR.
constexpr int k_bisection_max_iters = 30;
/// @brief Допуск по ширине интервала при бисекции RR.
constexpr double k_bisection_width_tol = 1e-7;
/// @brief Порог модуля V в Ньютоне RR по доле пара: если |V| не больше этого
///     значения, в знаменатель относительного критерия по значениям на
///     соседних итерациях подставляется 1.0.
constexpr double k_newton_rel_abs_ref_floor =
    4.0 * std::numeric_limits<double>::epsilon();

/// @brief Относительный допуск приращения в методе Ньютона: macheps^(2/3).
inline double rr_flash_newton_steptol()
{
    return std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0);
}
//*****************************************************************************


/// @brief Решение RR методом Ньютона (steptol, лимит итераций).
/// @param vapor_fraction_start Начальное приближение к V.
/// @param equation Уравнение RR(V).
/// @return (V_last, V_prev, iteration_count): последняя и предыдущая итерации,
///     счётчик шагов Ньютона.
std::tuple<double, double, int> solve_rr_newton(
    double vapor_fraction_start,
    const rr_am_equation_t& equation)
{
    const double steptol = rr_flash_newton_steptol();
    double V2 = vapor_fraction_start;
    double V1 = 2.0 * V2;
    int i = 1;

    for (;;) {
        const double ref = (std::abs(V1) > k_newton_rel_abs_ref_floor)
            ? std::abs(V1) : 1.0;
        if (std::abs(V2 - V1) / ref <= steptol) {
            break;
        }
        V1 = V2;
        const double df = equation.derivative(V1);
        const double f = equation.residual(V1);

        if (!std::isfinite(df) || df == 0.0 || !std::isfinite(f)) {
            break;
        }
        V2 = V1 - f / df;
        ++i;
        if (i > k_newton_max_iters) {
            break;
        }
    }
    return { V2, V1, i };
}
//*****************************************************************************


/// @brief Решение RR бисекцией на [Vl, Vr].
/// @param Vl Левая граница интервала локализации корня.
/// @param Vr Правая граница интервала локализации корня.
/// @param equation Уравнение RR(V).
/// @return (status, V, iterations): status 1 - успех; -1 - нет смены знака на
///     концах; -2 - лимит итераций.
std::tuple<int, double, int> solve_rr_bisection(
    double Vl, double Vr,
    const rr_am_equation_t& equation)
{
    double V1 = Vl;
    double V2 = Vr;
    double F1 = equation.residual(V1);
    double F2 = equation.residual(V2);

    if (F1 == 0.0 && std::isfinite(V1)) {
        return { 1, V1, 0 };
    }
    if (F2 == 0.0 && std::isfinite(V2)) {
        return { 1, V2, 0 };
    }
    if (!std::isfinite(F1) || !std::isfinite(F2) || F1 * F2 > 0.0) {
        return { -1, std::numeric_limits<double>::quiet_NaN(), 0 };
    }

    int num_iter = 0;
    while (std::abs(V2 - V1) > k_bisection_width_tol) {
        ++num_iter;
        if (num_iter > k_bisection_max_iters) {
            return { -2, std::numeric_limits<double>::quiet_NaN(), num_iter };
        }
        const double V3 = 0.5 * (V1 + V2);
        const double F3 = equation.residual(V3);
        if (!std::isfinite(F3)) {
            return { -2, std::numeric_limits<double>::quiet_NaN(), num_iter };
        }
        if (F3 == 0.0) {
            return { 1, V3, num_iter };
        }
        // Сохраняем интервал с разными знаками на концах: переносим пару (V, F)
        // на середину и обновляем тот конец, где знак совпал с F в середине.
        if (F1 * F3 < 0.0) {
            V2 = V3;
            F2 = F3;
        } else {
            V1 = V3;
            F1 = F3;
        }
    }
    return { 1, 0.5 * (V1 + V2), num_iter };
}

} // namespace


rr_am_equation_t::rr_am_equation_t(
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k)
    : feed_mole_fraction_m(feed_mole_fraction)
    , equilibrium_k_m(equilibrium_k)
{
}
//*****************************************************************************


double rr_am_equation_t::residual(double vapor_fraction) const
{
    const Eigen::Index n = feed_mole_fraction_m.size();
    const std::vector<double>& equilibrium_k = equilibrium_k_m;
    if (n <= 0 || static_cast<std::size_t>(n) != equilibrium_k.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double result = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) {
        const double z_i = feed_mole_fraction_m(i);
        const double k_i = equilibrium_k[static_cast<std::size_t>(i)];
        const double km1 = k_i - 1.0;
        const double den = 1.0 + vapor_fraction * km1;
        if (!std::isfinite(den) || den == 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        result += z_i * km1 / den;
    }
    return result;
}
//*****************************************************************************


double rr_am_equation_t::derivative(double vapor_fraction) const
{
    const Eigen::Index n = feed_mole_fraction_m.size();
    const std::vector<double>& equilibrium_k = equilibrium_k_m;
    if (n <= 0 || static_cast<std::size_t>(n) != equilibrium_k.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double result = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) {
        const double z_i = feed_mole_fraction_m(i);
        const double k_i = equilibrium_k[static_cast<std::size_t>(i)];
        const double km1 = k_i - 1.0;
        const double den = 1.0 + vapor_fraction * km1;
        if (!std::isfinite(den) || den == 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const double den2 = den * den;
        result -= z_i * (km1 * km1) / den2;
    }
    return result;
}
//*****************************************************************************


rachford_rice_am_t::rachford_rice_am_t(
    std::pair<double, double> vapor_bounds,
    std::pair<double, double> bisection_bounds,
    double vapor_initial,
    const rr_am_equation_t& equation)
    : vapor_bounds_m(vapor_bounds)
    , bisection_bounds_m(bisection_bounds)
    , vapor_initial_m(vapor_initial)
    , equation_m(equation)
{
}
//*****************************************************************************


rachford_rice_am_result_t rachford_rice_am_t::solve()
{
    const double nan_v = std::numeric_limits<double>::quiet_NaN();
    rachford_rice_am_result_t out{};

    // Решение уравнения Рэчфорда-Райса методом Ньютона
    // (solve_rr_newton: итерации по V, на выходе новое V, предыдущее V, число
    // итераций).
    const auto [rr_vapor_fraction_newton, rr_vapor_fraction_old_unused, newton_iters] =
        solve_rr_newton(vapor_initial_m, equation_m);
    (void)rr_vapor_fraction_old_unused;
    double rr_vapor_fraction = rr_vapor_fraction_newton;

    // Отказ от Ньютона: слишком много итераций или корень вне интервала по
    // паровой доле.
    if (newton_iters > k_newton_max_iters
        || rr_vapor_fraction < vapor_bounds_m.first
        || rr_vapor_fraction > vapor_bounds_m.second) {

        // Бисекция невязки RR на интервале, заданном вызывающим кодом
        // (деградация после Ньютона).
        const auto [bstatus, V_bisect, b_iters] =
            solve_rr_bisection(
                bisection_bounds_m.first,
                bisection_bounds_m.second,
                equation_m);
        if (bstatus == -1) {
            out.ok = false;
            out.vapor_split = nan_v;
            return out;
        }
        if (bstatus == -2 || b_iters > k_bisection_max_iters) {
            out.ok = false;
            out.vapor_split = nan_v;
            return out;
        }
        rr_vapor_fraction = V_bisect;
    }

    // Ограничение паровой мольной доли физическим отрезком [0, 1].
    if (rr_vapor_fraction < 0.0) {
        rr_vapor_fraction = 0.0;
    }
    if (rr_vapor_fraction > 1.0) {
        rr_vapor_fraction = 1.0;
    }

    out.ok = true;
    out.vapor_split = rr_vapor_fraction;
    return out;
}

}

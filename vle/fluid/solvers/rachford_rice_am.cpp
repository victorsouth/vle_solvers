/// @file rachford_rice_am.cpp
/// @brief Реализация RR: границы V, невязка, Ньютон, бисекция; объявления -
///     в `rachford_rice_am.h`. Перенос процедур pt_flash Maple (F_RR,
///     F_RR_Newton, границы V и др.).

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
/// @brief Сторожевые значения при поиске границ V по K_i (RR).
constexpr double k_vb_sentinel_left = -1.0e12;
/// @brief Правая граница-сторож при поиске границ V по K_i (RR).
constexpr double k_vb_sentinel_right = 1.0e12;
/// @brief Относительный зазор от асимптот по паровой доле при поджатии концов
///     интервала RR (бисекция, отладочная сетка по V): нижняя граница сдвига
///     `v_min + |v_min| * eps`, верхняя `v_max * (1 - eps)`.
constexpr double k_bisection_inset_rel = 1e-12;

/// @brief Вычисляет величину V_i = 1/(1 - K_i), используемую при анализе
///     интервала решения уравнения Рейчфорда-Райса.
/// @param equilibrium_k Константа равновесия K_i.
/// @return Значение V_i; при K_i = 1 возвращает NaN, поскольку выражение
///     становится сингулярным.
double v_from_equilibrium_k(double equilibrium_k)
{
    const double den = 1.0 - equilibrium_k;
    if (den == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return 1.0 / den;
}
//*****************************************************************************


/// @brief Формирует массив V_i для всех компонентов.
/// @param k_values Вектор констант равновесия K_i.
/// @param v_out Выходной массив V_i.
/// @return false, если K_i или вычисленный V_i некорректен (не конечен).
bool fill_v_from_k_row(const std::vector<double>& k_values, std::vector<double>& v_out)
{
    v_out.clear();
    v_out.reserve(k_values.size());

    for (double k : k_values) {
        // Проверка корректности K_i
        if (!std::isfinite(k)) {
            return false;
        }

        // Преобразование K_i -> V_i
        const double v = v_from_equilibrium_k(k);
        if (!std::isfinite(v)) {
            return false;
        }

        v_out.push_back(v);
    }
    return true;
}
//*****************************************************************************


/// @brief Определяет стартовые границы интервала решения RR в общем случае,
/// когда в смеси присутствуют компоненты как с K_i < 1, так и с K_i > 1.
/// @return (Vleft, Vright). Если подходящих V_i нет, возвращаются сторожевые значения.
std::pair<double, double> get_start_vborders_1(const std::vector<double>& k_values)
{
    std::vector<double> v_row;
    if (!fill_v_from_k_row(k_values, v_row)) {
        return { std::numeric_limits<double>::quiet_NaN()
            , std::numeric_limits<double>::quiet_NaN() };
    }

    // максимальное среди отрицательных
    double v_left = k_vb_sentinel_left;
    // минимальное среди положительных
    double v_right = k_vb_sentinel_right;

    // Поиск максимального V_i < 0
    for (double v : v_row) {
        if (v > v_left && v < 0.0) {
            v_left = v;
        }
    }

    // Поиск минимального V_i > 0
    for (double v : v_row) {
        if (v < v_right && v > 0.0) {
            v_right = v;
        }
    }

    return { v_left, v_right };
}
//*****************************************************************************


/// @brief Границы RR, когда все K_i не меньше 1: две наибольшие V_i после сортировки.
/// @param k_values Константы равновесия K_i по компонентам (не менее двух).
/// @return Пара границ или (NaN, NaN), если данных недостаточно или не число.
std::pair<double, double> get_start_vborders_2(const std::vector<double>& k_values)
{
    std::vector<double> v_row;
    if (!fill_v_from_k_row(k_values, v_row) || v_row.size() < 2) {
        return { std::numeric_limits<double>::quiet_NaN()
            , std::numeric_limits<double>::quiet_NaN() };
    }
    std::sort(v_row.begin(), v_row.end());
    const size_t n = v_row.size();
    return { v_row[n - 2], v_row[n - 1] };
}
//*****************************************************************************


/// @brief Границы RR, когда все K_i не больше 1: две наименьшие V_i после сортировки.
/// @param k_values Константы равновесия K_i по компонентам (не менее двух).
/// @return Пара границ или (NaN, NaN), если данных недостаточно или не число.
std::pair<double, double> get_start_vborders_3(const std::vector<double>& k_values)
{
    std::vector<double> v_row;
    if (!fill_v_from_k_row(k_values, v_row) || v_row.size() < 2) {
        return { std::numeric_limits<double>::quiet_NaN()
            , std::numeric_limits<double>::quiet_NaN() };
    }
    std::sort(v_row.begin(), v_row.end());
    return { v_row[0], v_row[1] };
}
//*****************************************************************************


/// @brief Относительный допуск приращения в методе Ньютона: macheps^(2/3).
inline double rr_flash_newton_steptol()
{
    return std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0);
}
//*****************************************************************************


/// @brief Единая проверка z/K на входе: непустой feed_mole_fraction и тот
///     же размер у equilibrium_k. Иначе 0.
inline std::size_t check_feed_k_component_count(
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k) noexcept
{
    const Eigen::Index n = feed_mole_fraction.size();
    return (n > 0 && static_cast<std::size_t>(n) == equilibrium_k.size())
        ? static_cast<std::size_t>(n) : 0;
}
//*****************************************************************************


/// @brief Невязка RR(V). Число компонентов - feed_mole_fraction.size();
///     размер equilibrium_k тот же (гарантирует вызывающий код).
double rr_router_get_F_RR_impl(double vapor_fraction,
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k)
{
    const Eigen::Index n = feed_mole_fraction.size();
    double result = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) {
        const double z_i = feed_mole_fraction(i);
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


/// @brief Производная d/dV невязки RR. Число компонентов -
///     feed_mole_fraction.size(); размер equilibrium_k тот же.
double rr_router_get_DF_RR_impl(double vapor_fraction,
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k)
{
    const Eigen::Index n = feed_mole_fraction.size();
    double result = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) {
        const double z_i = feed_mole_fraction(i);
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


/// @brief Решение RR методом Ньютона (steptol, лимит итераций).
/// @param vapor_fraction_start Начальное приближение к V.
/// @param feed_mole_fraction Мольные доли подачи z_i.
/// @param equilibrium_k Константы равновесия K_i.
/// @return (V_last, V_prev, iteration_count): последняя и предыдущая итерации,
///     счётчик шагов Ньютона.
std::tuple<double, double, int> rr_router_solve_RR_newton(
    double vapor_fraction_start,
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k)
{
    if (check_feed_k_component_count(feed_mole_fraction, equilibrium_k) == 0) {
        return { std::numeric_limits<double>::quiet_NaN()
            , std::numeric_limits<double>::quiet_NaN(), 0 };
    }

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
        const double df = rr_router_get_DF_RR_impl(
            V1, feed_mole_fraction, equilibrium_k);
        const double f = rr_router_get_F_RR_impl(
            V1, feed_mole_fraction, equilibrium_k);

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
/// @param feed_mole_fraction Мольные доли подачи z_i.
/// @param equilibrium_k Константы равновесия K_i.
/// @return (status, V, iterations): status 1 - успех; -1 - нет смены знака на концах;
///     -2 - лимит итераций.
std::tuple<int, double, int> rr_router_solve_RR_bisection(
    double Vl, double Vr,
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k)
{
    if (check_feed_k_component_count(feed_mole_fraction, equilibrium_k) == 0) {
        return { -1, std::numeric_limits<double>::quiet_NaN(), 0 };
    }

    double V1 = Vl;
    double V2 = Vr;
    double F1 = rr_router_get_F_RR_impl(V1, feed_mole_fraction, equilibrium_k);
    double F2 = rr_router_get_F_RR_impl(V2, feed_mole_fraction, equilibrium_k);

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
        const double F3 = rr_router_get_F_RR_impl(V3, feed_mole_fraction, equilibrium_k);
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


std::pair<double, double> rachford_rice_am_t::get_start_vborders(
    const std::vector<double>& k_values)
{
    if (k_values.empty()) {
        return { std::numeric_limits<double>::quiet_NaN()
            , std::numeric_limits<double>::quiet_NaN() };
    }

    // Проверка: все K_i >= 1?
    bool all_k_ge_1 = true;
    for (double k : k_values) {
        if (k < 1.0) {
            all_k_ge_1 = false;
            break;
        }
    }
    if (all_k_ge_1) {
        return get_start_vborders_2(k_values);
    }

    // Проверка: все K_i <= 1?
    bool all_k_le_1 = true;
    for (double k : k_values) {
        if (k > 1.0) {
            all_k_le_1 = false;
            break;
        }
    }
    if (all_k_le_1) {
        return get_start_vborders_3(k_values);
    }

    // Общий случай
    return get_start_vborders_1(k_values);
}
//*****************************************************************************


double rachford_rice_am_t::residual(double vapor_fraction,
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k)
{
    const std::size_t n = check_feed_k_component_count(
        feed_mole_fraction, equilibrium_k);
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return rr_router_get_F_RR_impl(
        vapor_fraction, feed_mole_fraction, equilibrium_k);
}
//*****************************************************************************


double rachford_rice_am_t::derivative(double vapor_fraction,
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    const std::vector<double>& equilibrium_k)
{
    const std::size_t n = check_feed_k_component_count(
        feed_mole_fraction, equilibrium_k);
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return rr_router_get_DF_RR_impl(
        vapor_fraction, feed_mole_fraction, equilibrium_k);
}
//*****************************************************************************


rachford_rice_am_t::rachford_rice_am_t(
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
    std::vector<double> equilibrium_k,
    std::pair<double, double> vapor_bounds,
    double vapor_initial)
    : feed_mole_fraction_m(feed_mole_fraction)
    , equilibrium_k_m(std::move(equilibrium_k))
    , vapor_bounds_m(vapor_bounds)
    , vapor_initial_m(vapor_initial)
{
}
//*****************************************************************************


rachford_rice_am_result_t rachford_rice_am_t::solve()
{
    const double nan_v = std::numeric_limits<double>::quiet_NaN();
    rachford_rice_am_result_t out{};

    // Решение уравнения Рэчфорда-Райса методом Ньютона
    // (rr_router_solve_RR_newton: итерации по V, на выходе новое V,
    // предыдущее V, число итераций).
    const auto [rr_vapor_fraction_newton, rr_vapor_fraction_old_unused, newton_iters] =
        rr_router_solve_RR_newton(
            vapor_initial_m, feed_mole_fraction_m, equilibrium_k_m);
    (void)rr_vapor_fraction_old_unused;
    double rr_vapor_fraction = rr_vapor_fraction_newton;

    // Отказ от Ньютона: слишком много итераций или корень вне интервала по
    // паровой доле.
    if (newton_iters > k_newton_max_iters
        || rr_vapor_fraction < vapor_bounds_m.first
        || rr_vapor_fraction > vapor_bounds_m.second) {

        // Повторное построение интервала по паровой доле и бисекция невязки
        // RR (деградация после Ньютона).
        const std::pair<double, double> b =
            get_start_vborders(equilibrium_k_m);
        const double Vl =
            b.first + std::abs(b.first) * k_bisection_inset_rel;
        const double Vr = b.second * (1.0 - k_bisection_inset_rel);
        const auto [bstatus, V_bisect, b_iters] =
            rr_router_solve_RR_bisection(
                Vl, Vr, feed_mole_fraction_m, equilibrium_k_m);
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

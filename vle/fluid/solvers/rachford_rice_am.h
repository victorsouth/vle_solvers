#pragma once

/// @file rachford_rice_am.h
/// @brief API решения RR при фиксированных K_i (Ньютон с бисекцией).
///     Соответствие алгоритму pt_flash Maple (F_RR, F_RR_Newton и др.).

namespace vlelib {
;

/// @brief Результат решения уравнения Рэчфорда-Райса.
struct rachford_rice_am_result_t {
    /// @brief true при успешной паровой доле.
    bool ok = false;
    /// @brief Паровая мольная доля V.
    double vapor_split = std::numeric_limits<double>::quiet_NaN();
};

/// @brief Решение RR при фиксированных K_i: Ньютон с деградацией к бисекции.
class rachford_rice_am_t {
    /// @brief Мольные доли подачи z_i; ссылка на вектор вызывающего кода,
    ///     объект-носитель должен существовать до уничтожения этого экземпляра
    ///     и завершения `solve`.
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction_m;
    /// @brief Константы равновесия K_i, заданные на вызов `solve`.
    std::vector<double> equilibrium_k_m;
    /// @brief Интервал паровой доли (Vleft, Vright) для проверки результата
    ///     Ньютона перед переходом к бисекции.
    std::pair<double, double> vapor_bounds_m;
    /// @brief Начальное приближение паровой доли V для итераций Ньютона.
    double vapor_initial_m;

public:
    /// @brief Задаёт подачу z_i, K_i, границы V и начальное приближение для
    ///     `solve`.
    /// @param feed_mole_fraction Мольные доли подачи z_i; вектор-носитель не
    ///     уничтожать до уничтожения этого объекта и завершения `solve`.
    /// @param equilibrium_k Константы равновесия K_i.
    /// @param vapor_bounds Пара (Vleft, Vright) интервала по паровой доле.
    /// @param vapor_initial Начальное приближение к V для Ньютона.
    rachford_rice_am_t(Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
        std::vector<double> equilibrium_k,
        std::pair<double, double> vapor_bounds,
        double vapor_initial);

    /// @brief Решает RR на интервале: Ньютон с переходом на бисекцию при
    ///     неудаче.
    /// @return ok=false и NaN в vapor_split при неудаче Ньютона или бисекции;
    ///     ok=true и V на отрезке [0, 1] после ограничения при успехе.
    rachford_rice_am_result_t solve();

    /// @brief Невязка уравнения Рэчфорда-Райса: сумма z_i (K_i-1) / (1 + V
    ///     (K_i-1)).
    /// @param vapor_fraction Мольная доля отгона V (доля пара в подводимой
    ///     смеси).
    /// @param feed_mole_fraction Мольные доли подачи z_i.
    /// @param equilibrium_k Константы равновесия K_i в том же порядке.
    /// @return Значение RR(V); NaN при несогласованных размерах или нуле
    ///     знаменателя слагаемого.
    static double residual(double vapor_fraction,
        Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
        const std::vector<double>& equilibrium_k);

    /// @brief Производная d/dV функции Рэчфорда-Райса.
    /// @param vapor_fraction Мольная доля отгона V.
    /// @param feed_mole_fraction Мольные доли подачи z_i.
    /// @param equilibrium_k Константы равновесия K_i.
    /// @return Производная; NaN при ошибке слагаемого или нуле знаменателя.
    static double derivative(double vapor_fraction,
        Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
        const std::vector<double>& equilibrium_k);

    /// @brief Пара границ интервала для мольной доли пара (Rachford-Rice), из
    ///     K_i.
    /// @param k_values Константы равновесия по компонентам.
    /// @return (Vleft, Vright).
    static std::pair<double, double> get_start_vborders(
        const std::vector<double>& k_values);
};

}

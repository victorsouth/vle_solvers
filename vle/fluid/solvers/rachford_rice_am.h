#pragma once

/// @file rachford_rice_am.h
/// @brief RR при фиксированных K_i: невязка, производная, численный solve
///     (Ньютон с бисекцией). Перенос pt_flash Maple - в @file tdext PR router.

namespace vlelib {
;

/// @brief Результат решения уравнения Рэчфорда-Райса.
struct rachford_rice_am_result_t {
    /// @brief true при успешной паровой доле.
    bool ok = false;
    /// @brief Паровая мольная доля V.
    double vapor_split = std::numeric_limits<double>::quiet_NaN();
};

/// @brief RR(V) при фиксированных подаче z_i и K_i (контракт PR router).
class rr_am_equation_t {
    /// @brief Мольные доли подачи z_i.
    Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction_m;
    /// @brief Константы равновесия K_i; живут до уничтожения объекта.
    const std::vector<double>& equilibrium_k_m;

public:
    /// @brief Связывает подачу и K_i на время вычисления F(V) и F'(V).
    /// @param feed_mole_fraction Мольные доли подачи z_i.
    /// @param equilibrium_k Константы равновесия K_i; живут до уничтожения
    ///     `rr_am_equation_t`.
    rr_am_equation_t(
        Eigen::Ref<const Eigen::VectorXd> feed_mole_fraction,
        const std::vector<double>& equilibrium_k);

    /// @brief Невязка RR(V): сумма z_i (K_i-1) / (1 + V (K_i-1)); та же
    ///     формула, что у `rr_equation`, но обход **всех** i с проверкой
    ///     знаменателя до суммирования (контракт PR router). При z_i=0 вклад
    ///     нулевой, пока знаменатель конечен; на полюсе V=1/(1-K_i) - NaN на
    ///     всей сумме (в `rr_equation` слагаемое с z_i=0 не вычисляется).
    /// @param vapor_fraction Мольная доля отгона V.
    /// @return Значение RR(V); NaN при несогласованных размерах или нуле
    ///     знаменателя любого слагаемого.
    double residual(double vapor_fraction) const;

    /// @brief Производная d/dV невязки RR(V); тот же контракт знаменателя, что у
    ///     `residual`.
    /// @param vapor_fraction Мольная доля отгона V.
    /// @return Производная; NaN при ошибке слагаемого или нуле знаменателя.
    double derivative(double vapor_fraction) const;
};

/// @brief Решение RR: Ньютон с деградацией к бисекции для `rr_am_equation_t`.
class rachford_rice_am_t {
    /// @brief Интервал паровой доли для проверки результата Ньютона.
    std::pair<double, double> vapor_bounds_m;
    /// @brief Интервал бисекции при деградации (концы с поджатием снаружи).
    std::pair<double, double> bisection_bounds_m;
    /// @brief Начальное приближение паровой доли V для итераций Ньютона.
    double vapor_initial_m;
    /// @brief Уравнение RR; живёт до завершения `solve`.
    const rr_am_equation_t& equation_m;

public:
    /// @brief Задаёт границы, начальное V и уравнение для `solve`.
    /// @param vapor_bounds Интервал проверки корня Ньютона (Vleft, Vright).
    /// @param bisection_bounds Интервал бисекции при деградации после Ньютона.
    /// @param vapor_initial Начальное приближение к V для Ньютона.
    /// @param equation Уравнение RR; не уничтожать до завершения `solve`.
    rachford_rice_am_t(std::pair<double, double> vapor_bounds,
        std::pair<double, double> bisection_bounds,
        double vapor_initial,
        const rr_am_equation_t& equation);

    /// @brief Решает RR: Ньютон, при неудаче бисекция на bisection_bounds.
    /// @return ok=false и NaN в vapor_split при неудаче Ньютона или бисекции;
    ///     ok=true и V на отрезке [0, 1] после ограничения при успехе.
    rachford_rice_am_result_t solve();
};

}

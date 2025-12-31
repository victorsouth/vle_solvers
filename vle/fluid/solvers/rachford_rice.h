#pragma once

namespace vlelib {
;

/// @brief Доли компонента в паре и жидкости
struct component_phase_concentrations_t {
    /// @brief Доля компонента в жидкости
    double liquid{std::numeric_limits<double>::quiet_NaN()};
    /// @brief Доля компонента в паре
    double vapor{ std::numeric_limits<double>::quiet_NaN() };
};

/// @brief Результаты расчета уравнения Речфорда-Райса
struct rachford_rice_result_t {
    /// @brief Состав жидкой фазы
    Eigen::VectorXd x;
    /// @brief Состав паровой фазы
    Eigen::VectorXd y;
    /// @brief Мольная доля отгона
    double vapor_split;
};


/// @brief Уравнение Речфорда-Райса
/// https://chemicals.readthedocs.io/chemicals.rachford_rice.html#two-phase-implementations
/// @param split Доля отгона
/// @param molar_fraction Мольный химсостав
/// @param K_values Константы равновесия
/// @return Невязка по уравнению
inline double rr_equation(double split, const Eigen::VectorXd& molar_fraction, const Eigen::VectorXd& K_values)
{
    double result = 0;
    for (int index = 0; index < molar_fraction.size(); ++index) {
        double Ki = K_values(index);
        double zi = molar_fraction(index);
        if (zi != 0) {
            result += zi * (Ki - 1.) / (1. + split * (Ki - 1.));
        }
    }
    return result;
}

/// @brief Первая производная уравнения Речфорда-Райса
/// https://chemicals.readthedocs.io/chemicals.rachford_rice.html#two-phase-implementations
/// @param split 
/// @param molar_fraction 
/// @param K_values 
/// @return 
double rr_derivative(double split, const Eigen::VectorXd& molar_fraction, const Eigen::VectorXd& K_values);

inline std::pair<Eigen::VectorXd, Eigen::VectorXd> rr_calc_compositions(
    double split, const Eigen::VectorXd& molar_fraction, const Eigen::VectorXd& K_values, bool normalize_compositions = false)
{
    Eigen::VectorXd x(molar_fraction.size());
    Eigen::VectorXd y(molar_fraction.size());
    for (int index = 0; index < molar_fraction.size(); ++index) {
        double Ki = K_values(index);
        double zi = molar_fraction(index);
        if (zi != 0) {
            x(index) = zi / (1. + split * (Ki - 1.));
            y(index) = Ki * x(index);
        }
        else {
            x(index) = 0;
            y(index) = 0;
        }

    }
    if (normalize_compositions) {
        normalize_concentration(x);
        normalize_concentration(y);
    }

    return std::make_pair(std::move(x), std::move(y));
}

// TODO: Удалить эту бестолочь
double calc_omega_extrapolation(const fluid_t* fluid, double P, double T);

/// @brief Усовершенствованная реализация Речфорда-Райса
/// Перевод на fixed-солверы
/// Аналитический расчет производной
/// Аккуратная работа с источниками литературы
/// Использование в задачах TV-flash, UV-flash
class rachford_rice2_t : public fixed_system_t<1> {
private:
    /// @brief Параметры потока
    const fluid_t* fluid;
    /// @brief Константы равновесия, предподсчитанные по начальным термобарическим-условиям
    Eigen::VectorXd K_values;
    /// @brief Точка начала кипения
    double pressure;
    /// @brief Точка начала конденсации
    double temperature;
    /// @brief Начальное значение доли отгона
    double vapor_fraction_initial;
public:
    /// @brief Запоминает флюид, термобарические условия. Предподсчитывает K_values
    rachford_rice2_t(const fluid_t* fluid, double pressure, double temperature,
        double vapor_fraction_initial = std::numeric_limits<double>::quiet_NaN());
    /// @brief Уравнение Речфорда-Райса 
    /// @param split Доля газа
    /// @return Невязка
    virtual double residuals(const double& split) override;
    /// @brief Функция возвращает максимальное и минимальное значение констант равновесия для данной смеси
    /// без учета компонентов, имеющих 0-вые концентрации
    /// @return пара [минимальная К, максимальная К]
    std::pair<double, double> get_k_boundaries() const;
    /// @brief Вычисляет границы omega между ближайшими разрывами 
    /// К разрывам подходим с некоторым эпислоном
    std::pair<double, double> get_omega_boundaries() const;
    /// @brief Производная уравнения Речфорда-Райса 
    /// @param split Доля газа
    /// @return Невязка
    virtual double jacobian_dense(const double& split) override;
public:
    /// @brief Формирует результат расчета для заданной доли отгона
    /// @param vapor_split Доля отгона
    /// @return Результат Речфорда-Райса
    rachford_rice_result_t build_result(double vapor_split) const;
    /// @brief Решение уравнения Речфорда-Райса
    fixed_solver_result_t<1> solve(fixed_solver_result_analysis_t<1>* solver_analysis = nullptr);
};

}
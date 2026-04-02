#pragma once
#ifndef __ANTOINE_MODEL__
#define __ANTOINE_MODEL__

/// @brief Коэффициенты уравнения Антуана. Нулевой коэффициент - А и т.д.
typedef std::array<double, 6> antoine_coefficients_t;

/// @brief Тип формулы давления насыщенных паров
/// (номера нужны для импорта из БД!!!)
enum class antoine_formula {
    /// @brief десятичный логарифм, бар, кельвины
    LgBarKelvin = 0,
    /// @brief десятичный логарифм, бар, цельсии
    LgBarCelcium = 2,
    /// @brief десятичный логарифм, мм. рт. ст.цельсии
    LgMmHgCelcium = 3,
    /// @brief натуральный логарифм, мм. рт. ст., кельвины
    LnMmHgKelvin = 1,
    /// @brief натуральный логарифм, кПа, кельвины
    ExtLnKPaKelvin = 4
};

/// @brief Модель Антуана для давления насыщенных паров
struct antoine_model_t {
    /// @brief размерности уравнения Антуана
    antoine_formula formula{ antoine_formula::LnMmHgKelvin };
    /// @brief коэффициенты в уравнении Антуана
    antoine_coefficients_t antoine_coefficients{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    /// @brief Нижняя граница применимости урав. Антуана
    double min_bound{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Верхняя граница применимости урав. Антуана
    double max_bound{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief экстраполирующий коэффициент в формуле Антуана
    double extrapolation_coefficient{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief расчет давления насыщенных паров
    double get_saturated_pressure(double temperature) const;
    /// @brief расчет температуры, при которой давление насыщенных паров равно заданному
    double get_temperature_for_saturated_pressure(double saturated_pressure) const;
};

#endif

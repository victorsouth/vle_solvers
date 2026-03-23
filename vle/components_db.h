#pragma once
#ifndef __COMPONENTS_DB__
#define __COMPONENTS_DB__

#include <limits>
#include <array>
#include <map>
#include <set>
#include <unordered_map>
#include <mutex>
#include <optional>
#include <utility>
#include "vle_solvers.h"
//using namespace std;

/// @brief Используемые единицы количества вещества (мольные, массовые)
enum class AmountType { Molar, Mass };

/// @brief Коэффициенты уравнения Антуана. Нулевой коэффициент - А и т.д.
typedef std::array<double, 6> antoine_coefficients_t;

/// @brief Аппроксимационный полином теплоемкости
typedef std::vector<double> heat_capacity_coefficients_t;

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
    antoine_formula formula{ antoine_formula::LnMmHgKelvin};
    /// @brief коэффициенты в уравнении Антуана
    antoine_coefficients_t antoine_coefficients{ std::numeric_limits<double>::quiet_NaN() , std::numeric_limits<double>::quiet_NaN() , std::numeric_limits<double>::quiet_NaN() };
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

/// @brief Коэффициенты термодинамических функций (теплоемкость Cp, энтропия, энтальпия)
struct thermodynamic_functions_coefficients_t {
    /// @brief полиномиальная зависимость теплоемкости Cp от температуры
    heat_capacity_coefficients_t heat_capacity;
    /// @brief для формулы энтальпии
    double enthalpy;
    /// @brief для формулы энтропии
    double entropy;
    // далее фикс для энтропии (что за "фикс"?)
    /// @brief Должен использоваться дефолтный конструктор.
    thermodynamic_functions_coefficients_t() = default;
    /// @brief Должен использоваться дефолтный конструктор копии.
    thermodynamic_functions_coefficients_t(const thermodynamic_functions_coefficients_t&) = default;
};

struct component_properties_t;

/// @brief термодинамические функции (теплоемкость Cp, энтропия, энтальпия)
class thermodynamic_functions_t 
    : public fixed_solvers::ranged_function_t<thermodynamic_functions_coefficients_t>
{
    friend component_properties_t;
public:
    /// @brief Возвращает размерную идеальногазовую теплоемкость индивидуального компонента.
    double get_Cp_molar(double temperature) const;
    /// @brief Возвращает размерную идеальногазовую энтальпию индивидуального компонента.
    double get_enthalpy_gas_molar(double temperature) const;
    /// @brief Возвращает размерную идеальногазовую стандартную энтропию индивидуального компонента.
    double get_entropy_molar(double temperature) const;
    /// @brief Должен использоваться дефолтный конструктор для инициализации. (там где нет RAII)
    thermodynamic_functions_t() = default; 
    /// @brief Должен использоваться дефолтный конструктор копии.
    thermodynamic_functions_t(const std::vector<fixed_solvers::function_range_t<thermodynamic_functions_coefficients_t>>& ranges);
};

/// @brief Корреляция Ван-Вельцена для вязкости
struct van_velzen_viscosity_correlation
{
    /// @brief коэффициент эмпирической зависимости am[20] -1 при отсутствии в данных вместо мусора
    double B{std::numeric_limits<double>::quiet_NaN()};
    /// @brief коэффициент эмпирической зависимости am[21] -1 при отсутствии в данных вместо мусора
    double T0{std::numeric_limits<double>::quiet_NaN()}; 
};

/// @brief Параметры чистого вещества
struct component_properties_t {
    /// @brief название (формула)
    std::wstring name; 
    /// @brief название (название)
    std::wstring component_name;
    /// @brief название (формула)
    std::wstring CASno;
    /// @brief молярная масса
    double molar_mass; 
    /// @brief плотность жидкости при 20 град
    double density_liquid_20; 
    /// @brief коэффициент сжимаемости жидкости
    double elastic_modulus;
    /// @brief температура кипения при нормальных условиях
    double normal_boiling_temperature; 
    /// @brief критическая температура
    double critical_temperature; 
    /// @brief критическое давление
    double critical_pressure; 
    /// @brief критический молярный объем
    double critical_molarvolume; 
    /// @brief фактор ацентричности Питцера
    double acentric_factor; 
    /// @brief теплота конденсации
    double condensation_heat_molar; 
    /// @brief газокинетический диаметр am[18] -1 при отсутствии в данных вместо мусора
    double gas_kinetic_diameter{std::numeric_limits<double>::quiet_NaN()};
    /// @brief равновесная энергия      am[19] -1 при отсутствии в данных вместо мусора
    double equilibrium_energy{std::numeric_limits<double>::quiet_NaN()};

    /// @brief Корреляция Ван-Вельцена для вязкости am[20,21]
    van_velzen_viscosity_correlation viscosity_correlation;

    /// @brief модель давления насыщенных паров Антуана
    antoine_model_t antoine_model;
    /// @brief термодинамические функции компонента в газовом состоянии (теплоемкость Cp, энтропия, энтальпия)
    thermodynamic_functions_t functions;
    /// @brief коэффициенты теплоемкости жидкой фазы 
    /// (теплоемкость мольная, проверено по воде и википедии 27.09.2022)
    fixed_solvers::ranged_polynom_t<heat_capacity_coefficients_t> heat_capacity_liquid;

    /// @brief удельная энтальпия вещества в жидком состоянии
    template <AmountType amount_type>
    double get_enthalpy_liquid(double pressure, double temperature) const;

    /// @brief удельная массовая энтальпия вещества в газообразном состоянии
    template <AmountType amount_type>
    double get_enthalpy_gas(double temperature) const
    {
        if constexpr (amount_type == AmountType::Mass) {
            return functions.get_enthalpy_gas_molar(temperature) / molar_mass;
        }
        else {
            return functions.get_enthalpy_gas_molar(temperature);
        }
    }

    /// @brief удельная внутренняя энергия вещества в газовом фазовом состоянии
    template <AmountType amount_type>
    double get_inner_energy_gas(double temperature) const;

    /// @brief Удельная внутренняя энергия вещества в жидком состоянии
    template <AmountType amount_type>
    double get_inner_energy_liquid(double temperature) const;

        /// @brief удельная мольная теплоемкость вещества в газообразном состоянии
    double get_Cp_gas_molar(double temperature) const;
    /// @brief удельная массовая теплоемкость вещества в газообразном состоянии
    double get_Cp_gas_mass(double temperature) const;
    /// @brief удельная мольная энтропия вещества в газообразном состоянии
    double get_entropy_gas_molar(double temperature) const;
    /// @brief удельная массовая энтропия вещества в газообразном состоянии
    double get_entropy_gas_mass(double temperature) const;
    /// @brief расчет давления насыщенных паров по формуле экстраполяции
    double get_saturated_pressure(double temperature) const;
    /// @brief Расчет производной давления насыщенных паров, 
    /// численный расчет вызывает get_saturated_pressure
    double get_saturated_pressure_derivative(double temperature) const;

    /// @brief расчет давления насыщенных паров по формуле экстраполяции
    double get_saturated_pressure_extrapolation(double temperature) const;
    /// @brief Оценка коэффициента экстраполяции давления насыщенных паров
    double estimate_antoine_extrapolation_coeff() const;
};


/// @brief Запись бинарного коэффициента взаимодействия (BIP).
/// Содержит CAS-номера двух компонентов и коэффициент k_ij.
struct bip_record_t {
    /// @brief первый CAS-номер записи
    std::wstring cas1;
    /// @brief второй CAS-номер записи
    std::wstring cas2;
    /// @brief бинарный коэффициент
    double bip_value;
};

/// @brief База данных компонентов, индексированная по CAS-номеру.
using components_database_t = std::unordered_map<std::wstring, component_properties_t> ;

/// @brief Набор бинарных коэффициентов взаимодействия.
using bip_records_t = std::vector<bip_record_t>;


/// @brief Для использования const std::set<std::wstring>& pair как ключа в unordered_map
struct wstring_pair_set_hash_t {
    size_t operator()(const std::set<std::wstring>& pair) const noexcept
    {
        if (pair.size() != 2) {
            // при указании гарантии noexcept?
            throw std::runtime_error("Wrong pair set size");
        }
        const size_t h1 = std::hash<std::wstring>{}(*pair.begin());
        const size_t h2 = std::hash<std::wstring>{}(*pair.rbegin());
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

/// @brief База данных термодинамических свойств компонентов.
/// Содержит CAS-базу, отображение формула->CAS и бинарные коэффициенты взаимодействия.
/// После инициализации является неизменяемой.
class thermo_db_t {
public:
    /// @brief Инициализация базы данных по глобальной базе компонентов (по формулам).
    /// Формирует CAS-базу, отображение формул и загружает бинарные коэффициенты.
    thermo_db_t();

    /// @brief Инициализация базы данных из сериализованной строки по CAS-номерам.
    /// Используется для альтернативного пути загрузки.
    thermo_db_t(const std::string& thermo_db);
    /// @brief Возвращает базу данных компонентов, индексированную по CAS-номеру.
    const components_database_t& get_component_casno_database() const;
    /// @brief Возвращает свойства компонента по химической формуле.
    /// Если формула отсутствует, возвращает nullptr.
    const component_properties_t& get_component_by_formula(const std::wstring& formula) const;
    /// @brief Возвращает свойства компонента по CAS-номеру.
    /// Если компонент отсутствует, возвращает nullptr.
    const component_properties_t& get_component_by_casno(const std::wstring& casno) const;
    /// @brief Возвращает свойства компонента по химической формуле.
    /// Если формула отсутствует, возвращает nullptr.
    component_properties_t& get_component_by_formula(const std::wstring& formula);
    /// @brief Возвращает свойства компонента по CAS-номеру.
    /// Если компонент отсутствует, возвращает nullptr.
    component_properties_t& get_component_by_casno(const std::wstring& casno);
    /// @brief Возвращает CAS-номер компонента по химической формуле.
    /// Если формула отсутствует или неоднозначна, кидает исключение
    const std::wstring& get_casno_by_formula(const std::wstring& formula) const;
    /// @brief Возвращает бинарный коэффициент взаимодействия k_ij по паре формул.
    /// Если коэффициент отсутствует, возвращает std::nullopt.
    double get_bip_pair_formula(const std::wstring& formula1, const std::wstring& formula2) const;
    /// @brief Возвращает бинарный коэффициент взаимодействия k_ij по паре CAS-номеров.
    /// Если коэффициент отсутствует, возвращает std::nullopt.
    double get_bip_pair_casno(const std::wstring& cas1, const std::wstring& cas2) const;

private:
    /// @brief Расчёт коэффициента экстраполяции модели Антуана для всех компонентов.
    void init_extrapolation_coeff();
    /// @brief Инициализация бинарных коэффициентов взаимодействия из глобальной базы
    void init_bips(const bip_records_t& bip_records);
private:
    /// @brief CAS-база компонентов
    components_database_t components;
    /// @brief Формула -> CAS
    std::unordered_multimap<std::wstring, std::wstring> formula2cas_mapping_var1;
    /// @brief Бинарные коэффициенты взаимодействия
    std::unordered_map<std::set<std::wstring>, double, wstring_pair_set_hash_t> bips;
};

/// @brief Сериализованная база данных компонентов, индексированная по химическим формулам.
extern const char* thermo_db_serialized_by_formula;

/// @brief Глобальная база данных компонентов, индексированная по химическим формулам.
//extern const components_database_t components_database_by_formula;

/// @brief Глобальный набор бинарных коэффициентов взаимодействия.
extern const bip_records_t bip_records_global;

/// @brief Глобальная термодинамическая база данных (CAS-база).
extern const thermo_db_t components_database;

#endif

/// @brief удельная внутренняя энергия вещества в газовом фазовом состоянии


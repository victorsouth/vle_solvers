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
#include <ranges>
#include "vle_solvers.h"

// Был объявлен в #include "fluid\fluid_base.h"

/// @brief Состав потока
struct fluid_stubdata_t {
    /// @brief Список чистых компонентов из БД
    std::vector<std::wstring> component_list;
    /// @brief Мольный состав смеси
    std::vector<double> molar_fraction;

#ifdef VLELIB_SERIALIZATION_SUPPORT
    /// @brief Для сериализации
    friend class boost::serialization::access;
    /// @brief Интрузивная сериализация/десериализация
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_NVP(component_list);
        ar& BOOST_SERIALIZATION_NVP(molar_fraction);
    }
#endif
};


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
using components_database_t = std::unordered_map<std::wstring, component_properties_t>;

/// @brief Набор бинарных коэффициентов взаимодействия.
using bip_records_t = std::vector<bip_record_t>;


/// @brief Для использования const std::set<std::wstring>& pair как ключа в unordered_map
struct wstring_pair_set_hash_t {
    /// @brief Собственно, хэш
    size_t operator()(const std::set<std::wstring>& pair) const 
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
    /// @brief Основной конструктор, инициализация чистой БД и псевдо БД
    thermo_db_t(const components_database_t& pure_db, 
        const bip_records_t& bip_db,
        const components_database_t& pseudo_db);
    /// @brief Делегирующий конструктор, создает БД чистых по глобальной базе, 
    /// дополняет ее переданной БД псевдокомпонентов
    thermo_db_t(const components_database_t& pseudo_db);
    /// @brief Инициализация базы данных по глобальной базе чистых компонентов (по формулам).
    /// Формирует CAS-базу, отображение формул и загружает бинарные коэффициенты.
    thermo_db_t();
public:
    /// @brief Возвращает базу данных компонентов, индексированную по CAS-номеру.
    const components_database_t& get_component_casno_database() const;
    /// @brief Возвращает свойства компонента по химической формуле.
    /// Если формула отсутствует, возвращает nullptr.
    const component_properties_t& get_component_by_formula(const std::wstring& formula) const;
    /// @brief Возвращает свойства компонента по CAS-номеру.
    /// Если компонент отсутствует, возвращает nullptr.
    const component_properties_t& get_component_by_casno(const std::wstring& casno) const;
    /// @brief Возвращает CAS-номер компонента по химической формуле.
    /// Если формула отсутствует или неоднозначна, кидает исключение
    const std::wstring& get_casno_by_formula(const std::wstring& formula) const;
    /// @brief Возвращает бинарный коэффициент взаимодействия k_ij по паре формул.
    /// Если коэффициент отсутствует, возвращает std::nullopt.
    double get_bip_pair_formula(const std::wstring& formula1, const std::wstring& formula2) const;
    /// @brief Возвращает бинарный коэффициент взаимодействия k_ij по паре CAS-номеров.
    /// Если коэффициент отсутствует, возвращает std::nullopt.
    double get_bip_pair_casno(const std::wstring& cas1, const std::wstring& cas2) const;

public:
    /// @brief Функция создает поток флюида с заданным компонентным составом и мольными долями
    /// @tparam Fluid
    /// @param component_formulas Названия компонентов
    /// @param molar_fractions Мольные доли компонентов
    template <typename Fluid>
    inline std::unique_ptr<Fluid> create_fluid(const std::vector<std::wstring>& component_formulas
        , const std::vector<double>& molar_fractions = std::vector<double>()) const;

    template <typename Fluid>
    inline std::unique_ptr<Fluid> create_fluid(const fluid_stubdata_t& fluid_data) const;

    /// @brief Создаёт копию объекта Fluid с возможной заменой компонентного состава
    /// ***************************************************************************
    /// @detail Функция дублирует объект жидкости Fluid, копируя все его внутренние
    /// параметры, включая концентрации, бинарные коэффициенты и используемую
    /// термодинамическую базу данных.
    /// Если список dest_components_casno пустой, создаётся точная копия исходного
    /// объекта Fluid.
    /// Если список dest_components_casno задан, функция пытается сопоставить
    /// компоненты исходного объекта с компонентами, указанными в списке CAS‑номеров.
    /// Порядок компонентов в dest_components_casno определяет порядок компонентов
    /// в результирующем объекте Fluid.
    /// Если какой‑либо компонент из dest_components_casno отсутствует в базе данных
    /// исходного Fluid, выбрасывается исключение std::logic_error.
    /// @tparam Fluid — тип создаваемого объекта жидкости
    /// @param src_fluid_ptr — указатель на исходный объект Fluid, который требуется
    /// дублировать. Не должен быть nullptr.
    /// @param dest_components_casno — список CAS‑номеров компонентов, определяющий
    /// желаемый порядок и состав компонентов в результирующем объекте Fluid.
    /// Если список пустой, используется состав исходного Fluid.
    /// @return std::unique_ptr<Fluid> — умный указатель на созданный объект Fluid.
    template <typename Fluid>
    static std::unique_ptr<Fluid> duplicate_fluid(const Fluid* src_fluid_ptr
        , const std::vector<std::wstring>& dest_components_casno = {});
private:
    /// @brief Расчёт коэффициента экстраполяции модели Антуана для всех компонентов.
    void init_extrapolation_coeff();
    /// @brief Инициализация бинарных коэффициентов взаимодействия из глобальной базы
    void init_bips(const bip_records_t& bip_records);
private:
    /// @brief CAS-база компонентов
    components_database_t cas_components;
    /// @brief Формула -> CAS
    std::unordered_multimap<std::wstring, std::wstring> formula2cas_mapping;
    /// @brief Бинарные коэффициенты взаимодействия
    std::unordered_map<std::set<std::wstring>, double, wstring_pair_set_hash_t> bips;
};
//*****************************************************************************



/// @brief Сериализованная база данных компонентов, индексированная по химическим формулам.
//extern const char* thermo_db_serialized_by_formula;
extern const char* get_thermo_db_serialized_by_formula();

/// @brief Глобальный набор бинарных коэффициентов взаимодействия.
extern const bip_records_t& get_bip_records_global();

/// @brief Глобальная термодинамическая база данных (CAS-база).
extern const thermo_db_t components_database;
//*****************************************************************************



template <typename Fluid>
inline std::unique_ptr<Fluid> thermo_db_t::create_fluid(
    const std::vector<std::wstring>& component_formulas
    , const std::vector<double>& molar_fractions) const
{
    std::vector<const component_properties_t*> components_local;
    components_local.reserve(component_formulas.size());

    for (const std::wstring& formula : component_formulas) {

        // 1) Попытка найти компонент в БД
        try {
            const std::wstring& cas = get_casno_by_formula(formula);
            const auto& component_properties = cas_components.at(cas);
            components_local.emplace_back(&component_properties);
            continue;
        }
        catch (const std::exception&) {
            // формула не найдена в БД
        }

        // 2) Ошибка
        std::stringstream msg;
        msg << "Component does not exist in thermoDB: "
            << fixed_solvers::wide2string(formula);
        throw std::logic_error(msg.str());

    }

    if (molar_fractions.empty()) {
        return std::make_unique<Fluid>(components_local);
    }

    Eigen::VectorXd fractions = Eigen::VectorXd::Map(
        molar_fractions.data(),
        molar_fractions.size()
    );

    return std::make_unique<Fluid>(components_local, fractions);
}
//*****************************************************************************



//*****************************************************************************



//TODO !рефактор!
/// @brief Функция создает поток из состава потока
/// @tparam Fluid Тип выходного потока
/// @param fluid_data Состав потока
/// @return Уникальный указатель на поток
template <typename Fluid>
inline std::unique_ptr<Fluid> thermo_db_t::create_fluid(const fluid_stubdata_t & fluid_data) const
{
    return create_fluid<Fluid>(fluid_data.component_list, fluid_data.molar_fraction);
}
//*****************************************************************************



//*****************************************************************************




template <typename Fluid>
std::unique_ptr<Fluid> thermo_db_t::duplicate_fluid(
    const Fluid* src_fluid_ptr, const std::vector<std::wstring>& dest_components_casno)
{
    if (!src_fluid_ptr) {
        throw std::logic_error("Source fluid pointer is null");
    }

    // если список пустой, просто копируем
    size_t dest_size = dest_components_casno.size();
    if (!dest_size) {
        return std::make_unique<Fluid>(*src_fluid_ptr);
    }

    const std::vector<const component_properties_t*>& src_comp_props = src_fluid_ptr->get_components();
    size_t src_size = src_comp_props.size();
    if (dest_size > src_size) {
        throw std::logic_error("Destination has more components than source");
    }

    // для ускорения поиска создаем отображение CASno -> индекс в src_fluid_ptr
    std::unordered_map<std::wstring, size_t> casno_to_src_index;
    casno_to_src_index.reserve(src_size);
    for (size_t src_index = 0; src_index < src_size; ++src_index) {
        casno_to_src_index[src_comp_props[src_index]->CASno] = src_index;
    }

    // нужно перенести:
    // 1. свойства компонентов
    std::vector<const component_properties_t*> dest_comp_props(dest_size, nullptr);
    std::vector<size_t> map_dest_to_src(dest_size);

    for (size_t dest_i = 0; dest_i < dest_size; ++dest_i) {
        // дубликатов в списке dest_components_casno быть не должно, так что можно не проверять
        const auto& dest_casno = dest_components_casno[dest_i];
        auto it = casno_to_src_index.find(dest_casno);
        if (it == casno_to_src_index.end()) {
            // компонент не найден
            std::stringstream msg;
            msg << "Component with CAS " << fixed_solvers::wide2string(dest_casno)
                << " is not found in source fluid";
            throw std::logic_error(msg.str());
        }
        const size_t src_i = it->second;
        dest_comp_props[dest_i] = src_comp_props[src_i];
        map_dest_to_src[dest_i] = src_i;
    }

    // 2. мольные доли. Они могут меняться, но mutex в самом геттере
    const Eigen::VectorXd& src_molar_fractions = src_fluid_ptr->get_molar_fraction();
    Eigen::VectorXd dest_molar_fractions(dest_size);
    for (size_t dest_i = 0; dest_i < dest_size; dest_i++) {
        const size_t src_i = map_dest_to_src[dest_i];
        dest_molar_fractions(dest_i) = src_molar_fractions(src_i);
    }

    // 3. матрицу бинарных коэффициентов взаимодействия
    const Eigen::MatrixXd& src_bips_matrix = src_fluid_ptr->get_binary_coeffs_ref();
    // Если матрица пустая, то и копировать нечего
    if( !src_bips_matrix.size() ) {
        return std::make_unique<Fluid>(dest_comp_props, dest_molar_fractions);
    }

    Eigen::MatrixXd dest_bips(dest_size, dest_size);
    for (size_t dest_i = 0; dest_i < dest_size; dest_i++) {
        for (size_t dest_j = 0; dest_j < dest_size; dest_j++) {
            const size_t src_i = map_dest_to_src[dest_i];
            const size_t src_j = map_dest_to_src[dest_j];
            dest_bips(dest_i, dest_j) = src_bips_matrix(src_i, src_j);
        }
    }

    return std::make_unique<Fluid>(dest_comp_props, dest_molar_fractions, dest_bips);
}
//*****************************************************************************



#endif


#pragma once
#ifndef __THERMO_DB__
#define __THERMO_DB__

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
    /// Если формула отсутствует, кидает исключение.
    const component_properties_t& get_component_by_formula(const std::wstring& formula) const;
    /// @brief Возвращает свойства компонента по CAS-номеру.
    /// Если компонент отсутствует, кидает исключение.
    const component_properties_t& get_component_by_casno(const std::wstring& casno) const;
    /// @brief Возвращает CAS-номер компонента по химической формуле.
    /// Если формула отсутствует или неоднозначна, кидает исключение
    const std::wstring& get_casno_by_formula(const std::wstring& formula) const;
    /// @brief Возвращает CAS-номера компонентов по списку химических формул.
    /// Если хотя бы одна формула отсутствует или неоднозначна, кидает исключение.
    std::vector<std::wstring> get_casno_by_formulas(const std::vector<std::wstring>& formulas) const;
    /// @brief Возвращает свойства компонентов по списку идентификаторов (CAS или формула).
    /// Если хотя бы один компонент отсутствует или неоднозначен, кидает исключение.
    std::vector<const component_properties_t*> get_components(
        const std::vector<std::wstring>& component_list) const;
    /// @brief Возвращает бинарный коэффициент взаимодействия k_ij по паре формул.
    /// Если коэффициент отсутствует, возвращает NaN.
    double get_bip_pair_formula(const std::wstring& formula1, const std::wstring& formula2) const;
    /// @brief Возвращает бинарный коэффициент взаимодействия k_ij по паре CAS-номеров.
    /// Если коэффициент отсутствует, возвращает NaN.
    double get_bip_pair_casno(const std::wstring& cas1, const std::wstring& cas2) const;
    /// @brief Возвращает базу бинарных коэффициентов взаимодействия.
    const bip_database_t& get_bip_db() const;

public:

    /// @brief Создаёт флюид с пересчётом матрицы бинарных коэффициентов взаимодействия.
    template <typename Fluid>
    inline std::unique_ptr<Fluid> create_fluid(const std::vector<std::wstring>& component_list
        , const std::vector<double>& molar_fractions = std::vector<double>()
        , const bip_estimation_plan_t& bip_recalc_plan = bip_estimation_plan_t()
    ) const
    {
        // Заложим на будущее
        // static_assert(
        //     std::is_same_v<Fluid, vlelib::fluid_peng_robinson_t>,
        //     "thermo_db_t::create_fluid: BIP recalculation is supported only for fluid_peng_robinson_t"
        //     );

        std::vector<const component_properties_t*> components_local = get_components(component_list);

        std::vector<std::wstring> components_casno_list;
        components_casno_list.reserve(component_list.size());
        for (const auto& comp_cas : components_local) {
            components_casno_list.push_back(comp_cas->CASno);
        }

        Eigen::MatrixXd binary_coeffs_local;
        if (!bip_recalc_plan.empty()) 
        {
            binary_coeffs_local = estimate_bip_matrix(
                components_casno_list, components_local, bip_db, bip_recalc_plan);
        }

        if (!molar_fractions.empty()) {
            Eigen::VectorXd fractions = Eigen::VectorXd::Map(
                molar_fractions.data(),
                molar_fractions.size()
            );
            return std::make_unique<Fluid>(components_local, fractions, binary_coeffs_local);
        }
        else {
            return std::make_unique<Fluid>(components_local, Eigen::VectorXd(), binary_coeffs_local);
        }

    }

    /// @brief Создаёт флюид из ранее сохраненной стабдаты
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
    components_database_t cas_components_db;
    /// @brief Формула -> CAS
    std::unordered_multimap<std::wstring, std::wstring> formula2cas_mapping;
    /// @brief Бинарные коэффициенты взаимодействия
    bip_database_t bip_db;
};
//*****************************************************************************

/// @brief Сериализованная база данных компонентов, индексированная по химическим формулам.
//extern const char* thermo_db_serialized_by_formula;
extern const char* get_thermo_db_serialized_by_formula();

/// @brief Сериализованная в JSON база бинарных коэффициентов взаимодействия.
extern const char* get_components_db_data_BIP();

/// @brief Глобальная термодинамическая база данных (CAS-база).
extern const thermo_db_t components_database;
//*****************************************************************************

;
//*****************************************************************************



//*****************************************************************************

//TODO !рефактор!
/// @brief Функция создает поток из состава потока
/// @tparam Fluid Тип выходного потока
/// @param fluid_data Состав потока
/// @return Уникальный указатель на поток
template <typename Fluid>
inline std::unique_ptr<Fluid> thermo_db_t::create_fluid(const fluid_stubdata_t& fluid_data) const
{
    return create_fluid<Fluid>(fluid_data.component_list, fluid_data.molar_fraction);
}
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
    if (!src_bips_matrix.size()) {
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

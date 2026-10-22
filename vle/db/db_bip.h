#pragma once
#ifndef __COMPONENTS_DB_BIP__
#define __COMPONENTS_DB_BIP__

/// @brief Запись бинарного коэффициента взаимодействия (BIP).
/// Содержит CAS-номера двух компонентов и коэффициент k_ij.
struct bip_record_t {
    /// @brief первый CAS-номер записи
    std::wstring cas1;
    /// @brief второй CAS-номер записи
    std::wstring cas2;
    /// @brief бинарный коэффициент
    double bip_value = std::numeric_limits<double>::quiet_NaN();
};


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


/// @brief Тип корреляции для пересчёта бинарных коэффициентов.
enum class bip_correlation_t {
    Nishiumi,
    ChuehPrausnitz,
    Gao
};

/// @brief Способ задания бинарного коэффициента
enum class bip_estimation_rule_t {
    /// @brief установить фиксированное значение для всех пар в диапазоне
    set_value,  
    /// @brief использовать из БД, если есть, иначе задать нулевым
    use_db_only,    
    /// @brief пересчитать по указанной корреляции 
    use_correlation_only,
    /// @brief использовать из БД, если есть, иначе по корреляции
    use_db_or_correlation
};

/// @brief Правило пересчёта BIP для заданной пары веществ
struct bip_estimation_plan_entry_t {
    /// @brief Пара компонентов, для которых применяется правило. Указаны CAS-номера компонентов.
    std::set<std::wstring> cas_pair;
    /// @brief Подход к расчету
    bip_estimation_rule_t rule;
    union {
        /// @brief Выбранная корреляция для пересчёта BIP.
        bip_correlation_t correlation;
        /// @brief Фиксированное значение BIP (для set_value), иначе NaN.
        double value = std::numeric_limits<double>::quiet_NaN();
    };

};

using bip_estimation_plan_t = std::vector<bip_estimation_plan_entry_t>;


/// @brief Хеш-таблица BIP по множеству формул(?) компонентов.
using bip_database_t = std::unordered_map<std::set<std::wstring>, double, wstring_pair_set_hash_t>;

struct component_properties_t;

/// @brief Расчет парного коэффициента по указанной корреляции для пары компонентов.
double estimate_bip(
    const component_properties_t& component1, const component_properties_t& component2,
    bip_correlation_t bip_correlation);


/// @brief Корреляция Чуэ–Праусница для оценки коэффициента парного взаимодействия.
/// American Institute of Chemical Engineers Journal, Vol. 13, No. 6, 1967,
/// pp. 1099–1113.
double estimate_bip_ChuehPrausnitz(
    const component_properties_t& component1, const component_properties_t& component2);

/// @brief Корреляция Гао для оценки коэффициента парного взаимодействия.
/// Fluid Phase Equilibria, 74 (1992) 85–93.
double estimate_bip_Gao(
    const component_properties_t& component1, const component_properties_t& component2);


/// @brief Генерирует план расчета всех пар BIP по заданной корреляции
bip_estimation_plan_t generate_bip_estimation_plan_with_correlation(
    const std::vector<std::wstring>& component_list, bip_correlation_t correlation);

/// @brief Генерирует план расчета всех пар BIP по правилу use_db_only.
bip_estimation_plan_t generate_bip_estimation_plan_use_db_only(
    const std::vector<std::wstring>& component_list);


/// @brief Расчет матрицы бинарных коэффициентов взаимодействия для 
/// заданного состава компонентов и плана пересчёта BIP.
Eigen::MatrixXd estimate_bip_matrix(
    const std::vector<std::wstring>& components_casno_list, 
    const std::vector<const component_properties_t*>& components, 
    const bip_database_t& bip_db, 
    const bip_estimation_plan_t& bip_estimation_plan);

//#ifndef __COMPONENTS_DB_BIP__
#endif


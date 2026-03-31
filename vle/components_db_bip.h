#pragma once
#ifndef __COMPONENTS_DB_BIP__
#define __COMPONENTS_DB_BIP__

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
enum bip_correlation_t {
    Nothing,
    Nishiumi,
    ChuehPrausnitz,
    Gao,
    SetValue
};

/// @brief Правило пересчёта BIP для заданного набора индексов.
struct bip_estimation_rule_t {
    /// @brief Индексы компонентов, для которых применяется корреляция.
    std::vector<std::pair<size_t, size_t>> indexes_ranges;

    /// @brief Выбранная корреляция для пересчёта BIP.
    bip_correlation_t correlation;

    /// @brief Фиксированное значение BIP (для SetValue), иначе NaN.
    double value = std::numeric_limits<double>::quiet_NaN();
};

/// @brief План пересчёта BIP как набор правил.
using bip_recalc_plan_t = std::vector<bip_estimation_rule_t>;

/// @brief Хеш-таблица BIP по множеству формул(?) компонентов.
using bips_hashedmap_t = std::unordered_map<std::set<std::wstring>, double, wstring_pair_set_hash_t>;


//#ifndef __COMPONENTS_DB_BIP__
#endif


#include "../vle_solvers.h"


bips_hashedmap_t get_bips_for_component_list(
    const std::vector<std::wstring>& components_casno_list
    , const bips_hashedmap_t& bips)
{
    bips_hashedmap_t result;

    // Быстрый поиск: превращаем список в set
    std::set<std::wstring> components_set(
        components_casno_list.begin(),
        components_casno_list.end()
    );

    // Перебираем все BIP в базе
    for (const auto& [pair_set, bip_value] : bips) {

        // pair_set — это set из двух CAS-номеров
        // Проверяем, что оба CAS входят в список компонентов
        bool all_included = true;
        for (const auto& cas : pair_set) {
            if (components_set.count(cas) == 0) {
                all_included = false;
                break;
            }
        }

        if (all_included) {
            result.emplace(pair_set, bip_value);
        }
    }

    return result;
}
//*****************************************************************************



/// @brief Корреляция Чуэ–Праусница (AIChE Journal, 1967).
inline double correlation_ChuehPrausnitz(
    const std::vector<const component_properties_t*>& components, size_t i, size_t j);

/// @brief Корреляция Гао (Fluid Phase Equilibria, 1992).
inline double correlation_Gao(
    const std::vector<const component_properties_t*>& components, size_t i, size_t j);
//*****************************************************************************




Eigen::MatrixXd estimate_BIP_formulas(
    const std::vector<std::wstring>& components_casno_list
    , const std::vector<const component_properties_t*>& components
    , const bips_hashedmap_t& bips
    , const bip_recalc_plan_t& bip_recalc_plan)
{
    const size_t N = components_casno_list.size();

    // 1. Проверяем, что нет правил с корреляцией Nishiumi — она не реализована
    for (const auto& rule : bip_recalc_plan) {
        if (rule.correlation == bip_correlation_t::Nishiumi) {
            throw std::runtime_error(
                "estimate_BIP_formulas: correlation Nishiumi is not implemented"
            );
        }
    }

    for (const auto& rule : bip_recalc_plan) {
        for (auto [i, j] : rule.indexes_ranges) {
            if (i >= N || j >= N) {
                throw std::runtime_error("estimate_BIP_formulas: index out of range");
            }
        }
    }

    // 2. Проверяем, покрывают ли правила SetValue весь диапазон индексов
    {
        std::set<size_t> covered;
        double set_value = std::numeric_limits<double>::quiet_NaN();
        bool has_setvalue = false;

        for (const auto& rule : bip_recalc_plan) {
            if (rule.correlation == bip_correlation_t::SetValue) {

                if (!std::isnan(rule.value)) {
                    if (!has_setvalue) {
                        set_value = rule.value;
                        has_setvalue = true;
                    }
                    else if (set_value != rule.value) {
                        // значит есть разные значения, будем использовать полный 
                        // цикл создания матрицы
                        break;
                    }
                }

                for (auto [i, j] : rule.indexes_ranges) {
                    for (size_t ins_index = i;ins_index <= j;++ins_index) {
                        covered.insert(ins_index);
                    }
                }
            }
        }

        // Если покрыт весь диапазон 0..N-1
        if (has_setvalue && covered.size() == N) {
            Eigen::MatrixXd M = Eigen::MatrixXd::Constant(N, N, set_value);
            return M;
        }
    }

    // 3. Сужаем BIP по умолчанию
    bips_hashedmap_t bips_for_component_list =
        get_bips_for_component_list(components_casno_list, bips);

    // 4. Строим матрицу BIP по умолчанию
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i == j) {
                continue;
            }

            std::set<std::wstring> key = {
                components_casno_list[i]
                , components_casno_list[j] };

            auto it = bips_for_component_list.find(key);
            if (it != bips_for_component_list.end()) {
                M(i, j) = it->second;
            }
        }
    }

    // 5. Применяем правила пересчёта
    for (const auto& rule : bip_recalc_plan) {

        switch (rule.correlation) {

        case bip_correlation_t::Nothing:
            // ничего не делаем
            break;

        case bip_correlation_t::SetValue:
            for (auto [i, j] : rule.indexes_ranges) {
                for (size_t ins_index1 = i;ins_index1 <= j;++ins_index1) {
                    for (size_t ins_index2 = i;ins_index2 <= j;++ins_index2) {
                        if (ins_index1 == ins_index2) {
                            continue;
                        }
                        M(ins_index1, ins_index2) = rule.value;
                        M(ins_index2, ins_index1) = rule.value;
                    }
                }
            }
            break;

        case bip_correlation_t::ChuehPrausnitz:
            for (auto [i, j] : rule.indexes_ranges) {
                for (size_t ins_index1 = i;ins_index1 <= j;++ins_index1) {
                    for (size_t ins_index2 = i;ins_index2 <= j;++ins_index2) {
                        M(ins_index1, ins_index2) =
                            correlation_ChuehPrausnitz(components, ins_index1, ins_index2);
                        M(ins_index2, ins_index1) = M(ins_index1, ins_index2);
                    }
                }
            }
            break;

        case bip_correlation_t::Gao:
            for (auto [i, j] : rule.indexes_ranges) {
                for (size_t ins_index1 = i;ins_index1 <= j;++ins_index1) {
                    for (size_t ins_index2 = i;ins_index2 <= j;++ins_index2) {
                        M(ins_index1, ins_index2) =
                            correlation_Gao(components, ins_index1, ins_index2);
                        M(ins_index2, ins_index1) = M(ins_index1, ins_index2);
                    }
                }
            }
            break;


        default:
            throw std::runtime_error("estimate_BIP_formulas: unknown correlation");
        }
    }

    return M;
}


/// @brief Корреляция Чуэ–Праусница (AIChE Journal, 1967).
inline double correlation_ChuehPrausnitz(
    const std::vector<const component_properties_t*>& components,
    size_t i, size_t j)
{
    // Внимание! Размерность?
    const double Vc1 = components[i]->critical_molarvolume;
    const double Vc2 = components[j]->critical_molarvolume;

    const double A = 1.0;
    const double B = 3.0;

    const double term = 2.0 * std::pow(Vc1 * Vc2, 1.0 / 6.0)
        / (std::pow(Vc1, 1.0 / 3.0) + std::pow(Vc2, 1.0 / 3.0));

    const double kij = A * (1.0 - std::pow(term, B));
    return kij;
}


/// @brief Корреляция Гао (Fluid Phase Equilibria, 1992).
inline double correlation_Gao(
    const std::vector<const component_properties_t*>& components,
    size_t i, size_t j)
{
    const double Tc1 = components[i]->critical_temperature;
    const double Tc2 = components[j]->critical_temperature;

    // Zc в статье Гао — константа 0.3074...
    constexpr double Zc1 = 0.30740130869870384801;
    constexpr double Zc2 = 0.30740130869870384801;

    const double term = 2.0 * std::sqrt(Tc1 * Tc2) / (Tc1 + Tc2);
    const double exponent = (Zc1 + Zc2) / 2.0;

    const double kij = 1.0 - std::pow(term, exponent);
    return kij;
}

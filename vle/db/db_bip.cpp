#include "../vle_solvers.h"

bip_database_t get_bips_subset(
    const std::vector<std::wstring>& casno_subset
    , const bip_database_t& bips)
{
    bip_database_t result;

    // Быстрый поиск: превращаем список в set
    std::set<std::wstring> components_set(
        casno_subset.begin(),
        casno_subset.end()
    );

    // Перебираем все BIP в базе
    for (const auto& [pair_set, bip_value] : bips) {

        // pair_set — это set из двух CAS-номеров
        // Проверяем, что оба CAS входят в список компонентов
        bool all_included = true;
        for (const std::wstring& cas : pair_set) {
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

double estimate_bip(const component_properties_t& component1, const component_properties_t& component2, bip_correlation_t bip_correlation)
{
    switch (bip_correlation) {
    case bip_correlation_t::ChuehPrausnitz:
        return estimate_bip_ChuehPrausnitz(component1, component2);
    case bip_correlation_t::Gao:
        return estimate_bip_Gao(component1, component2);
    }

    throw std::runtime_error("estimate_bip: unknown correlation");
}

double estimate_bip_ChuehPrausnitz(
    const component_properties_t& component1,
    const component_properties_t& component2)
{
    // Внимание! Размерность?
    const double Vc1 = component1.critical_molarvolume;
    const double Vc2 = component2.critical_molarvolume;

    const double A = 1.0;
    const double B = 3.0;

    const double term = 2.0 * std::pow(Vc1 * Vc2, 1.0 / 6.0)
        / (std::pow(Vc1, 1.0 / 3.0) + std::pow(Vc2, 1.0 / 3.0));

    const double kij = A * (1.0 - std::pow(term, B));
    return kij;
}


double estimate_bip_Gao(
    const component_properties_t& component1,
    const component_properties_t& component2)
{
    const double Tc1 = component1.critical_temperature;
    const double Tc2 = component2.critical_temperature;

    // Zc в статье Гао — константа 0.3074...
    constexpr double Zc1 = 0.30740130869870384801;
    constexpr double Zc2 = 0.30740130869870384801;

    const double term = 2.0 * std::sqrt(Tc1 * Tc2) / (Tc1 + Tc2);
    const double exponent = (Zc1 + Zc2) / 2.0;

    const double kij = 1.0 - std::pow(term, exponent);
    return kij;
}

Eigen::MatrixXd estimate_bip_matrix(const std::vector<std::wstring>& components_casno_list, 
    const std::vector<const component_properties_t*>& components, 
    const bip_database_t& bip_db, const bip_estimation_plan_t& bip_estimation_plan)
{
    const size_t N = components_casno_list.size();

    Eigen::MatrixXd bip_matrix = Eigen::MatrixXd::Zero(N, N);
    std::unordered_map<std::wstring, std::size_t> local_db;
    for (std::size_t index = 0; index < components.size(); ++index)
    {
        const auto& casno = components[index]->CASno;
        local_db[casno] = index;
    }

    for (const bip_estimation_plan_entry_t& info : bip_estimation_plan) {
        if (info.cas_pair.size() != 2) {
            throw std::runtime_error("estimate_bip_matrix: cas_pair must contain exactly 2 CAS numbers");
        }
        std::size_t i = local_db.at(*info.cas_pair.begin());
        std::size_t j = local_db.at(*info.cas_pair.rbegin());
        const component_properties_t* comp1 = components[i];
        const component_properties_t* comp2 = components[j];

        switch (info.rule)
        {
        case bip_estimation_rule_t::set_value:
            bip_matrix(i, j) = bip_matrix(j, i) = info.value;
            break;
        case bip_estimation_rule_t::use_correlation_only:
            bip_matrix(i, j) = bip_matrix(j, i) =
                estimate_bip(*comp1, *comp2, info.correlation);
            break;
        case bip_estimation_rule_t::use_db_only:
            // BIP взять из БД, если есть, иначе останется нулевым
            if (bip_db.contains(info.cas_pair)) {
                bip_matrix(i, j) = bip_matrix(j, i) = bip_db.at(info.cas_pair);
            }
            break;
        case bip_estimation_rule_t::use_db_or_correlation:
            if (bip_db.contains(info.cas_pair)) {
                bip_matrix(i, j) = bip_matrix(j, i) =
                    bip_db.at(info.cas_pair);
            }
            else {
                bip_matrix(i, j) = bip_matrix(j, i) =
                    estimate_bip(*comp1, *comp2, info.correlation);

            }
            break;
        default:
            throw std::runtime_error("estimate_bip_matrix: unknown estimation rule");
        }
    }

    return bip_matrix;
}


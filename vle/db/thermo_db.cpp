#include "../vle_solvers.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2026_02_02.h"


/// @brief Вызываем конструктор по умолчанию - инициализация БД и BIP из константных JSON
const thermo_db_t components_database; 


void thermo_db_t::init_extrapolation_coeff() {
    for (auto& [name, data] : cas_components_db)
    {
        double my_estimation = data.estimate_antoine_extrapolation_coeff();
        data.antoine_model.extrapolation_coefficient = my_estimation;
    }
}

void thermo_db_t::init_bips(const bip_records_t& bip_records)
{
    for (const auto& rec : bip_records) {
        std::set<std::wstring> key = { rec.cas1, rec.cas2 };
        bip_db[key] = rec.bip_value;
    }
}

thermo_db_t::thermo_db_t(const components_database_t& pure_db, const bip_records_t& bip_db, 
    const components_database_t& pseudo_db)
{
    // собираем из текущий базы компонентов (которая по формулам)
    for (const auto& [_, component] : pure_db) {
        const auto& casno = component.CASno;
        cas_components_db[casno] = component;
        formula2cas_mapping.insert({ component.name, casno });
    }

    // для псевдокомпонентов casno == formula
    for (const auto& [pseudo_name, properties] : pseudo_db) {
        if (cas_components_db.count(pseudo_name) != 0)
            throw std::runtime_error("Pseudocomonent name duplicates with pure");
        cas_components_db[pseudo_name] = properties;
        cas_components_db[pseudo_name].CASno = pseudo_name;
        cas_components_db[pseudo_name].component_name = pseudo_name;
        formula2cas_mapping.emplace(pseudo_name, pseudo_name);
    }
    // так как BIP нулевые для псевдокомпонентов и их взаимодействия с остальными 
    // компонентами, то не вставляем их в bips

    init_bips(bip_db);
    init_extrapolation_coeff();

}

thermo_db_t::thermo_db_t(const components_database_t& pseudo_db)
    : thermo_db_t(serializer_2026_02_02::
        deserialize_from_string<components_database_t>(std::string(
            get_thermo_db_serialized_by_formula())),
        serializer_2026_02_02::deserialize_BIP_from_string(
            get_components_db_data_BIP()),
        pseudo_db)
{

}

thermo_db_t::thermo_db_t()
    : thermo_db_t(serializer_2026_02_02::
        deserialize_from_string<components_database_t>(std::string(
            get_thermo_db_serialized_by_formula())),
        serializer_2026_02_02::deserialize_BIP_from_string(
            get_components_db_data_BIP()),
        components_database_t())
{
}

const components_database_t& thermo_db_t::get_component_casno_database() const
{
    return cas_components_db;
}

const component_properties_t& thermo_db_t::get_component_by_formula(const std::wstring& formula) const
{
    return get_component_by_casno(get_casno_by_formula(formula));
}

const component_properties_t& thermo_db_t::get_component_by_casno(const std::wstring& casno) const
{
    std::size_t ncomp = cas_components_db.count(casno);
    if (!ncomp) {
        throw std::runtime_error("CASno not found");
    }
    return cas_components_db.at(casno);
}

const std::wstring& thermo_db_t::get_casno_by_formula(const std::wstring& formula) const
{
    if (formula.empty()) {
        throw std::runtime_error("Cannot find formula");
    }
    std::size_t ncomps = formula2cas_mapping.count(formula);
    if (ncomps == 0) {
        throw std::runtime_error("Cannot find formula");
    }
    if (ncomps > 1) {
        throw std::runtime_error("Fatal error! There are several components for chemical formula ");
    }
    else {
        auto iter = formula2cas_mapping.find(formula);
        return iter->second;
    };
}

std::vector<std::wstring> thermo_db_t::get_casno_by_formulas(const std::vector<std::wstring>& formulas) const
{
    std::vector<std::wstring> casnos;
    casnos.reserve(formulas.size());

    for (const auto& formula : formulas) {
        casnos.push_back(get_casno_by_formula(formula));
    }

    return casnos;
}

std::vector<const component_properties_t*> thermo_db_t::get_components(
    const std::vector<std::wstring>& component_list) const
{
    std::vector<const component_properties_t*> components;
    components.reserve(component_list.size());

    for (const auto& component_id : component_list) {
        if (cas_components_db.count(component_id) == 1) {
            // Трактуем component_id как CAS. Дублей CAS нет, поэтому проверяем только на count == 1
            const component_properties_t& properties = cas_components_db.at(component_id);
            components.emplace_back(&properties);
        }
        else {
            // Не нашли component_id среди CAS-номеров,
            // Трактуем component_id как формулу, по которой пробуем получить CAS
            const std::wstring& cas = get_casno_by_formula(component_id);
            const component_properties_t& properties = cas_components_db.at(cas); // здесь будет ошибка, если не найдем cas
            components.emplace_back(&properties);
        }
    }

    return components;
}



double thermo_db_t::get_bip_pair_formula(const std::wstring& formula1, const std::wstring& formula2) const
{
    std::wstring cas1 = get_casno_by_formula(formula1);
    std::wstring cas2 = get_casno_by_formula(formula2);
    return get_bip_pair_casno(cas1, cas2);
}



double thermo_db_t::get_bip_pair_casno(const std::wstring& cas1, const std::wstring& cas2) const
{
    const auto key = std::set<std::wstring>{ cas1, cas2 };
    auto iter = bip_db.find(key);
    if (iter == bip_db.end()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return iter->second;
}

const bip_database_t& thermo_db_t::get_bip_db() const { 
    return bip_db; 
}



#include "../vle_solvers.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2026_02_02.h"


static thermo_db_t make_db_old()
{
    // база данных собранная по химической формуле
    thermo_db_t db = thermo_db_t();
    return db;
}


static thermo_db_t make_db_new()
{
    // база данных собранная по CAS-номеру
    //thermo_db_t db = thermo_db_t(std::string(thermo_db_serialized_by_casno));
    //return db;
    return {};
}


const thermo_db_t components_database = make_db_old();


thermo_db_t::thermo_db_t(const std::string& thermo_db_by_casno)
{
    throw("thermo_db_t::thermo_db_t(const std::string& thermo_db_by_casno) not realized");
    components = serializer_2026_02_02::deserialize_from_string<components_database_t>(
        std::string(thermo_db_by_casno));
    formula2cas_mapping_var1.clear();
    for (const auto& component_by_casno : components) {
        const auto& chemical_formula = component_by_casno.second.name;
        formula2cas_mapping_var1.insert({ chemical_formula, component_by_casno.first});
    }
    init_bips(get_bip_records_global());
    init_extrapolation_coeff();
}



thermo_db_t::thermo_db_t()
{
    // так как собираем здесь поэлементно, то очищаем.
    components.clear();
    formula2cas_mapping_var1.clear();
    
    components_database_t db_by_formula = serializer_2026_02_02::
        deserialize_from_string<components_database_t>(std::string(
            get_thermo_db_serialized_by_formula()));

    // собираем из текущий базы компонентов (которая по формулам)
    for (const auto& [formula, component] : db_by_formula) {
        const auto& casno = component.CASno;
        components[casno] = component;
        components[casno].name = formula;
        formula2cas_mapping_var1.insert({ formula, casno });
    }

    init_bips(get_bip_records_global());
    init_extrapolation_coeff();
}



void thermo_db_t::init_extrapolation_coeff()
{
    for (auto& [name, data] : components)
    {
        double my_estimation = data.estimate_antoine_extrapolation_coeff();
        data.antoine_model.extrapolation_coefficient = my_estimation;
    }
}



void thermo_db_t::init_bips(const bip_records_t& bip_records)
{
    for (const auto& rec : bip_records) {
        std::set<std::wstring> key = { rec.cas1, rec.cas2 };
        bips[key] = rec.bip_value;
    }
}

const components_database_t& thermo_db_t::get_component_casno_database() const
{
    return components;
}

const component_properties_t& thermo_db_t::get_component_by_formula(const std::wstring& formula) const
{
    return get_component_by_casno(get_casno_by_formula(formula));
}

const component_properties_t& thermo_db_t::get_component_by_casno(const std::wstring& casno) const
{
    std::size_t ncomp = components.count(casno);
    if (!ncomp) {
        throw std::runtime_error("CASno not found");
    }
    return components.at(casno);
}

const std::wstring& thermo_db_t::get_casno_by_formula(const std::wstring& formula) const
{
    if (formula.empty()) {
        throw std::runtime_error("Cannot find formula");
    }
    std::size_t ncomps = formula2cas_mapping_var1.count(formula);
    if (ncomps == 0) {
        throw std::runtime_error("Cannot find formula");
    }
    if (ncomps > 1) {
        throw std::runtime_error("Fatal error! There are several components for chemical formula ");
    }
    else {
        auto iter = formula2cas_mapping_var1.find(formula);
        return iter->second;
    };
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
    auto iter = bips.find(key);
    if (iter == bips.end()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return iter->second;
}
//*****************************************************************************



void thermo_db_t::add_hypocomps(const components_database_t& oils_hypocomps_parameters)
{
    // для псевдокомпонентов casno == formula
    for (const auto& [casno, properties] : oils_hypocomps_parameters) {
        if (components.count(casno) == 0) {
            components[casno] = properties;
            components[casno].CASno = casno;
            components[casno].component_name = casno;
            formula2cas_mapping_var1.emplace(casno, casno);
        }
    }
    // так как BIP нулевые для псевдокомпонентов и их взаимодействия с остальными 
    // компонентами, то не вставляем их в bips
}
//*****************************************************************************



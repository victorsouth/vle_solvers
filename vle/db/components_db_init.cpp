#include "../vle_solvers.h"

#define BOOST_EXCEPTION_DISABLE
#include "tdb_serialize_2026_02_02.h"


/// @brief Вызываем конструктор по умолчанию - инициализация БД и BIP из константных JSON
const thermo_db_t components_database; 



void thermo_db_t::init_bips(const bip_records_t& bip_records)
{
    for (const auto& rec : bip_records) {
        std::set<std::wstring> key = { rec.cas1, rec.cas2 };
        bips[key] = rec.bip_value;
    }
}


thermo_db_t::thermo_db_t(const components_database_t& pseudo_db)
    : thermo_db_t(serializer_2026_02_02::
        deserialize_from_string<components_database_t>(std::string(
            get_thermo_db_serialized_by_formula())),
        get_bip_records_global(),
        pseudo_db)
{

}

thermo_db_t::thermo_db_t()
    : thermo_db_t(serializer_2026_02_02::
        deserialize_from_string<components_database_t>(std::string(
            get_thermo_db_serialized_by_formula())),
        get_bip_records_global(),
        components_database_t())
{
}

const components_database_t& thermo_db_t::get_component_casno_database() const
{
    return cas_components;
}

const component_properties_t& thermo_db_t::get_component_by_formula(const std::wstring& formula) const
{
    return get_component_by_casno(get_casno_by_formula(formula));
}

const component_properties_t& thermo_db_t::get_component_by_casno(const std::wstring& casno) const
{
    std::size_t ncomp = cas_components.count(casno);
    if (!ncomp) {
        throw std::runtime_error("CASno not found");
    }
    return cas_components.at(casno);
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


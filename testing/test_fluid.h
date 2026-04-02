#pragma once

/// @brief Валидирует теплоемкости воды по данным из открытых источников 
/// https://ru.wikipedia.org/wiki/Вода
/// Адаптирован из vlelib/testing/base/thermodyn/test_thermodyn_base.cpp
TEST(ThermoDB, ValidatesWaterCp_Wiki) {

    auto fluid = components_database.create_fluid<vlelib::fluid_rault_dalton_t>({L"H2O"});

    double error_border = 0.06;

    //Мол. теплоёмк. 	75,37 Дж/(моль·К)
    double molar_heat_capacity = 75.37;
    print_errors(L"Молярная теплоёмкость:", molar_heat_capacity,
        fluid->get_heat_capacity_molar(ATMOSPHERIC_PRESSURE, 
            vle_solvers::celcium2kelvin(20)), error_border);

    //Удельная теплоёмкость воды, кДж/(кг*К) не учтено в жидкой фазе
    std::map<double, double> mass_heat_capacity = {
        {0, 	4.218},
        {10, 	4.192},
        {20,	4.182},
        {40,	4.178},
        {60,    4.184},
        {80,	4.196},
        {100,	4.216}
    };

    for(auto [temperature, heat_capacity] : mass_heat_capacity){
        double heat_capacity_calc = fluid->get_heat_capacity_mass(
            ATMOSPHERIC_PRESSURE, vle_solvers::celcium2kelvin(temperature));
        print_errors(L"Удельная теплоёмкость:", heat_capacity * 1000, heat_capacity_calc, error_border);
    }
}

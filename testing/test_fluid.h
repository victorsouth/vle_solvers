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

/// @brief Перевод массовых долей в мольные возвращает исходный мольный состав
TEST(FluidComposition, RestoresMolarFraction_WhenSetByMassFraction)
{
    // Arrange: смесь CH4/C2H6
    std::vector<std::wstring> components = { L"CH4", L"C2H6" };
    std::vector<double> molar_fraction = { 0.3, 0.7 };
    auto fluid = components_database.create_fluid<vlelib::fluid_rault_dalton_t>(components, molar_fraction);
    Eigen::VectorXd initial_molar_fraction = fluid->get_molar_fraction();
    Eigen::VectorXd mass_fraction = fluid->get_mass_fraction();
    // затираем исходный мольный состав
    fluid->set_molar_fraction(std::vector<double>{ 0.5, 0.5 });

    // Act: перевод массовых долей обратно в мольные
    fluid->set_mass_fraction(mass_fraction);

    // Assert: мольный состав совпадает с исходным
    EXPECT_LT((fluid->get_molar_fraction() - initial_molar_fraction).cwiseAbs().maxCoeff(), 1e-12);
}

#pragma once

/// @brief Тест расчета давления насыщенных паров для гликоля
TEST(ThermoEngine, ResearchGlycole) {


    const auto data = components_database.get_component_by_formula(L"C4H10O3");

    double t1 = data.get_saturated_pressure(vle_solvers::celcium2kelvin(-200));
    double t2 = data.get_saturated_pressure(vle_solvers::celcium2kelvin(-170));
    double t3 = data.get_saturated_pressure(vle_solvers::celcium2kelvin(-10));

    // Проверяем, что расчеты выполнены без ошибок
    ASSERT_FALSE(std::isnan(t1));
    ASSERT_FALSE(std::isnan(t2));
    ASSERT_FALSE(std::isnan(t3));
}

/// @brief Тест сравнения теплоемкости метана
TEST(ThermoDB, Compare1) {
    const auto component = components_database.get_component_by_formula(L"CH4");


    double density_std = 0.7168;

    double Cp = component.get_Cp_gas_mass(KELVIN_OFFSET);

    std::vector<std::wstring> component_list{ L"CH4" };
    auto fluid = components_database.create_fluid<vlelib::fluid_rault_dalton_t>(component_list);

    double Cpm = fluid->get_heat_capacity_mass(1e5, KELVIN_OFFSET);

    // Проверяем, что расчеты выполнены без ошибок
    ASSERT_FALSE(std::isnan(Cp));
    ASSERT_FALSE(std::isnan(Cpm));
}


/// @brief Пара "имя компонента + CAS номер" для использования как ключа.
struct name_casno {
    /// @brief Имя компонента.
    std::string name;
    /// @brief CAS номер компонента.
    std::string CasNo;
};


bool operator<(const struct name_casno& a, const struct name_casno& b) {
    return a.CasNo<b.CasNo;
}



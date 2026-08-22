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



/// @brief Тест для проверки экстраполяции давления паров
TEST(ComponentPropertiesTest, SaturatedPressureAtCriticalPoint) {
    // --- Arrange ---
    // Используем типичные параметры (например, для воды или метана)
    const auto& component = components_database.get_component_by_formula(L"CH4");

    // --- Act ---
    // Вызываем экстраполяцию ровно в критической точке
    double result = component.get_saturated_pressure_extrapolation(component.critical_temperature);

    // --- Assert ---
    // В критической точке P должно быть равно Pc, 
    EXPECT_NEAR(result, component.critical_pressure, 1e-7);
}


/// @brief Тест физической корректности влияния фактора ацентричности(omega).
TEST(ComponentPropertiesTest, AcentricFactorPhisicalLogicCheck) {
    // --- Arrange ---
    const auto& heavy_comp = components_database.get_component_by_formula(L"C2H5OH");
    // Температура ниже критической (T < Tc)
    double test_temp = 350.0;

    // --- Act ---
    // Физическая проверка: для ацентричных молекул (omega > 0)
    // давление паров должно быть НИЖЕ, чем у простых (omega = 0)
    double p_extrapolated = heavy_comp.get_saturated_pressure_extrapolation(test_temp);
    // Считаем значение без учета фактора ацентричности (как будто простая молекула)
    double alpha = heavy_comp.antoine_model.extrapolation_coefficient;
    double p_simple = heavy_comp.critical_pressure * std::exp(alpha *
        (1 - heavy_comp.critical_temperature / test_temp));

    // --- Assert ---
    // Если в коде используется стандартная логика (1 + omega), то p_extrapolated < p_simple
    // При вычитании (1 - omega) давление p_extrapolated будет аномально ВЫШЕ p_simpl
    EXPECT_LT(p_extrapolated, p_simple);
}

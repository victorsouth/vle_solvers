#pragma once

/// @brief Сравнивает получение BIP по формулам и по CAS-номерам.
TEST(BinaryCoefficients, BipsByFormulaAndCas_ReturnSameValue) {
    const auto methane = components_database.get_component_by_formula(L"CH4");
    const auto ethane = components_database.get_component_by_formula(L"C2H6");

    const double bip_by_formula = components_database.get_bip_pair_formula(L"CH4", L"C2H6");
    const double bip_by_cas = components_database.get_bip_pair_casno(methane.CASno, ethane.CASno);

    ASSERT_DOUBLE_EQ(bip_by_formula, bip_by_cas);
}

/// @brief Проверяет, что для отсутствующей пары CAS возвращается NaN.
TEST(BinaryCoefficients, BipsForUnknownCasPair_ReturnsNaN) {
    const double bip = components_database.get_bip_pair_casno(
        L"000-00-0", L"999-99-9");

    ASSERT_TRUE(std::isnan(bip));
}


/// @brief Сравнивает получение BIP по формулам и по CAS-номерам.
TEST(BinaryCoefficientsPlan, HandlesOptionDbOnly) {
    // Arrange
    std::vector<std::wstring> casno_list =
        components_database.get_casno_by_formulas({ L"CH4", L"C2H6" });
    std::vector<const component_properties_t*> components =
        components_database.get_components(casno_list);

    bip_estimation_plan_t plan;
    plan.push_back({
        .cas_pair = { casno_list[0], casno_list[1] },
        .rule = bip_estimation_rule_t::use_db_only,
        .correlation = bip_correlation_t::Gao
        });

    // Act
    const auto& bip_db = components_database.get_bip_db();
    Eigen::MatrixXd bip_matrix =
        estimate_bip_matrix(casno_list, components, bip_db, plan);

    // Assert - проверяем матрицу
    double etalon_bip = bip_db.at({ casno_list[0], casno_list[1] });

    ASSERT_DOUBLE_EQ(bip_matrix(0, 0), 0.0); // диагональные элементы...
    ASSERT_DOUBLE_EQ(bip_matrix(1, 1), 0.0); // ... должны быть нулевыми
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), bip_matrix(1, 0)); // симметрия

    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), etalon_bip); // элемент взят из БД

}

/// @brief Проверяет пересчёт BIP по корреляции Чуи–Праусница.
TEST(BinaryCoefficientsPlan, HandlesOptionChuehPrausnitz) {
    // Arrange
    std::vector<std::wstring> casno_list =
        components_database.get_casno_by_formulas({ L"CH4", L"C2H6" });
    std::vector<const component_properties_t*> components =
        components_database.get_components(casno_list);

    bip_estimation_plan_t plan;
    plan.push_back({
        .cas_pair = { casno_list[0], casno_list[1] },
        .rule = bip_estimation_rule_t::use_correlation_only,
        .correlation = bip_correlation_t::Chueh_Prausnitz
        });

    // Act
    const auto& bip_db = components_database.get_bip_db();
    Eigen::MatrixXd bip_matrix =
        estimate_bip_matrix(casno_list, components, bip_db, plan);

    // Assert
    const double etalon_bip = estimate_bip_ChuehPrausnitz(*components[0], *components[1]);

    ASSERT_DOUBLE_EQ(bip_matrix(0, 0), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(1, 1), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), bip_matrix(1, 0));
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), etalon_bip);
}

/// @brief Проверяет пересчёт BIP по корреляции Гао.
TEST(BinaryCoefficientsPlan, HandlesOptionGao) {
    // Arrange
    std::vector<std::wstring> casno_list =
        components_database.get_casno_by_formulas({ L"CH4", L"C2H6" });
    std::vector<const component_properties_t*> components =
        components_database.get_components(casno_list);

    bip_estimation_plan_t plan;
    plan.push_back({
        .cas_pair = { casno_list[0], casno_list[1] },
        .rule = bip_estimation_rule_t::use_correlation_only,
        .correlation = bip_correlation_t::Gao
        });

    // Act
    const auto& bip_db = components_database.get_bip_db();
    Eigen::MatrixXd bip_matrix =
        estimate_bip_matrix(casno_list, components, bip_db, plan);

    // Assert
    const double etalon_bip = estimate_bip_Gao(*components[0], *components[1]);

    ASSERT_DOUBLE_EQ(bip_matrix(0, 0), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(1, 1), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), bip_matrix(1, 0));
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), etalon_bip);
}

/// @brief Проверяет установку фиксированного ненулевого BIP через set_value.
TEST(BinaryCoefficientsPlan, HandlesOptionSetValue) {
    // Arrange
    std::vector<std::wstring> casno_list =
        components_database.get_casno_by_formulas({ L"CH4", L"C2H6" });
    std::vector<const component_properties_t*> components =
        components_database.get_components(casno_list);

    constexpr double custom_bip = 0.12345;

    bip_estimation_plan_t plan;
    plan.push_back({
        .cas_pair = { casno_list[0], casno_list[1] },
        .rule = bip_estimation_rule_t::set_value,
        .value = custom_bip
        });

    // Act
    const auto& bip_db = components_database.get_bip_db();
    Eigen::MatrixXd bip_matrix =
        estimate_bip_matrix(casno_list, components, bip_db, plan);

    // Assert
    ASSERT_DOUBLE_EQ(bip_matrix(0, 0), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(1, 1), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), bip_matrix(1, 0));
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), custom_bip);
}

/// @brief Проверяет fallback на корреляцию в use_db_or_correlation при пустой БД BIP.
TEST(BinaryCoefficientsPlan, HandlesOptionDbOrCorrelationWithEmptyDb) {
    // Arrange
    std::vector<std::wstring> casno_list =
        components_database.get_casno_by_formulas({ L"CH4", L"C2H6" });
    std::vector<const component_properties_t*> components =
        components_database.get_components(casno_list);

    bip_estimation_plan_t plan;
    plan.push_back({
        .cas_pair = { casno_list[0], casno_list[1] },
        .rule = bip_estimation_rule_t::use_db_or_correlation,
        .correlation = bip_correlation_t::Gao
        });

    // Act
    bip_database_t empty_bip_db; // Пустая БД: значение должно быть рассчитано по корреляции.
    Eigen::MatrixXd bip_matrix =
        estimate_bip_matrix(casno_list, components, empty_bip_db, plan);

    // Assert
    const double etalon_bip = estimate_bip_Gao(*components[0], *components[1]);

    ASSERT_DOUBLE_EQ(bip_matrix(0, 0), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(1, 1), 0.0);
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), bip_matrix(1, 0));
    ASSERT_DOUBLE_EQ(bip_matrix(0, 1), etalon_bip);
}

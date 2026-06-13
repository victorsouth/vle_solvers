#pragma once

/// @brief Создаёт rachford_rice2_t с заданным составом z и явным вектором K.
std::pair<std::unique_ptr<vlelib::fluid_rault_dalton_t>, vlelib::rachford_rice2_t>
make_rr_with_k_values(const std::vector<double>& z, const std::vector<double>& K)
{
    std::vector<std::wstring> components = { L"CH4", L"C2H6" };
    auto fluid = components_database.create_fluid<vlelib::fluid_rault_dalton_t>(components, z);
    Eigen::VectorXd K_values(2);
    K_values(0) = K[0];
    K_values(1) = K[1];
    vlelib::rachford_rice2_t rr(fluid.get(), K_values);
    return { std::move(fluid), std::move(rr) };
}

/// @brief Возвращает 0.0, если все K-значения меньше единицы.
TEST(RachfordRicePhysicalConstrained, ReturnsLiquid_WhenKValuesBelowOne)
{
    // Arrange: z=[0.5, 0.5], K=[0.5, 0.8], k_max < 1
    auto [fluid, rr] = make_rr_with_k_values({ 0.5, 0.5 }, { 0.5, 0.8 });

    // Act: быстрый путь без численного solve
    double V = rr.try_nonphysical_solve();

    // Assert: однофазная жидкость
    EXPECT_DOUBLE_EQ(V, 0.0);
}

/// @brief Возвращает 1.0, если все K-значения больше единицы.
TEST(RachfordRicePhysicalConstrained, ReturnsVapor_WhenKValuesAboveOne)
{
    // Arrange: z=[0.5, 0.5], K=[1.5, 2.0], k_min > 1
    auto [fluid, rr] = make_rr_with_k_values({ 0.5, 0.5 }, { 1.5, 2.0 });

    // Act: быстрый путь без численного solve
    double V = rr.try_nonphysical_solve();

    // Assert: однофазный пар
    EXPECT_DOUBLE_EQ(V, 1.0);
}

/// @brief Возвращает 1.0 при одинаковом положительном знаке граничных невязок.
TEST(RachfordRicePhysicalConstrained, ReturnsVapor_WhenBoundaryResidualsHaveSamePositiveSign)
{
    // Arrange: z=[0.1, 0.9], K=[0.2, 2.0]; F(0)=0.82>0, F(1)=0.05>0
    auto [fluid, rr] = make_rr_with_k_values({ 0.1, 0.9 }, { 0.2, 2.0 });
    const auto& z = fluid->get_molar_fraction();
    Eigen::VectorXd K_values(2);
    K_values << 0.2, 2.0;
    double f0 = vlelib::rr_equation(0.0, z, K_values);
    double f1 = vlelib::rr_equation(1.0, z, K_values);

    // Act: быстрый путь без численного solve
    double V = rr.try_nonphysical_solve();

    // Assert: F(0)>0, F(1)>0 — корня в (0, 1) нет, однофазный пар
    ASSERT_GT(f0, 0.0);
    ASSERT_GT(f1, 0.0);
    EXPECT_DOUBLE_EQ(V, 1.0);
}

/// @brief Возвращает 0.0 при одинаковом отрицательном знаке граничных невязок.
TEST(RachfordRicePhysicalConstrained, ReturnsLiquid_WhenBoundaryResidualsHaveSameNegativeSign)
{
    // Arrange: z=[0.9, 0.1], K=[0.5, 1.5]; F(0)=-0.4<0, F(1)<0
    auto [fluid, rr] = make_rr_with_k_values({ 0.9, 0.1 }, { 0.5, 1.5 });
    const auto& z = fluid->get_molar_fraction();
    Eigen::VectorXd K_values(2);
    K_values << 0.5, 1.5;
    double f0 = vlelib::rr_equation(0.0, z, K_values);
    double f1 = vlelib::rr_equation(1.0, z, K_values);

    // Act: быстрый путь без численного solve
    double V = rr.try_nonphysical_solve();

    // Assert: F(0)<0, F(1)<0 — корня в (0, 1) нет, однофазная жидкость
    ASSERT_LT(f0, 0.0);
    ASSERT_LT(f1, 0.0);
    EXPECT_DOUBLE_EQ(V, 0.0);
}

/// @brief Возвращает NaN, если на границах невязки разных знаков (двухфазный корень).
TEST(RachfordRicePhysicalConstrained, ReturnsNaN_WhenTwoPhaseRootExists)
{
    // Arrange: z=[0.3, 0.7], K=[0.3, 3.0]; F(0)>0, F(1)<0
    auto [fluid, rr] = make_rr_with_k_values({ 0.3, 0.7 }, { 0.3, 3.0 });
    const auto& z = fluid->get_molar_fraction();
    Eigen::VectorXd K_values(2);
    K_values << 0.3, 3.0;
    double f0 = vlelib::rr_equation(0.0, z, K_values);
    double f1 = vlelib::rr_equation(1.0, z, K_values);

    // Act: быстрый путь без численного solve
    double V = rr.try_nonphysical_solve();

    // Assert: F(0)>0, F(1)<0 — нужен численный solve для двухфазного корня
    ASSERT_GT(f0, 0.0);
    ASSERT_LT(f1, 0.0);
    EXPECT_TRUE(std::isnan(V));
}

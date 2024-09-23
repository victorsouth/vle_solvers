#pragma once

using namespace vle_solvers;

/// @brief ���� ��� ������� ����, ��� �������� � ����� ������
TEST(DataBase, ComponentPropertiesInitialization)
{
	// �������������� ���������
	wstring component_name = L"CO2";

	// ��������� ���������� �� ���������� �����������
	component_properties_t component_properties;

	// �������������� �������� ��� ��������� ����������
	component_properties = components_database.at(component_name);

	ASSERT_FALSE(std::isnan(component_properties.molar_mass));
}
// ����������� ���������� vle_solvers.h
#include "../src/vle_solvers.h"

// ���������� ����������� ������������ ����
using namespace vle_solvers;


int main()
{
	// ������� ����������, ������� ����� ��������� ���� ������ ����������� � �������������� ��
	components_database_t components_database;
	db_initializer_t db(components_database);

	// �������� ����� ���������, ��� �������� ���������� ���������� ��������
	wstring component = L"CO2";

	// ���������� �������� ����������
	component_properties_t component_properties = components_database.at(component);

	return 1;
}
// ����������� ���������� vle_solvers.h
#include "../src/vle_solvers.h"

// ���������� ����������� ������������ ����
using namespace vle_solvers;


int main()
{
	// �������� ����� ���������, ��� �������� ���������� ���������� ��������
	wstring component = L"CO2";

	// ���������� �������� ����������
	component_properties_t component_properties = components_database.at(component);

	return 1;
}
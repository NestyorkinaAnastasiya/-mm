#pragma once
/*grid.h*/
#define _CRT_SECURE_NO_WARNINGS
#pragma once
#include <stdio.h>
#include <vector>
#include <fstream>
#include <functional>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
using namespace std;

namespace grd
{
	//�������� �������
	struct Element
	{
		// ����� �������� ������� ��������
		int nDof;
		// ����
		int nodes[4];
		// ���������� ������ �������� ������� ��������
		vector <int> dof;
		int numberOfMaterial;

		Element& operator=(Element element)
		{
			for (int i = 0; i < 4; i++)
				nodes[i] = element.nodes[i];
			numberOfMaterial = element.numberOfMaterial;

			return *this;
		}
		// ����� ����������� ������ nodesNumber � ��������
		bool SearchNode(int nodesNumber);
		int GetGlobalDofNumber(int dofNumber);
	};
	ifstream& operator>>(std::ifstream& is, std::vector <Element>& elements);

	struct Point
	{
		double x;
		double y;
		Point() {};
		~Point() {};
		Point(double xx, double yy) {
			x = xx;
			y = yy;
		}

		bool operator==(Point point) {
			if (point.x == x && point.y == y) return true;
			else return false;
		}

	};
	ifstream& operator>>(std::ifstream& is, std::vector <Point>& points);

	//��������� ����� �������
	struct Grid
	{
		Grid();
		//������ �������� ���������
		vector <Element> elements;
		//������ �����
		vector <Point> nodes;
		//������ ����� � ������� �������� ���������
		vector <int> ku;
		//������������ ������ ���������, ���������� ���������� ����� �.�.
		//������ nodesNumber
		void SearchElements(int nodesNumber, vector<int>& elList);
		int SearchElement(double x, double y);
		double hx(int elementNumber);
		double hy(int elementNumber);
		void DataFromTelma();

		~Grid();
	};
	Grid grid;
}
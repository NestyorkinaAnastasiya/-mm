#include "grid.h"

namespace grd
{
	// �������� ����������� ������
	/*
	# inf2tr.dat (��������� ����):
	kuzlov - ����� ����� �����
	ktr    - ����� �������� ��������� (���������������)
	kt1    - ����� ����� � ������� �������� ���������
	# rz.dat (�������� ����):
	��������� i-� ������d (i=1..kuzlov):
	2*double (x,y),	��� x,y - (x,y)-���������� i-� �������
	# nvtr.dat (�������� ����):
	��������� i-� ������ (i=1..ktr):
	6*long (i1,i2,i3,i4,0,1), ��� i1,i2,i3,i4 - ���������� ������
	������ i-�� �������������� (����� �������, ������ �������,
	����� ������, ������ ������)
	# nvkat2d.dat (�������� ����):
	��������� i-� ������ (i=1..ktr):
	1*long (m), ��� m - ����� ��������� i-�� �������������� �
	������������ � ������� sreda � mu
	# l1.dat (�������� ����):
	��������� i-� ������ (i=1..kt1):
	1*long (k), ��� k - ���������� ����� i-� ������� � ������
	�������� ������� ��������
	*/

	// ���������� �����
	Grid::Grid() 
	{
		DataFromTelma();
		ifstream gridIn("grid.dt"), elementsIn("elements.dt");
		gridIn >> nodes;
		elementsIn >> elements;
	}

	void Grid::DataFromTelma()
	{
		int countOfNodes, countOfElements, countOfBC, l;
		FILE *fNodes, *fElements, *l1, *fp;
		ofstream gridOut("grid.dt"), elemOut("elements.dt");

		/* ������� fseek ���������� ��������� ������� � ������.
		������������� ���������� ��������� ��������� � �����,
		� ����� �������, ������� ������������ ����� ����������
		�������� � ��������� ���������.

		int fseek( FILE * filestream, long int offset, int origin );
		- filestream	��������� �� ������ ���� FILE, ���������������� �����.

		- offset		���������� ���� ��� ��������, ������������ ���������� ��������� ���������.

		- origin		������� ���������, ������������ ������� ����� ����������� ��������.
		����� ������� ������� ����� �� ��������� ��������, ����������� �
		������������ ����� <cstdio>:

		SEEK_SET	������ �����
		SEEK_CUR	������� ��������� �����
		SEEK_END	����� �����*/


		fopen_s(&fp, "inf2tr.dat", "r");
		char c = ' ';
		fseek(fp, 56, SEEK_SET);
		fscanf(fp, "%d", &countOfNodes);

		while (c != '=') fscanf(fp, "%c", &c);
		fscanf(fp, "%d", &countOfElements);

		c = ' ';
		while (c != '=') fscanf(fp, "%c", &c);
		fscanf(fp, "%d", &countOfBC);
		fclose(fp);
		gridOut << countOfNodes << endl;
		elemOut << countOfElements << endl;
		
		fopen_s(&fNodes, "rz.dat", "rb");
		for (int i = 0; i < countOfNodes; i++) {
			double xy[2];
			fread(xy, sizeof(double), 2, fNodes);
			gridOut << setprecision(21) << xy[0] << " " << setprecision(21) << xy[1] << endl;//?? iomanip
		}
		fclose(fNodes);

		fopen_s(&fElements, "nvtr.dat", "rb");
		fp = fopen("nvkat2d.dat", "rb");
		for (int i = 0; i < countOfElements; i++) {
			int nodes[6], material;
			fread(nodes, sizeof(int), 6, fElements);
			fread(&material, sizeof(int), 1, fp);

			// ��������� � ����� ���������� � 1
			for (int j = 0; j < 4; j++)
				nodes[j]-= 1;

			sort(nodes, nodes + 4);
			elemOut << material;
			for (int j = 0; j < 4; j++)	elemOut << " " << nodes[j];
			elemOut << endl;
		}
		fclose(fp);
		fclose(fElements);
		gridOut.close();
		elemOut.close();
		
		// ���������� ����� � �������� ���������
		fopen_s(&l1, "l1.dat", "rb");
		for (int i = 0; i < countOfBC; i++)	{
			fread(&l, sizeof(int), 1, l1);
			ku[i] = l - 1;
		}

		sort(ku.begin(), ku.end());
		fclose(l1);
	}

	Grid::~Grid() {}

	// ����� ����������� ������ �.�. � ��������
	bool Element::SearchNode(int nodesNumber)
	{
		for (int i = 0; i < 4; i++)
			if (nodesNumber == nodes[i]) return true;

		return false;
	}

	// ����� ����������� ������ �.�. � ��������
	int Element::GetGlobalDofNumber(int dofNumber)
	{
		// ��������� ����� �����, � ���-� ������������� ��
		int num1 = dofNumber / 4;  
		// ����� �� ����� ������� �� ����
		int num2 = dofNumber - num1 * 4; 
		return nodes[num1] * 4 + num2;
	}

	// ������������ ������ ������� ���������, ���������� ���������� ����� �.�.
	// ������ nodesNumber
	void Grid::SearchElements(int nodesNumber, vector <int> &elList)
	{
		int count;
		int size = elements.size();
		elList.reserve(4);

		count = 0;
		for (int i = 0; i < size && count < 4; i++)	{
			if (elements[i].SearchNode(nodesNumber)) {
				count++;
				elList.push_back(i);
			}
		}
	}
		
	int Grid::SearchElement(double x, double y)
	{
		double xLeft, xRight, yLow, yUp;
		int size = elements.size();

		for (int i = 0; i < size; i++) {
			xLeft = nodes[elements[i].nodes[0]].x;
			xRight = nodes[elements[i].nodes[1]].x;
			yLow = nodes[elements[i].nodes[0]].y;
			yUp = nodes[elements[i].nodes[3]].y;
			if (xLeft <= x && x <= xRight && yLow <= y && y <= yUp)
				return i;
		}
		return -1;
	}

	double Grid::hx(int elementNumber)
	{
		Element element = elements[elementNumber];
		return nodes[element.nodes[1]].x - nodes[element.nodes[0]].x;
	}

	double Grid::hy(int elementNumber)
	{
		Element element = elements[elementNumber];
		return nodes[element.nodes[2]].y - nodes[element.nodes[0]].y;
	}

	ifstream& operator>>(ifstream& is, vector <Element>& elements)
	{
		int tmp;
		Element element;

		is >> tmp;
		elements.reserve(tmp);

		for (int i = 0; i < tmp; i++) {
			element.nDof = 16;

			is >> element.numberOfMaterial;

			for (int j = 0; j < 4; j++)
				is >> element.nodes[j];
			element.dof.resize(element.nDof);
			for (int j = 0; j < element.nDof; j++)
				element.dof[j] = element.GetGlobalDofNumber(j);

			elements.push_back(element);
		}

		return is;
	}

	ifstream& operator>>(ifstream& is, vector <Point>& points)
	{
		int tmp;
		Point point_tmp;

		is >> tmp;

		points.reserve(tmp);
		for (int i = 0; i < tmp; i++) {
			is >> point_tmp.x >> point_tmp.y;
			points.push_back(point_tmp);
		}

		return is;
	}
}
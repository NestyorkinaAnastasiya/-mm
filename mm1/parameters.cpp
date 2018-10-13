#include "parameters.h"
// Папка с библиотекой
#include "rapidjson/include/rapidjson/document.h"
#include <fstream>
#include <map>

using namespace std;
using namespace rapidjson;

namespace parameters
{
	Parameters::~Parameters() {}

	Parameters::Parameters()
	{
		ifstream file("testing_parameters.json");
		string json, solverName;

		// Считали всё в одну строку
		while (!file.eof()) {
			string line;
			file >> line;
			json += line;
		}
		file.close();

		Document doc_testp;
		// The JSON is now parsed into document as a DOM tree:
		doc_testp.Parse(json.c_str());

		auto& j_testp = doc_testp["testp"];
		useLU = j_testp["useLU"].GetBool();
		solverName = j_testp["solver"].GetString();
		alpha = j_testp["alpha"].GetDouble();
		betta = j_testp["betta"].GetDouble();

		map<string, int> _map;

		_map["BiCGStab"] = 1;
		_map["GMRES"] = 2;
		_map["BCGandGMRES"] = 3;
		// Если не нашли, то
		if (_map.find(solverName) == _map.end()) solver = 2;
		else solver = _map[solverName];

		/* Параметры решателя ?*/
		file.open("solver.json");
		string json;

		while (!file.eof())
		{
			string line;
			file >> line;
			json += line;
		}

		file.close();

		Document doc_solver;
		doc_solver.Parse(json.c_str());

		auto& j_solver = doc_solver["solver"];
		epsilon = j_solver["epsilon"].GetDouble();
		gmresM = j_solver["gmres_m"].GetInt();
		maxCountOfIterations = j_solver["max_number_of_iterations"].GetInt();
	}
}
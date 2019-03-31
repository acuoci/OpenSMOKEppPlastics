/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2019  Alberto Cuoci                                      |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#pragma once

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

#include "ThermogravimetricAnalysis_Options.h"
#include "ExplicitOdeSolver_Parameters.h"

//#include "math/native-ode-solvers/parameters/OdeSolver_Parameters.h;"

namespace OpenSMOKE
{
	class Grammar_ThermogravimetricAnalysis : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PolymerType",
				OpenSMOKE::SINGLE_STRING,
				"Polymer type: polyethylene | polypropylene | polystyrene",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Pressure",
				OpenSMOKE::SINGLE_MEASURE,
				"Pressure",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Temperature",
				OpenSMOKE::SINGLE_MEASURE,
				"Temperature",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary defining the temperature profile",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureSlope",
				OpenSMOKE::SINGLE_MEASURE,
				"Temperature profile slope (example: 10 C/min)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EndTime",
				OpenSMOKE::SINGLE_MEASURE,
				"Total time for integration of equations (default: 10 h)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Options",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing additional options for solving the thermogravimetric analysis",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeParameters",
				OpenSMOKE::SINGLE_DICTIONARY,
				"Dictionary containing the numerical parameters for solving the stiff ODE system",
				false));
		}
	};

	enum PolymerType { POLYETHYLENE, POLYPROPYLENE, POLYSTYRENE };

	void ThermogravimetricAnalysisFromDictionary(	const std::string input_file_name,
													const std::string dictionary_name,
													double& P, double& tEnd,
													Eigen::VectorXd& profile_t, 
													Eigen::VectorXd& profile_T,
													PolymerType& polymer_type,
													OpenSMOKE::ThermogravimetricAnalysis_Options& thermo_options,
													OdeSMOKE::ExplicitOdeSolver_Parameters& ode_options)
	{
		// Define the dictionaries
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;

		OpenSMOKE::Grammar_ThermogravimetricAnalysis grammar_thermogravimetry;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(dictionary_name).SetGrammar(grammar_thermogravimetry);

		if (dictionaries(dictionary_name).CheckOption("@PolymerType") == true)
		{
			std::string polymer_name;
			dictionaries(dictionary_name).ReadString("@PolymerType", polymer_name);

			if (polymer_name == "polyethylene")
				polymer_type = POLYETHYLENE;
			else if (polymer_name == "polypropylene")
				polymer_type = POLYPROPYLENE;
			else if (polymer_name == "polystyrene")
				polymer_type = POLYSTYRENE;
			else
			{
				std::cout << "Fatal error: wrong @PolymerType" << std::endl;
				std::cout << "Available options: polyethylene | polypropylene | polystyrene" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		double T = 300.;
		if (dictionaries(dictionary_name).CheckOption("@Temperature") == true)
		{
			std::string units;
			dictionaries(dictionary_name).ReadMeasure("@Temperature", T, units);
			if (units == "K")
				T = T;
			else if (units == "C")
				T += 273.15;
			else
			{
				std::cout << "Fatal error: wrong units for temperature. Available units: K | C" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		P = 1.;
		if (dictionaries(dictionary_name).CheckOption("@Pressure") == true)
		{
			std::string units;
			dictionaries(dictionary_name).ReadMeasure("@Pressure", P, units);
			if (units == "atm")
				P = P;
			else if (units == "Pa")
				P /= 101325.;
			else if (units == "bar")
				P /= 1.01325;
			else
			{
				std::cout << "Fatal error: wrong units for pressure. Available units: atm | bar | Pa" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		tEnd = 10.*3600.;
		if (dictionaries(dictionary_name).CheckOption("@EndTime") == true)
		{
			std::string units;
			dictionaries(dictionary_name).ReadMeasure("@EndTime", tEnd, units);
			if (units == "s")
				tEnd = tEnd;
			else if (units == "min")
				tEnd *= 60.;
			else if (units == "h" || units == "hr")
				tEnd *= 3600.;
			else
			{
				std::cout << "Fatal error: wrong units for time. Available units: s | min | h" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		// Fixed temperature profile
		{
			if (dictionaries(dictionary_name).CheckOption("@TemperatureProfile") == true)
			{
				std::string name_of_temperature_profile_subdictionary;
				dictionaries(dictionary_name).ReadDictionary("@TemperatureProfile", name_of_temperature_profile_subdictionary);

				Eigen::VectorXd x, y;
				std::string x_variable, y_variable;
				GetXYProfileFromDictionary(dictionaries(name_of_temperature_profile_subdictionary), x, y, x_variable, y_variable);

				if (x_variable != "time")
					OpenSMOKE::FatalErrorMessage("The @TemperatureProfile must be defined versus time");
				if (y_variable != "temperature")
					OpenSMOKE::FatalErrorMessage("The @TemperatureProfile must be define the temperature profile");

				profile_t = x;
				profile_T = y;
			}
			else if (dictionaries(dictionary_name).CheckOption("@TemperatureSlope") == true)
			{
				double dT_over_dt;
				std::string units;
				dictionaries(dictionary_name).ReadMeasure("@TemperatureSlope", dT_over_dt, units);

				if (units == "C/s" || units == "K/s")
					dT_over_dt = dT_over_dt;
				else if (units == "C/min" || units == "K/min")
					dT_over_dt /= 60.;
				else if (units == "C/h" || units == "K/h")
					dT_over_dt /= 3600.;
				else
				{
					std::cout << "Fatal error: wrong units for temperature slope. Available units: C/s | K/s | C/min | K/min | C/h | K/h" << std::endl;
					std::cout << "Press enter to exit..." << std::endl;
					getchar();
					exit(-1);
				}

				const double sigma = 1.1;

				profile_t.resize(2);
				profile_t(0) = 0.;
				profile_t(1) = sigma * tEnd;

				profile_T.resize(2);
				profile_T(0) = T;
				profile_T(1) = T + dT_over_dt * profile_t(1);
			}
			else
			{
				const double sigma = 1.1;

				profile_t.resize(2);
				profile_t(0) = 0.;
				profile_t(1) = sigma * tEnd;

				profile_T.resize(2);
				profile_T(0) = T;
				profile_T(1) = T;
			}
		}

		// Options
		if (dictionaries(dictionary_name).CheckOption("@Options") == true)
		{
				std::string name_of_options_subdictionary;
				dictionaries(dictionary_name).ReadDictionary("@Options", name_of_options_subdictionary);
				thermo_options.SetupFromDictionary(dictionaries(name_of_options_subdictionary));

				if (!boost::filesystem::exists(thermo_options.output_path()))
					OpenSMOKE::CreateDirectory(thermo_options.output_path());
		}

		// ODE Parameters
		if (dictionaries(dictionary_name).CheckOption("@OdeParameters") == true)
		{
			std::string name_of_ode_parameters_subdictionary;
			dictionaries(dictionary_name).ReadDictionary("@OdeParameters", name_of_ode_parameters_subdictionary);
			ode_options.SetupFromDictionary(dictionaries(name_of_ode_parameters_subdictionary));
		}
	}
}
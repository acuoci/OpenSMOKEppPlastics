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
|   This file is part of OpenSMOKE++.                                     |
|                                                                         |
|   Copyright(C) 2019  Alberto Cuoci                                      |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#pragma once

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_ThermogravimetricAnalysis : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{

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
		}
	};


	void ThermogravimetricAnalysisFromDictionary(	const std::string input_file_name,
													const std::string dictionary_name,
													double& T, double& P, double& tEnd,
													Eigen::VectorXd& profile_t, 
													Eigen::VectorXd& profile_T		)
	{
		// Define the dictionaries
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;

		OpenSMOKE::Grammar_ThermogravimetricAnalysis grammar_thermogravimetry;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(dictionary_name).SetGrammar(grammar_thermogravimetry);


		T = 300.;
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
	}
}
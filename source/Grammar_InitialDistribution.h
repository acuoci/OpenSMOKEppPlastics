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

namespace OpenSMOKE
{

	class Grammar_InitialDistribution : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
				OpenSMOKE::SINGLE_STRING,
				"Initial distribution type: Schultz-Flory | Schultz | Dirac-Delta",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Threshold",
				OpenSMOKE::SINGLE_DOUBLE,
				"Threshold at which the distribution is cut (default: 1e-3)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StretchingFactor",
				OpenSMOKE::SINGLE_INT,
				"Stretching factor for building the initial distribution (default: 20)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MolecularWeightPolymer",
				OpenSMOKE::SINGLE_MEASURE,
				"Polymer molecular weight",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MolecularWeightMonomer",
				OpenSMOKE::SINGLE_MEASURE,
				"Monomer molecular weight (default: PS=104, PE=14. g/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Lumping",
				OpenSMOKE::SINGLE_BOOL,
				"Lumping on/off (default for PS: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LumpingStart",
				OpenSMOKE::SINGLE_INT,
				"Number of monomeric units from which lumping is applied (default: 100.)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@LumpingStep",
				OpenSMOKE::SINGLE_INT,
				"Number of monomeric units corresponding to the lumping step (default: 10.)",
				false));
		}
	};


	void InitialDistributionFromDictionary(	const std::string input_file_name,
											const std::string dictionary_name,
											Eigen::VectorXd& y, int& N,
											double& MW_polymer, double& MW_monomer, 
											bool& is_lumping_enabled, int& lumping_start, int& lumping_step, 
											const std::string& folder_name )
	{
		// Define the dictionaries
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;

		OpenSMOKE::Grammar_InitialDistribution grammar_distribution;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(dictionary_name).SetGrammar(grammar_distribution);

		opensmokepp::plastics::PolymerDistribution distribution_type = opensmokepp::plastics::SCHULTZ_FLORY;
		if (dictionaries(dictionary_name).CheckOption("@Type") == true)
		{
			std::string type;
			dictionaries(dictionary_name).ReadString("@Type", type);

			if (type == "Schultz-Flory")
				distribution_type = opensmokepp::plastics::SCHULTZ_FLORY;
			else if (type == "Schultz")
				distribution_type = opensmokepp::plastics::SCHULTZ;
			else if (type == "Dirac-Delta")
				distribution_type = opensmokepp::plastics::DIRAC_DELTA;
			else
			{
				std::cout << "Fatal error: wrong initial distribution type. Available types: Schultz-Flory | Schultz | Dirac-Delta" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		double threshold = 1.e-3;
		if (dictionaries(dictionary_name).CheckOption("@Threshold") == true)
			dictionaries(dictionary_name).ReadDouble("@Threshold", threshold);

		int stretching_factor = 20;
		if (dictionaries(dictionary_name).CheckOption("@StretchingFactor") == true)
			dictionaries(dictionary_name).ReadInt("@StretchingFactor", stretching_factor);

		is_lumping_enabled = true;
		if (dictionaries(dictionary_name).CheckOption("@Lumping") == true)
			dictionaries(dictionary_name).ReadBool("@Lumping", is_lumping_enabled);

		lumping_start = 100;
		if (dictionaries(dictionary_name).CheckOption("@LumpingStart") == true)
			dictionaries(dictionary_name).ReadInt("@LumpingStart", lumping_start);

		lumping_step = 10;
		if (dictionaries(dictionary_name).CheckOption("@LumpingStep") == true)
			dictionaries(dictionary_name).ReadInt("@LumpingStep", lumping_step);

		MW_polymer = 50000.;
		if (dictionaries(dictionary_name).CheckOption("@MolecularWeightPolymer") == true)
		{
			std::string units;
			dictionaries(dictionary_name).ReadMeasure("@MolecularWeightPolymer", MW_polymer, units);
			if (units != "g/mol" && units != "kg/kmol")
			{
				std::cout << "Fatal error: wrong units for molecular weight. Available units: g/mol | kg/kmol" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		MW_monomer = 104.;
		if (dictionaries(dictionary_name).CheckOption("@MolecularWeightMonomer") == true)
		{
			std::string units;
			dictionaries(dictionary_name).ReadMeasure("@MolecularWeightMonomer", MW_monomer, units);
			if (units != "g/mol" && units != "kg/kmol")
			{
				std::cout << "Fatal error: wrong units for molecular weight. Available units: g/mol | kg/kmol" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		// Create the initial distribution
		opensmokepp::plastics::InitialDistribution(distribution_type, MW_monomer, MW_polymer,
			lumping_start, lumping_step, threshold*100., stretching_factor, y, N, folder_name);

		// Adjust the initial distribution if lumping is enabled
		opensmokepp::plastics::LumpingSetup(is_lumping_enabled, lumping_start, lumping_step, MW_monomer, y, N, folder_name);
	}

	void InitialDistributionFromDictionary(const std::string input_file_name,
		const std::string dictionary_name,
		Eigen::VectorXd& y, int& N,
		double& MW_polymer, double& MW_monomer, double& wg, const std::string& folder_name)
	{
		// Define the dictionaries
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;

		OpenSMOKE::Grammar_InitialDistribution grammar_distribution;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(dictionary_name).SetGrammar(grammar_distribution);

		double threshold = 1.e-3;
		if (dictionaries(dictionary_name).CheckOption("@Threshold") == true)
			dictionaries(dictionary_name).ReadDouble("@Threshold", threshold);

		MW_polymer = 5000.;
		if (dictionaries(dictionary_name).CheckOption("@MolecularWeightPolymer") == true)
		{
			std::string units;
			dictionaries(dictionary_name).ReadMeasure("@MolecularWeightPolymer", MW_polymer, units);
			if (units != "g/mol" && units != "kg/kmol")
			{
				std::cout << "Fatal error: wrong units for molecular weight. Available units: g/mol | kg/kmol" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		MW_monomer = 14.;
		if (dictionaries(dictionary_name).CheckOption("@MolecularWeightMonomer") == true)
		{
			std::string units;
			dictionaries(dictionary_name).ReadMeasure("@MolecularWeightMonomer", MW_monomer, units);
			if (units != "g/mol" && units != "kg/kmol")
			{
				std::cout << "Fatal error: wrong units for molecular weight. Available units: g/mol | kg/kmol" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}
		}

		// Create the initial distribution
		opensmokepp::plastics::InitialDistribution(y, threshold, MW_monomer, MW_polymer, N, wg, folder_name);
	}

}
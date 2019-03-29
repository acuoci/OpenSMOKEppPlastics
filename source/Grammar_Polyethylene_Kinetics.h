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
#include "Utilities.h"

namespace OpenSMOKE
{

	class Grammar_Polyethylene_Kinetics : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@As",
				OpenSMOKE::SINGLE_MEASURE,
				"Random scission (PE -> Rp + Rsb): frequency factor (default: 7.94328e14 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Es",
				OpenSMOKE::SINGLE_MEASURE,
				"Random scission (PE -> Rp + Rsb): activation energy (default: 78000. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Aa",
				OpenSMOKE::SINGLE_MEASURE,
				"Allylic scission (PE -> Ra + Rsb): frequency factor (default: 4.00e13 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ea",
				OpenSMOKE::SINGLE_MEASURE,
				"Allylic scission (PE -> Ra + Rsb): activation energy (default: 68000. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@At",
				OpenSMOKE::SINGLE_MEASURE,
				"Termination (	Rsb + Rp -> PE  and  Rp + Rp -> PE): frequency factor (default: 2.5e7 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Et",
				OpenSMOKE::SINGLE_MEASURE,
				"Termination (	Rsb + Rp -> PE  and  Rp + Rp -> PE): activation energy (default: 6000. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Aaf",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstraction reaction (TODO): frequency factor (default: 4.74341649e8 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Eaf",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstraction reaction (TODO): activation energy (default: 12000. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Aab",
				OpenSMOKE::SINGLE_MEASURE,
				"H-reabstraction reaction (TODO): frequency factor (default: 4.74341649e8 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Eab",
				OpenSMOKE::SINGLE_MEASURE,
				"H-reabstraction reaction (TODO): activation energy (default: 13500. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Abeta",
				OpenSMOKE::SINGLE_MEASURE,
				"Beta-decomposition: Rt -> PE + Rsb: frequency factor (default: 1.25892541e14 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ebeta",
				OpenSMOKE::SINGLE_MEASURE,
				"Beta-decomposition: Rt -> PE + Rsb: activation energy (default: 30500.0 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Abb14",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting reaction (isomerization 1,4): TODO: frequency factor (default: 1.00e11 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ebb14",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting reaction (isomerization 1,4): TODO: activation energy (default: 20600. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Abb15",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting reaction (isomerization 1,5): TODO: frequency factor (default: 1.60e10 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ebb15",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting reaction (isomerization 1,6): TODO: activation energy (default: 14500. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Abb16",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting reaction (isomerization 1,6): TODO: frequency factor (default: 5.00e9 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ebb16",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting reaction (isomerization 1,6): TODO: activation energy (default: 14500. cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Cbutane",
				OpenSMOKE::SINGLE_DOUBLE,
				"Acceleration coefficient for butane (default: 2.)",
				false));
			
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Cpropane",
				OpenSMOKE::SINGLE_DOUBLE,
				"Acceleration coefficient for propane (default: 2.)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialAccelerationCorrection",
				OpenSMOKE::SINGLE_BOOL,
				"Initial acceleration due to ramifications and impurities (default: false)",
				false));
		}
	};


	void KineticsFromDictionary(const std::string input_file_name,
		const std::string dictionary_name,
		opensmokepp::plastics::PolyethyleneKinetics& kinetics)
	{
		// Define the dictionaries
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;

		OpenSMOKE::Grammar_Polyethylene_Kinetics grammar_kinetics;
		dictionaries.ReadDictionariesFromFile(input_file_name);
		dictionaries(dictionary_name).SetGrammar(grammar_kinetics);

		double value;
		std::string units;

		if (dictionaries(dictionary_name).CheckOption("@As") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@As", value, units);
			kinetics.SetAs(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Es") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Es", value, units);
			kinetics.SetEs(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@At") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@At", value, units);
			kinetics.SetAt(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Et") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Et", value, units);
			kinetics.SetEt(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Abeta") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Abeta", value, units);
			kinetics.SetAbeta(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ebeta") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ebeta", value, units);
			kinetics.SetEbeta(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Aaf") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Aaf", value, units);
			kinetics.SetAaf(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Eaf") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Eaf", value, units);
			kinetics.SetEaf(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Aab") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Aab", value, units);
			kinetics.SetAab(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Eab") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Eab", value, units);
			kinetics.SetEab(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Abb14") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Abb14", value, units);
			kinetics.SetAbb14(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ebb14") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ebb14", value, units);
			kinetics.SetEbb14(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Abb15") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Abb15", value, units);
			kinetics.SetAbb15(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ebb15") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ebb15", value, units);
			kinetics.SetEbb15(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Abb16") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Abb16", value, units);
			kinetics.SetAbb16(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ebb16") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ebb16", value, units);
			kinetics.SetEbb16(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Cbutane") == true)
		{
			dictionaries(dictionary_name).ReadDouble("@Cbutane", value);
			kinetics.SetCbutane(value);
		}

		if (dictionaries(dictionary_name).CheckOption("@Cpropane") == true)
		{
			dictionaries(dictionary_name).ReadDouble("@Cpropane", value);
			kinetics.SetCpropane(value);
		}

		if (dictionaries(dictionary_name).CheckOption("@InitialAccelerationCorrection") == true)
		{
			bool flag = false;
			dictionaries(dictionary_name).ReadBool("@InitialAccelerationCorrection", flag);
			kinetics.SetInitialAccelerationFactor(flag);
		}

	}
}
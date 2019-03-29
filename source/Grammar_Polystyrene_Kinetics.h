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
	
	class Grammar_Polystyrene_Kinetics : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@As",
				OpenSMOKE::SINGLE_MEASURE,
				"Random scission (PS -> Rp + Rsb): frequency factor (default: 5e13 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Es",
				OpenSMOKE::SINGLE_MEASURE,
				"Random scission (PS -> Rp + Rsb): activation energy (default: 67500 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Aa",
				OpenSMOKE::SINGLE_MEASURE,
				"Allylic scission (PS -> Ra + Rsb): frequency factor (default: 5e12 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ea",
				OpenSMOKE::SINGLE_MEASURE,
				"Allylic scission (PS -> Ra + Rsb): activation energy (default: 62500 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@At",
				OpenSMOKE::SINGLE_MEASURE,
				"Termination (	Rsb + Rp -> PS  and  Rp + Rp -> PS): frequency factor (default: 5e6 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Et",
				OpenSMOKE::SINGLE_MEASURE,
				"Termination (	Rsb + Rp -> PS  and  Rp + Rp -> PS): activation energy (default: 14000 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ar",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rp + PS -> PS + Rt: frequency factor (default: 5e7 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Er",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rp + PS -> PS + Rt: activation energy (default: 16500 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ar1",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rp + PS -> PS + Rt: frequency factor (default: 5e7 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Er1",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rp + PS -> PS + Rt: activation energy (default: 8500 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ar2",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rsb + PS -> PS + Rt: frequency factor (default: 5e7 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Er2",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rsb + PS -> PS + Rt: activation energy (default: 13500 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ar3",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rsb + PS -> PS + Rt: frequency factor (default: 5e7 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Er3",
				OpenSMOKE::SINGLE_MEASURE,
				"H-abstractions: Rsb + PS -> PS + Rt: activation energy (default: 16500 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Au",
				OpenSMOKE::SINGLE_MEASURE,
				"Unzipping: Rsb -> Rsb + S: frequency factor (default: 1e13 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Eu",
				OpenSMOKE::SINGLE_MEASURE,
				"Unzipping: Rsb -> Rsb + S: activation energy (default: 26000 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Abb",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting (intramolecular abs): Rsb -> Rt  and  Rp -> Rt: frequency factor (default: 1e9 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ebb",
				OpenSMOKE::SINGLE_MEASURE,
				"Back-biting (intramolecular abs): Rsb -> Rt  and  Rp -> Rt: activation energy (default: 16000 cal/mol)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Abeta",
				OpenSMOKE::SINGLE_MEASURE,
				"Beta-decomposition: Rt -> PS + Rsb: frequency factor (default: 1e13 mol,l,s)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Ebeta",
				OpenSMOKE::SINGLE_MEASURE,
				"Beta-decomposition: Rt -> PS + Rsb: activation energy (default: 27000 cal/mol)",
				false));


			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RightSideBetaScissions",
				OpenSMOKE::SINGLE_BOOL,
				"Correction for right-side beta scissions (default: false)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@BackBiting",
				OpenSMOKE::SINGLE_BOOL,
				"Back-biting reactions on/off (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RandomScissionEfficiency",
				OpenSMOKE::SINGLE_DOUBLE,
				"Effective fraction of radicals from random scission (default: 0.01)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@AbstractionReactionWeights",
				OpenSMOKE::VECTOR_DOUBLE,
				"Weights for effective kinetic constant of abstraction reactions (default: 0. 1. 0.)",
				false));

		}
	};

	


	void KineticsFromDictionary(	const std::string input_file_name,
									const std::string dictionary_name,
									opensmokepp::plastics::PolystyreneKinetics& kinetics )
	{
		// Define the dictionaries
		OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;

		OpenSMOKE::Grammar_Polystyrene_Kinetics grammar_kinetics;
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

		if (dictionaries(dictionary_name).CheckOption("@Aa") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Aa", value, units);
			kinetics.SetAa(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ea") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ea", value, units);
			kinetics.SetEa(opensmokepp::plastics::ConversionActivationEnergy(value, units));
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

		if (dictionaries(dictionary_name).CheckOption("@Ar") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ar", value, units);
			kinetics.SetAr(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Er") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Er", value, units);
			kinetics.SetEr(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ar1") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ar1", value, units);
			kinetics.SetAr1(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Er1") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Er1", value, units);
			kinetics.SetEr1(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ar2") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ar2", value, units);
			kinetics.SetAr2(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Er2") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Er2", value, units);
			kinetics.SetEr2(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ar3") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ar3", value, units);
			kinetics.SetAr3(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Er3") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Er3", value, units);
			kinetics.SetEr3(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Au") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Au", value, units);
			kinetics.SetAu(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Eu") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Eu", value, units);
			kinetics.SetEu(opensmokepp::plastics::ConversionActivationEnergy(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Abb") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Abb", value, units);
			kinetics.SetAbb(opensmokepp::plastics::ConversionFrequencyFactor(value, units));
		}

		if (dictionaries(dictionary_name).CheckOption("@Ebb") == true)
		{
			dictionaries(dictionary_name).ReadMeasure("@Ebb", value, units);
			kinetics.SetEbb(opensmokepp::plastics::ConversionActivationEnergy(value, units));
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


		if (dictionaries(dictionary_name).CheckOption("@RightSideBetaScissions") == true)
		{
			bool flag;
			dictionaries(dictionary_name).ReadBool("@RightSideBetaScissions", flag);
			kinetics.SetRightSideBetaScissions(flag);
		}

		if (dictionaries(dictionary_name).CheckOption("@BackBiting") == true)
		{
			bool flag;
			dictionaries(dictionary_name).ReadBool("@BackBiting", flag);
			kinetics.SetBackBiting(flag);
		}

		if (dictionaries(dictionary_name).CheckOption("@RandomScissionEfficiency") == true)
		{
			dictionaries(dictionary_name).ReadDouble("@RandomScissionEfficiency", value);
			kinetics.SetRandomScissionEfficiency(value);
		}

		if (dictionaries(dictionary_name).CheckOption("@AbstractionReactionWeights") == true)
		{
			std::vector<double> values;
			dictionaries(dictionary_name).ReadOption("@AbstractionReactionWeights", values);
			kinetics.SetAbstractionReactionWeights(values);
		}

	}
}
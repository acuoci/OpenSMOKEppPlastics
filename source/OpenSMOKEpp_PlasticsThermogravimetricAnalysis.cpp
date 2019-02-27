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

// Eigen C++ library
#include <Eigen/Dense>

// Utilities from OpenSMOKE++ framework
#include "idealreactors/utilities/Utilities"

// OpenSMOKE++ Plastics framework
#include "Utilities.h"
#include "PolystyreneKinetics.h"

// Grammar files
#include "Grammar_Kinetics.h"
#include "Grammar_InitialDistribution.h"
#include "Grammar_ThermogravimetricAnalysis.h"

// Folder containing data
std::string test_folder = "../../../../../data/";

int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_PlasticsThermogravimetricAnalysis", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	//unsigned int max_number_allowed_species = 100000;
	//OpenSMOKE::OpenSMOKE_CheckLicense(executable_folder, "PlugFlowReactor", max_number_allowed_species);

	std::string input_file_name = test_folder + "input.dic";

	std::string thermogravimetricanalysis_dictionary_name = "ThermogravimetricAnalysis";
	std::string initial_distribution_dictionary_name = "InitialDistribution";
	std::string kinetics_dictionary_name = "Kinetics";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKEpp_PlasticsThermogravimetricAnalysis");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"ThermogravimetricAnalysis\")");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return OPENSMOKE_SUCCESSFULL_EXIT;
			}

			if (vm.count("input"))
				input_file_name = vm["input"].as<std::string>();

			if (vm.count("dictionary"))
				thermogravimetricanalysis_dictionary_name = vm["dictionary"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}
	}

	// Thermogravimetric analysis from input file
	double T, P, tEnd;
	Eigen::VectorXd profile_t, profile_T;
	OpenSMOKE::ThermogravimetricAnalysisFromDictionary(input_file_name, thermogravimetricanalysis_dictionary_name, T, P, tEnd, profile_t, profile_T);
	OpenSMOKE::FixedProfile temperature_profile(profile_t.size(), profile_t.data(), profile_T.data());
	T = temperature_profile.Interpolate(0.);


	// Initial distribution from input file
	Eigen::VectorXd y;
	int N = 0;
	double MW_polymer, MW_monomer;
	bool is_lumping_enabled;
	int lumping_start, lumping_step;
	OpenSMOKE::InitialDistributionFromDictionary(input_file_name, initial_distribution_dictionary_name, 
		y, N, MW_polymer, MW_monomer, is_lumping_enabled, lumping_start, lumping_step);

	// Define polystyrene kinetics
	opensmokepp::plastics::PolystyreneKinetics kinetics;
	kinetics.SetVerbose(false);
	kinetics.SetMaxNumberOfUnits(N);
	kinetics.SetMolecularWeightPolymer(MW_polymer);
	kinetics.SetMolecularWeightMonomer(MW_monomer);
	kinetics.SetLumping(is_lumping_enabled);
	kinetics.SetLumpingStart(lumping_start);
	kinetics.SetLumpingStep(lumping_step);

	// Additional options for polystyrene kinetics
	OpenSMOKE::KineticsFromDictionary(input_file_name, kinetics_dictionary_name, kinetics);

	

	kinetics.UpdateSharedSpecies(T, P);
	kinetics.UpdateSharedSpeciesDistribution(T, P, y);
	std::cout << "Sum: " << y.sum() << " moles" << std::endl;


	const double rho = kinetics.LiquidDensity(T);
	const double Ctot = rho / kinetics.MolecularWeightPolymer() * 1.001059912;
	Eigen::VectorXd x = y / y.sum();				// mole fractions
	Eigen::VectorXd c = Ctot * x;					// concentrations (in mol/l or kmol/m3)

	kinetics.SetStatus(T, P, c);

	const double m_tot = 0.100;
	const double rho_gas = kinetics.SumGasMW(c);				// gaseous density [kg/m3]
	const double rho_liq = kinetics.SumLiquidMW(c);				// liquid density [kg/m3]
	const double m_gas = m_tot * rho_gas / (rho_gas + rho_liq); // gaseous mass [kg]
	const double m_liq = m_tot * rho_liq / (rho_gas + rho_liq); // liquid mass [kg]
	const double V_liq = m_tot / rho_liq;						// volume of liquid phase [m3]


	Eigen::VectorXd n = c * V_liq;			// moles (in kmol)
	n *= 1000.;								// moles (in mol)

	std::cout << rho_liq << " " << rho << " " << rho_liq - rho << std::endl;
	std::cout << m_gas << " " << m_liq << std::endl;
	std::cout << "Sum: " << y.sum() << " moles" << std::endl;
	std::cout << "VL:  " << V_liq << " m3" << std::endl;
	std::cout << "CL:  " << rho / kinetics.MolecularWeightPolymer() * 1.001059912 << " kmol/m3" << std::endl;
	std::cout << "MWliq: " << kinetics.SumLiquidMW(x) << std::endl;
	std::cout << "MWliq: " << kinetics.MolecularWeightPolymer() / 1.001059912 << std::endl;
	std::cout << "MWliq: " << kinetics.MolecularWeightPolymer() * (rho_liq / (rho_liq + rho_gas)) << std::endl;


	const double CCsr = 1. / (V_liq*1000.);	// concentration (in mol/l or kmol/m3)
	const double CCsa = 0. / (V_liq*1000.);	// concentration (in mol/l or kmol/m3)
	kinetics.SetCCBonds(CCsr, CCsa);
	kinetics.KineticConstants();
	kinetics.FormationRates();

	bool integrate = true;
	if (integrate == true)
	{
		std::ofstream fOut("History.out", std::ios::out);
		fOut.setf(std::ios::scientific);

		double t = 0.;
		const double dt = 1.e-4;
		const int nsteps = static_cast<int>(tEnd / dt);
		for (int k = 0; k <= nsteps; k++)
		{
			t += dt;

			const double nG = kinetics.SumGas(n);		// mol
			const double nL = kinetics.SumLiquid(n);	// mol
			const double mG = kinetics.SumGasMW(n);		// g
			const double mL = kinetics.SumLiquidMW(n);	// g
			const double VL = (mL / 1000.) / rho;		// m3

			x = n / n.sum();
			c = n / (VL*1000.);

			T = temperature_profile.Interpolate(t);
			kinetics.SetStatus(T, P, c);
			kinetics.UpdateCCBonds();
			kinetics.KineticConstants();
			kinetics.FormationRates();
			n += dt * kinetics.R()*VL*1000.;			// moles
			for (int i = 0; i < n.size(); i++)
				if (n(i) < 0.) n(i) = 0.;

			if (k % 100 == 1)
				std::cout << k << " " << t << " " << T << " " << mG << " " << mL << " " << mG + mL << std::endl;

			if (k % 100 == 1)
				fOut << k << " " << t << " " << mG << " " << mL << " " << mG + mL << std::endl;
		}

		fOut.close();
	}

	return 1;
}


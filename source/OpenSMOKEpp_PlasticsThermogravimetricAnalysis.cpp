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

// MKL Support Enabled
#define MKL_SUPPORT 1
#define EIGEN_USE_MKL_ALL
#define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_CXX17_OLD_ALLOCATOR_MEMBERS_DEPRECATION_WARNING
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING

// Eigen C++ library
#include <Eigen/Dense>

// Boost C++ Libraries
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

// Utilities from OpenSMOKE++ framework
#include "idealreactors/utilities/Utilities"

// OpenSMOKE++ Plastics framework
#include "Utilities.h"
#include "PolystyreneKinetics.h"
#include "PolyethyleneKinetics.h"
#include "PolypropyleneKinetics.h"
#include "ThermogravimetricAnalysis.h"

// Grammar files
#include "Grammar_Polyethylene_Kinetics.h"
#include "Grammar_Polypropylene_Kinetics.h"
#include "Grammar_Polystyrene_Kinetics.h"
#include "Grammar_InitialDistribution.h"

// OdeSMOKE Solvers
#include "ode_solver_virtual_class.h"
#include "runge_kutta_4th.h"
#include "runge_kutta_fehlberg.h"

// Pointers to kinetics
opensmokepp::plastics::PolystyreneKinetics*	ptPS;
opensmokepp::plastics::PolyethyleneKinetics* ptPE;
opensmokepp::plastics::PolypropyleneKinetics* ptPP;

// Folder containing data
std::string test_folder = "../../../../../data/";

#include "Polystyrene_ODESystem.h"
#include "Polyethylene_ODESystem.h"
#include "Polypropylene_ODESystem.h"

ThermogravimetricAnalysis* tg;

int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_PlasticsThermogravimetricAnalysis", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	//unsigned int max_number_allowed_species = 100000;
	//OpenSMOKE::OpenSMOKE_CheckLicense(executable_folder, "OpenSMOKEpp_PlasticsThermogravimetricAnalysis", max_number_allowed_species);

	std::string input_file_name = "input.dic";

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
	ThermogravimetricAnalysis thermo_gravimetric_analysis;
	tg = &thermo_gravimetric_analysis;
	tg->operator()(input_file_name, thermogravimetricanalysis_dictionary_name);

	if (tg->Polymer() == OpenSMOKE::POLYETHYLENE)
	{
		// Initial distribution from input file
		Eigen::VectorXd y(600001);
		double wg = 1.;
		int N = 0;
		double MW_polymer, MW_monomer;
		OpenSMOKE::InitialDistributionFromDictionary(input_file_name, initial_distribution_dictionary_name,
			y, N, MW_polymer, MW_monomer, wg, tg->options().output_path().string());


		tg->SetInitialMass(wg);

		// Define polystyrene kinetics
		opensmokepp::plastics::PolyethyleneKinetics kinetics;
		ptPE = &kinetics;
		kinetics.SetVerbose(false);
		kinetics.SetMaxNumberOfUnits(N);
		kinetics.SetMolecularWeightPolymer(MW_polymer);
		kinetics.SetMolecularWeightMonomer(MW_monomer);

		// Boiling temperature table
		std::cout << std::endl;
		std::cout << "------------------- -------------" << std::endl;
		std::cout << "   Boiling Temperature Table     " << std::endl;
		std::cout << "---------------------------------" << std::endl;
		std::cout << "   #LC        T(K)        T(C)   " << std::endl;
		std::cout << "---------------------------------" << std::endl;
		std::vector<double> list_boiling_temperature;
		for (int i = 1; i <= 100; i++)
		{
			const double Teb = kinetics.BoilingTemperature(i, tg->P(0.));
			list_boiling_temperature.push_back(Teb);
			std::cout << std::setw(6) << std::right << i;
			std::cout << std::setw(12) << std::right << std::setprecision(3) << std::fixed << Teb;
			std::cout << std::setw(12) << std::right << std::setprecision(3) << std::fixed << Teb - 273.15;
			std::cout << std::endl;
		}
		tg->SetBoilingTemperatureList(list_boiling_temperature);
		std::cout << "---------------------------------------------------------------" << std::endl;
		std::cout << std::endl;

		// Additional options for polystyrene kinetics
		OpenSMOKE::KineticsFromDictionary(input_file_name, kinetics_dictionary_name, kinetics);

		const int rho = kinetics.LiquidDensity(tg->T(0));
		const int LC = kinetics.MinNumberOfUnits(tg->T(0), tg->P(0));
		const double Teb = kinetics.BoilingTemperature(LC, tg->P(0));

		Eigen::VectorXd x = y / y.sum();
		const double Ctot = kinetics.LiquidDensity(tg->T(0)) / kinetics.MolecularWeightMonomer();	// (kmol/m3 or mol/l)
		Eigen::VectorXd c = Ctot * x;	// (kmol/m3 or mol/l)

		kinetics.SetMinNumberOfUnits(LC);

		Eigen::VectorXd n = y;

		const double rhoG = kinetics.SumGasMW(c);							// Gaseous density [kg/m3]
		const double rhoL = kinetics.LiquidDensity(tg->T(0));				// Liquid density [kg/m3]
		const double nG = kinetics.SumGas(n);								// Gaseous-phase moles
		const double nL = kinetics.SumLiquid(n);							// Liquid-phase moles [mol]
		const double mG = kinetics.SumGasMW(n);								// Gaseous-phase mass [g]
		const double mL = kinetics.SumLiquidMW(n);							// Liquid-phase mass [g]
		const double mTot = mG + mL;										// Total mass [g]
		const double VL = (mL / 1000.) / kinetics.LiquidDensity(tg->T(0));	// Liquid-phase volume [m3]
		const double Vtot = (mTot / 1000.) / rho;							// Total volume [m3]

		kinetics.SetStatus(tg->T(0), tg->P(0), c);
		kinetics.UpdateInitialAccelerationCoefficient(VL*1000., wg);
		kinetics.KineticConstants();
		kinetics.FormationRates();

		std::cout << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << " Summary                                                             " << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << " * Polymer MW [g/mol]:    " << MW_polymer << std::endl;
		std::cout << " * Monomer MW [g/mol]:    " << MW_monomer << std::endl;
		std::cout << " * Density [kg/m3]:       " << rho << std::endl;
		std::cout << " * Concentration [mol/l]: " << Ctot << std::endl;
		std::cout << " * Gas mass fraction:     " << rhoG / (rhoG + rhoL) << std::endl;
		std::cout << " * Liquid mass fraction:  " << rhoL / (rhoG + rhoL) << std::endl;
		std::cout << " * Total mass [g]:        " << mTot << std::endl;
		std::cout << " * Total moles [mol]:     " << y.sum() << std::endl;
		std::cout << " * Total volume [cm3]:    " << Vtot * 1.e6 << std::endl;
		std::cout << " * Temperature [K]:       " << tg->T(0) << std::endl;
		std::cout << " * Temperature [C]:       " << tg->T(0) - 273.15 << std::endl;
		std::cout << " * Pressure [atm]:        " << tg->P(0) << std::endl;
		std::cout << " * LC [-]:                " << LC << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;

		// Integration of ODE system
		if (tg->ode_options().type() == OdeSMOKE::ExplicitOdeSolver_Parameters::ODE_RUNGEKUTTA45)
		{
			typedef PE_ODESystemObject ode_system;
			typedef OdeSMOKE::OdeRungeKutta4th< ode_system, Eigen::VectorXd > Method_RungeKutta4thOrder;
			typedef OdeSMOKE::OdeRungeKuttaFamily< Eigen::VectorXd, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
			OdeSMOKE::ODESolverVirtualClass< Eigen::VectorXd, Solver_RungeKutta4thOrder > ode_solver_rk;

			ode_solver_rk.PrepareOutputFiles(tg->options().output_path());
			ode_solver_rk.SetThermogravimetricAnalysis(tg);
			ode_solver_rk.SetInitialConditions(n);
			ode_solver_rk.Solve(0, tg->FinalTime());
			n = ode_solver_rk.y();
		}
		else
		{
		}
	}

	else if (tg->Polymer() == OpenSMOKE::POLYPROPYLENE)
	{
		// Initial distribution from input file
		Eigen::VectorXd y(600001);
		double wg = 1.;
		int N = 0;
		double MW_polymer, MW_monomer;
		OpenSMOKE::InitialDistributionFromDictionary(input_file_name, initial_distribution_dictionary_name,
			y, N, MW_polymer, MW_monomer, wg, tg->options().output_path().string());


		tg->SetInitialMass(wg);

		// Define polypropylene kinetics
		opensmokepp::plastics::PolypropyleneKinetics kinetics;
		ptPP = &kinetics;
		kinetics.SetVerbose(false);
		kinetics.SetMaxNumberOfUnits(N);
		kinetics.SetMolecularWeightPolymer(MW_polymer);
		kinetics.SetMolecularWeightMonomer(MW_monomer);

		// Boiling temperature table
		std::cout << std::endl;
		std::cout << "------------------- -------------" << std::endl;
		std::cout << "   Boiling Temperature Table     " << std::endl;
		std::cout << "---------------------------------" << std::endl;
		std::cout << "   #LC        T(K)        T(C)   " << std::endl;
		std::cout << "---------------------------------" << std::endl;
		std::vector<double> list_boiling_temperature;
		for (int i = 1; i <= 100; i++)
		{
			const double Teb = kinetics.BoilingTemperature(i, tg->P(0.));
			list_boiling_temperature.push_back(Teb);
			std::cout << std::setw(6) << std::right << i;
			std::cout << std::setw(12) << std::right << std::setprecision(3) << std::fixed << Teb;
			std::cout << std::setw(12) << std::right << std::setprecision(3) << std::fixed << Teb - 273.15;
			std::cout << std::endl;
		}
		tg->SetBoilingTemperatureList(list_boiling_temperature);
		std::cout << "---------------------------------------------------------------" << std::endl;
		std::cout << std::endl;

		// Additional options for polystyrene kinetics
		OpenSMOKE::KineticsFromDictionary(input_file_name, kinetics_dictionary_name, kinetics);

		const int rho = kinetics.LiquidDensity(tg->T(0));
		const int LC = kinetics.MinNumberOfUnits(tg->T(0), tg->P(0));
		const double Teb = kinetics.BoilingTemperature(LC, tg->P(0));

		Eigen::VectorXd x = y / y.sum();
		const double Ctot = kinetics.LiquidDensity(tg->T(0)) / kinetics.MolecularWeightMonomer();	// (kmol/m3 or mol/l)
		Eigen::VectorXd c = Ctot * x;	// (kmol/m3 or mol/l)

		kinetics.SetMinNumberOfUnits(LC);

		Eigen::VectorXd n = y;

		const double rhoG = kinetics.SumGasMW(c);							// Gaseous density [kg/m3]
		const double rhoL = kinetics.LiquidDensity(tg->T(0));				// Liquid density [kg/m3]
		const double nG = kinetics.SumGas(n);								// Gaseous-phase moles
		const double nL = kinetics.SumLiquid(n);							// Liquid-phase moles [mol]
		const double mG = kinetics.SumGasMW(n);								// Gaseous-phase mass [g]
		const double mL = kinetics.SumLiquidMW(n);							// Liquid-phase mass [g]
		const double mTot = mG + mL;										// Total mass [g]
		const double VL = (mL / 1000.) / kinetics.LiquidDensity(tg->T(0));	// Liquid-phase volume [m3]
		const double Vtot = (mTot / 1000.) / rho;							// Total volume [m3]

		kinetics.SetStatus(tg->T(0), tg->P(0), c);
		kinetics.UpdateInitialAccelerationCoefficient(VL*1000., wg);
		kinetics.KineticConstants();
		kinetics.FormationRates();

		std::cout << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << " Summary                                                             " << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << " * Polymer MW [g/mol]:    " << MW_polymer << std::endl;
		std::cout << " * Monomer MW [g/mol]:    " << MW_monomer << std::endl;
		std::cout << " * Density [kg/m3]:       " << rho << std::endl;
		std::cout << " * Concentration [mol/l]: " << Ctot << std::endl;
		std::cout << " * Gas mass fraction:     " << rhoG / (rhoG + rhoL) << std::endl;
		std::cout << " * Liquid mass fraction:  " << rhoL / (rhoG + rhoL) << std::endl;
		std::cout << " * Total mass [g]:        " << mTot << std::endl;
		std::cout << " * Total moles [mol]:     " << y.sum() << std::endl;
		std::cout << " * Total volume [cm3]:    " << Vtot * 1.e6 << std::endl;
		std::cout << " * Temperature [K]:       " << tg->T(0) << std::endl;
		std::cout << " * Temperature [C]:       " << tg->T(0) - 273.15 << std::endl;
		std::cout << " * Pressure [atm]:        " << tg->P(0) << std::endl;
		std::cout << " * LC [-]:                " << LC << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;

		// Integration of ODE system
		if (tg->ode_options().type() == OdeSMOKE::ExplicitOdeSolver_Parameters::ODE_RUNGEKUTTA45)
		{
			typedef PP_ODESystemObject ode_system;
			typedef OdeSMOKE::OdeRungeKutta4th< ode_system, Eigen::VectorXd > Method_RungeKutta4thOrder;
			typedef OdeSMOKE::OdeRungeKuttaFamily< Eigen::VectorXd, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
			OdeSMOKE::ODESolverVirtualClass< Eigen::VectorXd, Solver_RungeKutta4thOrder > ode_solver_rk;

			ode_solver_rk.PrepareOutputFiles(tg->options().output_path());
			ode_solver_rk.SetThermogravimetricAnalysis(tg);
			ode_solver_rk.SetInitialConditions(n);
			ode_solver_rk.Solve(0, tg->FinalTime());
			n = ode_solver_rk.y();
		}
		else
		{
		}
	}

	else if (tg->Polymer() == OpenSMOKE::POLYSTYRENE)
	{
		// Initial distribution from input file
		Eigen::VectorXd y;
		int N = 0;
		double MW_polymer, MW_monomer;
		bool is_lumping_enabled;
		int lumping_start, lumping_step;
		OpenSMOKE::InitialDistributionFromDictionary(input_file_name, initial_distribution_dictionary_name,
							y, N, MW_polymer, MW_monomer, is_lumping_enabled, lumping_start, lumping_step, tg->options().output_path().string());

		// Set initial mass
		tg->SetInitialMass(100.);

		// Define polystyrene kinetics
		opensmokepp::plastics::PolystyreneKinetics kinetics;
		ptPS = &kinetics;
		kinetics.SetVerbose(false);
		kinetics.SetMaxNumberOfUnits(N);
		kinetics.SetMolecularWeightPolymer(MW_polymer);
		kinetics.SetMolecularWeightMonomer(MW_monomer);
		kinetics.SetLumping(is_lumping_enabled);
		kinetics.SetLumpingStart(lumping_start);
		kinetics.SetLumpingStep(lumping_step);

		// Boiling temperature table
		std::cout << std::endl;
		std::cout << "------------------- -------------" << std::endl;
		std::cout << "   Boiling Temperature Table     " << std::endl;
		std::cout << "---------------------------------" << std::endl;
		std::cout << "   #LC        T(K)        T(C)   " << std::endl;
		std::cout << "---------------------------------" << std::endl;
		std::vector<double> list_boiling_temperature;
		for (int i = 1; i <= 50; i++)
		{
			const double Teb = kinetics.BoilingTemperature(i, tg->P(0.));
			list_boiling_temperature.push_back(Teb);
			std::cout << std::setw(6) << std::right << i;
			std::cout << std::setw(12) << std::right << std::setprecision(3) << std::fixed << Teb;
			std::cout << std::setw(12) << std::right << std::setprecision(3) << std::fixed << Teb - 273.15;
			std::cout << std::endl;
		}
		tg->SetBoilingTemperatureList(list_boiling_temperature);
		std::cout << "---------------------------------------------------------------" << std::endl;
		std::cout << std::endl;


		// Additional options for polystyrene kinetics
		OpenSMOKE::KineticsFromDictionary(input_file_name, kinetics_dictionary_name, kinetics);

		// Complete the definition of initial distribution
		kinetics.UpdateSharedSpecies(tg->T(0), tg->P(0));
		kinetics.UpdateSharedSpeciesDistribution(tg->T(0), tg->P(0), y);

		// Calculate relevant data for thermogravimetric analysis
		Eigen::VectorXd x = y / y.sum();						// mole fractions
		const double rho = kinetics.LiquidDensity(tg->T(0));	// density [kg/m3]
		const double Ctot = rho / kinetics.SumMW(x);			// concentration (in mol/l or kmol/m3)
		Eigen::VectorXd c = Ctot * x;							// concentrations (in mol/l or kmol/m3)
		const double mTot = 0.100;								// total mass [kg]
		const double rhoG = kinetics.SumGasMW(c);				// gaseous density [kg/m3]
		const double rhoL = kinetics.SumLiquidMW(c);			// liquid density [kg/m3]
		const double mG = mTot * rhoG / (rhoG + rhoL);			// gaseous mass [kg]
		const double mL = mTot * rhoL / (rhoG + rhoL);			// liquid mass [kg]
		const double VL = mL / rho;								// volume of liquid phase [m3]
		const double Vtot = mTot / rho;							// volume of liquid phase [m3]

		std::cout << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << " Summary                                                             " << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << " * Polymer MW [g/mol]:    " << MW_polymer << std::endl;
		std::cout << " * Monomer MW [g/mol]:    " << MW_monomer << std::endl;
		std::cout << " * Density [kg/m3]:       " << rho << std::endl;
		std::cout << " * Concentration [mol/l]: " << Ctot << std::endl;
		std::cout << " * Gas mass fraction:     " << rhoG / (rhoG + rhoL) << std::endl;
		std::cout << " * Liquid mass fraction:  " << rhoL / (rhoG + rhoL) << std::endl;
		std::cout << " * Total mass [g]:        " << mTot * 1000. << std::endl;
		std::cout << " * Total moles [mol]:     " << y.sum() << std::endl;
		std::cout << " * Total volume [cm3]:    " << Vtot * 1.e6 << std::endl;
		std::cout << "---------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;


		const double CCsr = 1. / (VL*1000.);	// concentration (in mol/l or kmol/m3)
		const double CCsa = 0. / (VL*1000.);	// concentration (in mol/l or kmol/m3)

		kinetics.SetStatus(tg->T(0), tg->P(0), c);
		kinetics.SetCCBonds(CCsr, CCsa);
		kinetics.KineticConstants();
		kinetics.FormationRates();

		Eigen::VectorXd n = c * Vtot;			// moles (in kmol)
		n *= 1000.;

		if (tg->ode_options().type() == OdeSMOKE::ExplicitOdeSolver_Parameters::ODE_RUNGEKUTTA45)
		{
			typedef PS_ODESystemObject ode_system;
			typedef OdeSMOKE::OdeRungeKutta4th< ode_system, Eigen::VectorXd > Method_RungeKutta4thOrder;
			typedef OdeSMOKE::OdeRungeKuttaFamily< Eigen::VectorXd, Method_RungeKutta4thOrder> Solver_RungeKutta4thOrder;
			OdeSMOKE::ODESolverVirtualClass< Eigen::VectorXd, Solver_RungeKutta4thOrder > ode_solver_rk;

			ode_solver_rk.PrepareOutputFiles(tg->options().output_path());
			ode_solver_rk.SetThermogravimetricAnalysis(tg);
			ode_solver_rk.SetAbsoluteTolerances(tg->ode_options().absolute_tolerance());
			ode_solver_rk.SetRelativeTolerances(tg->ode_options().relative_tolerance());
			ode_solver_rk.SetInitialConditions(n);
			ode_solver_rk.Solve(0, tg->FinalTime());
			n = ode_solver_rk.y();
		}
	}

	return 1;
}



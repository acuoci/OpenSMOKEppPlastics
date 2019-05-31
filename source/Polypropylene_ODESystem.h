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

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "PolypropyleneKinetics.h"

extern opensmokepp::plastics::PolypropyleneKinetics* ptPP;

class PP_ODESystemObject
{
protected:

	unsigned int number_of_equations_;

public:

	PP_ODESystemObject()
	{
		number_of_equations_ = 3 * ptPP->MaxNumberOfUnits() + 34;
		step_ = 0;
	}

	void SetThermogravimetricAnalysis(ThermogravimetricAnalysis* tg)
	{
		tg_ = tg;
	}

	void PrepareOutputFiles(const boost::filesystem::path& folder_name)
	{
		boost::filesystem::path file_name_mass = folder_name / "GasDistributionMass.out";

		fDistMass_.open(file_name_mass.c_str(), std::ios::out);
		fDistMass_.setf(std::ios::scientific);
		fDistMass_ << std::setprecision(6);
		fDistMass_ << std::left;

		fDistMass_ << std::setw(7)  << "#(1)";
		fDistMass_ << std::setw(16) << "time[s](2)";
		fDistMass_ << std::setw(16) << "T[K](3)";
		fDistMass_ << std::setw(16) << "LC[-](4)";
		fDistMass_ << std::setw(16) << "dummy(5)";
		fDistMass_ << std::setw(16) << "res_liq(6)";
		fDistMass_ << std::setw(16) << "mass_liq[g](7)";
		fDistMass_ << std::setw(16) << "mass_gas[g](8)";
		fDistMass_ << std::setw(16) << "mass_tot[g](9)";
		fDistMass_ << std::setw(16) << "P_w(10)";
		fDistMass_ << std::setw(16) << "O_w(11)";
		fDistMass_ << std::setw(16) << "D_w(12)";
		fDistMass_ << std::endl;



		boost::filesystem::path file_name_moles = folder_name / "GasDistributionMoles.out";

		fDistMoles_.open(file_name_moles.c_str(), std::ios::out);
		fDistMoles_.setf(std::ios::scientific);
		fDistMoles_ << std::setprecision(6);
		fDistMoles_ << std::left;

		fDistMoles_ << std::setw(7) << "#(1)";
		fDistMoles_ << std::setw(16) << "time[s](2)";
		fDistMoles_ << std::setw(16) << "T[K](3)";
		fDistMoles_ << std::setw(16) << "LC[-](4)";
		fDistMoles_ << std::setw(16) << "dummy(5)";
		fDistMoles_ << std::setw(16) << "res_liq(6)";
		fDistMoles_ << std::setw(16) << "mol_liq[mol](7)";
		fDistMoles_ << std::setw(16) << "mol_gas[mol](8)";
		fDistMoles_ << std::setw(16) << "mol_tot[mol](9)";
		fDistMoles_ << std::setw(16) << "P_x(10)";
		fDistMoles_ << std::setw(16) << "O_x(11)";
		fDistMoles_ << std::setw(16) << "D_x(12)";
		fDistMoles_ << std::endl;
	}

	void CloseOutputFiles()
	{
		fDistMass_.close();
		fDistMoles_.close();
	}

	void GetFunctions(const double t, const Eigen::VectorXd& n, Eigen::VectorXd& dn_over_dx)
	{
		const double T = tg_->T(t);
		const double P = tg_->P(t);
		const double W = tg_->InitialMass();

		const double nG = ptPP->SumGas(n);							// mol
		const double nL = ptPP->SumLiquid(n);						// mol
		const double mG = ptPP->SumGasMW(n);						// g
		const double mL = ptPP->SumLiquidMW(n);						// g
		const double VL = (mL / 1000.) / ptPP->LiquidDensity(T);	// m3

		c_ = n / (VL*1000.);		// concentrations (in kmol/m3 or mol/l)

		ptPP->SetStatus(T, P, c_);
		ptPP->UpdateInitialAccelerationCoefficient(VL*1000., W);
		ptPP->KineticConstants();
		ptPP->FormationRates();

		// Mass conservation equations
		dn_over_dx = ptPP->R()*(VL*1000.);
	}

	void PrintStep(const double t, const Eigen::VectorXd& n, const Eigen::VectorXd& dn_over_dx)
	{
		step_++;

		const double T = tg_->T(t);

		//if (step_ % 10 == 1)
		{
			//	ptPP->UpdateGasDistribution();


			const double P = tg_->P(t);
			const double W = tg_->InitialMass();

			const int N = ptPP->MaxNumberOfUnits();		// number of monomeric units tracked
			const int LC = ptPP->MinNumberOfUnits();	// minimum number of units in liquid phase
			const double rho = ptPP->LiquidDensity(T);	// liquid-phase density [kg/m3]
			const double mG = ptPP->SumGasMW(n);		// gas-phase mass [g]
			const double mL = ptPP->SumLiquidMW(n);		// liquid-phase mass [g]
			const double nG = ptPP->SumGas(n);			// gas-phase number of moles [mol]
			const double nL = ptPP->SumLiquid(n);		// liquid-phase number of moles [mol]
			const double VL = mL / rho;					// liquid-phase volume [l]
			const double resL = mL / W;			     	// liquid-phase residual [-]

			// Concentrations
			c_ = n / VL;								// concentrations [mol/l]


			// Print on the screen
			{
				if (step_ % 500 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(7) << "#";
					std::cout << std::setw(16) << "time[s]";
					std::cout << std::setw(12) << "T[K]";
					std::cout << std::setw(12) << "LC[-]";
					std::cout << std::setw(12) << "res_liq";
					std::cout << std::setw(12) << "mass_liq[g]";
					std::cout << std::setw(12) << "mass_gas[g]";
					std::cout << std::setw(12) << "mass_tot[g]";
					std::cout << std::endl;
				}

				// On the screen
				std::cout << std::setw(7) << step_;
				std::cout << std::setw(16) << t;
				std::cout << std::setw(12) << T;
				std::cout << std::setw(12) << LC;
				std::cout << std::setw(12) << resL;
				std::cout << std::setw(12) << mL;
				std::cout << std::setw(12) << mG;
				std::cout << std::setw(12) << mL + mG;
				std::cout << std::endl;
			}

			double P_moles = 0.;
			double O_moles = 0.;
			double D_moles = 0.;
			ptPP->Sum(n, P_moles, O_moles, D_moles);
			double T_moles = P_moles + O_moles + D_moles;

			double P_mass = 0.;
			double O_mass = 0.;
			double D_mass = 0.;
			ptPP->SumMW(n, P_mass, O_mass, D_mass);
			double T_mass = P_mass + O_mass + D_mass;

			const double mol_par = P_moles / T_moles;
			const double mol_ole = O_moles / T_moles;
			const double mol_dio = D_moles / T_moles;

			const double mas_par = P_mass / T_mass;
			const double mas_ole = O_mass / T_mass;
			const double mas_dio = D_mass / T_mass;

			fDistMass_ << std::setw(7)  << step_;
			fDistMass_ << std::setw(16) << t;
			fDistMass_ << std::setw(16) << T;
			fDistMass_ << std::setw(16) << LC;
			fDistMass_ << std::setw(16) << 0.;
			fDistMass_ << std::setw(16) << resL;

			fDistMass_ << std::setw(16) << mL;
			fDistMass_ << std::setw(16) << mG;
			fDistMass_ << std::setw(16) << mL + mG;

			fDistMass_ << std::setw(16) << mas_par;
			fDistMass_ << std::setw(16) << mas_ole;
			fDistMass_ << std::setw(16) << mas_dio;
			fDistMass_ << std::endl;

			fDistMoles_ << std::setw(7)  << step_;
			fDistMoles_ << std::setw(16) << t;
			fDistMoles_ << std::setw(16) << T;
			fDistMoles_ << std::setw(16) << LC;
			fDistMoles_ << std::setw(16) << 0.;
			fDistMoles_ << std::setw(16) << resL;

			fDistMoles_ << std::setw(16) << nL;
			fDistMoles_ << std::setw(16) << nG;
			fDistMoles_ << std::setw(16) << nL + nG;

			fDistMoles_ << std::setw(16) << mol_par;
			fDistMoles_ << std::setw(16) << mol_ole;
			fDistMoles_ << std::setw(16) << mol_dio;
			fDistMoles_ << std::endl;
		}

		const int LCnew = tg_->SearchForLC(T);
		const int LCcurrent = ptPP->MinNumberOfUnits();
		if (LCnew != LCcurrent)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << " New boiling temperature (in K)                      " << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << "  * " << tg_->BoilingTemperature(LCcurrent - 1) << " -> " << tg_->BoilingTemperature(LCnew - 1) << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << std::endl;

			ptPP->SetMinNumberOfUnits(LCnew);
		}
	}

private:

	int step_;

	Eigen::VectorXd c_;

	ThermogravimetricAnalysis*						tg_;

	std::ofstream fDistMass_;
	std::ofstream fDistMoles_;
};

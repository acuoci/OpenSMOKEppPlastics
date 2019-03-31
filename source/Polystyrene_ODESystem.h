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

#include "PolystyreneKinetics.h"

extern opensmokepp::plastics::PolystyreneKinetics* ptPS;


class PS_ODESystemObject
{
protected:

	unsigned int number_of_equations_;

public:

	PS_ODESystemObject()
	{
		number_of_equations_ = 5 * ptPS->MaxNumberOfUnits() + 1;
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

		fDistMass_ << std::setw(7) << "#(1)";
		fDistMass_ << std::setw(16) << "time[s](2)";
		fDistMass_ << std::setw(16) << "T[K](3)";
		fDistMass_ << std::setw(16) << "LC[-](4)";
		fDistMass_ << std::setw(16) << "alpha(5)";
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
		fDistMoles_ << std::setw(16) << "alpha[-](5)";
		fDistMoles_ << std::setw(16) << "res_liq(6)";
		fDistMoles_ << std::setw(16) << "mol_liq[mol](7)";
		fDistMoles_ << std::setw(16) << "mol_gas[mol](8)";
		fDistMoles_ << std::setw(16) << "mol_tot[mol](9)";
		fDistMoles_ << std::setw(16) << "P_x(10)";
		fDistMoles_ << std::setw(16) << "O_x(11)";
		fDistMoles_ << std::setw(16) << "D_x(12)";
		fDistMoles_ << std::endl;
	}

	void GetFunctions(const double t, const Eigen::VectorXd& n, Eigen::VectorXd& dn_over_dx)
	{
		const double T = tg_->T(t);
		const double P = tg_->P(t);

		const double rho = ptPS->LiquidDensity(T);	// density [kg/m3]
		const double mL = ptPS->SumLiquidMW(n);		// liquid mass [g]
		const double VL = mL / rho;					// liquid volume [l]

		// Concentrations
		c_ = n / VL;								// concentrations [mol/l]

		// Calculation of formation rates
		ptPS->SetStatus(T, P, c_);
		ptPS->UpdateCCBonds();
		ptPS->KineticConstants();
		ptPS->FormationRates();

		// Mass conservation equations
		dn_over_dx = ptPS->R()*VL;
	}

	void PrintStep(const double t, Eigen::VectorXd& n, const Eigen::VectorXd& dn_over_dx)
	{
		step_++;

		const double T = tg_->T(t);
		const double P = tg_->P(t);

		//if (step_ % 10 == 1)
		{
			ptPS->UpdateGasDistribution();

			const double W = tg_->InitialMass();

			const int N = ptPS->MaxNumberOfUnits();
			const double alpha = ptPS->SplittingCoefficient(T, P);
			const int LC = ptPS->MinNumberOfUnits();				// minimum number of units in liquid phase
			const double rho = ptPS->LiquidDensity(T);				// density [kg/m3]
			const double mG = ptPS->SumGasMW(n);					// gas mass [g]
			const double mL = ptPS->SumLiquidMW(n);					// liquid mass [g]
			const double nG = ptPS->SumGas(n);						// gas moles [mol]
			const double nL = ptPS->SumLiquid(n);					// liquid moles [mol]
			const double VL = mL / rho;								// liquid volume [l]
			const double resL = mL / W;			     				// liquid-phase residual [-]

			// Concentrations
			c_ = n / VL;											// concentrations [mol/l]

			// Print on the screen
			{
				if (step_ % 500 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(7) << "#";
					std::cout << std::setw(16) << "time[s]";
					std::cout << std::setw(12) << "T[K]";
					std::cout << std::setw(12) << "LC[-]";
					std::cout << std::setw(12) << "alpha";
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
				std::cout << std::setw(12) << alpha;
				std::cout << std::setw(12) << resL;
				std::cout << std::setw(12) << mL;
				std::cout << std::setw(12) << mG;
				std::cout << std::setw(12) << mL + mG;
				std::cout << std::endl;
			}

			

			const double mas_par = ptPS->wpar() / ptPS->wtot();
			const double mas_ole = ptPS->wole() / ptPS->wtot();
			const double mas_dio = ptPS->wdio() / ptPS->wtot();
			const double mas_etilben = c_(0) * 106. / ptPS->wtot();
			const double mas_toluene = c_(0 + N) * 92. / ptPS->wtot();
			const double mas_13difprop = c_(1 + N) * 196. / ptPS->wtot();
			const double mas_ametstir = c_(0 + 3 * N) * 118. / ptPS->wtot();
			const double mas_stirene = c_(0 + 2 * N) * 104. / ptPS->wtot();
			const double mas_dimero = c_(1 + 2 * N) * 208. / ptPS->wtot();
			const double mas_trimero = c_(2 + 2 * N) * 312. / ptPS->wtot();

			fDistMass_ << std::setw(7)  << step_;
			fDistMass_ << std::setw(16) << t;
			fDistMass_ << std::setw(16) << T;
			fDistMass_ << std::setw(16) << LC;
			fDistMass_ << std::setw(16) << alpha;
			fDistMass_ << std::setw(16) << resL;

			fDistMass_ << std::setw(16) << mL;
			fDistMass_ << std::setw(16) << mG;
			fDistMass_ << std::setw(16) << mL+mG;

			fDistMass_ << std::setw(16) << mas_par;
			fDistMass_ << std::setw(16) << mas_ole;
			fDistMass_ << std::setw(16) << mas_dio;

		//	fDistMass_ << std::setw(16) << mas_etilben;
		//	fDistMass_ << std::setw(16) << mas_toluene;
		//	fDistMass_ << std::setw(16) << mas_13difprop;
		//	fDistMass_ << std::setw(16) << mas_ametstir;
		//	fDistMass_ << std::setw(16) << mas_stirene;
		//	fDistMass_ << std::setw(16) << mas_dimero;
		//	fDistMass_ << std::setw(16) << mas_trimero;

			fDistMass_ << std::endl;

			const double mol_par = ptPS->ypar() / ptPS->ytot();
			const double mol_ole = ptPS->yole() / ptPS->ytot();
			const double mol_dio = ptPS->ydio() / ptPS->ytot();
			const double mol_etilben = c_(0) / ptPS->ytot();
			const double mol_toluene = c_(0 + N) / ptPS->ytot();
			const double mol_13difprop = c_(1 + N) / ptPS->ytot();
			const double mol_ametstir = c_(0 + 3 * N) / ptPS->ytot();
			const double mol_stirene = c_(0 + 2 * N) / ptPS->ytot();
			const double mol_dimero = c_(1 + 2 * N) / ptPS->ytot();
			const double mol_trimero = c_(2 + 2 * N) / ptPS->ytot();

			fDistMoles_ << std::setw(7)  << step_;
			fDistMoles_ << std::setw(16) << t;
			fDistMoles_ << std::setw(16) << T;
			fDistMoles_ << std::setw(16) << LC;
			fDistMoles_ << std::setw(16) << alpha;
			fDistMoles_ << std::setw(16) << resL;

			fDistMoles_ << std::setw(16) << nL;
			fDistMoles_ << std::setw(16) << nG;
			fDistMoles_ << std::setw(16) << nL+nG;

			fDistMoles_ << std::setw(16) << mol_par;
			fDistMoles_ << std::setw(16) << mol_ole;
			fDistMoles_ << std::setw(16) << mol_dio;

		//	fDistMoles_ << std::setw(16) << mol_etilben;
		//	fDistMoles_ << std::setw(16) << mol_toluene;
		//	fDistMoles_ << std::setw(16) << mol_13difprop;
		//	fDistMoles_ << std::setw(16) << mol_ametstir;
		//	fDistMoles_ << std::setw(16) << mol_stirene;
		//	fDistMoles_ << std::setw(16) << mol_dimero;
		//	fDistMoles_ << std::setw(16) << mol_trimero;
			fDistMoles_ << std::endl;
		}

		const int LCnew = tg_->SearchForLC(T)-1;
		const int LCcurrent = ptPS->MinNumberOfUnits();
		if (LCnew != LCcurrent)
		{
			std::cout << LCnew << " " << LCcurrent << std::endl;

			std::cout << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << " New boiling temperature (in K)                      " << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << "  * " << tg_->BoilingTemperature(LCcurrent - 1) << " -> " << tg_->BoilingTemperature(LCnew - 1) << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << std::endl;

			// Complete the definition of initial distribution
			ptPS->UpdateSharedSpeciesDistribution(LCnew, n);
			ptPS->UpdateSharedSpecies(T, P);
		}
	}

private:

	int step_;

	Eigen::VectorXd c_;

	ThermogravimetricAnalysis*					tg_;

	std::ofstream fDistMass_;
	std::ofstream fDistMoles_;
};



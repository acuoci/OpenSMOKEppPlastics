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

#include "PolyethyleneKinetics.h"

extern opensmokepp::plastics::PolyethyleneKinetics* ptPE;

class PE_ODESystemObject
{
protected:

	unsigned int number_of_equations_;
	unsigned int max_index_distribution_;

public:

	PE_ODESystemObject()
	{
		number_of_equations_ = 3 * ptPE->MaxNumberOfUnits() + 34;
		step_ = 0;
		max_index_distribution_ = 100;
	}

	void SetThermogravimetricAnalysis(ThermogravimetricAnalysis* tg)
	{
		tg_ = tg;
	}

	void PrepareOutputFiles(const boost::filesystem::path& folder_name)
	{
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
			fDistMass_ << std::setw(16) << "dummy(5)";
			fDistMass_ << std::setw(16) << "res_liq(6)";
			fDistMass_ << std::setw(16) << "mass_liq[g](7)";
			fDistMass_ << std::setw(16) << "mass_gas[g](8)";
			fDistMass_ << std::setw(16) << "mass_tot[g](9)";
			fDistMass_ << std::setw(16) << "P_w(10)";
			fDistMass_ << std::setw(16) << "O_w(11)";
			fDistMass_ << std::setw(16) << "D_w(12)";
			fDistMass_ << std::endl;
		}

		{
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

		{
			boost::filesystem::path file_name_moles = folder_name / "GasDistributionMass.P.out";

			fDistMassP_.open(file_name_moles.c_str(), std::ios::out);
			fDistMassP_.setf(std::ios::scientific);
			fDistMassP_ << std::setprecision(6);
			fDistMassP_ << std::left;

			fDistMassP_ << std::setw(7) << "#(1)";
			fDistMassP_ << std::setw(16) << "time[s](2)";
			fDistMassP_ << std::setw(16) << "T[K](3)";
			fDistMassP_ << std::setw(16) << "LC[-](4)";
			fDistMassP_ << std::setw(16) << "dummy(5)";
			fDistMassP_ << std::setw(16) << "res_liq(6)";
			fDistMassP_ << std::setw(16) << "mass_liq[g](7)";
			fDistMassP_ << std::setw(16) << "mass_gas[g](8)";
			fDistMassP_ << std::setw(16) << "mass_tot[g](9)";

			unsigned int count = 10;
			for (unsigned int i = 0; i < max_index_distribution_; i++)
			{
				std::stringstream index; index << i + 1;
				std::stringstream counter; counter << count;
				std::string label = "m_P" + index.str() + "[g](" + counter.str() +")";
				fDistMassP_ << std::setw(16) << label;
				count++;
			}

			fDistMassP_ << std::endl;
		}

		{
			boost::filesystem::path file_name_moles = folder_name / "GasDistributionMass.O.out";

			fDistMassO_.open(file_name_moles.c_str(), std::ios::out);
			fDistMassO_.setf(std::ios::scientific);
			fDistMassO_ << std::setprecision(6);
			fDistMassO_ << std::left;

			fDistMassO_ << std::setw(7) << "#(1)";
			fDistMassO_ << std::setw(16) << "time[s](2)";
			fDistMassO_ << std::setw(16) << "T[K](3)";
			fDistMassO_ << std::setw(16) << "LC[-](4)";
			fDistMassO_ << std::setw(16) << "dummy(5)";
			fDistMassO_ << std::setw(16) << "res_liq(6)";
			fDistMassO_ << std::setw(16) << "mass_liq[g](7)";
			fDistMassO_ << std::setw(16) << "mass_gas[g](8)";
			fDistMassO_ << std::setw(16) << "mass_tot[g](9)";

			unsigned int count = 10;
			for (unsigned int i = 0; i < max_index_distribution_; i++)
			{
				std::stringstream index; index << i + 1;
				std::stringstream counter; counter << count;
				std::string label = "m_O" + index.str() + "[g](" + counter.str() + ")";
				fDistMassO_ << std::setw(16) << label;
				count++;
			}

			fDistMassO_ << std::endl;
		}

		{
			boost::filesystem::path file_name_moles = folder_name / "GasDistributionMass.D.out";

			fDistMassD_.open(file_name_moles.c_str(), std::ios::out);
			fDistMassD_.setf(std::ios::scientific);
			fDistMassD_ << std::setprecision(6);
			fDistMassD_ << std::left;

			fDistMassD_ << std::setw(7) << "#(1)";
			fDistMassD_ << std::setw(16) << "time[s](2)";
			fDistMassD_ << std::setw(16) << "T[K](3)";
			fDistMassD_ << std::setw(16) << "LC[-](4)";
			fDistMassD_ << std::setw(16) << "dummy(5)";
			fDistMassD_ << std::setw(16) << "res_liq(6)";
			fDistMassD_ << std::setw(16) << "mass_liq[g](7)";
			fDistMassD_ << std::setw(16) << "mass_gas[g](8)";
			fDistMassD_ << std::setw(16) << "mass_tot[g](9)";

			unsigned int count = 10;
			for (unsigned int i = 0; i < max_index_distribution_; i++)
			{
				std::stringstream index; index << i + 1;
				std::stringstream counter; counter << count;
				std::string label = "m_D" + index.str() + "[g](" + counter.str() + ")";
				fDistMassD_ << std::setw(16) << label;
				count++;
			}

			fDistMassD_ << std::endl;
		}
	}

	void CloseOutputFiles()
	{
		fDistMass_.close();
		fDistMoles_.close();
		fDistMassP_.close();
		fDistMassO_.close();
		fDistMassD_.close();

	}

	void GetFunctions(const double t, const Eigen::VectorXd& n, Eigen::VectorXd& dn_over_dx)
	{
		const double T = tg_->T(t);
		const double P = tg_->P(t);
		const double W = tg_->InitialMass();

		Eigen::VectorXd n_(n.size());
		for (unsigned int i = 0; i < n.size(); i++)
			n_(i) = std::max(0., n(i));

		const double nG = ptPE->SumGas(n_);							// mol
		const double nL = ptPE->SumLiquid(n_);						// mol
		const double mG = ptPE->SumGasMW(n_);						// g
		const double mL = ptPE->SumLiquidMW(n_);						// g
		const double VL = (mL / 1000.) / ptPE->LiquidDensity(T);	// m3

		c_ = n_ / (VL*1000.);		// concentrations (in kmol/m3 or mol/l)

		ptPE->SetStatus(T, P, c_);
		ptPE->UpdateInitialAccelerationCoefficient(VL*1000., W);
		ptPE->KineticConstants();
		ptPE->FormationRates();

		// Mass conservation equations
		dn_over_dx = ptPE->R()*(VL*1000.);
	}

	void PrintStep(const double t, const Eigen::VectorXd& n, const Eigen::VectorXd& dn_over_dx)
	{
		step_++;

		const double T = tg_->T(t);

		//if (step_ % 10 == 1)
		{
			const double P = tg_->P(t);
			const double W = tg_->InitialMass();

			const int N = ptPE->MaxNumberOfUnits();		// number of monomeric units tracked
			const int LC = ptPE->MinNumberOfUnits();	// minimum number of units in liquid phase
			const double rho = ptPE->LiquidDensity(T);	// liquid-phase density [kg/m3]
			const double mG = ptPE->SumGasMW(n);		// gas-phase mass [g]
			const double mL = ptPE->SumLiquidMW(n);		// liquid-phase mass [g]
			const double nG = ptPE->SumGas(n);			// gas-phase number of moles [mol]
			const double nL = ptPE->SumLiquid(n);		// liquid-phase number of moles [mol]
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
			ptPE->Sum(n, P_moles, O_moles, D_moles);
			double T_moles = P_moles + O_moles + D_moles;

			double P_mass = 0.;
			double O_mass = 0.;
			double D_mass = 0.;
			ptPE->SumMW(n, P_mass, O_mass, D_mass);
			double T_mass = P_mass + O_mass + D_mass;

			const double mol_par = P_moles / T_moles;
			const double mol_ole = O_moles / T_moles;
			const double mol_dio = D_moles / T_moles;

			const double mas_par = P_mass / T_mass;
			const double mas_ole = O_mass / T_mass;
			const double mas_dio = D_mass / T_mass;

			// Mass distribution of paraffins, olefins, and diolefins
			Eigen::VectorXd Ps;
			Eigen::VectorXd Os;
			Eigen::VectorXd Ds;
			ptPE->DistributionMW(n, Ps, Os, Ds);

			{
				fDistMass_ << std::setw(7) << step_;
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
			}

			{
				fDistMoles_ << std::setw(7) << step_;
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

			{
				fDistMassP_ << std::setw(7) << step_;
				fDistMassP_ << std::setw(16) << t;
				fDistMassP_ << std::setw(16) << T;
				fDistMassP_ << std::setw(16) << LC;
				fDistMassP_ << std::setw(16) << 0.;
				fDistMassP_ << std::setw(16) << resL;

				fDistMassP_ << std::setw(16) << mL;
				fDistMassP_ << std::setw(16) << mG;
				fDistMassP_ << std::setw(16) << mL + mG;

				const int LC = ptPE->MinNumberOfUnits();
	
				for (unsigned int i = 0; i < LC-1; i++)
					fDistMassP_ << std::setw(16) << Ps(i);

				for (unsigned int i = LC-1; i < max_index_distribution_; i++)
					fDistMassP_ << std::setw(16) << 0.;

				fDistMassP_ << std::endl;
			}

			{
				fDistMassO_ << std::setw(7) << step_;
				fDistMassO_ << std::setw(16) << t;
				fDistMassO_ << std::setw(16) << T;
				fDistMassO_ << std::setw(16) << LC;
				fDistMassO_ << std::setw(16) << 0.;
				fDistMassO_ << std::setw(16) << resL;

				fDistMassO_ << std::setw(16) << mL;
				fDistMassO_ << std::setw(16) << mG;
				fDistMassO_ << std::setw(16) << mL + mG;

				const int LC = ptPE->MinNumberOfUnits();

				for (unsigned int i = 0; i < LC - 1; i++)
					fDistMassO_ << std::setw(16) << Os(i);

				for (unsigned int i = LC - 1; i < max_index_distribution_; i++)
					fDistMassO_ << std::setw(16) << 0.;

				fDistMassO_ << std::endl;
			}

			{
				fDistMassD_ << std::setw(7) << step_;
				fDistMassD_ << std::setw(16) << t;
				fDistMassD_ << std::setw(16) << T;
				fDistMassD_ << std::setw(16) << LC;
				fDistMassD_ << std::setw(16) << 0.;
				fDistMassD_ << std::setw(16) << resL;

				fDistMassD_ << std::setw(16) << mL;
				fDistMassD_ << std::setw(16) << mG;
				fDistMassD_ << std::setw(16) << mL + mG;

				const int LC = ptPE->MinNumberOfUnits();

				for (unsigned int i = 0; i < LC - 1; i++)
					fDistMassD_ << std::setw(16) << Ds(i);

				for (unsigned int i = LC - 1; i < max_index_distribution_; i++)
					fDistMassD_ << std::setw(16) << 0.;

				fDistMassD_ << std::endl;
			}
		}

		const int LCnew = tg_->SearchForLC(T);
		const int LCcurrent = ptPE->MinNumberOfUnits();
		if (LCnew != LCcurrent)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << " New boiling temperature (in K)                      " << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << "  * " << tg_->BoilingTemperature(LCcurrent - 1) << " -> " << tg_->BoilingTemperature(LCnew - 1) << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << std::endl;

			ptPE->SetMinNumberOfUnits(LCnew);
		}
	}

private:

	int step_;

	Eigen::VectorXd c_;

	ThermogravimetricAnalysis*						tg_;

	std::ofstream fDistMass_;
	std::ofstream fDistMoles_;
	std::ofstream fDistMassP_;
	std::ofstream fDistMassO_;
	std::ofstream fDistMassD_;
};

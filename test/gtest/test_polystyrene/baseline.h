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

#include "Utilities.h"
#include "PolystyreneKinetics.h"

TEST(PS, baseline_01)
{
	// Initial distribution
	const opensmokepp::plastics::PolymerDistribution 
		distribution = opensmokepp::plastics::SCHULTZ_FLORY;	// distribution type
	const double MWE = 50000.;									// polymer molecular weight [g/mol]
	const double MW_monomer = 104.;								// monomer molecular weight [g/mol]
	const double epsilon = 0.1;									// tail-cut in percentage
	const int stretching_coefficient = 20;						// stretching coefficient

	// Lumping
	const bool is_lumping = true;		
	const int lumping_start = 100;
	const int lumping_step  = 10;

	// 
	const double eff = 100.;

	// Operating conditions
	const double T = 600. + 273.;		// temperature (in K)
	const double P = 1.;				// pressure (in atm)
	

	Eigen::VectorXd y;
	int N = 0;

	// Create the initial distribution
	opensmokepp::plastics::InitialDistribution(	distribution, MW_monomer, MWE, 
												lumping_start, lumping_step,
												epsilon, stretching_coefficient, 
												y, N);

	// Adjust the initial distribution according to lumping
	opensmokepp::plastics::LumpingSetup(is_lumping, lumping_start, lumping_step, MW_monomer, 
										y, N);


	opensmokepp::plastics::PolystyreneKinetics kinetics;
	kinetics.SetVerbose(false);
	kinetics.SetMaxNumberOfUnits(N);
	kinetics.UpdateSharedSpecies(T, P);
	kinetics.UpdateSharedSpeciesDistribution(T, P, y);

	const double rho = kinetics.LiquidDensity(T);
	const double Ctot = rho / MWE * 1.001059912;
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
		const double tf = 1;
		const double dt = 1.e-4;
		const int nsteps = static_cast<int>(tf/dt);
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

			kinetics.SetStatus(T, P, c);
			kinetics.UpdateCCBonds();
			kinetics.KineticConstants();
			kinetics.FormationRates();
			n += dt * kinetics.R()*VL*1000.;			// moles
			for (int i = 0; i < n.size(); i++)
				if (n(i) < 0.) n(i) = 0.;

			if (k % 100 == 1)
				std::cout << k << " " << t << " " << mG << " " << mL << " " << mG + mL << std::endl;

			fOut << k << " " << t << " " << mG << " " << mL << " " << mG + mL << std::endl;
		}

		fOut.close();
	}

	EXPECT_TRUE(true);
}


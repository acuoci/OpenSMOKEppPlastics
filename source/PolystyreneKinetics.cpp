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

#include "PolystyreneKinetics.h"

namespace opensmokepp::plastics
{

	double PolystyreneKinetics::Rgas_ = 1.987;

	PolystyreneKinetics::~PolystyreneKinetics()
	{
	}

	PolystyreneKinetics::PolystyreneKinetics()
	{
		// ------------------------------------------------------------------------------------------------
		// Default values from kinetics
		// ------------------------------------------------------------------------------------------------
		DefaultKineticConstants();

		// ------------------------------------------------------------------------------------------------
		// Default values from initial distribution
		// ------------------------------------------------------------------------------------------------
		MWm_ = 104.;
		MWp_ = 50000.;

		is_lumping_enabled_ = true;
		lumping_step_ = 10;
		lumping_start_ = 100;

		// ------------------------------------------------------------------------------------------------
		// Initialization
		// ------------------------------------------------------------------------------------------------
		N_ = 0;
		kd_ = 0;

		CCsr_ = 1.;
		CCsa_ = 0.;

		ypar_ = 0.;
		yole_ = 0.;
		ydio_ = 0.;
		ytot_ = 0.;

		wpar_ = 0.;
		wole_ = 0.;
		wdio_ = 0.;
		wtot_ = 0.;

		is_verbose_ = false;


		// ------------------------------------------------------------------------------------------------
		// Memory allocation
		// ------------------------------------------------------------------------------------------------
		
		radbbp1_.resize(4); radbbp1_.setZero();
		radbbp3_.resize(4); radbbp3_.setZero();
		radbbo1_.resize(4); radbbo1_.setZero();
		polbbp1_.resize(4); polbbp1_.setZero();
		polbbp3_.resize(4); polbbp3_.setZero();
		polbbo1_.resize(4); polbbo1_.setZero();
	}

	void PolystyreneKinetics::SetMaxNumberOfUnits(const int N)
	{
		N_ = N;
		R_.resize(5 * N + 1);
		R_.setZero();
	}

	void PolystyreneKinetics::SetMinNumberOfUnits(const int kd)
	{
		kd_ = kd;
	}

	void PolystyreneKinetics::SetRightSideBetaScissions(const bool flag)
	{
		is_delta_enabled_ = flag;
	}

	void PolystyreneKinetics::SetBackBiting(const bool flag)
	{
		is_backbiting_enabled_ = flag;
	}

	void PolystyreneKinetics::SetRandomScissionEfficiency(const double efficiency)
	{
		eff_ = 1. / efficiency;
	}

	void PolystyreneKinetics::SetAbstractionReactionWeights(const std::vector<double> weights)
	{
		fe1_ = weights[0];
		fe2_ = weights[1];
		fe3_ = weights[2];
	}

	void PolystyreneKinetics::UpdateSharedSpecies(const double T, const double P)
	{
		const int L = static_cast<int>(std::pow(T / 136.*(1. - std::log(P) / 10.5), 2.) / 8.);

		// Boiling temperature
		const double Tb1 = BoilingTemperature(L,P);
		const double Tb2 = BoilingTemperature(L+1, P);

		// Modulation interval (TODO)
		double deltem = 34.;
		if (P != 1.)
			deltem = -32. / std::log(P);

		int kd = L;

		if (std::fabs(T - Tb1) <= deltem)
			kd = L;
		else if (std::abs(T - Tb2) <= deltem)
			kd = L + 1;

		SetMinNumberOfUnits(kd);
	}

	double PolystyreneKinetics::Beta(const double C)
	{
		return kef_ * C / (kef_*C + ku_ + kb_);
	}

	double PolystyreneKinetics::BetaWithoutBackBiting(const double C)
	{
		return kef_ * C / (kef_*C + ku_);
	}

	double PolystyreneKinetics::BetaWithLumping(const double C)
	{
		return kef_ * C / (kef_*C + ku_ / lumping_step_ + kb_);
	}

	double PolystyreneKinetics::Gamma(const double C)
	{
		return kb_ / (kef_*C + ku_ + kb_);
	}

	double PolystyreneKinetics::GammaWithoutBackBiting(const double C)
	{
		return 0.;
	}

	double PolystyreneKinetics::GammaWithLumping(const double C)
	{
		return kb_ / (kef_*C + ku_ / lumping_step_ + kb_);
	}

	void PolystyreneKinetics::DefaultKineticConstants()
	{
		// ------------------------------------------------------------------------------------------------
		// 1. Initiation reactions
		// ------------------------------------------------------------------------------------------------

		As_ = 5.00e13;		// Random scission: PS -> Rp + Rsb
		Es_ = 67500.;		// Random scission: PS -> Rp + Rsb

		Aa_ = 5.00e12;		// Allylic scission: PS -> Ra + Rsb	
		Ea_ = 62500.;		// Allylic scission: PS -> Ra + Rsb	


	   // ------------------------------------------------------------------------------------------------
	   // 2. Termination reactions
	   // ------------------------------------------------------------------------------------------------

		At_ = 5.00e6;		// Termination:	Rsb + Rp -> PS  and  Rp + Rp -> PS
		Et_ = 14000.;		// Termination:	Rsb + Rp -> PS  and  Rp + Rp -> PS


		// ------------------------------------------------------------------------------------------------
		// 3. Primary benzyl radicals
		// ------------------------------------------------------------------------------------------------

		Aar_ = 5e7;			// H-abstractions: Rp + PS -> PS + Rt
		Ear_ = 16500.;		// H-abstractions: Rp + PS -> PS + Rt

		Aar1_ = 5e7;		// H-abstractions: Rp + PS -> PS + Rt
		Ear1_ = 8500.;		// H-abstractions: Rp + PS -> PS + Rt


		// ------------------------------------------------------------------------------------------------
		// 4. Secondary benzyl radicals
		// ------------------------------------------------------------------------------------------------

		Aar2_ = 5e7;		// H-abstractions: Rsb + PS -> PS + Rt
		Ear2_ = 13500.;		// H-abstractions: Rsb + PS -> PS + Rt

		Au_ = 1.0e13;		// Unzipping: Rsb -> Rsb + S
		Eu_ = 26000.;		// Unzipping: Rsb -> Rsb + S	

		Abb_ = 1.0e9;		// Back-biting (intramolecular abstractions): Rsb -> Rt  and  Rp -> Rt
		Ebb_ = 16000.;		// Back-biting (intramolecular abstractions): Rsb -> Rt  and  Rp -> Rt


		// ------------------------------------------------------------------------------------------------
		// 5. Tertiary benzyl radicals
		// ------------------------------------------------------------------------------------------------
		Abeta_ = 1.0e13;	// Beta-decomposition: Rt -> PS + Rsb	
		Ebeta_ = 27000.;	// Beta-decomposition: Rt -> PS + Rsb

		Aar3_ = 5e7;		// H-abstractions: Rt + PS -> PS + Rt'
		Ear3_ = 16500.;		// H-abstractions: Rt + PS -> PS + Rt'


		// ------------------------------------------------------------------------------------------------
		// Right-side beta scissions
		// ------------------------------------------------------------------------------------------------
		is_delta_enabled_ = false;
		delta_ = 0.;
		
		// ------------------------------------------------------------------------------------------------
		// Back-biting
		// ------------------------------------------------------------------------------------------------
		is_backbiting_enabled_ = true;

		// ------------------------------------------------------------------------------------------------
		// Efficiency for random scissions
		// ------------------------------------------------------------------------------------------------
		eff_ = 1./0.01;

		// ------------------------------------------------------------------------------------------------
		// Weights for H-abstraction reactions
		// ------------------------------------------------------------------------------------------------
		fe1_ = 0.;
		fe2_ = 1.;
		fe3_ = 0.;
	}

	double PolystyreneKinetics::LiquidDensity(const double T) const
	{
		return 1086.5 - 0.619*(T - 273.) + 1.36e-4*std::pow(T - 273., 2.);
	}

	double PolystyreneKinetics::BoilingTemperature(const double L, const double P) const
	{
		return 136.*std::sqrt(8. * L) / (1. - std::log(P) / 10.5);
	}

	double PolystyreneKinetics::SplittingCoefficient(const double T, const double P) const
	{
		const double Tb = BoilingTemperature(kd_, P);

		// Modulation interval (TODO)
		double deltem = 34.;
		if (P != 1.)
			deltem = -32. / std::log(P);

		double alpha = 0.5*(1. + (std::tanh((T - Tb) / deltem)) / std::tanh(1.));

		if (alpha >= 1.)      alpha = 1.;
		else if (alpha <= 0.) alpha = 0.;

		return alpha;
	}

	void PolystyreneKinetics::SetCCBonds(const double CCsr, const double CCsa)
	{
		CCsr_ = CCsr;
		CCsa_ = CCsa;
	}

	void PolystyreneKinetics::SetStatus(const double T, const double P, const double Ctot, const Eigen::VectorXd& yy)
	{
		T_ = T;
		P_ = P;
		Ctot_ = Ctot;
		y = Ctot_ * yy;
	}

	void PolystyreneKinetics::SetStatus(const double T, const double P, const Eigen::VectorXd& c)
	{
		T_ = T;
		P_ = P;
		Ctot_ = c.sum();
		y = c;
	}

	void PolystyreneKinetics::KineticConstants()
	{
		// -----------------------------------------------------------------------------------
		// Kinetic constants
		// -----------------------------------------------------------------------------------
		const double ks = As_ * std::exp(-Es_ / (Rgas_*T_));			// random scission
		const double ksa = Aa_ * std::exp(-Ea_ / (Rgas_*T_));			// allylic scission
		const double kt = At_ * T_ * std::exp(-Et_ / (Rgas_*T_));		// termination
		const double kar1 = Aar1_ * std::exp(-Ear1_ / (Rgas_*T_));		// abstraction 1->3b
		const double kar2 = Aar2_ * std::exp(-Ear2_ / (Rgas_*T_));		// abstraction 2b->3b
		const double kar3 = Aar3_ * std::exp(-Ear3_ / (Rgas_*T_));		// abstraction 3b->3b
		const double ker = Aar_ * std::exp(-Ear_ / (Rgas_*T_));			// re-abstraction 3b
		const double kbeta = Abeta_ * std::exp(-Ebeta_ / (Rgas_*T_));	// beta-scission
		ku_ = Au_ * std::exp(-Eu_ / (Rgas_*T_));						// unzipping
		const double kbb = Abb_ * std::exp(-Ebb_ / (Rgas_*T_));			// back-biting

		// Density of liquid phase
		const double rho = LiquidDensity(T_);

		// Concentration (in kmol/m3 or mol/l)
		const double C = rho / MWm_;

		// Averaged kinetic constant for intermolecular abstractions (kef)
		kef_ = (fe1_*kar1 + fe2_ * kar2 + fe3_ * kar3) / (fe1_ + fe2_ + fe3_);

		// Concentration of C-C bonds which can undergo random scission (CCsr_) or allyl scission (CCsa_)
		const double coefficient = (ks / eff_ * CCsr_ + ksa * CCsa_) / kt;
		if (coefficient < 0.)
			return;

		// Total concentration of pseudo radicals (R.)
		const double Rtot = std::sqrt(coefficient);

		// Equivalent rate constant of apparent propagation reactions involving Rt (kp)
		kp_ = kbeta * kef_ / (kbeta + C * ker) * Rtot;

		// Equivalent kinetic constant considering back-biting as abstraction mechanism (kb)
		kb_ = 0.;
		if (is_backbiting_enabled_ == true)
			kb_ = kbeta * kbb / (kbeta + C * ker);

		// Possible corrections
		delta_ = 0.;
		if (is_delta_enabled_ == true)
		{
			const double kdd1 = 1.0e12 *exp(-31000. / (Rgas_*T_));
			const double kdd2 = 1.0e14 *exp(-23000. / (Rgas_*T_));

			const double kp_d = kdd1 * kbb / (kdd1 + C * ker);
			const double kp_t = kdd2 * kbb / (kdd2 + C * ker);

			delta_ = kp_d / (kp_d + kp_t);
		}
	}

	void PolystyreneKinetics::FormationRates()
	{
		sum1_ = 0.;
		sum2_ = 0.;
		sum3_ = 0.;
		sum4_ = 0.;
		sum5_ = 0.;
		sum6_ = 0.;
		sum7_ = 0.;
		sum8_ = 0.;
		sum9_ = 0.;
		sum10_ = 0.;
		sum11_ = 0.;
		sum12_ = 0.;
		sum13_ = 0.;
		sum14_ = 0.;
		sum15_ = 0.;
		sum16_ = 0.;
		sum17_ = 0.;
		sum18_ = 0.;
		sum19_ = 0.;
		sum20_ = 0.;

		radp1_ = 0.;
		radp3_ = 0.;
		rado1_ = 0.;

		sum_monomer_ = 0.;
		sum_trimer_ = 0.;
		sum_13diphenylpropyl_ = 0.;

		radbbp1_.setZero();
		radbbp3_.setZero();
		radbbo1_.setZero();
		polbbp1_.setZero();
		polbbp3_.setZero();
		polbbo1_.setZero();

		R_.setZero();
		if (is_lumping_enabled_ == true)	LiquidFormationRatesWithLumping();
		else                                LiquidFormationRatesWithoutLumping();
		GaseousFormationRates();
	}

	void PolystyreneKinetics::LiquidFormationRatesWithoutLumping()
	{
		// Density of liquid phase
		const double rho = LiquidDensity(T_);

		// Concentration of tertiary Hs on which abstraction can occur (in mol/l or kmol/m3)
		const double C = rho / MWm_;

		// Fraction of secondary benzyl radicals Rsb abstracting tertiary Hs
		double beta = Beta(C);

		// Fraction of secondary benzyl radicals Rsb following back-biting
		double gamma = Gamma(C);

		// Fraction of secondary benzyl radicals Rsb following unzipping
		double lambda = 1. - beta - gamma;

		// TODO
		int indrad = 4;
		int indpol = 4;

		for (int j1 = N_ - 1; j1 >= kd_ + 1; j1--)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + N_ * 2;
			const int j4 = j1 + N_ * 3;
			const int j5 = j1 + N_ * 4;

			// Sums: j=(n+1,Inf)
			if (j1 <= N_ - 2)
			{
				sum3_ += y(j1 + 1 - 1);		// Sum(j=n+1,Inf)PIj
				sum6_ += y(j1 + 1 - 1);		// Sum(j=n+1,Inf)PIj
				sum7_ += y(j2 + 1 - 1);		// Sum(j=n+1,Inf)PIIIj
				sum9_ += y(j4 + 1 - 1);		// Sum(j=n+1,Inf)OIIj
				sum11_ += y(j1 + 1 - 1);	// Sum(j=n+1,Inf)PIj
				sum13_ += y(j3 + 1 - 1);	// Sum(j=n+1,Inf)OIj
				sum14_ += y(j4 + 1 - 1);	// Sum(j=n+1,Inf)OIIj
			}

			// Sums: j=(n+2,Inf)
			if (j1 <= N_ - 3)
			{
				sum1_ += y(j1 + 2 - 1);		// Sum(j=n+2,Inf)PIj
				sum2_ += y(j4 + 2 - 1);		// Sum(j=n+2,Inf)OIIj
				sum4_ += y(j2 + 2 - 1);		// Sum(j=n+2,Inf)PIIIj
				sum5_ += y(j3 + 2 - 1);		// Sum(j=n+2,Inf)OIj
				sum8_ += y(j3 + 2 - 1);		// Sum(j=n+2,Inf)OIj
				sum10_ += y(j5 + 2 - 1);	// Sum(j=n+2,Inf)DIIj
				sum12_ += y(j4 + 2 - 1);	// Sum(j=n+2,Inf)OIIj
				sum15_ += y(j5 + 2 - 1);	// Sum(j=n+2,Inf)DIIj
			}

			// Back-biting conditions
			if (j1 == N_ - 2)
			{
				indrad = 1;
				indpol = 1;
			}
			if (j1 == N_ - 3)
			{
				indrad = 2;
				indpol = 2;
			}
			if (j1 == N_ - 4)
			{
				indrad = 3;
				indpol = 1;
			}
			if (j1 <= N_ - 5)
			{
				if (indrad == 3)	indrad = 1;
				else				indrad += 1;

				if (indpol == 2)	indpol = 1;
				else                indpol += 1;
			}

			// Secondary radicals
			const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_ + radbbp1_(indrad - 1);
			const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_ + radbbp3_(indrad - 1);
			const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_ + radbbo1_(indrad - 1);

			// Paraffins (I): dPIn/dt = -kp*(n-1.5)*PIn + Beta*RPIn
			R_(j1 - 1) = -kp_ * (j1 - 1.5)*y(j1 - 1) + beta * RPI;

			// Paraffins (III): dPIIIn/dt = -kp*(n-2.0)*PIIn + Beta*RPIIIn
			R_(j2 - 1) = -kp_ * (j1 - 2.0)*y(j2 - 1) + beta * RPIII;

			// Olefins (I): dOIn/dt = -kp*(n-2.5)*OIn + 1/2*kp*Sum6 + kp*Sum7 + 1/2*Sum8 + Beta*ROIn + delta*gamma*RPIIIn+2
			R_(j3 - 1) = -kp_ * (j1 - 2.5)*y(j3 - 1) + 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_ + beta * ROI + polbbp3_(indpol - 1);

			// Olefins (II): dOIIn/dt = -kp*(n-2.0)*OIIn + 1/2*kp*Sum11 + 1/2*kp*Sum12 + delta*gamma*RPIn+2
			R_(j4 - 1) = -kp_ * (j1 - 2.0)*y(j4 - 1) + 0.5*kp_*sum11_ + 0.5*kp_*sum12_ + polbbp1_(indpol - 1);

			// Diolefins (II): dDIIn/dt = -kp*(n-3.0)*DIn + 1/2*kp*Sum13 + 1/2*kp*Sum14 + kp*Sum15 + delta*gamma*ROIn+2
			R_(j5 - 1) = -kp_ * (j1 - 3.0)*y(j5 - 1) + 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_ + polbbo1_(indpol - 1);

			// Unzipping: (1. - beta - gamma) is the unzipping fraction
			sum_monomer_ += (1. - beta - gamma)*(RPI + RPIII + ROI);
			radp1_ = (1. - beta - gamma)*RPI;
			radp3_ = (1. - beta - gamma)*RPIII;
			rado1_ = (1. - beta - gamma)*ROI;

			// Back-biting with left-side scission
			sum_trimer_ += (1. - delta_)*gamma*(RPI + RPIII + ROI);
			radbbp1_(indrad - 1) = (1. - delta_)*gamma*RPI;
			radbbp3_(indrad - 1) = (1. - delta_)*gamma*RPIII;
			radbbo1_(indrad - 1) = (1. - delta_)*gamma*ROI;

			// Back-biting with right-side scission
			sum_13diphenylpropyl_ += delta_ * gamma*(RPI + RPIII + ROI);
			polbbp1_(indpol - 1) = delta_ * gamma*RPI;
			polbbp3_(indpol - 1) = delta_ * gamma*RPIII;
			polbbo1_(indpol - 1) = delta_ * gamma*ROI;

			// Shortest species giving back-biting
			if (j1 == kd_ + 4)
			{
				beta = BetaWithoutBackBiting(C);
				gamma = GammaWithoutBackBiting(C);
				lambda = 1. - beta - gamma;
			}
		}

		// Add shortest liquid species
		sum1_ += y(kd_ + 2 - 1);				// Sum(j=n+2,Inf)PIj
		sum2_ += y(kd_ + 2 + 3 * N_ - 1);		// Sum(j=n+2,Inf)OIIj
		sum3_ += y(kd_ + 1 - 1);				// Sum(j=n+1,Inf)PIj
		sum4_ += y(kd_ + 2 + N_ - 1);			// Sum(j=n+2,Inf)PIIIj
		sum5_ += y(kd_ + 2 + 2 * N_ - 1);		// Sum(j=n+2,Inf)OIj
		sum6_ += y(kd_ + 1 - 1);				// Sum(j=n+1,Inf)PIj
		sum7_ += y(kd_ + 1 + N_ - 1);			// Sum(j=n+1,Inf)PIIIj
		sum8_ += y(kd_ + 2 + 2 * N_ - 1);		// Sum(j=n+2,Inf)OIj
		sum9_ += y(kd_ + 1 + 3 * N_ - 1);		// Sum(j=n+1,Inf)OIIj
		sum10_ += y(kd_ + 2 + 4 * N_ - 1);		// Sum(j=n+2,Inf)DIIj
		sum11_ += y(kd_ + 1 - 1);				// Sum(j=n+1,Inf)PIj
		sum12_ += y(kd_ + 2 + 3 * N_ - 1);		// Sum(j=n+2,Inf)OIIj
		sum13_ += y(kd_ + 1 + 2 * N_ - 1);		// Sum(j=n+1,Inf)DIj
		sum14_ += y(kd_ + 1 + 3 * N_ - 1);		// Sum(j=n+1,Inf)OIIj
		sum15_ += y(kd_ + 2 + 4 * N_ - 1);		// Sum(j=n+2,Inf)DIIj

		if (is_verbose_ == true)
		{
			std::ofstream fOut("TestWithoutLumping.out.class", std::ios::out);
			fOut.setf(std::ios::scientific);
			fOut.right; fOut.precision(9);
			for (int j1 = N_; j1 >= 1; j1--)
			{
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				fOut << std::setw(8) << j1; fOut << std::setw(18) << y(j1 - 1);  fOut << std::setw(18) << R_(j1 - 1) << std::endl;
				fOut << std::setw(8) << j2; fOut << std::setw(18) << y(j2 - 1);  fOut << std::setw(18) << R_(j2 - 1) << std::endl;
				fOut << std::setw(8) << j3; fOut << std::setw(18) << y(j3 - 1);  fOut << std::setw(18) << R_(j3 - 1) << std::endl;
				fOut << std::setw(8) << j4; fOut << std::setw(18) << y(j4 - 1);  fOut << std::setw(18) << R_(j4 - 1) << std::endl;
				fOut << std::setw(8) << j5; fOut << std::setw(18) << y(j5 - 1);  fOut << std::setw(18) << R_(j5 - 1) << std::endl;
			}
			fOut.close();
		}
	}

	void PolystyreneKinetics::LiquidFormationRatesWithLumping()
	{
		// Density of liquid phase
		const double rho = LiquidDensity(T_);

		// Concentration (in kmol/m3)
		const double C = rho / MWm_;

		{
			// Fraction of secondary benzyl radicals Rsb following H-abstraction
			double beta = BetaWithLumping(C);

			// Fraction of secondary benzyl radicals Rsb following back-biting
			double gamma = GammaWithLumping(C);

			// Fraction of secondary benzyl radicals Rsb following unzipping
			double lambda = 1. - beta - gamma;

			// Back-biting conditions
			const int indrad = 1;
			const int indpol = 1;

			for (int j1 = N_ - 1; j1 >= lumping_start_ + 1; j1--)
			{
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;
				const int jj = lumping_start_ + lumping_step_ * (j1 - lumping_start_);

				double c1 = 0.;
				double c2 = 0.;
				double c3 = 0.;
				double c4 = 0.;
				double c5 = 0.;

				if (j1 == N_ - 1)
				{
					c1 = 0.;
					c2 = 0.;
					c3 = (lumping_step_ - 1) / 2. - (lumping_step_ - 1) / lumping_step_;
					c4 = 0.;
					c5 = (lumping_step_ - 1) / 2.;
				}
				else if (j1 == N_ - 2)
				{
					c1 = 0.;
					c2 = lumping_step_ - 1. / lumping_step_ - ((lumping_step_ - 1.) / 2. - (lumping_step_ - 1.) / lumping_step_);
					c3 = (lumping_step_ - 1.) / 2. - (lumping_step_ - 1.) / lumping_step_;
					c4 = lumping_step_ - (lumping_step_ - 1.) / 2.;
					c5 = (lumping_step_ - 1.) / 2.;
				}
				else if (j1 <= N_ - 3)
				{
					c1 = lumping_step_ - (lumping_step_ - 1. / lumping_step_);
					c2 = lumping_step_ - 1. / lumping_step_ - ((lumping_step_ - 1.) / 2. - (lumping_step_ - 1.) / lumping_step_);
					c3 = (lumping_step_ - 1.) / 2. - (lumping_step_ - 1.) / lumping_step_;
					c4 = lumping_step_ - (lumping_step_ - 1.) / 2.;
					c5 = (lumping_step_ - 1.) / 2.;
				}

				sum1_ += c1 * y(j1 + 2 - 1) + c2 * y(j1 + 1 - 1) + c3 * y(j1 - 1);		// sums n + 2
				sum2_ += c1 * y(j4 + 2 - 1) + c2 * y(j4 + 1 - 1) + c3 * y(j4 - 1);
				sum4_ += c1 * y(j2 + 2 - 1) + c2 * y(j2 + 1 - 1) + c3 * y(j2 - 1);
				sum5_ += c1 * y(j3 + 2 - 1) + c2 * y(j3 + 1 - 1) + c3 * y(j3 - 1);
				sum8_ += c1 * y(j3 + 2 - 1) + c2 * y(j3 + 1 - 1) + c3 * y(j3 - 1);
				sum10_ += c1 * y(j5 + 2 - 1) + c2 * y(j5 + 1 - 1) + c3 * y(j5 - 1);
				sum12_ += c1 * y(j4 + 2 - 1) + c2 * y(j4 + 1 - 1) + c3 * y(j4 - 1);
				sum15_ += c1 * y(j5 + 2 - 1) + c2 * y(j5 + 1 - 1) + c3 * y(j5 - 1);

				sum3_ += c4 * y(j1 + 1 - 1) + c5 * y(j1 - 1);	// sums n+1
				sum6_ += c4 * y(j1 + 1 - 1) + c5 * y(j1 - 1);
				sum7_ += c4 * y(j2 + 1 - 1) + c5 * y(j2 - 1);
				sum9_ += c4 * y(j4 + 1 - 1) + c5 * y(j4 - 1);
				sum11_ += c4 * y(j1 + 1 - 1) + c5 * y(j1 - 1);
				sum13_ += c4 * y(j3 + 1 - 1) + c5 * y(j3 - 1);
				sum14_ += c4 * y(j4 + 1 - 1) + c5 * y(j4 - 1);

				sum16_ += y(j1 - 1);			// sums non-lumped terms
				sum17_ += y(j2 - 1);
				sum18_ += y(j3 - 1);
				sum19_ += y(j4 - 1);
				sum20_ += y(j5 - 1);


				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_ + radbbp1_(indrad - 1);
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_ + radbbp3_(indrad - 1);
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_ + radbbo1_(indrad - 1);


				R_(j1 - 1) = -kp_ * (jj - 1.5)*y(j1 - 1) + beta * RPI;
				R_(j2 - 1) = -kp_ * (jj - 2.0)*y(j2 - 1) + beta * RPIII;
				R_(j3 - 1) = -kp_ * (jj - 2.5)*y(j3 - 1) + 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_ + beta * ROI + polbbp3_(indpol - 1);
				R_(j4 - 1) = -kp_ * (jj - 2.0)*y(j4 - 1) + 0.5*kp_*sum11_ + 0.5*kp_*sum12_ + polbbp1_(indpol - 1);
				R_(j5 - 1) = -kp_ * (jj - 3.0)*y(j5 - 1) + 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_ + polbbo1_(indpol - 1);

				// Unzipping
				sum_monomer_ += lumping_step_ * (1. - beta - gamma)*(RPI + RPIII + ROI);
				radp1_ = (1. - beta - gamma)*RPI;
				radp3_ = (1. - beta - gamma)*RPIII;
				rado1_ = (1. - beta - gamma)*ROI;

				// Back-biting with left-scission
				sum_trimer_ += (lumping_step_ / 3.)*(1. - delta_)*gamma*(RPI + RPIII + ROI);
				radbbp1_(indrad - 1) = (1. - delta_)*gamma*RPI;
				radbbp3_(indrad - 1) = (1. - delta_)*gamma*RPIII;
				radbbo1_(indrad - 1) = (1. - delta_)*gamma*ROI;

				// Back-biting with right-scission
				sum_13diphenylpropyl_ += delta_ * gamma*(RPI + RPIII + ROI);
				polbbp1_(indpol - 1) = (delta_)*gamma*RPI;
				polbbp3_(indpol - 1) = (delta_)*gamma*RPIII;
				polbbo1_(indpol - 1) = (delta_)*gamma*ROI;

				//sum_monomer_ += lumping_step_*(1.-beta)*(RPI+RPIII+ROI)
				//radp1_  = (1.-beta)*RPI
				//radp3_  = (1.-beta)*RPIII
				//rado1_  = (1.-beta)*ROI
			}
		}

		{
			// Fraction of secondary benzyl radicals Rsb following H-abstraction
			const double beta = Beta(C);

			// Fraction of secondary benzyl radicals Rsb following back-biting
			const double gamma = Gamma(C);

			// Fraction of secondary benzyl radicals Rsb following unzipping
			const double lambda = 1. - beta - gamma;

			// Back-biting conditions
			const int indrad = 1;
			const int indpol = 1;

			{
				const int j1 = lumping_start_;
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;
				const int jj = lumping_start_ + lumping_step_ * (j1 - lumping_start_);
				const double c1 = 1. + (lumping_step_ - 1.) / 2.;
				const double c2 = 1. / lumping_step_;

				const double RPI = radp1_ + 0.5*kp_*(c1*sum16_ - c2 * y(j1 + 1 - 1)) +
					0.5*kp_*(c1*sum19_ - c2 * y(j4 + 1 - 1)) + radbbp1_(indrad - 1);
				const double RPIII = radp3_ + 0.5*kp_* c1*sum16_ + kp_ * (c1*sum17_ - c2 * y(j2 + 1 - 1)) +
					0.5*kp_*(c1*sum18_ - c2 * y(j3 + 1 - 1)) + radbbp3_(indrad - 1);
				const double ROI = rado1_ + 0.5*kp_*(c1*sum18_ - c2 * y(j3 + 1 - 1)) +
					0.5*kp_* c1*sum19_ + kp_ * (c1*sum20_ - c2 * y(j5 + 1 - 1)) + radbbo1_(indrad - 1);

				R_(j1 - 1) = -kp_ * (jj - 1.5)*y(j1 - 1) + beta * RPI;
				R_(j2 - 1) = -kp_ * (jj - 2.0)*y(j2 - 1) + beta * RPIII;
				R_(j3 - 1) = -kp_ * (jj - 2.5)*y(j3 - 1) + 0.5*kp_*c1*sum16_ + kp_ * c1*sum17_ +
					0.5*kp_*(c1*sum18_ - c2 * y(j3 + 1 - 1)) + beta * ROI + polbbp3_(indpol - 1);
				R_(j4 - 1) = -kp_ * (jj - 2.0)*y(j4 - 1) + 0.5*kp_*c1*sum16_ +
					0.5*kp_*(c1*sum19_ - c2 * y(j4 + 1 - 1)) + polbbp1_(indpol - 1);
				R_(j5 - 1) = -kp_ * (jj - 3.0)*y(j5 - 1) + 0.5*kp_*c1*sum18_ + 0.5*kp_* c1*sum19_ +
					kp_ * (c1*sum20_ - c2 * y(j5 + 1 - 1)) + polbbo1_(indpol - 1);

				// Unzipping
				sum_monomer_ += (1. - beta - gamma)*(RPI + RPIII + ROI);
				radp1_ = (1. - beta - gamma)*RPI;
				radp3_ = (1. - beta - gamma)*RPIII;
				rado1_ = (1. - beta - gamma)*ROI;

				// Back-biting with left-scission
				sum_trimer_ += (1. - delta_)*gamma*(RPI + RPIII + ROI);
				radbbp1_(indrad - 1) = (1. - delta_)*gamma*RPI;
				radbbp3_(indrad - 1) = (1. - delta_)*gamma*RPIII;
				radbbo1_(indrad - 1) = (1. - delta_)*gamma*ROI;

				// Back-biting with right-scission
				sum_13diphenylpropyl_ += delta_ * gamma*(RPI + RPIII + ROI);
				polbbp1_(indpol - 1) = delta_ * gamma*RPI;
				polbbp3_(indpol - 1) = delta_ * gamma*RPIII;
				polbbo1_(indpol - 1) = delta_ * gamma*ROI;

				// sum_monomer_ += (1. - beta)*(RPI + RPIII + ROI);
				// radp1_ = (1. - beta)*RPI;
				// radp3_ = (1. - beta)*RPIII;
				// rado1_ = (1. - beta)*ROI;
			}
		}

		{
			// Fraction of secondary benzyl radicals Rsb following H-abstraction
			double beta = Beta(C);

			// Fraction of secondary benzyl radicals Rsb following back-biting
			double gamma = Gamma(C);

			// Fraction of secondary benzyl radicals Rsb following unzipping
			double lambda = 1. - beta - gamma;

			// Back-biting conditions
			int indrad = 1;
			int indpol = 1;

			for (int j1 = lumping_start_ - 1; j1 >= kd_ + 1; j1--)
			{
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;
				const int jj = j1;

				// Back-biting conditions
				if (indrad == 3)	indrad = 1;
				else				indrad += 1;
				if (indpol == 2)	indpol = 1;
				else				indpol += 1;

				const double RPI = radp1_ + 0.5*kp_* sum16_ + 0.5*kp_* sum19_ + radbbp1_(indrad - 1);
				const double RPIII = radp3_ + 0.5*kp_*(sum16_ + y(j1 + 1 - 1)) + kp_ * sum17_ + 0.5*kp_* sum18_ + radbbp3_(indrad - 1);
				const double ROI = rado1_ + 0.5*kp_* sum18_ + 0.5*kp_*(sum19_ + y(j4 + 1 - 1)) + kp_ * sum20_ + radbbo1_(indrad - 1);

				R_(j1 - 1) = -kp_ * (jj - 1.5)*y(j1 - 1) + beta * RPI;
				R_(j2 - 1) = -kp_ * (jj - 2.0)*y(j2 - 1) + beta * RPIII;
				R_(j3 - 1) = -kp_ * (jj - 2.5)*y(j3 - 1) + 0.5*kp_*(sum16_ + y(j1 + 1 - 1)) + kp_ * (sum17_ + y(j2 + 1 - 1)) +
					0.5*kp_* sum18_ + beta * ROI + polbbp3_(indpol - 1);
				R_(j4 - 1) = -kp_ * (jj - 2.0)*y(j4 - 1) + 0.5*kp_*(sum16_ + y(j1 + 1 - 1)) + 0.5*kp_* sum19_ + polbbp1_(indpol - 1);
				R_(j5 - 1) = -kp_ * (jj - 3.0)*y(j5 - 1) + 0.5*kp_*(sum18_ + y(j3 + 1 - 1)) + 0.5*kp_*(sum19_ + y(j4 + 1 - 1)) + kp_ * sum20_ + polbbo1_(indpol - 1);

				// Unzipping
				sum_monomer_ += (1. - beta - gamma)*(RPI + RPIII + ROI);
				radp1_ = (1. - beta - gamma)*RPI;
				radp3_ = (1. - beta - gamma)*RPIII;
				rado1_ = (1. - beta - gamma)*ROI;

				// Back-biting with left-scission
				sum_trimer_ += (1.)*gamma*(RPI + RPIII + ROI);
				radbbp1_(indrad - 1) = (1.)*gamma*RPI;
				radbbp3_(indrad - 1) = (1.)*gamma*RPIII;
				radbbo1_(indrad - 1) = (1.)*gamma*ROI;

				// Back-biting with right-scission
				sum_13diphenylpropyl_ += delta_ * gamma*(RPI + RPIII + ROI);
				polbbp1_(indpol - 1) = delta_ * gamma*RPI;
				polbbp3_(indpol - 1) = delta_ * gamma*RPIII;
				polbbo1_(indpol - 1) = delta_ * gamma*ROI;

				// sum_monomer_ += (1. - beta)*(RPI + RPIII + ROI);
				// radp1_ = (1. - beta)*RPI;
				// radp3_ = (1. - beta)*RPIII;
				// rado1_ = (1. - beta)*ROI;

				sum16_ += y(j1 + 1 - 1);
				sum17_ += y(j2 + 1 - 1);
				sum18_ += y(j3 + 1 - 1);
				sum19_ += y(j4 + 1 - 1);
				sum20_ += y(j5 + 1 - 1);

				if (j1 == kd_ + 4)	// shortest species to give back-biting
				{
					// Fraction of secondary benzyl radicals Rsb following H-abstraction
					beta = BetaWithoutBackBiting(C);

					// Fraction of secondary benzyl radicals Rsb following back-biting
					gamma = GammaWithoutBackBiting(C);

					// Fraction of secondary benzyl radicals Rsb following unzipping
					lambda = 1. - beta - gamma;
				}
			}
		}

		sum1_ = sum16_;
		sum2_ = sum19_;
		sum3_ = sum16_ + y(kd_ + 1 - 1);
		sum4_ = sum17_;
		sum5_ = sum18_;
		sum6_ = sum16_ + y(kd_ + 1 - 1);
		sum7_ = sum17_ + y(kd_ + 1 + N_ - 1);
		sum8_ = sum18_;
		sum9_ = sum19_ + y(kd_ + 1 + 3 * N_ - 1);
		sum10_ = sum20_;
		sum11_ = sum16_ + y(kd_ + 1 - 1);
		sum12_ = sum19_;
		sum13_ = sum18_ + y(kd_ + 1 + 2 * N_ - 1);
		sum14_ = sum19_ + y(kd_ + 1 + 3 * N_ - 1);
		sum15_ = sum20_;

		if (is_verbose_ == true)
		{
			std::ofstream fOut("TestWithLumping.out.class", std::ios::out);
			fOut.setf(std::ios::scientific);
			fOut.right; fOut.precision(9);
			for (int j1 = N_; j1 >= 1; j1--)
			{
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				fOut << std::setw(8) << j1; fOut << std::setw(18) << y(j1 - 1);  fOut << std::setw(18) << R_(j1 - 1) << std::endl;
				fOut << std::setw(8) << j2; fOut << std::setw(18) << y(j2 - 1);  fOut << std::setw(18) << R_(j2 - 1) << std::endl;
				fOut << std::setw(8) << j3; fOut << std::setw(18) << y(j3 - 1);  fOut << std::setw(18) << R_(j3 - 1) << std::endl;
				fOut << std::setw(8) << j4; fOut << std::setw(18) << y(j4 - 1);  fOut << std::setw(18) << R_(j4 - 1) << std::endl;
				fOut << std::setw(8) << j5; fOut << std::setw(18) << y(j5 - 1);  fOut << std::setw(18) << R_(j5 - 1) << std::endl;
			}
			fOut.close();
		}
	}

	void PolystyreneKinetics::GaseousFormationRates()
	{
		// Density of liquid phase
		const double rho = LiquidDensity(T_);

		// Concentration (in kmol/m3)
		const double C = rho / MWm_;

		// Boiling temperature
		const double alpha = SplittingCoefficient(T_, P_);

		// Case: kd > 3
		if (kd_ > 3)
		{
			// Fraction of secondary benzyl radicals Rsb following H-abstraction
			double beta = BetaWithoutBackBiting(C);

			// Shared species
			{
				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_;

				// Shared species (kd/N): liquid side (N)
				{
					const int j1 = N_;
					const int j2 = j1 + N_;
					const int j3 = j1 + N_ * 2;
					const int j4 = j1 + N_ * 3;
					const int j5 = j1 + N_ * 4;

					R_(j1 - 1) = -kp_ * (kd_ - 1.5)*y(j1 - 1) + (1. - alpha)*beta*RPI;
					R_(j2 - 1) = -kp_ * (kd_ - 2.0)*y(j2 - 1) + (1. - alpha)*beta*RPIII;
					R_(j3 - 1) = -kp_ * (kd_ - 2.5)*y(j3 - 1) + (1. - alpha)*(0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_) + (1. - alpha)*beta*ROI;
					R_(j4 - 1) = -kp_ * (kd_ - 2.0)*y(j4 - 1) + (1. - alpha)*(0.5*kp_*sum11_ + 0.5*kp_*sum12_);
					R_(j5 - 1) = -kp_ * (kd_ - 3.0)*y(j5 - 1) + (1. - alpha)*(0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_);
				}

				// Shared species (kd/N): gaseous side (kd)
				{
					const int j1 = kd_;
					const int j2 = j1 + N_;
					const int j3 = j1 + N_ * 2;
					const int j4 = j1 + N_ * 3;
					const int j5 = j1 + N_ * 4;

					R_(j1 - 1) = alpha * beta*RPI;
					R_(j2 - 1) = alpha * beta*RPIII;
					R_(j3 - 1) = alpha * (0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_) + alpha * beta*ROI;
					R_(j4 - 1) = alpha * (0.5*kp_*sum11_ + 0.5*kp_*sum12_);
					R_(j5 - 1) = alpha * (0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_);
				}

				sum_monomer_ += (1. - beta)*(RPI + RPIII + ROI);
				radp1_ = (1. - beta)*RPI;
				radp3_ = (1. - beta)*RPIII;
				rado1_ = (1. - beta)*ROI;

				sum1_ += y(kd_ + 1 - 1);
				sum2_ += y(kd_ + 1 + 3 * N_ - 1);
				sum3_ += y(N_ - 1);
				sum4_ += y(kd_ + 1 + N_ - 1);
				sum5_ += y(kd_ + 1 + 2 * N_ - 1);
				sum6_ += y(N_ - 1);
				sum7_ += y(N_ + N_ - 1);
				sum8_ += y(kd_ + 1 + 2 * N_ - 1);
				sum9_ += y(N_ + 3 * N_ - 1);
				sum10_ += y(kd_ + 1 + 4 * N_ - 1);
				sum11_ += y(N_ - 1);
				sum12_ += y(kd_ + 1 + 3 * N_ - 1);
				sum13_ += y(N_ + 2 * N_ - 1);
				sum14_ += y(N_ + 3 * N_ - 1);
				sum15_ += y(kd_ + 1 + 4 * N_ - 1);
			}

			// Last species entirely in the gaseous phase (just before the shared species)
			{
				const int j1 = kd_ - 1;
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_;

				R_(j1 - 1) = beta * RPI;
				R_(j2 - 1) = beta * RPIII;
				R_(j3 - 1) = 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_ + beta * ROI;
				R_(j4 - 1) = 0.5*kp_*sum11_ + 0.5*kp_*sum12_;
				R_(j5 - 1) = 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_;

				sum_monomer_ += (1. - beta)*(RPI + RPIII + ROI);
				radp1_ = (1. - beta)*RPI;
				radp3_ = (1. - beta)*RPIII;
				rado1_ = (1. - beta)*ROI;

				sum1_ += y(N_ - 1);
				sum2_ += y(N_ + 3 * N_ - 1);
				sum4_ += y(N_ + N_ - 1);
				sum5_ += y(N_ + 2 * N_ - 1);
				sum8_ += y(N_ + 2 * N_ - 1);
				sum12_ += y(N_ + 3 * N_ - 1);
				sum10_ += y(N_ + 4 * N_ - 1);
				sum15_ += y(N_ + 4 * N_ - 1);
			}

			// Loop in the gaseous phase from 3 to kd-2 units (reverse order)
			for (int j1 = kd_ - 2; j1 >= 3; j1--)
			{
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_;

				R_(j1 - 1) = beta * RPI;
				R_(j2 - 1) = beta * RPIII;
				R_(j3 - 1) = 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_ + beta * ROI;
				R_(j4 - 1) = 0.5*kp_*sum11_ + 0.5*kp_*sum12_;
				R_(j5 - 1) = 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_;

				sum_monomer_ += (1. - beta)*(RPI + RPIII + ROI);
				radp1_ = (1. - beta)*RPI;
				radp3_ = (1. - beta)*RPIII;
				rado1_ = (1. - beta)*ROI;
			}

			// Gaseous species with 2 monomeric units 
			{
				const double beta1 = beta;

				const int j1 = 2;
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_;

				R_(j1 - 1) = beta * RPI;
				R_(j2 - 1) = beta * RPIII;
				R_(j3 - 1) = 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_ + beta1 * ROI;
				R_(j4 - 1) = 0.5*kp_*sum11_ + 0.5*kp_*sum12_;
				R_(j5 - 1) = 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_;

				sum_monomer_ += (1. - beta)*(RPI + RPIII) + (1. - beta1)*ROI;
				radp1_ = (1. - beta)*RPI;
				radp3_ = (1. - beta)*RPIII;
				rado1_ = (1. - beta1)*ROI;
			}

			// Gaseous species with 1 monomeric unit
			{
				const int j1 = 1;
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;

				R_(j1 - 1) = RPI;
				R_(j2 - 1) = RPIII;
				R_(j3 - 1) = 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_;
				R_(j4 - 1) = 0.5*kp_*sum11_ + 0.5*kp_*sum12_;
				R_(j5 - 1) = 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_;

				// Explicit correction on OI1: styrene can be obtained only by unzipping
				R_(j3 - 1) = sum_monomer_ + rado1_;

				// Explict correction on DII1: diolefin DII1 does not exist
				R_(j5 - 1) = 0.;

				// Additional term on OI3 (i.e. trimer): back-biting
				R_(j3 - 1 + 2) += sum_trimer_;

				// Additional term on PIII2: 1,3-diphenylpropyl: back-biting 
				R_(j2 - 1 + 1) += sum_13diphenylpropyl_;
			}
		}

		else if (kd_ == 3)
		{
			// Fraction of secondary benzyl radicals Rsb following H-abstraction
			double beta = BetaWithoutBackBiting(C);

			// Shared species
			{
				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;				//+ radbbp1(indrad-1)
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;	//+ radbbp3(indrad-1)
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_;		//+ radbbo1(indrad-1)

				// Species distributed in liquid phase
				{
					const int j1 = N_;
					const int j2 = j1 + N_;
					const int j3 = j1 + N_ * 2;
					const int j4 = j1 + N_ * 3;
					const int j5 = j1 + N_ * 4;

					R_(j1 - 1) = -kp_ * (kd_ - 1.5)*y(j1 - 1) + (1. - alpha)*beta*RPI;
					R_(j2 - 1) = -kp_ * (kd_ - 2.0)*y(j2 - 1) + (1. - alpha)*beta*RPIII;
					R_(j3 - 1) = -kp_ * (kd_ - 2.5)*y(j3 - 1) + (1. - alpha)*(0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_) + (1. - alpha)*beta*ROI;
					R_(j4 - 1) = -kp_ * (kd_ - 2.0)*y(j4 - 1) + (1. - alpha)*(0.5*kp_*sum11_ + 0.5*kp_*sum12_);
					R_(j5 - 1) = -kp_ * (kd_ - 3.0)*y(j5 - 1) + (1. - alpha)*(0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_);
				}

				// Species distributed in gaseous phase
				{
					const int j1 = kd_;
					const int j2 = j1 + N_;
					const int j3 = j1 + N_ * 2;
					const int j4 = j1 + N_ * 3;
					const int j5 = j1 + N_ * 4;

					R_(j1 - 1) = alpha * beta*RPI;
					R_(j2 - 1) = alpha * beta*RPIII;
					R_(j3 - 1) = alpha * (0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_) + alpha * beta*ROI;
					R_(j4 - 1) = alpha * (0.5*kp_*sum11_ + 0.5*kp_*sum12_);
					R_(j5 - 1) = alpha * (0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_);
				}

				sum_monomer_ += (1. - beta)*(RPI + RPIII + ROI);
				radp1_ = (1. - beta)*RPI;
				radp3_ = (1. - beta)*RPIII;
				rado1_ = (1. - beta)*ROI;

				sum1_ += y(kd_ + 1 - 1);
				sum2_ += y(kd_ + 1 + 3 * N_ - 1);
				sum3_ += y(N_ - 1);
				sum4_ += y(kd_ + 1 + N_ - 1);
				sum5_ += y(kd_ + 1 + 2 * N_ - 1);
				sum6_ += y(N_ - 1);
				sum7_ += y(N_ + N_ - 1);
				sum8_ += y(kd_ + 1 + 2 * N_ - 1);
				sum9_ += y(N_ + 3 * N_ - 1);
				sum10_ += y(kd_ + 1 + 4 * N_ - 1);
				sum11_ += y(N_ - 1);
				sum12_ += y(kd_ + 1 + 3 * N_ - 1);
				sum13_ += y(N_ + 2 * N_ - 1);
				sum14_ += y(N_ + 3 * N_ - 1);
				sum15_ += y(kd_ + 1 + 4 * N_ - 1);
			}

			{
				const double beta1 = beta;

				const int j1 = 2;
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_ + sum_13diphenylpropyl_;
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_;

				R_(j1 - 1) = beta * RPI;
				R_(j2 - 1) = beta * RPIII;
				R_(j3 - 1) = 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_ + beta1 * ROI;
				R_(j4 - 1) = 0.5*kp_*sum11_ + 0.5*kp_*sum12_;
				R_(j5 - 1) = 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_;

				sum_monomer_ += (1. - beta)*(RPI + RPIII) + (1. - beta1)*ROI;
				radp1_ = (1. - beta)*RPI;
				radp3_ = (1. - beta)*RPIII;
				rado1_ = (1. - beta1)*ROI;

				sum1_ += y(N_ - 1);
				sum2_ += y(N_ + 3 * N_ - 1);
				sum4_ += y(N_ + N_ - 1);
				sum5_ += y(N_ + 2 * N_ - 1);
				sum8_ += y(N_ + 2 * N_ - 1);
				sum12_ += y(N_ + 3 * N_ - 1);
			}

			{
				const int j1 = 1;
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;

				R_(j1 - 1) = RPI;
				R_(j2 - 1) = RPIII;
				R_(j3 - 1) = 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_;
				R_(j4 - 1) = 0.5*kp_*sum11_ + 0.5*kp_*sum12_;
				R_(j5 - 1) = 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_;

				// Explicit correction on OI1: styrene can be obtained only by unzipping
				R_(j3 - 1) = sum_monomer_ + rado1_;

				// Explict correction on DII1: diolefin DII1 does not exist
				R_(j5 - 1) = 0.;

				// Additional term on OI3 (i.e. trimer): back-biting
				R_(j3 - 1 + 2) += sum_trimer_;

				// Additional term on PIII2: 1,3-diphenylpropyl: back-biting 
				// R_(j2 - 1 + 1) += sum_13diphenylpropyl_;

				// Manipulation of equations for removing diolefin 3 (to be checked)
				R_(j5 - 1 + (N_ - 1)) += -kp_ * y(j5 - 1 + (N_ - 1));
				R_(j5 - 1 + 1) += kp_ * y(j5 - 1 + (N_ - 1));
				R_(j3 - 1) += kp_ * y(j5 - 1 + (N_ - 1));
			}
		}
		else if (kd_ == 2)
		{
			// Fraction of secondary benzyl radicals Rsb following H-abstraction
			double beta = BetaWithoutBackBiting(C);

			{
				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;
				const double ROI = kp_ * (0.5*sum8_ + 0.5*sum9_ + sum10_) + rado1_;

				{
					const int j1 = N_;
					const int j2 = j1 + N_;
					const int j3 = j1 + N_ * 2;
					const int j4 = j1 + N_ * 3;
					const int j5 = j1 + N_ * 4;

					const double coeff1 = std::max(kd_ - 1.5, 0.);
					const double coeff2 = std::max(kd_ - 2., 0.);
					const double coeff3 = std::max(kd_ - 2.5, 0.);
					const double coeff4 = std::max(kd_ - 2., 0.);
					const double coeff5 = std::max(kd_ - 3., 0.);

					R_(j1 - 1) = -kp_ * coeff1*y(j1 - 1) + (1. - alpha)*beta*RPI;
					R_(j2 - 1) = -kp_ * coeff2*y(j2 - 1) + (1. - alpha)*beta*RPIII;
					R_(j3 - 1) = -kp_ * coeff3*y(j3 - 1) + (1. - alpha)*(0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_) + (1. - alpha)*ROI;	// Those radicals do not unzip
					R_(j4 - 1) = -kp_ * coeff4*y(j4 - 1) + (1. - alpha)*(0.5*kp_*sum11_ + 0.5*kp_*sum12_);
					R_(j5 - 1) = -kp_ * coeff5*y(j5 - 1) + (1. - alpha)*(0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_);
				}

				{
					const int j1 = kd_;
					const int j2 = j1 + N_;
					const int j3 = j1 + N_ * 2;
					const int j4 = j1 + N_ * 3;
					const int j5 = j1 + N_ * 4;

					R_(j1 - 1) = alpha * beta*RPI;
					R_(j2 - 1) = alpha * beta*RPIII;
					R_(j3 - 1) = alpha * (0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_) + alpha * ROI; // Those radicals do not unzip
					R_(j4 - 1) = alpha * (0.5*kp_*sum11_ + 0.5*kp_*sum12_);
					R_(j5 - 1) = alpha * (0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_);

					sum_monomer_ += (1. - beta)*(RPI + RPIII + ROI);
					radp1_ = (1. - beta)*RPI;
					radp3_ = (1. - beta)*RPIII;
					rado1_ = 0.;

					sum1_ += y(kd_ + 1 - 1);
					sum2_ += y(kd_ + 1 + 3 * N_ - 1);
					sum3_ += y(N_ - 1);
					sum4_ += y(kd_ + 1 + N_ - 1);
					sum5_ += y(kd_ + 1 + 2 * N_ - 1);
					sum6_ += y(N_ - 1);
					sum7_ += y(N_ + N_ - 1);
					sum8_ += y(kd_ + 1 + 2 * N_ - 1);
					sum9_ += y(N_ + 3 * N_ - 1);
					sum10_ += y(kd_ + 1 + 4 * N_ - 1);
					sum11_ += y(N_ - 1);
					sum12_ += y(kd_ + 1 + 3 * N_ - 1);
					sum13_ += y(N_ + 2 * N_ - 1);
					sum14_ += y(N_ + 3 * N_ - 1);
					sum15_ += y(kd_ + 1 + 4 * N_ - 1);
				}
			}

			{
				const int j1 = 1;
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				const double RPI = kp_ * (0.5*sum1_ + 0.5*sum2_) + radp1_;
				const double RPIII = kp_ * (0.5*sum3_ + sum4_ + 0.5*sum5_) + radp3_;

				R_(j1 - 1) = RPI;
				R_(j2 - 1) = RPIII;
				R_(j3 - 1) = 0.5*kp_*sum6_ + kp_ * sum7_ + 0.5*kp_*sum8_;
				R_(j4 - 1) = 0.5*kp_*sum11_ + 0.5*kp_*sum12_;
				R_(j5 - 1) = 0.5*kp_*sum13_ + 0.5*kp_*sum14_ + kp_ * sum15_;

				// Those sums do not play any role (never used)
				// To be checked, it seems there is an error
				sum1_ += y(N_ - 1);
				sum2_ += y(N_ + 3 * N_ - 1);
				sum4_ += y(N_ + N_ - 1);
				sum5_ += y(N_ + 2 * N_ - 1);
				sum8_ += y(N_ + 2 * N_ - 1);
				sum12_ += y(N_ + 3 * N_ - 1);


				// Explicit correction on OI1: styrene can be obtained only by unzipping
				R_(j3 - 1) = sum_monomer_;

				// Explict correction on DII1: diolefin DII1 does not exist
				R_(j5 - 1) = 0.;

				// Additional term on OI3 (i.e. trimer): back-biting
				R_(j3 - 1 + 2) += sum_trimer_;

				// Additional term on PIII2: 1,3-diphenylpropyl: back-biting 
				R_(j2 - 1 + 1) += sum_13diphenylpropyl_;

				// Manipulation of equations for removing diolefin 3 (to be checked)
				R_(j5 - 1 + 2) += -kp_ * y(j5 - 1 + 2);
				R_(j5 - 1 + 1) += kp_ * y(j5 - 1 + 2);
				R_(j3 - 1) += kp_ * y(j5 - 1 + 2);
			}
		}

		if (is_verbose_ == true)
		{
			const double V = (0.09999873 / rho) * 1000.; // Volume (in l)

			std::ofstream fOut("Test1.out.ps.class", std::ios::out);
			fOut.setf(std::ios::scientific);
			fOut.right; fOut.precision(9);
			for (int j1 = N_; j1 >= 1; j1--)
			{
				const int j2 = j1 + N_;
				const int j3 = j1 + N_ * 2;
				const int j4 = j1 + N_ * 3;
				const int j5 = j1 + N_ * 4;

				fOut << std::setw(8) << j1; fOut << std::setw(18) << y(j1 - 1)*V;  fOut << std::setw(18) << V * R_(j1 - 1) << std::endl;
				fOut << std::setw(8) << j2; fOut << std::setw(18) << y(j2 - 1)*V;  fOut << std::setw(18) << V * R_(j2 - 1) << std::endl;
				fOut << std::setw(8) << j3; fOut << std::setw(18) << y(j3 - 1)*V;  fOut << std::setw(18) << V * R_(j3 - 1) << std::endl;
				fOut << std::setw(8) << j4; fOut << std::setw(18) << y(j4 - 1)*V;  fOut << std::setw(18) << V * R_(j4 - 1) << std::endl;
				fOut << std::setw(8) << j5; fOut << std::setw(18) << y(j5 - 1)*V;  fOut << std::setw(18) << V * R_(j5 - 1) << std::endl;
			}
			fOut.close();
		}
	}
	
	void PolystyreneKinetics::UpdateCCBonds()
	{
		double sumP = 0.;
		double sumO = 0.;
		double sumD = 0.;
		double sumaO = 0.;
		double sumaD = 0.;
		for (int j1 = kd_ + 1; j1 <= N_ - 1; j1++)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + N_ * 2;
			const int j4 = j1 + N_ * 3;
			const int j5 = j1 + N_ * 4;

			sumP += (y(j1 - 1) + y(j2 - 1))*(2.*j1 - 2.);
			sumO += (y(j3 - 1) + y(j4 - 1))*(2.*j1 - 4.);
			sumD += y(j5 - 1)*(2.*j1 - 6.);
			sumaO += y(j3 - 1) + y(j4 - 1);
			sumaD += y(j5 - 1)*2.;
		}

		CCsr_ = sumP + sumO + sumD;
		CCsa_ = sumaO + sumaD;

		if (CCsr_ < 0.) CCsr_ = 0.;
		if (CCsa_ < 0.) CCsa_ = 0.;
	}

	void PolystyreneKinetics::UpdateGasDistribution()
	{
		ysompar1_ = 0.;
		ysompar3_ = 0.;
		ysomole1_ = 0.;
		ysomole2_ = 0.;
		ysomdio2_ = 0.;
		for (int j = 1; j <= kd_; j++)
		{
			ysompar1_ += y(j - 1);
			ysompar3_ += y(j + N_ - 1);
			ysomole1_ += y(j + 2 * N_ - 1);
			ysomole2_ += y(j + 3 * N_ - 1);
			ysomdio2_ += y(j + 4 * N_ - 1);
		}
		ypar_ = ysompar1_ + ysompar3_;
		yole_ = ysomole1_ + ysomole2_;
		ydio_ = ysomdio2_;
		ytot_ = ypar_ + yole_ + ydio_;

		wpar1_ = 0.;
		wpar3_ = 0.;
		wole1_ = 0.;
		wole2_ = 0.;
		wdio2_ = 0.;
		for (int j = 1; j <= kd_; j++)
		{
			const double MW = MWm_ * j;

			wpar1_ += y(j - 1)*(MW + 2.);
			wpar3_ += y(j + N_ - 1)*(MW - 12.);
			wole1_ += y(j + 2 * N_ - 1)*(MW);
			wole2_ += y(j + 3 * N_ - 1)*(MW + 14.);
			wdio2_ += y(j + 4 * N_ - 1)*(MW + 12.);
		}
		wpar_ = wpar1_ + wpar3_;
		wole_ = wole1_ + wole2_;
		wdio_ = wdio2_;
		wtot_ = wpar_ + wole_ + wdio_;
	}

	void PolystyreneKinetics::UpdateSharedSpeciesDistribution(const double T, const double P, Eigen::VectorXd& v) const
	{
		const double alpha = SplittingCoefficient(T, P);

		v(kd_ - 1 - 1) += v(N_ - 1);
		v(kd_ - 1 + N_ - 1) += v(2 * N_ - 1);
		v(kd_ - 1 + 2 * N_ - 1) += v(3 * N_ - 1);
		v(kd_ - 1 + 3 * N_ - 1) += v(4 * N_ - 1);
		v(kd_ - 1 + 4 * N_ - 1) += v(5 * N_ - 1);

		v(N_ - 1) = (1. - alpha) *  v(kd_ - 1);
		v(2 * N_ - 1) = (1. - alpha) *  v(kd_ + N_ - 1);
		v(3 * N_ - 1) = (1. - alpha) *  v(kd_ + 2 * N_ - 1);
		v(4 * N_ - 1) = (1. - alpha) *  v(kd_ + 3 * N_ - 1);
		v(5 * N_ - 1) = (1. - alpha) *  v(kd_ + 4 * N_ - 1);

		v(kd_ - 1) *= alpha;
		v(kd_ + N_ - 1) *= alpha;
		v(kd_ + 2 * N_ - 1) *= alpha;
		v(kd_ + 3 * N_ - 1) *= alpha;
		v(kd_ + 4 * N_ - 1) *= alpha;
	}

	double PolystyreneKinetics::SumGasMW(const Eigen::VectorXd& v) const
	{
		double gas = 0.;
		for (int j1 = 1; j1 <= kd_ - 1; j1++)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + 2 * N_;
			const int j4 = j1 + 3 * N_;
			const int j5 = j1 + 4 * N_;
			const double MW = MWm_ * j1;
			gas +=	v(j1 - 1) * (MW + 2.) +
					v(j2 - 1) * (MW - 12.) +
					v(j3 - 1) * (MW)+
					v(j4 - 1) * (MW + 14.) +
					v(j5 - 1) * (MW + 12.);
		}

		// Largest chain in the gaseous phase
		const double MW = MWm_ * kd_;
		gas +=	v(kd_ - 1)*(MW + 2.) +
				v(kd_ + N_ - 1)*(MW - 12.) +
				v(kd_ + 2 * N_ - 1)*(MW)+
				v(kd_ + 3 * N_ - 1)*(MW + 14.) +
				v(kd_ + 4 * N_ - 1)*(MW + 12.);

		return gas;
	}

	double PolystyreneKinetics::SumLiquidMW(const Eigen::VectorXd& v) const
	{
		double liquid = 0.;
		for (int j1 = kd_ + 1; j1 <= lumping_start_ - 1; j1++)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + 2 * N_;
			const int j4 = j1 + 3 * N_;
			const int j5 = j1 + 4 * N_;

			const double MW = MWm_ * j1;
			liquid +=	v(j1 - 1) * (MW + 2.) +
						v(j2 - 1) * (MW - 12.) +
						v(j3 - 1) * (MW)+
						v(j4 - 1) * (MW + 14.) +
						v(j5 - 1) * (MW + 12.);
		}

		for (int j1 = lumping_start_; j1 <= N_ - 1; j1++)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + 2 * N_;
			const int j4 = j1 + 3 * N_;
			const int j5 = j1 + 4 * N_;

			int px = j1;
			if (is_lumping_enabled_ == true)
				px = lumping_start_ + lumping_step_ * (j1 - lumping_start_);

			const double MW = MWm_ * px;
			liquid +=	v(j1 - 1) * (MW + 2.) +
						v(j2 - 1) * (MW - 12.) +
						v(j3 - 1) * (MW)+
						v(j4 - 1) * (MW + 14.) +
						v(j5 - 1) * (MW + 12.);
		}

		// Be careful: the last position (N_) is occupied by the smallest liquid chain with kd_ units
		//             thus, the molecular weight is MWm_*kd_
		const double MW = MWm_ * kd_;
		liquid += v(N_ - 1)*(MW + 2.) +
			v(N_ + N_ - 1)*(MW - 12.) +
			v(N_ + 2 * N_ - 1)*(MW)+
			v(N_ + 3 * N_ - 1)*(MW + 14.) +
			v(N_ + 4 * N_ - 1)*(MW + 12.);

		return liquid;
	}

	double PolystyreneKinetics::SumGas(const Eigen::VectorXd& v) const
	{
		double gas = 0.;
		for (int j1 = 1; j1 <= kd_ - 1; j1++)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + 2 * N_;
			const int j4 = j1 + 3 * N_;
			const int j5 = j1 + 4 * N_;
			gas += v(j1 - 1) + v(j2 - 1) + v(j3 - 1) + v(j4 - 1) + v(j5 - 1);
		}

		// Largest chain in the gaseous phase
		gas += v(kd_ - 1) + v(kd_ + N_ - 1) + v(kd_ + 2 * N_ - 1) + v(kd_ + 3 * N_ - 1) + v(kd_ + 4 * N_ - 1);

		return gas;
	}

	double PolystyreneKinetics::SumLiquid(const Eigen::VectorXd& v) const
	{
		double liquid = 0.;
		for (int j1 = kd_ + 1; j1 <= lumping_start_ - 1; j1++)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + 2 * N_;
			const int j4 = j1 + 3 * N_;
			const int j5 = j1 + 4 * N_;
			liquid += v(j1 - 1) + v(j2 - 1) + v(j3 - 1) + v(j4 - 1) + v(j5 - 1);
		}

		for (int j1 = lumping_start_; j1 <= N_ - 1; j1++)
		{
			const int j2 = j1 + N_;
			const int j3 = j1 + 2 * N_;
			const int j4 = j1 + 3 * N_;
			const int j5 = j1 + 4 * N_;
			liquid += v(j1 - 1) + v(j2 - 1) + v(j3 - 1) + v(j4 - 1) + v(j5 - 1);
		}

		// Be careful: the last position (N_) is occupied by the smallest liquid chain with kd_ units
		//             thus, the molecular weight is MWm_*kd_
		liquid += v(N_ - 1) + v(N_ + N_ - 1) + v(N_ + 2 * N_ - 1) + v(N_ + 3 * N_ - 1) + v(N_ + 4 * N_ - 1);

		return liquid;
	}
	
}
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

#include "PolyethyleneKinetics.h"

namespace opensmokepp::plastics
{

	double PolyethyleneKinetics::Rgas_ = 1.987;

	PolyethyleneKinetics::~PolyethyleneKinetics()
	{
	}

	PolyethyleneKinetics::PolyethyleneKinetics()
	{
		// ------------------------------------------------------------------------------------------------
		// Default values from kinetics
		// ------------------------------------------------------------------------------------------------
		DefaultKineticConstants();

		// ------------------------------------------------------------------------------------------------
		// Default values from initial distribution
		// ------------------------------------------------------------------------------------------------
		MWm_ = 14.;
		MWp_ = 5000.;

		// ------------------------------------------------------------------------------------------------
		// Initialization
		// ------------------------------------------------------------------------------------------------
		NPA_ = 0;
		NOL_ = 0;
		NDO_ = 0;
		LC_ = 0;

		FR_ = 0.5;
		Cbut_ = 2.;
		Cprop_ = 2.;

		is_initial_acceleration_factor_ = false;
		initial_acceleration_coefficient_ = 1.;

		is_verbose_ = false;


		// ------------------------------------------------------------------------------------------------
		// Memory allocation
		// ------------------------------------------------------------------------------------------------
		k_.resize(9);		k_.setZero();
		CorrP_.resize(20);	CorrP_.setZero();
		CorrO_.resize(20);	CorrO_.setZero();
		sumY_.resize(22);	sumY_.setZero();
		sumY1_.resize(33);	sumY1_.setZero();
		sumY12_.resize(11);	sumY12_.setZero();
		sumY2_.resize(22);	sumY2_.setZero();
		alpha_.resize(11);	alpha_.setZero();
		M_.resize(17);		M_.setZero();
		LW_.resize(16);		LW_.setZero();
		N_.resize(19);		N_.setZero();

	}

	void PolyethyleneKinetics::SetMaxNumberOfUnits(const int N)
	{
		NPA_ = N;
		NOL_ = N;
		NDO_ = N;

		R_.resize(3*N + 1 + 33);
		R_.setZero();
	}

	void PolyethyleneKinetics::SetMinNumberOfUnits(const int LC)
	{
		LC_ = LC;
	}

	void PolyethyleneKinetics::SetInitialAccelerationFactor(const bool flag)
	{
		is_initial_acceleration_factor_ = flag;
	}

	void PolyethyleneKinetics::DefaultKineticConstants()
	{
		A_.resize(9);
		Beta_.resize(9);
		E_.resize(9);

		// ------------------------------------------------------------------------------------------------
		// 1. Initiation reactions
		// ------------------------------------------------------------------------------------------------
		A_(0) = 7.94328e14;		Beta_(0) = 0.0;		E_(0) = 78000.0;	// initiation: random scission
		A_(1) = 4.00e13;		Beta_(1) = 0.0;		E_(1) = 68000.0;	// initiation: allylic scission


	   // ------------------------------------------------------------------------------------------------
	   // 2. Termination reactions
	   // ------------------------------------------------------------------------------------------------
		A_(2) = 2.5e7;			Beta_(2) = 1.0;		E_(2) = 6000.0;		// termination


		// ------------------------------------------------------------------------------------------------
		// 3. Abstractions
		// ------------------------------------------------------------------------------------------------
		A_(3) = 4.74341649e8;	Beta_(3) = 0.0;		E_(3) = 12000.0;	// inter-molecuar abstraction
		A_(4) = 4.74341649e8;	Beta_(4) = 0.0;		E_(4) = 13500.0;	// re-abstraction


		// ------------------------------------------------------------------------------------------------
		// 4. Beta-decompositions
		// ------------------------------------------------------------------------------------------------
		A_(5) = 1.25892541e14;	Beta_(5) = 0.0;		E_(5) = 30500.0;	// beta-decomposition


		// ------------------------------------------------------------------------------------------------
		// 5. Back-biting (or isomerization)
		// ------------------------------------------------------------------------------------------------
		A_(6) = 1.00e11;		Beta_(6) = 0.0;		E_(6) = 20600.0;	// back-biting abstraction 1-4
		A_(7) = 1.60e10;		Beta_(7) = 0.0;		E_(7) = 14500.0;	// back-biting abstraction 1-5
		A_(8) = 5.00e9;			Beta_(8) = 0.0;		E_(8) = 14500.0;	// back-biting abstraction 1-6
	}

	double PolyethyleneKinetics::LiquidDensity(const double T) const
	{
		return 760.;
	}

	double PolyethyleneKinetics::BoilingTemperature(const double L, const double P) const
	{
		return 1.*std::sqrt((L + 1.) / 5.46e-5) / (1. - std::log(P) / 10.5);
	}

	int PolyethyleneKinetics::MinNumberOfUnits(const double T, const double P)
	{
		return static_cast<int>( 5.46e-5*std::pow(T, 2.)*std::pow(1. - std::log(P) / 10.5, 2.) );
	}

	void PolyethyleneKinetics::SetStatus(const double T, const double P, const double Ctot, const Eigen::VectorXd& x)
	{
		T_ = T;
		P_ = P;
		c_ = Ctot * x;
	}

	void PolyethyleneKinetics::SetStatus(const double T, const double P, const Eigen::VectorXd& c)
	{
		T_ = T;
		P_ = P;
		c_ = c;
	}

	void PolyethyleneKinetics::UpdateInitialAccelerationCoefficient(const double VL, const double WG)
	{
		if (is_initial_acceleration_factor_ == true)
		{
			double SOMMAR3 = 0.;
			for (int i = NPA_; i >= LC_; i--)
				SOMMAR3 += i * (c_(i - 1) + c_(i + NPA_ - 1) + c_(i + NPA_ + NOL_ - 1));

			const double csi = SOMMAR3 * MWm_ / WG * VL;
			const double ART = (0.9 - csi) / 0.05;

			initial_acceleration_coefficient_ = std::sqrt(1. + 19.5*(1. - std::tanh(ART)));
		}
	}

	void PolyethyleneKinetics::KineticConstants()
	{
		// Density of liquid phase (kg/m3)
		const double rho = LiquidDensity(T_);

		// Density of liquid phase (kmol/m3 or mol/l)
		const double Ctot = rho / MWm_;

		// Molar volume of flow unit (in m3/kmol or l/mol)
		const double Vs = LC_ / Ctot;	

		// Kinetic constants (iunits: mol,l,s)
		for (int i = 0; i < 9; i++)
			k_(i) = A_(i)*std::pow(T_, Beta_(i))*std::exp(-E_(i) / Rgas_ / T_);

		// Liquid phase
		double SOMMAR1 = 0.;
		double SOMMAR3 = 0.;
		for (int i = NPA_; i >= LC_; i--)
		{
			SOMMAR1 += (i - 3)*c_(i - 1) + (i - 5)*c_(i + NPA_ - 1) + (i - 7)*c_(i + NPA_ + NOL_ - 1);
			SOMMAR3 += i * (c_(i - 1) + c_(i + NPA_ - 1) + c_(i + NPA_ + NOL_ - 1));
		}

		// Liquid phase
		double SOMMAR2 = 0.;
		for (int i = NPA_ + NOL_; i >= NPA_ + LC_; i--)
			SOMMAR2 += c_(i - 1) + 2.*c_(i + NOL_ - 1);

		
		// Acceleration coefficient during the first phases of degradation
		// It accounts for the possible presence of ramifications or impurities
		// having a catalytic effect (i.e. reduction of activation energy)
		k_(0) *= initial_acceleration_coefficient_;
		
		// Termination reaction kinetic constant
		const double AKT = k_(2);

		// Chain initiation kinetic constant (to form 2 primary radicals)
		const double AKS = (k_(0)*SOMMAR1 + k_(1)*SOMMAR2) / SOMMAR3;

		// Radical concentration (in mol/l or kmol/m3)
		const double RAD = std::sqrt(2.*AKS/AKT * Ctot);

		// Apparent propagation reaction kinetic constant
		const double AKI = k_(5)*k_(3) / (k_(4)*Ctot + k_(5));

		// Propagation kinetic constant: kp=ki*[R]
		kp_ = RAD * AKI;

		// Auxiliary kinetic constants
		kd_ = 0.5*kp_*FR_*Cbut_;
		ka_ = 0.5*kp_*Cbut_;
	}

	void PolyethyleneKinetics::FormationRates()
	{
		// Density of liquid phase
		const double rho = LiquidDensity(T_);

		// Density of liquid phase
		const double Ctot = rho / MWm_;

		// Reset variables
		sum1_ = 0.;
		sum2_ = 0.;
		sum3_ = 0.;
		sum4_ = 0.;
		sum5_ = 0.;
		sum6_ = 0.;
		sum7_ = 0.;
		sumY_.setZero();
		sumY1_.setZero();
		sumY12_.setZero();
		sumY2_.setZero();

		// Calculate the corrections
		UpdateCorrections();

		// Calculates the back-biting reactions probabilities
		{
			Eigen::VectorXd AKES(2);
			for (int i = 1; i <= 2; i++)
				AKES(i - 1) = k_(i + 3 - 1);

			Eigen::VectorXd AKIS(3);
			for (int i = 1; i <= 3; i++)
				AKIS(i - 1) = k_(i + 6 - 1);

			Eigen::VectorXd p(3);
			Eigen::VectorXd ps(3);
			for (int i = 1; i <= 3; i++)
			{
				p(i - 1) = AKIS(i - 1) / (AKIS(i - 1) + AKES(1 - 1)*Ctot);
				ps(i - 1) = AKIS(i - 1) / (AKIS(i - 1) + AKES(2 - 1)*Ctot);
			}

			alpha_(0) = p(0);
			alpha_(1) = p(1);
			alpha_(2) = p(2);
			alpha_(3) = p(0)*ps(0);
			alpha_(4) = p(0)*ps(1) + p(1)*ps(0);
			alpha_(5) = p(0)*ps(2) + p(2)*ps(0) + p(1)*ps(1);
			alpha_(6) = p(1)*ps(2) + p(2)*ps(1) + 2.*p(0)*std::pow(ps(0), 2.);
			alpha_(7) = p(2)*ps(2) + p(0)*ps(0)*ps(1) +
				p(0)*ps(1)*ps(0) + 2.*p(1)*std::pow(ps(0), 2.);
			alpha_(8) = p(0)*ps(0)*ps(2) + p(0)*ps(2)*ps(0) + 2.*p(2) *
				std::pow(ps(0), 2.) + p(1)*ps(1)*ps(0) +
				2.*p(0)*std::pow(ps(1), 2.) + p(1)*ps(0)*ps(1);
			alpha_(9) = 6.*p(0)*std::pow(ps(0), 3.) + p(0)*ps(1)*ps(2) +
				p(0)*ps(2)*ps(1) + p(1)*ps(0)*ps(2) +
				p(1)*ps(2)*ps(0) + p(2)*ps(1)*ps(0) +
				p(2)*ps(0)*ps(1) + 2.*p(1)*std::pow(ps(1), 2.);
			alpha_(10) = 2.*p(0)*std::pow(ps(0), 2.) * ps(1) +
				2.*p(0)*ps(0)*ps(1)*ps(0) +
				2.*p(0)*ps(1)*std::pow(p(0), 2.) +
				6.*p(1)*std::pow(ps(0), 3.) + 2.*p(2)*std::pow(ps(1), 2.) +
				p(1)*ps(2)*ps(1) + p(1)*ps(1)*ps(2) +
				p(2)*ps(2)*ps(0) + p(2)*ps(0)*ps(2) +
				2.*p(0)*std::pow(ps(2), 2.);
		}

		ATOT_ = alpha_.sum();

		M_.setZero();
		LW_.setZero();
		N_.setZero();

		M_(0) = 1;
		N_(0) = 1;
		LW_(0) = 1;

		R_.setZero();
		LiquidFormationRates();
		GaseousFormationRates();
	}

	void PolyethyleneKinetics::LiquidFormationRates()
	{
		int IV = 1;

		for (int j = NPA_; j >= LC_; j--)
		{
			UpdateSums(j);

			const int j1 = j + NPA_;
			const int j2 = j + NPA_ * 2;

			if (j >= LC_)
			{
				if (j >= NPA_ - 1)
				{
					M_(IV - 1) = 1;
					N_(IV - 1) = 1;
					LW_(IV - 1) = 1;
				}
				else
				{
					IV++;

					if (IV > 19)
					{
						// do nothing
					}
					else
					{
						N_(IV - 1) = 1;
					}

					if (IV > 16)
					{
						// do nothing
					}
					else
					{
						LW_(IV - 1) = 1;
					}

					if (IV >= 5)
					{
						if (IV > 19)
						{
							// leave the if region
						}
						else
						{
							M_(IV - 2 - 1) = 1;
						}
					}
					else
					{
						M_(IV - 1) = 1;
					}
				}


				Eigen::VectorXd SC(7);
				SC.setZero();
				for (int ISC = 0; ISC < 11; ISC++)
				{
					SC(0) += alpha_(ISC)*sumY_(ISC);
					SC(1) += (alpha_(ISC)*(M_(ISC + 5)*Cbut_*c_(j1 + ISC + 8) + sumY_(ISC + 11)));
					SC(2) += alpha_(ISC)*sumY1_(ISC);
					SC(3) += (alpha_(ISC)*(N_(ISC + 4)*Cbut_*c_(j1 + ISC + 5) + sumY1_(ISC + 11) + sumY12_(ISC)));
					SC(4) += (alpha_(ISC)*(N_(ISC + 7)*Cbut_*c_(j2 + ISC + 8) + sumY1_(ISC + 22)));
					SC(5) += (alpha_(ISC)*(LW_(ISC + 4)*Cbut_*c_(j2 + ISC + 5) + sumY2_(ISC + 11)));
					SC(6) += alpha_(ISC)*sumY2_(ISC);
				}

				// Paraffins
				R_(j - 1) = kp_ * (-(j - 4)*c_(j - 1) + sum1_ * (1. - ATOT_)
					+ (1. - FR_)*(sum2_ + M_(3 - 1)*Cbut_*c_(j1 + 4 - 1))*(1. - ATOT_)
					+ 0.5*SC(1 - 1) + 0.5*FR_*SC(2 - 1));

				// Olefins
				R_(j1 - 1) = kp_ * (-((j - 6) + FR_ * Cprop_ + (1. - FR_)*Cbut_)*c_(j1 - 1)
					+ sum3_ + FR_ * (N_(2 - 1)*Cprop_*c_(j1 + 3 - 1) + sum4_)
					+ (N_(3 - 1)*Cbut_*c_(j2 + 4 - 1) + sum5_)*(1. - ATOT_)
					+ FR_ * (N_(2 - 1)*c_(j1 + 3 - 1) + sum4_)*(1. - ATOT_)
					+ 0.5*SC(3 - 1) + 0.5*FR_*SC(4 - 1) + 0.5*SC(5 - 1));

				// Diolefins
				R_(j2 - 1) = kp_ * (-((j - 8) + Cbut_ + Cprop_)*c_(j2 - 1) + (1. - FR_)*sum6_
					+ sum7_ + LW_(2 - 1)*Cprop_*c_(j2 + 3 - 1)
					+ 0.5*SC(6 - 1) + 0.5*FR_*SC(7 - 1));
			}

		}
	}

	void PolyethyleneKinetics::GaseousFormationRates()
	{
		const double ABB1 = alpha_(0);
		const double ABB2 = alpha_(1);
		const double ABB3 = alpha_(2);
		const double ABB4 = alpha_(3);
		const double ABB5 = alpha_(4);
		const double ABB6 = alpha_(5);
		const double ABB7 = alpha_(6);
		const double ABB8 = alpha_(7);
		const double ABB9 = alpha_(8);
		const double ABB10 = alpha_(9);
		const double ABB11 = alpha_(10);

		for (int J = LC_ - 1; J >= 2; J--)
		{
			UpdateSums(J);

			const double SY_1 = sumY_(0);
			const double SY_2 = sumY_(1);
			const double SY_3 = sumY_(2);
			const double SY_4 = sumY_(3);
			const double SY_5 = sumY_(4);
			const double SY_6 = sumY_(5);
			const double SY_7 = sumY_(6);
			const double SY_8 = sumY_(7);
			const double SY_9 = sumY_(8);
			const double SY_10 = sumY_(9);
			const double SY_11 = sumY_(10);
			const double SY_12 = sumY_(11);
			const double SY_13 = sumY_(12);
			const double SY_14 = sumY_(13);
			const double SY_15 = sumY_(14);
			const double SY_16 = sumY_(15);
			const double SY_17 = sumY_(16);
			const double SY_18 = sumY_(17);
			const double SY_19 = sumY_(18);
			const double SY_20 = sumY_(19);
			const double SY_21 = sumY_(20);
			const double SY_22 = sumY_(21);

			const double SY1_1 = sumY1_(0);
			const double SY1_2 = sumY1_(1);
			const double SY1_3 = sumY1_(2);
			const double SY1_4 = sumY1_(3);
			const double SY1_5 = sumY1_(4);
			const double SY1_6 = sumY1_(5);
			const double SY1_7 = sumY1_(6);
			const double SY1_8 = sumY1_(7);
			const double SY1_9 = sumY1_(8);
			const double SY1_10 = sumY1_(9);
			const double SY1_11 = sumY1_(10);
			const double SY1_12 = sumY1_(11);
			const double SY1_13 = sumY1_(12);
			const double SY1_14 = sumY1_(13);
			const double SY1_15 = sumY1_(14);
			const double SY1_16 = sumY1_(15);
			const double SY1_17 = sumY1_(16);
			const double SY1_18 = sumY1_(17);
			const double SY1_19 = sumY1_(18);
			const double SY1_20 = sumY1_(19);
			const double SY1_21 = sumY1_(20);
			const double SY1_22 = sumY1_(21);
			const double SY1_23 = sumY1_(22);
			const double SY1_24 = sumY1_(23);
			const double SY1_25 = sumY1_(24);
			const double SY1_26 = sumY1_(25);
			const double SY1_27 = sumY1_(26);
			const double SY1_28 = sumY1_(27);
			const double SY1_29 = sumY1_(28);
			const double SY1_30 = sumY1_(29);
			const double SY1_31 = sumY1_(30);
			const double SY1_32 = sumY1_(31);
			const double SY1_33 = sumY1_(32);

			const double SY12_1 = sumY12_(0);
			const double SY12_2 = sumY12_(1);
			const double SY12_3 = sumY12_(2);
			const double SY12_4 = sumY12_(3);
			const double SY12_5 = sumY12_(4);
			const double SY12_6 = sumY12_(5);
			const double SY12_7 = sumY12_(6);
			const double SY12_8 = sumY12_(7);
			const double SY12_9 = sumY12_(8);
			const double SY12_10 = sumY12_(9);
			const double SY12_11 = sumY12_(10);

			const double SY2_1 = sumY2_(0);
			const double SY2_2 = sumY2_(1);
			const double SY2_3 = sumY2_(2);
			const double SY2_4 = sumY2_(3);
			const double SY2_5 = sumY2_(4);
			const double SY2_6 = sumY2_(5);
			const double SY2_7 = sumY2_(6);
			const double SY2_8 = sumY2_(7);
			const double SY2_9 = sumY2_(8);
			const double SY2_10 = sumY2_(9);
			const double SY2_11 = sumY2_(10);
			const double SY2_12 = sumY2_(11);
			const double SY2_13 = sumY2_(12);
			const double SY2_14 = sumY2_(13);
			const double SY2_15 = sumY2_(14);
			const double SY2_16 = sumY2_(15);
			const double SY2_17 = sumY2_(16);
			const double SY2_18 = sumY2_(17);
			const double SY2_19 = sumY2_(18);
			const double SY2_20 = sumY2_(19);
			const double SY2_21 = sumY2_(20);
			const double SY2_22 = sumY2_(21);

			const int j0 = J - 1;
			const int j1 = J + NPA_ - 1;
			const int j2 = J + NPA_ * 2 - 1;

			if (J >= 16)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ATOT_))
					+ (1. - FR_)*(sum2_)*(1. - (ATOT_))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22)));

				if (J      >= LC_)	R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_)	R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4) * (1. - ATOT_);
				if (J +  9 >= LC_)	R_(j0) += kd_ * ABB1*c_(j1 + 9);
				if (J + 10 >= LC_)	R_(j0) += kd_ * ABB2*c_(j1 + 10);
				if (J + 11 >= LC_)	R_(j0) += kd_ * ABB3*c_(j1 + 11);
				if (J + 12 >= LC_)	R_(j0) += kd_ * ABB4*c_(j1 + 12);
				if (J + 13 >= LC_)	R_(j0) += kd_ * ABB5*c_(j1 + 13);
				if (J + 14 >= LC_)	R_(j0) += kd_ * ABB6*c_(j1 + 14);
				if (J + 15 >= LC_)	R_(j0) += kd_ * ABB7*c_(j1 + 15);
				if (J + 16 >= LC_)	R_(j0) += kd_ * ABB8*c_(j1 + 16);
				if (J + 17 >= LC_)	R_(j0) += kd_ * ABB9*c_(j1 + 17);
				if (J + 18 >= LC_)	R_(j0) += kd_ * ABB10*c_(j1 + 18);
				if (J + 19 >= LC_)	R_(j0) += kd_ * ABB11*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ATOT_))
					+ FR_ * (sum4_)*(1. - (ATOT_))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
						+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10
						+ ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33)));

				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) * Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - ATOT_));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4)*(1. - ATOT_);
				if (J +  6 >= LC_) R_(j1) += kd_ * ABB1*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += kd_ * ABB2*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += kd_ * ABB3*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += kd_ * ABB4*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += kd_ * ABB5*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += kd_ * ABB6*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += kd_ * ABB7*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += kd_ * ABB8*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += kd_ * ABB9*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += kd_ * ABB10*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += kd_ * ABB11*c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += ka_ * ABB1*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += ka_ * ABB2*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += ka_ * ABB3*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += ka_ * ABB4*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += ka_ * ABB5*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += ka_ * ABB6*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += ka_ * ABB7*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += ka_ * ABB8*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += ka_ * ABB9*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += ka_ * ABB10*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += ka_ * ABB11*c_(j2 + 19);

				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));

				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += ka_ * ABB1*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += ka_ * ABB2*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += ka_ * ABB3*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += ka_ * ABB4*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += ka_ * ABB5*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += ka_ * ABB6*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += ka_ * ABB7*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += ka_ * ABB8*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += ka_ * ABB9*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += ka_ * ABB10*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += ka_ * ABB11*c_(j2 + 16);
			}

			else if (J == 15)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ATOT_))
					+ (1. - FR_)*(sum2_)*(1. - (ATOT_))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22)));

				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4) * (1 - ATOT_);
				if (J +  9 >= LC_) R_(j0) += kd_ * ABB1*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += kd_ * ABB2*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += kd_ * ABB3*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += kd_ * ABB4*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += kd_ * ABB5*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += kd_ * ABB6*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += kd_ * ABB7*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += kd_ * ABB8*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += kd_ * ABB9*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += kd_ * ABB10*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += kd_ * ABB11*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9 + ABB10))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9 + ABB10))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
						+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32))
					+ ABB11 * CorrO_(j0));

				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) * Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ATOT_ - ABB11)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) * (1. - (ATOT_ - ABB11));
				if (J +  6 >= LC_) R_(j1) += kd_ * ABB1*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += kd_ * ABB2*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += kd_ * ABB3*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += kd_ * ABB4*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += kd_ * ABB5*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += kd_ * ABB6*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += kd_ * ABB7*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += kd_ * ABB8*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += kd_ * ABB9*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += kd_ * ABB10*c_(j1 + 15);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);

				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));

				if (J >= LC_) R_(j2) += -kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J + 3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J + 6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J + 7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J + 8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J + 9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);
			}

			else if (J == 14)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6
						+ ABB7 + ABB8 + ABB9 + ABB10))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9 + ABB10))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22)));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - (ATOT_ - ABB11));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
						+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB11 * (SY1_33))
					+ ABB10 * CorrO_(j0));

				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) *Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ATOT_ - ABB11 - ABB10)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) * (1. - (ATOT_ - ABB11 - ABB10));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 14);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ *c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);


				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));

				if (J      >= LC_) R_(j2) += - kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);
			}

			else if (J == 13)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6
						+ ABB7 + ABB8 + ABB9))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22)));

				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - (ATOT_ - ABB11 - ABB10));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
						+ ABB8 * SY1_8 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB9 * CorrO_(j0));

				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) * Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ATOT_ - ABB11 - ABB10 - ABB9)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) *(1. - (ATOT_ - ABB11 - ABB10 - ABB9));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ * c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ * c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ *c_(j1 + 13);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ *c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ * c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);


				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));

				if (J      >= LC_) R_(j2) += -kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);
			}

			else if (J == 12)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6
						+ ABB7 + ABB8))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7 + ABB8))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21))
					+ ABB11 * CorrP_(3 + J - 1));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - (ATOT_ - ABB11 - ABB10 - ABB9));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB8 * CorrO_(j0));

				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) *Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) * (1. - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ * c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ *c_(j1 + 12);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ * c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ *c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ *  c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);

				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));

				if (J      >= LC_) R_(j2) += -kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);
			}

			else if (J == 11)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6 + ABB7))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6 + ABB7))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB11 * (SY_22))
					+ ABB10 * CorrP_(3 + J - 1));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5
						+ ABB6))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5
						+ ABB6))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB8 * SY1_8
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB7 * CorrO_(j0));


				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) *Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8 - ABB7)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) *(1. - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8 - ABB7));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ * c_(j1 + 11);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ *c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ *c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ * c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ * c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);


				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));


				if (J      >= LC_) R_(j2) += -kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);
			}

			else if (J == 10)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5 + ABB6))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB9 * CorrP_(3 + J - 1));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8 - ABB7));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB7 * SY1_7 + ABB8 * SY1_8
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB6 * CorrO_(j0));


				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) *Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8 - ABB7 - ABB6)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) *(1. - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8 - ABB7 - ABB6));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ * c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ *c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ *c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ *c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ *c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);


				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));


				if (J      >= LC_) R_(j2) += -kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);
			}

			else if (J == 9)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
						+ ABB5))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB8 * CorrP_(3 + J - 1));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - (ATOT_ - ABB11 - ABB10 - ABB9 - ABB8 - ABB7 - ABB6));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB4 * SY1_4
						+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB5 * CorrO_(j0));

				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) *Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ABB1 + ABB2 + ABB3 + ABB4)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) *(1. - (ABB1 + ABB2 + ABB3 + ABB4));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ *c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ *c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ * c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ *c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ * c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ * c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);


				R_(j2) = kp_ * ((1. - FR_)*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));


				if (J      >= LC_) R_(j2) += -kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);

			}

			else if (J == 8)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB7 * CorrP_(3 + J - 1));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1. - (ABB1 + ABB2 + ABB3 + ABB4));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);

				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2 + ABB3))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2 + ABB3))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
						+ ABB5 * SY1_5
						+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB4 * CorrO_(j0));


				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) *Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ABB1 + ABB2 + ABB3)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) * (1. - (ABB1 + ABB2 + ABB3));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ *c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ *c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ *c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ * c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ * c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ * c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);




				R_(j2) = kp_ * (1. / 2.*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));


				if (J      >= LC_) R_(j2) += -kp_ * ((J - 8) + Cbut_ + Cprop_)*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);

			}

			else if (J == 7)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2 + ABB3))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2 + ABB3))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB7 * (+SY_18)
						+ ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB6 * CorrP_(3 + J - 1));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - (ABB1 + ABB2 + ABB3));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - (ABB1 + ABB2))
					+ FR_ * (sum4_)*(1. - (ABB1 + ABB2))
					+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB4 * SY1_4
						+ ABB5 * SY1_5
						+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB2 * (SY1_13 + SY12_2)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB3 * CorrO_(j0));


				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) * Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ABB1 + ABB2)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) *(1. - (ABB1 + ABB2));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ *c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ * c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ *c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ * c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ * c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ * c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);

				R_(j2) = kp_ * (1. / 2.*sum6_ + sum7_
					+ 0.5*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 0.5*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));

				if (J      >= LC_) R_(j2) += -kp_ * Cbut_*Cprop_*c_(j2);
				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 16);
			}

			else if (J == 6)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - (ABB1 + ABB2))
					+ (1. - FR_)*(sum2_)*(1. - (ABB1 + ABB2))
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB5 * CorrP_(3 + J - 1));


				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4) * (1. - (ABB1 + ABB2));
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ +
					FR_ * (sum4_)
					+(sum5_)*(1. - ABB1)
					+ FR_ * (sum4_)*(1. - ABB1)
					+ 0.5*(ABB1*SY1_1 + ABB3 * SY1_3 + ABB4 * SY1_4
						+ ABB5 * SY1_5
						+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB1*(SY1_12 + SY12_1)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB2 * CorrO_(j0));


				if (J      >= LC_) R_(j1) += -kp_ * ((J - 6) + FR_ * Cprop_ + (1. - FR_) *Cbut_)*c_(j1);
				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3)*(1. - (ABB1)));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4) *(1. - (ABB1));
				if (J +  6 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ *c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ *c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ * c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ * c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ * c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ * c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);

				R_(j2) = kp_ * (1. / 2.*sum6_ + sum7_
					+ 1. / (1. + Cprop_)*(ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ 1. / (1. + Cprop_)*FR_*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));


				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB1*Cbut_ *c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB2*Cbut_ *c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB3*Cbut_ * c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB4*Cbut_ * c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB5*Cbut_ *c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB6*Cbut_ * c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB7*Cbut_ * c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB8*Cbut_ *c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB9*Cbut_ *c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB10*Cbut_ * c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += 1. / (1. + Cprop_)*kp_*ABB11*Cbut_ * c_(j2 + 16);

			}

			else if (J == 5)
			{
				R_(j0) = kp_ * (
					sum1_*(1. - ABB1)
					+ (1. - FR_)*(sum2_)*(1. - ABB1)
					+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
						+ ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB3 * (SY_14)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB4 * CorrP_(3 + J - 1));

				if (J      >= LC_) R_(j0) += -kp_ * (J - 4)*c_(j0);
				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4)* (1 - ABB1);
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ + FR_ * (sum4_)
					+(sum5_)+FR_ * (sum4_)
					+0.5*(ABB2*SY1_2 + ABB3 * SY1_3 + ABB4 * SY1_4
						+ ABB5 * SY1_5
						+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
						+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
					+ 0.5*FR_*(ABB2*(SY1_13 + SY12_2)
						+ ABB3 * (SY1_14 + SY12_3)
						+ ABB4 * (SY1_15 + SY12_4)
						+ ABB5 * (SY1_16 + SY12_5)
						+ ABB6 * (SY1_17 + SY12_6)
						+ ABB7 * (SY1_18 + SY12_7)
						+ ABB8 * (SY1_19 + SY12_8)
						+ ABB9 * (SY1_20 + SY12_9)
						+ ABB10 * (SY1_21 + SY12_10)
						+ ABB11 * (SY1_22 + SY12_11))
					+ 0.5*(ABB2*(SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33))
					+ ABB1 * CorrO_(j0));

				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4);
				if (J +  7 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB6*Cbut_ * c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB7*Cbut_ *c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB8*Cbut_ * c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB9*Cbut_ *  c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB10*Cbut_ * c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*FR_*ABB11*Cbut_ *   c_(j1 + 16);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);


				R_(j2) = kp_ * (1. / 2.*sum6_ + sum7_
					+ (ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ FR_ * (ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));

				if (J +  3 >= LC_) R_(j2) += kp_ * Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += kp_ * ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += kp_ * ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += kp_ * ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += kp_ * ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += kp_ * ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += kp_ * ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += kp_ * ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += kp_ * ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += kp_ * ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += kp_ * ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += kp_ * ABB11*Cbut_*c_(j2 + 16);

			}

			else if (J == 4)
			{
				R_(j0) = kp_ * (sum1_ + (1. - FR_)*(sum2_)
					+0.5*(ABB1*SY_1 + ABB2 * SY_2
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB2 * (SY_13)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB3 * CorrP_(3 + J - 1));

				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4);
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);



				R_(j1) = kp_ * (sum3_ + FR_ * (sum4_)
					+(sum5_)+FR_ * (sum4_)
					+(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3 + ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
						+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10
						+ ABB11 * SY1_11)
					+ FR_ * (ABB1*(SY1_12 + 0.5*SY12_1)
						+ ABB2 * (SY1_13 + 0.5*SY12_2)
						+ ABB3 * (SY1_14 + 0.5*SY12_3)
						+ ABB4 * (SY1_15 + 0.5*SY12_4)
						+ ABB5 * (SY1_16 + 0.5*SY12_5)
						+ ABB6 * (SY1_17 + 0.5*SY12_6)
						+ ABB7 * (SY1_18 + 0.5*SY12_7)
						+ ABB8 * (SY1_19 + 0.5*SY12_8)
						+ ABB9 * (SY1_20 + 0.5*SY12_9)
						+ ABB10 * (SY1_21 + 0.5*SY12_10)
						+ ABB11 * (SY1_22 + 0.5*SY12_11))
					+ 0.5*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33)));


				if (J +  3 >= LC_) R_(j1) += kp_ * FR_*(Cprop_*c_(j1 + 3) + c_(j1 + 3));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*c_(j2 + 4);
				if (J +  6 >= LC_) R_(j1) += kp_ * FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += kp_ * FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += kp_ * FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += kp_ * FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += kp_ * FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += kp_ * FR_*ABB6*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += kp_ * FR_*ABB7*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += kp_ * FR_*ABB8*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += kp_ * FR_*ABB9*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += kp_ * FR_*ABB10*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += kp_ * FR_*ABB11*Cbut_*c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += 0.5*kp_*ABB1*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += 0.5*kp_*ABB2*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += 0.5*kp_*ABB3*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += 0.5*kp_*ABB4*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += 0.5*kp_*ABB5*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += 0.5*kp_*ABB6*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += 0.5*kp_*ABB7*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += 0.5*kp_*ABB8*Cbut_*c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += 0.5*kp_*ABB9*Cbut_*c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += 0.5*kp_*ABB10*Cbut_*c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += 0.5*kp_*ABB11*Cbut_*c_(j2 + 19);



				R_(j2) = kp_ * (1. / 2.*sum6_*Cbut_ + sum7_ * Cbut_
					+ (ABB1*(SY2_12)
						+ABB2 * (SY2_13)
						+ABB3 * (SY2_14)
						+ABB4 * (SY2_15)
						+ABB5 * (SY2_16)
						+ABB6 * (SY2_17)
						+ABB7 * (SY2_18)
						+ABB8 * (SY2_19)
						+ABB9 * (SY2_20)
						+ABB10 * (SY2_21)
						+ABB11 * (SY2_22))
					+ FR_ * (ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
						+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
						+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
						+ ABB10 * SY2_10 + ABB11 * SY2_11));


				if (J +  3 >= LC_) R_(j2) += kp_ * Cbut_*Cprop_*c_(j2 + 3);
				if (J +  6 >= LC_) R_(j2) += kp_ * ABB1*Cbut_*c_(j2 + 6);
				if (J +  7 >= LC_) R_(j2) += kp_ * ABB2*Cbut_*c_(j2 + 7);
				if (J +  8 >= LC_) R_(j2) += kp_ * ABB3*Cbut_*c_(j2 + 8);
				if (J +  9 >= LC_) R_(j2) += kp_ * ABB4*Cbut_*c_(j2 + 9);
				if (J + 10 >= LC_) R_(j2) += kp_ * ABB5*Cbut_*c_(j2 + 10);
				if (J + 11 >= LC_) R_(j2) += kp_ * ABB6*Cbut_*c_(j2 + 11);
				if (J + 12 >= LC_) R_(j2) += kp_ * ABB7*Cbut_*c_(j2 + 12);
				if (J + 13 >= LC_) R_(j2) += kp_ * ABB8*Cbut_*c_(j2 + 13);
				if (J + 14 >= LC_) R_(j2) += kp_ * ABB9*Cbut_*c_(j2 + 14);
				if (J + 15 >= LC_) R_(j2) += kp_ * ABB10*Cbut_*c_(j2 + 15);
				if (J + 16 >= LC_) R_(j2) += kp_ * ABB11*Cbut_*c_(j2 + 16);

			}

			else if (J == 3)
			{
				R_(j0) = kp_ * (sum1_ + (1. - FR_)*(sum2_)
					+0.5*(ABB1*SY_1 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(ABB1*(SY_12)
						+ABB3 * (SY_14)
						+ABB4 * (SY_15)
						+ABB5 * (SY_16)
						+ABB6 * (SY_17)
						+ABB7 * (SY_18)
						+ABB8 * (SY_19)
						+ABB9 * (SY_20)
						+ABB10 * (SY_21)
						+ABB11 * (SY_22))
					+ ABB2 * CorrP_(3 + J - 1));

				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4);
				if (J +  9 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB1*Cbut_*c_(j1 + 9);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);


				R_(j1) = kp_ * (sum3_ + FR_ * (Cprop_ + 1.)*sum4_ +
					(sum5_)*Cprop_
					+ (ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3 + ABB4 * SY1_4
						+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
						+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10
						+ ABB11 * SY1_11)
					+ FR_ * (ABB1*(SY1_12 + Cprop_ / (1. + Cprop_)*SY12_1)
						+ ABB2 * (SY1_13 + Cprop_ / (1. + Cprop_)*SY12_2)
						+ ABB3 * (SY1_14 + Cprop_ / (1. + Cprop_)*SY12_3)
						+ ABB4 * (SY1_15 + Cprop_ / (1. + Cprop_)*SY12_4)
						+ ABB5 * (SY1_16 + Cprop_ / (1. + Cprop_)*SY12_5)
						+ ABB6 * (SY1_17 + Cprop_ / (1. + Cprop_)*SY12_6)
						+ ABB7 * (SY1_18 + Cprop_ / (1. + Cprop_)*SY12_7)
						+ ABB8 * (SY1_19 + Cprop_ / (1. + Cprop_)*SY12_8)
						+ ABB9 * (SY1_20 + Cprop_ / (1. + Cprop_)*SY12_9)
						+ ABB10 * (SY1_21 + Cprop_ / (1. + Cprop_)*SY12_10)
						+ ABB11 * (SY1_22 + Cprop_ / (1. + Cprop_)*SY12_11))
					+ Cprop_ / (1. + Cprop_)*(ABB1*(SY1_23)
						+ABB2 * (SY1_24)
						+ABB3 * (SY1_25)
						+ABB4 * (SY1_26)
						+ABB5 * (SY1_27)
						+ABB6 * (SY1_28)
						+ABB7 * (SY1_29)
						+ABB8 * (SY1_30)
						+ABB9 * (SY1_31)
						+ABB10 * (SY1_32)
						+ABB11 * (SY1_33)));


				if (J +  3 >= LC_) R_(j1) += kp_ * (Cprop_*c_(j1 + 3));
				if (J +  4 >= LC_) R_(j1) += kp_ * Cbut_*Cprop_*c_(j2 + 4);
				if (J +  6 >= LC_) R_(j1) += kp_ * FR_*ABB1*Cbut_*c_(j1 + 6);
				if (J +  7 >= LC_) R_(j1) += kp_ * FR_*ABB2*Cbut_*c_(j1 + 7);
				if (J +  8 >= LC_) R_(j1) += kp_ * FR_*ABB3*Cbut_*c_(j1 + 8);
				if (J +  9 >= LC_) R_(j1) += kp_ * FR_*ABB4*Cbut_*c_(j1 + 9);
				if (J + 10 >= LC_) R_(j1) += kp_ * FR_*ABB5*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j1) += kp_ * FR_*ABB6*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j1) += kp_ * FR_*ABB7*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j1) += kp_ * FR_*ABB8*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j1) += kp_ * FR_*ABB9*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j1) += kp_ * FR_*ABB10*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j1) += kp_ * FR_*ABB11*Cbut_*c_(j1 + 16);
				if (J +  9 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB1*Cbut_* c_(j2 + 9);
				if (J + 10 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB2*Cbut_* c_(j2 + 10);
				if (J + 11 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB3*Cbut_ * c_(j2 + 11);
				if (J + 12 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB4*Cbut_* c_(j2 + 12);
				if (J + 13 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB5*Cbut_ * c_(j2 + 13);
				if (J + 14 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB6*Cbut_* c_(j2 + 14);
				if (J + 15 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB7*Cbut_* c_(j2 + 15);
				if (J + 16 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB8*Cbut_* c_(j2 + 16);
				if (J + 17 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB9*Cbut_* c_(j2 + 17);
				if (J + 18 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB10*Cbut_* c_(j2 + 18);
				if (J + 19 >= LC_) R_(j1) += Cprop_ / (1. + Cprop_)*kp_*ABB11*Cbut_* c_(j2 + 19);

			}

			else if (J == 2)
			{
				R_(j0) = kp_ * (sum1_ + (1. - FR_)*(sum2_)
					+0.5*(ABB2*SY_2 + ABB3 * SY_3
						+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
						+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
						+ ABB10 * SY_10 + ABB11 * SY_11)
					+ 0.5*FR_*(
									ABB2*SY_13
									+ABB3 * SY_14
									+ABB4 * SY_15
									+ABB5 * SY_16
									+ABB6 * SY_17
									+ABB7 * SY_18
									+ABB8 * SY_19
									+ABB9 * SY_20
									+ABB10 * SY_21
									+ABB11 * SY_22 )
					+ ABB1 * CorrP_(3 + J - 1));


				if (J +  4 >= LC_) R_(j0) += kp_ * (1. - FR_)*Cbut_*c_(j1 + 4);
				if (J + 10 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB2*Cbut_*c_(j1 + 10);
				if (J + 11 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB3*Cbut_*c_(j1 + 11);
				if (J + 12 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB4*Cbut_*c_(j1 + 12);
				if (J + 13 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB5*Cbut_*c_(j1 + 13);
				if (J + 14 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB6*Cbut_*c_(j1 + 14);
				if (J + 15 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB7*Cbut_*c_(j1 + 15);
				if (J + 16 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB8*Cbut_*c_(j1 + 16);
				if (J + 17 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB9*Cbut_*c_(j1 + 17);
				if (J + 18 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB10*Cbut_*c_(j1 + 18);
				if (J + 19 >= LC_) R_(j0) += 0.5*kp_*FR_*ABB11*Cbut_*c_(j1 + 19);

			}

		} 

		if (is_verbose_ == true)
		{
			const double V = 0.99 / LiquidDensity(T_) / 1.006674894;

			std::ofstream fLiq("Liq.class.out", std::ios::out);
			fLiq.setf(std::ios::scientific);
			for (int J = NPA_; J >= LC_; J--)
			{
				fLiq << std::setw(12) << std::right << J;
				fLiq << std::setw(24) << std::right << std::setprecision(15) << c_(J - 1)*V;
				fLiq << std::setw(24) << std::right << std::setprecision(15) << R_(J - 1)*V;
				fLiq << std::endl;
			}
			fLiq.close();
		}

		if (is_verbose_ == true)
		{
			const double V = 0.99 / LiquidDensity(T_) / 1.006674894;

			std::ofstream fGas("Gas.class.out", std::ios::out);
			fGas.setf(std::ios::scientific);
			for (int J = LC_ - 1; J >= 2; J--)
			{
				fGas << std::setw(12) << std::right << J;
				fGas << std::setw(24) << std::right << std::setprecision(15) << c_(J - 1)*V;
				fGas << std::setw(24) << std::right << std::setprecision(15) << R_(J - 1)*V;
				fGas << std::endl;
			}
			fGas.close();
		}
	}

	void PolyethyleneKinetics::UpdateCorrections()
	{
		double SOMXX = c_(NPA_ - 1);
		double SOMZZ1 = c_(2 * NPA_ - 1);
		double SOMZZ2 = c_(2 * NPA_ - 1);
		double SOMDD = c_(3 * NPA_ - 1);
		double ASOMXX = SOMXX;
		double ASOMZZ1 = SOMZZ1;
		double ASOMZZ2 = SOMZZ2;
		double ASOMDD = SOMDD;

		int IVAL = 15;

		for (int jj = NPA_ - 4; jj >= IVAL + 2; jj--)
		{
			if (jj + 3 >= LC_)
				ASOMXX += c_(jj + 3 - 1);
			SOMXX += ASOMXX;
		}

		for (int jj = NPA_ - 6; jj >= IVAL + 2; jj--)
		{
			if (jj + 5 >= LC_)
				ASOMZZ1 += c_(jj + 5 + NPA_ - 1);
			SOMZZ1 += ASOMZZ1;
		}

		for (int jj = NPA_ - 4; jj >= IVAL + 4; jj--)
		{
			if (jj + 3 >= LC_)
				ASOMZZ2 += c_(jj + 3 + NPA_);
			SOMZZ2 += ASOMZZ2;
		}

		for (int jj = NPA_ - 6; jj >= IVAL + 4; jj--)
		{
			if (jj + 5 >= LC_)
				ASOMDD += c_(2 * NPA_ + jj + 5 - 1);
			SOMDD += ASOMDD;
		}


		double SOMP1 = 0.0;
		for (int j = 2 * IVAL; j <= NPA_; j++)
		{
			if (j >= LC_)
				SOMP1 += c_(j - 1);
		}

		double SOMD1 = 0.0;
		for (int j = IVAL + 1; j <= NPA_; j++)
		{
			if (j >= LC_)
				SOMD1 += c_(j + 2 * NPA_ - 1);
		}

		double SOMD2 = 0.0;
		for (int j = 2 * IVAL + 5; j <= NPA_; j++)
		{
			if (j >= LC_)
				SOMD2 += c_(j + 2 * NPA_ - 1);
		}

		double SOMO1 = 0.0;
		for (int j = IVAL + 6; j <= NPA_; j++)
		{
			if (j >= LC_)
				SOMO1 += c_(j + NPA_ - 1);
		}

		double SOMO2 = 0.0;
		for (int j = 2 * IVAL + 3; j <= NPA_; j++)
		{
			if (j >= LC_)
				SOMO2 += c_(j + NPA_ - 1);
		}

		double SOMP2 = 0.0;
		for (int j = IVAL + 4; j <= NPA_; j++)
		{
			if (j >= LC_)
				SOMP2 += c_(j - 1);
		}


		CorrO_.setZero();
		CorrP_.setZero();

		CorrO_(IVAL - 1) = 0.5*(SOMXX + SOMP1) + 0.5*(SOMDD + Cbut_ * SOMD1) +
			Cprop_ / (1. + Cprop_)*SOMD1 + 0.5*SOMD2 +
			0.5*(0.5*(SOMZZ1 + Cbut_ * SOMO1) + 0.5*SOMO2) +
			0.5*(0.5*SOMZZ2 + Cprop_ / (1. + Cprop_)*SOMO1 + 0.5*SOMO2);


		if (IVAL + 7 >= LC_)
			CorrO_(IVAL - 1) += Cprop_ / (1. + Cprop_)*Cbut_*c_(2 * NPA_ + IVAL + 7 - 1);

		if (2 * IVAL + 4 >= LC_)
			CorrO_(IVAL - 1) += 0.5*Cbut_*c_(2 * NPA_ + 2 * IVAL + 4 - 1);

		if (2 * IVAL + 1 >= LC_)
			CorrO_(IVAL - 1) += 0.25*Cbut_*c_(2 * IVAL + 1 + NPA_ - 1);

		if (2 * IVAL + 2 >= LC_)
			CorrO_(IVAL - 1) += 0.25*c_(2 * IVAL + 2 + NPA_ - 1);

		CorrP_(IVAL - 1) = 0.5*SOMXX + (SOMP2 + SOMP2 + 0.5*SOMP1) +
			0.5*(SOMDD + Cbut_ * SOMD1) + 1. / (1. + Cprop_) * SOMD1 + SOMD1 + SOMD1 +
			0.5*(0.5*(SOMZZ1 + Cbut_ * SOMO1) + 0.5*SOMO2 + SOMO1 + SOMO1) +
			0.5*(0.5*SOMZZ2 + 1. / (1. + Cprop_)*SOMO1 + SOMO1 + SOMO1);

		if (IVAL + 3 >= LC_)
			CorrP_(IVAL - 1) += c_(IVAL + 3 - 1);

		if (IVAL + 7 >= LC_)
			CorrP_(IVAL - 1) += Cbut_ / (1. + Cprop_)*c_(IVAL + 7 + 2 * NPA_ - 1);

		if (IVAL + 6 >= LC_)
			CorrP_(IVAL - 1) += Cbut_ * c_(2 * NPA_ + IVAL + 6 - 1);

		if (IVAL + 7 >= LC_)
			CorrP_(IVAL - 1) += c_(2 * NPA_ + IVAL + 7 - 1);

		if (IVAL + 5 >= LC_)
			CorrP_(IVAL - 1) += Cbut_ * c_(2 * NPA_ + IVAL + 5 - 1);

		if (IVAL + 7 >= LC_)
			CorrP_(IVAL - 1) += c_(2 * NPA_ + IVAL + 7 - 1);

		if (IVAL + 6 >= LC_)
			CorrP_(IVAL - 1) += c_(2 * NPA_ + IVAL + 6 - 1);

		if (2 * IVAL + 1 >= LC_)
			CorrP_(IVAL - 1) += 0.25*Cbut_*c_(NPA_ + 2 * IVAL + 1 - 1);

		if (2 * IVAL + 2 >= LC_)
			CorrP_(IVAL - 1) += 0.25*c_(NPA_ + 2 * IVAL + 2 - 1);

		if (IVAL + 5 >= LC_)
			CorrP_(IVAL - 1) += 0.5*Cbut_*c_(NPA_ + IVAL + 5 - 1);

		if (IVAL + 4 >= LC_)
			CorrP_(IVAL - 1) += 0.5*Cbut_*c_(NPA_ + IVAL + 4 - 1);

		if (IVAL + 5 >= LC_)
			CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 5 - 1);

		if (IVAL + 5 >= LC_)
			CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 5 - 1);

		if (IVAL + 4 >= LC_)
			CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 4 - 1);

		if (IVAL + 5 >= LC_)
			CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 5 - 1);

		for (int IVAL = 14; IVAL >= 5; IVAL--)
		{
			if (2 * IVAL + 1 >= LC_)
				SOMP1 += c_(2 * IVAL + 1 - 1);

			if (2 * IVAL >= LC_)
				SOMP1 += c_(2 * IVAL - 1);

			if (IVAL + 4 >= LC_)
				SOMP2 += c_(IVAL + 4 - 1);

			if (IVAL + 8 >= LC_)
				SOMD1 += c_(2 * NPA_ + IVAL + 8 - 1);

			if (2 * IVAL + 6 >= LC_)
				SOMD2 += c_(2 * NPA_ + 2 * IVAL + 5 + 1 - 1);

			if (2 * IVAL + 5 >= LC_)
				SOMD2 += c_(2 * NPA_ + 2 * IVAL + 5 - 1);

			if (IVAL + 6 >= LC_)
				SOMO1 += c_(NPA_ + IVAL + 6 - 1);

			if (2 * IVAL + 4 >= LC_)
				SOMO2 += c_(NPA_ + 2 * IVAL + 3 + 1 - 1);

			if (2 * IVAL + 3 >= LC_)
				SOMO2 += c_(NPA_ + 2 * IVAL + 3 - 1);

			int JJ = IVAL + 2;

			if (JJ + 3 >= LC_)
				ASOMXX += c_(JJ + 3 - 1);
			SOMXX += ASOMXX;

			if (JJ + 5 >= LC_)
				ASOMZZ1 += c_(JJ + NPA_ + 5 - 1);
			SOMZZ1 += ASOMZZ1;

			JJ = IVAL + 4;

			if (JJ + 3 >= LC_)
				ASOMZZ2 += c_(JJ + NPA_ + 3 - 1);
			SOMZZ2 += ASOMZZ2;

			if (JJ + 5 >= LC_)
				ASOMDD += c_(JJ + 2 * NPA_ + 5 - 1);
			SOMDD += ASOMDD;

			CorrO_(IVAL - 1) = 0.5*(SOMXX + SOMP1) + 0.5*(SOMDD + Cbut_ * SOMD1) +
				Cprop_ / (1. + Cprop_)*SOMD1 + 0.5*SOMD2 + 0.5*(0.5*(SOMZZ1 + Cbut_ * SOMO1) +
					0.5*SOMO2) + 0.5*(0.5*SOMZZ2 + Cprop_ / (1. + Cprop_)*SOMO1 + 0.5*SOMO2);

			if (IVAL + 7 >= LC_)
				CorrO_(IVAL - 1) += Cprop_ / (1. + Cprop_)*Cbut_*c_(2 * NPA_ + IVAL + 7 - 1);

			if (2 * IVAL + 4 >= LC_)
				CorrO_(IVAL - 1) += 0.5*Cbut_*c_(2 * NPA_ + 2 * IVAL + 4 - 1);

			if (2 * IVAL + 1 >= LC_)
				CorrO_(IVAL - 1) += 0.25*Cbut_*c_(2 * IVAL + 1 + NPA_ - 1);


			if (2 * IVAL + 2 >= LC_)
				CorrO_(IVAL - 1) += 0.25*c_(2 * IVAL + 2 + NPA_ - 1);


			CorrP_(IVAL - 1) = 0.5*SOMXX + (SOMP2 + SOMP2 + 0.5*SOMP1) +
				0.5*(SOMDD + Cbut_ * SOMD1) + 1. / (1. + Cprop_) *
				SOMD1 + SOMD1 + SOMD1 +
				0.5*(0.5*(SOMZZ1 + Cbut_ * SOMO1) + 0.5*SOMO2 + SOMO1 + SOMO1) +
				0.5*(0.5*SOMZZ2 + 1. / (1. + Cprop_)*SOMO1 + SOMO1 + SOMO1);


			if (IVAL + 3 >= LC_)
				CorrP_(IVAL - 1) += c_(IVAL + 3 - 1);

			if (IVAL + 7 >= LC_)
				CorrP_(IVAL - 1) += Cbut_ / (1. + Cprop_)*c_(IVAL + 7 + 2 * NPA_ - 1);

			if (IVAL + 6 >= LC_)
				CorrP_(IVAL - 1) += Cbut_ * c_(2 * NPA_ + IVAL + 6 - 1);

			if (IVAL + 7 >= LC_)
				CorrP_(IVAL - 1) += c_(2 * NPA_ + IVAL + 7 - 1);

			if (IVAL + 5 >= LC_)
				CorrP_(IVAL - 1) += Cbut_ * c_(2 * NPA_ + IVAL + 5 - 1);

			if (IVAL + 7 >= LC_)
				CorrP_(IVAL - 1) += c_(2 * NPA_ + IVAL + 7 - 1);

			if (IVAL + 6 >= LC_)
				CorrP_(IVAL - 1) += c_(2 * NPA_ + IVAL + 6 - 1);

			if (2 * IVAL + 1 >= LC_)
				CorrP_(IVAL - 1) += 0.25*Cbut_*c_(NPA_ + 2 * IVAL + 1 - 1);

			if (2 * IVAL + 2 >= LC_)
				CorrP_(IVAL - 1) += 0.25*c_(NPA_ + 2 * IVAL + 2 - 1);

			if (IVAL + 5 >= LC_)
				CorrP_(IVAL - 1) += 0.5*Cbut_*c_(NPA_ + IVAL + 5 - 1);

			if (IVAL + 4 >= LC_)
				CorrP_(IVAL - 1) += 0.5*Cbut_*c_(NPA_ + IVAL + 4 - 1);

			if (IVAL + 5 >= LC_)
				CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 5 - 1);

			if (IVAL + 5 >= LC_)
				CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 5 - 1);

			if (IVAL + 4 >= LC_)
				CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 4 - 1);

			if (IVAL + 5 >= LC_)
				CorrP_(IVAL - 1) += 0.5*c_(NPA_ + IVAL + 5 - 1);

		}
	}

	void PolyethyleneKinetics::UpdateSums(const int J)
	{

		if ((NPA_ - 2 >= J) && (J + 2 >= LC_))
		{
			const int M0 = J + 1;
			const int M1 = M0 + NPA_;
			sum3_ += c_(M0);
			sum6_ += c_(M1);
		}

		if ((NPA_ - 3 >= J) && (J + 3 >= LC_))
		{
			const int M0 = J + 2;
			sum1_ += c_(M0);
		}

		if ((NPA_ - 4 >= J) && (J + 4 >= LC_))
		{
			const int M0 = J + 3;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sum4_ += c_(M1);
			sum7_ += c_(M2);
		}

		if ((NPA_ - 5 >= J) && (J + 5 >= LC_))
		{
			const int M0 = J + 4;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY1_(0) += c_(M0);
			sumY2_(0) += c_(M1);
			sum2_ += c_(M1);
			sum5_ += c_(M2);
		}

		if ((NPA_ - 6 >= J) && (J + 6 >= LC_))
		{
			const int M0 = J + 5;
			const int M1 = M0 + NPA_;
			sumY2_(1) += c_(M1);
			sumY1_(1) += c_(M0);
		}

		if ((NPA_ - 7 >= J) && (J + 7 >= LC_))
		{
			const int M0 = J + 6;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY1_(2) += c_(M0);
			sumY1_(11) += c_(M1);
			sumY2_(2) += c_(M1);
			sumY2_(11) += c_(M2);
		}

		if ((NPA_ - 8 >= J) && (J + 8 >= LC_))
		{
			const int M0 = J + 7;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(0) += c_(M0);
			sumY1_(3) += c_(M0);
			sumY2_(3) += c_(M1);
			sumY12_(0) += c_(M1);
			sumY1_(12) += c_(M1);
			sumY2_(12) += c_(M2);
		}

		if ((NPA_ - 9 >= J) && (J + 9 >= LC_))
		{
			const int M0 = J + 8;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY2_(13) += c_(M2);
			sumY1_(13) += c_(M1);
			sumY12_(1) += c_(M1);
			sumY_(1) += c_(M0);
			sumY1_(4) += c_(M0);
			sumY2_(4) += c_(M1);
		}

		if ((NPA_ - 10 >= J) && (J + 10 >= LC_))
		{
			const int M0 = J + 9;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(2) += c_(M0);
			sumY_(11) += c_(M1);
			sumY1_(14) += c_(M1);
			sumY12_(2) += c_(M1);
			sumY1_(22) += c_(M2);
			sumY2_(14) += c_(M2);
			sumY1_(5) += c_(M0);
			sumY2_(5) += c_(M1);
		}

		if ((NPA_ - 11 >= J) && (J + 11 >= LC_))
		{
			const int M0 = J + 10;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(3) += c_(M0);
			sumY12_(3) += c_(M1);
			sumY_(12) += c_(M1);
			sumY1_(23) += c_(M2);
			sumY2_(15) += c_(M2);
			sumY1_(15) += c_(M1);
			sumY1_(6) += c_(M0);
			sumY2_(6) += c_(M1);
		}

		if ((NPA_ - 12 >= J) && (J + 12 >= LC_))
		{
			const int M0 = J + 11;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(4) += c_(M0);
			sumY12_(4) += c_(M1);
			sumY_(13) += c_(M1);
			sumY1_(24) += c_(M2);
			sumY2_(16) += c_(M2);
			sumY1_(16) += c_(M1);
			sumY1_(7) += c_(M0);
			sumY2_(7) += c_(M1);
		}

		if ((NPA_ - 13 >= J) && (J + 13 >= LC_))
		{
			const int M0 = J + 12;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(5) += c_(M0);
			sumY12_(5) += c_(M1);
			sumY_(14) += c_(M1);
			sumY1_(25) += c_(M2);
			sumY2_(17) += c_(M2);
			sumY1_(17) += c_(M1);
			sumY1_(8) += c_(M0);
			sumY2_(8) += c_(M1);
		}

		if ((NPA_ - 14 >= J) && (J + 14 >= LC_))
		{
			const int M0 = J + 13;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY1_(26) += c_(M2);
			sumY_(15) += c_(M1);
			sumY_(6) += c_(M0);
			sumY12_(6) += c_(M1);
			sumY2_(18) += c_(M2);
			sumY1_(18) += c_(M1);
			sumY1_(9) += c_(M0);
			sumY2_(9) += c_(M1);
		}

		if ((NPA_ - 15 >= J) && (J + 15 >= LC_))
		{
			const int M0 = J + 14;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(16) += c_(M1);
			sumY1_(27) += c_(M2);
			sumY_(7) += c_(M0);
			sumY12_(7) += c_(M1);
			sumY2_(19) += c_(M2);
			sumY1_(19) += c_(M1);
			sumY1_(10) += c_(M0);
			sumY2_(10) += c_(M1);
		}

		if ((NPA_ - 16 >= J) && (J + 16 >= LC_))
		{
			const int M0 = J + 15;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(17) += c_(M1);
			sumY1_(28) += c_(M2);
			sumY_(8) += c_(M0);
			sumY12_(8) += c_(M1);
			sumY1_(20) += c_(M1);
			sumY2_(20) += c_(M2);
		}

		if ((NPA_ - 17 >= J) && (J + 17 >= LC_))
		{
			const int M0 = J + 16;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(18) += c_(M1);
			sumY1_(29) += c_(M2);
			sumY_(9) += c_(M0);
			sumY12_(9) += c_(M1);
			sumY2_(21) += c_(M2);
			sumY1_(21) += c_(M1);
		}

		if ((NPA_ - 18 >= J) && (J + 18 >= LC_))
		{
			const int M0 = J + 17;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(19) += c_(M1);
			sumY1_(30) += c_(M2);
			sumY_(10) += c_(M0);
			sumY12_(10) += c_(M1);
		}

		if ((NPA_ - 19 >= J) && (J + 19 >= LC_))
		{
			const int M0 = J + 18;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(20) += c_(M1);
			sumY1_(31) += c_(M2);
		}

		if ((NPA_ - 20 >= J) && (J + 20 >= LC_))
		{
			const int M0 = J + 19;
			const int M1 = M0 + NPA_;
			const int M2 = M0 + NPA_ * 2;
			sumY_(21) += c_(M1);
			sumY1_(32) += c_(M2);
		}
	}

	double PolyethyleneKinetics::Sum(const Eigen::VectorXd& v) const
	{
		return SumGas(v) + SumLiquid(v);
	}

	double PolyethyleneKinetics::SumMW(const Eigen::VectorXd& v) const
	{
		return SumGasMW(v) + SumLiquidMW(v);
	}

	void PolyethyleneKinetics::Sum(const Eigen::VectorXd& v, double& P, double& O, double& D) const
	{
		double Pg, Og, Dg;
		double Pl, Ol, Dl;
		SumGas(v, Pg, Og, Dg);
		SumLiquid(v, Pl, Ol, Dl);
		P = Pg + Pl;
		O = Og + Ol;
		D = Dg + Dl;
	}

	void PolyethyleneKinetics::SumMW(const Eigen::VectorXd& v, double& P, double& O, double& D) const
	{
		double Pg, Og, Dg;
		double Pl, Ol, Dl;
		SumGasMW(v, Pg, Og, Dg);
		SumLiquidMW(v, Pl, Ol, Dl);
		P = Pg + Pl;
		O = Og + Ol;
		D = Dg + Dl;
	}

	

	double PolyethyleneKinetics::SumGas(const Eigen::VectorXd& v) const
	{
		double sum = 0.;
		for (int i = 1; i <= LC_ - 1; i++)
		{
			sum += v(i - 1);
			sum += v(i - 1 + NPA_);
			sum += v(i - 1 + NPA_ + NOL_);
		}

		return sum;
	}

	double PolyethyleneKinetics::SumGasMW(const Eigen::VectorXd& v) const
	{
		double sum = 0.;
		for (int i = 1; i <= LC_ - 1; i++)
		{
			const double MW = i * MWm_;

			sum += v(i - 1)*(MW + 2.);
			sum += v(i - 1 + NPA_)*(MW);
			sum += v(i - 1 + NPA_ + NOL_)*(MW - 2.);
		}

		return sum;
	}

	void PolyethyleneKinetics::SumGas(const Eigen::VectorXd& v, double& P, double& O, double& D) const
	{
		P = 0.;
		O = 0.;
		D = 0.;

		for (int i = 1; i <= LC_ - 1; i++)
		{
			P += v(i - 1);
			O += v(i - 1 + NPA_);
			D += v(i - 1 + NPA_ + NOL_);
		}
	}

	void PolyethyleneKinetics::SumGasMW(const Eigen::VectorXd& v, double& P, double& O, double& D) const
	{
		P = 0.;
		O = 0.;
		D = 0.;

		for (int i = 1; i <= LC_ - 1; i++)
		{
			const double MW = i * MWm_;

			P += v(i - 1)*(MW + 2.);
			O += v(i - 1 + NPA_)*(MW);
			D += v(i - 1 + NPA_ + NOL_)*(MW - 2.);
		}
	}

	double PolyethyleneKinetics::SumLiquid(const Eigen::VectorXd& v) const 
	{
		double sum = 0.;
		for (int i = NPA_; i >= LC_; i--)
		{
			sum += v(i - 1);
			sum += v(i - 1 + NPA_);
			sum += v(i - 1 + NPA_ + NOL_);
		}

		return sum;
	}

	double PolyethyleneKinetics::SumLiquidMW(const Eigen::VectorXd& v) const
	{
		double sum = 0.;
		for (int i = NPA_; i >= LC_; i--)
		{
			const double MW = i * MWm_;

			sum += v(i - 1)*(MW + 2.);
			sum += v(i - 1 + NPA_)*(MW);
			sum += v(i - 1 + NPA_ + NOL_)*(MW - 2.);
		}

		return sum;
	}

	void PolyethyleneKinetics::SumLiquid(const Eigen::VectorXd& v, double& P, double& O, double& D) const
	{
		P = 0.;
		O = 0.;
		D = 0.;

		for (int i = NPA_; i >= LC_; i--)
		{
			P += v(i - 1);
			O += v(i - 1 + NPA_);
			D += v(i - 1 + NPA_ + NOL_);
		}
	}

	void PolyethyleneKinetics::SumLiquidMW(const Eigen::VectorXd& v, double& P, double& O, double& D) const
	{
		P = 0.;;
		O = 0.;;
		D = 0.;;
	
		double sum = 0.;
		for (int i = NPA_; i >= LC_; i--)
		{
			const double MW = i * MWm_;

			P += v(i - 1)*(MW + 2.);
			O += v(i - 1 + NPA_)*(MW);
			D += v(i - 1 + NPA_ + NOL_)*(MW - 2.);
		}
	}

	void PolyethyleneKinetics::Distribution(const Eigen::VectorXd& v, Eigen::VectorXd& P, Eigen::VectorXd& O, Eigen::VectorXd& D) const
	{
		P.resize(NPA_);
		O.resize(NOL_);
		D.resize(NDO_);

		for (int i = 0; i < NPA_; i++)
			P(i) = v(i);

		for (int i = 0; i < NOL_; i++)
			O(i) = v(i + NPA_);

		for (int i = 0; i < NDO_; i++)
			D(i) = v(i + NPA_ + NOL_);
	}

	void PolyethyleneKinetics::DistributionMW(const Eigen::VectorXd& v, Eigen::VectorXd& P, Eigen::VectorXd& O, Eigen::VectorXd& D) const
	{
		P.resize(NPA_);
		O.resize(NOL_);
		D.resize(NDO_);

		for (int i = 0; i < NPA_; i++)
			P(i) = v(i)*(MWm_*(i+1.)+2.);

		for (int i = 0; i < NOL_; i++)
			O(i) = v(i + NPA_)*(MWm_*(i + 1.));

		for (int i = 0; i < NDO_; i++)
			D(i) = v(i + NPA_ + NOL_)*(MWm_*(i + 1.) - 2.);
	}

}
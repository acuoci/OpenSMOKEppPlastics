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

#include "PolypropyleneKinetics.h"

namespace opensmokepp::plastics
{

	PolypropyleneKinetics::~PolypropyleneKinetics()
	{
	}

	PolypropyleneKinetics::PolypropyleneKinetics()
	{
		// ------------------------------------------------------------------------------------------------
		// Default values from kinetics
		// ------------------------------------------------------------------------------------------------
		DefaultKineticConstants();

		// ------------------------------------------------------------------------------------------------
		// Default values from initial distribution
		// ------------------------------------------------------------------------------------------------
		MWm_ = 42.;
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

	void PolypropyleneKinetics::DefaultKineticConstants()
	{
		A_.resize(9);
		Beta_.resize(9);
		E_.resize(9);

		// ------------------------------------------------------------------------------------------------
		// 1. Initiation reactions
		// ------------------------------------------------------------------------------------------------
		A_(0) = 7.94328e14;		Beta_(0) = 0.0;		E_(0) = 73000.0;	// initiation: random scission
		A_(1) = 4.00e13;		Beta_(1) = 0.0;		E_(1) = 65000.0;	// initiation: allylic scission


	   // ------------------------------------------------------------------------------------------------
	   // 2. Termination reactions
	   // ------------------------------------------------------------------------------------------------
		A_(2) = 2.5e7;			Beta_(2) = 1.0;		E_(2) = 6000.0;		// termination


		// ------------------------------------------------------------------------------------------------
		// 3. Abstractions
		// ------------------------------------------------------------------------------------------------
		A_(3) = 4.74341649e8;	Beta_(3) = 0.0;		E_(3) = 11000.0;	// inter-molecuar abstraction
		A_(4) = 4.74341649e8;	Beta_(4) = 0.0;		E_(4) = 13500.0;	// re-abstraction


		// ------------------------------------------------------------------------------------------------
		// 4. Beta-decompositions
		// ------------------------------------------------------------------------------------------------
		A_(5) = 1.25892541e14;	Beta_(5) = 0.0;		E_(5) = 30500.0;	// beta-decomposition


		// ------------------------------------------------------------------------------------------------
		// 5. Back-biting (or isomerization)
		// ------------------------------------------------------------------------------------------------
		A_(6) = 1.00e10;		Beta_(6) = 0.0;		E_(6) = 20600.0;	// back-biting abstraction 1-4
		A_(7) = 1.60e9;			Beta_(7) = 0.0;		E_(7) = 14500.0;	// back-biting abstraction 1-5
		A_(8) = 5.00e8;			Beta_(8) = 0.0;		E_(8) = 14500.0;	// back-biting abstraction 1-6
	}

	double PolypropyleneKinetics::LiquidDensity(const double T) const
	{
		return 760.;
	}

	double PolypropyleneKinetics::BoilingTemperature(const double L, const double P) const
	{
		return std::sqrt(3.*(L + 1.) / 5.46e-5) / (1. - std::log(P) / 10.5);
	}

	int PolypropyleneKinetics::MinNumberOfUnits(const double T, const double P)
	{
		return static_cast<int>( 1./3.*5.46e-5*std::pow(T, 2.)*std::pow(1. - std::log(P) / 10.5, 2.) );
	}

}
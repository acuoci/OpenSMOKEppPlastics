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

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace opensmokepp::plastics
{
	/// @brief Type of initial distribution of chains
	/// @details SCHULTZ is the result of a radicalic polimerization, while
	///          SCHULTZ-FLORY is the result of a ionic polimerization
	enum PolymerDistribution { DIRAC_DELTA, SCHULTZ, SCHULTZ_FLORY };

	double ConversionActivationEnergy(const double E, const std::string units);

	double ConversionFrequencyFactor(const double A, const std::string units);

	void InitialDistribution(const PolymerDistribution distribution,
		const double MW_monomer, const double MWE,
		const int lumping_start, const int lumping_step,
		const double epsilon, const int stretching_coefficient,
		Eigen::VectorXd& y, int& N, const std::string folder_name);

	void LumpingSetup(	const bool is_lumping, const int lumping_start, const int lumping_step,
						const double MW_monomer, Eigen::VectorXd& y, int& N, const std::string folder_name);

	void InitialDistribution(Eigen::VectorXd& y, const double epsi, const double MWm, const double MWp,
		int& N, double& wg, const std::string folder_name);

	int SearchForLC(std::vector<double>& list_boiling_temperature, const double T);


}


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

#ifndef OpenSMOKE_ThermogravimetricAnalysis_H
#define OpenSMOKE_ThermogravimetricAnalysis_H

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

// Grammar
#include "Grammar_ThermogravimetricAnalysis.h"

class ThermogravimetricAnalysis
{

public:

	ThermogravimetricAnalysis();

	void operator() (const std::string input_file_name, const std::string dictionary_name);

	double T(const double t) { return temperature_profile_.Interpolate(t); }

	double P(const double t) { return P_; }

	double FinalTime() const { return tEnd_; }

	double InitialMass() const { return initial_mass_; }

	bool is_verbose() const { return is_verbose_; }

	double BoilingTemperature(const int i) const { return list_boiling_temperature_[i]; }

	int SearchForLC(const double T);

	void SetInitialMass(const double initial_mass) { initial_mass_ = initial_mass; }

	void SetFinallTime(const double tEnd) { tEnd_ = tEnd; }

	void SetVerbose(const double is_verbose) { is_verbose_ = is_verbose; }

	void SetBoilingTemperatureList(const std::vector<double>& list_boiling_temperature) { list_boiling_temperature_ = list_boiling_temperature; }

	OpenSMOKE::PolymerType Polymer() const { return polymer_type_; }

	const OpenSMOKE::ThermogravimetricAnalysis_Options& options() const { return *options_; }

	const OdeSMOKE::ExplicitOdeSolver_Parameters& ode_options() const { return *ode_options_; }


private:

	bool is_verbose_;
	double tEnd_;
	double P_;
	double initial_mass_;
	OpenSMOKE::FixedProfile temperature_profile_;
	std::vector<double> list_boiling_temperature_;

	OpenSMOKE::PolymerType polymer_type_;

	OpenSMOKE::ThermogravimetricAnalysis_Options* options_;

	OdeSMOKE::ExplicitOdeSolver_Parameters* ode_options_;

private:

	void FatalErrorMessage(const std::string message);
};

#include "ThermogravimetricAnalysis.hpp"

#endif // OpenSMOKE_ThermogravimetricAnalysis_H

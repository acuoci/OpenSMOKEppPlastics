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

#ifndef OpenSMOKE_ExplicitOdeSolver_Parameters_H
#define	OpenSMOKE_ExplicitOdeSolver_Parameters_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <math/OpenSMOKEFunctions.h>
#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OdeSMOKE
{
	class ExplicitOdeSolver_Parameters
	{
	public:

		enum ODE_INTEGRATOR { ODE_RUNGEKUTTA45 };

	public:
	
		/**
		*@brief Default constructor
		*/
		ExplicitOdeSolver_Parameters();

		/**
		*@brief Reads the options from a file
		*@param dictionary dictionary from which the oprions are extracted
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		template <typename OdeSolver>
		void TransferDataFromOdeSolver(const OdeSolver& dae_solver, const double cpuTime);

		/**
		*@brief Prints the current status on a file
		*@param fOut file where to print
		*/
		void Status(std::ostream& fOut);

		/**
		*@brief Sets the type of integrator
		*@param type integrator type: ODE_INTEGRATOR_OPENSMOKEPP, ODE_INTEGRATOR_BZZODE, ODE_INTEGRATOR_CVODE, ODE_INTEGRATOR_DASPK
		*/
		void SetType(const ODE_INTEGRATOR type) { type_ = type;}

		/**
		*@brief Sets the relative tolerance
		*@param relative_tolerance the relative tolerance (same for all the variables)
		*/
		void SetRelativeTolerance(const double relative_tolerance) { relative_tolerance_ = relative_tolerance; }

		/**
		*@brief Sets the absolute tolerance
		*@param absolute_tolerance the absolute tolerance (same for all the variables)
		*/
		void SetAbsoluteTolerance(const double absolute_tolerance) { absolute_tolerance_ = absolute_tolerance; }

		/**
		*@brief Sets the minimum step which can be used by the integrator
		*@param minimum_step the minimum step which can be used by the integrator
		*/
		void SetMinimumStep(const double minimum_step) { minimum_step_ = minimum_step; }

		/**
		*@brief Sets the maximum step which can be used by the integrator
		*@param maximum_step the maximum step which can be used by the integrator
		*/
		void SetMaximumStep(const double maximum_step) { maximum_step_ = maximum_step; }

		/**
		*@brief Sets the initial step which can has to be used by the integrator
		*@param intitial_step the initial step which can has to be used by the integrator
		*/
		void SetInitialStep(const double initial_step) { initial_step_ = initial_step; }

		/**
		*@brief Sets the maximum number of steps which can be used by the integrator
		*@param maximum_number_of_steps the maximum number of steps which can be used by the integrator
		*/
		void SetMaximumNumberOfSteps(const int maximum_number_of_steps) { maximum_number_of_steps_ = maximum_number_of_steps; }

		/**
		*@brief Sets the maximum order adopted by the integrator
		*@param maximum_order the maximum order adopted by the integrator
		*/
		void SetMaximumOrder(const int maximum_order) { maximum_order_ = maximum_order; }

		/**
		*@brief Sets the minimum residuals for considering the system of equations in steady state conditions
		*@param minimum_ypthe minimum residuals for considering the system of equations in steady state conditions
		*/
		void SetMinimumMeanThreshold(const double minimum_yp) { minimum_yp_ = minimum_yp; }

		/**
		*@brief Sets the verbosity level
		*@param verbosity_level the verbosity level (0 means no output)
		*/
		void SetVerbosityLevel(const int verbosity_level) { verbosity_level_ = verbosity_level; }
		
		/**
		*@brief Returns the integrator type: ODE_INTEGRATOR_OPENSMOKEPP, ODE_INTEGRATOR_BZZODE, ODE_INTEGRATOR_IDA, ODE_INTEGRATOR_DASPK
		*/
		ODE_INTEGRATOR type() const { return type_; }

		/**
		*@brief Returns the relative tolerance
		*/
		double relative_tolerance() const { return relative_tolerance_; }

		/**
		*@brief Returns the absolute tolerance
		*/
		double absolute_tolerance() const { return absolute_tolerance_; }

		/**
		*@brief Returns the minimum step (default or defined by the user)
		*/
		double minimum_step() const { return minimum_step_; }

		/**
		*@brief Returns the maximum step (default or defined by the user)
		*/
		double maximum_step() const { return maximum_step_; }

		/**
		*@brief Returns the initial step (default or defined by the user)
		*/
		double initial_step() const { return initial_step_; }

		/**
		*@brief Returns the maximum order (default or defined by the user)
		*/
		int maximum_order() const { return maximum_order_; }

		/**
		*@brief Returns the maximum number of steps (default or defined by the user)
		*/
		int maximum_number_of_steps() const { return maximum_number_of_steps_; }

		/**
		*@brief Returns the minimum residuals for considering the system of equations in steady state conditions (default or defined by the user)
		*/
		double minimum_yp() const { return minimum_yp_; }

		/**
		*@brief Returns the verbosity level
		*/
		int verbosity_level() const { return verbosity_level_; }

	private:

		ODE_INTEGRATOR type_;									//!< ODE integrator type
		double relative_tolerance_;								//!< relative tolerance
		double absolute_tolerance_;								//!< absolute tolerance
		double minimum_step_;									//!< minimum step to be adopted during the integration
		double maximum_step_;									//!< maximum step to be adopted during the integration
		double initial_step_;									//!< initial step to be adopted 
		int maximum_order_;										//!< maximum order
		int maximum_number_of_steps_;							//!< maximum number of steps

		double minimum_yp_;										//!< minimum residuals for considering the system of equations in steady state conditions
		int verbosity_level_;									//!< verbosity level
	};
}

#include "ExplicitOdeSolver_Parameters.hpp"

#endif	/* OpenSMOKE_ExplicitOdeSolver_Parameters_H */


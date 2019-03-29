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

#include "Grammar_ExplicitOdeSolver_Parameters.h"

namespace OdeSMOKE
{
	ExplicitOdeSolver_Parameters::ExplicitOdeSolver_Parameters()
	{
		// Default values
		type_ = ODE_RUNGEKUTTA45;
		
		// Default values are adopted for tolerances
		relative_tolerance_ = 100.*OpenSMOKE::MachEpsFloat();
		absolute_tolerance_ = 1.e-10;
		
		// Default values are adopted for initial and maximu steps
		minimum_step_ = -1.;
		maximum_step_ = -1.;
		initial_step_ = -1.;

		// Numerical details
		maximum_number_of_steps_ = 5000000;
		maximum_order_ = 5;
		minimum_yp_ = 1.e-7;
		verbosity_level_ = 1;
	}

	void ExplicitOdeSolver_Parameters::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_ExplicitOdeSolver_Parameters grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@OdeSolver") == true)
		{
			std::string name;
			dictionary.ReadString("@OdeSolver", name);
			
			if (name == "RungeKutta45")
			{
				type_ = ODE_RUNGEKUTTA45;
			}

			else OpenSMOKE::FatalErrorMessage("Unknown ODE Solver: " + name);
		}

		if (type_ == ODE_RUNGEKUTTA45)
		{
		}
		

		if (dictionary.CheckOption("@RelativeTolerance") == true)
			dictionary.ReadDouble("@RelativeTolerance", relative_tolerance_);

		if (dictionary.CheckOption("@AbsoluteTolerance") == true)
			dictionary.ReadDouble("@AbsoluteTolerance", absolute_tolerance_);

		if (dictionary.CheckOption("@MaximumOrder") == true)
			dictionary.ReadInt("@MaximumOrder", maximum_order_);

		if (dictionary.CheckOption("@MaximumStep") == true)
			dictionary.ReadDouble("@MaximumStep", maximum_step_);

		if (dictionary.CheckOption("@MinimumStep") == true)
			dictionary.ReadDouble("@MinimumStep", minimum_step_);

		if (dictionary.CheckOption("@InitialStep") == true)
			dictionary.ReadDouble("@InitialStep", initial_step_);

		if (dictionary.CheckOption("@MaximumNumberOfSteps") == true)
			dictionary.ReadInt("@MaximumNumberOfSteps", maximum_number_of_steps_);

		if (dictionary.CheckOption("@MeanResidualThreshold") == true)
			dictionary.ReadDouble("@MeanResidualThreshold", minimum_yp_);

		if (dictionary.CheckOption("@VerbosityLevel") == true)
			dictionary.ReadInt("@VerbosityLevel", verbosity_level_);
	}
}

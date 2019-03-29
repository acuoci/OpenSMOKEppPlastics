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

ThermogravimetricAnalysis::ThermogravimetricAnalysis()
{
	tEnd_ = 1.;
	initial_mass_ = 1.;
	is_verbose_ = true;
	options_ = new OpenSMOKE::ThermogravimetricAnalysis_Options();
	ode_options_ = new OdeSMOKE::ExplicitOdeSolver_Parameters();
}

void ThermogravimetricAnalysis::operator() (const std::string input_file_name, const std::string dictionary_name)
{
	Eigen::VectorXd profile_t, profile_T;

	OpenSMOKE::ThermogravimetricAnalysisFromDictionary(input_file_name, dictionary_name, P_, tEnd_, profile_t, profile_T, polymer_type_, *options_, *ode_options_);
	temperature_profile_(profile_t.size(), profile_t.data(), profile_T.data());
}

int ThermogravimetricAnalysis::SearchForLC(const double T)
{
	for (int i = 0; i < list_boiling_temperature_.size(); i++)
		if (list_boiling_temperature_[i] > T)
			return (i + 1);

	const int n = list_boiling_temperature_.size() - 1;
	std::stringstream label; label << list_boiling_temperature_[n];
	FatalErrorMessage("Fatal error: Maximum allowed temperature is " + label.str());

	return 0;
}


void ThermogravimetricAnalysis::FatalErrorMessage(const std::string message)
{
	std::cout << "Press enter to exit...";
	getchar();
	exit(-1);
}

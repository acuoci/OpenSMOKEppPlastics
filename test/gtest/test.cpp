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

// Include GoogleTest C++
#include "pch.h"
#include "include/gtest_cout.h"

// Include standard C++ libraries
#include <fstream>

// Include Eigen C++
#include <Eigen/Dense>

// Include utilities
#include "include/gtest_utilities.h"

// Folder containing data
std::string test_folder = "../../../../../data/";

// Include tests
#include "test_polystyrene/baseline.h"

int main(int argc, char *argv[])
{
	::testing::InitGoogleTest(&argc, argv);

	if (argc >= 2)
	{
		test_folder = argv[1];
		const std::string message = "Folder containing data: " + test_folder + "\n";
		PRINTF(message.c_str());
	}

	return RUN_ALL_TESTS();
}
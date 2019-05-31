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

#include "PolyethyleneKinetics.h"

namespace opensmokepp::plastics
{

	//!  A class for numerical modeling of thermal degradation of plastics
	/*!
		 A class for numerical modeling of thermal degradation of plastics
	*/

	class PolypropyleneKinetics : public PolyethyleneKinetics
	{
	public:

		/// @brief Default constructor
		PolypropyleneKinetics();

		/// @brief Default destructor
		~PolypropyleneKinetics();

	public:

		/// @brief Returns the density of the liquid phase (in kg/m3)
		/// @details TODO
		/// @param[in] T temperature (in K)
		/// @return the liquid density (in kg/m3)
		double LiquidDensity(const double T) const;

		/// @brief Returns the boiling temperature (in K), according to the Clausis-Clapeyron and Trouton-Meissner equations 
		/// @details TODO
		/// @param[in] L lower limit of number of monomeric units corresponding to liquid-phase species
		/// @param[in] P pressure (in atm)
		/// @return the boiling temperature (in K)
		double BoilingTemperature(const double L, const double P) const;

		/// @brief Returns the lower limit of number of monomeric units corresponding to liquid-phase species
		/// @details TODO
		/// @param[in] T temperature (in K)
		/// @param[in] P pressure (in Pa)
		/// @return the lower limit of number of monomeric units corresponding to liquid-phase species
		int MinNumberOfUnits(const double T, const double P);

		/// @brief Returns the current lower limit of number of monomeric units corresponding to liquid-phase species
		/// @details TODO
		/// @return the current lower limit of number of monomeric units corresponding to liquid-phase species
		inline int MinNumberOfUnits() { return LC_; }

	private:

		/// @brief Default kinetic constants (from polystirene)
		void DefaultKineticConstants();

	};

}


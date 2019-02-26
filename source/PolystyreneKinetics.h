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

	//!  A class for numerical modeling of thermal degradation of plastics
	/*!
		 A class for numerical modeling of thermal degradation of plastics
	*/

	class PolystyreneKinetics
	{
	public:

		/// @brief Default constructor
		PolystyreneKinetics();

		/// @brief Default destructor
		~PolystyreneKinetics();

		/// @brief Returns the density of the liquid phase (in kg/m3)
		/// @details TODO
		/// @param[in] T temperature (in K)
		/// @return the liquid density (in kg/m3)
		double LiquidDensity(const double T);

		/// @brief Returns the boiling temperature (in K), according to the Clausis-Clapeyron and Trouton-Meissner equations 
		/// @details TODO
		/// @param[in] L lower limit of number of monomeric units corresponding to liquid-phase species
		/// @param[in] P pressure (in atm)
		/// @return the boiling temperature (in K)
		double BoilingTemperature(const double L, const double P);

		/// @brief Sets the maximum number of monomeric units to be tracked
		/// @details TODO
		/// @param[in] N maximum number of monomeric units
		void SetMaxNumberOfUnits(const int N);

		/// @brief Sets the minimum number of monomeric units in the liquid phase
		/// @details TODO
		/// @param[in] kd minimum number of monomeric units
		void SetMinNumberOfUnits(const int kd);

		double SplittingCoefficient(const double T, const double P);

		void CalculateSharedSpecies(const double T, const double P);

		void UpdateSharedSpeciesDistribution(const double T, const double P, Eigen::VectorXd& y);

		/// @brief Returns the density of gaseous phase (in kg/m3)
		/// @details TODO
		/// @return the density of gaseous phase (in kg/m3)
		double Gas();

		/// @brief Returns the density of liquid phase (in kg/m3)
		/// @details TODO
		/// @return the density of liquid phase (in kg/m3)
		double Residual();

		/// @brief Calculates the concentration (in mol/l or kmol/m3) of radicals which can undergo 
		///        random scission (CCsr) or allyl scission (CCsa)
		/// @details TODO
		void CCBonds();

		/// @brief Sets the initial concentration (in mol/l or kmol/m3) of radicals which can undergo 
		///        random scission (CCsr) or allyl scission (CCsa)
		/// @details TODO
		/// @param[in] CCsr concentration of radicals which can undergo random scission (in mol/l or kmol/m3)
		/// @param[in] CCar concentration of radicals which can undergo allyl scission (in mol/l or kmol/m3)
		void SetCCBonds(const double CCsr, const double CCsa);

		/// @brief Sets the current status for internal calculations
		/// @details TODO
		/// @param[in] T temperature (in K)
		/// @param[in] P pressure (in atm)
		/// @param[in] Ctot total concentration (in mol/l or kmol/m3)
		/// @param[in] x mole fractions
		void SetStatus(const double T, const double P, const double Ctot, const Eigen::VectorXd& x);
		
		/// @brief Sets the current status for internal calculations
		/// @details TODO
		/// @param[in] T temperature (in K)
		/// @param[in] P pressure (in atm)
		/// @param[in] c concentrations (in mol/l or kmol/m3)
		void SetStatus(const double T, const double P, const Eigen::VectorXd& c);

		/// @brief Updates the kinetic constants (in mol,l,s) according to the internal status
		/// @details TODO
		void KineticConstants();

		/// @brief Updates the formation rates according to the internal status (in mol/l/s or kmol/m3/s)
		/// @details TODO
		void FormationRates();

		void UpdateGasDistribution();

		/// @brief Returns the current formation rates (in mol/l/s or kmol/m3/s)
		/// @return formation rates (in mol/l/s or kmol/m3/s)
		const Eigen::VectorXd& R() const { return R_; }

		void SetVerbose(const bool flag) { is_verbose_ = flag; }

		/// @brief Returns the total number of moles (in mol) of paraffins in gaseous phase
		/// @return total number of moles (in mol) of paraffins in gaseous phase
		double ypar() const { return ypar_; }

		/// @brief Returns the total number of moles (in mol) of olefins in gaseous phase
		/// @return total number of moles (in mol) of olefins in gaseous phase
		double yole() const { return yole_; }

		/// @brief Returns the total number of moles (in mol) of olefins in gaseous phase
		/// @return total number of moles (in mol) of diolefins in gaseous phase
		double ydio() const { return ydio_; }

		/// @brief Returns the total number of moles (in mol) in gaseous phase
		/// @return total number of moles (in mol) in gaseous phase
		double ytot() const { return ytot_; }

		double wpar() const { return wpar_; }
		double wole() const { return wole_; }
		double wdio() const { return wdio_; }
		double wtot() const { return wtot_; }

		double CLiquid();

		double SumGas(const Eigen::VectorXd& y);

		double SumLiquid(const Eigen::VectorXd& y);

		double SumGasMW(const Eigen::VectorXd& y);

		double SumLiquidMW(const Eigen::VectorXd& y);

	private:

		/// @brief Default kinetic constants (from polystirene)
		void DefaultKineticConstants();

		/// @brief Calculates the formation rates of species in liquid phase
		///        in case lumping is disabled
		void LiquidFormationRatesWithoutLumping();

		/// @brief Calculates the formation rates of species in liquid phase 
		///        in case lumping is enabled
		void LiquidFormationRatesWithLumping();

		/// @brief Calculates the formation rates of species in gaseous phase
		void GaseousFormationRates();

		/// @brief Returns the fraction of secondary benzyl radicals Rsb abstracting tertiary Hs
		/// @param[in] C concentration of tertiary Hs on which abstraction can occur (in mol/l or kmol/m3)
		/// @return fraction of secondary benzyl radicals Rsb abstracting tertiary Hs
		inline double Beta(const double C);

		/// @brief Returns the fraction of secondary benzyl radicals Rsb abstracting tertiary Hs
		/// @param[in] C concentration of tertiary Hs on which abstraction can occur (in mol/l or kmol/m3)
		/// @return fraction of secondary benzyl radicals Rsb abstracting tertiary Hs
		inline double BetaWithoutBackBiting(const double C);

		/// @brief Returns the fraction of secondary benzyl radicals Rsb abstracting tertiary Hs
		/// @param[in] C concentration of tertiary Hs on which abstraction can occur (in mol/l or kmol/m3)
		/// @return fraction of secondary benzyl radicals Rsb abstracting tertiary Hs
		inline double BetaWithLumping(const double C);

		/// @brief Returns the fraction of secondary benzyl radicals Rsb following back-biting
		/// @param[in] C concentration of tertiary Hs on which abstraction can occur (in mol/l or kmol/m3)
		/// @return fraction of secondary benzyl radicals Rsb following back-biting
		inline double Gamma(const double C);

		/// @brief Returns the fraction of secondary benzyl radicals Rsb following back-biting
		/// @param[in] C concentration of tertiary Hs on which abstraction can occur (in mol/l or kmol/m3)
		/// @return fraction of secondary benzyl radicals Rsb following back-biting
		inline double GammaWithoutBackBiting(const double C);

		/// @brief Returns the fraction of secondary benzyl radicals Rsb following back-biting
		/// @param[in] C concentration of tertiary Hs on which abstraction can occur (in mol/l or kmol/m3)
		/// @return fraction of secondary benzyl radicals Rsb following back-biting
		inline double GammaWithLumping(const double C);


	private:

		/// @section kinetic_parameters Kinetic parameters
		/// @brief   Kinetic parameters of polystyrene reactions
		/// @details Frequency factors are always in: mol, l, s
		///          Activation energies are always in: cal/mol

		double As_;		/*!< Random scission frequency factor: PS -> Rp + Rsb */
		double Es_;		/*!< Random scission activation energy: PS -> Rp + Rsb */

		double Aa_;		/*!< Allylic scission frequency factor: PS -> Ra + Rsb */
		double Ea_;		/*!< Allylic scission activation energy: PS -> Ra + Rsb */

		double At_;		/*!< Termination frequency factor:	Rsb + Rp -> PS  and  Rp + Rp -> PS */
		double Et_;		/*!< Termination activation energy:	Rsb + Rp -> PS  and  Rp + Rp -> PS */

		double Aar_;	/*!< H-abstractions frequency factor: Rp + PS -> PS + Rt */
		double Ear_;	/*!< H-abstractions activation energy: Rp + PS -> PS + Rt */

		double Aar1_;	/*!< H-abstractions frequency factor: Rp + PS -> PS + Rt */
		double Ear1_;	/*!< H-abstractions activation energy: Rp + PS -> PS + Rt */

		double Aar2_;	/*!< H-abstractions frequency factor: Rsb + PS -> PS + Rt */
		double Ear2_;	/*!< H-abstractions activation energy: Rsb + PS -> PS + Rt */

		double Au_;		/*!< Unzipping frequency factor: Rsb -> Rsb + S */
		double Eu_;		/*!< Unzipping activation energy: Rsb -> Rsb + S */

		double Abb_;	/*!< Back-biting (intramolecular abstractions) frequency factor: Rsb -> Rt  and  Rp -> Rt */
		double Ebb_;	/*!< Back-biting (intramolecular abstractions) activation energy: Rsb -> Rt  and  Rp -> Rt */

		double Abeta_;	/*!< Beta-decomposition frequency factor: Rt -> PS + Rsb */
		double Ebeta_;	/*!< Beta-decomposition activation energy: Rt -> PS + Rsb */

		double Aar3_;	/*!< H-abstractions frequency factor: Rt + PS -> PS + Rt' */
		double Ear3_;	/*!< H-abstractions activation energy: Rt + PS -> PS + Rt' */

		double fe1_;	/*!< Weight for effective kinetic constant of abstraction reactions */
		double fe2_;	/*!< Weight for effective kinetic constant of abstraction reactions */
		double fe3_;	/*!< Weight for effective kinetic constant of abstraction reactions */

		double kef_;	/*!< inter-molecular abstraction (effective) */
		double kp_;		/*!< propagation (apparent) */
		double ku_;		/*!< unzipping */
		double kb_;		/*!< equivalent kinetic constant considering back-biting as abstraction mechanism */

		Eigen::VectorXd radbbp1_;		/*!< TODO */
		Eigen::VectorXd radbbp3_;		/*!< TODO */
		Eigen::VectorXd radbbo1_;		/*!< TODO */
		Eigen::VectorXd polbbp1_;		/*!< TODO */
		Eigen::VectorXd polbbp3_;		/*!< TODO */
		Eigen::VectorXd polbbo1_;		/*!< TODO */

		bool is_delta_enabled_;			/*!< if true, correction for beta-scissions on the right side is enabled */
		double delta_;					/*!< correction factor for beta-scissions on the right side */
		

		double MWm_;					/*!< molecular weight of monomer (in g/mol or kg/kmol) */

		
		bool is_backbiting_enabled_;	/*!< if true, back-biting reactions are enabled */
		

		double eff_;					/*!< effective fraction of radicals from random scission */

		int N_;							/*!< maximum number of monomeric units to be tracked */
		int kd_;						/*!< number of monomeric units in the shared (gas/liquid) species */

		bool is_lumping_enabled_;		/*!< if true, lumping is enabled */
		int lumping_start_;				/*!< number of monomeric units at which lumping starts */
		int lumping_step_;				/*!< number of monomeric units defining the lumping interval */

		double sum3_;					/*!< Sum(j=n+1,Inf)PIj */
		double sum1_;					/*!< Sum(j=n+2,Inf)PIj */	
		double sum6_;					/*!< Sum(j=n+1,Inf)PIj equivalent to sum3_ */
		double sum11_;					/*!< Sum(j=n+1,Inf)PIj equivalent to sum3_ */

		double sum7_;					/*!< Sum(j=n+1,Inf)PIIIj */
		double sum4_;					/*!< Sum(j=n+2,Inf)PIIIj */

		double sum13_;					/*!< Sum(j=n+1,Inf)OIj */
		double sum5_;					/*!< Sum(j=n+2,Inf)OIj */
		double sum8_;					/*!< Sum(j=n+2,Inf)OIj equivalent to sum5_ */
		
		double sum9_;					/*!< Sum(j=n+1,Inf)OIIj */
		double sum2_;					/*!< Sum(j=n+2,Inf)OIIj equivalent to sum2_ */
		double sum12_;					/*!< Sum(j=n+2,Inf)OIIj */
		double sum14_;					/*!< Sum(j=n+1,Inf)OIIj equivalent to sum9_ */

		double sum10_;					/*!< Sum(j=n+2,Inf)DIIj */
		double sum15_;					/*!< Sum(j=n+2,Inf)DIIj equivalent to sum10_ */

		double sum16_;					/*!< Additional term associated to PI family */
		double sum17_;					/*!< Additional term associated to PIII family */
		double sum18_;					/*!< Additional term associated to OI family */
		double sum19_;					/*!< Additional term associated to OII family */
		double sum20_;					/*!< Additional term associated to DII family */

		double radp1_;					/*!< (1-Beta-Gamma)*RPI */
		double radp3_;					/*!< (1-Beta-Gamma)*RPIII */
		double rado1_;					/*!< (1-Beta-Gamma)*ROI */

		double sum_monomer_;			/*!< formation rate of TODO (in TODO) */
		double sum_trimer_;				/*!< formation rate of trimer OI3 (in TODO) */
		double sum_13diphenylpropyl_;	/*!< formation rate of 1,3-diphenylpropyl PIII2 (in TODO) */	

		Eigen::VectorXd R_;	/*!< formation rates (in mol/l/s or kmol/m3/s) */

		double CCsr_;		/*!< concentration (in mol/l or kmol/m3) of radicals which can undergo random scission */
		double CCsa_;		/*!< concentration (in mol/l or kmol/m3) of radicals which can undergo allyl scission */

		double ysompar1_;	/*!< total number of moles (in mol) of paraffins (type I) in gaseous phase */
		double ysompar3_;	/*!< total number of moles (in mol) of paraffins (type III) in gaseous phase */
		double ysomole1_;	/*!< total number of moles (in mol) of olefins (type I) in gaseous phase */
		double ysomole2_;	/*!< total number of moles (in mol) of olefins (type II) in gaseous phase */
		double ysomdio2_;	/*!< total number of moles (in mol) of diolefins (type II) in gaseous phase */
		double ypar_;		/*!< total number of moles (in mol) of paraffins in gaseous phase */
		double yole_;		/*!< total number of moles (in mol) of olefins in gaseous phase */
		double ydio_;		/*!< total number of moles (in mol) of diolefins in gaseous phase */
		double ytot_;		/*!< total number of moles (in mol) in gaseous phase */

		double wpar1_;		/*!< total mass (in g) of paraffins (type I) in gaseous phase */
		double wpar3_;		/*!< total mass (in g) of paraffins (type III) in gaseous phase */
		double wole1_;		/*!< total mass (in g) of olefins (type I) in gaseous phase */
		double wole2_;		/*!< total mass (in g) of olefins (type II) in gaseous phase */
		double wdio2_;		/*!< total mass (in g) of diolefins (type II) in gaseous phase */
		double wpar_;		/*!< total mass (in g) of paraffins in gaseous phase */
		double wole_;		/*!< total mass (in g) of olefins in gaseous phase */
		double wdio_;		/*!< total mass (in g) of diolefins in gaseous phase */
		double wtot_;		/*!< total mass (in g) in gaseous phase */

		Eigen::VectorXd y;		/*!< concentrations (in mol/l or kmol/m3) */
		double T_;				/*!< temperature (in K) */
		double P_;				/*!< pressure (in atm) */
		double Ctot_;			/*!< total concentration (in mol/l or kmol/m3) */

		bool is_verbose_;		/*!< if true, additional data are written on files */

	private:

		static double Rgas_;	/*!< Ideal gas universal constant R (in cal/mol/K) */
	};

}


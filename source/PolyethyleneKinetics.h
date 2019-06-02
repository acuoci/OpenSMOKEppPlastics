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

	class PolyethyleneKinetics
	{
	public:

		/// @brief Default constructor
		PolyethyleneKinetics();

		/// @brief Default destructor
		~PolyethyleneKinetics();

	public:

		/// @section setting_functions Setting functions
		/// @brief   Functions to set the internal parameters/options

		/// @brief Sets the maximum number of monomeric units to be tracked
		/// @details TODO
		/// @param[in] N maximum number of monomeric units
		void SetMaxNumberOfUnits(const int N);

		/// @brief Sets the minimum number of monomeric units in the liquid phase
		/// @details TODO
		/// @param[in] kd minimum number of monomeric units
		void SetMinNumberOfUnits(const int kd);

		/// @brief Sets the molecular weight of the polymer
		/// @details TODO
		/// @param[in] MW polymer molecular weight (in g/mol or kg/kmol)
		void SetMolecularWeightPolymer(const double MW) { MWp_ = MW; }

		/// @brief Sets the molecular weight of the monomer
		/// @details TODO
		/// @param[in] MW monomer molecular weight (in g/mol or kg/kmol)
		void SetMolecularWeightMonomer(const double MW) { MWm_ = MW; }


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


		/// @brief Sets the frequency factor for random-scission reactions
		/// @param[in] As the frequency factor (units in mol,l,s)
		void SetAs(const double As) { A_(0) = As; }

		/// @brief Sets the frequency factor for allylic-scission reactions
		/// @param[in] Aa the frequency factor (units in mol,l,s)
		void SetAa(const double Aa) { A_(1) = Aa; }

		/// @brief Sets the frequency factor for termination reactions
		/// @param[in] At the frequency factor (units in mol,l,s)
		void SetAt(const double At) { A_(2) = At; }

		/// @brief Sets the frequency factor for H-abstraction reactions
		/// @param[in] Ar the frequency factor (units in mol,l,s)
		void SetAaf(const double Aaf) { A_(3) = Aaf; }

		/// @brief Sets the frequency factor for H-reabstraction reactions
		/// @param[in] Ar1 the frequency factor (units in mol,l,s)
		void SetAab(const double Aab) { A_(4) = Aab; }

		/// @brief Sets the frequency factor for beta-scission reactions
		/// @param[in] Abeta the frequency factor (units in mol,l,s)
		void SetAbeta(const double Abeta) { A_(5) = Abeta; }

		/// @brief Sets the frequency factor for back-biting reactions (isomerization 1-4)
		/// @param[in] Abb14 the frequency factor (units in mol,l,s)
		void SetAbb14(const double Abb14) { A_(6) = Abb14; }

		/// @brief Sets the frequency factor for back-biting reactions (isomerization 1-5)
		/// @param[in] Abb15 the frequency factor (units in mol,l,s)
		void SetAbb15(const double Abb15) { A_(7) = Abb15; }

		/// @brief Sets the frequency factor for back-biting reactions (isomerization 1-6)
		/// @param[in] Abb16 the frequency factor (units in mol,l,s)
		void SetAbb16(const double Abb16) { A_(8) = Abb16; }


		/// @brief Sets the temperature exponent for termination reactions
		/// @param[in] Betat the temperature exponent
		void SetBetat(const double Betat) { Beta_(2) = Betat; }


		/// @brief Sets the activation energy for random-scission reactions
		/// @param[in] Es the activation energy (in cal/mol)
		void SetEs(const double Es) { E_(0) = Es; }

		/// @brief Sets the activation energy for allylic-scission reactions
		/// @param[in] Ea the activation energy (in cal/mol)
		void SetEa(const double Ea) { E_(1) = Ea; }

		/// @brief Sets the activation energy for termination reactions
		/// @param[in] Eat the activation energy (in cal/mol)
		void SetEt(const double Et) { E_(2) = Et; }

		/// @brief Sets the activation energy for H-abstraction reactions
		/// @param[in] Eaf the activation energy (in cal/mol)
		void SetEaf(const double Eaf) { E_(3) = Eaf; }

		/// @brief Sets the activation energy for H-reabstraction reactions
		/// @param[in] Eab the activation energy (in cal/mol)
		void SetEab(const double Eab) { E_(4) = Eab; }

		/// @brief Sets the activation energy for beta-scission reactions
		/// @param[in] Ebeta the activation energy (in cal/mol)
		void SetEbeta(const double Ebeta) { E_(5) = Ebeta; }

		/// @brief Sets the activation energy for back-biting reactions (isomerization 1-4)
		/// @param[in] Ebb14 the activation energy (in cal/mol)
		void SetEbb14(const double Ebb14) { E_(6) = Ebb14; }

		/// @brief Sets the activation energy for back-biting reactions (isomerization 1-5)
		/// @param[in] Ebb15 the activation energy (in cal/mol)
		void SetEbb15(const double Ebb15) { E_(7) = Ebb15; }

		/// @brief Sets the activation energy for back-biting reactions (isomerization 1-6)
		/// @param[in] Ebb16 the activation energy (in cal/mol)
		void SetEbb16(const double Ebb16) { E_(8) = Ebb16; }

		/// @brief Sets acceleration factor for butane
		/// @param[in] Cbut acceleration factor for butane
		void SetCbutane(const double Cbut) { Cbut_ = Cbut; }

		/// @brief Sets acceleration factor for propane
		/// @param[in] Cprop acceleration factor for propane
		void SetCpropane(const double Cprop) { Cprop_ = Cprop; }


		/// @brief Sets the verbosity level
		/// @param[in] flag verbosity level
		void SetVerbose(const bool flag) { is_verbose_ = flag; }

		/// @brief Sets the initial acceleration factor
		/// @param[in] if true, initial acceleration factor is enabled
		/// @param[in] initial weight (in g)
		void SetInitialAccelerationFactor(const bool flag);


	public:

		/// @section access_functions Access functions
		/// @brief   Functions to access the internal parameters/options

		/// @brief Returns the current formation rates (in mol/l/s or kmol/m3/s)
		/// @return formation rates (in mol/l/s or kmol/m3/s)
		inline const Eigen::VectorXd& R() const { return R_; }
		

		/// @brief Returns the frequency factor for random-scission reactions
		/// @return the frequency factor (units in mol,l,s)
		inline double As() const { return A_(0); }

		/// @brief Returns the frequency factor for allylic-scission reactions
		/// @return the frequency factor (units in mol,l,s)
		inline double Aa() const { return A_(1); }

		/// @brief Returns the frequency factor for termination reactions
		/// @return the frequency factor (units in mol,l,s)
		inline double At() const { return A_(2); }

		/// @brief Returns the frequency factor for H-abstraction reactions
		/// @return the frequency factor (units in mol,l,s)
		inline double Aaf() const { return A_(3); }

		/// @brief Returns the frequency factor for H-reabstraction reactions
		/// @return the frequency factor (units in mol,l,s)
		inline double Aab() const { return A_(4); }

		/// @brief Returns the frequency factor for beta-scission reactions
		/// @return the frequency factor (units in mol,l,s)
		inline double Abeta() const { return A_(5); }

		/// @brief Returns the frequency factor for back-biting reactions (isomerization 1-4)
		/// @return the frequency factor (units in mol,l,s)
		inline double Abb14() const { return A_(6); }

		/// @brief Returns the frequency factor for back-biting reactions (isomerization 1-5)
		/// @return the frequency factor (units in mol,l,s)
		inline double Abb15() const { return A_(7); }

		/// @brief Returns the frequency factor for back-biting reactions (isomerization 1-6)
		/// @return the frequency factor (units in mol,l,s)
		inline double Abb16() const { return A_(8); }


		/// @brief Returns the activation energy for random-scission reactions
		/// @return the activation energy (in cal/mol)
		inline double Es() const { return E_(0); }

		/// @brief Returns the activation energy for allylic-scission reactions
		/// @return the activation energy (in cal/mol)
		inline double Ea() const { return E_(1); }

		/// @brief Returns the activation energy for termination reactions
		/// @return the activation energy (in cal/mol)
		inline double Et() const { return E_(2); }

		/// @brief Returns the activation energy for H-abstraction reactions
		/// @return the activation energy (in cal/mol)
		inline double Eaf() const { return E_(3); }

		/// @brief Returns the activation energy for H-reabstraction reactions
		/// @return the activation energy (in cal/mol)
		inline double Eab() const { return E_(4); }

		/// @brief Returns the activation energy for beta-scission reactions
		/// @return the activation energy (in cal/mol)
		inline double Ebeta() const { return E_(5); }

		/// @brief Returns the activation energy for back-biting reactions (isomerization 1-4)
		/// @return the activation energy (in cal/mol)
		inline double Ebb14() const { return E_(6); }

		/// @brief Returns the activation energy for back-biting reactions (isomerization 1-5)
		/// @return the activation energy (in cal/mol)
		inline double Ebb15() const { return E_(7); }

		/// @brief Returns the activation energy for back-biting reactions (isomerization 1-6)
		/// @return the activation energy (in cal/mol)
		inline double Ebb16() const { return E_(8); }


		/// @brief Returns the maximum number of monomeric units to be tracked
		/// @details TODO
		/// @return the maximum number of monomeric units to be tracked
		inline int MaxNumberOfUnits() const { return NPA_; }

		/// @brief Returns the polymer molecular weight (in g/mol or kg/kmol)
		/// @return the polymer molecular weight (in g/mol or kg/kmol)
		inline double MolecularWeightPolymer() const { return MWp_; }

		/// @brief Returns the monomer molecular weight (in g/mol or kg/kmol)
		/// @return the monomer molecular weight (in g/mol or kg/kmol)
		inline double MolecularWeightMonomer() const { return MWm_; }


	public:

		/// @section todo todo
		/// @brief   todo todo

		/// @brief Returns the sum of elements of input vector belonging to the gaseous phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @return the requested sum
		double SumGas(const Eigen::VectorXd& v) const;

		/// @brief Returns the sum of elements of input vector belonging to the gaseous phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @param[out] P the sum contribution associated to paraffins
		/// @param[out] O the sum contribution associated to olefins
		/// @param[out] D the sum contribution associated to diolefins
		void SumGas(const Eigen::VectorXd& v, double& P, double& O, double& D) const;

		/// @brief Returns the sum of elements of input vector belonging to the liquid phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @return the requested sum
		double SumLiquid(const Eigen::VectorXd& v) const;

		/// @brief Returns the sum of elements of input vector belonging to the liquid phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @param[out] P the sum contribution associated to paraffins
		/// @param[out] O the sum contribution associated to olefins
		/// @param[out] D the sum contribution associated to diolefins
		void SumLiquid(const Eigen::VectorXd& v, double& P, double& O, double& D) const;

		/// @brief Returns the sum (weighted on the molecular weights) of elements of input vector belonging to the gaseous phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @return the requested weighted sum
		double SumGasMW(const Eigen::VectorXd& v) const;

		/// @brief Returns the sum (weighted on the molecular weights) of elements of input vector belonging to the gaseous phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @param[out] P the sum contribution associated to paraffins
		/// @param[out] O the sum contribution associated to olefins
		/// @param[out] D the sum contribution associated to diolefins
		void SumGasMW(const Eigen::VectorXd& v, double& P, double& O, double& D) const;

		/// @brief Returns the sum (weighted on the molecular weights) of elements of input vector belonging to the liquid phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @return the requested weighted sum
		double SumLiquidMW(const Eigen::VectorXd& v) const;

		/// @brief Returns the sum (weighted on the molecular weights) of elements of input vector belonging to the liquid phase
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @param[out] P the sum contribution associated to paraffins
		/// @param[out] O the sum contribution associated to olefins
		/// @param[out] D the sum contribution associated to diolefins
		void SumLiquidMW(const Eigen::VectorXd& v, double& P, double& O, double& D) const;

		/// @brief Returns the sum (weighted on the molecular weights) of elements of input vector
		/// @details TODO
		/// @param[in] v vector of elements to be processed
		/// @return the weighted sum
		double Sum(const Eigen::VectorXd& v) const;
		double SumMW(const Eigen::VectorXd& v) const;

		void Sum(const Eigen::VectorXd& v, double& P, double& O, double& D) const;

		void SumMW(const Eigen::VectorXd& v, double& P, double& O, double& D) const;

		void Distribution(const Eigen::VectorXd& v, Eigen::VectorXd& P, Eigen::VectorXd& O, Eigen::VectorXd& D) const;
		void DistributionMW(const Eigen::VectorXd& v, Eigen::VectorXd& P, Eigen::VectorXd& O, Eigen::VectorXd& D) const;

		
		


	public:

		/// @section todo todo
		/// @brief   todo todo

		/// @brief Updates the initial acceleration coefficient
		/// @details TODO
		/// @param[in] VL liquid-phase volume (in l)
		/// @param[in] WG initial mass (in g)
		void UpdateInitialAccelerationCoefficient(const double VL, const double WG);

		/// @brief Updates the kinetic constants (in mol,l,s) according to the internal status
		/// @details TODO
		void KineticConstants();

		/// @brief Updates the formation rates according to the internal status (in mol/l/s or kmol/m3/s)
		/// @details TODO
		void FormationRates();


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

		/// @brief Calculates the formation rates of species in liquid phase
		void LiquidFormationRates();

		/// @brief Calculates the formation rates of species in gaseous phase
		void GaseousFormationRates();

		/// @brief Calculates the correction factors
		void UpdateCorrections();

		/// @brief Calculates the sums needed for building the formation rates
		void UpdateSums(const int J);


	protected:

		/// @section kinetic_parameters Kinetic parameters
		/// @brief   Kinetic parameters of polystyrene reactions
		/// @details Frequency factors are always in: mol, l, s
		///          Activation energies are always in: cal/mol

		double MWp_;		/*!< molecular weight of polymer (in g/mol or kg/kmol) */
		double MWm_;		/*!< molecular weight of monomer (in g/mol or kg/kmol) */

		int LC_;			/*!< number of monomeric units in the shared (gas/liquid) species */

		int NPA_;			/*!< total number of paraffins */
		int NOL_;			/*!< total number of olefins */
		int NDO_;			/*!< total number of diolefins */

		Eigen::VectorXd R_;	/*!< formation rates (in mol/l/s or kmol/m3/s) */


		Eigen::VectorXd c_;	/*!< concentrations (in mol/l or kmol/m3) */
		double T_;			/*!< temperature (in K) */
		double P_;			/*!< pressure (in atm) */

		bool is_verbose_;	/*!< if true, additional data are written on files */

		Eigen::VectorXd A_;			/*!< frequency factors (in mol,l,s or kmol,m3,s) */
		Eigen::VectorXd Beta_;		/*!< temperature exponent in kinetic constant */
		Eigen::VectorXd E_;			/*!< activation energy (in cal/mol) */
		Eigen::VectorXd k_;			/*!< kinetic constants (in mol,l,s or kmol,m3,s) */

		double FR_;
		double Cbut_;	/*!< Butane acceleration factor for possible resonances (Nava & Fabini, 1996) */
		double Cprop_;	/*!< Propane acceleration factors for possible resonances (Nava & Fabini, 1996) */

		bool is_initial_acceleration_factor_;		/*!< initial acceleration factor due to ramifications and impurities */
		double initial_acceleration_coefficient_;	/*!< initial acceleration coefficient */

		Eigen::VectorXd CorrO_;		/*!< corrections on olefins */
		Eigen::VectorXd CorrP_;		/*!< corrections on paraffins */

		Eigen::VectorXd sumY_;		/*!< internal vectors */
		Eigen::VectorXd sumY1_;		/*!< internal vectors */
		Eigen::VectorXd sumY12_;	/*!< internal vectors */
		Eigen::VectorXd sumY2_;		/*!< internal vectors */

		double ATOT_;				/*!< total sum of probabilities for back-biting */
		Eigen::VectorXd alpha_;		/*!< probabilities for back-biting */

		Eigen::VectorXi M_;			/*!< internal vectors */
		Eigen::VectorXi LW_;		/*!< internal vectors */
		Eigen::VectorXi N_;			/*!< internal vectors */

		double sum1_;				/*!< sum on paraffins */
		double sum2_;				/*!< sum on olefins */
		double sum3_;				/*!< sum on paraffins */
		double sum4_;				/*!< sum on olefins */
		double sum5_;				/*!< sum on diolefins */
		double sum6_;				/*!< sum on olefins */
		double sum7_;				/*!< sum on diolefins */

		
		double kp_;					/*!< Propagation constant kp=ki*[R] */
		double kd_;					/*!< Auxiliary kinetic constant */
		double ka_;					/*!< Auxiliary kinetic constant */

	protected:

		static double Rgas_;		/*!< Ideal gas universal constant R (in cal/mol/K) */
	};

}


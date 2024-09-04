// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__LOWEFIT__H
#define __BAT__LOWEFIT__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>

// This is a lowEfit header file.
// Model source code is located in file GammaFit/lowEfit.cxx

// ---------------------------------------------------------
class lowEfit : public BCModel {

    public:

        // Constructor
        lowEfit(const std::string& name);

        // Destructor
        ~lowEfit();

        // Overload LogLikelihood to implement model
        double LogLikelihood(const std::vector<double>& pars);

        // Overload LogAprioriProbability if not using built-in 1D priors
        // double LogAPrioriProbability(const std::vector<double> & pars);

        // Define the model function
        double FitFunction(const std::vector<double>& pars, double x, double sigma);
        // Provide the data (bin centers, observed counts, sigmas)
        void SetData(const std::map<std::string, std::vector<double>>& x, const std::map<std::string, std::vector<double>>& y,
                     const std::map<std::string, double>& resolutions);

        // Overload CalculateObservables if using observables
        // void CalculateObservables(const std::vector<double> & pars);

    private:

        double E_peak = 575;

        std::map<std::string, double> sigma_map;

        std::map<std::string, std::vector<double>> BinHeight_map;
        std::map<std::string, std::vector<double>> BinCenter_map;

};
// ---------------------------------------------------------

void readCSV (const std::string& filename,
              std::map<std::string, std::vector<double>>& bin_value_map, 
              std::map<std::string, std::vector<double>>& bin_center_map,
              std::map<std::string, double>& sigma_map);

#endif

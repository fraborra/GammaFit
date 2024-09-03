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
class lowEfit : public BCModel
{

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
    double FitFunction(const std::vector<double>& params, double x, double sigma);
    // Provide the data (bin centers, observed counts, sigmas)
    void SetData(const std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& y,
                 const std::vector<double>& resolutions);

    // Overload CalculateObservables if using observables
    // void CalculateObservables(const std::vector<double> & pars);

private:

    double E_peak = 575;

    unsigned int Nge = 1;
    std::vector<double> sigmas;

    std::vector<std::vector<double>> BinHeight {Nge};
    std::vector<std::vector<double>> BinCenter {Nge};

};
// ---------------------------------------------------------

#endif

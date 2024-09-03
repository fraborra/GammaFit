// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "lowEfit.h"
#include <BAT/BCMath.h>

// ---------------------------------------------------------
lowEfit::lowEfit(const std::string& name) : BCModel(name) {
    // Define parameters here in the constructor.
    // Also define their priors, if using built-in priors.
    // For example:
    // AddParameter("mu", -2, 1, "#mu", "[GeV]");
    // GetParameters.Back().SetPrior(new BCGaussianPrior(-1, 0.25));

    AddParameter("a", -10, 0);   // Quadratic coefficient
    AddParameter("b", -10, 10);   // Linear coefficient
    AddParameter("c", 0, 10);   // Constant term
    AddParameter("A", 0, 100);    // Amplitude of the Gaussian

    // Define observables here, too. For example:
    // AddObservable("mu_squared", 1, 4, "#mu^{2}", "[GeV^{2}]");
}

// ---------------------------------------------------------
lowEfit::~lowEfit() {
    // destructor
}

// ---------------------------------------------------------
double lowEfit::LogLikelihood(const std::vector<double>& pars) {
    // return the log of the conditional probability p(data|pars).
    // This is where you define your model.
    // BCMath contains many functions you will find helpful.

    // store our log-likelihood as we loop through bins
    double LL = 0.;
    // loop over detectors
    for (unsigned int j = 0; j <= Nge; ++j) {

        double sigma = sigmas[j];
        unsigned int NBins = BinHeight[j].size();

        // loop over bins of our data
        for (unsigned int i = 0; i <= NBins; ++i) {
            // retrieve observed number of events
            double x = BinHeight[j][i];
            // retrieve bin center
            double m = BinCenter[j][i];
            // calculate expected number of events, using ROOT Gaus function
            double nu = FitFunction(pars, m, sigma);
            // add to log-likelihood sum
            LL += BCMath::LogPoisson(x, nu);
        }
    } 
    // return log-likelihood
    return LL; 
}


// ---------------------------------------------------------
// double lowEfit::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     double LL = 0.;
//     double s = pars[0], b = pars[1], A = pars[2], B = pars[3], C = pars[4];
//     double bmax = h; // bin height
    
//     double s_max = 10;
    
//     if (s < 0 || s > s_max) {
//         LL += log(0.0);
//     } else {
//         LL += log(1.0/s_max);
//     }
//     if (b < 0 || b > bmax || A > 0 || C < 0) {
//         return -std::numeric_limits<double>::infinity();
//     }
//     /*
//     if (bmax * bmax * A + bmax * B + C < 0) {
//         return -std::numeric_limits<double>::infinity();
//     }
    
//     if ((A / 3) * std::pow(bmax, 3) + (B / 2) * std::pow(bmax, 2) + C * bmax != 1) {
//         return -std::numeric_limits<double>::infinity();
//     }
//     */
    
//     LP += std::log(b * b * A + b * B + C);
//     return LP;
// }

// ---------------------------------------------------------
// void lowEfit::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }


// Define the model function f(x)
double lowEfit::FitFunction(const std::vector<double>& params, double x, double sigma) {
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double A = params[3];
    // Quadratic term + Gaussian (with known mu (E_peak) and sigma)
    return a * x * x + b * x + c + A * std::exp(-0.5 * std::pow((x - E_peak) / sigma, 2));
}

// Provide the data (bin centers, observed counts, sigmas)
void lowEfit::SetData(const std::vector<std::vector<double>>& x, const std::vector<std::vector<double>>& y,
                      const std::vector<double>& resolutions) {
    BinHeight = y;
    BinCenter = x;
    sigmas = resolutions;
}


// Read the data (bin centers, observed counts, sigmas, ges)
void readCSV (){
    // to be implemented
}
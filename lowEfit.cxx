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

    AddParameter("a", -100, 100);   // Quadratic coefficient
    GetParameters().Back().SetPriorConstant();
    GetParameter("a").Fix(0); // Fixing to a linear

    AddParameter("b", -100, 100);   // Linear coefficient
    GetParameters().Back().SetPriorConstant();
    AddParameter("c", -10, 100);   // Constant term
    GetParameters().Back().SetPriorConstant();
    AddParameter("A", 0, 40);    // Amplitude of the Gaussian
    GetParameters().Back().SetPriorConstant();

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

    int exit = 0;

    // loop over detectors
    for (const auto& entry : sigma_map) {

        if (exit > 0) {
            break;
        }

        std::string detector = entry.first;

        double sigma = entry.second;

        std::vector<double> BinHeight = BinHeight_map[detector];
        std::vector<double> BinCenter = BinCenter_map[detector];
        unsigned int NBins = BinHeight.size();

        // loop over bins of our data
        for (unsigned int i = 200; i <= NBins; ++i) {
            // retrieve observed number of events
            double x = BinHeight[i];
            // retrieve bin center
            double m = BinCenter[i];
            // calculate expected number of events, using quadratic + gaussian function
            double nu = FitFunction(pars, m, sigma);
            // add to log-likelihood sum
            if (nu <= 0){
                LL += log(0.0);
            } else {
                LL += BCMath::LogPoisson(x, nu);
            }
        }

        exit += 1 ;
    } 
    // return log-likelihood
    return LL; 
}


// ---------------------------------------------------------
// double lowEfit::LogAPrioriProbability(const std::vector<double>& pars)
// {
// }

// ---------------------------------------------------------
// void lowEfit::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }


// Define the model function f(x)
double lowEfit::FitFunction(const std::vector<double>& pars, double x, double sigma) {
    double a = pars[0];
    double b = pars[1];
    double c = pars[2];
    double A = pars[3];
    // Quadratic term + Gaussian (with known mu (E_peak) and sigma)
    return a * x * x + b * x + c + A * std::exp(-0.5 * std::pow((x - E_peak) / sigma, 2));
}

// Provide the data (bin centers, observed counts, sigmas)
void lowEfit::SetData(const std::map<std::string, std::vector<double>>& x, const std::map<std::string, std::vector<double>>& y,
                     const std::map<std::string, double>& resolutions) {
    BinHeight_map = y;
    BinCenter_map = x;
    sigma_map = resolutions;

    // // DEBUG
    // for (const auto& entry : sigma_map) {
    //     std::string detector = entry.first;
    //     double sigma = entry.second;

    //     std::vector<double> BinHeight = BinHeight_map[detector];
    //     std::vector<double> BinCenter = BinCenter_map[detector];
    //     unsigned int NBins = BinHeight.size();

    //     std::cout << "======================" << std::endl;
    //     std::cout << "Detector: " << detector << " sigma: " << sigma << std::endl;
    //     // loop over bins of our data
    //     for (unsigned int i = 0; i <= NBins; ++i) {
    //         std::cout << "i: " << i << ", bin value: " << BinHeight[i] << ", bin center " << BinCenter[i] << std::endl;
    //     }
    //     std::cout << "======================" << std::endl;
    // }


}


// Read the data (bin centers, observed counts, sigmas, ges)
void readCSV (const std::string& filename,
              std::map<std::string, std::vector<double>>& bin_value_map, 
              std::map<std::string, std::vector<double>>& bin_center_map,
              std::map<std::string, double>& sigma_map){

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error in file opening!" << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Read the header line (if present)

    std::map<std::string, std::vector<double>> bin_value_map_tmp;
    std::map<std::string, std::vector<double>> bin_center_map_tmp;
    std::map<std::string, double> sigma_map_tmp;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string bin_center_str, bin_value_str, detector, sigma_str;

        // Leggi i valori da una riga del CSV
        std::getline(ss, bin_center_str, ',');
        std::getline(ss, bin_value_str, ',');
        std::getline(ss, detector, ',');
        std::getline(ss, sigma_str, ',');

        double bin_value = std::stod(bin_value_str);
        double bin_center = std::stod(bin_center_str);
        double sigma = std::stod(sigma_str);

        // Inserisci i valori nelle mappe in base al detector
        bin_value_map_tmp[detector].push_back(bin_value);
        bin_center_map_tmp[detector].push_back(bin_center);
        sigma_map_tmp[detector] = sigma;
    }

    bin_value_map = bin_value_map_tmp;
    bin_center_map = bin_center_map_tmp;
    sigma_map = sigma_map_tmp;

    file.close();

    // to be implemented
}
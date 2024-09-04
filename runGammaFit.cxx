// ***************************************************************
// This file was created using the bat-project script
// for project GammaFit.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>

#include "lowEfit.h"

int main()
{
    int Nch = 5;           // number of parallel MCMC chains
    int NIter = 1*100000;     // number of step per chain

    // input variables
    std::string input_filename = "detectors_data.csv";
    std::map<std::string, std::vector<double>> bin_value_map;
    std::map<std::string, std::vector<double>> bin_center_map;
    std::map<std::string, double> sigma_map;

    // read data
    readCSV(input_filename, bin_value_map, bin_center_map, sigma_map);

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new lowEfit object
    lowEfit m("QuadraticAndGaussian");
    
    // pass data
    m.SetData(bin_center_map, bin_value_map, sigma_map);

    // set marginalization method
    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    // Setting MC run iterations and number of parallel chains
    m.SetNIterationsRun(NIter);
    m.SetNChains(Nch);
    
    
    // Setting initial position for the parameters ## NEED TO BE checked
    std::vector<double> x0;
    x0.push_back(0);       // a
    x0.push_back(-0.05);   // b
    x0.push_back(40);      // c
    x0.push_back(5);       // A
    m.SetInitialPositions(x0);
    m.SetInitialPositionAttemptLimit(1000);


    BCLog::OutSummary("Model" + m.GetSafeName() + "created");

    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    // m.Normalize();

    // Write Markov Chain to a ROOT file as a TTree
    m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll();

    // run mode finding; by default using Minuit
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print summary plots
    // m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
    // m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
    // m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();

    // Check if the pre run has converged:
    int status = m.GetNIterationsConvergenceGlobal();
    std::cout << "Prerun status: " << status << std::endl;

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}

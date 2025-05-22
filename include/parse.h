#ifndef TYPES_H
#define TYPES_H

#include <iostream>
#include <filesystem>
#include <complex>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std::complex_literals;

using Clock   = std::chrono::high_resolution_clock;

using cdouble = std::complex<double>;
using Vec     = Eigen::VectorXd;
using Mat     = Eigen::MatrixXd;

using cVec    = Eigen::VectorXcd;
using cMat    = Eigen::MatrixXcd;
using cDMat   = Eigen::Matrix<cdouble, Eigen::Dynamic, Eigen::Dynamic>;

using SpMat   = Eigen::SparseMatrix<double>;
using cSpMat  = Eigen::SparseMatrix<std::complex<double>>;

using Triplet = Eigen::Triplet<std::complex<double>>;

enum class StateType
{
    Ground,
    Excited,
    Dicke, 
    Product, 
    Thermal, 
    CoherentSpin
};

enum class IntegratorType
{
    Euler,
    RK4, 
    RK5, 
    RK8, 
    BackwardEuler,
    CrankNicolson 
};

typedef struct 
{
    int N;
    int outs;
    int mtime;

    double hthr;
    double sthr;
    double omega;
    double beta;
    double theta;
    double phi;
    double gamma; // Control param gamma_p/gamma_m
    double h;
    double t_i;
    double t_f;

    bool fast;
    bool states;
    
    std::string outf;
    std::string jfile;

    StateType state;
    IntegratorType integrator;

} Params;
 

Params defaults()
{
    Params pars;
    pars.N     = 100;
    pars.outs  = 50;
    pars.mtime = 0;

    pars.hthr  = 0.1f;
    pars.sthr  = 0.01f;
    pars.omega = 1.0f;
    pars.beta  = 1.0f;
    pars.theta = 3.14/2.0;
    pars.phi   = 0.0f;
    pars.gamma = 0.8;

    pars.h   = 0.01;
    pars.t_i = 0.01;
    pars.t_f = 10.0f;

    pars.fast   = false;
    pars.states = false;

    pars.outf = "output";
    pars.jfile = "";

    pars.state = StateType::Product;
    pars.integrator = IntegratorType::RK4;

    return pars;
}

StateType parseState(const std::string& name)
{
    if (name == "Ground")  return StateType::Ground;
    if (name == "Excited") return StateType::Excited;
    if (name == "Dicke")   return StateType::Dicke;
    if (name == "Product") return StateType::Product;
    if (name == "Thermal") return StateType::Thermal;
    if (name == "CSS")     return StateType::CoherentSpin;

    throw std::invalid_argument("Unknown state: " + name);
}

IntegratorType parseIntegrator(const std::string& name)
{
    if (name == "Euler") return IntegratorType::Euler;
    if (name == "RK4")   return IntegratorType::RK4;
    if (name == "RK5")   return IntegratorType::RK5;
    if (name == "RK8")   return IntegratorType::RK8;
    if (name == "BE")    return IntegratorType::BackwardEuler;
    if (name == "CN")    return IntegratorType::CrankNicolson;

    throw std::invalid_argument("Unknown integrator: " + name);
}

void printHelp()
{
    std::cout << "\n";
    std::cout << "Parameters available: \n" << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Number of spins            --N              <int>    (Default: 100)"                         << std::endl;
    std::cout << "Find steady state only     --fast                                  "                         << std::endl;
    std::cout << "                                                                   "                         << std::endl; 
    std::cout << "Initial state              --state          <str>    Options: Ground, Excited"               << std::endl;
    std::cout << "                                                              Dicke, Product,"              << std::endl; 
    std::cout << "                                                              CSS, Thermal)  "              << std::endl; 
    std::cout << "                                                                   "                         << std::endl; 
    std::cout << "Time integrator            --integrator     <str>    Options: Explicit: Euler, RK4, RK5"     << std::endl; 
    std::cout << "                                                               Implicit: BE, CN"             << std::endl; 
    std::cout << "                                                     (Default: RK4)"                         << std::endl; 
    std::cout << "                                                                   "                         << std::endl; 
    std::cout << "Initial time               --ti             <dbl>    (Default: 0.01)"                        << std::endl; 
    std::cout << "Final time                 --tf             <dbl>    (Default: 1.0)"                         << std::endl; 
    std::cout << "Number of measurements     --meas           <int>    (Default: 50)   "                       << std::endl; 
    std::cout << "Set measurements times     --mtime          <int>    Options: 0: Linear"                     << std::endl; 
    std::cout << "                                                              1: Log (default)"              << std::endl; 
    std::cout << "                                                                   "                         << std::endl; 
    std::cout << "Timestep threshold         --thr            <dbl>    (Default: 0.5)"                         << std::endl; 
    std::cout << "Slope threshold            --sthr           <dbl>    (Default: 0.01)"                        << std::endl; 
    std::cout << "                                                                   "                         << std::endl; 
    std::cout << "Omega (energy unit)        --omega          <dbl>    (Default: 1.0)"                         << std::endl; 
    std::cout << "Inverse Temp (omega^-1)    --beta           <dbl>    (Default: 1.0)"                         << std::endl; 
    std::cout << "Control parameter (gp/gm)  --gamma          <dbl>    (Default: 0.9)"                         << std::endl; 
    std::cout << "                                                 "                                           << std::endl; 
    std::cout << "Output folder              --out            <str>"                                           << std::endl; 
    std::cout << "Save state matrix          --states         <str>"                                           << std::endl; 
    std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
}

int parseArgs(int argc, char* argv[], Params* pars)
{
    for (int i = 1; i < argc; ++i)
    {
        if      (!strcmp(argv[i], "--N") && i+1 < argc)       { pars->N       = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--meas") && i+1 < argc)    { pars->outs    = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--mtime") && i+1 < argc)   { pars->mtime   = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--thr") && i+1 < argc)     { pars->hthr    = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--sthr") && i+1 < argc)    { pars->sthr    = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--h") && i+1 < argc)       { pars->h       = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--ti") && i+1 < argc)      { pars->t_i     = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--tf") && i+1 < argc)      { pars->t_f     = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--omega") && i+1 < argc)   { pars->omega   = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--beta") && i+1 < argc)    { pars->beta    = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--theta") && i+1 < argc)   { pars->theta   = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--phi") && i+1 < argc)     { pars->phi     = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--gamma") && i+1 < argc)   { pars->gamma   = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--out") && i+1 < argc)     { pars->outf    = argv[++i];       }

        else if (!strcmp(argv[i], "--state") && i + 1 < argc)
        {
            std::string state_str = argv[++i];
            try{
                pars->state = parseState(state_str);
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << "\n";
                return 1;
            }
        }

        else if (!strcmp(argv[i], "--integrator") && i + 1 < argc)
        {
            std::string integrator_str = argv[++i];
            try{
                pars->integrator = parseIntegrator(integrator_str);
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << "\n";
                return 1;
            }
        }

        else if (!strcmp(argv[i], "--info"))
        {
            printHelp();
            exit(0);
        }

        else if (!strcmp(argv[i], "--fast"))
        {
            pars->fast=true;
        }
        
        else if (!strcmp(argv[i], "--states"))
        {
            pars->states=true;
        }

    }

    return 0;
}

std::string setDir(std::string base, int N, std::string state, std::string inte, double thr, double gnet)
{
    std::stringstream dir1;
    std::stringstream dir2;

    dir1 << base << "/N" << N;
    dir2 << dir1.str() << "/" << state << "_" << inte << "_h" 
         << thr << "_g" << gnet << ".txt";
    
    std::filesystem::create_directories(dir1.str());
    
    return dir2.str();
}


void printParams(Params* pars)
{
    std::string st_name  = "";
    std::string int_name = "";

    if (pars->state == StateType::Product)
        st_name = "Product";
    else if (pars->state == StateType::Ground)
        st_name = "Ground";
    else if (pars->state == StateType::Excited)
        st_name = "Excited";
    else if (pars->state == StateType::Dicke)
        st_name = "Dicke";
    else if (pars->state == StateType::Thermal)
        st_name = "Thermal";
    else if (pars->state == StateType::CoherentSpin)
        st_name = "Coherent Spin";

    if (pars->integrator == IntegratorType::Euler)
        int_name = "Euler";
    else if (pars->integrator == IntegratorType::RK5)
        int_name = "RK5";
    else if (pars->integrator == IntegratorType::RK8)
        int_name = "RK8";
    else if (pars->integrator == IntegratorType::BackwardEuler)
        int_name = "Backward Euler";
    else if (pars->integrator == IntegratorType::CrankNicolson)
        int_name = "Crank-Nicolson";
    else
        int_name = "RK4";

    std::string filename = setDir(pars->outf, pars->N, st_name, int_name, 
                                  pars->hthr, pars->gamma); 

    if (std::filesystem::exists(filename))
        std::filesystem::remove(filename);
    
    if (pars->states)
    {
        std::stringstream ss;
        ss << pars->outf << "/states";
        std::filesystem::create_directory(ss.str());
    }

    pars->jfile = filename;
    
    std::cout << "\nEigen is using " << Eigen::nbThreads() << " threads." << std::endl;
    std::cout << "                                "   << std::endl;
    std::cout << "                                "   << std::endl;
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "N              = " << pars->N          << std::endl; 
    std::cout << "state          = " << st_name          << std::endl; 
    std::cout << "integrator     = " << int_name         << std::endl; 
    std::cout << "ti             = " << pars->t_i        << std::endl; 
    std::cout << "tf             = " << pars->t_f        << std::endl; 
    std::cout << "omega          = " << pars->omega      << std::endl; 
    std::cout << "beta           = " << pars->beta       << std::endl; 
    std::cout << "theta          = " << pars->theta      << std::endl; 
    std::cout << "phi            = " << pars->phi        << std::endl; 
    std::cout << "gamma+/gamma_- = " << pars->gamma      << std::endl; 
    std::cout << "# meas         = " << pars->outs       << std::endl; 
    std::cout << "--------------------------------"   << std::endl;
    std::cout << "                                "   << std::endl;
    std::cout << "                                "   << std::endl;

}

void printBanner() {
    std::cout << "\033[1;36m";  // Light green
    std::cout << R"(
 _____                  _   _       
|  _  |                | \ | |      
| | | |_ __   ___ _ __ |  \| |_   _ 
| | | | '_ \ / _ \ '_ \| . ` | | | |
\ \_/ / |_) |  __/ | | | |\  | |_| |
 \___/| .__/ \___|_| |_\_| \_/\__,_|
      | |                           
      |_|                                       
    )" << std::endl;
    std::cout << "\033[0m";
}


Params init(int argc, char* argv[])
{
    Params params = defaults();
    parseArgs(argc, argv, &params); 
    printBanner();
    printParams(&params);
    return params;
}


#endif 

#ifndef _STATE_H__
#define _STATE_H__

#include "parse.h"


bool isMatrixFinite(const cMat& mat)
{
    return mat.allFinite(); 
}

// m = -N/2
cMat GroundState(int N) {
    cMat X0 = cMat::Zero(N + 1, N + 1);
    X0(N, N) = 1.0;
    return X0;
}

// m = +N/2 
cMat ExcitedState(int N) {
    cMat X0 = cMat::Zero(N + 1, N + 1);
    X0(0, 0) = 1.0;
    return X0;
}

// m = 0 
cMat DickeState(int N)
{
    if (N % 2 != 0)
    {
        std::cout << "Can't use odd N!" << std::endl;
        exit(0);
    }

    cMat X0 = cMat::Zero(N + 1, N + 1);
    X0(N / 2, N / 2) = cdouble(1.0, 0.0);
    return X0;
}


cMat ProductState(int N)
{
    int dim = N + 1;
    double S = N / 2.0;

    cMat X0 = cMat::Zero(dim, dim);

    std::vector<double> log_fact_table(N + 1);
    log_fact_table[0] = 0.0;
    for (int i = 1; i <= N; ++i)
        log_fact_table[i] = log_fact_table[i - 1] + std::log(i);

    double logNfact = log_fact_table[N];
    double log_pref = -N * std::log(2.0); 

    for (int i = 0; i < dim; ++i)
    {
        double m_i = S - i;
        int idx_i1 = static_cast<int>(S + m_i);
        int idx_i2 = static_cast<int>(S - m_i);
        double log_fac_i = logNfact - log_fact_table[idx_i1] - log_fact_table[idx_i2];

        for (int j = 0; j < dim; ++j)
        {
            double m_j = S - j;
            int idx_j1 = static_cast<int>(S + m_j);
            int idx_j2 = static_cast<int>(S - m_j);
            double log_fac_j = logNfact - log_fact_table[idx_j1] - log_fact_table[idx_j2];

            double log_val = log_pref + 0.5 * (log_fac_i + log_fac_j);
            X0(i, j) = std::exp(log_val);  
        }
    }

    if (!isMatrixFinite(X0))
    {
        std::cerr << "\nProduct state contains non-finite entries, shutting down!\n" << std::endl;
        exit(1);
    }
    return X0;
}


cMat CoherentSpinState(int N, double theta, double phi)
{
    int dim = N + 1;
    cVec coeffs(dim);

    std::vector<double> log_fact(N + 1, 0.0);
    for (int i = 1; i <= N; ++i)
        log_fact[i] = log_fact[i - 1] + std::log(i);

    for (int i = 0; i < dim; ++i)
    {
        int k = i;
        double log_binom = log_fact[N] - log_fact[k] - log_fact[N - k];

        double mag = std::exp(0.5 * log_binom +
                              (N - k) * std::log(std::cos(theta / 2.0)) +
                              k * std::log(std::sin(theta / 2.0)));

        cdouble phase = std::polar(1.0, k * phi);  // e^{i k phi}
        coeffs[i] = mag * phase;
    }

    coeffs.normalize();
    cMat rho = coeffs * coeffs.adjoint();

    if (!isMatrixFinite(rho))
    {
        std::cerr << "\nCSS state contains non-finite entries, shutting down!\n" << std::endl;
        exit(1);
    }

    return rho;
}

double computePolarisation(int N, double beta, double omega_0)
{
    int dim = N + 1;
    double S = N / 2.0;
    double Z = 0.0;
    double P = 0.0;

    for (int i = 0; i < dim; ++i)
    {
        double m = S - i;
        double E = omega_0 * m;
        double w = std::exp(-beta * E);
        Z += w;
        P += m * w;
    }

    double pol = std::abs(P/Z/S);
    return pol;
}


double findBeta(int N, double omega_0, double target_pol,
                               double tol = 1e-6, int max_iter = 100)
{
    if (target_pol >= 1.0) {
        std::cerr << "Warning: target polarisation = 1 not achievable for finite beta. "
                     "Returning large beta instead.\n";
        return 10; 
    }

    double beta_low = 0.0;
    double beta_high = 10.0;

    auto f = [&](double beta) {
        return computePolarisation(N, beta, omega_0) - target_pol;
    };

    while (f(beta_high) < 0 && beta_high < 10) {
        beta_high *= 2;
    }

    double f_low = f(beta_low);
    double f_high = f(beta_high);

    if (f_low > 0 || f_high < 0)
    {
        std::cerr << "Target polarisation is outside achievable range (0 <= p < 1).\n";
        exit(1);
    }

    for (int i = 0; i < max_iter; ++i)
    {
        double beta_mid = 0.5 * (beta_low + beta_high);
        double f_mid = f(beta_mid);

        if (std::abs(f_mid) < tol)
            return beta_mid;

        if (f_mid > 0)
            beta_high = beta_mid;
        else
            beta_low = beta_mid;
    }

    std::cerr << "Warning: beta not found within tolerance\n";
    return 0.5 * (beta_low + beta_high);
}

cMat ThermalState(int N, double pol_target, double omega_0)
{
    double beta = findBeta(N, omega_0, pol_target);
    std::cout << "\n\n[Thermal State] Found beta = " << beta << std::endl;

    int dim = N + 1;
    double S = N / 2.0;
    cMat X0 = cMat::Zero(dim, dim);

    double Z = 0.0;
    std::vector<double> weights(dim);

    for (int i = 0; i < dim; ++i)
    {
        double m = S - i;
        double E = omega_0 * m;
        double w = std::exp(-beta * E);
        weights[i] = w;
        Z += w;
    }

    double P = 0.0;
    for (int i = 0; i < dim; ++i)
    {
        double m = S - i;
        double prob = weights[i] / Z;
        X0(i, i) = prob;
        P += m * prob;
    }

    double pol = std::abs(P/S);
    std::cout << "[Thermal State] Target polarisation = "   << pol_target << std::endl;
    std::cout << "[Thermal State] Achieved polarisation = " << pol << std::endl;
    std::cout << "\n";

    return X0;
}


void writeState(const cMat& vstate, const std::string& fname)
{
    std::ofstream outFile(fname);
    if (!outFile.is_open())
    {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    size_t rows = vstate.rows(), cols = vstate.cols();
    outFile << "# Matrix dimensions: " << rows << " x " << cols << "\n\n";

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            const cdouble& element = vstate(i, j);
            outFile << element.real();
            if (j < cols - 1) outFile << "\t";
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Initial state matrix saved to " << fname << std::endl;
}


cMat initState(Params* pars, bool save = false) 
{
    int       N    = pars->N;
    StateType type = pars->state;
    
    cMat state;
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Initialising state for " << N << " spins  ... ";

    switch (type)
    {
        case StateType::Product:
            state = ProductState(N);
            if (save) { writeState(state, "ProductState.txt"); }
            break;

        case StateType::Dicke:
            state = DickeState(N);
            if (save) { writeState(state, "DickeState.txt"); }
            break;
            
        case StateType::Ground:
            state = GroundState(N);
            if (save) { writeState(state, "GroundState.txt"); }
            break;
        
        case StateType::Excited:
            state = ExcitedState(N);
            if (save) { writeState(state, "ExcitedState.txt"); }
            break;

        case StateType::Thermal:
        {
            double pol  = pars->pol;
            double om0  = pars->omega;
            state = ThermalState(N, pol, om0);
            if (save) { writeState(state, "ThermalState.txt"); }
            break;
        }
            
        case StateType::CoherentSpin:
        {
            double theta = pars->theta;
            double phi   = pars->phi;
            state = CoherentSpinState(N, theta, phi);
            if (save) { writeState(state, "CoherentSpinState.txt"); }
            break;
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::string time = formatDuration(dur);
    std::cout << "took " << time << std::endl;
    return state;
}

#endif

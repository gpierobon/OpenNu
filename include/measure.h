#ifndef _MEASURE_H__
#define _MEASURE_H__

#include "parse.h"
#include "utils.h"
#include "spinops.h"

void cacheL(const cSpMat& mat, const std::string& fname)
{
    std::ofstream file(fname);
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (cSpMat::InnerIterator it(mat, k); it; ++it)
        {
            file << it.row() << "," << it.col() << "," << it.value().real() << "," 
                 << it.value().imag() << "\n";
        }
    }
    file.close();
}

void cacheX(Params* pars, int meas, const cMat& vstate)
{
    std::stringstream ff;
    ff << pars->outf << "/states" << "/state" << "_" << std::setw(3) << std::setfill('0') << meas << ".txt"; 
    std::ofstream outFile(ff.str());
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
}

void cacheTraces(Params* pars, double& t, double* js)
{
    std::ofstream file(pars->jfile, std::ios::app);

    if (file.is_open())
    {
        file << t << " " << js[0] << " " << js[1] << " " << js[2]
                  << " " << js[3] << " " << js[4] << " " << js[5]
                  << " " << js[6] << " " << js[7] << " " << js[8]
                  << "\n";
        file.close();
    }
    else
        std::cerr << "Error: could not open file for writing.\n";

}


void steadyState(Params* pars, double& jz, double& sx)
{
    std::cout << "Computing steady state values of Jz, Sx... " << std::endl;
    
    int N = pars->N;
    int dim = N + 1;
    double S = N / 2.0;
    
    double ga = pars->gamma;

    std::vector<double> x(N + 1, 0.0);
    std::vector<double> main_diag(N + 1, 0.0);
    std::vector<double> lower_diag(N, 0.0);
    std::vector<double> upper_diag(N, 0.0);

    x[0] = 1.0;

    for (int i = 0; i <= N; ++i)
    {
        if (i < N)
        {
            lower_diag[i] = ga * (i + 1) * (N - i);
            main_diag[i]  = -ga * (i + 1) * (N - i) - i * (N - i + 1);
            upper_diag[i] =  (i + 1) * (N - i);
        }
        else
            main_diag[i] = -ga * (i + 1) * (N - i) - i * (N - i + 1);
    }

    for (int i = 1; i <= N; ++i)
    {
        if (i == 1)
            x[i] = -x[0] * main_diag[0] / upper_diag[0];
        else if (i == N)
            x[i] = -lower_diag[i - 1] * x[i - 1] / main_diag[i];
        else
            x[i] = -(x[i - 2] * lower_diag[i - 2] + x[i - 1] * main_diag[i - 1]) / upper_diag[i - 1];
    }

    // Normalize
    double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
    for (auto& xi : x)
        xi /= sum_x;
    
    Vec m_vals(dim);
    for (int i = 0; i < dim; ++i) {
        m_vals(i) = S - i;
    }

    double statjz = 0.0;
    double statsx = 0.0;
    std::vector<double> sigxp(N + 1, 0.0);

    for (int i = 0; i <= N; ++i)
    {
        statjz += x[i] * m_vals[i];
        sigxp[i] = 0.25 * x[i] *
                     ((S-m_vals(i)) * (S+m_vals(i)+1) + (S+m_vals(i)) * (S-m_vals(i)+1));
        statsx += sigxp[i];
    } 
    
    jz = statjz;
    sx = std::sqrt(statsx);
    std::cout << "\n";
    std::cout << "Jz = " << jz << std::endl;
    std::cout << "Sx = " << sx << std::endl;
}


void getTraces(const cMat& mat, const SpinOps* sOps, double* js)
{
    js[0] = (mat * sOps->Sx).trace().real();                   
    js[1] = (mat * sOps->Sy).trace().real();
    js[2] = (mat * sOps->Sz).trace().real();

    js[3] = (mat * sOps->Sx * sOps->Sx).trace().real();
    js[4] = (mat * sOps->Sy * sOps->Sy).trace().real();
    js[5] = (mat * sOps->Sz * sOps->Sz).trace().real();
    
    js[6] = (mat * sOps->Sp * sOps->Sm).trace().real();
    js[7] = (mat * sOps->Sm * sOps->Sp).trace().real();
    js[8] = js[3] + js[4] + js[5];      
}

void measure(Params* pars, const cMat& mstate, SpinOps* S, double t, size_t i, 
             size_t num_steps, std::chrono::duration<double> dur, int meas_count)
{
    int N = pars->N;
    double js[9];
    getTraces(mstate, S, js);
    
    double p = std::abs(js[2]/std::sqrt(js[8]));
    double jpm = js[6]/((N+1)*(N+1)/4);

    printStatus2(i, num_steps, t, p, jpm, dur, meas_count, pars->outs);
    cacheTraces(pars, t, js); 
    if (pars->states)
        cacheX(pars, meas_count, mstate);
    return;
}


#endif

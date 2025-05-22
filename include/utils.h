#ifndef _UTILS_H__
#define _UTILS_H__

#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>
#include "parse.h"

cVec Mat2Vec(const cMat& mat)
{
    cVec vec = Eigen::Map<const cVec>(mat.data(), mat.size());
    return vec;
}

cMat Vec2Mat(const cVec& vec, int rows, int cols)
{
    assert(vec.size() == rows * cols);
    return Eigen::Map<const cMat>(vec.data(), rows, cols);
}


void SpMatPrint(const SpMat& mat) {
    int rows = mat.rows();
    int cols = mat.cols();

    std::cout << "\n\n";
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (mat.coeff(i, j) != 0.0)
                std::cout << std::setprecision(3) << mat.coeff(i, j) << " ";
            else
                std::cout << " 0.000 ";
        }
        std::cout << std::endl;
    }
}

void cSpMatPrint(const cSpMat& mat) {
    int rows = mat.rows();
    int cols = mat.cols();

    std::cout << "\n\n";
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (mat.coeff(i, j) != cdouble(0.0, 0.0))
                std::cout << std::setprecision(3) << mat.coeff(i, j) << " ";
            else
                std::cout << " 0.000 ";
        }
        std::cout << std::endl;
    }
}

std::vector<size_t> generateMeasList(Params* pars, double dt, size_t num_steps)
{
    int t0 = pars->t_i;
    size_t n_bins = pars->outs;
    std::vector<size_t> bins;


    if (pars->mtime == 0)
    {
        double t_max = t0 + num_steps * dt;
        for (size_t k = 0; k < n_bins; ++k)
        {
            double t = t0 + dt; 
            double t_bin = t + (t_max - t) * k / (n_bins - 1);
            size_t i_bin = static_cast<size_t>((t_bin - t0) / dt);
            bins.push_back(i_bin);
        }
    }
    else
    {
        double log_t_max = std::log10(t0 + num_steps * dt);
        for (size_t k = 0; k < n_bins; ++k)
        {
            double log_t = std::log10(t0 + dt); 
            double t_bin = std::pow(10.0, log_t + (log_t_max - log_t) * k / (n_bins - 1));
            size_t i_bin = static_cast<size_t>((t_bin - t0) / dt);
            bins.push_back(i_bin);
        }
    }
    return bins;
}


std::string formatDuration(std::chrono::duration<double> dur) {
    using namespace std::chrono;

    double seconds = dur.count();
    std::ostringstream out;
    out << std::fixed << std::setprecision(1);

    if (seconds < 1.0)
    {
        int ms = static_cast<int>(seconds * 1000);
        out << ms << "ms";
    }
    else if (seconds < 60.0)
    {
        out << seconds << "s";
    }
    else if (seconds < 3600.0)
    {
        int mins = static_cast<int>(seconds) / 60;
        double rem = seconds - mins * 60;
        out << mins << "m " << rem << "s";
    }
    else
    {
        int hrs  = static_cast<int>(seconds) / 3600;
        int mins = (static_cast<int>(seconds) % 3600) / 60;
        double rem = seconds - hrs * 3600 - mins * 60;
        out << hrs << "h " << mins << "m " << rem << "s";
    }

    return out.str();
}


void printStatus2(int step, size_t num_steps, double t, double p, double jpm,
                 std::chrono::duration<double> dur, int meas, int num_meas)
{
    double percent = 100.0 * step / num_steps;
    std::string time_str = formatDuration(dur);

    std::cout << std::fixed    <<  std::setprecision(1)
              << std::setw(6)  << "Meas #"    << std::setw(4) << meas 
              << std::fixed    <<  std::setprecision(5)
              << std::setw(5)  << "    |    t:"     << " " << std::setw(6) << t
              << std::fixed    <<  std::setprecision(3)           
              << std::setw(6)  << "    |    P:"    << " " << std::setw(7) << p
              << std::setw(6)  << "    |    J+-:"  << " " << std::setw(7) << jpm 
              << std::setw(4)  << "    |  "        << std::setw(6) << std::setprecision(1) << percent << "%"
              << std::setw(12) << "Walltime:"  << " " << std::setw(12) << time_str << std::endl;
}



#endif

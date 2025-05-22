#include <iostream>
#include <stdexcept>

#include "nulindblad.h"

#define EIGEN_USE_OPENMP

int main(int argc, char* argv[])
{
    Params pars = init(argc, argv);
    Clock::time_point start = Clock::now();

    if (pars.fast)
    {
        double statJz, statSx;
        steadyState(&pars, statJz, statSx);
        exit(0);
    }

    // Initialise state and Liouvillian
    cMat mstate = initState(&pars);
    cVec vstate = Mat2Vec(mstate);
    int mdim    = std::sqrt(vstate.size());
    cSpMat L    = Liouvillian(&pars);
    
    SpinOps spOps(pars.N);
    spOps.build();
    
    // Compute largest eigenvalue and set timestep
    double rad       = spectral_radius(L);
    double thr       = pars.hthr;
    double h_check   = thr / rad;
    double dt        = std::min(pars.h, h_check);
    size_t num_steps = static_cast<size_t>((pars.t_f - pars.t_i) / dt);

    std::cout << "\nTime step set to: dt = " << h_check << std::endl;
    std::cout << "Running a loop with " << num_steps << " steps ... \n" << std::endl;

    std::vector<size_t> times = generateMeasList(&pars, dt, num_steps);
    double t0 = pars.t_i;
    int meas_count = 0;

    // First measurement 
    Clock::time_point curr = Clock::now();
    auto ti   = curr - start;
    measure(&pars, mstate, &spOps, t0, 0, num_steps, ti, meas_count);

    // Evolve 
    for (size_t i = 0; i < num_steps; ++i)
    {
        double t = t0 + i * dt;

        evolve(pars.integrator, vstate, L, dt, t);
                                                            
        if (std::binary_search(times.begin(), times.end(), i))
        {
            meas_count += 1;
            Clock::time_point curr = Clock::now();
            auto dur = curr - start;
            cMat mstate = Vec2Mat(vstate, mdim, mdim);
            measure(&pars, mstate, &spOps, t, i, num_steps, dur, meas_count);
        }
    }

    Clock::time_point end = Clock::now();
    auto wall = end - start;
    std::string time_str = formatDuration(wall);
    std::cout << "\nFinished: it took " << time_str << std::endl;  

    std::cout << "\nOutput saved to " << pars.jfile << std::endl;
    std::cout << "(t, Jx, Jy, Jz, Jx^2, Jy^2, Jz^2, Jpm, Jmp, J)" << std::endl;
    
    return 0;
}

#ifndef _EVOLVE_H__
#define _EVOLVE_H__

#include "parse.h"

cVec odefun(double t, const cVec& X_vec, const cSpMat& lind)
{
    return lind * X_vec;
}

void Euler_step(cVec& X, const cSpMat& lind, double h)
{
    X += h * (lind * X);
}

void RK4_step(cVec& X, const cSpMat& lind, double t, double h)
{
    cVec k1 = odefun(t, X, lind);
    cVec k2 = odefun(t + 0.5 * h, X + 0.5 * h * k1, lind);
    cVec k3 = odefun(t + 0.5 * h, X + 0.5 * h * k2, lind);
    cVec k4 = odefun(t + h, X + h * k3, lind);

    X += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

    
// RK Cash-Karp method
void RK5_step(cVec& X, const cSpMat& lind, double t, double h)
{
    static const double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;
    static const double b21 = 0.2;
    static const double b31 = 3.0 / 40.0,  b32 = 9.0 / 40.0;
    static const double b41 = 0.3,        b42 = -0.9,       b43 = 1.2;
    static const double b51 = -11.0 / 54.0, b52 = 2.5,      b53 = -70.0 / 27.0, b54 = 35.0 / 27.0;
    static const double b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
                         b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0;

    static const double c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0,
                         c6 = 512.0 / 1771.0;

    cVec k1 = odefun(t, X, lind);
    cVec k2 = odefun(t + a2 * h, X + h * b21 * k1, lind);
    cVec k3 = odefun(t + a3 * h, X + h * (b31 * k1 + b32 * k2), lind);
    cVec k4 = odefun(t + a4 * h, X + h * (b41 * k1 + b42 * k2 + b43 * k3), lind);
    cVec k5 = odefun(t + a5 * h, X + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4), lind);
    cVec k6 = odefun(t + a6 * h, X + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5), lind);

    X += h * (c1 * k1 + c3 * k3 + c4 * k4 + c6 * k6);
}
    

// Dormand-Prince 8th-order method
void RK8_step(cVec& X, const cSpMat& lind, double t, double h)
{
    static const double a[12][12] =
    {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {3.0 / 10.0, 3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {4.0 / 5.0, 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {8.0 / 9.0, 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };

    static const double b[12] =
    {
        5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0,
        -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    };

    //static const double b_err[12] =
    //{
    //    35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0,
    //    -2187.0 / 6784.0, 11.0 / 84.0, 0.0, 0.0,
    //    0.0, 0.0, 0.0, 0.0
    //};

    cVec k[12];
    k[0] = odefun(t, X, lind);

    for (int i = 1; i < 12; ++i)
    {
        cVec sum = cVec::Zero(X.size());

        for (int j = 0; j < i; ++j)
            sum += a[i][j] * k[j];
        
        k[i] = odefun(t + h * (i / 12.0), X + h * sum, lind);
    }

    cVec X_new = X;
    for (int i = 0; i < 12; ++i)
        X_new += h * b[i] * k[i];

    //cVec err = 0.0;
    //for (int i = 0; i < 12; ++i)
    //    err += h * b_err[i] * k[i];

    X = X_new;
}



void backward_euler_step(cVec& X, const cSpMat& L, double h)
{
    const int dim = L.rows();
    cSpMat I(dim, dim);
    I.setIdentity(); 

    cSpMat A = I - h * L;

    Eigen::SparseLU<cSpMat> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed in Backward Euler step");
    }

    X = solver.solve(X);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed in Backward Euler step");
    }
}

void crank_nicolson_step(cVec& X, const cSpMat& L, double h)
{
    const int dim = L.rows();
    cSpMat I(dim, dim);
    I.setIdentity(); 

    cSpMat A = I - (h / 2.0) * L;
    cSpMat B = I + (h / 2.0) * L;
    cVec rhs = B * X;

    Eigen::SparseLU<cSpMat> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Decomposition failed in Crank-Nicolson step");
    }

    X = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed in Crank-Nicolson step");
    }
    
}

void evolve(IntegratorType type, cVec &X, const cSpMat& L, double h, double t)
{
    switch (type)
    {
        case IntegratorType::Euler:
            Euler_step(X, L, h);
            break;

        case IntegratorType::RK4:
            RK4_step(X, L, t, h);
            break;

        case IntegratorType::RK5:
            RK5_step(X, L, t, h);
            break;
        
        case IntegratorType::RK8:
            RK8_step(X, L, t, h);
            break;
       
        case IntegratorType::BackwardEuler:
            backward_euler_step(X, L, h);
            break;

        case IntegratorType::CrankNicolson:
            crank_nicolson_step(X, L, h);
            break;
    }

}

#endif

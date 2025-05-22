#ifndef SPINOPS_H
#define SPINOPS_H

#include "parse.h"

void buildSz(int N, SpMat& Sz)
{
    int   dim = N+1;
    double S  = 0.5*N;

    for (int j = 0; j < dim; ++j)
        Sz.insert(j,j) = S-j;

    Sz.makeCompressed();
}


void buildSpm(int N, cSpMat& Sp, cSpMat& Sm)
{
    int    dim = N+1;
    double S   = 0.5*N;

    std::vector<double> m_vals(dim);
    for (int i = 0; i < dim; ++i)
        m_vals[i] = S - i;

    for (int j = 0; j < dim - 1; ++j)
    {
        double val = std::sqrt((S - m_vals[j + 1]) * (S + m_vals[j + 1] + 1));

        if (val != 0.0)
        {
            Sp.insert(j, j + 1) = cdouble(val, 0.0);
            Sm.insert(j + 1, j) = cdouble(val, 0.0);
        }
    }

    Sp.makeCompressed();
    Sm.makeCompressed();
    //printSparseMatrix(Sp);
    //printSparseMatrix(Sm);
}

void buildSxy(int N, const cSpMat& Sp, const cSpMat& Sm, cSpMat& Sx, cSpMat& Sy)
{
    int dim = N + 1;

    Sx = cSpMat(dim, dim);
    Sy = cSpMat(dim, dim);

    Sx = 0.5 * (Sp + Sm);

    Sy = cdouble(0, -0.5) * (Sp - Sm);

    Sx.makeCompressed();
    Sy.makeCompressed();
}


struct SpinOps
{
    int N;
    double S;
    int dim;

    SpMat  Sz;
    cSpMat Sx, Sy, Sp, Sm;

    SpinOps(int N_) : N(N_), S(0.5 * N_), dim(N_ + 1)
    {
        Sx = cSpMat(dim, dim);
        Sy = cSpMat(dim, dim);
        Sz =  SpMat(dim, dim);
        Sp = cSpMat(dim, dim);
        Sm = cSpMat(dim, dim);
    }

    void build()
    {
        buildSz (N, Sz);    
        buildSpm(N, Sp, Sm); 
        buildSxy(N, Sp, Sm, Sx, Sy);     
    }

};

#endif

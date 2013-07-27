#include <iostream>

#include "thermoelasticity1dp.h"

int main()
{
    DiffusionData t_data;
    ElasticityData e_data;

    t_data.alpha = 1.6e-6Q;

    t_data.k = 401.q;

    e_data.Ax = 1.0e-4Q;
    e_data.E = 1.1e+11Q;

    tFloat q_ = 5.0e+4Q;

    t_data.heatSource = new PolynomialConstant(-q_/t_data.k);

    tFloat t0 = 20.0q;
    tFloat tl = 30.0q;

    Boundary1D left(0.0q, Dirichlet, t0);
    Boundary1D right(1.0q, Dirichlet, tl);

    Thermoelasticity1Dp mesh(20001, left, right, e_data, t_data);

    PolynomialQuadratic AS_t(t0, tl-t0+0.5q*q_/t_data.k, -0.5q*q_/t_data.k);

    PolynomialCubic AS_u(0.0q,
                         t_data.alpha/6.0q*(-3.0q*(tl-t0+0.5q*q_/t_data.k)+q_/t_data.k),
                         t_data.alpha/2.0q*(tl-t0+0.5q*q_/t_data.k),
                         -t_data.alpha*q_/(t_data.k*6.0q)
                         );

    PolynomialQuadratic AS_e(t_data.alpha/6.0q*(-3.0q*(tl-t0+0.5q*q_/t_data.k)+q_/t_data.k),
                         t_data.alpha*(tl-t0+0.5q*q_/t_data.k),
                         -t_data.alpha*q_/(t_data.k*2.0q)
                         );

    PolynomialQuadratic AS_s(e_data.E*(t_data.alpha/6.0q*(-3.0q*(tl-t0+0.5q*q_/t_data.k)+q_/t_data.k)),
                         e_data.E*(t_data.alpha*(tl-t0+0.5q*q_/t_data.k) - t_data.alpha*(tl-t0+0.5q*q_/t_data.k)),
                         e_data.E*(-t_data.alpha*q_/(t_data.k*2.0q) - t_data.alpha*(-0.5q*q_/t_data.k))
                         );

    tFloat AS_f = e_data.E*t_data.alpha*e_data.Ax/6.0q*(3.0q*(t0-tl) - 0.5q*q_/t_data.k);

    mesh.solver();

    mesh.printSolution(AS_t, AS_u, AS_e, AS_s, AS_f);

    mesh.thermo->plotSolution(AS_t);
    mesh.plotSolution(AS_u);
    mesh.plotSecondarySolutions(AS_e, AS_s);

    std::cout<<std::endl<<std::endl<<std::endl;
}


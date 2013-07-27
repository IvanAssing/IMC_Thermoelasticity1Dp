#ifndef THERMOELASTICITY1DP_H
#define THERMOELASTICITY1DP_H

#include "diffusion1dp.h"


// Classe para coleta de dados dos problemas
class ElasticityData
{
    public:
        tFloat Ax;
        tFloat E;

        ElasticityData(){
            Ax = 0.0q;
            E = 0.0q;
        }
};


class Thermoelasticity1Dp
{
    public:
        Boundary1D boundaryLeft, boundaryRight;
        Diffusion1Dp *thermo;
        ElasticityData eData;
        TDMA equationsSystem;
        tFloat l, h;
        tInteger n;

        tFloat *strain, *stress;

        Thermoelasticity1Dp();
        Thermoelasticity1Dp(tInteger nodes, Boundary1D left, Boundary1D right, ElasticityData eData, DiffusionData dData);

        void plotSolution(Functor1D &analyticalSolution);
        void printSolution(Functor1D &aT, Functor1D &aU, Functor1D &aE, Functor1D &aS, tFloat aF);
        void plotSecondarySolutions(Functor1D &AS_Strain, Functor1D &AS_Stress);
        void solver(void);
};

#endif // THERMOELASTICITY1DP_H

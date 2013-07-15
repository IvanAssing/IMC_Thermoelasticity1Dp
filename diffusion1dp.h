#ifndef DIFFUSION1DP_H
#define DIFFUSION1DP_H

#include "imc_dfm.h"
#include "tdma.h"

#include "functor1d.h"
#include "boundary.h"


// Classe para coleta de dados dos problemas
class DiffusionData
{
    public:
        tFloat k;
        tFloat H;
        tFloat Tinf;
        tFloat Tb;
        tFloat Ab;
        tFloat P;
        tFloat mi;
        tFloat C;
        Functor1D *heatSource;
        tFloat alpha;


        DiffusionData(){
            k = 0.0q;
            H = 0.0q;
            Tinf = 0.0q;
            Tb = 0.0q;
            Ab = 0.0q;
            P = 0.0q;
            mi = 0.0q;
            C = 0.0q;
            alpha = 0.0q;
            heatSource = new PolynomialConstant(0.0q);
        }
};


// Tipos de problemas
enum DiffusionProblem{
    ParedePlana,
    Aleta,
    QML
};

// Classe para resolução dos problemas de Difusão 1D Permanente
class Diffusion1Dp
{
    public:
        Boundary1D boundaryLeft, boundaryRight;
        DiffusionData data;
        DiffusionProblem problemType;
        tFloat l, h;
        tInteger n;

    public:
        TDMA equationsSystem;
        tFloat  averageValue, maxValue, heatFlowLeft, heatFlowRight;

        Diffusion1Dp(DiffusionProblem type, tInteger nodes, Boundary1D left, Boundary1D right, DiffusionData data);

        void plotSolution(Functor1D &analyticalSolution);
        void printSolution(Functor1D &analyticalSolution);
        void printSecondaryResults(tFloat AS_average = 0.0q, tFloat AS_left = 0.0q, tFloat AS_right = 0.0q,
                                   tFloat AS_left_1 = 0.0q, tFloat AS_right_1 = 0.0q);

        void solver(void);

};



#endif // DIFFUSION1DP_H

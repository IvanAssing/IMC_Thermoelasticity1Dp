#ifndef TDMA_H
#define TDMA_H

#include <stdio.h>
#include "imc_dfm.h"

#define NCOEF 4


// Classe para resolução de Sistema de Equações Lineares tipo TDMA
class TDMA
{
    private:
        tFloat **equation;
        tInteger neq;
        tInteger nEqMax;
        tFloat *P,*Q,*T;

    public:
        TDMA(void):nEqMax(0){}
        TDMA(tInteger nEquations);
        void addEquation(tFloat ap, tFloat aw, tFloat ae, tFloat bp);
        void operator ()(tFloat ap, tFloat aw, tFloat ae, tFloat bp);

        void solver(void);
        void printCoefficients(void);
        void printSolution(void);
        void setMaxEquations(tInteger nEquations);

        inline tFloat getP(tInteger i) { return P[i];}
        inline tFloat getQ(tInteger i) { return Q[i];}
        inline tFloat getT(tInteger i) { return T[i];}
        inline void setT(tInteger i,tFloat Tp) { T[i] = Tp;}
        inline tInteger getNeq(void) { return neq;}

        inline void emptyEquations(void) { neq = 0;}
        inline tFloat getAw(tInteger i) { return equation[i][1];}
        inline tFloat getAp(tInteger i) { return equation[i][0];}
        inline tFloat getAe(tInteger i) { return equation[i][2];}
        inline tFloat getBp(tInteger i) { return equation[i][3];}

    public:
        virtual ~TDMA();
};

#endif // TDMA_H


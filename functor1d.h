#ifndef FUNCTOR1D_H
#define FUNCTOR1D_H

#include "imc_dfm.h"

// Objeto Funcão (Super Classe Abstrata)
class Functor1D
{
    public:
        Functor1D();
        virtual tFloat operator()(tFloat x) = 0;
};

// Objeto Polinômio Constante
class PolynomialConstant : public Functor1D
{
    public:
        tFloat a0;
        PolynomialConstant(tFloat a0);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 1
class PolynomialLinear : public Functor1D
{
    public:
        tFloat a0, a1;
        PolynomialLinear(tFloat a0, tFloat a1);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 2
class PolynomialQuadratic : public Functor1D
{
    public:
        tFloat a0, a1, a2;
        PolynomialQuadratic(tFloat a0, tFloat a1, tFloat a2);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 3
class PolynomialCubic : public Functor1D
{
    public:
        tFloat a0, a1, a2, a3;
        PolynomialCubic(tFloat a0, tFloat a1, tFloat a2, tFloat a3);
        virtual tFloat operator()(tFloat x);
};

// Objeto Polinômio de Grau 4
class PolynomialQuartic : public Functor1D
{
    public:
        tFloat a0, a1, a2, a3, a4;
        PolynomialQuartic(tFloat a0, tFloat a1, tFloat a2, tFloat a3, tFloat a4);
        virtual tFloat operator()(tFloat x);
};

// Objeto Função para Temperatura Analítica em Aleta
class SpecialFunctorTA : public Functor1D
{
    public:
        tFloat k;
        tFloat H;
        tFloat Tinf;
        tFloat Tb;
        tFloat Ab;
        tFloat P;
        tFloat L;
        tFloat m;

        SpecialFunctorTA(tFloat k,
                         tFloat H,
                         tFloat Tinf,
                         tFloat Tb,
                         tFloat Ab,
                         tFloat P,
                         tFloat L);
        virtual tFloat operator()(tFloat x);
};

// Objeto Função para Fluxo de Calor Analítico em Aleta
class SpecialFunctorQA : public Functor1D
{
    public:
        tFloat k;
        tFloat H;
        tFloat Tinf;
        tFloat Tb;
        tFloat Ab;
        tFloat P;
        tFloat L;
        tFloat m;

        SpecialFunctorQA(tFloat k,
                         tFloat H,
                         tFloat Tinf,
                         tFloat Tb,
                         tFloat Ab,
                         tFloat P,
                         tFloat L);
        virtual tFloat operator()(tFloat x);
};

#endif // FUNCTOR1D_H

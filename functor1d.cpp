#include "functor1d.h"

Functor1D::Functor1D()
{
}

PolynomialConstant::PolynomialConstant(tFloat _a0)
    :a0(_a0)
{

}

tFloat PolynomialConstant::operator ()(tFloat x)
{
    return a0;
}

PolynomialLinear::PolynomialLinear(tFloat _a0, tFloat _a1)
    :a0(_a0), a1(_a1)
{

}

tFloat PolynomialLinear::operator ()(tFloat x)
{
    return a1*x+a0;
}

PolynomialQuadratic::PolynomialQuadratic(tFloat _a0, tFloat _a1, tFloat _a2)
    :a0(_a0), a1(_a1), a2(_a2)
{

}

tFloat PolynomialQuadratic::operator ()(tFloat x)
{
    return (a2*x+a1)*x+a0;
}

PolynomialCubic::PolynomialCubic(tFloat _a0, tFloat _a1, tFloat _a2, tFloat _a3)
    :a0(_a0), a1(_a1), a2(_a2), a3(_a3)
{

}

tFloat PolynomialCubic::operator ()(tFloat x)
{
    return ((a3*x+a2)*x+a1)*x+a0;
}

PolynomialQuartic::PolynomialQuartic(tFloat _a0, tFloat _a1, tFloat _a2, tFloat _a3, tFloat _a4)
    :a0(_a0), a1(_a1), a2(_a2), a3(_a3), a4(_a4)
{

}

tFloat PolynomialQuartic::operator ()(tFloat x)
{
    return (((a4*x+a3)*x+a2)*x+a1)*x+a0;
}


SpecialFunctorTA::SpecialFunctorTA(tFloat _k,
                                   tFloat _H,
                                   tFloat _Tinf,
                                   tFloat _Tb,
                                   tFloat _Ab,
                                   tFloat _P,
                                   tFloat _L)
    :k(_k), H(_H), Tinf(_Tinf), Tb(_Tb), Ab(_Ab), P(_P), L(_L)
{
#ifndef QUAD_PRECISION
    m = sqrt(_H*_P/(_k*_Ab));
#else
    m = sqrtq(_H*_P/(_k*_Ab));
#endif
}

tFloat SpecialFunctorTA::operator ()(tFloat x)
{
#ifndef QUAD_PRECISION
    return (Tinf + (Tb - Tinf)*((cosh(m*(L-x))+(H/(m*k))*sinh(m*(L-x)))/(cosh(m*L)+(H/(m*k))*sinh(m*L))));
#else
    return (Tinf + (Tb - Tinf)*((coshq(m*(L-x))+(H/(m*k))*sinhq(m*(L-x)))/(coshq(m*L)+(H/(m*k))*sinhq(m*L))));
#endif
}

SpecialFunctorQA::SpecialFunctorQA(tFloat _k,
                                   tFloat _H,
                                   tFloat _Tinf,
                                   tFloat _Tb,
                                   tFloat _Ab,
                                   tFloat _P,
                                   tFloat _L)
    :k(_k), H(_H), Tinf(_Tinf), Tb(_Tb), Ab(_Ab), P(_P), L(_L)
{
#ifndef QUAD_PRECISION
    m = sqrt(_H*_P/(_k*_Ab));
#else
    m = sqrtq(_H*_P/(_k*_Ab));
#endif
}

tFloat SpecialFunctorQA::operator ()(tFloat x)
{
#ifndef QUAD_PRECISION
    return (sqrt(H*P*k*Ab)*(Tb - Tinf)*((sinh(m*L)+(H/(m*k))*cosh(m*L))/(cosh(m*L)+(H/(m*k))*sinh(m*L))));
#else
    return (sqrtq(H*P*k*Ab)*(Tb - Tinf)*((sinhq(m*L)+(H/(m*k))*coshq(m*L))/(coshq(m*L)+(H/(m*k))*sinhq(m*L))));
#endif
}


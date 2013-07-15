#include "boundary.h"

Boundary1D::Boundary1D()
{
}

Boundary1D::Boundary1D(tFloat _x, BoundaryCondition _type, tFloat _bcValue)
{
    x = _x;
    type = _type;
    bcValue = _bcValue;
}

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "imc_dfm.h"


// Tipos de condição de contorno
enum BoundaryCondition{
    Dirichlet,
    Neumann,
    Robin
};

// Objeto para definir parametros de contorno
class Boundary1D
{

    public:
        tFloat x, bcValue;
        BoundaryCondition type;
        Boundary1D();
        Boundary1D(tFloat x, BoundaryCondition type = Dirichlet, tFloat bcValue = 0.0);
};

#endif // BOUNDARY_H

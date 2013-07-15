
#include "tdma.h"

// Construtor do objeto TDMA
TDMA::TDMA(tInteger n):nEqMax(n+2)
{
    // Alocação de memória
    // Aloca matriz para os coeficientes ap,aw,ae,bp
    equation = new tFloat*[nEqMax];
    for(tInteger p=0;p<nEqMax;p++)
        equation[p] = new tFloat[NCOEF];

    neq = 0; // Nenhuma equação adicionada

    // Aloca vetores
    P = new tFloat[nEqMax];
    Q = new tFloat[nEqMax];
    T = new tFloat[nEqMax];
}

// Redimensiona tamanho esperado para o sistema
void TDMA::setMaxEquations(tInteger n)
{
    if(nEqMax){
        delete P;
        delete Q;
        delete T;

        for(tInteger p=0;p<nEqMax;p++)
            delete equation[p];
        delete equation;
    }

    nEqMax = n;

    equation = new tFloat*[nEqMax];
    for(tInteger p=0;p<nEqMax;p++)
        equation[p] = new tFloat[NCOEF];

    neq = 0;

    P = new tFloat[nEqMax];
    Q = new tFloat[nEqMax];
    T = new tFloat[nEqMax];
}


// Método para adicionar equação
void TDMA::addEquation(tFloat ap,tFloat aw,tFloat ae,tFloat bp)
{
    // Verifica excesso de equações
    if(neq+1>nEqMax)
    {
        printf("Erro, número máximo de equações alcançado!!!!");
        return;
    }

    ++neq; // Incrementa contador

    // Adiciona coeficientes
    equation[neq-1][0] = ap;
    equation[neq-1][1] = aw;
    equation[neq-1][2] = ae;
    equation[neq-1][3] = bp;

}


// Operador para adicionar equação
void TDMA::operator()(tFloat ap,tFloat aw,tFloat ae,tFloat bp)
{
    this->addEquation(ap,aw,ae,bp);
}


// Imprimir coeficientes
void TDMA::printCoefficients(void)
{
    std::cout<<std::endl<<std::endl;
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
    std::cout<<std::endl<<std::setw(5)<<std::right<<"p";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"ap";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"aw";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"ae";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"bp";
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');

    std::cout.precision(OUT_FLOAT_PRECISION);
    std::cout<<std::scientific;


#ifdef QUAD_PRECISION
    char str[1000];
    for(tInteger i=0;i<neq;i++){
        std::cout<<std::endl;
        std::cout<<std::setw(5)<<i;
        quadmath_snprintf(str, 1000, Q_FORMAT,equation[i][0]);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,equation[i][1]);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,equation[i][2]);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,equation[i][3]);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
    }
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
#else
    for(tInteger i=0;i<neq;i++){
        std::cout<<std::endl;
        std::cout<<std::setw(5)<<i;
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<equation[i][0];
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<equation[i][1];
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<equation[i][2];
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<equation[i][3];
    }
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
#endif
}


// Solver do sistema TDMA
void TDMA::solver(void)
{
    // Passo 1
    P[0] = equation[0][2]/equation[0][0]; // Equação (7)
    Q[0] = equation[0][3]/equation[0][0]; // Equação (8)

    // Passo 2
    for(tInteger i=1;i<neq;i++)
    {
        P[i] = equation[i][2]/(equation[i][0]-equation[i][1]*P[i-1]); // Equação (5)
        Q[i] = (equation[i][3]+equation[i][1]*Q[i-1])/(equation[i][0]-equation[i][1]*P[i-1]); // Equação (6)
    }

    // Passo 3
    T[neq-1] = Q[neq-1]; // Equação (9)

    // Passo 4
    for(tInteger i=neq-2;i>=0;i--)
        T[i] = P[i]*T[i+1] + Q[i]; // Equação (2)
}


// Imprime solução na tela
void TDMA::printSolution(void)
{
    std::cout<<std::endl;
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+3*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
    std::cout<<std::endl<<std::setw(5)<<std::right<<"p";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"Pp";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"Qp";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"Tp";
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+3*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');

    std::cout.precision(OUT_FLOAT_PRECISION);
    std::cout<<std::scientific;

#ifdef QUAD_PRECISION
    char str[1000];
    for(tInteger i=0;i<neq;i++){
        std::cout<<std::endl;
        std::cout<<std::setw(5)<<i;
        quadmath_snprintf(str, 1000, Q_FORMAT,P[i]);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,Q[i]);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,T[i]);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
    }
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+3*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
#else
    for(tInteger i=0;i<neq;i++){
        std::cout<<std::endl;
        std::cout<<std::setw(5)<<i;
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<P[i];
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<Q[i];
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<T[i];
    }
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+3*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
#endif
}


// Libera memória alocada pelo objeto
TDMA::~TDMA()
{
    if(nEqMax){
        delete P;
        delete Q;
        delete T;

        for(tInteger p=0;p<nEqMax;p++)
            delete equation[p];
        delete equation;
    }
}



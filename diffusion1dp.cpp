#include "diffusion1dp.h"

Diffusion1Dp::Diffusion1Dp(DiffusionProblem _problemType, tInteger _nodes, Boundary1D _left,
                           Boundary1D _right, DiffusionData _data)
{
    problemType = _problemType;
    boundaryLeft = _left;
    boundaryRight = _right;
    data = _data;
    n = _nodes;

    // Calcula L e h
    l = boundaryRight.x - boundaryLeft.x;
    h = l/(n - 1);

    if(boundaryLeft.type != Dirichlet && boundaryRight.type != Dirichlet){
        std::cout<<"Problema inconsistente, solução indeterminada!";
        return;
    }

}


void Diffusion1Dp::solver(void)
{
    // Configura o número de equações do solver TDMA
    equationsSystem.setMaxEquations(n);

    // DISCRETIZAÇÃO PARA O PROBLEMA DE DIFUSÃO DE CALOR EM PAREDA PLANA
    if(problemType == ParedePlana){

        // Aplica a condição de contorno esquerdo
        if(boundaryLeft.type == Dirichlet)
            equationsSystem(1.0q, 0.0q, 0.0q, boundaryLeft.bcValue); // Se Dirichlet
        else
            equationsSystem(1.0q, 0.0q, 1.0q, boundaryLeft.bcValue*h); // Se Neumann

        // Aplica discretização de d2T/dx2 = S(x), com CDS-2
        for(tInteger i=1; i<n-1; i++){
            tFloat xp = boundaryLeft.x+i*h;
            equationsSystem(2.0q, 1.0q, 1.0q, -h*h*data.heatSource->operator ()(xp));
        }

        // Aplica a condição de contorno direito
        if(boundaryRight.type == Dirichlet)
            equationsSystem(1.0q, 0.0q, 0.0q, boundaryRight.bcValue); // Se Dirichlet
        else
            equationsSystem(1.0q, 1.0q, 0.0q, -boundaryRight.bcValue*h); // Se Neumann
    }

    // DISCRETIZAÇÃO PARA O PROBLEMA DE DIFUSÃO DE CALOR ALETA
    if(problemType == Aleta){

        // Aplica a condição de contorno esquerdo
        if(boundaryLeft.type == Dirichlet)
            equationsSystem(1.0q, 0.0q, 0.0q, data.Tb); // Se Dirichlet
        else
            equationsSystem(data.k-data.H*h, 0.0q, data.k, -h*data.H*data.Tinf); // Se Robin

        tFloat m2_h2 = (data.H*data.P*h*h)/(data.k*data.Ab);
        tFloat ap = 2.0q+m2_h2;
        tFloat bp = m2_h2*data.Tinf;

        // Aplica discretização de d2T/dx2 = S(T), com CDS-2
        for(tInteger i=1; i<n-1; i++)
            equationsSystem(ap, 1.0q, 1.0q, bp);

        // Aplica a condição de contorno direito
        if(boundaryRight.type == Dirichlet)
            equationsSystem(1.0q, 0.0q, 0.0q, data.Tb); // Se Dirichlet
        else
            equationsSystem(data.k+data.H*h, data.k, 0.0q, h*data.H*data.Tinf); // Se Robin
    }

    // DISCRETIZAÇÃO PARA O PROBLEMA DE DIFUSÃO DE QUANTIDADE DE MOVIMENTO LINEAR
    if(problemType == QML){

        // Aplica a condição de contorno esquerdo
        equationsSystem(1.0q, 0.0q, 0.0q, 0.0q); // Dirichlet

        tFloat bp = -data.C*h*h/data.mi;

        // Aplica discretização de mi.d2u/dy2 = C, com CDS-2
        for(tInteger i=1; i<n-1; i++)
            equationsSystem(2.0q, 1.0q, 1.0q, bp);

        // Aplica a condição de contorno direito
        equationsSystem(1.0q, 0.0q, 0.0q, 0.0q); // Dirichlet
    }


    // Resolve o sistema de equações lineares
    equationsSystem.solver();

    // Calcula a temperatura média
    averageValue = 0.0q;
    for(tInteger i=1; i<n; i++)
        averageValue += (equationsSystem.getT(i-1)+equationsSystem.getT(i));
    averageValue *= 0.5q*h/l;

    // Calcula o fluxo de calor
    heatFlowLeft = -data.k/h*(equationsSystem.getT(1)-equationsSystem.getT(0));
    heatFlowRight = -data.k/h*(equationsSystem.getT(n-1)-equationsSystem.getT(n-2));
    if(problemType == Aleta)
        heatFlowLeft *= data.Ab;

    // Calcula Valor máximo
    maxValue = equationsSystem.getT(0);
    for(tInteger i=1; i<n; i++)
        maxValue = equationsSystem.getT(i)>maxValue ? equationsSystem.getT(i):maxValue;

}

void Diffusion1Dp::plotSolution(Functor1D &analyticalSolution)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "image.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<n; i++)
        file1<<static_cast<double>(boundaryLeft.x+i*h)<<"\t"<<
               static_cast<double>(equationsSystem.getT(i))<<std::endl;
    file1.close();

    // Solução Analítica
    std::ofstream file2(dat2_filename.c_str());
    tFloat h_10 = l/(10*n-1); // Aumenta o número de pontos em 10X
    for(tInteger i=0; i<10*n; i++){
        tFloat xp = boundaryLeft.x+i*h_10;
        file2<<static_cast<double>(xp)<<"\t"<<static_cast<double>(analyticalSolution(xp))<<std::endl;
    }
    file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n";

    if(problemType == ParedePlana)
        file3 << "set title \""<<STR_PAREDE_PLANA<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
                 "set xlabel 'x'\n"
                 "set ylabel 'Temperatura'\n";
    else if(problemType == Aleta)
        file3 << "set title \""<<STR_ALETA<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
                 "set xlabel 'x'\n"
                 "set ylabel 'Temperatura'\n";
    else
        file3 << "set title \""<<STR_QML<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
                 "set xlabel 'y'\n"
                 "set ylabel 'Velocidade'\n";
    file3 <<

             "plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "
             "'" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}

// Imprime a solução do problema
void Diffusion1Dp::printSolution(Functor1D &analyticalSolution)
{
    std::cout<<std::endl;
    if(problemType == ParedePlana)
        std::cout<<STR_PAREDE_PLANA;
    else if(problemType == Aleta)
        std::cout<<STR_ALETA;
    else
        std::cout<<STR_QML;
    std::cout<<std::endl;
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
    std::cout<<std::endl<<std::setw(5)<<std::right<<"p";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"Xp";
    std::cout<<std::setw(OUT_FLOAT_WIDTH+1)<<std::right<<"Yp (Numérica)";
    std::cout<<std::setw(OUT_FLOAT_WIDTH+1)<<std::right<<"Yp (Analítica)";
    std::cout<<std::setw(OUT_FLOAT_WIDTH)<<std::right<<"Erro";
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');

    std::cout.precision(OUT_FLOAT_PRECISION);
    std::cout<<std::scientific;

#ifdef QUAD_PRECISION
    char str[1000];
    for(tInteger i=0;i<n;i++){
        tFloat xp = boundaryLeft.x+i*h;
        std::cout<<std::endl;
        std::cout<<std::setw(5)<<i;
        quadmath_snprintf(str, 1000, Q_FORMAT,xp);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,equationsSystem.getT(i));
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,analyticalSolution(xp));
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,analyticalSolution(xp)-equationsSystem.getT(i));
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<str;
    }
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
#else
    for(tInteger i=0;i<n;i++){
        tFloat xp = boundaryLeft.x+i*h;
        std::cout<<std::endl;
        std::cout<<std::setw(5)<<i;
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<xp;
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<equationsSystem.getT(i);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<analyticalSolution(xp);
        std::cout<<std::setw(OUT_FLOAT_WIDTH)<<analyticalSolution(xp)-equationsSystem.getT(i);
    }
    std::cout<<std::endl<<std::setfill('-')<<std::setw(5+4*OUT_FLOAT_WIDTH)<<""<<std::setfill(' ');
#endif
}


// Imprime os resultados das variáveis secundárias
void Diffusion1Dp::printSecondaryResults(tFloat AS_average, tFloat AS_left, tFloat AS_right, tFloat AS_left_1, tFloat AS_right_1)
{
    std::cout<<std::endl;
    std::cout<<std::endl<<std::setfill('-')<<std::setw(30+3*(OUT_FLOAT_WIDTH+5))<<""<<std::setfill(' ');
    std::cout<<std::endl<<std::setw(30)<<"";
    std::cout<<std::setw(OUT_FLOAT_WIDTH+5+3)<<"Solução Numérica";
    std::cout<<std::setw(OUT_FLOAT_WIDTH+5+3)<<"Solução Analítica";
    std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<"Erro";
    std::cout<<std::endl<<std::setfill('-')<<std::setw(30+3*(OUT_FLOAT_WIDTH+5))<<""<<std::setfill(' ');

#ifdef QUAD_PRECISION
    char str[1000];
    if(problemType == ParedePlana){
        std::cout<<std::endl<<std::setw(30+1)<<"Temperatura Média";
        quadmath_snprintf(str, 1000, Q_FORMAT,averageValue);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_average);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_average - averageValue);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        std::cout<<std::endl<<std::setw(30)<<"Temperatura em X = X0";
        quadmath_snprintf(str, 1000, Q_FORMAT,equationsSystem.getT(0));
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left_1);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left_1 - equationsSystem.getT(0));
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        std::cout<<std::endl<<std::setw(30)<<"Temperatura em X = XL";
        quadmath_snprintf(str, 1000, Q_FORMAT,equationsSystem.getT(n-1));
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_right_1);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_right_1 - equationsSystem.getT(n-1));
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        std::cout<<std::endl<<std::setw(30)<<"Fluxo de Calor em X = X0";
        quadmath_snprintf(str, 1000, Q_FORMAT,heatFlowLeft);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left - heatFlowLeft);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        std::cout<<std::endl<<std::setw(30)<<"Fluxo de Calor em X = XL";
        quadmath_snprintf(str, 1000, Q_FORMAT,heatFlowRight);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_right);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_right - heatFlowRight);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
    }
    else if(problemType == Aleta){
        std::cout<<std::endl<<std::setw(30)<<"Fluxo de Calor";
        quadmath_snprintf(str, 1000, Q_FORMAT,heatFlowLeft);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left - heatFlowLeft);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
    }
    else{
        std::cout<<std::endl<<std::setw(30+1)<<"Velocidade Média";
        quadmath_snprintf(str, 1000, Q_FORMAT,averageValue);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_average);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_average - averageValue);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        std::cout<<std::endl<<std::setw(30)<<"Velocidade Máxima";
        quadmath_snprintf(str, 1000, Q_FORMAT,maxValue);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
        quadmath_snprintf(str, 1000, Q_FORMAT,AS_left - maxValue);
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<str;
    }

#else
    if(problemType == ParedePlana){
        std::cout<<std::endl<<std::setw(30+1)<<"Temperatura Média";
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<averageValue;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_average;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_average - averageValue;
        std::cout<<std::endl<<std::setw(30)<<"Fluxo de Calor em X = X0";
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<heatFlowLeft;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_left;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_left - heatFlowLeft;
        std::cout<<std::endl<<std::setw(30)<<"Fluxo de Calor em X = XL";
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<heatFlowRight;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_right;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_right - heatFlowRight;
    }
    else if(problemType == Aleta){
        std::cout<<std::endl<<std::setw(30)<<"Fluxo de Calor";
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<heatFlowLeft;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_left;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_left - heatFlowLeft;
    }
    else{
        std::cout<<std::endl<<std::setw(30+1)<<"Velocidade Média";
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<averageValue;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_average;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_average - averageValue;
        std::cout<<std::endl<<std::setw(30)<<"Velocidade Máxima";
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<maxValue;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_left;
        std::cout<<std::setw(OUT_FLOAT_WIDTH+5)<<AS_left - maxValue;
    }

#endif
    std::cout<<std::endl<<std::setfill('-')<<std::setw(30+3*(OUT_FLOAT_WIDTH+5))<<""<<std::setfill(' ');
}


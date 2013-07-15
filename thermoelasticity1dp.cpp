#include "thermoelasticity1dp.h"

Thermoelasticity1Dp::Thermoelasticity1Dp()
{
}

Thermoelasticity1Dp::Thermoelasticity1Dp(tInteger _nodes, Boundary1D _left, Boundary1D _right,
                                         ElasticityData _eData, DiffusionData _dData)
{
    boundaryLeft = _left;
    boundaryRight = _right;
    eData = _eData;
    n = _nodes;

    // Calcula L e h
    l = boundaryRight.x - boundaryLeft.x;
    h = l/(n - 1);

    thermo = new Diffusion1Dp(ParedePlana, _nodes, _left, _right, _dData);
}

void Thermoelasticity1Dp::solver(void)
{
    thermo->solver();

    equationsSystem.setMaxEquations(n);

    // DISCRETIZAÇÃO
    equationsSystem(1.0q, 0.0q, 0.0q, 0.0q); // Se Dirichlet

    // Aplica discretização de d2T/dx2 = S(x), com CDS-2
    for(tInteger i=1; i<n-1; i++)
        equationsSystem(2.0q, 1.0q, 1.0q, -0.5q*h*thermo->data.alpha*
                        (thermo->equationsSystem.getT(i+1) - thermo->equationsSystem.getT(i-1)));

    equationsSystem(1.0q, 0.0q, 0.0q, 0.0q); // Se Dirichlet

    // Resolve o sistema de equações lineares
    equationsSystem.solver();


    // Calculo da deformação
    strain = new tFloat[n];
    strain[0] = (-3.0q*equationsSystem.getT(0)+4.0q*equationsSystem.getT(1)-equationsSystem.getT(2))/(2.0q*h);
    for(tInteger i=1; i<n-1; i++)
        strain[i] = (equationsSystem.getT(i+1)-equationsSystem.getT(i-1))/(2.0q*h);
    strain[n-1] = (3.0q*equationsSystem.getT(n-1)-4.0q*equationsSystem.getT(n-2)+equationsSystem.getT(n-3))/(2.0q*h);

    // Calculo da tensão
    stress = new tFloat[n];
    for(tInteger i=0; i<n; i++){
        stress[i] = eData.E*(strain[i]-thermo->data.alpha*(thermo->equationsSystem.getT(i)-thermo->boundaryLeft.bcValue));
        std::cout<<"\n"<<static_cast<double>(stress[i]);
    }

}

void Thermoelasticity1Dp::plotSolution(Functor1D &analyticalSolution)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "thermoelasticity1dp_solution.png";
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

    file3 << "set title \""<<STR_THERMOELASTICITY<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set xlabel 'x'\n"
             "set ylabel 'Deslocamento'\n";

    file3 <<

             "plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "
             "'" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}


void Thermoelasticity1Dp::plotSecondarySolutions(Functor1D &AS_Strain, Functor1D &AS_Stress)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "thermoelasticity1dp_solution.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<n; i++)
        file1<<static_cast<double>(boundaryLeft.x+i*h)<<"\t"<<
               static_cast<double>(strain[i])<<std::endl;
    file1.close();

    // Solução Analítica
    std::ofstream file2(dat2_filename.c_str());
    tFloat h_10 = l/(10*n-1); // Aumenta o número de pontos em 10X
    for(tInteger i=0; i<10*n; i++){
        tFloat xp = boundaryLeft.x+i*h_10;
        file2<<static_cast<double>(xp)<<"\t"<<static_cast<double>(AS_Strain(xp))<<std::endl;
    }
    file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n";

    file3 << "set title \""<<STR_THERMOELASTICITY<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set xlabel 'x'\n"
             "set ylabel 'Deformação'\n";

    file3 <<

             "plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "
             "'" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());


    // Solução numérica
    std::ofstream file10(dat1_filename.c_str());
    for(tInteger i=0; i<n; i++)
        file10<<static_cast<double>(boundaryLeft.x+i*h)<<"\t"<<
               static_cast<double>(stress[i])<<std::endl;
    file10.close();



    // Solução Analítica
    std::ofstream file20(dat2_filename.c_str());
    //tFloat h_10 = l/(10*n-1); // Aumenta o número de pontos em 10X
    for(tInteger i=0; i<10*n; i++){
        tFloat xp = boundaryLeft.x+i*h_10;
        file20<<static_cast<double>(xp)<<"\t"<<static_cast<double>(AS_Stress(xp))<<std::endl;
    }
    file20.close();

    std::ofstream file30(cmd_filename.c_str());
    file30 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n";

    file30 << "set title \""<<STR_THERMOELASTICITY<<"\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set xlabel 'x'\n"
             "set ylabel 'Tensão'\n";

    file30 <<

             /*"plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "*/
             "plot '" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file30.close();

    //const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    //const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());

}

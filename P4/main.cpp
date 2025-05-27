#include <iostream>
#include <filesystem>
#include "bc_solver.hpp"
//#include "jacobi.hpp"
//#define USE_JACOBI

int main() {
   std::filesystem::path myPath=getenv("HOME");
   myPath/="Data/test03";

#ifdef USE_JACOBI
    Jacobi JJ;
    JJ.C=200;
    JJ.resize(12000);
    JJ.Init();
    JJ.Inv();
    try{
    JJ.Save(myPath/"C200N20.JACOBI");
    }catch(const std::exception& e)
    {
        std::cerr << "Jacobiのやろうなんかえらーしおったで. " << e.what() << std::endl;
    }
#else
    BC_Solver solver;
    solver.C=200.;
    solver.N=40;
    solver.solve();
    solver.save(myPath / "C200_N40.txt");
    solver.N = 20;
    solver.solve();
    solver.save(myPath / "C200_N20.txt");
#endif
    return 0;
}



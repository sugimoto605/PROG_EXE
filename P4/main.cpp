#include <iostream>
#include <filesystem>
#include "bc_solver.hpp"
#include "jacobi.hpp"
#include <map>
#define USE_JACOBI

int main() {
    std::filesystem::path myPath=getenv("HOME");
    myPath/="Data/test03";
#ifdef USE_JACOBI
    Jacobi JJ;
    JJ.C=0.2;
    JJ.resize(200);
    JJ.Init();
    JJ.filename = myPath / "C02_N200.dat";
    JJ.SaveTimes.insert(1000);
    JJ.SaveTimes.insert(2000);
    JJ.SaveTimes.insert(5000);
    JJ.SaveTimes.insert(10000);
    JJ.SaveTimes.insert(20000);
    JJ.SaveTimes.insert(50000);
    JJ.SaveTimes.insert(100000);
    if (!JJ.Solve(myPath/"C02_N200.log"))
        std::cout << "ERROR!! OOPS!! DID NOT CONVERGE!!" << std::endl;
    else
    {
        JJ.Save(myPath/"C02_N200.jac");
        JJ.Inv();        //逆行列で
        JJ.Save(myPath/"C02_N200.sol");     //こっちが正解
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



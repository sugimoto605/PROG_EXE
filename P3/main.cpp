#include <iostream>
#include "ode_euler.hpp"
int main() {
    std::cout << "Hello, World!" << std::endl;
    EulerSolver solver(1.0,0.0);
    solver.dx=0.02;
    int Nmax=1.0/solver.dx;
    for(int i=0;i<Nmax;i++){
        solver.step(); // 1ステップ進める
    }
    solver.print(); // 結果を表示
    return 0;
}
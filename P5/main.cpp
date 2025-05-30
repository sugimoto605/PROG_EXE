#include <iostream>
#include <iomanip>
#include <cmath>
#include "IntFunc.hpp"
#include "intde.hpp"
int main() {
    IntDE myIntegrator;
    double err,err_req=1e-6;
    size_t count=0;
    myIntegrator.Initialize(err_req); // 初期化
    auto v=myIntegrator.Integrate([&count](double x, void *param) {
        auto P=static_cast<size_t*>(param);
        count++;
        return 1.0 / std::sqrt(x * (1.0 - x)); 
    }, 0.0, 1.0, err);
    std::cout << "Integral value: " << std::setprecision(12) << v << ", Error: " << err << " requested: " << err_req << " by " << count << " steps"<< std::endl;
    // IntFunc OBJ;
    // std::cout << OBJ.Integrate([](double x) {
    //         return 1.0/std::sqrt(x*(1.0-x)); 
    //     }, 0.0, 1.0,1e-5) << std::endl;
    return 0;
}



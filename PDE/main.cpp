#include <iostream>
#include "CONV.hpp"
int main()
{
    size_t nx=100; // 空間分割数
    double dt=0.01; // 時間刻み幅
    CONV mySolver(nx,dt,"Data/pde100.data");
    return 0;
}
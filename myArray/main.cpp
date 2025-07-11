#include <iostream>
#include "../include/vtk_write.hpp"
int main()
{
    vtk_write Writer;
    Writer.ox << -3,-3, 0;
    Writer.dx << 0.1, 0.1, 0.;
    Writer.nx << 60, 60, 1;
    std::filesystem::path myHome = getenv("HOME");
    Writer.filename=(myHome/"Data")/"test05";
    Writer.filename/="LEGACY_2D.vtk";
    Writer.write_legacy_vtk2d([&Writer](int i, int j) -> double
    {
        double x=Writer.ox[0]+Writer.dx[0]*i;
        double y=Writer.ox[1]+Writer.dx[1]*j;
        return std::exp(-x*x-y*y);
    });
    return 0;
}
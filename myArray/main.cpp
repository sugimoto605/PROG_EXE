#include <iostream>
#include "../include/vtk_myIO.hpp"
int main()
{
    vtk_myIO Writer;
    Writer.ox << -3, -3, 0;
    Writer.dx << 0.01, 0.01, 0.;
    Writer.nx << 600, 600, 1;
    std::filesystem::path myHome = getenv("HOME");
    auto myFolder = (myHome / "Data") / "test05";
    auto FUNC = [&Writer](int i, int j) -> double
    {
        double x = Writer.ox[0] + Writer.dx[0] * i;
        double y = Writer.ox[1] + Writer.dx[1] * j;
        return std::exp(-x * x - y * y);
    };
    Writer.filename = myFolder / "LEGACY_2D.vtk";
    Writer.write_legacy_vti2d(FUNC);
    Writer.filename = myFolder / "XML_2D.vti";
    Writer.write_xml_vti2d(FUNC);
    return 0;
}
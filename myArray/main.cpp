// main.cpp
#include <iostream>
#include "../include/vtk_myIO.hpp"
#include "../include/pvdwrite.hpp"
int T2()
{
    std::filesystem::path myHome = getenv("HOME");
    auto myFolder = (myHome / "Data") / "test06";
    PVDwrite Indexer(myFolder/"test.pvd");
    vtk_myIO Writer;
    Writer.ox << -3, -3, 0;
    Writer.dx << 0.01, 0.01, 0.;
    Writer.nx << 600, 600, 1;
    double time=0.0, dt=0.01;
    auto FUNC = [&time,&Writer](int i, int j) -> double
    {
        double x = Writer.ox[0] + Writer.dx[0] * i + 0.5*std::cos(2*M_PI*time);
        double y = Writer.ox[1] + Writer.dx[1] * j + 0.5*std::sin(2*M_PI*time);
        double rate= 1.0 + 0.5 * std::sin(8 * M_PI * time);
        return std::exp(-x * x - y * y)*rate;
    };
    for(int nt=0; nt<100; nt++)
    {
        time = nt * dt;
        std::ostringstream oss;
        oss << "XML_2D_" << std::setw(5) << std::setfill('0') << nt << ".vti";
        Writer.filename = (myFolder / "Contents") / oss.str();
        Writer.write_xml_vti2d(FUNC);
        Indexer.Append(Writer.filename, time);
    }
    Indexer.Write();
    std::cout << "PVD file written: " << Writer.filename << std::endl;
    return 0;
}

int T1()
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
int main()
{
    T2();
    return 0;
}

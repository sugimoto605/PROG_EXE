#pragma once
#include <iostream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <fstream>
class CONV
{
    double dt,CFL;
    double time;
    class Point {
    public:
        double X;
        double U_VAL;
        double U_PRE;
    };
    std::vector<Point> Data;
    std::ofstream ofs;
public:
    CONV(const size_t nx, const double dt,const std::function<double(double)> U0) : dt(dt)
    {
        double xmax = 1.0;
        time = 0.0;            // 初期時間
        double dx = xmax / nx; // 空間刻み幅の計算
        CFL = dt / dx;     // CFL条件の計算
        Data.resize(nx + 1); // 空間分割数に基づいて初期化
        for (size_t i = 0; i < Data.size(); ++i)
            Data[i].U_PRE = Data[i].U_VAL = U0(Data[i].X=i*dx); // 初期条件の適用
        std::cout << "Starting CONV solver with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    ~CONV(){
        if (ofs.is_open()) ofs.close();
    }
    bool Write(const std::string &filename)
    {
        if (ofs) throw std::runtime_error("File stream already open.");
        std::filesystem::path datafile = getenv("HOME");
        datafile /= filename;
        std::filesystem::create_directories(datafile.parent_path());
        ofs.open(datafile);
        if (!ofs) throw std::ios_base::failure("Failed to open file: " + datafile.string());
        return Write();
    }
    bool Write()
    {
        try
        {
            ofs << std::endl << "# time= " << time << std::endl;
            for (auto &v : Data)
                ofs << v.X << " " << v.U_VAL << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error writing to file: " << e.what() << std::endl;
            return false;
        }
        return true;
    }
    void Step(const double U_W)
    {
        // PDEの解法ロジックをここに実装
        std::cout << "Solving PDE..." << std::endl;
    }
};
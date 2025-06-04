#pragma once
#include <iostream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <fstream>
class CONV
{
    double dt, CFL;
    double time;
    class Point
    {
    public:
        size_t id;
        double X;
        double U_VAL; // U_n
        double U_PRE; // U_{n-1}
    };
    std::vector<Point> Data;
    std::ofstream ofs;

public:
    CONV(const size_t& nx, const double& dt, const std::function<double(double)>& U0 = [](double x)
                                           { return 0.; }) : dt(dt)
    {
        double xmax = 10.0;
        time = 0.0;            // 初期時間
        double dx = xmax / nx; // 空間刻み幅の計算
        CFL = dt / dx;         // CFL条件の計算
        Data.resize(nx + 1);   // 空間分割数に基づいて初期化
        for (size_t i = 0; i < Data.size(); ++i) Data[i].id = i; // 各点のIDを設定
        for (auto& P:Data) P.U_PRE = P.U_VAL = U0(P.X = P.id * dx); // 初期条件の適用
        std::cout << "Starting CONV solver with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    bool Write(const std::string &filename)
    {
        if (ofs) ofs.close(); // 既存のファイルストリームを閉じる
        std::filesystem::path datafile = getenv("HOME");
        datafile /= filename;
        std::filesystem::create_directories(datafile.parent_path());
        ofs.open(datafile);
        if (!ofs)
            throw std::ios_base::failure("Failed to open file: " + datafile.string());
        return Write();
    }
    bool Write()
    {
        try
        {
            ofs << std::endl
                << "# time= " << time << " CFL= " << CFL << std::endl;

            for (auto &v : Data)
                ofs << v.X << " " << v.U_VAL << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error writing to file: " << e.what() << std::endl;
            return false;
        }
        std::cout << "time[" << time << "]written to file successfully." << std::endl;
        return true;
    }
    void Step(const double U_W=1.0)
    {
        time+=dt;
        Data[0].U_VAL = U_W;
        for(size_t i=1;i<Data.size();++i)
            Data[i].U_VAL=(1.-CFL)*Data[i].U_PRE + CFL*Data[i-1].U_PRE;
        for(auto& v:Data) v.U_PRE=v.U_VAL;
    }
    void Step_Implicit(const double U_W=1.0)
    {
        time+=dt;
        Data[0].U_VAL = U_W;
        for(size_t i=1;i<Data.size();++i)
            Data[i].U_VAL = (Data[i].U_PRE+CFL*Data[i-1].U_VAL)/(1.+CFL);
        for(auto& v:Data) v.U_PRE=v.U_VAL;
    }
    double get_time() const { return time; }
};
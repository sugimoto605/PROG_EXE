// pCONV.hpp    Solve advection equation using explicit method. Periodic problem
#pragma once
#include <iostream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <fstream>
#include "pvector.hpp"

class pCONV
{
    double xmax = 1.0;
    double dt, CFL;
    double time;
    int nt = 0; // 時間ステップのカウンタ
    class Point
    {
    public:
        Point* pre = nullptr;   // ポインタを追加して前のポイントへのリンクを作成
        double X;
        pvector<double> U;      // U_n
        Point() : U(3) {}       // デフォルトコンストラクタでサイズ3のpvectorを初期化
        auto &operator[](int index) { return U[index]; }
        const auto &operator[](int index) const { return U[index]; }
    };
    pvector<Point> Data;
    std::ofstream ofs;

public:
    pCONV(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                             { return 0.; }) : dt(dt)
    {
        double C=1.0; // 速度
        time = 0.0;             // 初期時間
        double dx = xmax / nx;  // 空間刻み幅の計算
        CFL = C* dt / dx;          // CFL条件の計算
        Data.resize(nx);        // 空間分割数に基づいて初期化
        for (size_t i = 0; auto &P : Data)
        {
            P.pre = &Data[i - 1];
            P[nt] = P[nt - 1] = P[nt - 2] = U0(P.X = i * dx); // 初期条件の適用
            i++;
        }
        std::cout << "Starting pCONV solver with nx=" << Data.size() << ", dt=" << dt << std::endl;
        // for(auto& P:Data)
        // {
        //     std::cout << std::setprecision(5) << std::fixed << std::showpos 
        //     << P.X << " : " << P[nt] << " " << (*P.pre)[nt] << std::endl;
        // }
    }
    bool Write(const std::string &filename)
    {
        if (ofs) ofs.close(); // 既存のファイルストリームを閉じる
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
            ofs << std::endl
                << "# nt= " << nt << " time= " << time << " CFL= " << CFL << std::endl;
            for (auto &v : Data) ofs << v.X << " " << v[nt] << std::endl;
            ofs << Data.front().X+xmax << " " << Data.front()[nt] << std::endl; // 最初のデータを追加
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error writing to file: " << e.what() << std::endl;
            return false;
        }
        std::cout << "time[" << time << "]written to file successfully." << std::endl;
        return true;
    }
    void StepLQ()   //  1次精度, explicit method
    {
        time = (++nt) * dt;
        for (auto &P : Data)
            P[nt] = (1. - CFL) * P[nt - 1] + CFL * (*P.pre)[nt - 1];
    }
    void StepHQ()
    {
        time = (++nt) * dt;
        for (auto &P : Data)
        {
            auto& Pn1=*P.pre;   //左隣
            auto& Pn2=*Pn1.pre; //左左隣
            // 時間1次, 空間2次 これはできる
            // ∂u/∂t = (u[i][nt]-u[i][nt-1])/dt
            // ∂u/∂x = (1.5*u[i][nt-1]-2.0*u[i-1][nt-1]+0.5*u[i-2][nt-1])/dx
            // (u[i][nt]-u[i][nt-1])/dt + C*(1.5*u[i][nt-1]-2.0*u[i-1][nt-1]+0.5*u[i-2][nt-1])/dx=0
            // u[i][nt]-u[i][nt-1] + CFL*(1.5*u[i][nt-1]-2.0*u[i-1][nt-1]+0.5*u[i-2][nt-1])=0
            // u[i][nt]=(1-1.5*CFL)*u[i][nt-1]+2*CFL*u[i-1][nt-1]-0.5*CFL*u[i-2][nt-1]
            P[nt] = (1.-1.5*CFL)*P[nt - 1] +2.*CFL*Pn1[nt-1]-0.5*CFL* Pn2[nt - 1];
        }
    }
    // 時間2次, 空間2次 V1 @[i][nt-1] これでは，うまくいかない
    // ∂u/∂t = (u[i][nt]-u[i][nt-2])/(2*dt)
    // ∂u/∂x = (1.5*u[i][nt-1]-2.0*u[i-1][nt-1]+0.5*u[i-2][nt-1])/dx
    // (u[i][nt]-u[i][nt-2])/2 + (dt*C/dx)(1.5*u[i][nt-1]-2.0*u[i-1][nt-1]+0.5*u[i-2][nt-1])=0
    // u[i][nt]=u[i][nt-2]-3*CFL*u[i][nt-1]+4*CFL*u[i-1][nt-1]-CFL*u[i-2][nt-1]
    // P[nt] = P[nt - 2] - 3. * CFL * P[nt - 1] + 4. * CFL * Pn1[nt - 1] - CFL * Pn2[nt - 1];
    // 時間2次, 空間2次  @[i][nt-2] これでもダメ
    // ∂u/∂t = (-1.5*u[i][nt-2]+2.0*u[i][nt-1]-0.5*u[i][nt])/dt
    // ∂u/∂x = (+1.5*u[i][nt-2]-2.0*u[i-1][nt-2]+0.5*u[i-2][nt-2])/dx
    // (-1.5*u[i][nt-2]+2.0*u[i][nt-1]-0.5*u[i][nt])/dt+C*(+1.5*u[i][nt-2]-2.0*u[i-1][nt-2]+0.5*u[i-2][nt-2])/dx=0
    // -3*u[i][nt-2]+4*u[i][nt-1]-u[i][nt]+CFL*(+3*u[i][nt-2]-4*u[i-1][nt-2]+u[i-2][nt-2])=0
    // u[i][nt]=4*u[i][nt-1]-3*u[i][nt-2]+3*CFL*u[i][nt-2]+CFL*u[i-2][nt-2]-4*CFL*u[i-1][nt-2]
    // P[nt]=4*P[nt-1]-3*P[nt-2]+3*CFL*P[nt-2]+CFL*Pn2[nt-2]-4*CFL*Pn1[nt-2];
    double get_time() const { return time; }
};
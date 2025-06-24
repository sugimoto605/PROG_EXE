// bSolver1DBase.hpp    Solve some equation using some method. Bounded region
#pragma once
#include <iostream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <fstream>
#include "../include/pvector.hpp"
template <typename T>
class bSolver1DBase
{
protected:
    double C = 1.0, K=1.0; // Cは速度、Kは拡散係数
    double CFL, DN; // CFL条件と拡散数
    double dt, dx, time, xmax = 1.0;
    int nt = 0; // 時間ステップのカウンタ
    class Point
    {
    public:
        Point *pre = nullptr; // ポインタを追加して前のポイントへのリンクを作成
        Point *next = nullptr;
        double X;
        pvector<T> U;     // U_n
        Point() : U(3) {} // デフォルトコンストラクタでサイズ3のpvectorを初期化
        auto &operator[](int index) { return U[index]; }
        const auto &operator[](int index) const { return U[index]; }
    };
    std::vector<Point> Data;
    std::ofstream ofs;
    double LBC_value = 0.0; // 下端の境界条件の値   境界条件がtに依存するケースは除外
    double RBC_value = 0.0; // 上端の境界条件の値   境界条件がtに依存するケースは除外
public:
    bSolver1DBase(const size_t &nx, const double &dt, const std::function<double(double)> &U0) : dt(dt)
    {
        time = 0.0;         // 初期時間
        dx = xmax / nx;     // 空間刻み幅の計算
        CFL = C * dt / dx;  // CFL条件の計算
        DN = K * dt /dx/dx; // 拡散数の計算
        Data.resize(nx+1);  // 空間分割数に基づいて初期化
        size_t i=0;
        for (auto &P : Data)
        {
            P.pre = &Data[i - 1];
            P.next = &Data[i + 1];
            P[nt] = P[nt - 1] = P[nt - 2] = U0(P.X = i * dx); // 初期条件の適用
            i++;
        }
        std::cout << "Booting bSolver1DBase with C=" << C << " K=" << K << " nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    double& LBC() { return LBC_value; } // 下端の境界条件
    double& RBC() { return RBC_value; } // 上端の境界条件
    virtual std::string what() const
    {
        return "nSolver1DBase: nx=" + std::to_string(Data.size()) + ", dt=" + std::to_string(dt) +
               ", CFL=" + std::to_string(CFL) + " , DN=" + std::to_string(DN)  + ", time=" + std::to_string(time);
    }
    bool Write(const std::string &filename)
    {
        if (ofs)
            ofs.close(); // 既存のファイルストリームを閉じる
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
            ofs << std::endl << "#" << what() << " nt= " << nt << std::endl;
            for (auto &v : Data)
                ofs << v.X << " " << v[nt] << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error writing to file: " << e.what() << std::endl;
            return false;
        }
        std::cout << "time[" << time << "]written to file successfully." << std::endl;
        return true;
    }
    virtual void Step() = 0;
    virtual void Initialize(void *parm = nullptr) {};
    double get_time() const { return time; }
};
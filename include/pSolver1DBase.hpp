// pSolver1DBase.hpp    Solve some equation using some method. Periodic problem
#pragma once
#include <iostream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <fstream>
#include "../include/pvector.hpp"

template <typename T>
class pSolver1DBase
{
    protected:
        double C = 1.0, CFL;
        double K=1.0, DN;
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
    pvector<Point> Data;
    std::ofstream ofs;
public:
    pSolver1DBase(const size_t &nx, const double &dt, const std::function<double(double)> &U0) : dt(dt)
    {
        time = 0.0;            // 初期時間
        dx = xmax / nx;        // 空間刻み幅の計算
        CFL = C * dt / dx;     // CFL条件の計算
        Data.resize(nx);       // 空間分割数に基づいて初期化
        for (size_t i = 0; auto &P : Data)
        {
            P.pre = &Data[i - 1];
            P.next= &Data[i + 1];
            P[nt] = P[nt - 1] = P[nt - 2] = U0(P.X = i * dx); // 初期条件の適用
            i++;
        }
        std::cout << "Booting pSolver1DBase with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    std::string what() const
    {
        return "nx=" + std::to_string(Data.size()) + ", dt=" + std::to_string(dt) +
               ", CFL=" + std::to_string(CFL) + ", DN=" + std::to_string(DN) + ", time=" + std::to_string(time);
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
            ofs << std::endl
                << "# nt= " << nt << " " << what() << std::endl;
            for (auto &v : Data)
                ofs << v.X << " " << v[nt] << std::endl;
            ofs << Data.front().X + xmax << " " << Data.front()[nt] << std::endl; // 最初のデータを追加
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error writing to file: " << e.what() << std::endl;
            return false;
        }
        std::cout << "time[" << time << "]written to file successfully." << std::endl;
        return true;
    }
    virtual void Step()=0;
    virtual void Initialize(void *parm=nullptr){};
    double get_time() const { return time; }
    double get_dx() const { return dx; }
};
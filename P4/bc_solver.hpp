#pragma once // ifdef...とかしなくても, 現代では，これで行けるらしいな
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class BC_Solver
{
public:
    double C = 1.0;
    int N = 100;
    double Xmin = 0.0;
    double Xmax = 1.0;
    double Y0 = 0.0; // 左端の境界条件
    double Y1 = 1.0; // 右端の境界条件
    std::vector<double> y;
    double dx;
    // 解く
    void solve()
    {
        dx = (Xmax - Xmin) / N;
        std::vector<double> b;
        b.resize(N - 1); // 方程式の数だけ確保
        y.resize(N + 1); // 境界条件で２個追加
        class ROW
        {
        public:
            double L, C, R; // 左列,中列,右列
        };
        std::vector<ROW> A(N - 1); // 要素数を最初から定義しても良い
        for (auto &row : A)
        {
            row.L = -1.;
            row.C = 2. - C * dx * dx;
            row.R = -1.;
        }
        // 境界条件
        for (auto &v : b)
            v = 0.0;
        b[0] = Y0;
        b[N - 2] = Y1;
        // 前進消去
        for (int i = 1; i < N - 1; ++i)
        {
            auto m = A[i].L / A[i - 1].C;
            A[i].C -= m * A[i - 1].R;
            b[i] -= m * b[i - 1];
        }
        // 後退代入
        y[N] = Y1;
        for (int i = N - 1; i > 0; --i)
            y[i] = (b[i - 1] - ((i < (N - 1)) ? A[i - 1].R * y[i + 1] : 0.0)) / A[i - 1].C;
        y[0] = Y0;
    }
    // ファイル出力
    void save(std::filesystem::path filename)
    {
        std::filesystem::create_directories(filename.parent_path());
        std::ofstream ofs(filename);
        for (int i = 0; i <= N; ++i)
        {
            ofs << i * dx << " " << y[i] << "\n";
        }
        ofs.close();
        std::cout << "Saved: " << filename << std::endl;
    }
};
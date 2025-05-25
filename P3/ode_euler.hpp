#pragma once
#include <iostream>
// y' = -y をオイラー法で解くクラス
class EulerSolver {
    double y;
    double x;
public:
    double dx;
    // コンストラクタ
    // y0: 初期値, t0: 開始時刻, h: 刻み幅
    EulerSolver(){
        std::cout << "Starting solver-A" << std::endl;
    }
    EulerSolver(double y0, double x0):y(y0),x(x0),dx(0.1)
    {
        //x=x0;y=y0;
        dx=0.1;
        std::cout << "Starting solver-B" << std::endl;
    }
    // 1ステップ進める
    void step() {
        y += dx * (-y);
        x += dx;
    }
    void print() {
        double y_honma=exp(-x);
        std::cout << "x=" << x << ", y=" << y << " error= " << y-y_honma <<  std::endl;
    }
    // y' = -y をオイラー法で解くクラス

//     // nステップ進めて結果を返す
//     std::vector<double> solve(int n) {
//         std::vector<double> ys;
//         ys.push_back(y_);
//         for (int i = 0; i < n; ++i) {
//             step();
//             ys.push_back(y_);
//         }
//         return ys;
//     }

//     // 現在のt, yを取得
//     double t() const { return t_; }
//     double y() const { return y_; }

// private:
//     double y_;
//     double t_;
//     double h_;
};
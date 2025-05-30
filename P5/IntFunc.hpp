//  int_func: 与えられた関数を積分
#pragma once // ifdef...とかしなくても, 現代では，これで行けるらしいな
#include <iostream>
class IntFunc
{
public:
    IntFunc(){};
    double Get(std::function<double(double)> F, double x)
    {
        return F(x);
    }
    double Integrate(std::function<double(double)> f, double a, double b, double eps = 1e-6) const {
        double sum = 0.0;
        int n = 100; // 初期分割数
        double h = (b - a) / n;
        while (true) {
            double new_sum = 0.0;
            for (int i = 0; i < n; ++i) {
                new_sum += f(a + i * h + h / 2.0); // 中点法
            }
            new_sum *= h;
            std::cout << "n= " << n << ", sum= " << new_sum << std::endl; // デバッグ用出力
            if (std::fabs(new_sum - sum) < eps) {
                return new_sum; // 精度が満たされたら終了
            }
            sum = new_sum;
            n *= 2; // 分割数を倍にして再計算
            if (n>1000000000) // もし分割数が大きすぎる場合は終了
                throw std::runtime_error("Too many iterations, consider increasing the initial n or decreasing eps.");
            h = (b - a) / n;
        }
    }
};

    // // 積分区間 [a, b]
    // double a = 0.0;
    // double b = 1.0;

    // // 積分の精度
    // double eps = 1e-6;

    // // 関数 f(x) を定義
    // virtual double f(double x) const = 0;

    // // 数値積分を行う
    // double integrate() const {
    //     double sum = 0.0;
    //     int n = 100; // 初期分割数
    //     double h = (b - a) / n;

    //     while (true) {
    //         double new_sum = 0.0;
    //         for (int i = 0; i < n; ++i) {
    //             new_sum += f(a + i * h + h / 2.0); // 中点法
    //         }
    //         new_sum *= h;

    //         if (std::fabs(new_sum - sum) < eps) {
    //             return new_sum; // 精度が満たされたら終了
    //         }

    //         sum = new_sum;
    //         n *= 2; // 分割数を倍にして再計算
    //         h = (b - a) / n;
    //     }
    // }

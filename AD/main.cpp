// main.cpp
#include "heat1D_explicit.hpp"
#include "heat1D_implicit.hpp"
#include <set>
int main()
{
    size_t nx = 160;         // 空間分割数 100(text), 40, 20, 8
    double dt = 0.25;       // 時間刻み幅 DN=100 に固定しとこ
    size_t ntmax = 5. / dt;  // 最大時間ステップ数
    size_t ntsave = 0.5 / dt; // データ保存間隔
    // std::set<int> tset;
    // tset.insert(round(0.005/dt));tset.insert(round(0.010/dt));tset.insert(round(0.020/dt));
    // tset.insert(round(0.050/dt));tset.insert(round(0.10/dt));tset.insert(round(0.20/dt));
    // tset.insert(round(0.50/dt));tset.insert(round(1/dt));tset.insert(round(2/dt));
    // tset.insert(round(5/dt));
    auto I_SIN = [](double x)
    { return std::pow(std::sin(2.0 * M_PI * x), 5); }; // 初期条件として正弦波^5を設定
    auto I_LEVEQUE = [](double x)
    {
        double b = 200;
        double xc = x - std::floor(x + 0.2) - 0.3;
        double y = exp(-b * xc * xc);
        if ((xc > 0.3) && (xc < 0.5))
            y += 1;
        return y;
    }; // Lebequeのテキストの初期関数
    auto I_DC = [](double x)
    { return (x > 0.8) ? 0 : std::max(2. * x - 0.6, 0.); };
    //--------------------------------------------------
    double K = 1e-2;
    std::string filename = "Data/test04/ADi_160_K1e-2.dat";
    // Heat1DExplicit mySolver(nx, dt);
    Heat1DImplicit mySolver(nx, dt); // 初期条件を指定してインスタンス化
    mySolver.LBC() = 1.0; // 下端の境界条件
    //--------------------------------------------------
    mySolver.Initialize(&K);
    mySolver.Write(filename);
    for (size_t nt = 0; nt < ntmax; nt++)
    {
        mySolver.Step(); // 時間ステップを実行
       if ((nt + 1) % ntsave == 0)
        // if (tset.contains(nt+1))
            mySolver.Write();
    }
    mySolver.Write(filename+".last");
    return 0;
}
// main.cpp
#include "heat1D_explicit.hpp"
#include "heat1D_implicit.hpp"
#include "heat1D_sakurai.hpp"
#include <set>
int main()
{
    size_t nx = 80;         // 空間分割数 100(text), 40, 20, 8
    double dt = 0.25;       // 時間刻み幅 DN=100
    size_t ntmax = 20. / dt;  // 最大時間ステップ数
    size_t ntsave = 1 / dt; // データ保存間隔
    // std::set<int> tset;
    // tset.insert(round(0.005/dt));tset.insert(round(0.010/dt));tset.insert(round(0.020/dt));
    // tset.insert(round(0.050/dt));tset.insert(round(0.10/dt));tset.insert(round(0.20/dt));
    // tset.insert(round(0.50/dt));tset.insert(round(1/dt));tset.insert(round(2/dt));
    // tset.insert(round(5/dt));
    auto I_SIN = [](double x)
    { return std::pow(std::sin(2.0 * M_PI * x), 5); }; // 初期条件として正弦波^5を設定
    auto I_X = [](double x)
    { return x; }; // 初期条件としてy=xを設定
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
    double C = 1;
    double K = 0.01;
    auto parm= std::pair<double, double>(C, K); // パラメータのペアを作成
    std::string filename = "Data/test04/ADBs_80_K0.01.dat";
    // Heat1DExplicit mySolver(nx, dt);
    // Heat1DImplicit mySolver(nx, dt);
    Heat1D_Sakurai mySolver(nx, dt);
    mySolver.LBC_func = [](double x){return 1.0;}; // 下端の境界条件
    mySolver.dLBC_func = [](double t){return 0.0;};// 下端の境界条件の微分
    mySolver.dRBC_func = [](double t){return 0.0;}; // 上端の境界条件の微分
    // mySolver.I_func = I_LEVEQUE; // 初期条件は u(x) = x
    //--------------------------------------------------
    mySolver.Initialize(&parm); // 初期化
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
// main.cpp
#include "AD_P1D.hpp"
int main()
{
    size_t nx = 400;          // 空間分割数 100(text), 40, 20, 8
    double dt = 0.02/2;         // 時間刻み幅 DN=100 に固定しとこ
    size_t ntmax = 1 / dt;   // 最大時間ステップ数
    size_t ntsave = 0.2 / dt; // データ保存間隔
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
    AD_P1D mySolver(nx, dt, I_LEVEQUE); // 初期条件を指定してインスタンス化
    double dx=mySolver.get_dx(); // 空間刻み幅を取得
    double K = dx*dx*1; // 拡散係数
    std::string filename = "Data/test04/ADiP_400_Kdx2_CFL4.dat";
    //--------------------------------------------------
    mySolver.Initialize(&K);
    mySolver.Write(filename);
    for (size_t nt = 0; nt < ntmax; nt++)
    {
        mySolver.Step(); // 時間ステップを実行
        if ((nt + 1) % ntsave == 0) mySolver.Write();
    }
    mySolver.Write(filename + ".last");
    return 0;
}
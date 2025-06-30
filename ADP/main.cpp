// main.cpp
#include "AD_P1D.hpp"
#include "AD_P1D_sakurai.hpp"
int main()
{
    size_t nx = 400;          // 空間分割数 100(text), 40, 20, 8
    double dt = 80./nx;      // 時間刻み幅 CFL=0.8 に固定しとこ
    size_t ntmax = 5 / dt;   // 最大時間ステップ数
    size_t ntsave = 1 / dt; // データ保存間隔
    std::pair<double, double> P = {1.0, 1.0}; // 速度と拡散係数をペアで渡す
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
    auto dI_LEVEQUE = [&P](double x)
    {
        //  u_t = -C u_x + K u_xx
        double b = 200;
        double xc = x - std::floor(x + 0.2) - 0.3;
        double dy = -2.*b*xc * exp(-b * xc * xc);
        double ddy = -2. * b * exp(-b * xc * xc) + 4.*b*b*xc*xc *exp(-b*xc*xc);
        return -P.first*dy+P.second*ddy;
    };
    auto I_DC = [](double x)
    { return (x > 0.8) ? 0 : std::max(2. * x - 0.6, 0.); };
    //--------------------------------------------------
    AD_P1D_Sakurai mySolver(nx, dt);    // 初期条件を指定してインスタンス化
    double dx=mySolver.get_dx();        // 空間刻み幅を取得
    P.second = 0.02*dx;                    // 拡散係数 0.25でよさげ. 0, 0.25, 0.5, 1くらいで    
    std::string filename = "Data/test04/ADs_400_Kdx002_CFL80.dat";
    //--------------------------------------------------
    mySolver.I_func = I_LEVEQUE;
    mySolver.dI_func = dI_LEVEQUE;
    mySolver.Initialize(&P);
    mySolver.Write(filename);
    for (size_t nt = 0; nt < ntmax; nt++)
    {
        mySolver.Step(); // 時間ステップを実行
        if ((nt + 1) % ntsave == 0) mySolver.Write();
    }
    mySolver.Write(filename + ".last");
    return 0;
}
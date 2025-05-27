#pragma once // ifdef...とかしなくても, 現代では，これで行けるらしいな
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

class Jacobi
{
private:
    Eigen::MatrixXd myMAT; // Coefficient matrix
    Eigen::VectorXd myVEC; // Right-hand side vector
    Eigen::VectorXd mySOL; // Solution vector
public:
    double C = 1.0;
    double Xmin = 0.0;
    double Xmax = 1.0;
    double Y0 = 0.0; // 左端の境界条件
    double Y1 = 1.0; // 右端の境界条件
    double dx;
    Jacobi()
    {
        std::cout << "Jacobi solver initialized." << std::endl;
    };
    void Init()
    {
        dx = (Xmax - Xmin) / (myVEC.size()+1);
        // Initialize the matrix and vector
        for (int i = 0; i < myMAT.rows(); ++i)
        {
            myMAT(i, i) = 2.0 - C * dx * dx; // Diagonal elements
            if (i > 0)
                myMAT(i, i - 1) = -1.0; // Lower diagonal
            if (i < myMAT.rows() - 1)
                myMAT(i, i + 1) = -1.0; // Upper diagonal
        }
        myVEC[0] = Y0;
        myVEC[myVEC.size()-1] = Y1;
    }
    void Inv()
    {
        auto SOL = myMAT.inverse() * myVEC;
        mySOL[0]=Y0;
        for(int i=0; i<SOL.size(); ++i)
        {
            mySOL[i+1]=SOL[i];
        }
        mySOL[mySOL.size()-1]=Y1;
    }
    // ファイル出力
    void Save(std::filesystem::path filename)
    {
        // このプログラムだと，filenameの途中に，未作成のフォルダーがあるとエラーが発生するね
        std::ofstream ofs(filename);
        for (int i = 0; i < mySOL.size(); ++i)
        {
            ofs << i * dx << " " << mySOL[i] << "\n";
        }
        ofs.close();
        std::cout << "Saved: " << filename << std::endl;
    }
    void Solve()
    {
        // Jacobi method implementation
        int n = myMAT.rows();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n);
        double tolerance = 1e-10;
        int max_iterations = 1000;

        for (int k = 0; k < max_iterations; ++k)
        {
            for (int i = 0; i < n; ++i)
            {
                double sum = 0.0;
                for (int j = 0; j < n; ++j)
                {
                    if (i != j)
                    {
                        sum += myMAT(i, j) * x(j);
                    }
                }
                x_new(i) = (myVEC(i) - sum) / myMAT(i, i);
            }
            if ((x_new - x).norm() < tolerance)
            {
                break;
            }
            x = x_new;
        }
    }
    void resize(size_t n)
    {
        myMAT.resize(n-1, n-1);
        myMAT.setZero();
        myVEC.resize(n-1);
        myVEC.setZero();
        mySOL.resize(n + 1);
        mySOL.setZero();
    }
};
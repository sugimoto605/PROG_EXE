// Jacobi method: This project implements the Jacobi method for solving linear equations
// Branch Main
#pragma once // ifdef...とかしなくても, 現代では，これで行けるらしいな
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <unordered_set>

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
    std::unordered_set<size_t> SaveTimes;
    std::filesystem::path filename;
     Jacobi()
    {
        SaveTimes.clear();
        filename = "FuckingTaco.data"; // Default filename
        std::cout << "Jacobi solver initialized." << std::endl;
    };
    void Init()
    {
        dx = (Xmax - Xmin) / (myVEC.size() + 1);
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
        myVEC[myVEC.size() - 1] = Y1;
        mySOL[0] = Y0;
        mySOL[mySOL.size() - 1] = Y1;
    }
    void Inv()
    {
        auto SOL = myMAT.inverse() * myVEC;
        for (int i = 0; i < SOL.size(); ++i)
        {
            mySOL[i + 1] = SOL[i];
        }
    }
    // ファイル出力
    void Save(std::filesystem::path filename="SAVEFILE_TACO.taco")
    {
        std::filesystem::create_directories(filename.parent_path());
        std::ofstream ofs(filename);
        for (int i = 0; i < mySOL.size(); ++i)
        {
            ofs << i * dx << " " << mySOL[i] << "\n";
        }
        ofs.close();
        std::cout << "Saved: " << filename << std::endl;
    }
    bool Solve(size_t max_iter = 10000)
    {
        size_t message_iter = 1000;
        double tol = 1e-4;
        Eigen::VectorXd x_old = myVEC;
        Eigen::VectorXd x_new = myVEC;
        x_old.setZero();
        x_new.setZero();
        // 初期化
        size_t iter = 0;
        double err = 1.0;

        std::filesystem::create_directories(filename.parent_path());
        std::ofstream ofs(filename);
        while (iter < max_iter)
        {
            iter++;
            for (int i = 0; i < myMAT.rows(); ++i)
            {
                x_new[i] = 0.0; // 初期化
                for (int j = 0; j < myMAT.rows(); ++j)
                    if (j != i) // 対角要素を除く行列積
                        x_new[i] += myMAT(i, j) * x_old[j];
                // D * x_new = b - (A-D)*x_old   対角行列だから割れば良い
                x_new[i] = (myVEC[i] - x_new[i]) / myMAT(i, i);
            }
            if ((err = (x_new - x_old).norm()) < tol)
                break;
            x_old = x_new;
            // if (!(iter%message_iter))
            //     std::cout << "Iteration " << iter << ": Error = " << err << std::endl;
            if (SaveTimes.contains(iter))
            {
                std::cout << "Iteration " << iter << ": Error = " << err ;
                std::cout << " Save to [" << filename.string() << "]" << std::endl;
                ofs << std::endl
                    << "# iter = " << iter << std::endl;
                ofs << 0 * dx << " " << Y0 << std::endl;
                for (size_t i = 0; i < x_new.size(); ++i)
                    ofs << (i + 1) * dx << " " << x_new[i] << std::endl;
                ofs << x_new.size() * dx << " " << Y1 << std::endl;
            }
            // if (iter % message_iter == 0)
            //     std::cout << "Iteration " << iter << ": Error = " << err << std::endl;
        }
        if (iter == max_iter)
            std::cout << "Timeout. Error= " << err << std::endl;
        else
            std::cout << "Converged in " << iter << " iterations." << std::endl;
        for (int i = 0; i < myMAT.rows(); ++i)
            mySOL[i + 1] = x_new[i];
        return iter != max_iter;
    }
    void resize(size_t n)
    {
        myMAT.resize(n - 1, n - 1);
        myMAT.setZero();
        myVEC.resize(n - 1);
        myVEC.setZero();
        mySOL.resize(n + 1);
        mySOL.setZero();
    }
    void PrintSaveTimes()
    {
        std::cout << "Save times: ";
        for (auto& time:SaveTimes)
            std::cout << time << " ";
        std::cout << std::endl;
    }
};
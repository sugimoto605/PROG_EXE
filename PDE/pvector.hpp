#pragma once
#include <iostream>
#include <vector>
#include <cmath>

class pvector
{
    std::vector<double> data;
    auto pindex(int index) // インデックスを正の値に変換するヘルパー関数
    const {
        int n = data.size();
        return ((index % n) + n) % n;
    }
public:
    pvector(size_t size = 0) : data(size) {}    
    auto &operator[](int index)
    {return data[pindex(index)];}
    const auto &operator[](int index) const
    {return data[pindex(index)];}
    auto size() const
    {return data.size();}
    void resize(size_t newSize)
    {data.resize(newSize);}
    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }
};
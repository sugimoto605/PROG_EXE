// myArray.hpp
#pragma once
#include <cassert>
#include <vector>
template <typename T>
class myArray
{
    std::vector<size_t> _sizes;
    std::vector<T> _data;
public:
    // 一般の関数
    const size_t size() // myArray.size()でデータの総数
    const {
        size_t sz = 1;
        for (auto k : _sizes)
            sz *= k;
        return sz;
    }
    double *data() // myArray.data()で1次元配列のデータにアクセスできる
    {
        return _data.data(); // std::vectorでは, data()で1次元配列のデータにアクセスできる
    }
    const double *data() // myArray.data()で1次元配列のデータにアクセスできる
    const {
        return _data.data(); // std::vectorでは, data()で1次元配列のデータにアクセスできる
    }
    double &operator[](const int k) // myArray[k]がk番目のデータ
    {
        return _data[k];
    }
    const double &operator[](const int k) // myArray[k]がk番目のデータ
    const {
        return _data[k];
    }
    auto begin() { return _data.begin(); }
    auto end() { return _data.end(); }
    const auto begin() const { return _data.begin(); }
    const auto end() const { return _data.end(); }
    auto &front() { return _data.front(); }
    auto &back() { return _data.back(); }
    const auto &front() const { return _data.front(); }
    const auto &back() const { return _data.back(); }
    // 1D配列
    myArray(const size_t n0) {resize(n0);} // コンストラクタ0:(N)を呼び出す
    void resize(const size_t n0)
    {
        _sizes.resize(1); // 1次元配列にする
        _sizes[0] = n0;
        _data.resize(size());
    }
    const size_t serial(const int i0)
    const {
        assert(i0 >= 0);
        assert(i0 < _sizes[0]);
        return i0;
    }
    double &operator()(const int i0) // オペレータ0:(i)を呼び出す
    {
        assert(_sizes.size()==1);
        return _data[serial(i0)];
    }
    const double &operator()(const int i0) // オペレータ0:(i)を呼び出す
    const {
        assert(_sizes.size() == 1);
        return _data[serial(i0)];
    }
    // 2D配列
    myArray(const size_t n0, const size_t n1) // コンストラクタ1:(N,M)を呼び出す
    {
        _sizes.resize(2); // 2次元配列にする
        _sizes[0] = n0;
        _sizes[1] = n1;
        _data.resize(size());
    }
    const size_t serial(const int i0, const int i1)
    const {
        assert(i1 >= 0);
        assert(i1 < _sizes[1]);
        return serial(i0) + _sizes[0] * i1;
    }
    double &operator()(const int i0, const int i1) // オペレータ1:(i,j)を呼び出す
    {
        assert(_sizes.size() == 2);
        return _data[serial(i0, i1)];
    }
    double &operator()(const int i0, const int i1) // オペレータ1:(i,j)を呼び出す
    const {
        assert(_sizes.size() == 2);
        return _data[serial(i0, i1)];
    }
    // 3次元配列
    myArray(const int n0, const int n1, const int n2) // コンストラクタ2:(N,M,O)を呼び出す
    {
        _sizes.resize(3); // 3次元配列にする
        _sizes[0] = n0;
        _sizes[1] = n1;
        _sizes[2] = n2;
        _data.resize(size());
    }
    const size_t serial(const int i0, const int i1, const int i2)
    const {
        assert(i2 >= 0);
        assert(i2 < _sizes[2]);
        return serial(i0, i1) + _sizes[0] * _sizes[1] * i2;
    }
    double &operator()(const int i0, const int i1, const int i2) // オペレータ2:(i,j,k)を呼び出す
    {
        assert(_sizes.size() == 3);
        return _data[serial(i0, i1, i2)];
    }
    const double &operator()(const int i0, const int i1, const int i2) // オペレータ2:(i,j,k)を呼び出す
    const {
        assert(_sizes.size() == 3);
        return _data[serial(i0, i1, i2)];
    }
};
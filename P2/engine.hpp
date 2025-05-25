#pragma once
#include <iostream>
class ENGINE{
public:
    int i;      //回転数
    std::string name;//エンジン名
    ENGINE()
    {
        i = 0;
        name = "noname-engine";
        std::cout << name << " created" << std::endl;
    };
    void print()
    {
        std::cout << name << " rotation=" << i << std::endl;
    };
};

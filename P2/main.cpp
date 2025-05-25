#include <iostream>
#include "engine.hpp"   //仕様書を読み込め
int main() {
    ENGINE E1,E2;
    E1.name="engine-1";
    E2.name="engine-2";
    E1.i=1000;
    E2.i=2000;
    E1.print();
    E2.print();
    return 0;
}
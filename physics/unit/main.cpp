#include <iostream>

#include "ScalarUnit.hpp"

int main (int argc, char *argv[]) {


Area a = 1.5_km2;

AccelerationComponent A = 6_kmps2;
Mass M = 5_mg;

//ForceComponent F;
TorqueComponent F = M*A;


::std::cout << F << ::std::endl;


}

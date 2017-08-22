//***************************************************************************************
//  This program is unit test of PS::F64vec class functions.
//***************************************************************************************

#include <particle_simulator.hpp>

int main(int argc, char* argv[]) {

    //--- test for Vecter3<>
    //------ constructor (x, y, z)
    PS::F64vec v;
    PS::F64vec v0 = PS::F64vec(2.0);
    PS::F64vec v1 = PS::F64vec(1.0, 0.0, -3.0);

    std::cout << " PS::F64vec constructor" << std::endl
              << "   v  =" << v.x  << ", " << v.y  << ", " << v.z  << std::endl
              << "   v0 =" << v0.x << ", " << v0.y << ", " << v0.z << std::endl
              << "   v1 =" << v1.x << ", " << v1.y << ", " << v1.z << std::endl << std::endl;

    //------ assign
    v = v1;
    std::cout << " PS::F64vec assign v = v1" << std::endl
              << "   v  =" << v.x << ", " << v.y << ", " << v.z << std::endl << std::endl;

    //------ operator +,-
    v = v0 + v1;
    std::cout << " PS::F64vec operator +" << std::endl
              << "   v0 + v1  =" << v.x << ", " << v.y << ", " << v.z << std::endl << std::endl;

    v = v0 - v1;
    std::cout << " PS::F64vec operator -" << std::endl
              << "   v0 - v1  =" << v.x << ", " << v.y << ", " << v.z << std::endl << std::endl;

    v -= v1;
    std::cout << " PS::F64vec operator -=" << std::endl
              << "   v -= v1  :" << v.x << ", " << v.y << ", " << v.z << std::endl << std::endl;

    v += v1;
    std::cout << " PS::F64vec operator +=" << std::endl
              << "   v += v1  :" << v.x << ", " << v.y << ", " << v.z << std::endl << std::endl;

    //------ vec*scalar
    v = v1*(-10.0);
    std::cout << " PS::F64vec vec*scalar" << std::endl
              << "   v1*(-10.0) =" << v.x << ", " << v.y << ", " << v.z << std::endl << std::endl;

    v *= -2.0;
    std::cout << " PS::F64vec vec*=scalar" << std::endl
              << "   v *= (-2.0) :" << v.x << ", " << v.y << ", " << v.z << std::endl << std::endl;

    //------ inner product
    PS::F64 s = v0*v1;
    std::cout << " PS::F64vec inner product (*)" << std::endl
              << "   v0*v1 =" << s << std::endl << std::endl;

    //------ cross product
    v0 = PS::F64vec(1.0, 0.0, 0.0);
    v1 = PS::F64vec(0.0, 1.0, 0.0);
    v = v0^v1;
    std::cout << " PS::F64vec cross product (^)" << std::endl
              << "   v0    =" << v0.x << ", " << v0.y << ", " << v0.z << std::endl
              << "   v1    =" << v1.x << ", " << v1.y << ", " << v1.z << std::endl
              << "   v0^v1 =" << v.x  << ", " << v.y  << ", " << v.z  << std::endl << std::endl;



    //--- test for MatrixSym3<>
    //------ constructor (xx, yy, zz, xy, xz, yz)
    PS::F64mat m;
    PS::F64mat m0 = PS::F64mat(1.0);
    PS::F64mat m1 = PS::F64mat(3.0, 2.0, 1.0, 0.0, -1.0, -2.0);

    //------ operator +,-

    //------ trace (xx + yy + zz)
}

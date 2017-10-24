//***************************************************************************************
//  This program is unit test of PS::F64vec class functions.
//***************************************************************************************

#include <particle_simulator.hpp>
#include <vec_ext.hpp>

int main(int argc, char* argv[]) {

    //--- test for Vecter3<>
    //------ constructor (x, y, z)
    PS::F64vec v;
    PS::F64vec v0 = PS::F64vec(2.0);
    PS::F64vec v1 = PS::F64vec(1.0, 0.0, -3.0);

    std::cout << " PS::F64vec constructor" << "\n"
              << "   v  = " << v  << "\n"
              << "   v0 = " << v0 << "\n"
              << "   v1 = " << v1 << "\n" << std::endl;

    //------ assign
    v = v1;
    std::cout << " PS::F64vec assign v = v1" << "\n"
              << "   v = " << v << "\n" << std::endl;

    //------ operator +,-
    v = v0 + v1;
    std::cout << " PS::F64vec operator +" << "\n"
              << "   v0 + v1 = " << v << "\n" << std::endl;

    v = v0 - v1;
    std::cout << " PS::F64vec operator -" << "\n"
              << "   v0 - v1 = " << v << "\n" << std::endl;

    v -= v1;
    std::cout << " PS::F64vec operator -=" << "\n"
              << "   v -= v1 : " << v << "\n" << std::endl;

    v += v1;
    std::cout << " PS::F64vec operator +=" << "\n"
              << "   v += v1 : " << v << "\n" << std::endl;

    //------ vec*scalar
    v = v1*(-10.0);
    std::cout << " PS::F64vec vec*scalar" << "\n"
              << "   v1*(-10.0) = " << v << "\n" << std::endl;

    v *= -2.0;
    std::cout << " PS::F64vec vec*=scalar" << "\n"
              << "   v *= (-2.0) : " << v << "\n" << std::endl;

    //------ inner product
    PS::F64 s = v0*v1;
    std::cout << " PS::F64vec inner product (*)" << "\n"
              << "   v0*v1 = " << s << "\n" << std::endl;

    s = VEC_EXT::dot(v0, v1);
    std::cout << " syntactic sugar of inner product." << "\n"
              << "   VEC_EXT::dot(v0, v1) = " << s << "\n" << std::endl;

    //------ cross product
    v0 = PS::F64vec(1.0, 0.0, 0.0);
    v1 = PS::F64vec(0.0, 1.0, 0.0);
    v = v0^v1;
    std::cout << " PS::F64vec cross product (^)" << "\n"
              << "   v0    = " << v0 << "\n"
              << "   v1    = " << v1 << "\n"
              << "   v0^v1 = " << v  << "\n" << std::endl;

    v = VEC_EXT::cross(v0, v1);
    std::cout << " syntactic sugar of cross product." << "\n"
              << "   VEC_EXT::cross(v0, v1) = " << v << "\n" << std::endl;

    //------ rotation
    const PS::F64 pi  = 3.141592653589793;
    PS::F64       rad;

    rad = 1.0/4.0*pi;
    v0 = PS::F64vec(0.0, 1.0, 0.0);
    v  = VEC_EXT::rot_x(v0, rad);
    std::cout << " rotation on axis X." << "\n"
              << "   v0  = " << v0  << "\n"
              << "   rad = " << rad << "\n"
              << "   VEC_EXT::rot_x(v0, rad) = " << v << "\n" << std::endl;

    v0 = PS::F64vec(1.0, 0.0, 0.0);
    v  = VEC_EXT::rot_y(v0, rad);
    std::cout << " rotation on axis Y." << "\n"
              << "   v0  = " << v0  << "\n"
              << "   rad = " << rad << "\n"
              << "   VEC_EXT::rot_y(v0, rad) = " << v << "\n" << std::endl;

    v0 = PS::F64vec(1.0, 0.0, 0.0);
    v  = VEC_EXT::rot_z(v0, rad);
    std::cout << " rotation on axis Z." << "\n"
              << "   v0  = " << v0  << "\n"
              << "   rad = " << rad << "\n"
              << "   VEC_EXT::rot_z(v0, rad) = " << v << "\n" << std::endl;

    //------ vector norm
    std::cout << " vector norm." << "\n"
              << "   v = " << v << "\n"
              << "   VEC_EXT::norm(v) = " << VEC_EXT::norm(v) << "\n" << std::endl;

    //--- test for MatrixSym3<>
    //------ constructor (xx, yy, zz, xy, xz, yz)
    //PS::F64mat m;
    //PS::F64mat m0 = PS::F64mat(1.0);
    //PS::F64mat m1 = PS::F64mat(3.0, 2.0, 1.0, 0.0, -1.0, -2.0);

    //------ operator +,-

    //------ trace (xx + yy + zz)
}

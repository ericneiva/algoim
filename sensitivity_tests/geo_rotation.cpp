// Examples to demonstrate Algoim's methods for computing high-order accurate quadrature schemes
// for implicitly defined domains in hyperrectangles. The file contains a single main() routine;
// compile it as you would for any .cpp file with a main() entry point.

#include <fstream>
#include <algoim_quad.hpp>

template<int N>
struct Ellipsoid
{
    Algoim::Real ang = 0.0;

    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        T s =  x(0)*cos(ang) + x(1)*sin(ang);
        T t = -x(0)*sin(ang) + x(1)*cos(ang);
        if (N == 2)
            return s*s + 4.0*t*t - 1.0;
        else
            return s*s + 4.0*t*t + 9.0*x(2)*x(2) - 1.0;
    }

    template<typename T>
    blitz::TinyVector<T,N> grad(const blitz::TinyVector<T,N>& x) const
    {
        T a = T(cos(ang)); // T is not generally of type Algoim::Real
        T b = T(sin(ang)); // T is not generally of type Algoim::Real
        T ds = ( 2.0*(a*a+4.0*b*b)*x(0)-6.0*a*b*x(1) );
        T dt = ( 2.0*(4.0*a*a+b*b)*x(1)-6.0*a*b*x(0) );
        if (N == 2)
            return blitz::TinyVector<T,N>(ds, dt);
        else
            return blitz::TinyVector<T,N>(ds, dt, 18.0*x(2));
    }
};

template<int N>
Ellipsoid<N> ellipsoidGen(Algoim::Real ang){
    Ellipsoid<N> e;
    e.ang = ang;
    return e;
};

// Surface area of a 3D ellipsoid, computed via the cells of a Cartesian grid
Algoim::Real ellipsoid_surface_area(int n, int order, Algoim::Real l, std::vector<Algoim::Real> o, Algoim::Real t)
{
    // std::cout << "Surface area of a 3D ellipsoid, computed via the cells of a " << n << " by " << n << " by " << n << " Cartesian grid:\n";
    Algoim::Real dx = l / n;
    Ellipsoid<3> phi = ellipsoidGen<3>(t);
    Algoim::Real volume = 0.0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k)
    {
        blitz::TinyVector<Algoim::Real,3> xmin = {o[0] + i*dx, o[1] + j*dx, o[2] + k*dx};
        blitz::TinyVector<Algoim::Real,3> xmax = {o[0] + i*dx + dx, o[1] + j*dx + dx, o[2] + k*dx + dx};
        volume += Algoim::quadGen<3>(phi, Algoim::BoundingBox<Algoim::Real,3>(xmin, xmax), 3, -1, order).sumWeights();
    }
    return volume;
};

// Volume of a 3D ellipsoid, computed via the cells of a Cartesian grid
Algoim::Real ellipsoid_volume(int n, int order, Algoim::Real l, std::vector<Algoim::Real> o, Algoim::Real t)
{
    // std::cout << "Volume of a 3D ellipsoid, computed via the cells of a " << n << " by " << n << " by " << n << " Cartesian grid:\n";
    Algoim::Real dx = l / n;
    Ellipsoid<3> phi = ellipsoidGen<3>(t);
    Algoim::Real volume = 0.0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k)
    {
        blitz::TinyVector<Algoim::Real,3> xmin = {o[0] + i*dx, o[1] + j*dx, o[2] + k*dx};
        blitz::TinyVector<Algoim::Real,3> xmax = {o[0] + i*dx + dx, o[1] + j*dx + dx, o[2] + k*dx + dx};
        volume += Algoim::quadGen<3>(phi, Algoim::BoundingBox<Algoim::Real,3>(xmin, xmax), -1, -1, order).sumWeights();
    }
    return volume;
};

int main(int argc, char* argv[])
{
    // std::cout << "Algoim Examples - High-order quadrature algorithms for implicitly defined domains\n\n";
    // std::cout << "h-convergence test on ellipse (2D) and ellipsoid (3D)\n\n";
    std::cout << std::fixed << std::setprecision(64);

    std::vector< std::vector<Algoim::Real> > err2D;
    std::vector< std::vector<Algoim::Real> > err3D;

    std::vector<Algoim::Real> coord;

    std::vector<Algoim::Real> c_err2D;
    std::vector<Algoim::Real> c_err3D;

    Algoim::Real exa3D = (qd_real::_2pi/qd_real(9.0));
    // Computed using ArbNumerics.jl. Strings must be used to assign to a quad-double.
    Algoim::Real sur3D = qd_real("4.4008095646649703416002003892297059434836743233771458003566869");

    Algoim::Real l = 2.2;
    Algoim::Real t = 0.0;
    Algoim::Real dt = (qd_real::_2pi/qd_real(100.0));

    std::vector<Algoim::Real> o = { -1.1, -1.1, -1.1 };

    Algoim::Real err = 0;

    int n = 32;
    int nods = 3;

    std::vector<int> ods = { 4, 7, 10 };

    for (int j = 0; j < nods; ++j)
    {

        int order = ods[j];
        int count = 0;
        // std::cout << "Solving for order " << order << "...\n";

        while (count < 101)
        {

            err = abs(sur3D-ellipsoid_surface_area(n,order,l,o,t));
            c_err2D.push_back(err);

            err = abs(exa3D-ellipsoid_volume(n,order,l,o,t));
            c_err3D.push_back(err);

            if (j == 0)
                coord.push_back(t);
            
            t += dt;
            ++count;

        }

        err2D.push_back(c_err2D);
        err3D.push_back(c_err3D);

        c_err2D.clear();
        c_err3D.clear();

        t = 0.0;

    }

    for (int i = 0; i < coord.size(); ++i)
    {
        std::cout << coord[i] << " ";
        for (int j = 0; j < nods; ++j)
        {
            std::cout << err2D[j][i] << " " << err3D[j][i] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
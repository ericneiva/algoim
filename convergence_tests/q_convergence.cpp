// Examples to demonstrate Algoim's methods for computing high-order accurate quadrature schemes
// for implicitly defined domains in hyperrectangles. The file contains a single main() routine;
// compile it as you would for any .cpp file with a main() entry point.

#include <fstream>
#include <algoim_quad.hpp>

// Define PI
#if defined(ALGOIM_QDREAL)
inline static qd_real pi() { return qd_real::_pi; }
#elif defined(ALGOIM_DDREAL)
inline static dd_real pi() { return dd_real::_pi; }
#else
inline static double pi() { return std::atan(1)*4; }
#endif

template<int N>
struct Ellipsoid
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return x(0)*x(0) + 4.0*x(1)*x(1) - 1.0;
        else
            return x(0)*x(0) + 4.0*x(1)*x(1) + 9.0*x(2)*x(2) - 1.0;
    }

    template<typename T>
    blitz::TinyVector<T,N> grad(const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return blitz::TinyVector<T,N>(2.0*x(0), 8.0*x(1));
        else
            return blitz::TinyVector<T,N>(2.0*x(0), 8.0*x(1), 18.0*x(2));
    }
};

// Area of a 2D ellipse, computed via the cells of a Cartesian grid
Algoim::Real ellipse_area(int n, int order)
{
    // std::cout << "Area of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
    Algoim::Real dx = 2.2 / n;
    Ellipsoid<2> phi;
    Algoim::Real area = 0.0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
    {
        blitz::TinyVector<Algoim::Real,2> xmin = {-1.1 + i*dx, -1.1 + j*dx};
        blitz::TinyVector<Algoim::Real,2> xmax = {-1.1 + i*dx + dx, -1.1 + j*dx + dx};
        area += Algoim::quadGen<2>(phi, Algoim::BoundingBox<Algoim::Real,2>(xmin, xmax), -1, -1, order).sumWeights();
    }
    return area;
};

// Volume of a 3D ellipsoid, computed via the cells of a Cartesian grid
Algoim::Real ellipsoid_volume(int n, int order)
{
    // std::cout << "Volume of a 3D ellipsoid, computed via the cells of a " << n << " by " << n << " by " << n << " Cartesian grid:\n";
    Algoim::Real dx = 2.2 / n;
    Ellipsoid<3> phi;
    Algoim::Real volume = 0.0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k)
    {
        blitz::TinyVector<Algoim::Real,3> xmin = {-1.1 + i*dx, -1.1 + j*dx, -1.1 + k*dx};
        blitz::TinyVector<Algoim::Real,3> xmax = {-1.1 + i*dx + dx, -1.1 + j*dx + dx, -1.1 + k*dx + dx};
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

    std::vector<Algoim::Real> c_err2D;
    std::vector<Algoim::Real> c_err3D;

    Algoim::Real exa2D = qd_real::_pi2;
    Algoim::Real exa3D = (qd_real::_2pi/qd_real(9.0));

    Algoim::Real err = 0;

    int ods = 10;
    int ns = 3;
    std::vector<int> n{ 4, 32, 128 };

    for (int j = 0; j < ns; ++j)
    {

        for (int i = 1; i <= ods; ++i)
        {

            // std::cout << "Solving for order " << i << "...\n";

            err = abs(exa2D-ellipse_area(n[j],i));
            c_err2D.push_back(err);

            err = abs(exa3D-ellipsoid_volume(n[j],i));
            c_err3D.push_back(err);

        }

        err2D.push_back(c_err2D);
        err3D.push_back(c_err3D);

        c_err2D.clear();
        c_err3D.clear();

    }

    for (int i = 0; i < ods; ++i)
    {
        std::cout << i+1 << " ";
        for (int j = 0; j < ns; ++j)
        {
            std::cout << err2D[j][i] << " " << err3D[j][i] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
// Examples to demonstrate Algoim's methods for computing high-order accurate quadrature schemes
// for implicitly defined domains in hyperrectangles. The file contains a single main() routine;
// compile it as you would for any .cpp file with a main() entry point.

#include <fstream>
#include <algoim_quad.hpp>

const Algoim::Real L = 4.25; 

template<int N>
struct Gyroid
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return cos(x(0))*sin(x(1));
        else
            return cos(x(0))*sin(x(1)) + cos(x(1))*sin(x(2)) + cos(x(2))*sin(x(0));
    }

    template<typename T>
    blitz::TinyVector<T,N> grad(const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return blitz::TinyVector<T,N>( -sin(x(0))*sin(x(1)), cos(x(0))*cos(x(1)) );
        else
            return blitz::TinyVector<T,N>( -sin(x(0))*sin(x(1)) + cos(x(2))*cos(x(0)), \
                                            cos(x(0))*cos(x(1)) - sin(x(1))*sin(x(2)), \
                                            cos(x(1))*cos(x(2)) - sin(x(2))*sin(x(0)) );
    }
};

template<int N>
struct Plane
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        return x(0) - x(1);
    }

    template<typename T>
    blitz::TinyVector<T,N> grad(const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return blitz::TinyVector<T,N>( T(1.0), T(-1.0) );
        else
            return blitz::TinyVector<T,N>( T(1.0), T(-1.0), T(0.0) );
    }
};

template<int N>
struct ExampleFunction
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return log( (1/(L*L)) * ( x(0)*x(0) + x(1)*x(1) ) + 3.0/8.0 );
        else
            return log( (1/(L*L)) * ( x(0)*x(0) + x(1)*x(1) + x(2)*x(2) ) + 3.0/8.0 );
    }
};

template<int N>
struct QuadraticFunction
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return x(0)*x(0) + x(1)*x(1);
        else
            return x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
    }
};

template<int N>
struct QuarticFunction
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return x(0)*x(0)*x(0)*x(0) + x(1)*x(1)*x(1)*x(1);
        else
            return x(0)*x(0)*x(0)*x(0) + x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2);
    }
};

Algoim::Real integrate(int n, int order)
{
    Algoim::Real dx = 2*L / n;
    Gyroid<3> phi;
    ExampleFunction<3> fun;
    Algoim::Real result = 0.0;
    Algoim::QuadratureRule<3> q;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n/2; ++k)
    {
        blitz::TinyVector<Algoim::Real,3> xmin = {-L + i*dx, -L + j*dx, -0.5*L + k*dx};
        blitz::TinyVector<Algoim::Real,3> xmax = {-L + i*dx + dx, -L + j*dx + dx, -0.5*L + k*dx + dx};
        q = Algoim::quadGen<3>(phi, Algoim::BoundingBox<Algoim::Real,3>(xmin, xmax), 3, -1, order);
        result += q(fun);
    }
    return result;
};

Algoim::Real integrate_exact(int n, int order)
{
    Algoim::Real dx = 2*L / n;
    Plane<3> phi;
/*     QuadraticFunction<3> fun; */
    QuarticFunction<3> fun;
    Algoim::Real result = 0.0;
    Algoim::QuadratureRule<3> q;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n/2; ++k)
    {
        blitz::TinyVector<Algoim::Real,3> xmin = {-L + i*dx, -L + j*dx, -0.5*L + k*dx};
        blitz::TinyVector<Algoim::Real,3> xmax = {-L + i*dx + dx, -L + j*dx + dx, -0.5*L + k*dx + dx};
        q = Algoim::quadGen<3>(phi, Algoim::BoundingBox<Algoim::Real,3>(xmin, xmax), 3, -1, order);
        result += q(fun);
    }
    return result;
};

int main(int argc, char* argv[])
{
    std::cout << std::fixed << std::setprecision(64);

    std::vector< std::vector<Algoim::Real> > err_array;
    std::vector<Algoim::Real> c_err_array;
/*     Algoim::Real exact = qd_real("6.897665194490618059924850963768989519102402631696"); Section 4.2 */
/*     Algoim::Real exact = qd_real("692.08904849392541936112799180498941339971229318373875"); // Quadratic */
    Algoim::Real exact = qd_real("6875.4721411318403379657058935876917037427668125972046456831796882454"); // Quartic
    Algoim::Real err = 0;

    int pts = 5;
    int nods = 5;
    std::vector<int> ods{ 1, 2, 4, 6, 10 };

    for (int j = 0; j < nods; ++j)
    {

        int order = ods[j];

        for (int i = 0; i < pts; ++i)
        {

            int n = pow(2,i+1);

            err = abs(exact-integrate_exact(n,order));
            c_err_array.push_back(err);
            
        }

        err_array.push_back(c_err_array);
        c_err_array.clear();

    }

    for (int i = 0; i < pts; ++i)
    {
        std::cout << pow(2,i+1) << " ";
        for (int j = 0; j < nods; ++j)
        {
            std::cout << err_array[j][i] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
//NA 2017/2018: Zadaća 5, Zadatak 1
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <limits>
#include <algorithm>
#include <complex>
#include <stdexcept>
#include <utility>

enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};
const double inf = std::numeric_limits<double>::infinity();
const double epsilon = std::numeric_limits<double>::epsilon();

template<typename FunTip>
bool BracketRoot(FunTip f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4)
{
    if(hinit <= 0 || hmax <= 0 || lambda <= 0) throw std::domain_error("Invalid parameters");
    double temp_A = x0, f1 = f(temp_A), h= hinit;
    while(std::abs(h) < hmax) {
        double temp_b = temp_A + h;
        double f2 = f(temp_b);
        if(f1 * f2 <= 0) {
            if(temp_A > temp_b) std::swap(temp_A,temp_b);
            a = temp_A;
            b = temp_b;
            return true;
        }
        h *= lambda;
        temp_A = temp_b;
        f1 = f2;
    }
    h = hinit;
    temp_A = x0;
    f1 = f(temp_A);
    while(std::abs(h) < hmax) {
        double temp_b = temp_A - h;
        double f2 = f(temp_b);
        if(f1 * f2 <= 0) {
            if(temp_A > temp_b) std::swap(temp_A,temp_b);
            a = temp_A;
            b = temp_b;
            return true;
        }
        h *= lambda;
        temp_b = temp_A;
        f2 = f1;
    }
    return false;
}

template<typename FunTip>
double RegulaFalsiSolve(FunTip f, double a, double b, RegulaFalsiMode mode= Slavic,double eps = 1e-10,int maxiter = 100)
{
    if(f(a) * f(b) > 0) throw std::range_error("Root must be bracketed");
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    double f1, f2, c, c_old;
    int br;

    if(mode == Unmodified) {
        f1 = f(a);
        f2 = f(b);
        c = a;
        c_old = b;
        br = 0;
        while(std::abs(c - c_old) > eps) {
            if(br == maxiter) throw std::logic_error("Given accuracy has not achieved");
            c_old = c;
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3 = f(c);
            if(std::abs(f3) < eps) return c;
            if(f1 * f3 < 0) {
                b = a;
                f2 = f1;
            }
            a = c;
            f1 = f3;
            br++;
        }
        return c;
    }
    else if(mode == Illinois) {
        f1 = f(a);
        f2 = f(b);
        c = a;
        c_old = b;
        br = 0;
        while(std::abs(c - c_old)>eps) {
            if(br == maxiter) throw std::logic_error("Given accuracy has not achieved");
            c_old = c;
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3 = f(c);
            if(std::abs(f3) < eps) return c;
            if(f1 * f3 < 0) {
                b = a;
                f2 = f1;
            }
            else f2 /= 2;
            a = c;
            f1 = f3;
            br++;
        }
        return c;
    }
    else if(mode == Slavic) {
        std::function<double(double)> help;
        help = [](double x) { return x/(1 + std::abs(x)); };
        f1 = help(f(a));
        f2 = help(f(b));
        c = a;
        c_old = b;
        br = 0;
        while(std::abs(c - c_old) > eps) {
            if(br == maxiter)
                throw std::logic_error("Given accuracy has not achieved");
            c_old = c;
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3 = help(f(c));
            if(std::abs(f3)<eps) return c;
            if(f1 * f3 < 0) {
                b = c;
                f2 = f3;
            }
            else {
                a = c;
                f1 = f3;
            }
            br++;
        }
        return c;
    }
    else if(mode == IllinoisSlavic) {
        std::function<double(double)> help;
        help = [](double x) { return x/(1 + std::abs(x)); };
         f1 = help(f(a));
        f2 = help(f(b));
        c = a;
        c_old = b;
        br = 0;
        while(std::abs(c - c_old) > eps) {
            if(br == maxiter)
                throw std::logic_error("Given accuracy has not achieved");
            c_old = c;
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3 = help(f(c));
            if(std::abs(f3)<eps) return c;
            if(f1 * f3 < 0) {
                b = a;
                f2 = f1;
            }
            else f2/=2;

            a = c;
            f1 = f3;
            br++;
        }
        return c;
    }
}

   template<typename FunTip>
   double RiddersSolve(FunTip f, double a, double b, double eps = 1e-10,int maxiter = 100)
   {
       if(f(a) * f(b) > 0) throw std::range_error("Root must be bracketed");
        if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");

            double f1=f(a), f2 = f(b);
            int br(0);
            while(std::abs(b-a) >= eps) {
                if(br == maxiter) throw std::logic_error("Given accuracy has not achieved");

                    double c = (a+b) / 2, f3 = f(c);
                    if(std::abs(f3) < eps) return c;
                    int sign;
                    if((f1 - f2) > 0) sign = 1;
                    else if((f1 - f2) < 0) sign  = -1;
                    else sign = 0;
                    double d = c + (f3 * (c - a) * sign) / std::sqrt(f3 * f3 - f1 * f2);
                    double f4 = f(d);
                    if(std::abs(f4) < eps) return d;
                    if(f3 * f4 <= 0) {
                        a = c;
                        b = d;
                        f1 = f3;
                        f2 = f4;
                    }
                    else if(f1 * f4 <= 0) {
                        b = d;
                        f2 = f4;
                    }
                    else {
                        a = d;
                        f1 = f4;
                    }
                    br++;
            }
       return (a + b) / 2;
   }

template <typename FunTip1, typename FunTip2>
double NewtonRaphsonSolve(FunTip1 f, FunTip2 fprim, double x0, double eps = 1e-10,double damping = 0, int maxiter = 100)
{
    // strana 311.
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    else if( damping < 0 || damping >= 1) throw std::domain_error("Invalid parameters");
    double dX(inf), v(f(x0)), d(fprim(x0));
    int br(0);

    while(std::fabs(dX) > eps) {

        if(std::abs(fprim(x0)) < eps || br == maxiter || !std::isfinite(x0)) throw std::logic_error("Convergence has not achieved");
        if(std::fabs(v) <= eps) return x0;
        dX = v / d;
        double w= v;
        v = f(x0 - dX);
        d = fprim(x0 - dX);
        while(std::abs(v) > std::abs(w) || !std::isfinite(v) || d == 0)
        {
            dX *= damping;
            v = f(x0 - dX);
            d = fprim(x0 - dX);
        }
        x0 -= dX;
        br++;
    }
    return x0;
}

bool operator==(std::complex<double>c1, std::complex<double> c2) {
    return(std::abs(c1.real() - c2.real()) < epsilon && std::abs(c1.imag() - c2.imag()) < epsilon);
}

std::complex<double>operator*(std::complex<double>c1, std::complex<double>c2) {
    return {c1.real()*c2.real()-c1.imag()*c2.imag(), c1.real()*c2.imag() + c1.imag()*c2.real()};
}

std::complex<double>operator *(double x, std::complex<double>c) {
    std::complex<double> temp(x,0);
    return temp * c;
}

std::complex<double>operator *(std::complex<double>c,double x) {
    return x * c;
}

std::pair<std::complex<double>,bool> Laguerre(std::vector<double>p, int n, std::complex<double>x, double eps, int maxiter)
{
    std::complex<double> dX(inf);
    int k = 1;
    while(std::abs(dX)> eps && k < maxiter) {
        std::complex<double> f = p[n];
        std::complex<double>d = 0;
        std::complex<double> s = 0;
        for(int i(n-1); i >= 0; i--) {
            s = s * x + 2 * d;
            d = d * x + f;
            f = f * x + p[i];
        }
        if(f == 0) return {x,true};
        std::complex<double> r = std::sqrt((n - 1) * ((n - 1) * d * d - n * f * s));
        if(std::abs(d + r) > std::abs(d - r)) dX = n * f / (d + r);
        else dX =  n * f / (d - r);
        x -= dX;
        k++;
    }
    if(std::abs(dX) <= eps) return {x,true};
    return {x,false};
}

std::pair<std::complex<double>,bool> Laguerre(std::vector<std::complex<double>>p, int n, std::complex<double>x, double eps, int maxiter)
{
    std::complex<double> dX(inf);
    int k = 1;
    while(std::abs(dX)> eps && k < maxiter) {
        std::complex<double> f = p[n];
        std::complex<double>d = 0;
        std::complex<double> s = 0;
        for(int i(n-1); i >= 0; i--) {
            s = s * x + 2 * d;
            d = d * x + f;
            f = f * x + p[i];
        }
        if(f == 0) return {x,true};
        std::complex<double> r = std::sqrt((n - 1) * ((n - 1) * d * d - n * f * s));
        if(std::abs(d + r) > std::abs(d - r)) dX = n * f / (d + r);
        else dX =  n * f / (d - r);
        x -= dX;
        k++;
    }
    if(std::abs(dX) <= eps) return {x,true};
    return {x,false};
}

std::vector<std::complex<double>> PolyRoots(std::vector<std::complex<double>> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10)
{
    if(eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    int n = coefficients.size() - 1, i(n), it(0);
    std::vector<std::complex<double>> a;
    while(i >= 1) {
        if(it == maxiters) throw std::logic_error("Convergence has not achieved");
        int t = 1;
        bool var(false);
        std::complex<double> x;
        while(!var && (t < maxtrials)) {
            x = {0.0};
            std::pair<std::complex<double>,bool> pair = Laguerre(coefficients,i,x,eps,maxiters);
            var = pair.second;
            x = pair.first;
            t++;
        }
        if(!var) throw std::logic_error("Convergence has not achieved");
        if(std::abs(x.imag()) <= eps) x = x.real();
        a.push_back(x);

        std::complex<double> v = coefficients[i];

        for(int j = i - 1; j >= 0; j--) {
            std::complex<double> w = coefficients[j];
            coefficients[j] = v;
            v = w + x * v;
        }
        i--;
        it++;
    }

    std::sort(a.begin(), a.end(), [] (std::complex<double> x, std::complex<double> y) {
        if (x.real() < y.real()) return true;
        else if (x.real() > y.real()) return false;
        return x.imag() < y.imag();
    });
    return a;
}

std::vector<std::complex<double>> PolyRoots(std::vector<double> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10)
{
    if(eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    int n = coefficients.size() - 1, i(n), it(0);
    std::vector<std::complex<double>> a(n + 1);
    while(i >= 1) {
        if(it == maxiters) throw std::logic_error("Convergence has not achieved");
        int t = 1;
        bool var(false);
        std::complex<double> x;
        while(!var && (t < maxtrials)) {
            x = {1,1};
            std::pair<std::complex<double>,bool> pair = Laguerre(coefficients,i,x,eps,maxiters);
            var = pair.second;
            x = pair.first;
            t++;
        }
        if(!var) throw std::logic_error("Convergence has not achieved");
        if(std::abs(x.imag()) <= eps) {
            x = x.real();
            a[i] = x;
            double v = coefficients[i];
            for(int j = i - 1; j >= 0; j--) {
                double w = coefficients[j];
                coefficients[j] = v;
                v = w + x.real() * v;
            }
            i--;
        }
        else {
            a[i] = x;
            a[i - 1] = std::conj(x);
            double alfa = 2 * x.real(), beta = std::abs(x);
            beta *= beta;
            double u = coefficients[i];
            double v = coefficients[i - 1] + alfa * u;
            for(int j = i - 2; j >= 0; j--) {
                double w = coefficients[j];
                coefficients[j] = u;
                u = v;
                v = w + alfa * v - beta * coefficients[j];
            }
            i -= 2;
        }
        it++;
    }
    a.erase(a.begin());
    std::sort(a.begin(),a.end(),[](std::complex<double>c1, std::complex<double>c2) {
        if(std::abs(c1.real() - c2.real()) > epsilon) return c1.real() < c2.real();
        return c1.imag() < c2.imag();
    });
    return a;
}


void testBracketRoot() {
    double a(0),b(0);
    try {
        if(BracketRoot([](double x) { return (x-3)/3; },0.1,a,b)) std::cout << "[" << a << "," << b << "]" << std::endl;
        else std::cout << "Nema nule" << std::endl;
        if(BracketRoot([](double x) { return x*x+1; },0.1,a,b)) std::cout << "[" << a << "," << b << "]" << std::endl;
        else std::cout << "Nema nule" << std::endl;
    }
    catch(std::domain_error e) {
        std::cout << e.what() << std::endl;
    }
    catch(...) {
        std::cout << "Unknown error" << std::endl;
    }
}

void testRegulaFalsi() {
    auto f1=[](double x) { return 0.05*(std::exp(10*(x-3))-1); };
    std::cout << "Slavic: " <<RegulaFalsiSolve(f1,1,4) << "\nIllinois: " << RegulaFalsiSolve(f1,0,4,Illinois) << "\nIllinoisSlavic: " <<RegulaFalsiSolve(f1,0,4,IllinoisSlavic) << std::endl;
}

void testRidders() {
    auto f1=[](double x) { return 0.05*(std::exp(10*(x-3))-1); };
    std::cout << "Ridders: " <<RiddersSolve(f1,1,4) << std::endl;
}

void testNewtonRaphson() {
    auto f=[](double x) { return x*x+5*x+6; };
    auto fprim=[](double x) { return 2*x+5; };
    std::cout << "Newton-Raphson: " << NewtonRaphsonSolve(f,fprim,1) << std::endl;
}

void testPolyRoots() {
    std::cout << "PolyRoots (kompleksni koeficijenti): ";
    std::vector<std::complex<double>> comp{{-42, -9}, {29, -6}, {-8, 3}, 1};
    std::vector<std::complex<double>> x{{1, 2}, 3, {4, -5}};
    std::vector<std::complex<double>> z1 = PolyRoots(comp);
    for(int i = 0; i < z1.size(); i++)
        std::cout << (std::abs(z1[i]-x[i]) < 1e-8) << " ";

    std::cout << std::endl;

    std::cout << "PolyRoots (realni koeficijenti): ";
    std::vector<double> real{6, 11, 6, 1};
    std::vector<std::complex<double>> y{-3, -2, -1};
    std::vector<std::complex<double>> z2 = PolyRoots(real);
    for(int i = 0; i < z2.size(); i++)
        std::cout << (std::abs(z2[i]-y[i]) < 1e-8) << " ";
    std::cout << std::endl;
}

int main ()
{
    testBracketRoot();
	std::cout << std::endl;
	testRegulaFalsi();
	std::cout << std::endl;
	testRidders();
	std::cout << std::endl;
	testNewtonRaphson();
	std::cout << std::endl;
	testPolyRoots();
	std::cout << std::endl;

	return 0;
}

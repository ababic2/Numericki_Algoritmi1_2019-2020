//NA 2017/2018: Zadaća 4, Zadatak 1
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <limits>
const double pi = 4 * atan(1);

class ChebyshevApproximation
{
    // std::function<double(double)> f;
    double xmin;
    double xmax;
    int n;
    int m;
    std::vector<double> c;
    ChebyshevApproximation(std::vector<double> c, double xmin, double xmax)
        : c(c), xmin(xmin), xmax(xmax), m(c.size() - 1) {}
public:

    template <typename FunTip>
    ChebyshevApproximation(FunTip f, double xmin, double xmax, int n) : xmin(xmin), xmax(xmax), n(n)
    {

        if(xmin >= xmax || n < 1) throw std::domain_error("Bad parameters");
        c.resize(n + 1);
        m = n;
        std::vector<double>w;
        for(int i = 0; i <= 4 * n + 3; i++) {
            double k = cos( pi * i / (2*n + 2));
            w.push_back(k);
        }

        std::vector<double>v;
        for (int i = 0; i <= n; i++) {
            double k = (xmin + xmax + (xmax - xmin)*w[2*i + 1]) / 2;
            v.push_back(f(k));
        }

        for(int k = 0; k <= n; k++) {
            double s = 0;
            for(int i  = 0; i <= n; i++) {
                int j = (k * (2 * i + 1)) % (4 * n + 4);
                s = s + v[i] * w[j];
            }
            c[k] = 2 * s / (n + 1);
        }
    }
    void set_m(int m);
    void trunc(double eps);
    double operator()(double x) const;
    double derivative(double x) const;
    ChebyshevApproximation derivative() const;
    ChebyshevApproximation antiderivative() const;
    double integrate(double a, double b) const;
    double integrate() const;
};

double ChebyshevApproximation::integrate(double a, double b) const
{
    if( a < xmin || a > xmax || b < xmin || b > xmax) throw std::domain_error("Bad interval");
    ChebyshevApproximation F (antiderivative());
    return F(b) - F(a);
}

double ChebyshevApproximation::integrate() const
{
    double s = 0;
    for(int k = 1; k <= (m + 1) / 2 ; k++)
        s += 2 * c[2 * k] / ( 1 - 4 * k * k);

    s *= (xmax - xmin) / 2;
    s += c[0] * (xmax - xmin) / 2;
    return s;
}

ChebyshevApproximation ChebyshevApproximation::antiderivative() const
{

    std::vector<double> koef(m + 2);
    koef[0] = 0;
    for(int k = 1; k <= m + 1; k++)
        koef[k] = (xmax - xmin) / (4 *k) * (c[k - 1] - (k < m?c[k + 1]:0));

    return ChebyshevApproximation(koef, xmin, xmax);
}

ChebyshevApproximation ChebyshevApproximation:: derivative() const
{
    double mi = 4./(xmax - xmin);
    std::vector<double>koef(c.size());
    koef[m - 1] = mi * m * c[m];
    koef[m - 2] = mi * (m - 1) * c[m - 1];
    for(int k = m - 3; k >= 0; k--)
        koef[k] = koef[k + 2] + mi * (k + 1) * c[k + 1];
    return ChebyshevApproximation(koef, xmin, xmax);
}

void ChebyshevApproximation::trunc(double eps)
{
    if( eps < 0 ) throw std::domain_error("Bad tolerance");
    for(int i = m; i >= 0; i--) {
        if(i < 1 ) throw std::domain_error("Bad tolerance");
        else if(fabs(c[i]) > eps) {
            m = i;
            break;
        }
    }
}

double ChebyshevApproximation::operator()(double x) const
{
    if(x > xmax || x < xmin) throw std::domain_error("Bad argument");
    double t = (2 * x - xmin - xmax) / (xmax - xmin);
    double p(1), q(t), s(c[0]/2 + c[1] * t);

    for(int k = 2; k <= m; k++) {
        double r = 2 * t * q - p;
        s += c[k] * r;
        p = q;
        q = r;
    }
    return s;
}

double ChebyshevApproximation::derivative(double x) const
{
    if(x > xmax || x < xmin) throw std::domain_error("Bad argument");
    double t = (2 * x - xmin - xmax) / (xmax - xmin);
    double p(1), q(4 * t), s(c[1] + q * c[2]),r;

    for(int k = 3; k <= m; k++) {
        r = k * (2 * t * q / ( k - 1) - p/(k - 2));
        s += c[k] * r;
        p = q;
        q = r;
    }
    return 2 * s  /(xmax - xmin);
}

void ChebyshevApproximation::set_m(int k)
{
    if(!(k > 1 && k < n)) throw std::domain_error("Bad order");
    this->m = m;
}


int main ()
{
    //TESTIRANJE IZUZETAKA

    try {
        ChebyshevApproximation ch([](double x) {
            return x;
        }, 0, -1, 10);
    } catch(std::domain_error e) {
        std::cout << "OK ";
        std::cout << e.what() << std::endl;
    }
    try {
        ChebyshevApproximation ch([](double x) {
            return x;
        }, 0, 1, 0);
    } catch(std::domain_error e) {
        std::cout << "OK ";
        std::cout << e.what() << std::endl;
    }

    ChebyshevApproximation ch([](double x) {
        return 3 * x * x * x + 2 * x * x + 5;
    }, 0, 10, 10);
    try {
        ch.set_m(11);
    } catch(std::domain_error e) {
        std::cout << "OK ";
        std::cout << e.what() << std::endl;
    }
    ch.set_m(9);
    try {
        std::cout << ch(-1) << std::endl;
    } catch(std::domain_error e) {
        std::cout << "OK ";
        std::cout << e.what() << std::endl;
    }
    try {
        ch.trunc(-5);
    } catch(std::domain_error e) {
        std::cout << "OK ";
        std::cout << e.what() << std::endl;
    }
    try {
        ChebyshevApproximation ch3([](double x) {
            return sin(x);
        }, 0, pi, 16);
        std::cout << ch3.integrate(-pi / 2, pi) << std::endl;
    } catch(std::domain_error e) {
        std::cout << "OK ";
        std::cout << e.what() << std::endl;
    }
    {
        ChebyshevApproximation ch([](double x) {
            return 3 * x * x * x + 2 * x * x + 5;
        }, 0, 10, 10);
        ch.set_m(9);
        std::cout << ch(1) << std::endl;
        std::cout << ch.derivative(1) << std::endl;
        std::cout << ch.derivative()(1) << std::endl;
        ChebyshevApproximation ch2([](double x) {
            return 1. / x;
        }, 1, 2.71, 16);
        ch2.trunc(0.0001);

        std::cout << ch2(1) << std::endl;
        std::cout << ch2(1) << std::endl;
        std::cout << ch2.antiderivative()(2.71) - ch2.antiderivative()(1) << std::endl;
        std::cout << ch2.integrate(1, 2.71) << std::endl;
        std::cout << ch2.integrate() << std::endl;

        ChebyshevApproximation ch3([](double x) {
            return sin(x);
        }, 0, pi, 16);

        std::cout << ch3(pi / 2) << std::endl;
        std::cout << ch3.antiderivative()(pi / 2) - ch3.antiderivative()(0) << std::endl;
        std::cout << ch3.integrate(pi / 4, pi / 2) << std::endl;
        std::cout << ch3.integrate() << std::endl;
    }
    return 0;
}

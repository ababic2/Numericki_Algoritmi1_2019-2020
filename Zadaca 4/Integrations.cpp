//NA 2017/2018: Zadaća 4, Zadatak 2
#include <iostream>
#include <cmath>
#include <vector>

const double pi = 4 * atan(1);

template <typename FunType>
std::pair<double, bool>RombergIntegration(FunType f, double a, double b,double eps=1e-8, int nmax = 1000000, int nmin = 50) {
    if(eps < 0 || nmax < 0 || nmin < 0 || nmax < nmin)
        throw std::domain_error("Bad parameter");
    int N = 2;
    double h = (b - a) / N;
    double s = (f(a) + f(b)) / 2;
    double Iold = s;
    std::vector<double>I;
    for(int i = 1; N <= nmax; i++) {
        for(int j = 1; j <= N / 2; j++)
            s += f(a + (2 * j - 1) * h);
        I.push_back(h * s);
        double p = 4;
        for(int k = I.size() - 2; k >= 0 ; k--) {
            I[k] = (p * I[k + 1] - I[k]) / (p - 1);
            p *= 4;
        }
        if(N >= nmin && std::fabs(I[0] - Iold) <= eps) return {I[0],true};
        Iold = I[0];
        h /= 2;
        N *= 2;
    }
    return {Iold, false};
}

template <typename FunType>
std::pair<double, bool>TanhSinhIntegration(FunType f, double a,double b,double eps=1e-8,int nmax = 1000000, int nmin = 20, double range = 3.5) {
    if(eps < 0 || nmin < 0 || range < 0 || nmax < nmin)
        throw std::domain_error("Bad parameter");

    int predznak = 1;
    if(a > b){
        predznak *= -1;
        std::swap(a,b);
    }
    double N(2), h(2 * range / N), p((b + a) / 2), q((b - a) / 2), s(0);
    double Iold(s);
    while(N < nmax) {
        for(int i = 1; i <= N /2; i++) {
            double t = -range + (2 * i - 1) * h;
            double u = pi * std::sinh(t) / 2;
            double v = f(p + q * std::tanh(u));
            if(std::isfinite(v))
                s += q * pi * std::cosh(t) * v / (2*std::cosh(u)*std::cosh(u));
        }
        double I = h * s;
        if(N >= nmin && std::fabs(I - Iold) <= eps) return {predznak*I, true};
        Iold = I;
        N*=2;
        h/=2;
    }
    return {Iold*predznak,false};
}

template <typename FunType>
std::pair<double,bool> AdaptiveAux(FunType f, double a, double b, double eps, double f1, double f2, double f3, double R) {
    if(!std::isfinite(f1)) f1 = 0;
    if(!std::isfinite(f2)) f2 = 0;
    if(!std::isfinite(f3)) f3 = 0;
    double c((a + b) / 2);
    double I1 = (b - a) * (f1 + 4 * f3 + f2) / 6;
    double f4 = f((a + c) / 2);
    double f5 = f((c + b) / 2);
    if(!std::isfinite(f4)) f4 = 0;
    if(!std::isfinite(f5)) f5 = 0;
    double I2 = (b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12;
    if(std::fabs(I1 - I2) <= eps) return {I2, true};
    if(R <= 0) return {I2, false};
    return AdaptiveAux(f, a, c, eps, f1, f3, f4, R-1) + AdaptiveAux(f, c, b, eps, f3, f2, f5, R-1);
}

std::pair<double,bool>operator+ (std::pair<double,bool>a, std::pair<double,bool>b) {
    return {a.first + b.first, a.second && b.second};
}

template <typename FunType>
std::pair<double,bool>AdaptiveIntegration(FunType f, double a, double b, double eps=1e-10, int maxdepth = 30, int nmin = 1) {
    if(eps < 0 || maxdepth < 0 || nmin < 0) throw std::domain_error("Bad parameter");
    double h((b - a)/nmin);
    std::pair<double,bool> s{0, true};
    for(int i = 1; i <= nmin; i++) {
        s = s + AdaptiveAux(f, a, a + h, eps, f(a), f(a + h), f((a + h) / 2), maxdepth);
        a += h;
    }
    return s;
}

int main ()
{
    {   //test f(x) = sin(x)
        auto sinf = [](double x) { return std::sin(x); };
        auto aig = AdaptiveIntegration(sinf, 0, pi);
        std::cout << aig.first << " " << aig.second << std::endl;;

        auto aig2 = TanhSinhIntegration(sinf, 0, pi);
        std::cout << aig2.first << " " << aig2.second << std::endl;

        auto aig3 = RombergIntegration(sinf, 0, pi);
        std::cout << aig3.first << " " << aig3.second << std::endl;;

    }

    {
        auto p = RombergIntegration([](double x) { return sin(x); }, 0, pi / 2);
        std::cout << p.first << " " << p.second << std::endl;

        p = RombergIntegration([](double x) { return exp(x); }, 0, 1);
        std::cout << p.first << " " << p.second << std::endl;

        p = RombergIntegration([](double x) { return 3 * x * x; }, 1, 3);
        std::cout << p.first << " " << p.second << std::endl;

        p = AdaptiveIntegration([](double x) { return sin(x); }, 0, pi / 2);
        std::cout << p.first << " " << p.second << std::endl;

        p = AdaptiveIntegration([](double x) { return exp(x); }, 0, 1);
        std::cout << p.first << " " << p.second << std::endl;

        p = AdaptiveIntegration([](double x) { return 3 * x * x; }, 1, 3);
        std::cout << p.first << " " << p.second << std::endl;

        p = TanhSinhIntegration([](double x) { return sin(x); }, 0, pi / 2);
        std::cout << p.first << " " << p.second << std::endl;

        p = TanhSinhIntegration([](double x) { return exp(x); }, 0, 1);
        std::cout << p.first << " " << p.second << std::endl;

        p = TanhSinhIntegration([](double x) { return 3 * x * x; }, 1, 3);
        std::cout << p.first << " " << p.second << std::endl;
    }
    {
        auto f = [](double x) { return x == 0 ? 0 : 1 / std::sqrt(x); };
        std::pair<double,bool> r1 = RombergIntegration(f, 0, 1);
        std::cout << r1.first << " " << r1.second << std::endl;

        auto f2 = [](double x) { return x == 0 ? 0 : 1 / std::sqrt(x); };
        r1 = AdaptiveIntegration(f2, 0, 1);
        std::cout << r1.first << " " << r1.second << std::endl;

        auto f3 = [](double x) { return x == 0 ? 0 : 1 / std::sqrt(x); };
        r1 = TanhSinhIntegration(f3, 0, 1);
        std::cout << r1.first << " " << r1.second << std::endl;
    }

	return 0;
}

#include<iostream>
#include<cmath>
#include<vector>
#include<algorithm>
#include <functional>
#define EPSILON 0.0001

template<typename FunType>
std::pair<double,bool>Limit(FunType f, double x0, double h = 0, double eps = 1e-8, double nmax = 20) {
    if(eps <= 0|| !(nmax > 3 && nmax <= 30)) throw std::domain_error("Invalid parameters");
    if(h < 0.00001) h = 0.001 * std::max(1.,std::fabs(x0)); // h ne nula !
    double yOld = std::numeric_limits<int>::max();

    std::vector<double> y(nmax);

    for(int i = 0; i < nmax; i++) {
        y.at(i) = f(x0 + h);
        double p = 2;
        for(int k = i - 1; k >= 0; k--) {
            y.at(k) = (p * y.at(k + 1) - y.at(k)) / (p - 1);
            p *= 2;
        }
        if(std::fabs(y.at(0) - yOld) < eps) return std::make_pair(y.at(0),true);
        yOld = y.at(0);
        h /= 2.;
    }
    return std::make_pair(y.at(0),false);
}

int main() {

  // test 1   // neki klasicni primjeri
    auto limes = Limit([] (double x) {return (x*x-9)/(x-3); }, 3);
    if (fabs(limes.first - 6) < EPSILON && limes.second == 1) std::cout << "1. OK\n";
    else std::cout << "1. NOT OK\n";

    // test 2
    limes = Limit([] (double x) {return (x-4)/(sqrt(x)-2);}, 4);
    if (fabs(limes.first - 4) < EPSILON && limes.second == 1) std::cout << "2. OK\n";
    else std::cout << "2. NOT OK\n";

    // test 3 // nepreciznost
    limes = Limit([] (double x) {return x/(x*x);}, 0.000001);
    if (limes.first < 1e+06 && limes.second == 0) std::cout << "3. OK\n";
    else std::cout << "3. NOT OK\n";

    // test 4   lijevi i desni limesi
    limes = Limit([] (double x) {return x/(x*x);}, -0.000001);
    if (limes.first > -1e+06 && limes.second == 0) std::cout << "4. OK\n";
    else std::cout << "4. NOT OK\n";

    // test 5
    limes = Limit([] (double x) {return std::abs(x)/x;}, 0);
    if (abs(limes.first - 1.) < 1e-5 && limes.second == 1) std::cout << "5. OK\n";
    else std::cout << "5. NOT OK\n";


    // test 6
    limes = Limit([] (double x) {return std::abs(sin(x))/sin(x);}, 0);
    if (abs(limes.first - 1) < EPSILON && limes.second == 1) std::cout << "6. OK\n";
    else std::cout << "7. NOT OK\n";

    // test 7 jos nepreciznosti
    limes = Limit([] (double x) {return (x*x - 2*x -3)/(x-1);}, 1);
    if (abs(limes.first + -4.1943+0.9) < EPSILON && limes.second == 0) std::cout << "7. OK\n";
    else std::cout << "7. NOT OK\n";


    // // testovi sa izuzecima
    try {
        limes = Limit([] (double x) {return x*x + x;}, 2, 0.01, 0);
        std::cout << "Izutetak 1. NOT OK\n";
    } catch(std::domain_error e) { std::cout << "Izutetak 1. OK\n"; }

    try {
        limes = Limit([] (double x) {return x*x + x;}, 2, 0.01, 1e-8, 2);
        std::cout << "Izutetak 2. NOT OK\n";
    } catch(std::domain_error e) { std::cout << "Izutetak 2. OK\n"; }

    try {
        limes = Limit([] (double x) {return x*x + x;}, 2, 0.01, 0, 31);
        std::cout << "Izutetak 3. NOT OK\n";
    } catch(std::domain_error e) { std::cout << "Izutetak 3. OK\n"; }


    try {
        limes = Limit([] (double x) {return x*x + x;}, 2, 0.01, -1);
        std::cout << "Izutetak 4. NOT OK\n";
    } catch(std::domain_error e) { std::cout << "Izutetak 4. OK\n"; }

    return 0;
}

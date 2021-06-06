//NA 2017/2018: Zadaća 5, Zadatak 2
#include <iostream>
#include <cmath>

template<typename FunTip>
double FindMinimum(FunTip f, double x0, double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
    if(eps<=0 || hinit<=0 || hmax<=0 || lambda<=0) throw std::domain_error("Invalid parameters");
    double a = x0 - hinit;
    double b = x0 - hinit;
    double c = x0;
    bool found(false);
    while(std::abs(hinit) < hmax) {
        if(f(c + hinit) < f(c)) {
            b = c + hinit;
            a = c - hinit;
        }
        else if(f(c - hinit) < f(c)) {
            b = c - hinit;
            a = b - hinit;
        }
        else {
            a = c - hinit;
            b = c + hinit;
            found = true;
            break;
        }
        c = b;
        hinit *= lambda;
    }
    // Zlantni presjek
    if(!found) throw std::logic_error("Minimum has not found");
    double golden = (1 + std::sqrt(5)) / 2;
    double d;
    if(std::abs(c - a) < std::abs(b - c)) d = b- (b - c)/golden;
    else {
        d = c;
        c = a + (c - a) / golden;
    }

    double u = f(c), v = f(d);
    while(std::abs(b - a) > eps) {
        if(u < v) {
            b = d;
            d = c;
            c = a + (c - a) / golden;
            v = u;
            u = f(c);
        }
        else {
            a = c;
            c = d;
            d = b - (b -d) / golden;
            u = v;
            v = f(d);
        }
    }
    return ( a + b) / 2;
}

void testFindMinimum() {
    std::cout << "FindMinimum: ";
    std::cout << FindMinimum([](double x) { return 1 + (x - 5) * (x - 5); }, 20);
    std::cout << std::endl;
}

int main ()
{
    // TESTOVI
    testFindMinimum();
	std::cout << std::endl;
    //y=x^2-6x-15 minimum=3
    std::cout << FindMinimum( [] (double x) {
        return x*x - 6*x - 15;
    }, 0);  //pocetna tacka lijevo od minimuma
    std::cout<<"  ";
    std::cout << FindMinimum( [] (double x) {
        return x*x - 6*x - 15;
    }, 6);  //pocetna tacka desno od minimuma
    std::cout << std::endl;

    //funkcija nema minimum
    try {
        std::cout << FindMinimum([](double x) {
            return 3*x;
        }, 1);
    } catch (std::logic_error e) {
        std::cout<<e.what()<<std::endl;
    }
    try {
        std::cout << FindMinimum([](double x) {
            return 1/x;
        }, 1);
    } catch (std::logic_error e) {
        std::cout<<e.what()<<std::endl;
    }

    //izuzeci
    try {
        FindMinimum([](double x) {
            return x;
        }, 1, -1, 1, 1, 1);
    } catch(std::domain_error e) {
        std::cout << "'" << e.what() << "'" << std::endl;
    }
    try {
        FindMinimum([](double x) {
            return x;
        }, 1, 1, -1, 1, 1);
    } catch (std::domain_error e) {
         std::cout << "'" << e.what() << "'" << std::endl;
    }
    try {
        FindMinimum([](double x) {
            return x;
        }, 1, 1, 1, -1, 1);
    } catch (std::domain_error e) {
         std::cout << "'" << e.what() << "'" << std::endl;
    }
    try {
        FindMinimum([](double x) {
            return x;
        }, 1, 1, 1, 1, -1);
    } catch (std::domain_error e) {
         std::cout << "'" << e.what() << "'" << std::endl;
    }
	return 0;
}

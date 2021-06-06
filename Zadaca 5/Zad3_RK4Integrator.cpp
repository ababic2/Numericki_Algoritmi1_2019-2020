//NA 2017/2018: Zadaća 5, Zadatak 3
#include <iostream>
#include <cmath>
#include <vector>

template<typename FunTip>
double RK4Step(FunTip f, double x, double y, double h) {
    double k1 = f(x,y);
    double k2 = f(x + h/2, y + h * k1 / 2);
    double k3 = f(x + h/2, y + h * k2 / 2);
    double k4 = f(x + h, y + h * k3);
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

template<typename FunTip>
std::vector<std::pair<double,double>> RK4Integrator(FunTip f, double x0, double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false)
{
    if((h > 0 && xmax < x0) || (h <= 0 && xmax > x0)) return {{x0,y0}};
    std::vector<std::pair<double,double>>koef;
    if(!adaptive) {
        double x = x0, y = y0;
        if(h > 0) {
            while(x <= xmax + eps) {
                koef.push_back({x,y});
                double k1 = f(x,y);
                double k2 = f(x + h/2, y + h * k1 / 2);
                double k3 = f(x + h / 2, y + h * k2 / 2);
                double k4 = f(x + h, y + h * k3);
                y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                x += h;
            }
        }
        else {
            while(x >= xmax - eps) {
                koef.push_back({x,y});
                double k1 = f(x,y);
                double k2 = f(x + h/2, y + h * k1 / 2);
                double k3 = f(x + h / 2, y + h * k2 / 2);
                double k4 = f(x + h, y + h * k3);
                y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                x += h;
            }
        }
    }
    else {
        // alforitam strana 411.
        if(h > 0) {
            double x = x0, y = y0;
            koef.push_back({x,y});
            while(x <= xmax + eps) {
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f, x + h/2,u, h/2);
                double w = RK4Step(f,x,y,h);
                double delta = std::abs(w - v) / h;
                if(delta <= eps) {
                    x += h;
                    y = v;
                    koef.push_back({x,y});
                }
                h = h * std::min(5.0,0.9 * std::pow(eps/delta,1/4.));
            }
            if(std::abs(koef[koef.size() - 1].first - xmax) > eps) {
                koef.erase(koef.begin());
                h = xmax - koef[koef.size() - 1].first;
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f, x + h / 2, u , h/2);
                double w = RK4Step(f,x,y,h);
                koef[koef.size() - 1] = {xmax,v};
            }
        }
        else {
            double x = x0, y = y0;
            koef.push_back({x,y});
            while(x >= xmax - eps) {
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f, x + h / 2, u , h/2);
                double w = RK4Step(f,x,y,h);
                double delta = std::abs(w-v) / ((-1) * h);
                if(delta <= eps) {
                    x += h;
                    y = v;
                    koef.push_back({x,y});
                }
                h *= std::min(5.0,0.9 * std::pow(eps/delta,1/4.));
            }
            if(std::abs(koef[koef.size()-1].first - xmax) > eps) {
                koef.erase(koef.begin());
                h = xmax - koef[koef.size() - 1].first;
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f,x + h/2, u, h/2);
                double w = RK4Step(f,x,y,h);
                koef[koef.size() - 1] = {xmax,v};
            }
        }
    }
    return koef;
}


void testRK4() {
    //Diferencijalna jednačina y'=x-y+3

    //bez adaptacije, pozitivni korak
    try {
        std::vector<std::pair<double,double>> v = RK4Integrator([](double x, double y) { return x-y+3; }, 0, 1, 1.5, 0.1);
        auto f = [](double x) { return -1*std::exp((-1)*x)+x+2; };
        for(int i = 0; i < v.size(); i++) {
            double x = v[i].first;
            std::cout << x << " " << v[i].second << " " << f(x) << std::endl;
        }
    }
    catch(...) {
        std::cout << "Error";
    }
    std::cout << std::endl;
    //bez adaptacije, negativni korak
    try {
        std::vector<std::pair<double,double>> v = RK4Integrator([](double x, double y) { return x-y+3; }, 0, 1, -1.5, -0.1);
        auto f = [](double x) { return -1*std::exp((-1)*x)+x+2; };
        for(int i = 0; i < v.size(); i++) {
            double x = v[i].first;
            std::cout << x << " " << v[i].second << " " << f(x) << std::endl;
        }
    }
    catch(...) {
        std::cout << "Error";
    }
    std::cout << std::endl;

    //sa adaptacijom, pozitivan korak
    try {
        std::vector<std::pair<double,double>> v = RK4Integrator([](double x, double y) { return x-y+3; }, 0, 1, 1.5, 0.1, 1e-8, true);
        auto f = [](double x) { return -1*std::exp((-1)*x)+x+2; };
        for(int i = 0; i < v.size(); i++) {
            double x = v[i].first;
            std::cout << x << " " << v[i].second << " " << f(x) << std::endl;
        }
    }
    catch(...) {
        std::cout << "Error";
    }
    std::cout << std::endl;

    //sa adaptacijom, negativan korak
    try {
        std::vector<std::pair<double,double>> v = RK4Integrator([](double x, double y) { return x-y+3; }, 0, 1, -1.5, -0.1, 1e-8, true);
        auto f = [](double x) { return -1*std::exp((-1)*x)+x+2; };
        for(int i = 0; i < v.size(); i++) {
            double x = v[i].first;
            std::cout << x << " " << v[i].second << " " << f(x) << std::endl;
        }
    }
    catch(...) {
        std::cout << "Error";
    }
    std::cout << std::endl;
}

int main ()
{
    //TESTOVI

    testRK4();
	std::cout << std::endl;

    //y'=3x+2y
    std::cout<<"Jednacina y'=3x+2y, bez adaptacije: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return 3 * x + 2*y;
        },
        0, 1, 1.2, 0.01);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 1.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;
    std::cout<<"Jednacina y'=3x+2y, sa adaptacijom: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return 3 * x + 2*y;
        },
        0, 1, 1.2, 0.01, true);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 1.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;

    std::cout<<"Jednacina y'=3x+2y, bez adaptacije, unazad: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return 3 * x + 2*y;
        },
        0, 1, -1.2, -0.01);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 1.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;


    std::cout<<"Jednacina y'=3x+2y, sa adaptacijom, unazad: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return 3 * x + 2*y;
        },
        0, 1, -1.2, -0.01,true);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 1.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;






//y'=x+2y+1
    std::cout<<"Jednacina y'=x+2y+1, bez adaptacije: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return  x + 2*y + 1;
        },
        0, 1, 1.1, 0.01);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 0.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;

    std::cout<<"Jednacina y'=x+2y+1, sa adaptacijom: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return  x + 2*y + 1;
        },
        0, 1, 1.1, 0.01, true);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 0.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;

    std::cout<<"Jednacina y'=x+2y+1, bez adaptacije, unazad: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return  x + 2*y + 1;
        },
        0, 1, -1.1, -0.01);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 0.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;

    std::cout<<"Jednacina y'=x+2y+1, sa adaptacijom, unazad: "<<std::endl;
    {
        auto rk4 = RK4Integrator([](double x, double y) {
            return  x + 2*y + 1;
        },
        0, 1, -1.1, -0.01, true);
        auto rj = [](double x) {
            return 1.75* std::exp(2*x) - 0.5*x - 0.75;
        };
        for(int i = 0; i < rk4.size(); i++) {
            double x = rk4[i].first;
            std::cout <<"x: "<< x << "  y: " << rk4[i].second << "  ispravno y: " << rj(x) << std::endl;
        }
    }
    std::cout<<std::endl;

	return 0;
}

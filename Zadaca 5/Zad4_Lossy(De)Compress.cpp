#include<iostream>
#include <complex>
#include <vector>
#include <cmath>

// u periodu kada je pao c9, napisano

void FFT(std::vector<double> &x, std::vector<std::complex<double>> &y, int vel, int s = 0, int d = 0, int t = 1) {
    const double pi = 4 * atan(1);
    if(vel == 1) y[d] = x[s];
    else {
        FFT(x,y,vel/2,s,d,t*2);
        FFT(x,y,vel/2,d+vel/2,t*2);
        std::complex<double> mi = 1;
        std::complex<double> w = std::complex<double>(std::cos(2*pi/vel),std::sin(2*pi/vel));
        std::complex<double> temp{1,0};
        std::complex<double> v = temp / w;
        for(int i = d; i< (d + vel/2 - 1); i++) {
            std::complex<double> c = y[i];
            std::complex<double>c2 = mi * y[i + vel/2];
            y[i] = c + c2;
            y[i + vel/2] = c - c2;
            mi *= v;
        }
    }
}


std::vector<double> LossyCompress(std::vector<double> data, int new_size) {

    if (new_size <= 1 || new_size > data.size()) throw std::range_error("Bad new size");
    else if (!( (data.size() & (data.size() - 1)) == 0)) throw std::range_error("Data size must be a power of two");

    const double pi = 4*atan(1);
    std::vector<double> v1(data.size());
    std::vector<std::complex<double>> v2(data.size());

    for (int i = 0; i < data.size()/2; i++) v1[i] = data[i*2];
    for (int i = data.size()/2; i < data.size(); i++)   v1[i] = data[2*(data.size()-i) - 1];

    FFT(v1, v2, v1.size());
    v2[new_size - 1] = data.size();

    std::vector<double> sekvenca(v2.size());
    for (int i = 0; i < v2.size()-1; i++) {
        std::complex<double> w = std::complex<double>(std::cos(pi/v2.size()), std::sin(pi/v2.size()));
        w = std::pow(w, (-1*i)/2.);
        sekvenca[i] = (w*v2[i]).real();
    }
    return sekvenca;
}

std::vector<double> LossyDecompress(std::vector<double> compressed) {

    int vel = compressed[compressed.size() - 1];
    if (vel < 1 || vel < compressed.size()) throw std::logic_error("Bad compressed sequence");
    else if (!( (compressed.size() & (compressed.size() - 1)) == 0)) throw std::range_error("Data size must be a power of two");

    for (int i = 0; i < vel; i++) compressed[i] = 0;    // popunimo je nulama do n-1
    return compressed;
}
int main() {
    try {
    LossyCompress({1, 6, 3, 5}, -5);
}
catch(std::range_error e) {
    std::cout << "'" << e.what() << "'" << std::endl;
}
try {
    LossyCompress({1, 2, 3, 4}, 6);
}
catch(std::range_error e) {
    std::cout << "'" << e.what() << "'" << std::endl;
}
try {
    LossyCompress({1, 2, 3, 4, 5}, 2);
}
catch(std::range_error e) {
    std::cout << "'" << e.what() << "'" << std::endl;
}
try {
    LossyDecompress({0, 0, 0, 0, 3});
}
catch(std::logic_error e) {
    std::cout << "'" << e.what() << "'" << std::endl;
}
try {
    LossyDecompress({0, 0, 0, 0, 0, 0, 4});
}
catch(std::logic_error e) {
    std::cout << "'" << e.what() << "'" << std::endl;
}
try {
    LossyDecompress({0, 0, 0, 0, 1.5});
}
catch(std::logic_error e) {
    std::cout << "'" << e.what() << "'";
}


    return 0;
}

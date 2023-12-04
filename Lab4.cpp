#include <iostream>
#include <iomanip>
#include <array>
#include <string>
#include <string.h>
#include <vector>
#include <cmath>
#include <math.h>
#include <functional>

using namespace std;

const long double a = 0.000000000000001;
const long double b = 1;
const long double eps = 1e-6;

long double f(long double x) {
    long double result = x / log(1 + x);
    return result;
}

double trapezoidal() {
    long double r = 1;
    int n = 1;
    long double I2n = 0;
    while (abs(r) > eps) {
        long double In = I2n;
        long double h = (b - a) / n;
        I2n = (f(b) - f(a)) / 2;
        for (int i = 0; i < n - 1; ++i) {
            I2n += f(a + (i + 1) * h);
        }
        I2n *= h;
        r = (I2n - In) / 3;
        n *= 2;
    }
    return I2n;
}

double parabolic() {
    long double r = 1;
    int n = 1;
    long double I2n = 0;
    while (abs(r) > eps) {
        long double In = I2n;
        long double h = (b - a) / n;
        I2n = f(b) - f(a);
        for (int i = 0; i < n / 2 - 1; ++i) {
            I2n += 2 * f(a + 2 * (i + 1) * h);
        }
        for (int i = 0; i < n / 2; ++i) {
            I2n += 4 * f(a + (2 * i + 1) * h);
        }
        I2n = I2n * h / 3;
        r = (I2n - In) / 15;
        n *= 2;
    }
    return I2n;
}


int main() {
    double result = 1.2292741343616127837489278679213389019648421242534218341935682723331333647814503772818287462097859244;
    double A = trapezoidal();
    double B = parabolic();

    std::cout << std::setprecision(17);

    std::cout << "Trapezoidal:" << '\n' << A << '\n' << "Difference:" << '\n' << A - result << '\n' << '\n';

    std::cout << "Parabolic:" << '\n' << B << '\n' << "Difference:" << '\n' << B - result << '\n' << '\n';

    return 0;
}

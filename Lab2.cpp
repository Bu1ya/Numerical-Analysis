#include <iostream>
#include <iomanip>
#include <array>
#include <string>
#include <string.h>
#include <vector>
#include <cmath>
#include <math.h>


using namespace std;


const long double e = 2.7182818284590452;

const long double A = 0; // [a;b]
const long double B = 1;

const long double Alpha = 0.2;
const int n = 20;

vector<long double> generateNodes(string type) {

    vector<long double> nodes;

    if (type == "aUsingDefaultNodes" || type == "aUsingDisplacedNodes") {
        for (int i = 0; i <= n; i++) {
            nodes.push_back(A + i * ((B - A) / n));
        }
    }
    else if (type == "bUsingDefaultNodes" || type == "bUsingDisplacedNodes") {
        for (int i = 0; i <= n; i++) {
            nodes.push_back((A + B - ((B - A) * cos(((2 * i + 1) * atan(1) * 4) / (2 * (n + 1))))) / 2);
        }
    }
    else {
        cout << "Invalid type";
    }
    return nodes;
}

long double function(long double x) {
    long double result = pow(e, -1 * x) * sin(x) + pow(x, 3);
    return result;
}

vector<long double> tridiagonalMatrixAlgorithm(vector<long double> a, vector<long double> b,
                                               vector<long double> c, vector<long double> d) {
    int count = b.size();
    a.insert(a.begin(), 0);
    c.push_back(0);
    vector<long double> alpha(count, 0), beta(count, 0), x(count, 0);

    //forward
    for (int i = 1; i < count; i++) {
        alpha[i] = c[i] / (b[i] - a[i] * alpha[i - 1]);
        beta[i] = (d[i] - a[i] * beta[i - 1]) / (b[i] - a[i] * alpha[i - 1]);
    }

    //backward
    x[count - 1] = beta[count - 1];
    for (int i = count - 2; i >= 0; i--) {
        x[i] = beta[i] - alpha[i] * x[i + 1];
    }

    return x;
}

long double splinePolynomial(long double m0, long double m1, long double x0, long double x1, long double x) {
    return (m0 * (pow((x1 - x), 3) / (6 * (x1 - x0))) +
        m1 * (pow((x - x0), 3) / (6 * (x1 - x0))) +
        (function(x0) - (m0 * (pow((x1 - x0), 2) / 6))) * ((x1 - x)/(x1 - x0)) +
        (function(x1) - (m1 * (pow((x1 - x0), 2) / 6))) * ((x - x0) / (x1 - x0))
        );
}

long double cubicSpline(vector<long double> nodes, string type) {
    
    vector<long double> a, b, c, d, m;

    vector<long double> displacedNodes = nodes;

    for (int i = 2; i < n; i++) {
        a.push_back((nodes[i] - nodes[i - 1]) / 6);
        c.push_back((nodes[i] - nodes[i - 1]) / 6);
    }
    for (int i = 2; i <= n; i++) {
        b.push_back(((nodes[i] - nodes[i - 1]) + (nodes[i - 1] - nodes[i - 2])) / 6);
    }
    for (int i = 2; i <= n; i++) {
        d.push_back(( - (function(nodes[i - 1]) - function(nodes[i - 2])) / (nodes[i - 1] - nodes[i - 2])) + 
                       ((function(nodes[i]) - function(nodes[i - 1])) / (nodes[i] - nodes[i - 1])));
    }

    m = tridiagonalMatrixAlgorithm(a, b, c, d);
    m.insert(m.begin(), 0);
    m.push_back(0);


    if (type == "aUsingDisplacedNodes" || type == "bUsingDisplacedNodes") {
        for (int i = 0; i < displacedNodes.size(); i++) {
            displacedNodes[i] += Alpha * ((B - A) / n);
        }
    }


    for (int i = 1; i < displacedNodes.size() - 1; i++) {
        long double polynomial = splinePolynomial(m[i - 1], m[i], nodes[i - 1], nodes[i], displacedNodes[i]);
        cout << fixed << " | " << setw(21) << displacedNodes[i] << " | "
            << setw(20) << function(displacedNodes[i]) << " | " 
            << setw(20) << polynomial << " | "
            << setw(21) << function(displacedNodes[i]) - polynomial << " | " << '\n';
    }


    return 0;
}

void printHeader() {
    cout << '\n';
    cout << " | " << "          x          " << " | "
        << "        f(x)        " << " | "
        << "        S(x)        " << " | "
        << "     |f(x)-S(x)|     " << " | " << '\n';
    cout << " |-----------------------|----------------------|----------------------|-----------------------|" << '\n';
}

int main()
{
    std::cout << std::setprecision(17);
    string aDef = "aUsingDefaultNodes",
        aDis = "aUsingDisplacedNodes",
        bDef = "bUsingDefaultNodes",
        bDis = "bUsingDisplacedNodes";
    cout << "\nResult Using Displaced Linear Nodes Generation:\n";
    printHeader();
    cubicSpline(generateNodes(aDis), aDis);
    cout << "\nResult Using Displaced Cosinus Nodes Generation:\n";
    printHeader();
    cubicSpline(generateNodes(bDis), bDis);
}
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

const long double a = -1;
const long double b = 2;
const long double alpha = 0.2;
const int n = 20;


vector<long double> generateNodes(string type) {
    
    vector<long double> nodes;
    
    if (type == "aUsingDefaultNodes") {
        for (int i = 0; i <= n; i++) {
            nodes.push_back(a + i * ((b - a) / n));
        }
    }
    else if (type == "bUsingDefaultNodes") {
        for (int i = 0; i <= n; i++) {
            nodes.push_back((a + b - ((b - a) * cos(((2 * i + 1) * atan(1)*4) / (2 * (n + 1)))))/2);
        }
    }
    else if (type == "aUsingDisplacedNodes") {
        for (int i = 0; i <= n; i++) {
            nodes.push_back(a + i * ((b - a) / n));
        }
    }
    else if (type == "bUsingDisplacedNodes") {
        for (int i = 0; i <= n; i++) {
            nodes.push_back((a + b - ((b - a) * cos(((2 * i + 1) * atan(1) * 4) / (2 * (n + 1))))) / 2);
        }
    }
    else {
        cout << "Invalid type";
    }
    return nodes;
}

long double function(long double x) {
    long double result = pow(e, -2 * x) * sin(x) + pow(x, 3);
    return result;
}

long double dividedDifference(int l, int r, vector<long double> nodes) {
    long double result = 0;

    int tempIndex = l;

    for (l; l <= r; l++) {
        long double numerator = function(nodes[l]), denominator = 1;
        for (int j = tempIndex; j <= r; j++) {
            if (j == l) {
                continue;
            }
            denominator *= (nodes[l] - nodes[j]);
        }
        result += numerator / denominator;
    }

    return result;
}

long double lagrangePolynomial(long double x, vector<long double> nodes, string mode) {
    long double result = 0;

    if (mode == "aUsingDisplacedNodes" || mode == "bUsingDisplacedNodes") {
        //for (int i = 0; i <= n; i++) {
        //    nodes[i] += alpha * ((b - a) / n);
        //}
        x += alpha * ((b - a) / n);
    }

    for (int i = 0; i <= n; i++) {
        long double tempResult = function(nodes[i]);
        for (int j = 0; j <= n; j++) {
            if (i == j) {
                continue;
            }
            tempResult *= ((x - nodes[j]) / (nodes[i] - nodes[j]));
        }
        result += tempResult;
    }

    return result;
}

long double newtonForwardPolynomial(long double x, vector<long double> nodes, string mode) {
    long double result = function(nodes[0]);

    vector<long double> nodesBefore = nodes;

    if (mode == "aUsingDisplacedNodes" || mode == "bUsingDisplacedNodes") {
        for (int i = 0; i <= n; i++) {
            nodes[i] += alpha * ((b - a) / n);
        }
        x += alpha * ((b - a) / n);
    }

    for (int i = 1; i <= n; i++) {
        long double tempResult = dividedDifference(0, i, nodesBefore);
        for (int j = 0; j <= i - 1; j++) {
            tempResult *= (x - nodes[j]);
        }
        result += tempResult;
    }

    return result;
}

long double newtonBackwardsPolynomial(long double x, vector<long double> nodes, string mode) {
    long double result = function(nodes[n]);

    vector<long double> nodesBefore = nodes;

    if (mode == "aUsingDisplacedNodes" || mode == "bUsingDisplacedNodes") {
        for (int i = 0; i <= n; i++) {
            nodes[i] += alpha * ((b - a) / n);
        }
        x += alpha * ((b - a) / n);
    }

    for (int i = 1; i <= n; i++) {
        long double tempResult = dividedDifference(n-i, n, nodesBefore);
        for (int j = 0; j <= i - 1; j++) {
            tempResult *= (x - nodes[n - j]);
        }
        result += tempResult;
    }

    return result;
}

void printHeader() {
    cout << " ***********************************************************************************************" << '\n';
    cout << " | " << "          x          " << " | "
        << "        f(x)        " << " | "
        << "        P(x)        " << " | "
        << "     |f(x)-P(x)|     " << " | " << '\n';
    cout << " ***********************************************************************************************" << '\n';
}

void printLagrangePolynomialResultUsingDefaultLinearNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Lagrange Polynomial Result Using Default Linear Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << lagrangePolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - lagrangePolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printLagrangePolynomialResultUsingDisplacedLinearNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Lagrange Polynomial Result Using Displaced Linear Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i] + alpha * ((b - a) / n)) << " | "
            << setw(20) << lagrangePolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i] + alpha * ((b - a) / n)) - lagrangePolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printLagrangePolynomialResultUsingDefaultCosinusNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Lagrange Polynomial Result Using Default Cosinus Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << lagrangePolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - lagrangePolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printLagrangePolynomialResultUsingDisplacedCosinusNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Lagrange Polynomial Result Using Displaced Cosinus Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i] + alpha * ((b - a) / n)) << " | "
            << setw(20) << lagrangePolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i] + alpha * ((b - a) / n)) - lagrangePolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonForwardPolynomialResultUsingDefaultLinearNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Forward Polynomial Result Using Default Linear Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonForwardPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonForwardPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonForwardPolynomialResultUsingDisplacedLinearNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Forward Polynomial Result Using Displaced Linear Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonForwardPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonForwardPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonForwardPolynomialResultUsingDefaultCosinusNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Forward Polynomial Result Using Default Cosinus Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= 20; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonForwardPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonForwardPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonForwardPolynomialResultUsingDisplacedCosinusNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Forward Polynomial Result Using Displaced Cosinus Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonForwardPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonForwardPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonBackwardsPolynomialResultUsingDefaultLinearNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Backwards Polynomial Result Using Default Linear Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonBackwardsPolynomialResultUsingDisplacedLinearNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Backwards Polynomial Result Using Displaced Linear Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonBackwardsPolynomialResultUsingDefaultCosinusNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Backwards Polynomial Result Using Default Cosinus Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

void printNewtonBackwardsPolynomialResultUsingDisplacedCosinusNodesGeneration(string mode) {
    vector<long double> nodes = generateNodes(mode);
    cout << " Newton Backwards Polynomial Result Using Displaced Cosinus Nodes Generation:" << '\n';
    printHeader();
    for (int i = 0; i <= n; i++) {
        cout << fixed << " | " << setw(21) << nodes[i] << " | "
            << setw(20) << function(nodes[i]) << " | "
            << setw(20) << newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | "
            << setw(21) << function(nodes[i]) - newtonBackwardsPolynomial(nodes[i], nodes, mode) << " | " << '\n';
    }
    cout << " ***********************************************************************************************" << '\n';
}

int main()
{
    std::cout << std::setprecision(17);
    printLagrangePolynomialResultUsingDefaultLinearNodesGeneration("aUsingDefaultNodes");
    printLagrangePolynomialResultUsingDisplacedLinearNodesGeneration("aUsingDisplacedNodes");
    printLagrangePolynomialResultUsingDefaultCosinusNodesGeneration("bUsingDefaultNodes");
    printLagrangePolynomialResultUsingDisplacedCosinusNodesGeneration("bUsingDisplacedNodes");

    printNewtonForwardPolynomialResultUsingDefaultLinearNodesGeneration("aUsingDefaultNodes");
    printNewtonForwardPolynomialResultUsingDisplacedLinearNodesGeneration("aUsingDisplacedNodes");
    printNewtonForwardPolynomialResultUsingDefaultCosinusNodesGeneration("bUsingDefaultNodes");
    printNewtonForwardPolynomialResultUsingDisplacedCosinusNodesGeneration("bUsingDisplacedNodes");

    printNewtonBackwardsPolynomialResultUsingDefaultLinearNodesGeneration("aUsingDefaultNodes");
    printNewtonBackwardsPolynomialResultUsingDisplacedLinearNodesGeneration("aUsingDisplacedNodes");
    printNewtonBackwardsPolynomialResultUsingDefaultCosinusNodesGeneration("bUsingDefaultNodes");
    printNewtonBackwardsPolynomialResultUsingDisplacedCosinusNodesGeneration("bUsingDisplacedNodes");
}
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;
const string FILENAME = "cubicSpline.txt";

double func(double x) {
    return log10(x);
}


struct SplineSegment {
    double a, b;
    double A, B;
    double x_start, x_end;
};

// s'(x)
double splineFirstDerivative(double x, const SplineSegment& segment) {
    if (x < segment.x_start || x > segment.x_end) {
        cerr << "Error: x is out of range for the spline segment." << endl;
        return 0.0;
    }

    double term1 = -3 * segment.a * pow(segment.x_end - x, 2);
    double term2 = 3 * segment.b * pow(x - segment.x_start, 2);
    double term3 = -segment.A;
    double term4 = +segment.B;
    //cout << "Derivative1 of x = " << x << " : " << term1 << " + " << term2 << " + " << term3 << " + " << term4 << endl;
    return term1 + term2 + term3 + term4;
}

// s''(x)
double splineSecondDerivative(double x, const SplineSegment& segment) {
    if (x < segment.x_start || x > segment.x_end) {
        cerr << "Error: x is out of range for the spline segment." << endl;
        return 0.0;
    }

    double term1 = 6 * segment.a * (segment.x_end - x);
    double term2 = 6 * segment.b * (x - segment.x_start);

    return term1 + term2;
}

void printMatrix(const vector<vector<double>>& matrix, const string& name) {
    cout << name << ":\n";
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << setw(10) << setprecision(4) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printVector(const vector<double>& vec, const string& name) {
    cout << name << ":\n";
    for (double val : vec) {
        cout << setw(10) << setprecision(4) << val << " ";
    }
    cout << endl << endl;
}

double calculateSplineValue(double x, const vector<double>& xs, const vector<double>& ms, const vector<double>& A, const vector<double>& B) {
    int n = xs.size();

    int i = 0;
    while (i < n - 1 && x > xs[i + 1]) {
        ++i;
    }

    if (i >= n - 1 || x < xs[0]) {
        cerr << "Error: x = " << x << " is outside the spline range [" << xs[0] << ", " << xs[n - 1] << "].\n";
        return 0.0;
    }

    double h = xs[i + 1] - xs[i];

    double term1 = (ms[i] / (6.0 * h)) * pow(xs[i + 1] - x, 3);
    double term2 = (ms[i + 1] / (6.0 * h)) * pow(x - xs[i], 3);
    double term3 = (A[i] / h) * (xs[i + 1] - x);
    double term4 = (B[i] / h) * (x - xs[i]);

    return term1 + term2 + term3 + term4;
}

int main() {
    double x_min, x_max;
    x_min = 1.; x_max = 100.;

    vector<double> x, f;
    for (double i = x_min; i <= x_max - 1; i += 7) {
        x.push_back(i);
        f.push_back(func(i));
    }
    x.push_back(x_max);
    f.push_back(func(x_max));
    int n = x.size();

    cout << setw(10) << "x" << setw(15) << "f(x)" << endl;
    for (const auto& value : x) {
            std::cout << setw(10)  << value <<  setw(10) << " " << func(value) ;
        std::cout << std::endl;
    }

    // h_i
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }
    printVector(h, "h (steps)");

    // C (n-2 x n-2)
    vector<vector<double>> C(n - 2, vector<double>(n - 2, 0.0));
    for (int i = 0; i < n - 2; ++i) {
        if (i > 0) C[i][i - 1] = h[i] / 6.0;
        C[i][i] = (h[i] + h[i + 1]) / 3.0;
        if (i < n - 3) C[i][i + 1] = h[i + 1] / 6.0;
    }
    printMatrix(C, "C (tridiagonal matrix)");

    // H (n-2 x n)
    vector<vector<double>> H(n - 2, vector<double>(n, 0.0));
    for (int i = 0; i < n - 2; ++i) {
        H[i][i] = 1.0 / h[i];
        H[i][i + 1] = -(1.0 / h[i] + 1.0 / h[i + 1]);
        H[i][i + 2] = 1.0 / h[i + 1];
    }
    printMatrix(H, "H matrix");

    //  b = H * f
    vector<double> b(n - 2, 0.0);
    for (int i = 0; i < n - 2; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i] += H[i][j] * f[j];
        }
    }
    printVector(b, "b");

    //  C * m = b
    vector<double> m(n, 0.0);
    vector<double> alpha(n - 2, 0.0), beta(n - 2, 0.0);

    //cout << "\nDirect pass:\n";
    for (int i = 1; i < n - 2; ++i) {
        double t = C[i][i - 1] / C[i - 1][i - 1];
        C[i][i] -= t * C[i - 1][i];
        b[i] -= t * b[i - 1];

        /*cout << "Step " << i << ":\n";
        printMatrix(C, "C (after modification)");
        printVector(b, "b (after modification)");*/
    }
    //cout << "\nReverse pass:\n";
    m[n - 2] = b[n - 3] / C[n - 3][n - 3];

    for (int i = n - 3; i >= 0; --i) {
        m[i + 1] = (b[i] - C[i][i + 1] * m[i + 2]) / C[i][i];
        //cout << "m[" << i + 1 << "] = " << m[i+1] << endl;
    }
    printVector(m, "m (second derivatives)");

    vector<double> A(n - 1), B(n - 1);
    for (int i = 1; i < n; ++i) {
        A[i - 1] = f[i - 1] - (m[i - 1] * h[i - 1] * h[i - 1]) / 6.0;
        B[i - 1] = f[i] - (m[i] * h[i - 1] * h[i - 1]) / 6.0;
    }
    printVector(A, "A coefficients");
    printVector(B, "B coefficients");

    //cout << fixed << setprecision(4);
    cout << "s(x) = {\n";
    for (int i = 0; i < n - 1; ++i) {
        cout << "    ";
        cout << m[i] << "* (" << x[i + 1] << " - x)^3 /  (6 * " << h[i] << ")+ "
             << "(" << m[i + 1] << "* (x - " << x[i] << ")^3  / (6 * "<< h[i] << ") + "
             << A[i] << " * (" << x[i + 1] << " - x) / " << h[i] << " + "
             << B[i] << " * (x - " << x[i] << ") / " << h[i];
        cout << ",    if " << x[i] << " <= x <= " << x[i + 1] << ";\n";
    }
    cout << "};\n";

    cout << "s(x) = {\n";
    for (int i = 0; i < n - 1; ++i) {
        double mi = m[i] / (6 * h[i]);
        double mi1 = m[i + 1] / (6 * h[i]);
        cout << "    "
        <<  mi* pow(x[i + 1],3) << " - " << mi*3*pow(x[i + 1],2) << "*x + " << mi*3*x[i + 1] << "*x^2 - "<< mi <<"x^3) + "
             << mi1 << "x^3 - "<< 3*mi1 * x[i] << "x^2 + "<< mi1 * 3* x[i] << "x^2 - " << mi1 * pow(x[i],3) << " + "
             << (A[i] / h[i]) * x[i + 1] << " - "<< (A[i] / h[i]) << "x + "
             << (B[i] / h[i]) << " x - " <<(B[i] / h[i]) * x[i] <<
                 ",    if " << x[i] << " <= x <= " << x[i + 1] << ";\n";
    }
    cout << "};\n";

    vector<SplineSegment> spline_segments(n-1);

    for(int i = 0; i < n-1; ++i) {
        spline_segments[i]={m[i]/(6*h[i]),m[i + 1] / (6 * h[i]), A[i] / h[i], B[i]/h[1], x[i], x[i+1]};
    }

    const int plot_points = 25;
    const double step = (x_max - x_min) / (plot_points);

    ofstream outFile("polynomial_plot_data.txt");
    if (outFile.is_open()) {
        for (double p = x_min; p <= x_max; p += step) {
            double y = calculateSplineValue(p, x, m, A, B);
            double log = func(p);
            outFile << p << " " << y << " " << log;

            for(const SplineSegment& segment : spline_segments) {
                if( p >= segment.x_start && p <=segment.x_end) {
                    double firstDerivative = splineFirstDerivative(p, segment);
                    double secondDerivative = splineSecondDerivative(p, segment);
                    outFile << " " << firstDerivative << " " << secondDerivative << endl;
                }
            }
        }


        outFile.close();
        cout << "Data for graph has been saved to 'polynomial_plot_data.txt'." << endl;
    } else {
        cerr << "Error opening file to write" << endl;
    }

    return 0;
}

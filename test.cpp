#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

// Функція для виводу матриці
void printMatrix(const vector<vector<double>>& matrix, const string& name) {
    cout << name << ":\n";
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << setw(5) << setprecision(4) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// Функція для виводу вектора
void printVector(const vector<double>& vec, const string& name) {
    cout << name << ":\n";
    for (double val : vec) {
        cout << setw(5) << setprecision(4) << val << " ";
    }
    cout << endl << endl;
}

int main() {
    // Вузли x і значення функції f(x)
    /*vector<double> x = {1, 8, 15, 22, 29, 36, 43, 50, 57, 64, 71, 78, 85, 92, 99};
    vector<double> f = {0, 0.9031, 1.1761, 1.3424, 1.4624, 1.5563, 1.6335, 1.6989, 1.7559, 1.8062, 1.8513, 1.8921, 1.9294, 1.9638, 1.9956};*/
    vector<double> x,f;
    /*
    for (double i = 0.0; i<=3.0; i++) {
        x.push_back(i);
    }
    f.push_back(0.0);
    f.push_back(0.5);
    f.push_back(2.0);
    f.push_back(1.5);
    */
    x = {-2,1,3};
    f = {3, 1, 0};
    int n = x.size();

    // Кроки h_i
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }
    printVector(h, "h (steps)");

    // Формування тридіагональної матриці C і вектора b
    vector<vector<double>> C(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);

    for (int i = 1; i < n - 1; ++i) {
        C[i][i - 1] = h[i - 1];
        C[i][i] = 2 * (h[i - 1] + h[i]);
        C[i][i + 1] = h[i];
        b[i] = 6 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
    }

    // Граничні умови (природний сплайн)
    C[0][0] = 1;
    C[n - 1][n - 1] = 1;

    printMatrix(C, "C (tridiagonal matrix)");
    printVector(b, "b (right-hand side)");

    // Розв'язання системи C * m = b (метод прогонки)
    vector<double> m(n, 0.0); // Результуючий вектор m
    vector<double> alpha(n, 0.0), beta(n, 0.0);

    cout << "Direct pass:\n";
    for (int i = 1; i < n; ++i) {
        double t = C[i][i - 1] / C[i - 1][i - 1];
        C[i][i] -= t * C[i - 1][i];
        b[i] -= t * b[i - 1];
        cout << "Step " << i << ":\n";
        printMatrix(C, "C (after modification)");
        printVector(b, "b (after modification)");
    }

    // Зворотній хід
    cout << "\nReverse pass:\n";
    m[n - 1] = b[n - 1] / C[n - 1][n - 1];
    cout << "m[" << n - 1 << "] = " << m[n - 1] << endl;
    for (int i = n - 2; i >= 0; --i) {
        m[i] = (b[i] - C[i][i + 1] * m[i + 1]) / C[i][i];
        cout << "m[" << i << "] = " << m[i] << endl;
    }

    printVector(m, "m (second derivatives)");

    // Обчислення коефіцієнтів A і B
    vector<double> A(n - 1), B(n - 1);
    for (int i = 1; i < n; ++i) {
        A[i - 1] = f[i - 1] - (m[i - 1] * h[i - 1] * h[i - 1]) / 6.0;
        B[i - 1] = f[i] - (m[i] * h[i - 1] * h[i - 1]) / 6.0;
    }
    printVector(A, "A coefficients");
    printVector(B, "B coefficients");

    cout << fixed << setprecision(4);
    cout << "s(x) = {\n";
    for (int i = 0; i < n - 1; ++i) {
        cout << "    ";
        cout << m[i]/(6 *h[i]) << "* (" << x[i+1] << " - x)^3 + "
             << "(" << m[i+1]/(6*h[i])<< "(x - " << x[i] << ")^3 + "
             << A[i]/h[i] << " * (" << x[i+1] << " - x) + "
             << B[i]/ h[i] << " * (x - " << x[i] << ")";
        cout << ",    if " << x[i] << " <= x <= " << x[i+1] << ";\n";
    }
    cout << "};\n";
    return 0;
}

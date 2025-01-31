#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Функція для розрахунку логарифма
double func(double x) {
    return log10(x);
}

// Коефіцієнти сплайна
struct Spline {
    double a, b, c, d, x;
};
vector<double> solveTridiagonal(const vector<double>& h, const vector<double>& alpha, int n) {
    vector<double> c(n + 1, 0); // Масив коефіцієнтів c
    vector<double> l(n + 1), mu(n + 1), z(n + 1);

    // Природні умови: c0 = cn = 0
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n; i++) {
        l[i] = 2 * (h[i - 1] + h[i]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    // Зворотний хід
    for (int j = n - 1; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j + 1];
    }

    return c;
}

// Розрахунок кубічного сплайна
vector<Spline> cubicSpline(const vector<double>& x, const vector<double>& y) {
    int n = x.size() - 1;
    vector<double> h(n), alpha(n+1,0), b(n), d(n);
    vector<Spline> splines(n);

    // Обчислення h[i]
    for (int i = 0; i < n; i++) {
        h[i] = x[i + 1] - x[i];
        std::cout << h[i] << std::endl;
    }
    // Формування правої частини alpha
    for (int i = 1; i < n; i++) {
        alpha[i] = (6/ h[i]) * (y[i + 1] - y[i]) - (6 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    auto c = solveTridiagonal(h, alpha, n);
    // Вивід матриці СЛАР
    cout << "СЛАР для коефіцієнтів c[i]:\n";
    for (int i = 1; i < n; i++) {
        cout << "h[" << i - 1 << "] = " << h[i - 1] << ", h[" << i << "] = " << h[i]
             << ", alpha[" << i << "] = " << alpha[i] << endl;
    }

    for (int i = 0; i < n; i++) {
        double value = (c[i + 1] - c[i]) / (h[i]);
        d[i] = value;
    }
    for (int i = 0; i < n; i++) {
        double value = (y[i + 1] - y[i]) / h[i] - (h[i]*h[i])/6 * d[i] + h[i]*c[i+1]/2;
        std::cout << "(" << y[i+1] << " - " << y[i] << ") / " << h[i] << " - (" << h[i] <<")^2 / 6 * " << d[i] << " + " << h[i] << " * " << c[i+1] <<" / 2 = " << h[i]*c[i+1] << "/2 - " <<  (h[i]*h[i])* d[i] << "/6 + " <<  (y[i + 1] - y[i]) / h[i] << std::endl;
        b[i] = value;
    }
    for (int i = 0; i < n; i++) {
        std::cout << "b[" << i << "] = " << b[i] << std::endl;
        splines[i] = {y[i+1], b[i], c[i], d[i], x[i]};
    }

    return splines;
}

// Значення сплайна в точці
double evaluateSpline(const vector<Spline>& splines, double x) {
    for (const auto& s : splines) {
        if (x >= s.x && x <= s.x + 7) {
            double dx = x - s.x;
            return s.a + s.b * dx + s.c * dx * dx + s.d * dx * dx * dx;
        }
    }
    return 0; // Зовні інтервалу
}

int main() {
    // Вузли з кроком 7
    vector<double> x, y;
    /*
    for (double i = 1; i <= 99; i += 7) {
        x.push_back(i);
        y.push_back(func(i));
    }
    */
    for (double i = 0.0; i<=3.0; i++) {
        x.push_back(i);
    }
    y.push_back(0.0);
    y.push_back(0.5);
    y.push_back(2.0);
    y.push_back(1.5);


    // Кубічний сплайн
    auto splines = cubicSpline(x, y);

    // Вивід коефіцієнтів сплайна
    cout << "\nКоефіцієнти сплайна:\n";
    for (const auto& s : splines) {
        cout << "x = " << s.x << ": a = " << s.a << ", b = " << s.b
             << ", c = " << s.c << ", d = " << s.d << endl;
    }

    // Обчислення та вивід графіка
    /*
    cout << "\nГрафік сплайна та оригінальної функції:\n";
    for (double xi = 1; xi <= 99; xi += 0.5) {
        double s_val = evaluateSpline(splines, xi);
        cout << "x = " << xi << ", S(x) = " << s_val
             << ", log10(x) = " << func(xi) << endl;
    }
    */

    return 0;
}

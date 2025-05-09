#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// Kita normalisasi tahun dengan mengurangi 2000 agar tidak terlalu besar
const int YEAR_SHIFT = 2000;

// Fungsi untuk menghitung regresi polinomial derajat 3
tuple<double, double, double, double> regresiPolinomial(const vector<int>& x_raw, const vector<double>& y) {
    int n = x_raw.size();
    vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = x_raw[i] - YEAR_SHIFT;

    double sum_x[7] = {0};
    double sum_xy[4] = {0};

    for (int i = 0; i < n; ++i) {
        double xi = x[i];
        double yi = y[i];
        for (int j = 0; j <= 6; ++j)
            sum_x[j] += pow(xi, j);
        for (int j = 0; j <= 3; ++j)
            sum_xy[j] += yi * pow(xi, j);
    }

    double A[4][5] = {
        {sum_x[0], sum_x[1], sum_x[2], sum_x[3], sum_xy[0]},
        {sum_x[1], sum_x[2], sum_x[3], sum_x[4], sum_xy[1]},
        {sum_x[2], sum_x[3], sum_x[4], sum_x[5], sum_xy[2]},
        {sum_x[3], sum_x[4], sum_x[5], sum_x[6], sum_xy[3]},
    };

    // Eliminasi Gauss Jordan
    for (int i = 0; i < 4; ++i) {
        double pivot = A[i][i];
        for (int j = 0; j < 5; ++j)
            A[i][j] /= pivot;
        for (int k = 0; k < 4; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < 5; ++j)
                    A[k][j] -= factor * A[i][j];
            }
        }
    }
    return {A[0][4], A[1][4], A[2][4], A[3][4]};
}

double evaluasiPolinomial(double x_raw, double a, double b, double c, double d) {
    double x = x_raw - YEAR_SHIFT;
    return a + b * x + c * x * x + d * x * x * x;
}

int main() {
    vector<int> tahun_pop = {1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969,
        1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985,
        1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001,
        2002, 2003, 2004, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2017, 2018, 2019, 2020, 2021, 2022, 2023};

    vector<double> populasi = {88296070, 90791249, 93375850, 96051424, 98833749, 101365130, 103792754, 106526393,
        109450006, 112517639, 115657495, 118833702, 122039841, 125288515, 128555045, 131843848, 135173655,
        138533541, 141953163, 145434834, 148950540, 152485035, 156052152, 159651381, 163251124, 166776185,
        170175065, 173511154, 176855065, 180201630, 183501098, 186778238, 190043744, 193305168, 196591828,
        199888057, 203204348, 206536095, 209826788, 213004668, 216077790, 219097902, 222088495, 225048008,
        227926649, 237062337, 240157903, 243220028, 246305322, 249470032, 252698525, 255852467, 258877399,
        267346658, 269951846, 272489381, 274814866, 276758053, 278830529, 281190067};

    vector<int> tahun_net = {1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003,
        2004, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2017, 2018, 2019, 2020, 2021, 2022, 2023};

    vector<double> persentase = {0.001059744, 0.026109477, 0.056623989, 0.194910264, 0.255306646, 0.444415936,
        0.925564, 2.01861, 2.13414, 2.38702, 2.60029, 5.78627, 7.91748, 6.92, 10.92, 12.28,
        14.52, 14.94, 17.1432, 32.3358, 39.9046, 47.6906, 53.7265, 62.1045, 66.4846, 69.2084};

    // Regresi populasi dan internet
    auto [a_p, b_p, c_p, d_p] = regresiPolinomial(tahun_pop, populasi);
    auto [a_i, b_i, c_i, d_i] = regresiPolinomial(tahun_net, persentase);

    vector<int> estimasi_tahun = {2005, 2006, 2015, 2016, 2030, 2035};

    cout << fixed << setprecision(2);
    cout << "\n=== Estimasi Populasi dan Internet ===\n";
    for (int t : estimasi_tahun) {
        double pop = evaluasiPolinomial(t, a_p, b_p, c_p, d_p);
        double pct = evaluasiPolinomial(t, a_i, b_i, c_i, d_i);

        // Tambahkan batas atas realistis untuk presentase jika diinginkan:
        // pct = min(pct, 100.0);

        double pengguna = pop * pct / 100.0;

        cout << "Tahun " << t << ":\n";
        cout << "  Populasi: " << (long long)pop << "\n";
        cout << "  % Internet: " << pct << "%\n";
        cout << "  Jumlah Pengguna Internet: " << (long long)pengguna << "\n\n";
    }

    // Catatan penting untuk laporan:
    cout << "\nCatatan: Model regresi polinomial tidak memiliki batas atas sehingga pada tahun-tahun jauh seperti 2030 dan 2035,\n"
            "persentase pengguna internet bisa lebih dari 100%. Ini mencerminkan keterbatasan model dan menunjukkan bahwa\n"
            "untuk prediksi jangka panjang lebih baik menggunakan model yang memiliki batas atas seperti regresi logistik.\n";

    // Cetak persamaan
    cout << "\nPersamaan Polinomial Populasi:\n";
    cout << "y = " << a_p << " + " << b_p << "(x-" << YEAR_SHIFT << ") + " << c_p << "(x-" << YEAR_SHIFT << ")^2 + " << d_p << "(x-" << YEAR_SHIFT << ")^3\n";
    cout << "\nPersamaan Polinomial % Internet:\n";
    cout << "y = " << a_i << " + " << b_i << "(x-" << YEAR_SHIFT << ") + " << c_i << "(x-" << YEAR_SHIFT << ")^2 + " << d_i << "(x-" << YEAR_SHIFT << ")^3\n";

    return 0;
}
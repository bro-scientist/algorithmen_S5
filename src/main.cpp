#include "algorithms.cpp"

using namespace std;

float func_to_integrate(float x) {
    return pow(x, 2);
}

float func_to_dgl(float x, float y) {
    return 2 * x * y;
}

int main()
{
    // serie 1 - root & primes

#if false
        int a, n;
        cout << "A:";
        cin >> a;
        cout << endl << "N:";
        cin >> n;
        cout << endl << n << ". Wurzel von " << a << " ist: " << nth_root(n, a) << endl;
#endif

#if false
        int n;
        cout << "primes to:";
        cin >> n;
        int* numbers = primes_to_n(n);

        cout << "-Primzahlen bis " << n << "-" << endl;
        for(int i = 0; i < n; i++)
        {
            if(numbers[i]) cout << numbers[i] << endl;
        }
#endif

#if false
        cout << "divisor: 432" << endl;
        cout << "divisor: " << biggest_common_divisor(38016, 2160) << endl;
        cout << "divisor: " << biggest_common_divisor(42768, 2534112) << endl;
#endif

    // serie 2 - pi & karazuba

#if false
        cout << "pi:     " << pi(256) << endl;
        cout << "better: " << better_pi(INT8_MAX) << endl;
        cout << "        " << better_pi(INT16_MAX) << endl;
        cout << "        " << better_pi(INT32_MAX) << endl;
#endif

#if false
        cout << "karazuba: 7006652" << endl;
        cout << "karazuba: " << karazuba_pos_even(24, 18) << endl;
#endif

#if false
        cout << "multiply: 432" << endl;
        cout << "multiply: " << multiply(24, 18) << endl;
#endif

    // serie 3 - differentiation & integral

#if false
    cout << "differentiate: 8" << endl;
    cout << "differentiate: " << differentiate(&func_to_integrate, 4, 0.01f) << endl;
    cout << "differentiate2: 2" << endl;
    cout << "differentiate2: " << differentiate2(&func_to_integrate, 400, 1.0f) << endl;
#endif

#if false
    float a = 15.0f, b = 47.0f;
    float deltaX = 0.8f, intervals = 8;

    float trapez = integrate_by_trapezoid(&func_to_integrate, a, b, deltaX);
    float shorter = integrate_shorter(&func_to_integrate, a, b, deltaX);
    float simpson = integrate_simpson(&func_to_integrate, a, b, intervals);

    cout << "integrate_by_trapez: " << trapez << endl;
    cout << "integrate_shorter: " << shorter << endl;
    cout << "integrate_simpson: " << simpson << endl;

    assert (round(trapez) == round(shorter));
#endif

    // serie 4 - regression & dings

    float temperatures[] = {
            1905, 7.7F,
            1915, 8.0F,
            1925, 7.9F,
            1935, 8.1F,
            1945, 8.3F,
            1955, 8.1F,
            1965, 7.9F,
            1975, 8.3F,
            1985, 8.5F,
            1995, 9.0F,
            2005, 9.2F
    };

#if false
    cout << "regression_get: ~y = -15.365 + 0.012x" << endl;
    regression_pair pair = regression_get(temperatures, 11);
    cout << "regression_get: ~y = " << pair.a << " + " << pair.b << "x" << endl;
    cout << "regression_predict: 2010: 8.9" << endl;
    cout << "regression_predict: 2010: " << regression_predict(temperatures, 11, 2010) << endl;
#endif

    // serie 5 - dgl

#if false
    float x = 4.0F, y = 1.0F, deltaX = 0.01, limit = 4.0F;
    cout << "         euler_cauchy: " << euler_cauchy(func_to_dgl, x, y, deltaX, limit) << endl;
    cout << "         euler_better: " << euler_better(func_to_dgl, x, y, deltaX, limit) << endl;
    cout << "runge_kutta_2nd_order: " << runge_kutta_2nd_order(func_to_dgl, x, y, deltaX, limit) << endl;
    cout << "runge_kutta_4th_order: " << runge_kutta_4th_order(func_to_dgl, x, y, deltaX, limit) << endl;
#endif

    // serie 6.1 - lgs

    const int size0 = 3;
    int input0[size0 * size0] = {
            3, 2, 1,
            6, 6, 3,
            9,10, 6
    };
    int b0[size0] = {2, 3, 5};

    const int size1 = 5;
    int input1[size1 * size1] = {
            4, 0, 3, 2, 7,
            8, 2, 9, 2, 1,
            1, 3, 4, 1, 3,
            2, 2, 4, 6, 8,
            6, 8, 8, 5, 0
    };
    int b1[size1] = {3, 1, 3, 7, 4};

#if false
    cout << "lr_zerlegung: " << endl;
    // pivotize_rows(input0, b0, size0);
    lr_zerlegung(input0, b0, size0);

    // pivotize_rows(input1, b1, size1);
    lr_zerlegung(input1, b1, size1);
#endif

    // serie 6.2 - lgs 2

    const int size2 = 3;
    int input2[size2 * size2] = {
            3, 1, 0,
            1, 3, 1,
            0, 1, 3
    };
    int b2[size2] = {1, 5, 7};
    double x[size2] = {0}, x_old[size2] = {0};

#if false
    jacobi_lgs(input2, b2, size2, 0.01, x, x_old);
#endif

#if false
    gau3seidel_lgs(input2, b2, size2, 0.01, x, x_old);
#endif

    return 0;
}
#include "algorithms.cpp"

using namespace std;

float func_to_integrate(float x) {
    return pow(x, 2);
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

#if true
    cout << "regression_get: ~y = -15.365 + 0.012x" << endl;
    regression_pair pair = regression_get(temperatures, 11);
    cout << "regression_get: ~y = " << pair.a << " + " << pair.b << "x" << endl;
    cout << "regression_predict: 2010: 8.9" << endl;
    cout << "regression_predict: 2010: " << regression_predict(temperatures, 11, 2010) << endl;
#endif

    return 0;
}
#include <iostream>

static int nth_root(int n, double a)
{
    int oldx = 10;
    int newx = 1;
    double epsilon = 0.001;

    while (abs(oldx - newx) > epsilon) {
        oldx = newx;
        newx = ((n - 1) * pow(oldx, n) + a) / (n * pow(oldx, n - 1));
    }
    return newx;
}

static int* primes_to_n(int n)
{
    int *numbers = new int[n-1];
    for (int i = 2; i < n+1; i++) {
        numbers[i-2] = i;
    }

    for (int i = 2; i * i < n; i++) {
        int k = 2, ki = k*i;
        while (ki < n) {
            numbers[ki] = 0;
            ki = ++k*i;
        }
    }

    return numbers;
}
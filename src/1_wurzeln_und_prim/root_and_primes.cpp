#include <iostream>

static void arr_pop(int array[], int size, int pos)
{
    int length = size/sizeof(int);
    for (int i = pos; i < length; i++)
    {
        array[i-1] = array[i];
    }
}

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
    int *numbers = new int[n];
    for (int i = 2; i < n; i++) {
        numbers[i] = i;
    }

    for (int i = 2; i * i < n; i++) {
        // for (int k = 2; k * i < n; k++) numbers[k * i] = 0;
        int k = 2, ki = k*i;
        while (ki < n) {
            numbers[ki] = 0; // TODO replace with pop
            ki = ++k*i;
        }
    }

    return numbers;
}
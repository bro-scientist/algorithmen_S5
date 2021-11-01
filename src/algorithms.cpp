#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

// serie 1 - root & primes

static int nth_root(int n, double a)
{
    int oldx = 10;
    int newx = 1;
    double epsilon = 0.001;

    while (abs(oldx - newx) > epsilon)
    {
        oldx = newx;
        newx = ((n - 1) * pow(oldx, n) + a) / (n * pow(oldx, n - 1));
    }
    return newx;
}

static int* primes_to_n(int n)
{
    int *numbers = new int[n-1];
    for (int i = 2; i < n+1; i++)
    {
        numbers[i-2] = i;
    }

    for (int i = 2; i * i < n; i++)
    {
        int k = 2, ki = k*i;
        while (ki < n)
        {
            numbers[ki] = 0;
            ki = ++k*i;
        }
    }

    return numbers;
}

static uint64_t biggest_common_divisor(uint64_t a, uint64_t b)
{

    while (b)
    {
        if (a < b)
        {
            uint64_t swap = a;
            a = b;
            b = swap;
        }
        uint64_t diff = a-b;
        a = b;
        b = diff;
    }
    return a;
}

// serie 2 - pi & karazuba

static float pi(int edges)
{
    float sn = 2.0;

    unsigned n = 1;
    while (n < edges)
    {
        sn = sqrt(2-2*sqrt(1-(pow(sn,2)/4)));
        n <<= 1; // bit shift multiplication
    }
    return sn*n;
}

static float better_pi(int edges)
{
    float sn = 2.0;

    unsigned n = 1;
    while (n < edges)
    {
        sn = sqrt(pow(sn,2)/(2+2*sqrt(1-pow(sn,2)/4)));
        n <<= 1; // bit shift multiplication
    }
    return sn*n;
}

static int highest_bit(unsigned long long number)
{
    int result = -1;
    while(number >> ++result) {}
    return --result;
}

static int karazuba_pos_even(int x, int y)
{
    int b = 10;
    int n = (int)fmax(highest_bit(x), highest_bit(y));

    if (n<4) return x*y;

    int k = n/2 + n%2;
    int bK = pow(b, k);

    int P = x/bK;
    int Q = x%P;
    int E = Q-P;

    int T = y/bK;
    int S = y%T;
    int F = T-S;

    int U = karazuba_pos_even(P, S);
    int V = karazuba_pos_even(abs(E), abs(F));
    int W = karazuba_pos_even(Q, T);

    return U*(int)pow(b,2)+(U+W-E*F*V)*(bK)+W;
}

long long multiply(long long x, long long y)
{
    int xLength = highest_bit(x);
    int yLength = highest_bit(y);

    // the bigger of the two lengths
    int N = (int)(fmax(xLength, yLength));

    // if the max length is small it's faster to just flat out multiply the two nums
    if (N < 10)
        return x * y;

    //max length divided and rounded up
    N = (N/2) + (N%2);

    long long multiplier = pow(10, N);

    long long b = x/multiplier;
    long long a = x - (b * multiplier);
    long long d = y / multiplier;
    long long c = y - (d * N);

    long long z0 = multiply(a,c);
    long long z1 = multiply(a + b, c + d);
    long long z2 = multiply(b, d);


    return z0 + ((z1 - z0 - z2) * multiplier) + (z2 * (long long)(pow(10, 2 * N)));

}

// serie 3 - differentiation & integral

float differentiate(float (*func)(float), float x, float h)
{
    // return (func(x + h) - func(x)) / h;
    return (func(x + h) - func(x - h)) / (h * 2);
}

float differentiate2(float (*func)(float), float x, float h)
{
    return (func(x + h) + func(x - h) - 2 * func(x)) / pow(h, 2);
}

float integrate_by_trapezoid(float (*func)(float), float a, float b, float deltaX)
{
    float stripes = (b - a) / deltaX;
    float sum = 0;
    for(float i = 0; i < stripes; i++)
    {
        float val = func(a + i * deltaX) + func(a + (i + 1) * deltaX);
        sum += val;
    }
    return sum * deltaX / 2;
}

float integrate_shorter(float (*func)(float), float a, float b, float deltaX)
{
    float sum = (func(a) + func(b)) / 2;
    a += deltaX;
    while(a < b - deltaX)
    {
        sum += func(a);
        a += deltaX;
    }
    return sum * deltaX;
}

float integrate_simpson(float (*func)(float), float a, float b, int intervals)
{
    assert(!(intervals % 2));

    float deltaX = 1.0f / (float)intervals;

    float sum = (func(a) + func(b));

    a += deltaX;
    int i = 1;
    while(a < b - deltaX)
    {
        sum += (float)((i++ % 2) ? 4 : 2) * func(a);
        a += deltaX;
    }
    return sum * deltaX / 3;
}

// serie 4 - regression & dings

struct regression_pair
{
    double a,b;
};

regression_pair regression_get(float values[], int length)
{
    double x_ = 0, y_ = 0, a , b, sum_xi_yi = 0, sum_xi_sq = 0;

    for (int i = 0; i < length; i++)
    {
        x_ += values[i << 1];
        y_ += values[(i << 1) + 1];
        sum_xi_yi += values[i << 1] * values[(i << 1) + 1];
        sum_xi_sq += pow(values[i << 1], 2);
    }

    x_ /= length;
    y_ /= length;
    b = (sum_xi_yi - length * x_ * y_) / (sum_xi_sq - length * pow(x_, 2));
    a = y_ - b * x_;

    return {a,b};
}

float regression_predict(float values[], int length, float x) {
    regression_pair params = regression_get(values, length);
    return params.a + params.b * x;
}
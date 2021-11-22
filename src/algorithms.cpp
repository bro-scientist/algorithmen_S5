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

// serie 5 - dgl

float euler_cauchy(float (*func)(float, float), float x, float y, float deltaX, float limit)
{
    limit /= deltaX;
    for(int i = 0; i < limit; i++)
    {
        y = func(x, y) * deltaX + y;
        x += deltaX;
        // std::cout << i << ": " << "|" << y << std::endl;
    }
    return y;
}

float euler_better(float (*func)(float, float), float x, float y, float deltaX, float limit)
{
    limit /= deltaX;
    for(int i = 0; i < limit; i++)
    {
        y = func(x + deltaX / 2, y + (deltaX * func(x, y)) / 2) * deltaX + y;
        x += deltaX;
        // std::cout << i << ": " << "|" << y << std::endl;
    }
    return y;
}

float runge_kutta_2nd_order(float (*func)(float, float), float x, float y, float s, float limit)
{
    float k1, k2;
    limit /= s;
    for(int i = 0; i < limit; i++)
    {
        k1 = func(x, y);
        x += s;
        k2 = func(x, y + s * k1);
        y += (s * (k1 + k2) / 2);
        // std::cout << i << ": " << "|" << y << std::endl;
    }
    return y;
}

float runge_kutta_4th_order(float (*func)(float, float), float x, float y, float s, float limit)
{
    float k1, k2, k3, k4;
    limit /= s;
    for(int i = 0; i < limit; i++)
    {
        k1 = func(x, y);
        k2 = func(x + (s / 2), y + s * k1 / 2);
        k3 = func(x + (s / 2), y + s * k2 / 2);
        x += s;
        k4 = func(x, y + s * k3);
        y += s * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        // std::cout << i << ": " << x << "|" << y << std::endl;
    }
    return y;
}

// serie 6.1 - lgs

void print_arr(int* arr, int size)
{
    for (int x = 0; x < size; x++)
    {
        for (int y = 0; y < size; y++)
        {
            cout << " " << arr[x * size + y];
        }
        cout <<  endl;
    }
    cout <<  endl;
}

void lr_zerlegung(int* input, int* b, int size)
{
    // x = v  ;  y = ->

    int* left;
    left = new int[size * size];

    for (int x = 0; x < size; x++)
    {
        for (int y = 0; y < size; y++)
        {
            left[x * size + y] = x == y;
        }
    }

    for (int i = 0; i < size - 1; i++)
    {
        for (int x = i + 1; x < size; x++)
        {
            left[x * size + i] = input[x * size + i] / input[i * size + i];
            for (int y = 0; y < size; y++)
            {
                input[x * size + y] -= left[x * size + i] * input[i * size + y];
            }
        }
    }

    cout << "left:" << endl;
    print_arr(left, size);

    cout << "right:" << endl;
    print_arr(input, size);

    cout << "result:" << endl;
    int* y = new int[size];

    y[0] = b[0];
    y[1] = b[1] - left[1 * size + 0] * y[0];
    y[2] = b[2] - left[2 * size + 0] * y[0] - left[2 * size + 1] * y[1];
    // TODO

    for (int x = 0; x < size; x++)
    {
        cout << y[x] << endl;
    }
}

void pivotize_rows(int* input, int* b, int size)
{
    // x = ->  ;  y = v

    for (int x = 0; x < size; x++)
    {
        int temp = INT32_MIN, row = 0;
        for (int y = x; y < size; y++)
        {
            if (input[y * size + x] > temp)
            {
                temp = input[y * size + x];
                row = y;
            }
        }

        temp = b[x];
        b[x] = b[row];
        b[row] = temp;

        for (int xx = 0; xx < size; xx++)
        {
            temp = input[x * size + xx];
            input[x * size + xx] = input[row * size + xx];
            input[row * size + xx] = temp;
        }
    }
}

void pivotize_cols(int* input, int size)
{
    // x = ->  ;  y = v

    for (int y = 0; y < size; y++)
    {
        int temp = INT32_MIN, col = 0;
        for (int x = y; x < size; x++)
        {
            if (input[y * size + x] > temp)
            {
                temp = input[y * size + x];
                col = x;
            }
        }

        for (int yy = 0; yy < size; yy++)
        {
            temp = input[yy * size + y];
            input[yy * size + y] = input[yy * size + col];
            input[yy * size + col] = temp;
        }
    }
}

// serie 6.2 - lgs 2

void jacobi_or_gau3seidel_lgs(int* input, int* b, int size, double epsilon, double* x, double* x_old, bool overwrite)
{
    double error;

    int iter = 0;
    do
    {
        error = 0;
        cout << ++iter << ". interation ";
        for (int i = 0; i < size; i++)
        {
            double sum = 0;
            for (int j = 0; j < size; j++)
            {
                if (i != j) sum += input[i * size + j] * (overwrite ? x[j] : x_old[j]);
            }
            x[i] = (b[i] - sum) / input[i * size + i];
            cout << "x" << i << ": " << x[i] << ", ";
        }
        cout << endl;

        for (int k = 0; k < size; k++)
        {
            error += abs(x_old[k] - x[k]);
            x_old[k] = x[k];
        }
    }
    while (error > epsilon);

    cout << endl << "result: " << endl;
    for (int r = 0; r < size; r++)
    {
        cout << "x" << r << ": " << round(x[r]) << ", " << endl;
    }
}

void jacobi_lgs(int* input, int* b, int size, double epsilon, double* x, double* x_old)
{
    jacobi_or_gau3seidel_lgs(input, b, size, epsilon, x, x_old, false);
}

void gau3seidel_lgs(int* input, int* b, int size, double epsilon, double* x, double* x_old)
{
    jacobi_or_gau3seidel_lgs(input, b, size, epsilon, x, x_old, true);
}

// serie 6.3
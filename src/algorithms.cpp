#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <set>

#define INFINITY INT32_MAX

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
        n *= 2;
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
        n *= 2;
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

void print_arr_2d_char_offset(int* arr, int sizeX, int sizeY, int offset)
{
    for (int x = 0; x < sizeX; x++)
    {
        for (int y = 0; y < sizeY; y++)
        {
            int value = arr[x * sizeX + y];
            cout << (value != -1 ? (char)(value + offset) : '#') << "|";
        }
        // move cursor back and replace last char with space
        cout << '\b' << ' ';
        cout << endl;
    }
}

void print_arr_2d(int* arr, int sizeX, int sizeY)
{
    for (int x = 0; x < sizeX; x++)
    {
        for (int y = 0; y < sizeY; y++)
        {
            int value = arr[x * sizeX + y];
            if (value < INFINITY) cout << value << '|';
            else cout << '\236' << '|';
        }
        // move cursor back and replace last char with space
        cout << '\b' << ' ';
        cout <<  endl;
    }
}

void print_arr_2d(int* arr, int size)
{
    print_arr_2d(arr, size, size);
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
    print_arr_2d(left, size);

    cout << "right:" << endl;
    print_arr_2d(input, size);

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

void abstiegsverfahren_dumm(float func(float, float), float step, float d0(float), float d1(float))
{
    int dim = 2;
    float diff = 1, epsilon = 0.01;
    float x_last[2], x[2] = {1, 1};

    while (diff > epsilon)
    {
        /*
        for (int i = 0; i < dim; i++)
        {
            x_last[i] = x[i];
            x[i] = x_last[i] + step *
        }
        */

        x_last[0] = x[0];
        x_last[1] = x[1];

        x[0] = x_last[0] + step * d0(x[0]);
        x[1] = x_last[1] + step * d1(x[1]);

        std::cout << "x: (" << x[0] << " , " << x[1] << ") " << std::endl;

        diff = abs(func(x_last[0], x_last[1]) - func(x[0], x[1]));
        std::cout << diff << std::endl;
    }
}


void abstiegsverfahren(float func(float[]), float step, float (*d[])(float), int dimensions, float x_start[], float x_last[])
{
    float diff = 1, epsilon = 0.01;

    while (diff > epsilon)
    {
        std::cout << "x: (";
        for (int i = 0; i < dimensions; i++)
        {
            x_last[i] = x_start[i];
            x_start[i] = x_last[i] + step * d[i](x_start[i]);
            std::cout << x_start[i] << ",";
        }

        diff = abs(func(x_last) - func(x_start));
        // delete last comma
        std::cout << "\x08) diff: " << diff << std::endl;
    }
}

// serie 6.3





// serie 10

int simulated_annealing(double epsilon, double t, int size, int A[], int queue[], int best[])
{
    double delta;
    int shortest = INFINITY;

    srand(time(nullptr));

    while (t > epsilon) {

        // calculate length of current queue
        int length = 0;
        for (int i = 0; i < size; i++)
        {
            length += A[queue[i] * size + queue[i + 1]];
        }
        delta = abs(length - shortest);

        // keep new queue as best if round-trip is shorter OR random value is smaller than p
        if (length < shortest || (double)(rand() % 1000) / 1000 < exp(-(delta/t)))
        {
            shortest = length;
            for (int i = 0; i < size+1; i++) best[i] = queue[i];

            cout << "new best: ";
            for (int i = 0; i < size+1; i++)
            {
                cout << best[i];
            }
            cout << " - " << shortest << endl;
        }

        // reverse part of the current queue
        int from = rand() % 6, to = rand() % 6;
        while(from >= to) from = rand() % 6, to = rand() % 6;
        while (from < to)
        {
            int mem = queue[from];
            queue[from++] = queue[to];
            queue[to--] = mem;
        }

        // match beginning and ending of euler-round-trip
        queue[6] = queue[0];

        // "cool down"
        t /= 2;
    }

    // output final result and return
    cout << "best:     ";
    for (int i = 0; i < size+1; i++)
    {
        cout << best[i];
    }
    cout << " - " << shortest << endl;

    return shortest;
}

// serie 11

// TODO: visualizeMST();
void prim_MST(int A[], const int size)
{
    int mstLength = 0;
    vector<int> nodesInMst = *new vector<int>();

    nodesInMst.push_back(0);
    while (nodesInMst.size() < size) // can be a for loop with size-1 iterations
    {
        int smallestValue = INFINITY, from, to;
        for (int node : nodesInMst)
        {
            for (int i = 0; i < size; i++)
            {
                if (std::find(nodesInMst.begin(), nodesInMst.end(), i) != nodesInMst.end()) continue;
                int currentValue = A[size * node + i];
                if (currentValue < smallestValue)
                {
                    smallestValue = currentValue;
                    from = node;
                    to = i;
                }
            }
        }
        nodesInMst.push_back(to);
        mstLength += smallestValue;
        cout << endl << "link: " << (char)(from + 65) << "->" << (char)(to + 65) << "  | nodesInMst:";
        for (int node : nodesInMst)
        {
            cout << " " << (char)(node + 65);
        }
    }
    cout << endl << "length: " << mstLength;
}

void dijkstra(int A[], const int size, int startNode, int distances[], int predecessors[], const int maxVal)
{
    int char_offset = 65;
    int v = 0, w;
    set<int> B = {};

    // initKW
    for (int i = 0; i < size; i++)
    {
        distances[i] = maxVal;
        predecessors[i] = -1;
    }
    distances[startNode] = 0;

    B.insert(startNode);

    while (!B.empty())
    {
        cout << "Distanzen:   ";
        print_arr_2d(distances, 1, size);
        cout << "Vorgaenger:  ";
        print_arr_2d_char_offset(predecessors, 1, size, char_offset);

        cout << "B: ";
        int dist = maxVal;
        for (int value : B)
        {
            cout << (char)(value + char_offset) << ',';
            if (dist > distances[value])
            {
                v = value;
                dist = distances[value];
            }
        }
        B.erase(v);
        cout << '\b' << ' ';
        cout << endl << "min: " << (char)(v + char_offset) << endl<< "________________________________________" << endl;

        for (w = 0; w < size; w++)
        {
            if (A[v * size + w] < maxVal)
            {
                if (distances[w] == maxVal)
                    B.insert(w);

                // verkuerze
                int newDist = distances[v] + A[v * size + w];
                if (newDist < distances[w])
                {
                    distances[w] = newDist;
                    predecessors[w] = v;
                }
            }
        }
    }
    cout << "Distanzen:   ";
    print_arr_2d(distances, 1, size);
    cout << "Vorgaenger:  ";
    print_arr_2d_char_offset(predecessors, 1, size, char_offset);
}

// serie 12

int getSteps(char c, char* steps, unsigned const int length)
{
    if (!c) return 0;
    for (int i = 0; i < length; i++)
    {
        if (steps[i] == c)
        {
            // cout << "found:" << c << " at " << i << endl;
            return i + 1;
        }
    }
    return 0;
}

unsigned int findText_broyer_moore_horspool(char* text, unsigned const int textLength, char* pattern, unsigned const int patternLength, char* steps)
{
    // run through pattern to initialize steps for char occurrences
    for (int i = patternLength-2; i >= 0; i--) {
        int found = getSteps(pattern[i], steps, patternLength);
        if (!found) steps[patternLength-2-i] = pattern[i];
    }
    if (!getSteps(pattern[patternLength-1], steps, patternLength))
        steps[patternLength-1] = pattern[patternLength-1];

    for (int k = 0; k < patternLength; k++)
        cout << (steps[k] ? steps[k] : '-') << " ";

    cout << endl << pattern << endl;

    unsigned int textPos = patternLength-1;
    while (textPos < textLength)
    {
        for (int diff = 0; diff < patternLength; diff++)
        {
            unsigned int tempPatternPos = patternLength-1-diff, tempTextPos = textPos-diff;
            if (pattern[tempPatternPos] != text[tempTextPos])
            {
                // go forward n steps (if zero is returned, the distance is actually 6)
                unsigned int moveBy = getSteps(text[textPos], steps, patternLength);
                moveBy = moveBy ? moveBy : 6;
                textPos += moveBy;
                break;
            }
            else if (!tempPatternPos)
                return tempTextPos;
        }
    }
    return 0;
}

// serie 14

int det(int x1, int y1, int x2, int y2)
{
    int result = x1*y2 - x2*y1;
    return !result ? 0 : result > 0 ? 1 : -1;
}

int det(int* points, int p, int q, int i)
{
    return det(points[q+0]-points[p+0], points[q+1]-points[p+1], points[i+0]-points[p+0], points[i+1]-points[p+1]);
}

vector<int> jarvis_march_convex_hull(int* points, int size)
{
    const int x = 0;
    const int y = 1;

    // output
    vector<int> convexHull = *new vector<int>;

    // point indices
    int l = 0, p = 0, q = 0;

    // compare y value of each point to find bottom left
    for (int i = 2; i < size*2; i+=2)
        if (points[i+y] < points[l+y] || points[i+y] == points[l+y] && points[i+x] < points[l+x])
            l = i;

    p = l;
    convexHull.push_back(p);
    cout << p/2 << " added." << endl;
    do
    {
        q = p+2;
        for (int i = 0; i < size*2; i+=2)
        {
            // only use when there are a lot of points
            // if (std::find(convexHull.begin(), convexHull.end(), i) != convexHull.end()) continue;
            if (q == i || p == i) continue;

            cout << " det: " << p/2 << "->" << q/2 << "  " << p/2 << "->" << i/2 << " = " << det(points, p, q, i) << endl;
            int determinant = det(points, p, q, i);

            // check if vector p->i is right from p->q or they have the same direction but p->i is longer
            if (determinant < 0 || !determinant && points[q+x]-points[p+x] < points[i+x]-points[p+x])
                q = i;
        }
        convexHull.push_back(q);
        p = q;
        cout << q/2 << " added." << endl;
    }
    while (p != l);

    return convexHull;
}
#include "algorithms.cpp"

using namespace std;

float func_to_integrate(float x) {
    return pow(x, 2);
}

float func_to_dgl(float x, float y) {
    return 2 * x * y;
}

float func_to_find_min(float x[]) {
    return pow(x[0], 2) + pow(x[1], 2);
}

float derivative_x(float x) {
    return -2 * x;
}

float derivative_y(float y) {
    return -2 * y;
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

#if false
    const int dimensions_ = 2;
    float x_last[dimensions_], x_start[dimensions_] = {1, 1};

    float (*derivatives[dimensions_])(float) = {derivative_x, derivative_y};

    abstiegsverfahren(func_to_find_min, 0.3F, derivatives, dimensions_, x_start, x_last);
#endif

#if false

    double epsilon = 0.00001, t = 10000;
    int size = 6, A[] = {
            1000, 3, 5, 6, 7, 4,
            3, 1000, 6, 1, 4, 9,
            5, 6, 1000, 4, 10, 8,
            6, 1, 4 ,1000, 7, 3,
            7, 4, 10, 7, 1000, 5,
            4, 9, 8, 3, 5, 1000
    };

    int queue[] = {0, 1, 2, 3, 4, 5, 0},
        best[] = {0, 1, 2, 3, 4, 5, 0};

    int shortest = simulated_annealing(epsilon, t, size, A, queue, best);

    cout << "nodes:    ";
    for (int i = 0; i < size+1; i++)
    {
        cout << best[i] + 1;
    }
#endif

#if false
    const int size = 6;
    int maxVal = INFINITY, A[] = {
            maxVal, 3, maxVal, 6, maxVal, 3,
            3, maxVal, 2, 5, 5, 1,
            maxVal, 2, maxVal, 1, 9, 6,
            6, 5, 1, maxVal, 2, 4,
            maxVal, 5, 9, 2, maxVal, 1,
            3, 1, 6, 4, 1, maxVal
    };

    prim_MST(A, size);
#endif

    // dijkstra

#if false
    const int size = 7, maxVal = INFINITY;
    int A[] = {
            maxVal, 4     , maxVal, maxVal, maxVal, 10    , 5     ,
            4     , maxVal, 7     , maxVal, maxVal, maxVal, 10    ,
            maxVal, 7     , maxVal, 12    , maxVal, maxVal, 1     ,
            maxVal, maxVal, 12    , maxVal, 4     , maxVal, maxVal,
            maxVal, maxVal, maxVal, 4     , maxVal, 3     , 8     ,
            10    , maxVal, maxVal, maxVal, 3     , maxVal, 4     ,
            5     , 10    , 1     , maxVal, 8     , 4     , maxVal
    };
    int distances[size] = { 0 }, predecessors[size] = { 0 };

    // print_arr_2d(A, size);

    dijkstra(A, size, 0, distances, predecessors, maxVal);
#endif

#if false

    char* file = "C:\\Users\\jacxt\\CLionProjects\\algorithmen_S5\\src\\heuhaufen.txt";
    int* size;
    FILE* fp;
    long lSize;
    char* buffer;
    fp = fopen(file, "rb");
    if(!fp) perror(file),exit(1);

    fseek(fp, 0L, SEEK_END);
    lSize = ftell(fp);
    rewind(fp);

    /* allocate memory for entire content */
    buffer = (char*)calloc(1, lSize+1);
    if(!buffer) fclose(fp),fputs("memory alloc fails",stderr),exit(1);

    /* copy the file into the buffer */
    if( 1!=fread(buffer, lSize, 1, fp))
        fclose(fp),free(buffer),fputs("entire read fails",stderr),exit(1);

    /* do your work here, buffer is a string contains the whole text */
    size = (int*)lSize;
    fclose(fp);



    char text[] = "Eine Nadel im Heuhaufen finden", pattern[] = "banane", pattern2[] = "finden", pattern3[] = "Heuhau";
    unsigned const int textLength = 30, patternLength = 6;
    char steps[patternLength] = { 0 }, steps2[patternLength] = { 0 }, step3[patternLength] = { 0 };

    unsigned int pos1 = findText_broyer_moore_horspool(text, textLength, pattern, patternLength, steps);
    if (pos1)
        cout << ">" << pattern << "< gefunden von " << pos1 << " bis " << pos1+patternLength << "in:\n " << text << endl << endl;
    else
        cout << ">"<< pattern << "< nicht gefunden in:\n " << text << endl<< endl;

    unsigned int pos2 = findText_broyer_moore_horspool(text, textLength, pattern2, patternLength, steps2);
    cout << ">" << pattern2 << "< gefunden von " << pos2 << " bis " << pos2+patternLength << " in:\n " << text << endl << endl;

    unsigned int pos3 = findText_broyer_moore_horspool(buffer, (unsigned int)size, pattern3, patternLength, step3);
    if (pos3)
        cout << ">" << pattern3 << "< gefunden von " << pos3 << " bis " << pos3+patternLength << " in:\n " << buffer << endl << endl;
    else
        cout << ">"<< pattern << "< nicht gefunden in " << buffer << endl<< endl;
#endif

#if false
    const int size = 8;
    int points[size * 2] = {
            2, 0,
            4, 4,
            0, 3,
            2, 9,
            6, 6,
            5, 1,
            8, 2,
            9, 7
    };

    vector<int> convexHull = jarvis_march_convex_hull(points, size);

    cout << "convex hull:" << endl;
    for (int member : convexHull)
        cout << " " << member/2;
#endif

}
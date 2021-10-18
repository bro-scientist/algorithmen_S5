#include "algorithms.cpp"

using namespace std;

int main()
{
    // serie 1 - root & primes

    if(false)
    {
        int a, n;
        cout << "A:";
        cin >> a;
        cout << endl << "N:";
        cin >> n;
        cout << endl << n << ". Wurzel von " << a << " ist: " << nth_root(n, a) << endl;
    }
    if (false)
    {
        int n;
        cout << "primes to:";
        cin >> n;
        int* numbers = primes_to_n(n);

        cout << "-Primzahlen bis " << n << "-" << endl;
        for(int i = 0; i < n; i++)
        {
            if(numbers[i]) cout << numbers[i] << endl;
        }
    }

    // serie 2 - pi & karazuba

    if(true)
    {
        cout << "pi:     " << pi(1024) << endl;
        cout << "better: " << better_pi(INT8_MAX) << endl;
        cout << "        " << better_pi(INT16_MAX) << endl;
        cout << "        " << better_pi(INT32_MAX) << endl;
    }

    return 0;
}
#include "1_wurzeln_und_prim/root_and_primes.cpp"

using namespace std;

int main()
{
    if(false)
    {
        int a, n;
        cout << "A:";
        cin >> a;
        cout << endl << "N:";
        cin >> n;
        cout << endl << n << ". Wurzel von " << a << " ist: " << nth_root(n, a) << endl;
    }
    if (true)
    {
        int n;
        cout << "primes to:";
        cin >> n;
        int* numbers = primes_to_n(n);

        cout << "-Primzahlen bis " << n << "-" << endl;
        for(int i = 0; i < n; i++)
        {
            cout << numbers[i] << endl;
        }
    }

    return 0;
}
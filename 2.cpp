#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>

using namespace std;

// void sort(int n, int m, double** a) {
//     double* new_a = new double[n*m];
//     int k = 0;
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < m; j++) {
//             new_a[k++] = a[i][j];
//         }
//     }
//     for (int i = 0; i < n*m-1; i++) {
//         for (int j = 0; j < n*m-i-1; j++) {
//             if (new_a[j] > new_a[j+1]) {
//                 double save = new_a[j];
//                 new_a[j] = new_a[j+1];
//                 new_a[j + 1] = save;
//             }
//         }
//     }
//     k = 0;
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < m; j++) {
//             a[i][j] = new_a[k++];
//         }
//     }
//     delete[] new_a;
// }

double median(int n, int m, double** a) {
    int size = n * m;

    double* new_a = new double[n * m];
    int k = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            new_a[k++] = a[i][j];
        }
    }

    double med = 0;
    if (size % 2 == 0) {
        med = (new_a[size / 2] + new_a[size / 2 - 1]) / 2;
    } else {
        med = new_a[size / 2];
    }
    delete[] new_a;

    return med;
}

// void cumul_weights(int n, int m, double** a) {
//     double tot = 0;
//     double* cw = new double [n*m];

//     sort(n, m, a);
    
//     double* new_a = new double[n*m];
//     int k = 0;
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < m; j++) {
//             new_a[k++] = a[i][j];
//         }
//     }

//     for (int i = 0; i < n*m; i++) {
//         tot += new_a[i];
//         cw[i] = tot;
//     }

//     double W = tot/2;
// }

void sort(int n, int m, double** a) {
    double** transposed = new double*[m];
    for (int i = 0; i < m; i++) {
        transposed[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            transposed[j][i] = a[i][j];
        }
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n - 1; j++) {
            for (int k = 0; k < n - j - 1; k++) {
                if (transposed[i][k] > transposed[i][k + 1]) {
                    swap(transposed[i][k], transposed[i][k + 1]);
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            a[i][j] = transposed[j][i];
        }
    }

    for (int i = 0; i < m; i++) {
        delete[] transposed[i];
    }
    delete[] transposed;
}

double compute_weighted_median(double** a, int size) {
    sort(size, 2, a);

    double total_weight = 0.0;
    for (int i = 0; i < size; i++) {
        total_weight += a[i][1];
    }

    double cumulative_weight = 0.0;

    for (int i = 0; i < size; i++) {
        cumulative_weight += a[i][1];
        if (cumulative_weight >= 0.5) {
            return a[i][0];
        }
    }

    return a[size - 1][0];
}

void w_med(void* inn, void* inoutt, int* len, MPI_Datatype* dt) {
    double* in = static_cast<double*>(inn);
    double* inout = static_cast<double*>(inoutt);

    int total_len = *len * 2;
    double* combined = new double[total_len];

    for (int i = 0; i < *len; i++) {
        combined[i] = in[i];
    }
    for (int i = 0; i < *len; i++) {
        combined[*len + i] = inout[i];
    }

    for (int i = 0; i < total_len - 1; i++) {
        for (int j = 0; j < total_len - i - 1; j++) {
            if (combined[j] > combined[j + 1]) {
                swap(combined[j], combined[j + 1]);
            }
        }
    }

    if (total_len % 2 == 0) {
        *inout = (combined[total_len / 2 - 1] + combined[total_len / 2]) / 2.0;
    } else {
        *inout = combined[total_len / 2];
    }

    delete[] combined;
}

int main() {
    int ProcRank, ProcNum;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    const int n = 5, m = 2;

    double** x = new double*[n];
    for (int i = 0; i < n; i++) {
        x[i] = new double[m];
    }
    srand(time(0) + ProcRank);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 1; j++) {
            x[i][j] = rand() % 10;
        }
    }

    double** w = new double*[n];
    for (int i = 0; i < n; i++) {
        w[i] = new double[m];
    }
    srand(time(0) + ProcRank);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 1; j++) {
            w[i][j] = static_cast<double>(rand()) / (RAND_MAX * 10);
        }
    }

    double** array = new double*[n];
    for (int i = 0; i < n; i++) {
        array[i] = new double[2];
        array[i][0] = x[i][0];
        array[i][1] = w[i][0];
    }

    // cout << "Process " << ProcRank << " data:" << endl;
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < m; j++) {
    //         cout << array[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    double x_med_per_pr = compute_weighted_median(array, n);

    cout << "Process " << ProcRank << ": weighted median = " << x_med_per_pr << endl;

    for (int i = 0; i < n; i++) {
        delete[] array[i];
    }
    delete[] array;

    double median;
    double start_time = MPI_Wtime();

    MPI_Op op;
    MPI_Op_create(&w_med, true, &op);

    MPI_Reduce(&x_med_per_pr, &median, 1, MPI_DOUBLE, op, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    if (ProcRank == 0) {
        cout << "Global median: " << median << endl;
        cout << "Time taken: " << end_time - start_time << " seconds" << endl;
    }

    MPI_Op_free(&op);
    MPI_Finalize();

    return 0;
}

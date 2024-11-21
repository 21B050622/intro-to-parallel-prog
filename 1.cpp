#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std; 

void weighted_median_reduce(void* inn, void* inoutt, int* len, MPI_Datatype* datatype) {
    double* in = static_cast<double*>(inn);
    double* inout = static_cast<double*>(inoutt);

    vector<pair<double, double>> all_data;
    for (int i = 0; i < *len; i += 2) {
        all_data.emplace_back(in[i], in[i + 1]);     // invec
        all_data.emplace_back(inout[i], inout[i + 1]); // inoutvec
    }

    sort(all_data.begin(), all_data.end());

    double total_weight = 0.0;
    for (const auto& pair : all_data) {
        total_weight += pair.second; 
    }

    double half_weight = total_weight / 2.0;
    double cumulative_weight = 0.0;
    double median_value = 0.0;

    for (const auto& pair : all_data) {
        cumulative_weight += pair.second;
        if (cumulative_weight >= half_weight) {
        // if (cumulative_weight >= 0.5) {
            median_value = pair.first; 
            break;
        }
    }

    for (int i = 0; i < *len; i += 2) {
        inout[i] = median_value;    
        inout[i + 1] = total_weight; 
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double local_data[2] = {static_cast<double>(rank), static_cast<double>(rank)/1999}; // (значение, вес)

    double global_result[2]; 

    cout << "Process " << rank << ": Value = " << local_data[0] << ", Weight = " << local_data[1] << endl;

    // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Op weighted_median_op;
    MPI_Op_create(&weighted_median_reduce, 1, &weighted_median_op);

    MPI_Reduce(local_data, global_result, 2, MPI_DOUBLE, weighted_median_op, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "\n=== Final Result ===" << endl;
        cout << "Weighted Median: " << global_result[0] << "\n";
        cout << "Total Weight: " << global_result[1] << "\n";
    }

    MPI_Op_free(&weighted_median_op);

    MPI_Finalize();
    return 0;
}


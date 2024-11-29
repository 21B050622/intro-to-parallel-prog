#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

const double eps = 1e-6;

int main()
{
    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 180, m = 180;
    double dx = 1.0/n, dy = 1.0/m;

    int k = sqrt(size);
    int dims[2] = {k, k};
    int periods[2] = {0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    int up, down, left, right;
    MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
    MPI_Cart_shift(cart_comm, 1, 1, &left, &right);

    // int x_start = n/k * (rank/k);
    // int y_start = m/k * (rank%k);

    int x_start = (n / k) * coords[0];
    int x_next = (m / k) * (coords[0] + 1);
    int y_start = (n / k) * coords[1];
    int y_next = (m / k) * (coords[1] + 1);

    int r = x_next - x_start;
    int c = y_next - y_start;

    vector<vector<double>> u(r+2, vector<double>(c+2));

    for (int i = 1; i <= r; ++i) {
        for (int j = 1; j <= c; ++j) {
            // grid
            double x = (x_start + i - 1) * dx;
            double y = (y_start + j - 1) * dy;

            // b.c
            if (x_start + i - 1 == 0)
                u[i][j] = 0;
            if (x_start + i - 1 == n - 1)
                u[i][j] = 2*y*y + y;
            if (y_start + j - 1 == 0)
                u[i][j] = 0;
            if (y_start + j - 1 == m - 1)
                u[i][j] = 0.5*x*x + 0.5*x;
        }
    }

    double max_diff;
    double start_time = MPI_Wtime(); 

    do {
        max_diff = 0;

        if (up != MPI_PROC_NULL) {
            double r1 = 0;
            MPI_Sendrecv(&u[x_start+1][1], c, MPI_DOUBLE, up, 0, &u[x_start][1], c, MPI_DOUBLE, up, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        if (down != MPI_PROC_NULL) {
            MPI_Sendrecv(&u[x_next-2][1], c, MPI_DOUBLE, down, 0, &u[x_next-1][1], c, MPI_DOUBLE, down, 0, cart_comm, MPI_STATUS_IGNORE);
        }
        
        if (left != MPI_PROC_NULL) {
            vector<double> l_c(r), recv_l_c(r);
            for (int i = 0; i < r; i++)
                l_c[i] = u[i+1][y_start+1];

            MPI_Sendrecv(l_c.data(), r, MPI_DOUBLE, left, 0, recv_l_c.data(), r, MPI_DOUBLE, left, 0, cart_comm, MPI_STATUS_IGNORE);

            for (int i = 0; i < r; i++)
                u[i + 1][y_start] = recv_l_c[i];
        }
        if (right != MPI_PROC_NULL) {
            vector<double> r_c(r), recv_r_c(r);
            for (int i = 0; i < r; i++)
                r_c[i] = u[i+1][y_next-2];

            MPI_Sendrecv(r_c.data(), r, MPI_DOUBLE, right, 0, recv_r_c.data(), r, MPI_DOUBLE, right, 0, cart_comm, MPI_STATUS_IGNORE);

            for (int i = 0; i < r; i++)
                u[i+1][y_next-1] = recv_r_c[i];
        }

        // if (up != MPI_PROC_NULL) {
        //     double r1 = 0;
        //     MPI_Sendrecv(&u[1][1], c, MPI_DOUBLE, up, 0, &u[0][1], c, MPI_DOUBLE, up, 0, cart_comm, MPI_STATUS_IGNORE);
        // }
        // if (down != MPI_PROC_NULL) {
        //     MPI_Sendrecv(&u[r][1], c, MPI_DOUBLE, down, 0, &u[r+1][1], c, MPI_DOUBLE, down, 0, cart_comm, MPI_STATUS_IGNORE);
        // }
        
        // if (left != MPI_PROC_NULL) {
        //     vector<double> l_c(r), recv_l_c(r);
        //     for (int i = 0; i < r; i++)
        //         l_c[i] = u[i+1][1];

        //     MPI_Sendrecv(l_c.data(), r, MPI_DOUBLE, left, 0, recv_l_c.data(), r, MPI_DOUBLE, left, 0, cart_comm, MPI_STATUS_IGNORE);

        //     for (int i = 0; i < r; i++)
        //         u[i+1][0] = recv_l_c[i];
        // }
        // if (right != MPI_PROC_NULL) {
        //     vector<double> r_c(r), recv_r_c(r);
        //     for (int i = 0; i < r; i++)
        //         r_c[i] = u[i+1][c-2];

        //     MPI_Sendrecv(r_c.data(), r, MPI_DOUBLE, right, 0, recv_r_c.data(), r, MPI_DOUBLE, right, 0, cart_comm, MPI_STATUS_IGNORE);

        //     for (int i = 0; i < r; i++)
        //         u[i+1][c-1] = recv_r_c[i];
        // }

        for (int i = 1; i <= r; i++) {
            for (int j = 1; j <= c; j++) {
                if ((x_start + i - 1 == 0 || x_start + i - 1 == n - 1) || (y_start + j - 1 == 0 || y_start + j - 1 == m - 1)) {
                    continue;
                }
                double old_u = u[i][j];
                u[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1])/4;

                double x = (x_start + i - 1) * dx;
                double y = (y_start + j - 1) * dy;

                if (x_start + i - 1 == 0)
                    u[i][j] = 0;
                if (x_start + i - 1 == n - 1)
                    u[i][j] = 2*y*y + y;
                if (y_start + j - 1 == 0)
                    u[i][j] = 0;
                if (y_start + j - 1 == m - 1)
                    u[i][j] = 0.5*x*x + 0.5*x;

                max_diff = max(max_diff, fabs(u[i][j] - old_u));
            }
        }

        double global_max_diff;
        MPI_Allreduce(&max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, cart_comm);
        max_diff = global_max_diff;

    } while (max_diff > eps);

    double end_time = MPI_Wtime(); 
    double time = end_time - start_time;

    vector<double> local_flattened(r*c);
    for (int i = 1; i <= r; ++i) {
        for (int j = 1; j <= c; ++j) {
            local_flattened[(i-1)*c + (j-1)] = u[i][j];
        }
    }

    vector<double> global_flattened;
    if (rank == 0) {
        global_flattened.resize(n*m);
    }

    MPI_Gather(local_flattened.data(), r*c, MPI_DOUBLE, global_flattened.data(), r*c, MPI_DOUBLE, 0, cart_comm);

    if (rank == 0) {
        ofstream outfile("p.txt");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                outfile << global_flattened[i*m + j] << " ";
            }
            outfile << endl;
        }
        outfile.close();
        cout << "time " << time << endl;
    }

    MPI_Finalize();
    return 0;
}


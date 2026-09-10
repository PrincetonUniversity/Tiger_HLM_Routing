#include <mpi.h>

#include "build_info.hpp"
#include "omp_info.hpp"
#include "model_setup.hpp"
#include "routing.hpp"
#include "end_info.hpp"

int main(int argc, char* argv[])
{
    // Threads call MPI only from the region outside the parallel constructs, but the
    // runtime is asked for FUNNELED explicitly rather than left to chance.
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int rank = 0, n_ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    if (rank == 0) printBuildInfo(); // Print build information
    setupOpenMP(); // Every rank: it sets the thread count, not just prints
    ModelSetup setup = setupModel(argv[1]); // Set up the model using the provided configuration file
    applyOmpSchedule(setup.config.omp_schedule, setup.config.omp_chunk); // Schedule for the link loop
    runRouting(setup, rank, n_ranks); // Run the routing process
    if (rank == 0) printEndInfo(); // Print end information
    MPI_Finalize();
    return 0;
}
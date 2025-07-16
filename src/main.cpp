#include "build_info.hpp"
#include "omp_info.hpp"
#include "model_setup.hpp"
#include "routing.hpp"
#include "end_info.hpp"

int main(int argc, char* argv[])
{   
    printBuildInfo(); // Print build information
    setupOpenMP(); // Set up OpenMP
    ModelSetup setup = setupModel(argv[1]); // Set up the model using the provided configuration file
    runRouting(setup); // Run the routing process
    printEndInfo(); // Print end information
    return 0;
}
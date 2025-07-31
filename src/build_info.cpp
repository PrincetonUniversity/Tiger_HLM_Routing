// build_info.cpp
#include <iostream>
#include "build_info.hpp"

/**
 * @brief Prints the build information of the Tiger HLM routing application.
 */
void printBuildInfo() {
    std::cout << "___________________BUILD INFO_____________________ \n" << std::endl;
    std::cout << "Tiger HLM routing version 1.0.0" << std::endl;
    std::cout << "Tiger HLM routing build date: " << __DATE__ << " " << __TIME__ << std::endl;
    std::cout << "Tiger HLM routing build by: am2192" << std::endl;
    std::cout << "__________________________________________________ \n" << std::endl;
}

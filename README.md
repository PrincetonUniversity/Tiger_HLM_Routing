# Tiger-HLM Routing

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


Tiger-HLM Routing software for converting hillslope generative runoff to streamflow with nonlinear routing equation. The complementary GPU based runoff generation software can be [here](https://github.com/PrincetonUniversity/Tiger_HLM_GPU). This software can be used as a stand alone as well. 


Please read [our wiki page](https://github.com/PrincetonUniversity/Tiger_HLM_Routing/wiki) to get started.


## Directory Structure

The organizational structure of the Tiger HLM Routing is as follows:

```text
├── bin                         # Compiled binaries directory
│   └── routing                 # Final executable for the routing model
│
├── build                       # Intermediate build files (object files)
│
├── compile.sh                  # Shell script to compile the project manually
│
├── data                        # Input data folder (contains example)
│   ├── hydrography             # Input hydrography files for example
│   │   ├── hillslopes.gpkg  
│   │   ├── links.gpkg
│   │   └── params_table.csv
│   ├── runoff                  # Input runoff files for example
│   │   ├── total_runoff_2017.nc  
│   │   └── total_runoff_2018.nc 
│   ├── config.yaml             # Model configuration file
│   ├── example.ipynb           # Notebook containing information for example
│   └── params.csv  # Network structure and routing parameters CSV
│
├── Makefile                    # Build automation instructions using `make`
│
└── src                         # Source code directory
    ├── build_info.cpp          # Stores Git/build metadata for inclusion in binary
    ├── build_info.hpp          # Header for build_info.cpp
    ├── end_info.cpp            # Handles logging and final reporting
    ├── end_info.hpp            # Header for end_info.cpp
    ├── I_O                     # Input/Output module source files
    │   ├── config_loader.cpp   # Loads YAML config file
    │   ├── config_loader.hpp
    │   ├── inputs.cpp          # Reads NetCDF, parameter CSVs, etc.
    │   ├── inputs.hpp
    │   ├── node_info.cpp       # Processes node/stream metadata
    │   ├── node_info.hpp
    │   ├── output_series.cpp   # Handles time series NetCDF output
    │   └── output_series.hpp
    ├── main.cpp                # Entry point for the routing model
    ├── models
    │   └── RHS.hpp             # Right-hand side of ODE system for routing equations
    ├── model_setup.cpp         # Initializes model state, time, and parameters
    ├── model_setup.hpp
    ├── omp_info.cpp            # Detects and logs OpenMP thread info
    ├── omp_info.hpp
    ├── routing.cpp             # Core routing equation logic
    ├── routing.hpp
    └── utils
        ├── time.cpp            # Time utilities (e.g., timestep parsing, conversions)
        └── time.hpp
```


## Citation
To cite this software in your publication, please use the following BibTeX (to be updated upon paper acceptance) to refer to the code's [method paper](empty):
```
@article{,
	doi = {},
	url = {},
	year = ,
	month = ,
	publisher = {},
	volume = {},
	number = {},
	pages = {}
	author = {},
	title = {},
	journal = {},
}
```

Finally, we will have DOIs for each released version on Zenodo. This approach promotes computational reproducibility by allowing you to specify the exact version of the code used to generate the results presented in your publication. A working zenodo badge will be added above once the first version is released. 

```
@software{,
  author       = {Tiger HLM development team},
  title        = {},
  month        = ,
  year         = ,
  publisher    = {},
  version      = {},
  doi          = {},
  url          = {}
}
```
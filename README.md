# EPS Constraint for the Chesapeake Bay Watershed

## Overview
This project implements a hybrid optimization strategy for the Chesapeake Bay watershed. It combines a point-based interior-point optimization (IPOPT) method and an $\epsilon$-constraint method to enhance the efficiency of multi-objective evolutionary algorithms. This approach systematically explores different constraints to generate a collection of viable solutions for environmental management.

## Description
The $\epsilon$-constraint method modifies the original single-objective problem into a series of single-objective problems by adjusting the constraint on the second objective function, $f_2$. This is done to produce a series of solutions that contribute to forming a Pareto front. The method is defined by the following optimization formulation:
Minimize $f_1(\mathbf{x})$, 
Subject to  $f_2(\mathbf{x}) \leq \epsilon_k f_2^{\text{base}},$
     $\mathbf{x} \in \mathbf{X},$
$$

where $f_2^{\text{base}}$ is the baseline nitrogen loading and $\mathbf{X}$ represents the feasible set of decision variables. The method involves multiple executions for different values of $\epsilon_k$, ranging from 1.0 to 0.70 in decreasing steps of 0.03, to achieve a comprehensive set of solutions.

## Prerequisites
Before using this project, ensure that [Ipopt](https://coin-or.github.io/Ipopt/index.html) is installed, as it is crucial for the initial seeding of the optimization process.


Follow instructions: https://coin-or.github.io/Ipopt/index.html

### Debian Linux:
```
apt install gcc g++ gfortran git patch wget pkg-config liblapack-dev libmetis-dev gawk
```


### Mac OS:
```
brew update
brew install bash gcc
brew link --overwrite gcc
brew install pkg-config
brew install metis
```

### ASL
Execute the following commands:

```
git clone https://github.com/coin-or-tools/ThirdParty-ASL.git
cd ThirdParty-ASL
./get.ASL
./configure
make
sudo make install

cd ../
```
### HSL
Obtain the HSL code in http://hsl.rl.ac.uk/ipopt.
```
git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
cd ThirdParty-HSL
copy HSL package downloaded from http://hsl.rl.ac.uk/ipopt and rename it as coinhsl
./configure
make
sudo make install
```
### MacOSX does not link the libhsl.
```
cd /usr/local/lib
ln -s libcoinhsl.2.dylib libhsl.dylib
```

### Mumps

https://github.com/coin-or-tools/ThirdParty-Mumps
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git

cd ThidParty-Mumps/
./configure ADD_FCFLAGS=-fallow-argument-mismatch
make
make install

##### Compile Ipopt
```
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
./configure CXX=g++--9 CC=gcc-9 F77=gfortran-9
make
sudo make install

# Compile
make
```


## Installation

Follow the standard building process as described in previous examples:

1. Clone the repository:
   ```bash
   git clone https://github.com/gtoscano/eps_cnstr.git
   ```
2. Navigate to the project directory:
   ```bash
   cd eps_cnstr
   ```
3. Follow the build instructions specific to your system, generally similar to:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ```

## Usage
To run the optimization, execute the built application from the command line, adjusting parameters as necessary for your specific analysis needs.

## Contributing
Contributions to enhance or expand the methodology are welcome. Please fork the repository and submit a pull request with your proposed changes.

## License
This project is available under the MIT License. See the LICENSE file for more details.


## Acknowledgments
This work relies on the foundational algorithms provided by the Ipopt project, a component of the COIN-OR initiative.

```


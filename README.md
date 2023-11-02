# Install Ipopt

Follow instructions: https://coin-or.github.io/Ipopt/index.html


## Linux:
```
wajig install gcc g++ gfortran git patch wget pkg-config liblapack-dev libmetis-dev gawk
```


## Mac OS:
```
brew update
brew install bash gcc
brew link --overwrite gcc
brew install pkg-config
brew install metis
```

## Instalation process
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
#### MacOSX does not link the libhsl.
```
cd /usr/local/lib
ln -s libcoinhsl.2.dylib libhsl.dylib
```

## Mumps

https://github.com/coin-or-tools/ThirdParty-Mumps
git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git

cd ThidParty-Mumps/
./configure ADD_FCFLAGS=-fallow-argument-mismatch
make
make install

# Compile Ipopt
```
git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
./configure CXX=g++--9 CC=gcc-9 F77=gfortran-9
make
sudo make install

# Compile
make
```

# Run
./epa prefix_files load_pct
```
./epa Bradford_PA .9
```


# Parallel machine scheduling with decision diagrams

In this repository, you will find the implementation used to produce the results
in the following paper:

We present a new flow-based formulation for identical parallel machine
scheduling with a regular objective function and no idle time. We use a decision
diagram containing all the possible sequences of jobs that follow some ordering
rules to construct the new formulation. These rules do not exclude the optimal
solution of a given instance but constrain the optimal solution to some
canonical form. To define these rules, we need to partition the planning horizon
into non-uniform periods. The novel formulation will have numerous variables and
constraints. Hence, we apply a Dantzig-Wolfe decomposition to compute the linear
programming relaxation of the new flow-based formulation in a reasonable amount
of time.

Moreover, we will see that the resulting lower bound will be stronger than the
lower bound provided by the linear relaxation of the classical time-indexed
formulation. We use a Branch-and-Price framework to solve the new formulation.
Several instances from the literature will be solved for the first time.​

## Instructions

In order to reproduce the results, you need to perform the following steps after
downloading the package.  These instructions will get you a copy of the project
up and running on your local machine for development and testing purposes.

### Prerequisites
* **CMake v3.21+** - found at [https://cmake.org/](https://cmake.org/)

* **C++ Compiler** - needs to support at least the **C++20** standard, i.e.
  *MSVC* at least version 142, *GCC* at least version 10 for linux based systems

* **Vcpkg** - C/C++ dependency manager from Microsoft. For more information on
  how to install and to use vcpkg, we refer to [https://vcpkg.io](https://vcpkg.io).
  Do not forget to integrate vcpkg with your environment, i.e. add the vcpkg binary
  directory to your environment variable path. For example, on linux based systems
  you can add the following lines to your .bashrc file or just execute them in the terminal:

```bash
export VCPKG_ROOT=<path to vcpkg>
```
  CMakelists.txt will automatically detect vcpkg and will install the required dependencies, but you need export the environment variable VCPKG_ROOT.

* **Gurobi** - We use the Gurobi optimizer to compute the linear programming
  (LP) relaxations. You can download Gurobi from the company's
  [website](https://www.gurobi.com/). Follow the installation instructions that
  are provide by gurobi. Per Operating system you will need to adjust
  environment variables such that cmake can find gurobi on your computer. On
  windows for example the environment variables are automatically set.

* **Osi** - Also [Osi](https://github.com/coin-or/Osi) is used. Osi (Open Solver
  Interface) provides an abstract class to generic linear LP, together with
  derived classes for specific solvers. In theory, we can use an arbitrary LP
  solver (Gurobi, CPLEX, XPress, Soplex, ...) to solve the LP relaxation, but
  this is not implemented yet. For now, we can only use gurobi. To install Osi,
  we use coinbrew which can be found
  [here](https://coin-or.github.io/coinbrew/). Please follow the instructions
  [on](https://coin-or.github.io/coinbrew/) to install on your system. Please
  install the Osi in a directory called ThirdParty, i.e. apply coinbrew in the
  directory ThirdParty:

#### Linux

```bash
mkdir -p ThirdParty && cd ThirdParty
wget -O coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Cgl@master
./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-release" --with-gurobi-lflags="-L$GUROBI_HOME/lib -l${GUROBI_VERSION} -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none
./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-debug" --reconfigure --enable-debug --with-gurobi-lflags="-L$GUROBI_HOME/lib -l${GUROBI_VERSION} -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none
``` 

Do not forget to adjust the environment variables such that cmake can find the Gurobi and Osi libraries. For example, on linux based systems you can add the following lines to your .bashrc file or just execute them in the terminal:

```bash
export GUROBI_HOME="/opt/gurobi912/linux64"
export GUROBI_VERSION="gurobi91"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
export COIN_OR_X64_LINUX_RELEASE_DIR="/home/daniel/ThirdParty/coin-or-x64-linux-release"
export COIN_OR_X64_LINUX_DEBUG_DIR="/home/daniel/ThirdParty/coin-or-x64-linux-debug"
```

#### Windows (Not supported, but you can try to compile it yourself)
To compile Osi on Windows based systems you should install
[Msys2](https://www.msys2.org/) first. Follow the next instructions instructions
to setup the developers environment in order to compile Osi with visual studio.
Use powershell to setup your environment correctly. First add the binaries of
msys2 to the environment variable path:

```powershell
$env:Path = "C:\msys64\usr\bin\;$env:Path"
```

Next define the following function in Powershell such that the compiler
executables of visual studio are recognized by your powershell environment.

```powershell
function Invoke-BatchFile
{
   param([string]$Path, [string]$Parameters)

   $tempFile = [IO.Path]::GetTempFileName()

   ## Store the output of cmd.exe.  We also ask cmd.exe to output
   ## the environment table after the batch file completes
   cmd.exe /c " `"$Path`" $Parameters && set > `"$tempFile`" "

   ## Go through the environment variables in the temp file.
   ## For each of them, set the variable in our local environment.
   Get-Content $tempFile | Foreach-Object {
       if ($_ -match "^(.*?)=(.*)$")
       {
           Set-Content "env:\$($matches[1])" $matches[2]
       }
   }

   Remove-Item $tempFile
}
```

You can now execute the following command in powershell. Here we use visual
studio code 2022, but after correct adjustments you will probably be also
to use earlier version of visual studio as well.

```powershell
Invoke-BatchFile 'C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat' x64
```

Start now msys2:

```powershell
bash

```
Install now some tools that are needed to compile Osi with visual studio:

```bash
pacman -S make wget tar patch dos2unix diffutils git svn pkg-config zip unzip
pacman -S mingw-w64-i686-toolchain mingw-w64-x86_64-toolchain
pacman -S mingw-w64-x86_64-lapack \
          mingw-w64-x86_64-winpthreads-git \
          mingw-w64-x86_64-readline \
          mingw-w64-x86_64-suitesparse \
          mingw-w64-x86_64-metis

```

Execute the following commands in bash with this git repository as the current working directory:

```bash
mkdir -p ThirdParty && cd ThirdParty
wget -O coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew build Cgl --prefix="coin-or-x64-MD" --tests=none --enable-msvc --build=x86_64-w64-mingw32 --enable-shared=MD --with-gurobi-lflags="-L/c/gurobi912/win64/lib -lgurobi91" --with-gurobi-cflags="-I/c/gurobi912/win64/include"
./coinbrew build Cgl --prefix="coin-or-x64-MDd" --tests=none --enable-debug --enable-msvc --build=x86_64-w64-mingw32 --enable-shared=MDd --with-gurobi-lflags="-L/c/gurobi912/win64/lib -lgurobi91" --with-gurobi-cflags="-I/c/gurobi912/win64/include"

```


### Building the project

To build the project:

#### Linux

```bash
cmake --preset gcc-11-release && cmake --build --preset build-gcc-11-release
cmake --preset gcc-11-debug && cmake --build --preset build-gcc-11-debug

```

#### Windows

```powershell
cmake --preset windows64-msvc-2022 && cmake --build --preset build-windows64-debug
cmake --preset windows64-msvc-2022-release && cmake --build --preset build-windows64-release

```

#### Docker
We also provide a docker file to build the project. To run the docker [container](https://hub.docker.com/repository/docker/danielkowalczyk/pm_linux_compute_env/general)

```bash
docker run -it --rm -v .:/workspaces/implementation -w /workspaces/implementation danielkowalczyk/pm_linux_compute_env:latest bash
cmake --preset gcc-11-release && cmake --build --preset build-gcc-11-release
```

Or you can also use the devcontainer provided in the repository. For more information on how to use devcontainers in VS code, please refer to the [documentation](https://code.visualstudio.com/docs/remote/containers).

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our how you can
become a contributor and the process for submitting pull requests to us.

## Authors

* **Daniel Kowalczyk** - [@danielkowalczyk](https://gitlab.kuleuven.be/u0056096)

## License

This project is licensed under the [Unlicense](https://unlicense.org/) - see the
[LICENSE](LICENSE) file for details

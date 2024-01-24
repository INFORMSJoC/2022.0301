#!/bin/bash

export CMAKE_VERSION="3.21.4"
export GRB_VERSION="9.1.2_linux64"
export GRB_SHORT_VERSION="9.1"
export GRB_VERSION_DIR="912"

mkdir -p ThirdParty
cd ThirdParty

pip install -U cmake

if [ ! -d "./gurobi" ]; then
	wget https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}.tar.gz
    tar -xvf gurobi${GRB_VERSION}.tar.gz
    rm gurobi${GRB_VERSION}.tar.gz
    mv gurobi${GRB_VERSION_DIR} gurobi
	rm -rf gurobi/linux64/docs
    export GUROBI_HOME=$(pwd)/gurobi/linux64
    ln -sf $GUROBI_HOME/lib/libgurobi_g++5.2.a $GUROBI_HOME/lib/libgurobi_c++.a
    export LD_LIBRARY_PATH=$GUROBI_HOME/lib
    ldconfig -v
else
	export GUROBI_HOME=$(pwd)/gurobi/linux64
fi

if [ ! -d "./coin-or-x64-linux-debug" ] || [ ! -d "coin-or-x64-linux-release" ]; then
	wget -O coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
	chmod u+x coinbrew
	./coinbrew fetch Cgl@master
	if [ ! -d "./coin-or-x64-linux-release" ]; then
		./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-release" --with-gurobi-lflags="-L$GUROBI_HOME/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none
	fi

	if [ ! -d "./coin-or-x64-linux-debug" ]; then
		./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-debug" --reconfigure --enable-debug --with-gurobi-lflags="-L$GUROBI_HOME/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none
	fi
fi

if [ ! -d "./vcpkg" ]; then
	git clone https://github.com/Microsoft/vcpkg.git vcpkg
	./vcpkg/bootstrap-vcpkg.sh
	./vcpkg/vcpkg integrate install
fi

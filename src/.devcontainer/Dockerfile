# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.209.3/containers/cpp/.devcontainer/base.Dockerfile
ARG VARIANT="jammy"
FROM mcr.microsoft.com/vscode/devcontainers/cpp:${VARIANT}

# [Optional] Uncomment this section to install additional packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends wget ca-certificates git curl unzip tar zip pkg-config gfortran ninja-build gpg-agent software-properties-common

# Add current cmake/ccmake, from Kitware
ARG CMAKE_URL="https://apt.kitware.com/ubuntu/"
ARG CMAKE_PKG="jammy"
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null \
    | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null && \
    apt-add-repository -y "deb ${CMAKE_URL} ${CMAKE_PKG} main" && \
    apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends cmake

# packages needed by coin-or
RUN apt-get -y install --no-install-recommends libgfortran-10-dev libblas-dev liblapack-dev libmetis-dev libnauty2-dev

ARG GRB_VERSION=9.1.2
ARG GRB_SHORT_VERSION=9.1

# install gurobi package and copy the files
WORKDIR /opt
RUN update-ca-certificates \
    && wget -v https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_linux64.tar.gz \
    && tar -xvf gurobi${GRB_VERSION}_linux64.tar.gz  \
    && rm -f gurobi${GRB_VERSION}_linux64.tar.gz \
    && mv -f gurobi* gurobi \
    && rm -rf gurobi/linux64/docs \
    && ln -sf /opt/gurobi/linux64/lib/libgurobi_g++5.2.a /opt/gurobi/linux64/lib/libgurobi_c++.a

ENV GUROBI_HOME /opt/gurobi/linux64
ENV PATH $PATH:$GUROBI_HOME/bin
ENV LD_LIBRARY_PATH $GUROBI_HOME/lib

# install coin-or packages into ThirdParty folder
WORKDIR /tmp
RUN mkdir coin-or

WORKDIR /tmp/coin-or
RUN wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
RUN chmod u+x coinbrew && \
    ./coinbrew build Cgl@master --with-gurobi-lflags="-L${GUROBI_HOME}/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I ${GUROBI_HOME}/include" --tests none --prefix=coin-or-x64-linux-release && \
    ./coinbrew build Cbc@master --enable-debug --with-gurobi-lflags="-L${GUROBI_HOME}/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I ${GUROBI_HOME}/include" --tests none --prefix=coin-or-x64-linux-debug && \
    cp -r /tmp/coin-or/ /opt/coin-or && \
    rm -rf /tmp/coin-or && \
    echo "/opt/coin-or/coin-or-x64-linux-release/lib" >> /etc/ld.so.conf.d/coin-or.conf && \
    ldconfig

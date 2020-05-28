FROM liambindle/penelope:0.1.0-ubuntu16.04-gcc7-netcdf4.5.0-netcdff4.4.4

# Make a directory to install GEOS-Chem to
RUN mkdir /opt/geos-chem && mkdir /opt/geos-chem/bin

# Copy the GEOS-Chem repository (".") to /gc-src
# This means this docker build command's context must be 
# GEOS-Chem's root source code directory
COPY . /gc-src
RUN cd /gc-src \
&&  mkdir build

# Commands to properly set up the environment inside the container
RUN echo "module load gcc/7" >> /init.rc \
&&  echo "spack load hdf5" >> /init.rc \
&&  echo "spack load netcdf" >> /init.rc \
&&  echo "spack load netcdf-fortran" >> /init.rc \
&&  echo "export PATH=$PATH:/opt/geos-chem/bin" >> /init.rc \
&&  echo "source /init.rc >> /etc/bash.bashrc"

# Make bash the default shell
SHELL ["/bin/bash", "-c"]

# Build Standard and copy the executable to /opt/geos-chem/bin
RUN cd /gc-src/build \
&&  cmake -DRUNDIR=IGNORE -DRUNDIR_SIM=standard .. \
&&  make -j install \
&&  cp geos /opt/geos-chem/bin/geos-chem-standard \
&& rm -rf /gc-src/build/*

# Build Tropchem and copy the executable to /opt/geos-chem/bin
RUN cd /gc-src/build \
&&  cmake -DRUNDIR=IGNORE -DRUNDIR_SIM=tropchem .. \
&&  make -j install \
&&  cp geos /opt/geos-chem/bin/geos-chem-tropchem \
&& rm -rf /gc-src/build/*

# Build SOA_SVPOA and copy the executable to /opt/geos-chem/bin
RUN cd /gc-src/build \
&&  cmake -DRUNDIR=IGNORE -DRUNDIR_SIM=complexSOA_SVPOA .. \
&&  make -j install \
&&  cp geos /opt/geos-chem/bin/geos-chem-soa_svpoa\
&& rm -rf /gc-src/build/*

RUN rm -rf /gc-src

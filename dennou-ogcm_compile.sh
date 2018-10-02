export FC=mpiifort
export FCFLAGS="-fPIC -cpp -DSJPACK -O3 -openmp -xHost -qopt-report -I${MKLROOT}/include/intel64/lp64"

#
HOST=
NETCDF_DIR=/ap/netcdf4/4.3.3.1-intel
NETCDFF_DIR=/ap/netcdf4-fortran/4.4.2-intel
INSTALL_DIR=/home/ykawai/lib/Dennou-OGCM/
GTOOL5LIB=/home/ykawai/lib/gtool5-serial/lib/libgtool5-serial.a
ISPACKLIB=/home/ykawai/lib/ispack/lib/libisp-avx.a
#ISPACKLIB=/home/ykawai/lib/ispack/lib/libisp-sse64.a
SPMLLIB=/home/ykawai/lib/spml/lib/libspml-omp.a
LAPACKLIB=no
export LDFLAGS="-L${MKLROOT}/lib/intel64"
export SYSLDLIBS=" -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm"


./configure \
    --with-netcdf=${NETCDF_DIR}/lib/libnetcdf.a \
    --with-netcdf-include=${NETCDF_DIR}/include/netcdf.h \
    --with-netcdff=${NETCDFF_DIR}/lib/libnetcdff.a \
    --with-gtool5=${GTOOL5LIB} \
    --with-spml=${SPMLLIB} \
    --with-ispack=${ISPACKLIB} \
    --with-dsogcm-mode=AXISYM \
    --enable-sjpack \
    --with-lapack=${LAPACKLIB} \
    --prefix=${INSTALL_DIR}

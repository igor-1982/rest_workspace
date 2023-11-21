Installation
------------------

* Prerequisites

   - libopenblas.so
     =========================================================
     -> git clone https://github.com/xianyi/OpenBLAS.git OpenBLAS
     -> cd OpenBLAS
     -> make
     -> cp libopenblas.so $REST_EXT_DIR/libpoenblas.so
     =========================================================
   - libcint.so
     =========================================================
     -> git clone https://github.com/sunqm/libcint.git libcint
     -> cd libcint
     -> mkdir build
     -> cmake -DWITH_RANGE_COULOMB=1 ..
     -> make
     -> cp libcint.so $REST_EXT_DIR/libcint.so
     =========================================================
   - libxc.so
     =========================================================
     -> based on your os, download the installation file from 
            https://www.tddft.org/programs/libxc/download/
     =========================================================
   - libhdf5.so
     =========================================================
     -> download the source code from 
          https://www.hdfgroup.org/downloads/hdf5
     -> tar -zcvf hdf5-*.tar.gz
     -> cd hdf
     -> source ./HDF5-*-Linux.sh
     -> cp HDF5-*/HDF_Group/HDF5/*/lib/libhdf5.so $REST_EXT_DIR/
     =========================================================
   - librest2fch.so
     =========================================================
     -> git clone https://gitlab.com/jeanwsr/MOKIT -b for-rest
     -> cd MOKIT/src
     -> make rest2fch -f Makefile.gnu_openblas
     -> cp MOKIT/mokit/lib/librest2fch.so $REST_EXT_DIR/
     =========================================================


* Build REST::

   1) cp Config.templet Config
   2) edit "Config" to make the prerequisite libraries aforementioned accessable 
      via some global variables heading with "REST". Please refer to the Config.templet 
      file for more details.
   3) bash Config
   4) Source $HOME/.bash_profile
   5-1) cargo build           (debug version: compile fast and run slowly)
   5-2) cargo build --release (release version: compile slowly and run fast)


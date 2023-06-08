Installation
------------------

* Prerequisites

   - libopenblas.so
   - libcint.so
   - libhdf5.so
   - libxc.so
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


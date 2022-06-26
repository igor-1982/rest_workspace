Installation
------------------

* Prerequisites

   - libopenblas.so
   - libcint.so
   - libhdf5.so
   - gcc version 8.0 or higher

* Build REST::

   1) cp Config.templet Config
   2) edit "Config" to make the prerequisite libraries aforementioned accessable
   3) bash Config
   4) Source $HOME/.bash_profile
   5-1) cargo build           (debug version: compile fast and run slowly)
   5-2) cargo build --release (release version: compile slowly and run fast)


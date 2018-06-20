# No module environment on the Mac
[[ $(uname) == "Darwin" ]] && return
# Initialise module environment if it is not
if [[ ! $(command -v module > /dev/null 2>&1) ]]; then
  . /usr/local/apps/module/init/bash
fi

module unload eccodes
module unload emos
module unload fftw
module unload libemos

module switch gnu clang
export FC=/usr/local/apps/gcc/6.3.0/bin/gfortran


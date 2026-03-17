#! /usr/bin/env bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if ! command -v fypp >/dev/null 2>&1
then
    echo "fypp could not be found"
    exit 1
fi

include_dir="${SCRIPT_DIR}/.."
for f in $(ls -1 ${SCRIPT_DIR}/../*.fypp); do
  fypp_file=$f
  fypp_file_basename=$(basename $f)
  F90_file="${fypp_file_basename%.*}.F90"
  if [[ ${fypp_file_basename} == "pluto_defs.fypp" ]]; then
    continue
  fi
  if [[ ${fypp_file_basename} == "pluto_module_memoryspace.fypp" ]]; then
    F90_file="pluto_module_host.F90"
    echo "+ fypp -Dmemoryspace=\"host\" ${fypp_file} ${F90_file}"
    fypp -Dmemoryspace=\"host\" -I${include_dir} ${fypp_file} ${F90_file}
    F90_file="pluto_module_device.F90"
    echo "+ fypp -Dmemoryspace=\"device\" ${fypp_file} ${F90_file}"
    fypp -Dmemoryspace=\"device\" -I${include_dir} ${fypp_file} ${F90_file}
  else
    echo "+ fypp ${fypp_file} ${F90_file}"
    fypp -I${include_dir} ${fypp_file} ${F90_file}
  fi
done

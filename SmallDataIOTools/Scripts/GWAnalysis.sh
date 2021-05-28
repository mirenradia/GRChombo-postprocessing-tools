#!/usr/bin/env bash

# Find a Chombo .ex binary within a provided directory and starting with a name
get_exec_path() {
  if [ "$#" != 2 ]; then
    echo "usage $0 <exec path> <exec prefix>";
    exit 1
  fi
  EXEC_PATH=$1;
  EXEC_PREFIX=$2;
  NUM_EXECS=$(find "$EXEC_PATH" -name "${EXEC_PREFIX}3d*.ex"| wc -l);
  if [ "$NUM_EXECS" != 1 ]; then
    echo "There are $NUM_EXECS executables in the directory with the prefix $EXEC_PREFIX";
    exit 1
  fi
  readlink -f "$(find "$EXEC_PATH" -name "${EXEC_PREFIX}3d*.ex")";
}

# Get directory with GRChombo postprocessing tools repository
GRCHOMBO_POSTPROCESSING_DIR=${GRCHOMBO_POSTPROCESSING_DIR-$HOME/GRChombo-postprocessing-tools};

# Set number of OpenMP threads
#OMP_NUM_THREADS=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}');
#export OMP_NUM_THREADS

# Find GWEnergyMomentum executable
GW_EXEC=$(get_exec_path "$GRCHOMBO_POSTPROCESSING_DIR/SmallDataIOTools/GWEnergyMomentum" GWEnergyMomentum);

# Set data subpath
DATA_SUBPATH=${PWD}/data

# Execute GWEnergyMomentum executable in current directory with no arguments
echo $GW_EXEC
$GW_EXEC -p ${DATA_SUBPATH}/Weyl4_extraction_

# Find IntegrateTime executable
INTEGRATETIME_EXEC=$(get_exec_path "$GRCHOMBO_POSTPROCESSING_DIR/SmallDataIOTools/IntegrateTime" IntegrateTime);

# Integrate power and momentum fluxes in time to get cumulative radiated quantities
$INTEGRATETIME_EXEC -o ${DATA_SUBPATH}/GW_energy.dat ${DATA_SUBPATH}/GW_power.dat
$INTEGRATETIME_EXEC -o ${DATA_SUBPATH}/GW_momentum.dat ${DATA_SUBPATH}/GW_momentum_flux.dat
$INTEGRATETIME_EXEC -o ${DATA_SUBPATH}/GW_angular_momentum.dat ${DATA_SUBPATH}/GW_angular_momentum_flux.dat

# Add radii comment line to time integrated files
RADII_LINE_1D=$(sed -n '2p' ${DATA_SUBPATH}/GW_power.dat)
RADII_LINE_3D=$(sed -n '2p' ${DATA_SUBPATH}/GW_momentum_flux.dat)
sed -i "2i${RADII_LINE_1D}" ${DATA_SUBPATH}/GW_energy.dat
sed -i "2i${RADII_LINE_3D}" ${DATA_SUBPATH}/GW_momentum.dat
sed -i "2i${RADII_LINE_3D}" ${DATA_SUBPATH}/GW_angular_momentum.dat

# Find Norm executable
NORM_EXEC=$(get_exec_path "$GRCHOMBO_POSTPROCESSING_DIR/SmallDataIOTools/Norm" Norm)

#Take norms of momentum and angular momentum
$NORM_EXEC -o ${DATA_SUBPATH}/GW_momentum_norm.dat ${DATA_SUBPATH}/GW_momentum.dat
$NORM_EXEC -o ${DATA_SUBPATH}/GW_angular_momentum_norm.dat ${DATA_SUBPATH}/GW_angular_momentum.dat

# Add radii comment lines to normed files
sed -i "2i${RADII_LINE_1D}" ${DATA_SUBPATH}/GW_momentum_norm.dat
sed -i "2i${RADII_LINE_1D}" ${DATA_SUBPATH}/GW_angular_momentum_norm.dat

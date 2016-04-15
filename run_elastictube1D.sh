#!/bin/sh
# set -e
BASE=$PWD
cd $BASE
# ---------------------------------------- PARAMETERS --------------------------------------------------------

# 1d tube parameters
N=100

# run paramters
FLUIDCORES=3
SOLIDCORES=3

# coupling parameters
PPNAME=v-imvj_RS-SVD
CP=serial-implicit
PP=IQN-IMVJ
EXTRAPOLATION=2

# QN parameters
REUSED=0
FILTER=QR2
FILTER_EPS=1e-2
PRECOND=residual-sum

RS=RS-SVD
R_RS=8
TRUNC_EPS=1e-4
CHUNKSIZE=8


COPY=1
ONLY_POSTPROC=0
POSTPROC=0

#NOW="$(date +'%Y-%m-%d')"

# ------------------------------------------------------------------------------------------------------------
D1=results/${PPNAME}/                #${NOW}_FSI-${N}-${NCOARSE}
# ------------------------------------------------------------------------------------------------------------

#DIRNAME=v-imvj_RS-SLIDE_M-${CHUNKSIZE}_reuse-${REUSE}_${FILTER}-${FILTER_EPS}
#DESCR=v-imvj_RS-SLIDE_M-${CHUNKSIZE}_reuse-${REUSE}_${FILTER}-${FILTER_EPS}

cp ConfigurationFiles/precice-config-parallel.xml preCICE.xml
FILE=preCICE.xml


sed -i s/timesteps-reused\ value=\"[0-9]*\"/timesteps-reused\ value=\"${REUSED}\"/g ${FILE}                # set reuse
sed -i s/extrapolation-order\ value=\"[0-9]*\"/extrapolation-order\ value=\"${EXTRAPOLATION}\"/g ${FILE}   # set extrpol. order
sed -i s/post-processing:[A-Za-z-]*/post-processing:${PP}/g ${FILE}                                        # set QN method
sed -i s/coupling-scheme:[A-Za-z-]*/coupling-scheme:${CP}/g ${FILE}                                        # set coupling scheme
sed -i -e s/type=\"[A-Z0-9]*\"\ limit=\"[0-9e-]*\"/type=\"${FILTER}\"\ limit=\"${FILTER_EPS}\"/g ${FILE}   # set filter method


echo "Start Simulation run"

# COMPUTATION, PARAMETER STUDY
if [ ${ONLY_POSTPROC} = 0 ]; then
    #sed -i s/extrapolation-order\ value=\"[0-9]*\"/extrapolation-order\ value=\"${EXTRAPOLATION}\"/g ${FILE}

for CHUNKSIZE in 4 # 4 8 16 32
do 
 for TRUNC_EPS in 1e-4 # 1e-4 1e-3 1e-2 1e-1 5e-1
 do
    for KAPPA in 10 # 10  100   1000
    do
        for TAU in 0.001 # 0.1 0.01 0.001
        do
            echo "\n ############################### \n"
            echo " run 1d elastictube with N="${N}", tau="${TAU}", kappa="${KAPPA}
            echo " coupling-scheme: "${CP}", post-processing: "${PP}
            echo " reuse: "${REUSED}
            if [ ${RS} = "RS-SVD" ]; then
              echo "RS-SVD: chunksize: "${CHUNKSIZE}", trunc: "${TRUNC_EPS}
            fi
            echo " filter: "${FILTER}", eps: "${FILTER_EPS}
            echo "\n ###############################"
 
            # clean up
            ./clean

            if [ ${FLUIDCORES} = 1 ]; then
              ./FluidSolver ${FILE} $N ${TAU} ${KAPPA}  > log.fluid 2>&1 &
            else
              mpirun -np ${FLUIDCORES} ./FluidSolver ${FILE} $N ${TAU} ${KAPPA}  > log.fluid 2>&1 &
            fi
            if [ ${SOLIDCORES} = 0 ]; then
              ./StructureSolver ${FILE} $N  > log.structure 2>&1
            else
              mpirun -np ${SOLIDCORES} ./StructureSolver ${FILE} $N  > log.structure 2>&1
            fi

            DEST_DIR=results/${PPNAME}/\[${M}_${TRUNC_EPS}_${REUSED}\]  
            DESCR=${PPNAME}\[${M}_${TRUNC_EPS}_${REUSED}\]
            
            # make result dir
            if [ ${COPY} = 1 ]; then
              if [ ! -d ${D1} ]; then
                 mkdir ${D1}
              fi
              if [ ! -d ${DEST_DIR} ]; then
                mkdir ${DEST_DIR}
              fi
           fi


            if [ ${COPY} = 1 ]; then      
              cp iterations-STRUCTURE_1D.txt ${DEST_DIR}/iter_FSI-${N}_${DESCRIPTION}_[${N}_${TAU}_${KAPPA}].txt
              cp convergence-STRUCTURE_1D.txt ${DEST_DIR}/conv_FSI-${N}_${DESCRIPTION}_[${N}_${TAU}_${KAPPA}].txt
              cp postProcessingInfo.txt ${DEST_DIR}/ppInfo_FSI-${N}_${DESCRIPTION}_[${N}_${TAU}_${KAPPA}].txt
              cp preCICE.xml ${DEST_DIR}/
              cp EventTimings.log ${DEST_DIR}/EventTimings_FSI-${N}_${DESCRIPTION}_[${N}_${TAU}_${KAPPA}].log
              cp log.fluid ${DEST_DIR}/
              cp log.solid ${DEST_DIR}/
            fi
        done
    done
 done
done
fi


# POST-PROCESSING OF OUTPUT DATA
if [ ${POSTPROC} = 1 ]; then
    
    python script_postprocessing_iter.py ${DEST_DIR}

    cp iterations_FSI-${N}_*.dat ${DEST_DIR}/
    cp log.pp ${DEST_DIR}/log.pp

    echo "\n -- pgfplots iteration table ---\n\n"
    cat iterations_FSI-${N}_*.dat
    echo "\n -- post procesing log file  ---\n\n"
    cat log.pp

    rm iterations_FSI-${N}_*.dat
    rm log.pp
fi

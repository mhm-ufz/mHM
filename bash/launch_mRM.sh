#!/bin/bash
#
# Evaluate Regionalization for MultiScale Routing Model (mRM)

# author: Stephan Thober
# created: 23.02.2016
# modified: Jun 2016, ST - added German and European basins

# set -x

# ---------------------------------------------------------------------------------------------------------------
# Start
#
set -e
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

function usage () {
    printf "${pprog} [-h] [-a False] [-b basins] [-c do_opti] [-d base_dir] [-e False] [-i in_dir] [-m mRM_exec] [-n namelist_dir] [-p python_dir] [-r res] [-t target_name]\n"
    printf "Run mRM with specific regionalization for different resolutions.\n"
    printf "The following settings can be given\n"
    printf "\n"
    printf "Options\n"
    printf "    -h     Prints this help screen.\n"
    printf "    -a     perform ALL mRM runs, only performed when 'True' is given\n"
    printf "           (Default: True)\n"
    printf "    -b     Do basin dependent namelists (Default: False).\n"
    printf "           The basin-dependent namelists have a fixed naming convention\n"
    printf "           and are located in the folder:\n"
    printf "           (/home/thober/projects/2016/routing/paper/data/mRM_default_params/)\n"
    printf "    -c     Perform basin calibration, otherwise forward run.\n"
    printf "           (Default: True).\n"
    printf "    -d     base directory where outputs are written\n"
    printf "           (Default: /work/thober/mRM)\n"
    printf "    -e     optimize ALL European basins, only performed when 'True' is given\n"
    printf "           (Default: True)\n"
    printf "    -i     input directory for mRM_exec [-m] and file_nml [-n]\n"
    printf "           filenames are assumed to be mhm, mhm_parameter.nml and run_times.sh, respectively\n"
    printf "           (Default: '')\n"
    printf "    -m     specifies the mRM to be used.\n"
    printf "           This file must be located in the base directory.\n"
    printf "           (Default: mrm_r2254)\n"
    printf "    -n     specify the mhm_parameter.nml file\n"
    printf "           (Default: mhm_parameter.nml)\n"
    printf "    -p     directory where python scripts are located\n"
    printf "           for calculation of bar plot and correlation plot\n"
    printf "           (Default: /home/thober/projects/2016/routing/paper/figures/).\n"
    printf "    -r     resolutions to evaluate\n"
    printf "           (Default: 0.5 0.25 0.125)\n"
    printf "    -t     Target directory that will be create as subdirectory\n"
    printf "           in base_dir (Default: '')\n"
    printf "    -z     run time file containing hash values for run times\n"
    printf "           See default table in this script for format.\n"
    printf "           (Default: '')\n"
    printf "\n"
    printf "Example\n"
    printf "    ${pprog} -b /work/thober/mRM -m mrm_for_test -n /work/thober/mRM/tennessee_savannah -p /usr/lib/local/\n"
}

# -----------------------------------------------------------------------------
# FUNCTION DEFINITIONS --------------------------------------------------------
# -----------------------------------------------------------------------------
function create_work () {
    if [[ -d ${work} ]] ; then rm -r  ${work}; fi
    mkdir -p ${work}
    # (1.4) copy namelists ----------------------------------------------------
    cp ${file_nml}           ${work}/mhm_parameter.nml
    cp mhm_template.nml      ${work}/mhm.nml
    cp mrm_outputs.nml       ${work}/mrm_outputs.nml
    cp submit_mRM.sh         ${work}/submit_mRM.sh
    cp run_python_scripts.sh ${work}/run_python_scripts.sh
    # model
    cp ${mRM_exec}           ${work}/mrm
    chmod +x ${work}/mrm
}

# colum where gauge id is located is input to this function
function change_namelists () {
    gauge_col=${1}
    is_EU=${2}
    is_ger=${3}
    datainput=${database}${basin_name}
    datatotalrun=${datainput}/total_runoff/
    datalatlon=${datainput}/latlon/latlon.nc
    # gauging station id
    ngauge=1
    gaugeid=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f ${gauge_col} -d " ")
    # get periods (separate evaluation for german and europe basins)
    if [[ ${is_ger} == 'True' ]]; then
        # get dates
        dStart=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f 13 -d " ")
        yStart=$(echo ${dStart} | cut -d '-' -f 1)
        mStart=$(echo ${dStart} | cut -d '-' -f 2)
        dStart=$(echo ${dStart} | cut -d '-' -f 3)
        dEnd=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f 14 -d " ")
        yEnd=$(echo ${dEnd} | cut -d '-' -f 1)
        mEnd=$(echo ${dEnd} | cut -d '-' -f 2)
        dEnd=$(echo ${dEnd} | cut -d '-' -f 3)
        # check whether there are six years of warmingDays avaialable
        set +e
        (( yDiff = ${yEnd} - ${yStart} ))
        if [[ ${yDiff} -gt 15 ]]; then
            (( yStart = ${yEnd} - 10 ))
        fi
        # set warming days to zero
        warmingDays=0
        set -e
    else
        yStart=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f $(( ${gauge_col} + 1 )) -d " ")
        mStart=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f $(( ${gauge_col} + 2 )) -d " ")
        dStart=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f $(( ${gauge_col} + 3 )) -d " ")
        yEnd=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f $(( ${gauge_col} + 4 )) -d " ")
        mEnd=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f $(( ${gauge_col} + 5 )) -d " ")
        dEnd=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f $(( ${gauge_col} + 6 )) -d " ")
        meteo_yStart=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f 2 -d " ")
        # check whether there are six years of warmingDays avaialable
        set +e
        (( yDiff = ${yStart} - ${meteo_yStart} ))
        if [[ ${yDiff} -gt 4 ]]; then
            warmingDays=1826
        else
            (( warmingDays = 365 * ${yDiff} ))
        fi
        set -e
    fi
    # German basins have more complicated path structure
    if [[ ${is_ger} == 'True' ]]; then
        datainput=${database}sub_${gaugeid}
        datatotalrun=/data/stohyd/mHM_project/germany/basin/small/sub_${gaugeid}/thober/total_runoff/
        datalatlon=${datainput}/latlon/latlon_v5.5.nc
    fi
    # apply changes to mhm.nml
    sed -e "s|STcoordSysST|${coordSys}|" \
        -e "s|STresHydroST|${resHydro}|" -e "s|STresRoutST|${rr}|" \
        -e "s|STdataInputST|${datainput}|" -e "s|SToutST|${work}|" \
        -e "s|STlatlonST|${datalatlon}|" -e "s|STtotalRunST|${datatotalrun}|" \
        -e "s|STfileLCST|${fileLC}|" -e "s|STgaugeIdST|${gaugeid}|" \
        -e "s|STyStartST|${yStart}|" -e "s|STmStartST|${mStart}|" \
        -e "s|STdStartST|${mStart}|" -e "s|STyEndST|${yEnd}|" \
        -e "s|STmEndST|${mEnd}|" -e "s|STdEndST|${dEnd}|" \
        -e "s|SToptiST|${do_opti}|" \
	mhm.nml > mhm.${pid}.nml
    mv mhm.${pid}.nml mhm.nml
        # -e "s|STwarmingDaysST|${warmingDays}|" \
    # add _nc for EU basins
    if [[ ${is_EU} == 'True' ]]; then
        sed -e "s|${datainput}/morph/|${datainput}/morph_v2/|" \
	    mhm.nml > mhm.${pid}.nml
        mv mhm.${pid}.nml mhm.nml
    fi    
    #
    if [[ ${do_opti} == 'True' ]]; then
        # (1.7) create subdirectories for optimization --------------------
        for ii in $(seq ${Nrealization}); do
            # create directory
            run_dir='run_'${ii}/
            mkdir ${run_dir}
            # symlink files
            ln -s ${work}mhm_parameter.nml ${run_dir}
            cp ${work}mhm.nml ${run_dir}
            ln -s ${work}mrm_outputs.nml ${run_dir}
            ln -s ${work}mrm ${run_dir}
            # change run directory
            sed -e "s|\(dir_Out.*${work}\)\(.*\)|\1${run_dir}\2|" \
	        ${run_dir}mhm.nml > mhm.${pid}.nml
            mv mhm.${pid}.nml ${run_dir}mhm.nml
        done
    fi
    # (1.8) submit script ---------------------------------------------
    sed -e "s|STjobnameST|${basin_name%_*}_${rr}|" submit_mRM.sh > sub.${pid}.sh
    mv sub.${pid}.sh submit_mRM.sh
}

deg_to_m() {
    case $1 in
        '0.5') echo '16000';;
        '0.25') echo '8000';;
        '0.125') echo '4000';;
        *) echo $1;;
    esac
}      

deg_to_m_EU() {
    case $1 in
        '0.5') echo '48000';;
        '0.25') echo '24000';;
        '0.125') echo '12000';;
        *) echo $1;;
    esac
}      
# -----------------------------------------------------------------------------
# FUNCTION DEFINITIONS --------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# switches for basins
do_ohio='False'
do_ger='True'
do_EU='False'
# switches for optimization
do_opti='True'
do_basin_nml='False'
do_mrm='True' # set to True if mRM should be run
# default setings
basin_nml_dir='/home/thober/projects/2016/routing/paper/data/mRM_default_params/'
base_dir='/work/thober/mRM/'
mRM_exec='/work/thober/mRM/mrm_r2254'
target_name=''
ptn_dir='/home/thober/projects/2016/routing/paper/figures/'
file_nml='mhm_parameter.nml' # parameter file
in_dir=''
while getopts ":h:a:b:c:d:e:i:m:n:p:r:t:" Option ; do
    case ${Option} in
        h) usage 1>&2; exit 0;;
        a) do_mrm="${OPTARG}";;
        b) do_basin_nml="${OPTARG}";;
        c) do_opti="${OPTARG}";;
        d) base_dir="${OPTARG}";;
        e) do_eu="${OPTARG}";;
        i) in_dir="${OPTARG}";;
        m) mRM_exec="${OPTARG}";;
        n) file_nml="${OPTARG}";;
        r) res="${OPTARG}";;
        t) target_name="${OPTARG}";;
        p) ptn_dir="${OPTARG}";;
        *) printf "Error ${pprog}: unimplemented option.\n\n";  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

# check whether in_dir has been given and reset variables
if [[ ${in_dir} != '' ]]; then
    mRM_exec=${in_dir}/mhm
    file_nml=${in_dir}/mhm_parameter.nml
fi

# determine opti_function to use
opti=$(grep 'opti_function =' ${file_nml%/*}/mhm.nml | cut -d '=' -f 2)
# set Number of realizations
Nrealization=10

# print startup
echo "*******************************"
echo "*** start up initialization ***"
echo "base_dir = "${base_dir}
echo "mRM_exec = "${mRM_exec}
echo "file_nml = "${file_nml}
echo "res = "${res}
echo "target_name = "${target_name}
echo "ptn_dir = "${ptn_dir}
echo "opti_function = "${opti}
echo "*** start up initialization ***"
echo "*******************************"

# consistency checks =========================================================
if [[ ${target_name} == '' ]]; then printf "Error ${pprog}: target_name cannot be empty.\n\n" ; exit 1 ; fi
if [[ ! -f ${mRM_exec} ]] ; then printf "Error ${pprog}: mRM_exec not found.\n\n" ; exit 1 ; fi
if [[ ! -f ${file_nml} ]] ; then printf "Error ${pprog}: ${file_nml} not found.\n\n" ; exit 1 ; fi
if [[ ${mRM_exec:0:1} != '/' ]]; then printf "Error ${pprog}: ${mRM_exec} is no absolute path.\n\n"; exit 1; fi
if [[ ${file_nml:0:1} != '/' ]]; then printf "Error ${pprog}: ${file_nml} is no absolute path.\n\n"; exit 1; fi
if [[ ! -f ${ptn_dir}'/plot_param_box.py' ]];then printf "Error ${pprog}: plot_param_box.py not found in ptn_dir.\n\n" ; exit 1 ; fi
if [[ ! -f ${ptn_dir}'/plot_param_corr.py' ]];then printf "Error ${pprog}: plot_param_corr.py not found in ptn_dir.\n\n" ; exit 1 ; fi
if [[ ! -f ${ptn_dir}'/ParamOpti_development.py' ]];then printf "Error ${pprog}: ParamOpti_development.py not found in ptn_dir.\n\n" ; exit 1 ; fi
if [[ ! -f ${ptn_dir}'/plot_Q_mRM.py' ]];then printf "Error ${pprog}: plot_Q_mRM.py not found in ptn_dir.\n\n" ; exit 1 ; fi


# create dictionary for translating resolution for neckar rockenau

# ============================================================================

# set target_dir
target_dir=${base_dir}'/'${target_name}
if [[ ! -d ${target_dir} ]]; then mkdir -p ${target_dir}; fi

# create configfile
config_file=${target_dir}'/Config.out'
cat > ${config_file} <<EOF
*** THIS IS THE CONFIGFILE OF EVAL_MRM_MPR.sh ***" 
run dir: ${PWD}
base_dir = ${base_dir} 
mRM_exec = ${mRM_exec} 
file_nml = ${file_nml} 
res = ${res} 
target_name = ${target_name} 
ptn_dir = ${ptn_dir} 
opti_function = ${opti} 
EOF

# create prefix for StdOut file
prefix=${target_name:0:1}${target_name#*_}

# create target dir
if [[ ! -d ${target_dir} ]]; then mkdir -p ${target_dir}; fi

# =============================================================================
# (1) PROCESS OHIO BASINS
# =============================================================================
if [[ ${do_ohio} == 'True' ]]; then
    # (0.) settings for US ----------------------------------------------------
    res='0.5 0.25 0.125 0.0625'
    basin_file='/data/stohyd/data/processed/USA/mopex/LUT_ohio_basins.txt'
    database='/data/stohyd/data/processed/USA/mopex/'
    resHydro='0.125'
    fileLC='lc_2001.asc'
    coordSys=1
    ncol=17 # number of columns in lut table
    # (1.1) determine number of basins ----------------------------------------
    job_ids='' # initialize job ids
    nbasin=$(($(cat ${basin_file} | wc -l)-1))
    # resolution loop
    for rr in ${res}; do
        echo 'resolution: '${rr}
        # basin loop
        for ((i=1; i<=${nbasin}; i++)) ; do
            # (1.2) get the basin name and verify whether it should be processed
            basin_name=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f 8 -d " ")
            echo -n "${basin_name} "
            # last colunm (FLAG) in lookup table is 0: exclude basin or 1: include basin
            flag=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f ${ncol} -d " ")
            if [[ ${flag} -eq 0 ]] ; then
	        echo "skipped."
	        continue
            fi
            # (1.3) create working directory ----------------------------------
            work=${target_dir}/${basin_name}/res_${rr}/
            if [[ ${do_basin_nml} == 'True' ]]; then
                file_nml=${basin_nml_dir}/${basin_name}_res_0.0625_run_5_FinalParam.nml
            fi
            create_work
            # (1.5) change dir into the work directory ------------------------
            cd ${work}
            # (1.6) change mhm namelist ---------------------------------------
            change_namelists 9 False False # nine is the column of the gauge id
            # (1.9) submit routing model --------------------------------------
            basin_list="${basin_list} ${basin_name}"
            if [[ ${do_opti} == 'True' ]]; then
                job_id=$(qsub -q ideve_120 -terse -t 1-${Nrealization} submit_mRM.sh)
            else
                job_id=$(qsub -q ideve_120 -terse submit_mRM.sh)
            fi
            # sleep for 5s such that qstat gets all jobs
            # printf "    waiting for 2s such that job is registered by qstat\n"
            # sleep 2
            # # submit python jobs for this resolution
            # qsub -hold_jid ${job_id%.*} run_python_scripts.sh -d ${target_dir}/ -b ${basin_name} -e ${Nrealization} -r ${rr} -p ${ptn_dir} -n ${file_nml} > /dev/null
            # append res_jid to basin_jids
            job_ids=${job_ids}${job_jid%.*},
            # (1.10) change back into original directory ----------------------
            cd ${isdir}
        done
        echo ''
    done
fi

# ============================================================================
# (2) PROCESS GERMAN BASINS
# ============================================================================
if [[ ${do_ger} == 'True' ]]; then
    # (0.) settings for German basins ----------------------------------------
    res='16000 8000' # '16000 8000 4000 2000 1000'
    # flag does not exist in german LUT, all basins are considered
    basin_file='/home/thober/projects/2016/routing/launch_scripts/LUT_german_basins_Matze.txt'
    database='/data/stohyd/data/processed/Germany/basins/'
    resHydro='4000'
    fileLC='lc_1990.asc'
    coordSys=0
    ncol=17 # number of columns in lut table
    # (1.1) determine number of basins ----------------------------------------
    job_ids='' # initialize job ids
    nbasin=$(($(cat ${basin_file} | wc -l)-1))
    # resolution loop
    for rr in ${res}; do
        rr=$(deg_to_m ${rr})
        echo 'resolution: '${rr}
        # basin loop
        for ((i=1; i<=${nbasin}; i++)) ; do
            # (1.2) get the basin name and verify whether it should be processed --
            basin_name=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f 18,19 -d " " | tr ' ' '_')
            basin_id=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f 1 -d " ")
            echo -n "${basin_name} "
            # (1.3) create working directory --------------------------------------
            work=${target_dir}/sub_${basin_id}/res_${rr}/
            if [[ ${do_basin_nml} == 'True' ]]; then
                file_nml=${basin_nml_dir}/sub_${basin_id}_res_1000_run_5_FinalParam.nml
            fi
            create_work
            # (1.5) change dir into the work directory ---------------------------
            cd ${work}
            # (1.6) change mhm namelist ------------------------------------------
            change_namelists 1 False True
            # (1.9) submit routing model --------------------------------------
            basin_list="${basin_list} ${basin_name}"
            if [[ ${do_opti} == 'True' ]]; then
                job_id=$(qsub -q ideve_120 -terse -t 1-${Nrealization} submit_mRM.sh)
            else
                job_id=$(qsub -q ideve_120 -terse submit_mRM.sh)
            fi
            # # sleep for 5s such that qstat gets all jobs
            # printf "    waiting for 2s such that job is registered by qstat\n"
            # sleep 2
            # # submit python jobs for this resolution
            # qsub -hold_jid ${job_id%.*} run_python_scripts.sh -d ${target_dir}/ -b ${basin_name} -e ${Nrealization} -r ${rr} -p ${ptn_dir} -n ${file_nml} > /dev/null
            # append res_jid to basin_jids
            job_ids=${job_ids}${job_jid%.*},
            # (1.10) change back into original directory ----------------------
            cd ${isdir}
        done
        echo ''
    done
fi

# ============================================================================
# (3) PROCESS EUROPEAN BASINS
# ============================================================================
if [[ ${do_EU} == 'True' ]]; then
    # (0.) settings for European basins --------------------------------------
    res='48000 24000 12000 6000 3000'
    basin_file='/home/thober/projects/2016/routing/launch_scripts/europe_basin_lut_JM_cross-valid.txt'
    database='/data/stohyd/data/processed/europe/sub_basins/'
    resHydro='24000'
    fileLC='lc_1990.asc'
    coordSys=0 
    ncol=37 # number of columns in lut table
    # (1.1) determine number of basins ----------------------------------------
    job_ids='' # initialize job ids
    nbasin=$(($(cat ${basin_file} | wc -l)-1))
    # resolution loop
    for rr in ${res}; do
        rr=$(deg_to_m_EU ${rr})
        echo 'resolution: '${rr}
        # basin loop
        for ((i=1; i<=${nbasin}; i++)) ; do
            # (1.2) get the basin name and verify whether it should be processed --
            basin_name=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f 8 -d " ")
            echo -n "${basin_name} "
            # last colunm (FLAG) in lookup table is 0: exclude basin or 1: include basin
            flag=$(cat ${basin_file} | sed "1d" | sed -n "${i}p" | cut -f ${ncol} -d " ")
            if [[ ${flag} -eq 0 ]] ; then
	        echo "skipped."
	        continue
            fi
            # (1.3) create working directory --------------------------------------
            work=${target_dir}/${basin_name}/res_${rr}/
            if [[ ${do_basin_nml} == 'True' ]]; then
                file_nml=${basin_nml_dir}/${basin_name}_res_48000_run_5_FinalParam.nml
            fi
            create_work
            # (1.5) change dir into the work directory ---------------------------
            cd ${work}
            # (1.6) change mhm namelist ------------------------------------------
            change_namelists 12 True False
            # (1.9) submit routing model --------------------------------------
            basin_list="${basin_list} ${basin_name}"
            if [[ ${do_opti} == 'True' ]]; then
                job_id=$(qsub -q ideve_120 -terse -t 1-${Nrealization} submit_mRM.sh)
            else
                job_id=$(qsub -q ideve_120 -terse submit_mRM.sh)
                # chmod +x submit_mRM.sh
                # ./submit_mRM.sh
                # time ./mrm > mrm.out
            fi
            # # sleep for 5s such that qstat gets all jobs
            # printf "    waiting for 2s such that job is registered by qstat\n"
            # sleep 2
            # # submit python jobs for this resolution
            # qsub -hold_jid ${job_id%.*} run_python_scripts.sh -d ${target_dir}/ -b ${basin_name} -e ${Nrealization} -r ${rr} -p ${ptn_dir} -n ${file_nml} > /dev/null
            # append res_jid to basin_jids
            job_ids=${job_ids}${job_jid%.*},
            # (1.10) change back into original directory ----------------------
            cd ${isdir}
        done
        echo ''
    done
fi
printf "Done!\n"
exit 0

# ============================================================================
# WHAT FOLLOWS, IS OLD CODE ==================================================
# ============================================================================

# # ============================================================================
# # SETUP HASH TABLE FOR RUN TIMES
# hash_index() {
#     case $1 in
#         'ohio_metropolis')
#                 case $2 in
#                      '0.5') return 1;;
#                      '0.25') return 2;;
#                      '0.125') return 3;;
#                 esac;;
#         'wabash_carmel')
#                 case $2 in
#                      '0.5') return 4;;
#                      '0.25') return 5;;
#                      '0.125') return 6;;
#                 esac;;
#         'ohio_louisville')
#                 case $2 in
#                      '0.5') return 7;;
#                      '0.25') return 8;;
#                      '0.125') return 9;;
#                 esac;;
#         'tennessee_savannah')
#                 case $2 in
#                      '0.5') return 10;;
#                      '0.25') return 11;;
#                      '0.125') return 12;;
#                 esac;;
#         'neckar_rockenau')
#                 case $2 in
#                      '16000') return 13;;
#                      '8000') return 14;;
#                      '4000') return 15;;
#                 esac;;
#     esac
# }

# if [[ ! -f ${rt_file} ]]; then
#     printf "Warning ${pprog}: ${rt_file} not found, using default table.\n"
#     run_times=("", # sort of like returning null/nil for a non existent key
#                "24:00:00" # ohio_metropolis 0.5
#                "48:00:00" # ohio_metropolis 0.25
#                "96:00:00" # ohio_metropolis 0.125
#                "01:30:00" # wabash_carmel 0.5
#                "03:00:00" # wabash_carmel 0.25
#                "06:00:00" # wabash_carmel 0.125
#                "02:00:00" # ohio_louisville 0.5
#                "04:00:00" # ohio_louisville 0.25
#                "08:00:00" # ohio_louisville 0.125
#                "01:00:00" # tennessee_savannah 0.5
#                "02:00:00" # tennessee_savannah 0.25
#                "04:00:00" # tennessee_savannah 0.125
#                "01:00:00" # neckar_rockenau 16000
#                "02:00:00" # neckar_rockenau 8000
#                "04:00:00" # neckar_rockenau 4000
#               );
# else
#     # include run times
#     source ${rt_file}
# fi
# # initialize job ids
# all_ids=""
# # loop over basins
# for bb in ${basins};do
#     printf "  processing basin: "${bb}"\n"
#     # set nml_file
#     nml_file=${base_dir}'/'${bb}'/namelists/mhm.nml'
#     # consistency checks
#     if [[ ! -f ${nml_file} ]];then printf "Error ${pprog}: mhm.nml not found in nml_dir.\n\n" ; exit 1 ; fi
#     # set basin_ids
#     basin_ids=""
#     # loop over resolutions
#     for rr in ${res}; do
#         # translate resolution if neckar_rockenau is evaluated
#         if [[ ${bb} == 'neckar_rockenau' ]]; then
#             rr=$(deg_to_m ${rr})
#             rr_name=${rr:0:2}
#         else
#             rr_name=${rr:2:3}
#         fi
#         #
#         printf "    processing resolution: "${rr}"\n"
#         if [[ ${do_mrm} == 'True' ]]; then
#             # set target_dir
#             run_dir=${target_dir}'/'${bb}'/res_'${rr}'/run_'
#             # create target dirs
#             for ii in $(seq ${Nrealization}); do
#                 tmp_dir=${run_dir}${ii}
#                 if [[ ${tmp_dir} == '' ]]; then printf "Error ${pprog}: ${BASE_DIR}${ii} is empty.\n\n" ; exit 1 ; fi
#                 if [[ ${tmp_dir} == '/' ]]; then printf "Error ${pprog}: ${BASE_DIR}${ii} is /.\n\n" ; exit 1 ; fi
#                 if [[ ${tmp_dir/*\/usr\/*/TRUE} == 'TRUE' ]]; then printf "Error ${pprog}: ${BASE_DIR}${ii} contains '/usr/'.\n\n" ; exit 1 ; fi
#                 if [[ -d ${tmp_dir} ]]; then rm -r ${tmp_dir}; fi
#                 mkdir -p ${tmp_dir}
#                 # copy files there
#                 ln -s ${file_nml} ${tmp_dir}/mhm_parameter.nml
#                 cp mrm_outputs.nml ${tmp_dir}
#                 # replace Routing Resolution in mhm.nml
#                 sed -e "s/resolution_Routing(1) = .*/resolution_Routing(1) = ${rr}/g" \
#                     -e "s/opti_function = .*/opti_function = ${opti}/g" ${nml_file} > ${tmp_dir}/mhm.nml
#                 # copy executable
#                 ln -s ${mRM_exec} ${tmp_dir}/mrm
#             done

#             # get run time
#             rt=$(hash_index ${bb} ${rr} || echo ${run_times[$?]})
#             # create name of StdOut File
#             gauge=${bb#*_}
#             # substitute run time and basin name in submit_mRM.sh
#             sed -i -e "s/#$ -l h_rt=.*/#$ -l h_rt=${rt}/g" -e "s/#$ -N .*/#$ -N ${prefix}_${gauge:0:1}_${rr_name}/g" submit_mRM.sh
#             # submit mRM optimization and obtain job id for that resolution
#             res_jid=$(qsub -terse -t 1-10 submit_mRM.sh ${run_dir})
#             # sleep for 5s such that qstat gets all jobs
#             printf "    waiting for 5s such that job is registered by qstat\n"
#             sleep 5
#             # submit python jobs for this resolution
#             qsub -hold_jid ${res_jid%.*} run_python_scripts.sh -d ${target_dir}/ -b ${bb} -r ${rr} -e ${Nrealization} -p ${ptn_dir} -o ${target_dir}'/'${bb}'/res_'${rr} -n ${file_nml} > /dev/null
#             # append res_jid to basin_jids
#             basin_ids=${basin_ids}${res_jid%.*},
#         else
#             ./run_python_scripts.sh -d ${target_dir}/ -b ${bb} -r ${rr} -e ${Nrealization} -p ${ptn_dir} -o ${target_dir}'/'${bb}'/res_'${rr} -n ${file_nml}
#         fi
#     # end resolution loop
#     done
#     if [[ ${do_mrm} == 'True' ]]; then
#         # submit python job for this basin, give -r flag to each resolutions
#         qsub -hold_jid ${basin_ids} run_python_scripts.sh -d ${target_dir}/ -b ${bb} -r $(echo ${res} | sed -e 's/ / -r /g') -e ${Nrealization} -p ${ptn_dir} -o ${target_dir}'/'${bb} -n ${file_nml} > /dev/null
#         # append basin ids to all ids
#         all_ids=${all_ids}${basin_ids}
#     else
#         ./run_python_scripts.sh -d ${target_dir}/ -b ${bb} -r $(echo ${res} | sed -e 's/ / -r /g') -e ${Nrealization} -p ${ptn_dir} -o ${target_dir}'/'${bb} -n ${file_nml}
#     fi
# # end basin loops
# done
# printf "submitted all mRM runs\n"
# if [[ ${do_mrm} == 'True' ]]; then
#     # submit python jobs for all basins and resolutions, give -b flag to each basin, give -r flag to each resolutions
#     qsub -hold_jid ${all_ids} run_python_scripts.sh -d ${target_dir}/ -b ${basins/ / -b } -r ${res/ / -r } -e ${Nrealization} -p ${ptn_dir} -o ${target_dir} -n ${file_nml}
# else
#     ./run_python_scripts.sh -d ${target_dir}/ -b $(echo ${basins} | sed -e 's/ / -b /g') -r $(echo ${res} | sed -e 's/ / -r /g') -e ${Nrealization} -p ${ptn_dir} -o ${target_dir} -n ${file_nml}
# fi    
# printf "Done!\n"

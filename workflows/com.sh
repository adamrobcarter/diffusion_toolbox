set -e

# bash workflows/com.sh sim_nohydro_0.1_L*_t1e4_unwrap

python -m MSD.calc_centre_of_mass_entire ${1}
python -m MSD.show centre_of_mass_entire_${1}
python -m MSD.calc ${1}
python -m MSD.show ${1}
python -m visualisation.Ds_msd_window ${1}
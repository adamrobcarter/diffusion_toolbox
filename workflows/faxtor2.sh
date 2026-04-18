set -e

python -m preprocessing.faxtor2 --sample=${1} --exp=${2} --meas=${3}
python -m preprocessing.remove_moving_avg faxtor2_${1}_${2}_${3}
python -m particle_detection.detect faxtor2_${1}_${2}_${3}_movavrem
python -m particle_linking.link faxtor2_${1}_${2}_${3}_movavrem
python -m particle_linking.show_movie faxtor2_${1}_${2}_${3}_movavrem
python -m stepsize.calc faxtor2_${1}_${2}_${3}_movavrem

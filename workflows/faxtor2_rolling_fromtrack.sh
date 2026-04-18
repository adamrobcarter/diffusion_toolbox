set -e

python -m particle_detection.detect faxtor2_${1}_${2}_${3}_movavrem_rolling5
python -m particle_detection.show_movie faxtor2_${1}_${2}_${3}_movavrem_rolling5
python -m particle_linking.link faxtor2_${1}_${2}_${3}_movavrem_rolling5
python -m particle_linking.show_movie faxtor2_${1}_${2}_${3}_movavrem_rolling5
python -m stepsize.calc faxtor2_${1}_${2}_${3}_movavrem_rolling5

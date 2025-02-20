set -e # stop on error
python -m particle_linking.link $1
python -m MSD.calc $1
python -m MSD.show $1
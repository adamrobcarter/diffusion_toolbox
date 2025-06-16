set -e # stop on error

python -m preprocessing.hpf $1
python -m preprocessing.remove_moving_avg ${1}_hpf
python -m preprocessing.to_tiff ${1}_hpf_movavrem
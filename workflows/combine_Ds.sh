set -e # stop on error
python -m box_counting.D_of_L $1
python -m box_counting.D_of_L ${1}_longer
python -m visualisation.merge_Ds $1 ${1}_longer
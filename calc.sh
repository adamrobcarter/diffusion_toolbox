python -m particle_linking.link "$1_div8"
python -m MSD.calc "$1_div8"
python -m isf.calc_f $1
python -m box_counting.count $1
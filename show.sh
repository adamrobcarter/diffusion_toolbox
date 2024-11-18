python -m MSD.show "$1_div8"
python -m isf.show_f $1
python -m box_counting.msd_single $1
python -m box_counting.D_of_L $1
python -m visualisation.Ds_overlapped $1
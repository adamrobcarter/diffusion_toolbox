set -e

# python -m particle_detection.singlet_to_old sim_nohydro_011_L320_test_singlet

python -m isf.calc_f_first sim_nohydro_011_L320_test_singlet
# python -m isf.calc_f_first sim_nohydro_011_L320_test_singlet_unmix
python -m isf.calc_f_first sim_nohydro_011_L320_test_mixt

python -m visualisation.Ds_overlapped_mult sim_nohydro_011_L320_test_unmix_longer sim_nohydro_011_L320_test_mixt
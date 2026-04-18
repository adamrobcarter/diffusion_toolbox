import common
import tifffile
import numpy as np
import argparse

parser = argparse.ArgumentParser()

def len4(s):
        if len(s) != 4:
                raise Exception("sample/exp/meas should have lenth 4")
        return s

parser.add_argument("--sample", type=len4)
parser.add_argument("--exp",    type=len4)
parser.add_argument("--meas",   type=len4)
args = parser.parse_args()

basepath = '/beamlines/bl31/projects/cycle2026-I/20250370300-lmichot/DATA/PROCESSED/'

dir = f"{basepath}SAMPLE_{args.sample}/EXPERIMENT_{args.exp}/MEASUREMENT_{args.meas}/PCO_EDGE/"

path = f"{dir}corrected_data.tiff"

stack = tifffile.imread(path)

stack = np.swapaxes(stack, 2, 1)

common.save_data(f"preprocessing/data/stack_faxtor2_{args.sample}_{args.exp}_{args.meas}.npz", stack=stack, pixel_size=0.65, time_step=1)

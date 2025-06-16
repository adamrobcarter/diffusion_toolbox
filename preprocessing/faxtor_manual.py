import common
import numpy as np
import tqdm

sample_files = [f'/data2/acarter/faxtor/SAMPLE_0025/EXPERIMENT_0000/MEASUREMENT_{str(n).rjust(4, "0")}/PCO_EDGE/PCO_EDGE_radiography_0000.tiff' for n in range(2, 79, 2)]
flat_files   = [f'/data2/acarter/faxtor/SAMPLE_0025/EXPERIMENT_0000/MEASUREMENT_{str(n).rjust(4, "0")}/PCO_EDGE/PCO_EDGE_radiography_0000.tiff' for n in range(3, 80, 2)]   

sample_tiffs = [common.load_tif(filename) for filename in tqdm.tqdm(sample_files)]
flat_tiffs   = [common.load_tif(filename) for filename in tqdm.tqdm(flat_files)]

samples = np.stack(sample_tiffs, axis=0)
flats   = np.stack(flat_tiffs,   axis=0)

stack = np.log( samples / flats )
print(common.nanfrac(stack), 'nanfrac')
assert np.isfinite(stack).all()

name = 'faxtor025b'
time_step = 16
kwargs = dict(
    NAME = 'Si 120um, PS 8um 2g/L injected'
)
common.save_data(f'preprocessing/data/stack_{name}.npz',       stack=stack,      pixel_size=0.65, time_step=time_step, **kwargs)
common.save_data(f'preprocessing/data/stack_{name}_small.npz', stack=stack[:50], pixel_size=0.65, time_step=time_step, **kwargs)
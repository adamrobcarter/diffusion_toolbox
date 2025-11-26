import common

if __name__ == '__main__':
    data = common.load_tif('raw_data/Project_Test4_ch00.tif')
    common.save_data('preprocessing/data/stack_hermann', stack=data, time_step=0.2, pixel_size=7.265625)
import countoscope as countoscope_new
import countoscope_old
import numpy as np
import common


system_size = np.array([30, 30, 30])

box_size=np.array([10, 10, 10])
results = countoscope_new.Countoscope(trajectory_file = path_to_data,
                    system_size=system_size,
                    box_size=box_size)
results.run()
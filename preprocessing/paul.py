import numpy as np
import common
import pandas

"""
LABEL,ID,TRACK_ID,QUALITY,POSITION_X,POSITION_Y,POSITION_Z,POSITION_T,FRAME,RADIUS,VISIBILITY,MANUAL_SPOT_COLOR,MEAN_INTENSITY_CH1,MEDIAN_INTENSITY_CH1,MIN_INTENSITY_CH1,MAX_INTENSITY_CH1,TOTAL_INTENSITY_CH1,STD_INTENSITY_CH1,CONTRAST_CH1,SNR_CH1
Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame,Radius,Visibility,Manual spot color,Mean intensity ch1,Median intensity ch1,Min intensity ch1,Max intensity ch1,Sum intensity ch1,Std intensity ch1,Contrast ch1,Signal/Noise ratio ch1
Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame,R,Visibility,Spot color,Mean ch1,Median ch1,Min ch1,Max ch1,Sum ch1,Std ch1,Ctrst ch1,SNR ch1
,,,(quality),(pixel),(pixel),(pixel),(frame),,(pixel),,,(counts),(counts),(counts),(counts),(counts),(counts),,
ID77827,77827,1,3.338761806488037,1323.3559827742383,1057.9013128675044,0.0,497.0,497,5.0,1,,135.83505154639175,129.0,113.0,169.0,13176.0,14.921796426404745,0.056441066945738146,0.972681546460471

"""

# data = np.loadtxt('raw_data/paul/DSCF5574_Active_Demo_spots.csv', skiprows=4, delimiter=',')

# particles = data[:, [4, 5, 6, 7, 1]]
# print(particles)

df = pandas.read_csv('raw_data/paul/DSCF5574_Active_Demo_spots.csv', skiprows=[1,2,3], header=0)

df = df.drop(columns=['RADIUS', 'VISIBILITY', 'MANUAL_SPOT_COLOR', 'MEAN_INTENSITY_CH1', 'MEDIAN_INTENSITY_CH1', 'MIN_INTENSITY_CH1', 'MAX_INTENSITY_CH1', 'TOTAL_INTENSITY_CH1', 'STD_INTENSITY_CH1', 'CONTRAST_CH1', 'SNR_CH1'])

print(df.describe())

data = df[['POSITION_X', 'POSITION_Y', 'FRAME', 'TRACK_ID']].to_numpy()
data[:, 2] -= data[:, 2].min()
data[:, 3] -= data[:, 3].min()

common.save_particles('paul0', data, time_step=1)
common.save_trajs    ('paul0', data, time_step=1)
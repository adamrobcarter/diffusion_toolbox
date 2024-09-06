import common
from preprocessing.marine_to_stack import do_file
          
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'PPR-')):
    do_file(i, file, '2_A', f's1, red {i}')
        
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'DDTPonly-')):
    do_file(i, file, '2_B', f's1, green {i}')
        
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'DDTP-PPR-')): # 2-color
    do_file(i, file, '2_C', f's1, red&green {i}')
        
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'MFLonly-')):
    do_file(i, file, '2_D', f's2, green {i}')
          
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'MFL-PPR-')): # 2-color
    do_file(i, file, '2_E', f's2, red&green {i}')
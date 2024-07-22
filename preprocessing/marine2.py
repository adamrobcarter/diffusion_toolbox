import common
from preprocessing.marine_to_stack import do_file
          
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'PPR-')):
    do_file(i, file, '2_A')
        
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'DDTPonly-')):
    do_file(i, file, '2_B')
        
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'DDTP-PPR-')): # 2-color
    do_file(i, file, '2_C')
        
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'MFLonly-')):
    do_file(i, file, '2_D')
          
for i, file in enumerate(common.get_directory_files('raw_data/marine2', 'czi', 'MFL-PPR-')): # 2-color
    do_file(i, file, '2_E')
import os
import common

if __name__ == '__main__':
    for folder in sorted(os.listdir('/media/com-psiche/Sans titre/psiche_export_1/')):
        try:
            with open(f'/media/com-psiche/Sans titre/psiche_export_1/{folder}/NAME', 'r') as namefile:
                name = namefile.read()
            print(folder, name)
        except FileNotFoundError:
            print(folder, common.name(folder))
import common
import matplotlib.pyplot as plt

for file in common.files_from_argv('greenkubo/data', 'Dcm_'):
    data = common.load(f'greenkubo/data/Dc_{file}.npz')

    Dcm = data['Dcm']
    print(Dcm)

    plt.plot(Dcm)
    plt.loglog()
    plt.savefig(f'greenkubo/figures_png/Dc_{file}.png')
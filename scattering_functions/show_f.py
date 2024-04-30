import common
import scattering

for file in sys.argv[1:]:
    
    for type_index, Ftype in enumerate(['f']):
    # for type_index, Ftype in enumerate(['Fs', 'f', 'DDM']):

        show_single_F_type(file, type_index, Ftype)

            
    plt.suptitle(fr'F or F_s (k, t), {file}')
    plt.tight_layout()

    common.save_fig(plt.gcf(), f'scattering_functions/figures_png/Fs_decay_t_{file}.png', dpi=300)
        
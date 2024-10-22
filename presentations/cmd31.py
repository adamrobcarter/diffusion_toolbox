

import MSD.show
MSD.show.go('eleanor0.01', SHOW_FIT=False, export_destination='/home/acarter/presentations/cmd31/figures')
MSD.show.go('eleanor0.01', SHOW_FIT=True)

import scattering_functions.show_S_of_k
scattering_functions.show_S_of_k.go('eleanor0.34')

import van_hove.show_G
van_hove.show_G.go('eleanor0.34')

import van_hove.show_g
van_hove.show_g.go('eleanor0.34')

import visualisation.Ds_hack

import scattering_functions.show_Fs_overlayed
scattering_functions.show_Fs_overlayed.go('eleanor0.01', SHOW_FIT=False)
scattering_functions.show_Fs_overlayed.go('eleanor0.01', SHOW_FIT=True)

# import visualisation.Ds_overlapped

import box_counting.showN1N2
box_counting.showN1N2.go('alice0.02_nodrift')
box_counting.showN1N2.go('alice0.02_drift0.03')

    """
    (['boxcounting_collective'             ], False, True, True, True),
    (['boxcounting_collective', 'MSD_short'], False, True, True, True),
    (['f_short', 'MSD_short'], True, True, True, True),
    (['f_short'             ], True, True, True, True),
    (['boxcounting_collective', 'f_short'], False, True, True, False),
    (['boxcounting_collective', 'f_short'], False, False, True, False),
    (['boxcounting_shorttime', 'MSD_short'], False, True, False, False),
    (['boxcounting_shorttime'             ], False, True, False, False),
    """
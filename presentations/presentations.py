import matplotlib.image
import matplotlib.offsetbox

def show_png_on_axis(ax, file, coords, size, color='black', data_coords=True):
    image = matplotlib.image.imread(file)
    imagebox = matplotlib.offsetbox.OffsetImage(image, zoom=size)   

    ab = matplotlib.offsetbox.AnnotationBbox(
        imagebox,
        coords,        # position in DATA units
        xycoords="data" if data_coords else "axes fraction", # <-- use data coordinates
        # frameon=False,
    )
    
    frame = ab.patch
    frame.set_edgecolor(color)     # frame/border color
    # frame.set_linewidth(2)         # thickness

    ax.add_artist(ab)
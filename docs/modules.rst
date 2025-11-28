API Reference
=============
The library is organised into various modules for processing and analysing videos and particle positions, listed below. 
Scripts are typically run as ``python -m <module>.<script> [dataset_name]`` from the toolbox root directory.
Do not include the ``.py`` suffix, note that the ``-m`` flag is required to run modules as scripts, and the module and script are separated by a dot not a slash.

Scripts look for data in the relevelant ``<module>/data`` folder, and save data in the data folder in their directory. 
Figure are saved in the relevent ``<module>/figures`` folder (older modules use ``figures_png``).

Data are saved in `npz format <https://numpy.org/doc/2.2/reference/generated/numpy.savez.html>`__.
This is a file that can contain multiple named numpy arrays.
There are three common formats:

- ``stack_<dataset_name>.npz`` contains a video - a stack of frames. Should contain the following objects, but may have more metadata:
  
  - ``stack``: 3D array of shape (num_frames, height, width)
  - ``time_step``: time between frames in seconds
  - ``pixel_size``: size of one pixel in microns

- ``particles_<dataset_name>.npz`` contains particle positions and properties. Should contain the following objects, but may have more metadata:
  
  - ``particles``: 2D array containing rows of (x, y, frame_number) for each detected particle
  - ``time_step``: time between frames in seconds

- ``trajs_<dataset_name>.npz`` contains linked particle positions. Has the same format as ``particles_`` files, but ``particles`` is rows of (x, y, frame_number, particle_id)

The modules are:

.. toctree::
   :maxdepth: 1
   :glob:

   api/*
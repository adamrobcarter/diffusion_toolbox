import common
import numpy as np
import matplotlib.pyplot as plt
import tqdm

for file in common.files_from_argv("particle_linking/data", "trajs"):
	data = common.load(f"particle_linking/data/trajs_{file}.npz")
	particles = data["particles"]

	timecol = common.get_particles_column("t", data)
	idcol   = common.get_particles_column("id", data)

	# sort by ID then time
	print('sorting')
	# particles = particles[np.lexsort((particles[:, time_column], particles[:, id_column]))]
	# chatgpt says the below is faster and more memory-efficient
	particles = particles[np.argsort(particles[:, timecol], kind='mergesort')]
	particles = particles[np.argsort(particles[:, idcol], kind='mergesort')]

	currentid=-1
	startrow=0

	xvar = []
	ymean = []

	for row in tqdm.trange(particles.shape[0]):
		if particles[row, idcol] == currentid:
			pass
		else:
			# we found a new particle, so first deal with the old
			x = particles[startrow:row-1, 0]
			y = particles[startrow:row-1, 1]
			#print(currentid, x.size)

			if x.size > 20:
				xvar.append(np.std(x))
				ymean.append((y[-1]-y[0])/y.size)

			# now move on
			currentid = particles[row, idcol]
			startrow = row

	fig, ax = plt.subplots()
	ax.scatter(xvar, ymean)
	ax.set_xlabel("x variance")
	ax.set_ylabel("y mean velocity")
	ax.set_xlim(0, 50)
	common.save_fig(fig, f"particle_linking/figures_png/yvelxspread_{file}.png")

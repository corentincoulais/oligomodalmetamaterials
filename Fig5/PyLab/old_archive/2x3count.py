"Aleksi Bossart, 18.4.2019"
"plot all essentially different 2x3 tilings"

import numpy as np
from matplotlib import pyplot as plt


"Mask to forbid NNNN interactions"

mask = np.array([[1,1,0,1,1,0],
                 [0,1,1,0,1,1],
                 [1,0,1,1,0,1],
                 [1,1,0,1,1,0],
                 [0,1,1,0,1,1],
                 [1,0,1,1,0,1]])


"Fixing the first tile leaves 1024 supercells"

supercells = 0 * np.stack([mask for _ in range(1024)])
supercells[:, 0, 0] += 1


"initialise all configurations"

s= 0

for i1 in np.nonzero(mask[:, 1])[0]:

	for i2 in np.nonzero(mask[:, 2])[0]:

		for i3 in np.nonzero(mask[:, 3])[0]:

			for i4 in np.nonzero(mask[:, 4])[0]:

				for i5 in np.nonzero(mask[:, 5])[0]:

					supercells[s, i1, 1] += 1
					supercells[s, i2, 2] += 1
					supercells[s, i3, 3] += 1
					supercells[s, i4, 4] += 1
					supercells[s, i5, 5] += 1

					s += 1


"remove open-ended diagonal edges"

for i in range(1024):
	supercells[i] = np.matmul(np.diag( (np.sum(supercells[i], axis=1) ) > 1 ).astype(int), supercells[i])


"remove duplicates"

supercells = np.unique(supercells, axis = 0)
print(np.shape(supercells))


"account for translation symmetries"

for i in range(int(np.shape(supercells)[0])):

	supercells[i] = supercells[i] * int(np.logical_not(np.array_equal(supercells[i, :, :3], np.roll(supercells[i, :, 3:], 3, axis = 0))))
	supercells[i] = supercells[i] * int(np.logical_not(np.array_equal(supercells[i, :, :3], np.roll(supercells[i, :, 3:], 3, axis = 0))))
	if np.logical_and(np.array_equal(supercells[i, :, (0,3)], np.roll(supercells[i, :, (1,4)], 2, axis = 1)), np.array_equal(np.roll(supercells[i, :, (1,4)], 2, axis = 1), np.roll(supercells[i, :, (2,5)], 2, axis = 1))):
		print(supercells[i])

"account for mirror symmetries"

"remove duplicates"

supercells = np.unique(supercells, axis = 0)
print(np.shape(supercells))


from __future__ import division
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

file_beg = 'out_p-'
max_trials = 5

n = range(10000, 1000001, 10000)
# times = np.array([])
p_times = []
lines = ['-b', '-r', '-g', '-k', '-c', '-m', '-y', '-0.75']

for p, c in zip(xrange(6, 49, 6), xrange(0,8)):
    with open(file_beg + '{}.txt'.format(p)) as f:
        trial = 0
        ave = 0
        for line in f:
            _ __ time s = line.split()
            if trial < max_trials:
                ave += float(time)
                trial += 1

            if trial == max_trials:
                ave /= max_trials
                p_times.append(ave)
                trial = 0
                ave = 0

    # times = np.append(times, p_times)

    plt.plot(n, p_times, c)

plt.xlabel('n')
plt.ylabel('Time (s)')
plt.savefig('n_vs_t_all-p.pdf')


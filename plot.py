import numpy as np
import json
from matplotlib import pyplot as plt

fig = plt.figure()
ax_conv = fig.add_subplot(4,1,1)
ax_vel = fig.add_subplot(4,1,2)
ax_acc = fig.add_subplot(4, 1, 3)


with open('world.json') as fh:
    world = json.load(fh)
    bx, by = np.array(world['bounds']).T
    ax_conv.plot(bx, by, 'k--', label="bounding box")

with open('result.json') as fh:
    result = json.load(fh)

    for thing in ('initialization', 'solution'):
        solution = result[thing]
        x, y = np.array(solution).reshape(-1, 2).T
        ax_conv.plot(x, y, 'o-', label=thing)
       
        # now "solution" points to the last one (the actual solution)
        vel = np.diff(np.array(solution).reshape(-1, 2), axis=0)
        acc = np.diff(vel, axis=0)
        acc_x, acc_y = acc.T
        ax_acc.plot(acc_x, label="acc x %s" % thing)
        ax_acc.plot(acc_y, label="acc y %s" % thing)
        vx, vy = vel.T
        ax_vel.plot(vy, label="vel x %s" % thing)
        ax_vel.plot(vx, label="vel y %s" % thing)

for ax in (ax_acc, ax_conv, ax_vel):
    ax.legend()
    
    
plt.show()
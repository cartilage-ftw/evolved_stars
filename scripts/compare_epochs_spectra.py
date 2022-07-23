import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif', 'serif':'Times'})

#h2o_ep1 = np.loadtxt('../data/ascii/epoch1/WHya_H2O.ascii')
h2o_ep1 = pd.read_csv('../data/ascii/epoch1/WHya_H2O_abs.ascii', names=['velocity','flux_density'])
#h2o_ep1_abs = pd.read_csv('../data/ascii/epoch1/WHya_H2O_abs.ascii', names=['velocity','flux_density'])
print(h2o_ep1.head())
h2o_ep2 = pd.read_csv('../data/ascii/epoch2/WHya_ep2_H2Oabs_ep2.ascii',
                        names=['velocity','flux_density'])
#h2o_ep2_abs = pd.read_csv('../data/ascii/epoch2/WHya_ep2_H2Oabs_ep2.ascii',
#                        names=['velocity','flux_density'])

fig, ax = plt.subplots(figsize=(6,6))

velocities = np.linspace(-5, 85, num=90)

#print("The velocities I'm using,", velocities)
ax.plot(h2o_ep1['velocity'], h2o_ep1['flux_density'], c='cornflowerblue', label='Nov 2015')
#ax.hist(h2o_ep1['flux_density'], histtype='step', label='Nov 2015')
ax.plot(h2o_ep2['velocity'], h2o_ep2['flux_density'], c='tab:pink', label='Nov 2017')
#ax.hist(h2o_ep2['flux_density'], histtype='step', label='Nov 2017')
# range=(min(h2o_ep1['velocity']), max(h2o_ep1['velocity'])),
ax.minorticks_on()

ax.tick_params(axis='both', which='major', direction='in', length=8, labelsize=14)
ax.tick_params(axis='both', which='minor', direction='in', length=4, labelsize=14)

ax.set_xlabel('Velocity (km/s)', fontsize=14)
ax.set_ylabel('Flux Density', fontsize=14)
plt.title('H_2O 331.1273 GHz', fontsize=14)
plt.axhline(y=0, xmin=0, xmax=1, ls='--', color='gray')
plt.axvline(x=39.5, ymin=0, ymax=1, ls='-', c='cyan')
plt.legend(fontsize=13)
plt.tight_layout()

plt.savefig('../figures/H2O_abs_spec_comparison.png', dpi=300)

plt.show()
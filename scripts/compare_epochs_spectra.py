import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif'})#, 'serif':'Times'

#h2o_ep1 = np.loadtxt('../data/ascii/epoch1/WHya_H2O.ascii')

def compare_spectral_lines(ep1_line_path, ep2_line_path, line_name=''):
    ep1_data = pd.read_csv(ep1_line_path, header=1)
    ep2_data = pd.read_csv(ep2_line_path, header=1)
    fig, ax = plt.subplots(figsize=(6,6))

    velocities = np.linspace(-5, 85, num=90)

    #print("The velocities I'm using,", velocities)
    
    ax.fill_between(ep1_data['Velocity'], ep1_data['Flux Density'], color='dodgerblue',
                zorder=2, alpha=0.6, step='pre', label='Nov 2015')
    ax.fill_between(ep2_data['Velocity'], ep2_data['Flux Density'],
                alpha=0.6, step='pre', color='gold', label='Nov 2017')

    #ax.step(ep1_data['Velocity'], ep1_data['Flux Density'], c='k', lw=0.25)
    ax.step(ep2_data['Velocity'], ep2_data['Flux Density'], c='r', lw=0.25, alpha=0.4)

    ax.minorticks_on()

    ax.tick_params(axis='both', which='major', direction='in', length=8, labelsize=14)
    ax.tick_params(axis='both', which='minor', direction='in', length=4, labelsize=14)

    ax.set_xlabel('Velocity (km/s)', fontsize=14)
    ax.set_ylabel('Flux Density', fontsize=14)
    plt.title(line_name, fontsize=14)
    plt.axhline(y=0, xmin=0, xmax=1, ls='--', color='gray')
    plt.axvline(x=39.5, ymin=0, ymax=1, ls='--', lw=1, c='dimgray', zorder=-1)

    plt.legend(fontsize=13)
    plt.tight_layout()

    f_name = ep1_line_path.replace('../data/ascii/epoch1/WHya_', '')[:-6]
    plt.savefig(f'../figures/{f_name}_spec_comparison.png', dpi=300)

    #plt.show()


if __name__ == '__main__':
    ep1_path = '../data/ascii/epoch1/'
    ep2_path = '../data/ascii/epoch2/'

    ep1_files = ['WHya_H2O.ascii', 'WHya_29SiO_v5.ascii', 'WHya_CO_v-1.ascii',
             'WHya_SiO_v2_larger_mask.ascii', 'WHya_SiO_v7.ascii', 
             # the low EP lines, I have commented out
             #'WHya_SiO_v-1.ascii', 'WHya_H13CN.ascii', 'WHya_SO_3sigma0.ascii'
             ]

    ep2_files = ['WHya_ep2_H2O.ascii', 'WHya_ep2_29SiO_v5.ascii', 'WHya_ep2_12CO_v-1.ascii',
            'WHya_ep2_SiO_v-2_larger_mask.ascii', 'WHya_ep2_SiO_v7.ascii',
            #'WHya_ep2_SiO_v-1.ascii', 'WHya_ep2_H13CN.ascii', 'WHya_ep2_SO_3sigma0.ascii'
            ]

    line_names = ['H$_2$O 331.1273 GHz', '$^{29}$SiO $v=5$', '$^{12}$CO $v=1$',
                'SiO $v=2$', 'SiO $v=7$',] #'SiO $v=1$', 'H$^{13}$CN', 'SO $3\Sigma=0$']

    for ep1_file, ep2_file, line_name in zip(ep1_files, ep2_files, line_names):
        compare_spectral_lines(ep1_path + ep1_file, ep2_path + ep2_file, line_name)
        # and the absorption regions
        ep1_abs = ep1_file[:-6] + '_abs.ascii'
        ep2_abs = ep2_file[:-6] + '_abs.ascii'
        compare_spectral_lines(ep1_path + ep1_abs, ep2_path + ep2_abs, line_name + ' absorption')
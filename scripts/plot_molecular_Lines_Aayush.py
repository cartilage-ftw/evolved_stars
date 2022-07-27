import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


"""
Author: Theo, Aayush
"""

plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif'})


def read_pf_table():
    """
    There's a table containing partition functions of molecules/ions for a range of temperatures.
    """

    path_to_table = '../data/part_func/table6.dat'
    df = pd.read_table(path_to_table, header=2, delim_whitespace=True, index_col=False)
    # a painful process of transposing and not screwing up the column names
    data = df.T
    temp_col = data.index.to_numpy().tolist()

    data.columns = data.iloc[0].copy()
    data = data.reset_index(drop=True).drop(0)

    data['T[K]'] = temp_col[1:]
    data.reset_index(drop=True, inplace=True)

    pf_data = data.apply(pd.to_numeric, errors='ignore')
    return pf_data

pf_table = read_pf_table()

# It turns out this looks like something easy to interpolate between
#plt.plot(pf_table['T[K]'], pf_table['H2'], color='gray', marker='o', mfc='magenta')
#plt.show()

def get_part_func(species, temp):
    """
    Get the partition function for an species, for an arbitrary temperature between the grids
    """
    return np.interp(temp, pf_table['T[K]'], pf_table[species])


def DopplerShift(vel,freq):

    c = 299792.458

    return (1 + vel/c)*freq


def Planck(T,f):

    """Calculates the Planck function in W m-2 sr-1 Hz-1.
        
        Args:
        T (float): black-body temperature.
        f (float or numpy array): frequency or frequencies where to evaluate the function.
        
        Returns:
        float: Values of the Planck function at the corresponding frequencies.
        """


    h = 6.62607004E-34 #J/Hz
    k = 1.380649E-23 #J/K
    c = 299792458.0

    return (2.0 * h * f**3/c/c)/(1.0-np.exp(-h*f/(k*T)))

class Line:

    """Line class: Each line to be considered for a given molecule must be defined using this class and included when the molecule is defined.
        """


    def __init__(self, ID, freq, A_ul, g_up, g_low, E_up, deltaVelAbs,deltaVelEm, dist):

        """Function for initializing the line class. The function does the following:
            - defines the frequency array over which the line will be calculated;
            - calculates the normalized line shape for the absorption and emission profiles assuming gaussian lines;
            - initializes the calculation of the optical depth, and the line emission and absorption using the variables that are available at this stage;
            
            
            Args:
            freq (float): frequency of the transition in Hertz.
            A_ul (float): Einstein A coefficient of the transition in s-1.
            g_up (float): degeneracy of the upper level of the transition.
            g_low (float): degeneracy of the lower level of the transition.
            E_up (float): energy of the upper level of the transition in Kelvin.
            deltaVelAbs (float): standard deviation of the gaussian line shape for the absorption component.
            deltaVelEm (float): standard deviation of the gaussian line shape for the emission component.
            dist (float): distance to the star
            """


        h = 6.62607004E-34 #J/Hz
        k = 1.380649E-23 #J/K

        self.ID = ID
        self.freq = freq #in Hertz
        self.E_up = E_up # in Kelvin
        self.A_ul = A_ul
        self.g_up = g_up
        self.g_low = g_low
        self.E_low = E_up-freq*h/k

        deltaNuAbs = self.freq*deltaVelAbs/(c)
        deltaNuEm = self.freq*deltaVelEm/(c)


        self.freqRange = np.arange(-100,100)/50.0
        self.freqRange = self.freqRange*6.0*deltaNuAbs + freq
        
        self.lineCenterAbs = 1.0/(deltaNuAbs * np.sqrt(2.0*np.pi))
        self.lineCenterEm = 1.0/(deltaNuEm * np.sqrt(2.0*np.pi))
        self.lineShapeEm = np.exp( - np.power(self.freqRange - freq,2.0) / (2.0*np.power(deltaNuEm,2.0) ) ) * self.lineCenterEm
        self.lineShapeAbs = np.exp( - np.power(self.freqRange - freq,2.0) / (2.0*np.power(deltaNuAbs,2.0) ) ) * self.lineCenterAbs
        self.tauNu = self.lineShapeAbs
        self.tauCenter = self.lineCenterAbs
        
        self.emission = self.lineShapeEm*self.A_ul*h*self.freq / (4.0*np.pi*dist*dist)
        self.absorption = self.lineShapeAbs
    
        self.nUp = 0.0
        self.nLow = 0.0

#        self.strength = self.A_ul*h*self.freq / (4.0*np.pi*dist*dist)
        
        self.d2 = dist*dist

        self.stellarIntensity = 0.0


class Molecule:

    """Molecule class: The molecule to be considered is defined by this class. The individual lines to be calculated are given as input.
        """

    def __init__(self, ID, partFunc, lines, Texc, Nmol, Rout, Rstar, Tstar, color):
        
        """Function for initializing the molecule class. The function does the following:
            - calculates the volume of gas that will be seen against the star;
            - calculates the average column density of molecules in the gas seen against the star;
            - calculates the volume of gas that will be seen around the star;
            - calculates the stellar emission;
            - calculates the populations for levels involved in the transitions considered using the excitation temperature provided;
            - updates the values of the emission and absortion for each line by adding the newly calculated quantities;
            - subtracts the stellar emission from the emission seen towards the star.
            
            Args:
            freq (float): frequency of the transition in Hertz.
            A_ul (float): Einstein A coefficient of the transition in s-1.
            g_up (float): degeneracy of the upper level of the transition.
            g_low (float): degeneracy of the lower level of the transition.
            E_up (float): energy of the upper level of the transition in Kelvin.
            deltaVelAbs (float): standard deviation of the gaussian line shape for the absorption component.
            deltaVelEm (float): standard deviation of the gaussian line shape for the emission component.
            dist (float): distance to the star
            """

        
        h = 6.62607004E-34
        k = 1.380649E-23 #J/K
        c = 299792458.0
        
        self.ID = ID
        self.Z = partFunc
        self.lines = lines
        self.Texc = Texc
        self.Nmol = Nmol    # In molecules per cubic meter
        self.color = color
        
        Rout_m = Rout*1.496E11      #Conversion from AU to m
        Rstar_m = Rstar*1.496E11    #Conversion from AU to m
        
        z_m      = np.sqrt(Rout_m * Rout_m - Rstar_m * Rstar_m)
        height_m = Rout_m - z_m
        self.volAbs = (np.pi*height_m*height_m/3.0)*(3.0*Rout_m - height_m)
        self.volAbs += (np.pi*Rstar_m*Rstar_m*z_m) - (2.0*np.pi*Rstar_m**3.0/3.0)
        self.averColDens = self.Nmol*self.volAbs/(np.pi*Rstar_m*Rstar_m) #Average column density in molecules per m^2 of the gas in front of the star

        self.volEm = (4.0*np.pi*Rout_m**3.0/3.0) - (4.0*np.pi*Rstar_m**3.0/3.0) - (2.0*self.volAbs)
                
        for line in self.lines:
            line.Istar = 2.0*line.freq*line.freq*k*Tstar/c/c #W m-2 sr-1 Hz-1
            
            line.nUp = line.g_up*np.exp(-line.E_up/self.Texc)/self.Z
            line.nLow = line.g_low*np.exp(-line.E_low/self.Texc)/self.Z
            line.emission = line.emission * self.Nmol*line.nUp * self.volEm * 1E26
            line.tauCenter = c*c/(8.0*np.pi*np.power(line.freq,2.0))*(line.g_up/line.g_low)*self.averColDens*line.nLow*line.A_ul*(1.0-np.exp(-h*line.freq/(k*Texc)))*line.tauCenter
            line.tauNu     = c*c/(8.0*np.pi*np.power(line.freq,2.0))*(line.g_up/line.g_low)*self.averColDens*line.nLow*line.A_ul*(1.0-np.exp(-h*line.freq/(k*Texc)))*line.tauNu

            line.absorption = Planck(Texc,line.freq)*(1.0 - np.exp(-line.tauNu)) + line.Istar*np.exp(-line.tauNu)
            line.absorption = line.absorption - line.Istar
            line.absorption = line.absorption * 1E26 * np.pi * Rstar_m * Rstar_m / line.d2


    def PrintValues(self):
    
        print("\nPrinting values of molecule: ",self.ID)
        print("Lines included: ",[line.ID for line in self.lines])
        print("Partition Function: ", self.Z)
        print("Excitation Temperature: ", self.Texc, " K")
        print("Density:", self.Nmol, " m-3")    # In molecules per cubic meter
        return ""


c = 299792.458
h = 6.62607004E-34
dist = 102 #in parsec
dist = dist*3.086E16 #conversion to meters
deltaVelAbs = 5.0 # in km/s
deltaVelEm = 5.3 # in km/s
Rstar = 2.0 #in au
vstar = 39.9
absorptionShift = -11 #km/s, Theo set it to 9 initially

Rstar_ep1_mas = 26.9 # in milliarcsecs, Vlemmings+17
Rstar_ep2_mas = 24.2
Tstar_ep1 = 2495
Tstar_ep2 = 2680

#convert dist to AU
Rstar_ep1 = Rstar_ep1_mas*102/1E3
Rstar_ep2 = Rstar_ep2_mas*102/1E3

print(Rstar_ep1)

SiOv2_8_7  = Line("SiOv2_8_7",342.50460700E9,0.00216616,17.0,15.0,3595.12278, deltaVelAbs,deltaVelEm, dist)
SiOv7_8_7  = Line("SiOv7_8_7",330.4775327E9,0.002080,17.0,15.0,12082.533, deltaVelAbs,deltaVelEm, dist)


#SiOv2_300K = Molecule("SiOv2",2240.7865,[SiOv2_8_7],1500.0,1.75E10,2.0,Rstar,2500.0,"magenta")
#SiOv2_1500K = Molecule("SiOv2",2098.81,[SiOv2_8_7],1500.0,1.75E10,2.0,Rstar,2500.0,"crimson")
#SiOv2_ep1 = Molecule("SiOv2",get_part_func('SiO', 770),[SiOv2_8_7],770,1.25E9,2.49*Rstar_ep1,Rstar_ep1,
#            Tstar_ep1,"dimgray")
#SiOv2_ep2_sameRout = Molecule("SiOv2", get_part_func('SiO', 1030), [SiOv2_8_7], 1030, 9E8, 2.49*Rstar_ep1, Rstar_ep2,
#                Tstar_ep2, 'red')
SiOv2_ep2_diffRout= Molecule("SiOv2", get_part_func('SiO', 1210), [SiOv2_8_7], 1210, 1E9, 2.48*Rstar_ep2, Rstar_ep2,
                Tstar_ep2, 'crimson')
#COv1_ep1 = Molecule("CO", get_part_func('CO', 770), )
#SiO_1000K = Molecule("SiO",1162.9451,[SiOv2_8_7,SiOv7_8_7],1000.0,3.0E10,2.0,Rstar,2500.0,"magenta")
#SiO_2000K = Molecule("SiO",3318.628,[SiOv2_8_7,SiOv7_8_7],2000.0,3.0E10,2.0,Rstar,2500.0,"magenta")
#SiO_3000K = Molecule("SiO",4481.5731,[SiOv2_8_7,SiOv7_8_7],3000.0,3.0E10,2.0,Rstar,2500.0,"magenta")

#molecules = [SiO_1000K]
#molecules = [SiO_2000K]
#molecules = [SiO_3000K]
molecules = [SiOv2_ep2_diffRout]

print('The partition function for CO predicted here is', get_part_func('CO', 2000))

ep1_line_path = '../data/ascii/epoch1/WHya_SiO_v2_larger_mask.ascii'
ep2_line_path = '../data/ascii/epoch2/WHya_ep2_SiO_v-2_larger_mask.ascii'

if __name__ == '__main__':
    # load ascii
    ep1_data = pd.read_csv(ep1_line_path, header=1)
    ep2_data = pd.read_csv(ep2_line_path, header=1)
    fig, ax = plt.subplots(figsize=(6,6.67))

    ep1_data['Frequency'] = DopplerShift(vstar, ep1_data['Frequency'])
    ep2_data['Frequency'] = DopplerShift(vstar, ep2_data['Frequency'])
    velocities = np.linspace(-5, 85, num=90)

    #print("The velocities I'm using,", velocities)


    ax.fill_between(ep1_data['Frequency']/1E9, ep1_data['Flux Density'], color='dodgerblue',
                zorder=2, alpha=0.6, step='pre', label='Nov 2015')
    ax.fill_between(ep2_data['Frequency']/1E9, ep2_data['Flux Density'],
                alpha=0.6, step='pre', color='gold', label='Nov 2017')

    #ax.step(ep1_data[freq_axis], ep1_data['Flux Density'], c='k', lw=0.25)
    ax.step(ep2_data['Frequency'], ep2_data['Flux Density'], c='r', lw=0.25, alpha=0.4, zorder=1)

    ax.minorticks_on()

    ax.tick_params(axis='both', which='major', direction='in', length=8, labelsize=14)
    ax.tick_params(axis='both', which='minor', direction='in', length=4, labelsize=14)

    ax.set_xlabel('Frequency (Hz)', fontsize=14)
    ax.set_ylabel('Flux Density', fontsize=14)
    plt.title('SiO $v=2$', fontsize=14)
    plt.axhline(y=0, xmin=0, xmax=1, ls='--', color='gray')
    #if freq_axis == 'Velocity':
    #    plt.axvline(x=39.5, ymin=0, ymax=1, ls='--', lw=1, c='dimgray', zorder=-1)
    #else:
    #   print('To mark the hrv=39.5 km/s line, convert that into frequency first.')

    #print(molecules[0].PrintValues())

    for m in molecules:
        m.PrintValues()
        for line in m.lines:
            #if ((line.freq/1E9-spec[i][0][0]) * (line.freq/1E9-spec[i][0][-1])) < 0.0:
            emissionPlot = np.interp((line.freqRange+line.freq*absorptionShift/299792.458), line.freqRange, line.emission)
            line_label='T$_{exc}=$' + str(round(m.Texc)) + ', $N_{mol}=$'+str(m.Nmol/1E9) + 'E9'
            ax.plot((line.freqRange+line.freq*absorptionShift/299792.458)/1E9, (line.absorption+emissionPlot)*1E3, color=m.color,
                        label=line_label)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=2, fontsize=13)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("Flux Dens. [mJy]")
    ax.set_xlim(min(ep1_data['Frequency'])/1E9, max(ep2_data['Frequency'])/1E9)
    plt.tight_layout()
    plt.savefig(f'../figures/match_with_synth_SiOv2.png', dpi=300)
    plt.show()
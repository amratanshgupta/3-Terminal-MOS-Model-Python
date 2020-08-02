import numpy as np
import matplotlib.pyplot as plt
import math as m
import random
from Mathematical_functions import diff
from Class_Surface_Potential_To_VCB_3Terminal import Surface_Potential_To_VCB
import matplotlib.colors as colors
colors_list = list(colors._colors_full_map.values())


class Total_charge_VCB(Surface_Potential_To_VCB):
    """Returns Inversion charge, Depletion Charge, Total Charge in 2-Terminal MOS"""
    vt = 25.8e-3
    q = 1.6e-19
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)
    ni = 1.5e10
    # Cox = epox / tox

    def __init__(self, tox, gamma, Vcb_max, Vgb=0, Na=1e17):
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)
        Surface_Potential_To_VCB.__init__(self, tox, gamma, Vcb_max, Vgb, Na=1e17)
        self.phi_f = vt * m.log(self.Na / Surface_Potential_To_VCB.ni)
        self.Na = Na
        self.Vcb_max = Vcb_max
        self.Vgb = Vgb
        self.Qtotal = []
        self.Qinv = []
        self.Qdep = []
        self.VCB_QTotal = {}
        self.VCB_Qdep = {}
        self.VCB_Qinv = {}
        self.VCB_Qinv_list = []

        self.Surface_Potential, self.phi_vcb_map, self.phi_regions = Surface_Potential_To_VCB.Psi_to_VCB(self)
        # print(self.phi_vcb_map)

    def Total_Charge(self):
        """Returns total semiconductor charge as a function of VCB"""
        self.Qtotal = []
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)
        for (k, t) in self.phi_vcb_map.items():
            first = m.sqrt(abs((vt * m.exp(-k / vt)) + (k - vt) + (m.exp((-(2 * self.phi_f + t) / vt)) * ((vt * m.exp(k / vt)) - k - vt))))
            if k >= 0:
                Q_total = -1 * first * (m.sqrt(2 * q * eps * self.Na))
                self.Qtotal.append(Q_total)
                mapped_VCB = self.phi_vcb_map[k]
                self.VCB_QTotal[mapped_VCB] = Q_total
                # self.VCB_QTotal[self.phi_vcb_map[k]] = Q_total

            elif k < 0:
                Q_total = first * (m.sqrt(2 * q * eps * self.Na))
                self.Qtotal.append(Q_total)
                mapped_VCB = self.phi_vcb_map[k]
                self.VCB_QTotal[mapped_VCB] = Q_total

        return(self.VCB_QTotal, self.Qtotal)

    def Qdepletion(self):
        """Returns Depletion charge as a function of VCB """
        self.Qdep = []
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)
        for k in self.phi_vcb_map.keys():
            first = m.sqrt(abs((vt * m.exp(-k / vt)) + (k - vt)))
            if k >= 0:
                Q_depletion = -1 * first * (m.sqrt(2 * q * eps * self.Na))
                self.Qdep.append(Q_depletion)
                mapped_VCB = self.phi_vcb_map[k]
                self.VCB_Qdep[mapped_VCB] = Q_depletion

            elif k < 0:
                Q_depletion = first * (m.sqrt(2 * q * eps * self.Na))
                self.Qdep.append(Q_depletion)
                mapped_VCB = self.phi_vcb_map[k]
                self.VCB_Qdep[mapped_VCB] = Q_depletion
        # print(f'Qdep {len(self.Qdep)}', f'{len(self.Surface_Potential_Qdep.keys())}')
        return(self.VCB_Qdep, self.Qdep)

    def Qinversion(self):
        """Returns Inversion charge as difference of Total charge and depletion charge as function of varying VCB """
        self.Qinv = []
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)

        (Qdep_dic, Qdep_list) = Total_charge_VCB.Qdepletion(self)
        (Qtotal_dic, Qtotal_list) = Total_charge_VCB.Total_Charge(self)

        for (i, j) in zip(Qtotal_list, Qdep_list):
            k = (i * 1e7) - (j * 1e7)
            Q_inver = k * 1e-7
            self.Qinv.append(Q_inver)

        for var in Qdep_dic.keys():
            self.VCB_Qinv_list.append(var)

        # print(f'Qinv: {len(self.Qinv)}', f'{len(self.Surface_Potential_Qinv_list)}')
        return(self.VCB_Qinv_list, self.Qinv)

    def Plot_Q_vs_VCB(self):
        vt = 25.8e-3
        common_plotting_factor = 2.2

        self.Qinv = np.array(self.Qinv)
        self.Qdep = np.array(self.Qdep)
        self.Qtotal = np.array(self.Qtotal)

        plt.plot(self.VCB_Qinv_list, self.Qinv * 1e6, 'r', label=r'$Q_{inv}$')
        plt.plot(self.VCB_Qdep.keys(), self.Qdep * 1e6, 'g', label=r'$Q_{dep}$')
        plt.plot(self.VCB_QTotal.keys(), self.Qtotal * 1e6, 'b', linewidth=2.5, linestyle=':', label=r'$Q_{sc}$')

        plt.title(r'Total Semiconductor Charge($Q_{SC}$) vs  Source-Body bias for 3-Terminal MOS')

        # Axis Label
        # plt.ylabel(r'Value [x 10^{-7}]')
        plt.ylabel(r'$Q_{Total}$  ($\mu C$)')
        plt.xlabel(r'$V_{CB}$ (V)')

        plt.legend()
        ## Putting Applied Gate-Body bias on plot
        plt.text(self.Vcb_max / 1.2, min(self.Qdep * 1e6) * 0.5, r'$V_{GB}$ =' + '{:.2f} V'.format(self.Vgb))

        plt.axvline(x=0, ymin=0.0, ymax=1.0, color='k')
        plt.axhline(y=0, xmin=0.0, xmax=1.0, color='k')

        # VCB Regions - Verticle lines : Weak/ Moderate/Strong Inversion
        plt.axvline(x=self.Vcb_regions[0], color='#968F8F', linestyle='--', alpha=0.75)
        plt.axvline(x=self.Vcb_regions[1], color='#968F8F', linestyle='--', alpha=0.75)
        plt.axvline(x=self.Vcb_regions[2], color='#968F8F', linestyle='--', alpha=0.75)

        #### Regions Text #####
        plt.annotate("Depletion", xy=((self.Vcb_regions[0] * 1.15, 0.85 * (min(self.Qdep * 1e6)))))
        plt.annotate("Weak\nInversion", xy=((self.Vcb_regions[1] * 1.05, 0.85 * (min(self.Qdep * 1e6)))))
        plt.annotate("Moderate\nInversion", xy=((self.Vcb_regions[2], 0.85 * (min(self.Qdep * 1e6)))))
        plt.annotate("Strong\nInversion", xy=((self.Vcb_regions[2] * 0.76, 0.85 * (min(self.Qdep * 1e6)))))

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    Na = 1e17
    q = 1.6e-19
    tox = 40e-7
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)
    Cox = epox / tox
    Vgb = 6
    Vcb_max = 3
    gamma = ((m.sqrt(2 * q * eps * Na)) / Cox)

    x = Total_charge_VCB(tox, gamma, Vcb_max, Vgb, Na)
    y = x.Total_Charge()
    z = x.Qinversion()
    v = x.Qdepletion()

    b = x.Plot_Q_vs_VCB()

import numpy as np
import matplotlib.pyplot as plt
import math as m
import random
from Mathematical_functions import diff
from Class_Surface_Potential_To_VGB_3Terminal import Surface_Potential_To_VGB
import matplotlib.colors as colors
colors_list = list(colors._colors_full_map.values())
# plt.style.use('seaborn-bright')
# print(plt.style.available)


class Total_charge_VGB_Surface_Potential(Surface_Potential_To_VGB):
    """Returns Inversion charge, Depletion Charge, Total Charge in 2-Terminal MOS"""
    vt = 25.8e-3
    q = 1.6e-19
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)
    ni = 1.5e10
    # Cox = epox / tox

    def __init__(self, tox, gamma, Vgb_max, Vcb=0, Na=1e17):
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)
        Surface_Potential_To_VGB.__init__(self, tox, gamma, Vgb_max, Vcb, Na=1e17)
        self.phi_f = vt * m.log(self.Na / Surface_Potential_To_VGB.ni)
        self.Na = Na
        self.Vcb = Vcb
        self.Qtotal = []
        self.Qinv = []
        self.Qdep = []
        self.Surface_Potential_QTotal = {}
        self.Surface_Potential_Qdep = {}
        self.Surface_Potential_Qinv = {}
        self.Surface_Potential_Qinv_list = []
        self.VGB_QToTal = {}
        self.VGB_Qdep = {}
        self.VGB_Qinv = {}

        self.Surface_Potential, self.phi_vgb_map, self.phi_regions = Surface_Potential_To_VGB.Psi_to_VGB(self)
        # print(self.phi_vgb_map)

    def Total_Charge(self):
        self.Qtotal = []
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)
        for k in self.phi_vgb_map.keys():
            first = m.sqrt(abs((vt * m.exp(-k / vt)) + (k - vt) + (m.exp((-(2 * self.phi_f + self.Vcb) / vt)) * ((vt * m.exp(k / vt)) - k - vt))))
            if k >= 0:
                Q_total = -1 * first * (m.sqrt(2 * q * eps * self.Na))
                self.Qtotal.append(Q_total)
                self.Surface_Potential_QTotal[k] = Q_total
                # self.VGB_QTotal[self.phi_vgb_map[k]] = Q_total

            elif k < 0:
                Q_total = first * (m.sqrt(2 * q * eps * self.Na))
                self.Qtotal.append(Q_total)
                self.Surface_Potential_QTotal[k] = Q_total

        return(self.Surface_Potential_QTotal, self.Qtotal)

    def Qdepletion(self):
        self.Qdep = []
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)
        for k in self.phi_vgb_map.keys():
            first = m.sqrt(abs((vt * m.exp(-k / vt)) + (k - vt)))
            if k >= 0:
                Q_depletion = -1 * first * (m.sqrt(2 * q * eps * self.Na))
                self.Qdep.append(Q_depletion)
                self.Surface_Potential_Qdep[k] = Q_depletion
                # self.VGB_QTotal[self.phi_vgb_map[k]] = Q_total

            elif k < 0:
                Q_depletion = first * (m.sqrt(2 * q * eps * self.Na))
                self.Qdep.append(Q_depletion)
                self.Surface_Potential_Qdep[k] = Q_depletion

        # print(f'Qdep {len(self.Qdep)}', f'{len(self.Surface_Potential_Qdep.keys())}')
        return(self.Surface_Potential_Qdep, self.Qdep)

    def Qinversion(self):
        self.Qinv = []
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)
        # for k in self.phi_vgb_map.keys():
        #     vt = 25.8e-3
        #     first = abs(m.sqrt(abs(((m.exp((-2 * self.phi_f) / vt) * ((vt * m.exp(k / vt)) - k - vt)) - vt))) - m.sqrt(abs(k - vt)))
        #     Q_inv = -1 * first * (m.sqrt(2 * q * eps * Na))
        #     self.Qinv.append(Q_inv)
        #     self.Surface_Potential_Qinv[k] = Q_inv

        # for i, j in zip(self.Surface_Potential_Qinv.keys(), self.Qinv):
        #     print(i, j, end='\n')
        # return(self.Surface_Potential_Qdep, self.Qinv)
        (Qdep_dic, Qdep_list) = Total_charge_VGB_Surface_Potential.Qdepletion(self)
        (Qtotal_dic, Qtotal_list) = Total_charge_VGB_Surface_Potential.Total_Charge(self)

        for (i, j) in zip(Qtotal_list, Qdep_list):
            k = (i * 1e7) - (j * 1e7)
            Q_inver = k * 1e-7
            self.Qinv.append(Q_inver)

        for var in Qdep_dic.keys():
            self.Surface_Potential_Qinv_list.append(var)

        # print(f'Qinv: {len(self.Qinv)}', f'{len(self.Surface_Potential_Qinv_list)}')
        return(self.Surface_Potential_Qinv_list, self.Qinv)

    def Plot_Q_vs_Surface_Potential(self):
        # Total_charge_VGB_Surface_Potential.Qdepletion(self)
        # Total_charge_VGB_Surface_Potential.Qinversion(self)

        # Plots for
        plt.plot(self.Surface_Potential_Qinv_list, self.Qinv, 'r', label=r'$Q_{inv}$')
        plt.plot(self.Surface_Potential_Qdep.keys(), self.Qdep, 'g', label=r'$Q_{dep}$')
        plt.plot(self.phi_vgb_map.keys(), self.Qtotal, 'k', linestyle='--', label=r'$Q_{sc}$')
        plt.title(r'Total Semiconductor Charge($Q_{SC}$) vs  Surface Potential for 3-Terminal MOS')

        plt.ylabel(r'$Q_{Total}$ (C)')
        plt.xlabel('Surface Potential (V)')

        plt.legend()

        plt.axvline(x=0, ymin=0.0, ymax=1.0, color='k')
        plt.axhline(y=0, xmin=0.0, xmax=1.0, color='k')
        # plt.plot(self.Surface_Potential_Qdep.keys(), self.Surface_Potential_Qdep.items(),'r')
        # plt.show()
        # plt.plot(self.Surface_Potential_Qinv.keys(), self.Surface_Potential_Qinv.items(),'g')
        plt.tight_layout()
        plt.show()

        # def Plot_Q_vs_VGB(self):
        #     plt.plot(self.phi_vgb_map.items(), self.Qtotal)
        #     plt.title('Total Semiconductor Charge vs  Gate-Body Voltage for 2 Terminal MOS')
        #     plt.plot(self.Surface_Potential_Qdep.keys(), self.Surface_Potential_Qdep.items()
        #     plt.plot(self.phi_vgb_map.items(), self.Qinv)

        #     plt.ylabel(r'$Q_{Total}$ (C)')
        #     plt.xlabel(r'$V_{GB}$ (V)')

        #     #######  Lines and Labels ######
        #     ##  X & Y Axis ###
        #     plt.axhline(y=0, xmin=0.0, xmax=1.0, color='k')
        #     plt.axvline(x=0, ymin=0.0, ymax=1.0, color='k')
        #     plt.show()


if __name__ == "__main__":
    Na = 1e17
    q = 1.6e-19
    tox = 40e-7
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)
    Cox = epox / tox
    Vcb = 1.0
    Vgb_max = 9
    gamma = ((m.sqrt(2 * q * eps * Na)) / Cox)

    x = Total_charge_VGB_Surface_Potential(tox, gamma, Vgb_max, Vcb, Na)
    y = x.Total_Charge()
    z = x.Qinversion()
    v = x.Qdepletion()

    b = x.Plot_Q_vs_Surface_Potential()
    # a = x.Plot_Q_vs_VGB()

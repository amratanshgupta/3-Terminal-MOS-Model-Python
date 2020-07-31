import numpy as np
import matplotlib.pyplot as plt
import math as m
import random
from Mathematical_functions import diff
from Class_Surface_Potential_To_VGB_3Terminal import Surface_Potential_To_VGB
from Class_Surface_Potential_VGB_To_Total_Charge_3Terminal import Total_charge_VGB_Surface_Potential
import matplotlib.colors as colors
colors_list = list(colors._colors_full_map.values())
# plt.style.use('seaborn-paper')
# print(plt.style.available)


class Capacitance_To_VGB(Total_charge_VGB_Surface_Potential):
    vt = 25.8e-3
    q = 1.6e-19
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)

    def __init__(self, tox, gamma, Vgb_max, Vcb=0, Na=1e17):
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)

        Total_charge_VGB_Surface_Potential.__init__(self, tox, gamma, Vgb_max, Vcb, Na=1e17)
        self.phi_f = vt * m.log(self.Na / Surface_Potential_To_VGB.ni)
        self.Na = Na
        self.tox = tox
        self.Vcb = Vcb
        self.Ctotal = []
        self.Cinv = []
        self.Cdep = []
        self.Cox = epox / self.tox
        self.Csemiconductor = []
        self.C_SC_total = []
        self.VGB_Ctotal = {}
        self.VGB_Cdep = {}
        self.VGB_Cinv = {}
        self.VGB_Cinv_list = []

    def Dep_Capacitance(self):
        self.Cdep = []
        self.VGB_Cdep = {}
        # print(self.phi_vgb_map.keys())
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)

        for k in self.phi_vgb_map.keys():
            if k > 0.075:
                # print(k)
                first = 1 / (m.sqrt(abs(k + vt * (m.exp((k - (2 * self.phi_f + self.Vcb)) / vt)))))
                Cdep = first * (m.sqrt(2 * q * eps * self.Na)) / 2
                self.Cdep.append(Cdep)
                local_VGB = self.phi_vgb_map[k]
                self.VGB_Cdep[local_VGB] = Cdep
        return(self.VGB_Cdep, self.Cdep)

    def Inv_Capacitance(self):
        self.Cinv = []
        self.Surface_Potential_Cinv = {}
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)

        for k in self.phi_vgb_map.keys():
            if k > 0.075:
                # print(k)
                first = 1 / (m.sqrt(abs(k + vt * (m.exp((k - (2 * self.phi_f + self.Vcb)) / vt)))))
                # print(first)
                C_inv = first * (m.sqrt(2 * q * eps * self.Na)) * m.exp((k - (2 * self.phi_f + self.Vcb)) / vt) / 2
                self.Cinv.append(C_inv)
                local_VGB = self.phi_vgb_map[k]
                self.VGB_Cinv[local_VGB] = C_inv
        return(self.VGB_Cinv, self.Cinv)

    def Semiconductor_Capacitance(self):
        VGB_Cinv, Cinv = Capacitance_To_VGB.Inv_Capacitance(self)
        VGB_Cdep, Cdep = Capacitance_To_VGB.Dep_Capacitance(self)
        self.Csemiconductor = []
        self.VGB_Csemiconductor = {}
        for i, j in zip(Cinv, Cdep):
            k = i + j
            self.Csemiconductor.append(k)

        for a, b in zip(VGB_Cdep.keys(), self.Csemiconductor):
            self.VGB_Csemiconductor[a] = b

        return(self.VGB_Csemiconductor, self.Csemiconductor)

    def Total_Capacitance(self):
        # print(self.Cox)
        VGB_Csemiconductor, Csemiconductor = Capacitance_To_VGB.Semiconductor_Capacitance(self)
        self.VGB_CTotal = {}
        self.Ctotal = []

        for i in Csemiconductor:
            Ctotal = (self.Cox * i) / (self.Cox + i)
            self.Ctotal.append(Ctotal)

        for a, b in zip(VGB_Csemiconductor.keys(), self.Ctotal):
            self.VGB_Ctotal[a] = b

        return(self.VGB_Ctotal, self.Ctotal)

    def Plot_C_vs_VGB(self):
        plt.semilogy(self.VGB_Cinv.keys(), self.Cinv, 'r', linestyle='--', label=r'$C_{inv}$')
        plt.semilogy(self.VGB_Cdep.keys(), self.Cdep, 'g', linestyle='--', label=r'$C_{dep}$')
        plt.semilogy(self.VGB_Csemiconductor.keys(), self.Csemiconductor, 'k', linestyle=':', label=r'$C_{semiconductor}$')
        plt.semilogy(self.VGB_Ctotal.keys(), self.Ctotal, 'b', linestyle='-.', label=r'$C_{Total}$')
        plt.axhline(y=self.Cox, color='#E111B7', linestyle=':', label=r'$C_{ox}$', alpha=0.85)
        plt.title(r'Capacitance($C_{}$) vs  $V_{GB}$ for 2-Terminal MOS')
        plt.ylabel('Capacitance (F) - (log Scale)')
        plt.xlabel(r'$V_{GB}$ (V)')
        plt.legend()
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    Na = 1e17
    q = 1.6e-19
    tox = 4e-7
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)
    Cox = epox / tox
    print(Cox)
    Vgb_max = 3.3
    Vcb = 1.0
    gamma = ((m.sqrt(2 * q * eps * Na)) / Cox)

    x = Capacitance_To_VGB(tox, gamma, Vgb_max, Vcb, Na)
    y = x.Dep_Capacitance()
    z = x.Inv_Capacitance()
    a = x.Semiconductor_Capacitance()
    b = x.Total_Capacitance()
    c = x.Plot_C_vs_VGB()

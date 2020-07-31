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


class Capacitance_To_Surface_Potential(Total_charge_VGB_Surface_Potential):
    vt = 25.8e-3
    ni = 1.5e10
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)
    ni = 1.5e10
    q = 1.6e-19

    def __init__(self, tox, gamma, Vgb_max, Vcb=0, Na=1e17):
        vt = 25.8e-3
        epox = 8.854e-14 * (3.9)
        eps = 8.854e-14 * (11.8)
        q = 1.6e-19
        Total_charge_VGB_Surface_Potential.__init__(self, tox, gamma, Vgb_max, Vcb, Na=1e17)
        self.phi_f = vt * m.log(self.Na / Surface_Potential_To_VGB.ni)
        self.Ctotal = []
        self.Cinv = []
        self.Cdep = []
        self.tox = tox
        self.Vcb = Vcb
        self.Cox = epox / self.tox
        self.Csemiconductor = []
        self.C_SC_total = []
        self.Surface_Potential_Ctotal = {}
        self.Surface_Potential_Cdep = {}
        self.Surface_Potential_Cinv = {}
        self.Surface_Potential_Cinv_list = []

    def Dep_Capacitance(self):
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        self.Cdep = []
        self.Surface_Potential_Cdep = {}
        # print(self.phi_vgb_map.keys())
        for k in self.phi_vgb_map.keys():
            if k > 0.05:
                first = 1 / (m.sqrt(abs(k + vt * (m.exp((k - (2 * self.phi_f + self.Vcb)) / vt)))))
                Cdep = first * (m.sqrt(2 * q * eps * self.Na)) / 2
                self.Cdep.append(Cdep)
                self.Surface_Potential_Cdep[k] = Cdep
        return(self.Surface_Potential_Cdep, self.Cdep)

    def Inv_Capacitance(self):
        self.Cinv = []
        self.Surface_Potential_Cinv = {}
        vt = 25.8e-3
        q = 1.6e-19
        eps = 8.854e-14 * (11.8)
        for k in self.phi_vgb_map.keys():
            if k > 0.05:
                # print(k)
                first = 1 / (m.sqrt(abs(k + vt * (m.exp((k - (2 * self.phi_f + self.Vcb)) / vt)))))
                # print(first)
                C_inv = first * (m.sqrt(2 * q * eps * self.Na)) * m.exp((k - (2 * self.phi_f + self.Vcb)) / vt) / 2
                self.Cinv.append(C_inv)
                self.Surface_Potential_Cinv[k] = C_inv

        return(self.Surface_Potential_Cinv, self.Cinv)

    def Semiconductor_Capacitance(self):
        Surface_Potential_Cinv, Cinv = Capacitance_To_Surface_Potential.Inv_Capacitance(self)
        Surface_Potential_Cdep, Cdep = Capacitance_To_Surface_Potential.Dep_Capacitance(self)
        self.Csemiconductor = []
        self.Surface_Potential_Csemiconductor = {}
        for i, j in zip(Cinv, Cdep):
            k = i + j
            self.Csemiconductor.append(k)

        for a, b in zip(Surface_Potential_Cdep.keys(), self.Csemiconductor):
            self.Surface_Potential_Csemiconductor[a] = b

        return(self.Surface_Potential_Csemiconductor, self.Csemiconductor)

    def Total_Capacitance(self):
        # self.Cox = Cox
        # print(self.Cox)
        Surface_Potential_Csemiconductor, Csemiconductor = Capacitance_To_Surface_Potential.Semiconductor_Capacitance(self)
        self.Surface_Potential_CTotal = {}
        self.Ctotal = []

        for i in Csemiconductor:
            Ctotal = (self.Cox * i) / (self.Cox + i)
            self.Ctotal.append(Ctotal)

        for a, b in zip(Surface_Potential_Csemiconductor.keys(), self.Ctotal):
            self.Surface_Potential_Ctotal[a] = b

        return(self.Surface_Potential_Ctotal, self.Ctotal)

    def Plot_C_vs_Surface_Potential(self):
        plt.semilogy(self.Surface_Potential_Cinv.keys(), self.Cinv, 'r', linestyle='--', label=r'$C_{inv}$')
        plt.semilogy(self.Surface_Potential_Cdep.keys(), self.Cdep, 'g', linestyle='--', label=r'$C_{dep}$')
        plt.semilogy(self.Surface_Potential_Csemiconductor.keys(), self.Csemiconductor, 'k', linestyle=':', label=r'$C_{semiconductor}$')
        plt.semilogy(self.Surface_Potential_Ctotal.keys(), self.Ctotal, 'b', linestyle='-.', label=r'$C_{ox}$')
        plt.axhline(y=self.Cox, color='#E111B7', linestyle='-.', alpha=0.85)
        plt.title(r'Capacitance($C_{}$) vs  Surface Potential for 2-Terminal MOS')
        plt.ylabel('Capacitance (F)')
        plt.xlabel('Surface Potential (V)')
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
    Vgb_max = 12
    Vcb = 1.0
    gamma = ((m.sqrt(2 * q * eps * Na)) / Cox)

    x = Capacitance_To_Surface_Potential(tox, gamma, Vgb_max, Vcb, Na)
    y = x.Dep_Capacitance()
    z = x.Inv_Capacitance()
    a = x.Semiconductor_Capacitance()
    b = x.Total_Capacitance()
    c = x.Plot_C_vs_Surface_Potential()

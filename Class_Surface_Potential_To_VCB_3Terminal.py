import numpy as np
import matplotlib.pyplot as plt
import math as m
import random
from Mathematical_functions import diff
import matplotlib.colors as colors
colors_list = list(colors._colors_full_map.values())
# plt.style.use('seaborn-paper')
# print(plt.style.available)


class Surface_Potential_To_VCB():
    vt = 25.8e-3
    phi = 0.05
    ni = 1.5e10
    vt = 25.8e-3
    phi = 0.05
    delta = 0.001
    vt = 25.8e-3
    vfb = 0.0
    ni = 1.5e10
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)

    def __init__(self, tox, gamma, Vcb_max, Vgb=0, Na=1e17):
        vt = 25.8e-3
        eps = 8.854e-14 * (11.8)
        epox = 8.854e-14 * (3.9)

        self.tox = tox
        self.Cox = epox / self.tox
        self.gamma = gamma
        self.Na = Na
        self.Vgb = Vgb
        self.Vcb_max = Vcb_max
        self.phi_list = []
        self.Vcb_list = []
        self.Vcb_regions = []
        self.phi_regions = []
        self.phi_vcb_map = {}
        self.count = 0
        self.phi_f = vt * m.log(self.Na / Surface_Potential_To_VCB.ni)

    def Psi_to_VCB(self):
        vt = 25.8e-3
        epox = 8.854e-14 * (3.9)
        eps = 8.854e-14 * (11.8)
        delta = 0.001
        phi = 0.05

        # print(self.phi_f + self.Vcb)
        self.Vcb = np.array(np.linspace(0, self.Vcb_max, 500))
        self.phi_regions = [self.phi_f, 2 * self.phi_f, (2 * self.phi_f) + (5 * vt)]
        # print(self.phi_regions)
        self.count = 0
        if (phi - vt) > 0:

            for j in self.Vcb:
                if self.count == 0:
                    phi = random.uniform((Surface_Potential_To_VCB.vfb * 0.1) + 1.25, (Surface_Potential_To_VCB.vfb * 0.1) + 1.35)
                    self.count = self.count + 1
                else:
                    phi = self.phi_list[-1]

                while 1:
                    first = m.sqrt(abs((vt * m.exp(-phi / vt) + (phi - vt) + (m.exp(-(2 * self.phi_f + j) / vt) * ((vt * m.exp(phi / vt)) - phi - vt)))))
                    # print(f'first: {first:.4f}')
                    second = m.sqrt(abs(vt * m.exp(-(phi + delta) / vt) + (phi + delta - vt) + m.exp(-(2 * self.phi_f + j) / vt) * ((vt * m.exp((phi + delta) / vt)) - (phi + delta) - vt)))
                    # print(f'second: {second:.4f}')

                    f1 = ((self.gamma * first) + phi + (Surface_Potential_To_VCB.vfb)) - self.Vgb
                    f2 = ((self.gamma * second) + (phi + delta) + (Surface_Potential_To_VCB.vfb)) - self.Vgb

                    test = phi - 0.1 * (f1 / diff(f2, f1, delta))
                    # print(f"Test Value: {test}")

                    if (abs((test - phi)) / phi) < 1e-2 and abs(f1) < 1e-3:
                        # print(f" ==  break_out === {break_out}")
                        break_out = break_out + 1

                        if break_out > 2:
                            break_out = 0
                            break

                    else:
                        phi = test
                        break_out = 0

                self.phi_list.append(test)
                self.Vcb_list.append(j)

        else:
            phi = phi + 0.001

        self.phi_list = np.array(self.phi_list)
        self.Vcb_list = np.array(self.Vcb_list)

        for i, j in zip(self.phi_list, self.Vcb_list):
            self.phi_vcb_map[i] = j

        self.Vu = ((- self.gamma / 2) + m.sqrt((self.gamma**2) / 4 + self.Vgb - Surface_Potential_To_VCB.vfb))**2 - self.phi_f
        self.Vcb_regions.append(self.Vu)
        self.Vw = ((- self.gamma / 2) + m.sqrt((self.gamma**2) / 4 + self.Vgb - Surface_Potential_To_VCB.vfb))**2 - (2 * self.phi_f)
        self.Vcb_regions.append(self.Vw)
        Vz = 0.5
        self.Vq = ((- self.gamma / 2) + m.sqrt((self.gamma**2) / 4 + self.Vgb - Surface_Potential_To_VCB.vfb - Vz))**2 - (2 * self.phi_f)
        self.Vcb_regions.append(self.Vq)

        # print(self.phi_vcb_map)
        # print(self.Vcb_regions)

        return (self.Vcb_regions, self.phi_vcb_map, self.phi_regions)
# # =====================================================================

    def Plot_Psi_Vs_VCB(self):
        vt = 25.8e-3
        common_plotting_factor = 2.2
        self.phi_sa = ((- self.gamma / 2) + m.sqrt((gamma**2) / 4 + self.Vgb - Surface_Potential_To_VCB.vfb))**2

        plt.plot(self.Vcb_list, self.phi_list)
        plt.title('Surface Potential  vs  Source-Body Voltage for 3-Terminal MOS')

        ###### Slope lines with VCB #############

        self.Vcb_plot = np.array(np.linspace(0, self.Vcb_max / 1.55, 500))
        plt.plot(self.Vcb_plot, self.Vcb_plot + 2 * self.phi_f, color='r', linestyle='--', alpha=0.5)
        plt.plot(self.Vcb_plot, self.Vcb_plot + self.phi_f, color='r', linestyle='--', alpha=0.5)

        #######  Lines and Labels ######
        ##  X & Y Axis ###
        plt.axhline(y=0, xmin=0.0, xmax=1.0, color='k')
        plt.axvline(x=0, ymin=0.0, ymax=1.0, color='k')

        ##### Axes Scale Limits ##
        plt.ylim(0,)
        # plt.xlim(0,)

        # Text mentioning the Depletion, Weak/Moderate/Strong Inversion regions
        plt.annotate("Depletion", xy=((self.Vcb_regions[0] + 1.5, (self.phi_f * common_plotting_factor * 2))))
        plt.annotate("Weak\nInversion", xy=((self.Vcb_regions[1] * 1.05, (self.phi_f * common_plotting_factor * 2))))
        plt.annotate("Moderate\nInversion", xy=((self.Vcb_regions[2], (2 * self.phi_f * common_plotting_factor))))
        plt.annotate("Strong\nInversion", xy=((self.Vcb_regions[2] - 0.76, (2 * self.phi_f * common_plotting_factor))))

        # # Horizontal lines at region boundaries
        plt.axhline(y=(2 * self.phi_f) + (5 * vt), color='g', linestyle='--', alpha=0.75)
        plt.axhline(y=(2 * self.phi_f), color='g', linestyle='--', alpha=0.75)
        plt.axhline(y=(self.phi_f), color='g', linestyle='--', alpha=0.75)

        # # Vertical lines at region boundaries
        plt.axvline(x=self.Vcb_regions[0], color='#968F8F', linestyle='--', alpha=0.75)
        plt.axvline(x=self.Vcb_regions[1], color='#968F8F', linestyle='--', alpha=0.75)
        plt.axvline(x=self.Vcb_regions[2], color='#968F8F', linestyle='--', alpha=0.75)

        # Printing the vlaue of boundaries of Strong/ Moderate / Weak Inversion in VCB terms
        # print(u"$V_{Q}$", end=' ')
        # print("= {}".format(self.Vcb_regions[0]))

        # print(r"$V_{W}$", end=' ')
        # print("= {}".format(self.Vcb_regions[1]))

        # print(r"$V_{U}$", end=' ')
        # print(" = {}".format(self.Vcb_regions[2]))

        #### Axis Labels ######
        plt.xlabel(r'$V_{CB}$ (V)')
        plt.ylabel('Surface Potential (V)')
        # plt.text(Vgb_regions[2] - vfb, phi_f + 0.1, r'$V_{FB}$ = - 1.0 V')

        #### Placing corresponding VL0, VM0, and VH0 at appropriate points on VGB-axis ####
        # plt.text(self.Vgb_regions[0] - Surface_Potential_To_VGB.vfb, 0.1, r'$V_{L0}$ =' + '{:.3f} V'.format(self.Vgb_regions[0]))
        # plt.text(self.Vgb_regions[1] - Surface_Potential_To_VGB.vfb, 0.1, r'$V_{M0}$ =' + '{:.3f} V'.format(self.Vgb_regions[1]))
        # plt.text(self.Vgb_regions[2] - Surface_Potential_To_VGB.vfb, 0.1, r'$V_{H0}$ =' + '{:.3f} V'.format(self.Vgb_regions[2]))

        # Placing VCB over the graph
        self.phi_sa = ((- self.gamma / 2) + m.sqrt((gamma**2) / 4 + self.Vgb - Surface_Potential_To_VCB.vfb))**2
        plt.text(self.Vcb_max / 1.2, self.phi_sa * common_plotting_factor * 0.5, r'$V_{GB}$ =' + '{:.2f} V'.format(self.Vgb))

        # Placing corresponding vlaues of phi, 2 phi, and (2 phi + 5vt) at appropriate points on surface potential-axis
        plt.text(0, self.phi_f - 0.1, r'$\phi_{F}$  = ' + '{:.3f} V'.format(self.phi_f))
        plt.text(0, (2 * self.phi_f) - 0.15, r'$2 \phi_{F}$  = ' + '{:.3f} V'.format(2 * self.phi_f))
        plt.text(0, (2 * self.phi_f + 5 * vt) - 0.1, r'$2 \phi_{F}$ + $5V_t$ = ' + '{:.3f} V'.format(2 * self.phi_f + 5 * vt))

       # plaicng the value of flat-band voltage(VFB) on the plot
        plt.text(self.Vcb_max - 1, self.phi_sa / 1.25, r'$V_{FB}$ =' + '{:.2f} V'.format(Surface_Potential_To_VCB.vfb))

        # plt.legend()
        plt.tight_layout()
        plt.show()
        return None


if __name__ == "__main__":
    tox = 40e-7
    Na = 1e17
    q = 1.6e-19
    eps = 8.854e-14 * (11.8)
    epox = 8.854e-14 * (3.9)
    Cox = epox / tox
    Vcb_max = 6
    Vgb = 8
    gamma = ((m.sqrt(2 * q * eps * Na)) / Cox)

    x = Surface_Potential_To_VCB(tox, gamma, Vcb_max, Vgb, Na)
    a, b, c = x.Psi_to_VCB()
    # print(a, b, c)
    # print(c)

    y = x.Plot_Psi_Vs_VCB()

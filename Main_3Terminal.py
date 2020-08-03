import numpy as np
import matplotlib.pyplot as plt
import math as m
from Mathematical_functions import diff
from Class_Surface_Potential_To_VGB_3Terminal import Surface_Potential_To_VGB
from Class_Surface_Potential_VGB_To_Total_Charge_3Terminal import Total_charge_VGB_Surface_Potential
from Class_Total_charge_VGB_3Terminal import Total_charge_VGB
from Class_Capacitance_To_VGB_3Terminal import Capacitance_To_VGB
from Class_Capacitance_To_Surface_Potential_3Terminal import Capacitance_To_Surface_Potential
from Class_Surface_Potential_To_VCB_3Terminal import Surface_Potential_To_VCB
from Class_Total_charge_VCB_3Terminal import Total_charge_VCB
import matplotlib.colors as colors
colors_list = list(colors._colors_full_map.values())

Na = 1e17
q = 1.6e-19
eps = 8.854e-14 * (11.8)
epox = 8.854e-14 * (3.9)

heading = """Amratansh Gupta developed this code for charge and capacitance modelling of 3-Terminal MOS Structure. Reuse or reproduction of the model without proper credits/citation is NOT allowed."""

print(heading.center(10))
print(f"""\nInstructions for using the model:\n
INPUTS -
1. Flat-band Voltage is implicitly set to 0.
2. Substrate doping ('Na') (enter in cm\u207B\u00B3).
3. Enter the value of Gate-oxide(tox) in nanometers.
4. Enter the value of maximum/fixed value of Gate-Body bias('VGB') in volt.
5. Enter the value of maximum/fixed value of Source-Body bias('VCB') in volt.\n

OUTPUTS - Enter the integer listed below in order to get the corresponding plot
1. Plot Surface Potential vs Gate-Body Voltage
2. Plot Charges - (Depletion, Inversion, and Total Semiconductor charge) vs Surface Potential
3. Plot Charges - (Depletion, Inversion, and Total Semiconductor charge) vs Gate-Body voltage
4. Plot Capacitance and its components vs Gate-Body Voltage
5. Plot Capacitance and its components vs Surface Potential
6. Plot Surface Potential vs Source-Body Voltage
7. Plot Charges - (Depletion, Inversion, and Total Semiconductor charge) vs Source-Body voltage

 """)

tox = float(input('Enter the oxide thickness (in nm): ')) * (1e-7)
Na = float(input('Enter the substrate dpoing: '))
Cox = epox / tox
# print(Cox)
gamma = ((m.sqrt(2 * q * eps * Na)) / Cox)
while True:
    try:
        plot_output = int(input('Enter the integer listed above in order to get the corresponding plot: '))
        if abs(plot_output) <= 7 and (plot_output) > 0:

            if abs(plot_output) == 1:
                Vgb_max = float(input('Enter the maximum value of VGB(in V): '))
                Vcb = float(input('Enter the fixed value of VCB(in V): '))
                x = Surface_Potential_To_VGB(tox, gamma, Vgb_max, Vcb, Na)
                y = x.Psi_to_VGB()
                print('Plotting')
                z = x.Plot_Psi_Vs_VGB()
                continue
            else:
                pass

            if abs(plot_output) == 2:
                Vgb_max = float(input('Enter the maximum value of VGB(in V): '))
                Vcb = float(input('Enter the fixed value of VCB(in V): '))
                x = Total_charge_VGB_Surface_Potential(tox, gamma, Vgb_max, Vcb, Na)
                y = x.Total_Charge()
                z = x.Qinversion()
                v = x.Qdepletion()
                print('Plotting')
                b = x.Plot_Q_vs_Surface_Potential()
                continue

            else:
                pass

            if abs(plot_output) == 3:
                Vgb_max = float(input('Enter the maximum value of VGB(in V): '))
                Vcb = float(input('Enter the fixed value of VCB(in V): '))
                x = Total_charge_VGB(tox, gamma, Vgb_max, Vcb, Na)
                y = x.Total_Charge()
                z = x.Qinversion()
                v = x.Qdepletion()
                print('Plotting')
                b = x.Plot_Q_vs_VGB()
                continue

            else:
                pass

            if abs(plot_output) == 4:
                Vgb_max = float(input('Enter the maximum value of VGB(in V): '))
                Vcb = float(input('Enter the fixed value of VCB(in V): '))
                x = Capacitance_To_VGB(tox, gamma, Vgb_max, Vcb, Na)
                y = x.Dep_Capacitance()
                z = x.Inv_Capacitance()
                a = x.Semiconductor_Capacitance()
                b = x.Total_Capacitance()
                print('Plotting')
                c = x.Plot_C_vs_VGB()
                continue
            else:
                pass

            if abs(plot_output) == 5:
                Vgb_max = float(input('Enter the maximum value of VGB(in V): '))
                Vcb = float(input('Enter the fixed value of VCB(in V): '))
                x = Capacitance_To_Surface_Potential(tox, gamma, Vgb_max, Vcb, Na)
                y = x.Dep_Capacitance()
                z = x.Inv_Capacitance()
                a = x.Semiconductor_Capacitance()
                b = x.Total_Capacitance()
                print('Plotting')
                c = x.Plot_C_vs_Surface_Potential()
                continue
            else:
                pass

            if abs(plot_output) == 6:
                Vgb = float(input('Enter the fixed value of VGB(in V): '))
                Vcb_max = float(input('Enter the maximum value of VCB(in V): '))
                x = Surface_Potential_To_VCB(tox, gamma, Vcb_max, Vgb, Na)
                y = x.Psi_to_VCB()
                print('Plotting')
                z = x.Plot_Psi_Vs_VCB()
                continue

            else:
                pass

            if abs(plot_output) == 7:
                Vcb_max = float(input('Enter the maximum value of VCB(in V): '))
                Vgb = float(input('Enter the fixed value of VGB(in V): '))
                x = Total_charge_VCB(tox, gamma, Vcb_max, Vgb, Na)
                y = x.Total_Charge()
                z = x.Qinversion()
                v = x.Qdepletion()
                print('Plotting')
                b = x.Plot_Q_vs_VCB()
                continue

            else:
                pass

        else:
            print('Input Integer exceeds the allowed range')

    except Exception as e:
        print(e)
        print('Invalid Input - Non-Integer Input.')
        continue

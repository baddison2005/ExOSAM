#Rossiter-McLaughlin Effect Analysis Python Wrapper for Exosam (RMLEAPWE).

#This python wrapper is used to determine the best fit parameters; lambda (the spin-orbit angle) and vsini (the
#rotational velocity of the host star), as well as their associated uncertainties by fitting for the Rossiter-McLaughlin
#effect velocity anomaly (RM effect). Please view the readme file for how to use ExOSAM and the steps taken to determine
#the best fit lambda and vsini and their associated uncertainties.



#---------------------------------------Import required Python packages.-----------------------------------------------#
import os, sys
import math
import random
import numpy as np
import pandas as pd
import pyexcel                        # import it to handle CSV files.
from tkinter import *
import tkinter as tk
import matplotlib.pyplot as plt
import subprocess
import fortranformat as ff
from astropy.stats import knuth_bin_width
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
#from joblib import Parallel, delayed
from tkinter import filedialog
sys.path.append(os.getcwd())
from RV_offset_RM_fit import RV_offset
from RV_RM_transit_curve import RV_model_array
from matplotlib.colors import LogNorm
from matplotlib.pyplot import *
import scipy.optimize as so
from scipy.ndimage.filters import gaussian_filter
import corner
from matplotlib.ticker import AutoMinorLocator




#-------------------------------------Constants and initial directory info---------------------------------------------#
Rss = 6.9634E8                        #Radius of our Sun (in meters).
AU = 1.4960E11                        #One astronomical unit in meters.
pi = math.pi
G = 6.6738E-11                        #Gravitation constant.
Mss = 1.9891E30                       #Mass of the sun in kilograms.
Time = 0                              #The time for the FOR loop.
Me = 5.9722E24                        #Mass of the Earth in kg.
Mj = 1.8986E27                        #Mass of Jupiter in kg.
Rj = 7.1492E7                         #Radius of Jupiter (in meters).
RE = 6.3781E6                         #Radius of Earth (in meters if planet is given in Earth radii).
Long_ascend_node = 90                 #Set to 90. Only measuring the relative or projected inclination of orbit.
day_sec = 86400.0

#Grab RM directory info and parameter input file.
root = tk.Tk()
root.withdraw()
RM_directory_info = filedialog.askopenfilename(parent=root, title='Please select the input parameter file.')

directory_input = RM_directory_info.rsplit('/', 1)[0] + '/'
parameter_filename = RM_directory_info.rsplit('/', 1)[1]
input_parameter_filename = RM_directory_info

#Read in the input model parameters values.
param = np.genfromtxt(input_parameter_filename, dtype='str')




#----------------------------------------------Custom Functions------------------------------------------------------#
def to_bool(value):
    """
       Converts 'something' to boolean. Raises exception for invalid formats
           Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
           Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    """
    if str(value).lower() in ("yes", "y", "true",  "t", "1"): return True
    if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): return False
    raise Exception('Invalid value for boolean conversion: ' + str(value))
#enddef

def isint(value):
    try:
        int(value)
        return True
    except:
        return False
#enddef

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level




#-----------------------Set the values for the parameters from the input file.-----------------------------------------#
directory_location = param[0][2]
Planet_name = param[1][2]
other_RV_files = param[2][2]
output_temp_data = param[3][2]
Jupiter_units = param[4][2]
Bessel_function_exit = float(param[5][2])
phase_shift = param[6][2]
subtract_JD_mid_time = param[7][2]
add_RV_zero_offset = param[8][2]
Number_orbits = float(param[9][2])
data_plot_model_interval = float(param[10][2])
Phase_angle_start = float(param[11][2])
Time_bef_aft = float(param[12][2])
Time_compare_vel = float(param[13][2])
my_data_plotting_symbols = param[14][2]
other_data_plotting_symbols = param[15][2]
Ms_solar_prior = float(param[16][2])
Ms_solar_end = float(param[17][2])
Ms_solar_begin = float(param[18][2])
Ms_solar_1sigerr = float(param[19][2])
Ms_solar_norm_prior = param[20][2]
Ms_solar_fixed = param[21][2]
Rs_solar_prior = float(param[22][2])
Rs_solar_end = float(param[23][2])
Rs_solar_begin = float(param[24][2])
Rs_solar_1sigerr = float(param[25][2])
Rs_solar_norm_prior = param[26][2]
Rs_solar_fixed = param[27][2]
Rp_Rs_ratio_prior = float(param[28][2])
Rp_Rs_ratio_end = float(param[29][2])
Rp_Rs_ratio_begin = float(param[30][2])
Rp_Rs_ratio_1sigerr = float(param[31][2])
Rp_Rs_ratio_norm_prior = param[32][2]
Rp_Rs_ratio_fixed = param[33][2]
use_Rp_Rs_ratio = param[34][2]

vmacro_prior = float(param[35][2])
vmacro_end = float(param[36][2])
vmacro_begin = float(param[37][2])
vmacro_1sigerr = float(param[38][2])
vmacro_norm_prior = param[39][2]
vmacro_fixed = param[40][2]

beta_rot = float(param[41][2])
M_pixels = int(param[42][2])
linear_quadratic = param[43][2]
u_linear = float(param[44][2])

q_1_prior = float(param[45][2])
q_1_end = float(param[46][2])
q_1_begin = float(param[47][2])
q_1_1sigerr = float(param[48][2])

q_2_prior = float(param[49][2])
q_2_end = float(param[50][2])
q_2_begin = float(param[51][2])
q_2_1sigerr = float(param[52][2])
q_norm_prior = param[53][2]
q_fixed = param[54][2]

Inc_prior = float(param[55][2])
Inc_end = float(param[56][2])
Inc_begin = float(param[57][2])
Inc_1sigerr = float(param[58][2])
Inc_norm_prior = param[59][2]
Inc_fixed = param[60][2]

Impact_prior = float(param[61][2])
Impact_end = float(param[62][2])
Impact_begin = float(param[63][2])
Impact_1sigerr = float(param[64][2])
Impact_norm_prior = param[65][2]
Impact_fixed = param[66][2]
use_Impact = param[67][2]

Rorb_Rs_prior = float(param[68][2])
Rorb_Rs_end = float(param[69][2])
Rorb_Rs_begin = float(param[70][2])
Rorb_Rs_1sigerr = float(param[71][2])
Rorb_Rs_norm_prior = param[72][2]
Rorb_Rs_fixed = param[73][2]
use_Rorb_Rs_ratio = param[74][2]

arg_periastron_prior = float(param[75][2])
arg_periastron_end = float(param[76][2])
arg_periastron_begin = float(param[77][2])
arg_periastron_1sigerr = float(param[78][2])
arg_periastron_norm_prior = param[79][2]
arg_periastron_fixed = param[80][2]
Ecc_prior = float(param[81][2])
Ecc_end = float(param[82][2])
Ecc_begin = float(param[83][2])
Ecc_1sigerr = float(param[84][2])
Ecc_norm_prior = param[85][2]
Ecc_fixed = param[86][2]
Mp_prior = float(param[87][2])
Mp_end = float(param[88][2])
Mp_begin = float(param[89][2])
Mp_1sigerr = float(param[90][2])
Mp_norm_prior = param[91][2]
Mp_fixed = param[92][2]
Rp_prior = float(param[93][2])
Rp_end = float(param[94][2])
Rp_begin = float(param[95][2])
Rp_1sigerr = float(param[96][2])
Rp_norm_prior = param[97][2]
Rp_fixed = param[98][2]
Orbital_period_prior = float(param[99][2])
Orbital_period_end = float(param[100][2])
Orbital_period_begin = float(param[101][2])
Orbital_period_1sigerr = float(param[102][2])
Orbital_period_norm_prior = param[103][2]
Orbital_period_fixed = param[104][2]
JD_time_mid_transit_prior = float(param[105][2])
JD_time_mid_transit_end = float(param[106][2])
JD_time_mid_transit_begin = float(param[107][2])
JD_time_mid_transit_1sigerr = float(param[108][2])
JD_time_mid_transit_norm_prior = param[109][2]
JD_time_mid_transit_fixed = param[110][2]
Albedo = float(param[111][2])
RV_zero_offset_prior = float(param[112][2])
RV_zero_offset_end = float(param[113][2])
RV_zero_offset_begin = float(param[114][2])
RV_zero_offset_interval = float(param[115][2])
RV_zero_offset_1sigerr = float(param[116][2])
RV_zero_offset_norm_prior = param[117][2]
RV_zero_offset_fixed = param[118][2]
RV_offset_datasets_prior = float(param[119][2])
RV_offset_datasets_end = float(param[120][2])
RV_offset_datasets_begin = float(param[121][2])
RV_offset_datasets_interval = float(param[122][2])
RV_offset_datasets_1sigerr = float(param[123][2])
RV_offset_datasets_norm_prior = param[124][2]
RV_offset_datasets_fixed = param[125][2]

K_amp_prior = float(param[126][2])
K_amp_end = float(param[127][2])
K_amp_begin = float(param[128][2])
K_amp_1sigerr = float(param[129][2])
K_amp_norm_prior = param[130][2]
K_amp_fixed = param[131][2]
use_K_amp = param[132][2]

vsini_prior = float(param[133][2])
vsini_end = float(param[134][2])
vsini_begin = float(param[135][2])
vsini_1sigerr = float(param[136][2])
vsini_norm_prior = param[137][2]
impose_prior_vsini = param[138][2]
vsini_fixed = param[139][2]
stellar_rotation_angle_prior = float(param[140][2])
stellar_rotation_angle_end = float(param[141][2])
stellar_rotation_angle_begin = float(param[142][2])
stellar_rotation_angle_1sigerr = float(param[143][2])
stellar_rotation_angle_norm_prior = param[144][2]
impose_prior_stellar_rotation_angle = param[145][2]
stellar_rotation_angle_fixed = param[146][2]
mcmc_accepted_iteration_size = int(param[147][2])
number_mcmc_walkers = int(param[148][2])
scale_factor = float(param[149][2])
chi_squared_change = float(param[150][2])
chi_squared_change_fit = float(param[151][2])
use_out_transit_rv_for_fit = param[152][2]
test_null = param[153][2]





#---------------------------------------------Set required working paths-----------------------------------------------#
save_data_directory = directory_location
output_data_directory = output_temp_data
RM_data_output = output_data_directory + 'data_array.txt'
other_data_output = output_data_directory + 'other_data_array.txt'
my_data_output = output_data_directory + 'my_data_array.txt'




#-------------------------------------------------Select RV data-------------------------------------------------------#
root = tk.Tk()
root.withdraw()
RM_data = filedialog.askopenfilename(initialdir = directory_location, parent=root, title='Choose your RV dataset')




#--------------------------------------------------Read in RV data-----------------------------------------------------#
template_RV_sheet = pyexcel.get_sheet(file_name=RM_data, name_columns_by_row=0)
num_my_rv = pyexcel.Sheet.number_of_rows(template_RV_sheet)
my_datafilelength = num_my_rv
Index = np.arange(num_my_rv)
template_RV = pd.DataFrame(index=Index, columns=['Time', 'RV', 'RV_error'])
template_RV['Time'] = template_RV_sheet.column[0]
template_RV['RV'] = template_RV_sheet.column[1]
template_RV['RV_error'] = template_RV_sheet.column[2]
template_RV = template_RV.astype(float)



#---------------------------------------------Read in other RV data----------------------------------------------------#
if other_RV_files == 'Y':
    root = tk.Tk()
    root.withdraw()
    other_data_file_name = filedialog.askopenfilename(initialdir = directory_location, parent=root,
                                                      title='Choose an other RV dataset')
    template_other1_sheet = pyexcel.get_sheet(file_name=other_data_file_name, name_columns_by_row=0)
    num_rv1 = pyexcel.Sheet.number_of_rows(template_other1_sheet)
    datafilelength_other1 = num_rv1
    Index_other = np.arange(num_rv1)
    template_other1 = pd.DataFrame(index=Index_other, columns=['Time', 'RV', 'RV_error'])
    template_other1['Time'] = template_other1_sheet.column[0]
    template_other1['RV'] = template_other1_sheet.column[1]
    template_other1['RV_error'] = template_other1_sheet.column[2]
    template_other1 = template_other1.astype(float)
#endif
else:
    num_rv1 = 0
#endelse




#---------------------------------------Planetary orbital calculations-------------------------------------------------#



Orbital_period_day = Orbital_period_prior
Orbital_period = Orbital_period_prior*3600.0*24.0
JD_time_mid_transit = JD_time_mid_transit_prior

Ms = Mss * Ms_solar_prior                      #Mass of the star (in kg).
if use_Impact == 'Y' and use_Rorb_Rs_ratio == 'Y':
    Inc_prior = math.acos((Impact_prior/(Rorb_Rs_prior))*((1.0 + (Ecc_prior*math.sin(arg_periastron_prior*
                (pi/180.0))))/(1.0 - Ecc_prior**2.0)))*(180.0/pi)
    Transit_length = (Orbital_period/pi) * math.asin((1.0/Rorb_Rs_prior)* \
                           (math.sqrt((1.0 + Rp_Rs_ratio_prior)**2.0 - Impact_prior**2.0)/
                            math.sin(Inc_prior*(pi/180.0)))) * (math.sqrt(1.0 - Ecc_prior**2.0)/
                            (1.0 + (Ecc_prior*math.sin(arg_periastron_prior*(pi/180.0)))))
#endif

if use_K_amp == 'Y' and Mp_fixed == 'Y' and Ecc_prior == 0:
    Mp = (Orbital_period/(2.0*pi*G))**(1.0/3.0)*((Ms**(2.0/3.0)*K_amp_prior)/math.sin(Inc_prior*(pi/180.0)))
else:
    if Jupiter_units == 'Y':
        Mp = Mp_prior * Mj
    else:
        Mp = Mp_prior * Me
    #endif
#endif

Rorb = ((Orbital_period ** (2.0) * G * (Ms + Mp)) / (4.0 * pi ** (2.0))) ** (1.0 / 3.0)  # Semi-major axis
Rorb_star = ((Orbital_period ** (2.0) * G * ((Mp ** (3.0)) / (Ms + Mp) ** (2.0))) / (4.0 * pi ** (2.0))) ** (1.0 / 3.0)  # Star Semi-major axis.

if use_Rorb_Rs_ratio == 'Y':
    Rs = Rorb/(Rorb_Rs_prior)
else:
    Rs = Rss * Rs_solar_prior
#endif

if use_Rp_Rs_ratio == 'Y':
    Rp = Rp_Rs_ratio_prior * Rs
else:
    if Jupiter_units == 'Y':
        Rp = Rp_prior * Rj
    else:
        Rp = Rp_prior * RE
    #endelse
#endif

RV = 0.0                                       #Set the RV to zero.
Rs2 = Rs**2.0                                   #Variable to speed up calculations.
Rp2 = Rp**2.0                                   #Variable to speed up calculations.

Ecc = Ecc_prior
Inc = Inc_prior

if use_K_amp == 'Y':
    RVamplitude = K_amp_prior
else:
    if Ecc == 0:
        #Maximum amplitude caused by the exoplanet in a circular orbit.
        #RVamplitude = (Mp/Ms)*sqrt((G*(Ms + Mp))/Rorb)*sin(Inc*(pi/180.0D0))
        RVamplitude = math.sqrt((G*(Mp**(3.0))*(math.sin(Inc*(pi/180.0)))**(3.0))/
                      (Rorb_star*math.sin(Inc*(pi/180.0))*(Ms + Mp)**(2.0)))
    else:
        #Maximum amplitude caused by the exoplanet in an eccentric orbit.
        #RVamplitude = (1.0D0/sqrt(1.0D0 - Ecc^2.0D0))*(Mp/Ms)*sqrt((G*(Ms + Mp))/Rorb)*sin(Inc*(pi/180.0D0))
        RVamplitude = math.sqrt((G*(Mp**(3.0))*(math.sin(Inc*(pi/180.0)))**(3.0))/
                      ((1.0 - Ecc**2.0)*Rorb_star*math.sin(Inc*(pi/180.0))*(Ms + Mp)**(2.0)))
    #endelse
#endelse




#--------------------------------------Apply any offsets to your RV data-----------------------------------------------#
if add_RV_zero_offset == 'Y':
    #Apply velocity offset to RV data.
    template_RV.loc[:, 'RV'] = template_RV.loc[:, 'RV'] + RV_zero_offset_prior
#endif

#If subtract_JD_mid_time is set to 'Y', then set the data time relative to JD mid transit time.
if subtract_JD_mid_time == 'Y':
    #Phase shift the data so that it fits in one orbital period, making sure that the phase is between 0.5 - 1.5.
    #Assuming time is HJD format.
    for i in range(0,my_datafilelength - 1):
        Time_minus_midtransit = template_RV.loc[i,'Time'] - JD_time_mid_transit
        if (Time_minus_midtransit/Orbital_period_day) - math.floor(Time_minus_midtransit/Orbital_period_day) < 0.5:
            phase = (Time_minus_midtransit/Orbital_period_day) - math.floor(Time_minus_midtransit/Orbital_period_day) \
                    + 1.0
        #endif
        else:
            phase = (Time_minus_midtransit/Orbital_period_day) - math.floor(Time_minus_midtransit/Orbital_period_day)
        #endelse
        template_RV.loc[i,'Time'] = (Orbital_period_day * phase) - Orbital_period_day
    #endfor
#endif

#Sort the RV array in time to feed into Fortran routine.
template_RV = template_RV.sort_values('Time', ascending=True)
#template_RV = template_RV.reset_index(drop=True)



#--------------------------------------Apply any offsets to other RV data----------------------------------------------#
if other_RV_files == 'Y':
    if add_RV_zero_offset == 'Y':
        # Apply velocity offset to RV data.
        template_other1.loc[:, 'RV'] = template_other1.loc[:, 'RV'] + RV_zero_offset_prior
    #endif
    if subtract_JD_mid_time == 'Y':
        # Phase shift the data so that it fits in one orbital period, making sure that the phase is between 0.5 - 1.5.
        # Assuming time is HJD format.
        for i in range(0, datafilelength_other1 - 1):
            Time_minus_midtransit = template_other1.loc[i, 'Time'] - JD_time_mid_transit
            if (Time_minus_midtransit / Orbital_period_day) - math.floor(
                Time_minus_midtransit / Orbital_period_day) < 0.5:
                phase = (Time_minus_midtransit / Orbital_period_day) - math.floor(
                         Time_minus_midtransit / Orbital_period_day) + 1.0
            #endif
            else:
                phase = (Time_minus_midtransit / Orbital_period_day) - math.floor(
                         Time_minus_midtransit / Orbital_period_day)
            #endelse
            template_other1.loc[i, 'Time'] = (Orbital_period_day * phase) - Orbital_period_day
        #endfor
    #endif

    #Sort the RV array in time to feed into Fortran routine.
    template_other1 = template_other1.sort_values('Time', ascending=True)
    #template_other1 = template_other1.reset_index(drop=True)
#endif

#---------------------------------------Read in plotting symbols for your RVs------------------------------------------#
my_data_plotting_symbols_array = pd.DataFrame(index=Index, columns=['Symbols'])

if len(my_data_plotting_symbols) >= 1:
    #Read in the plotting symbols.
    f = open(my_data_plotting_symbols, 'r')
    num_my_symbols = 0
    for line in f:
        line = line.strip()
        columns = line.split()
        if isint(columns[0]) == True:
            my_data_plotting_symbols_array.iloc[num_my_symbols] = int(columns[0])
        else:
            my_data_plotting_symbols_array.iloc[num_my_symbols] = columns[0]
        #endelse
        num_my_symbols = num_my_symbols + 1
    #endfor
    f.close()

    if num_my_symbols < my_datafilelength:
        #Less symbols than data points. Fill in the missing symbols using the last symbol in the array.
        for i in range(num_my_symbols, my_datafilelength - 1):
            my_data_plotting_symbols_array.iloc[i] = my_data_plotting_symbols_array.iloc[num_my_symbols - 1]
        # endfor
    #endif
#endif
else:
    my_data_plotting_symbols_array[:] = 16
#endelse

#rearrange order of plotting symbols to match the reordering done on RV dataset.
my_data_plotting_symbols_array = my_data_plotting_symbols_array.iloc[list(template_RV.index)]

template_RV = template_RV.reset_index(drop=True)
my_data_plotting_symbols_array = my_data_plotting_symbols_array.reset_index(drop=True)



#---------------------------------------Read in plotting symbols of other RVs------------------------------------------#
if other_RV_files == 'Y':
    other_data_plotting_symbols_array = pd.DataFrame(index=Index_other, columns=['Symbols'])

    if len(other_data_plotting_symbols) >= 1:
        #Read in the plotting symbols.
        f = open(other_data_plotting_symbols, 'r')
        num_other_symbols = 0
        for line in f:
            line = line.strip()
            columns = line.split()
            if isint(columns[0]) == True:
                other_data_plotting_symbols_array.iloc[num_other_symbols] = int(columns[0])
            else:
                other_data_plotting_symbols_array.iloc[num_other_symbols] = columns[0]
            #endelse
            num_other_symbols = num_other_symbols + 1
        # endfor
        f.close()

        if num_other_symbols < datafilelength_other1:
            #Less symbols than data points. Fill in the missing symbols using the last symbol in the array.
            for i in range(num_other_symbols, datafilelength_other1 - 1):
                other_data_plotting_symbols_array.iloc[i] = other_data_plotting_symbols_array.iloc[
                                                            num_other_symbols - 1]
            # endfor
        #endif
    #endif
    else:
        other_data_plotting_symbols_array[:] = 36
    # endelse

    # rearrange order of plotting symbols to match the reordering done on RV dataset.
    other_data_plotting_symbols_array = other_data_plotting_symbols_array.iloc[list(template_other1.index)]

    template_other1 = template_other1.reset_index(drop=True)
    other_data_plotting_symbols_array = other_data_plotting_symbols_array.reset_index(drop=True)
#endif
total_datalength = num_my_rv + num_rv1




#--------------------------------------Set the required units for simulation-------------------------------------------#
omega_arg_periastron_prior = arg_periastron_prior - 270.0 + Long_ascend_node         #Argument of the periastron as
                                                                                     #as measured from the transit mid
                                                                                     #mid time.
omega_arg_periastron_end = arg_periastron_end - 270.0 + Long_ascend_node             #Argument of the periastron.
omega_arg_periastron_begin = arg_periastron_begin - 270.0 + Long_ascend_node         #Argument of the periastron.

#omega_arg_periastron_prior = arg_periastron_prior + Long_ascend_node         #Argument of the periastron.
#omega_arg_periastron_end = arg_periastron_end + Long_ascend_node             #Argument of the periastron.
#omega_arg_periastron_begin = arg_periastron_begin + Long_ascend_node         #Argument of the periastron.

omega_arg_periastron = omega_arg_periastron_prior

True_anomaly_start = pi - omega_arg_periastron*(pi/180.0)

if True_anomaly_start >= pi:
    ecc_anomaly_start = 2.0*math.atan(math.tan(True_anomaly_start/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) + 2.0*pi
#endif
elif True_anomaly_start <= -pi:
    ecc_anomaly_start = 2.0*math.atan(math.tan(True_anomaly_start/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) - 2.0*pi
#endif
else:
    ecc_anomaly_start = 2.0*math.atan(math.tan(True_anomaly_start/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc)))
#endelse

#Calculate the amount of time (in seconds) after the passage of periastron has occurred when the planet is at mid
#transit. This is determined from the transit mid time and the argument of periastron.
if omega_arg_periastron > 180.0:
    Mean_anomaly_transit = (2.0*pi) - (omega_arg_periastron*(pi/180.0)) + pi
#endif
else:
    Mean_anomaly_transit = pi - (omega_arg_periastron*(pi/180.0))
#endelse

JD_time_peri = (JD_time_mid_transit*24.0*3600.0) - ((Mean_anomaly_transit*Orbital_period)/(2.0*pi))

time_peri_passage = (JD_time_mid_transit*24.0*3600.0) - JD_time_peri

print("time_peri_passage ", time_peri_passage)

if Ecc == 0:
    Time_start = ((ecc_anomaly_start*Orbital_period)/(2.0*pi)) + time_peri_passage
#endif
else:
    Time_start = (((ecc_anomaly_start - (Ecc*math.sin(ecc_anomaly_start)))*Orbital_period)/(2.0*pi)) + time_peri_passage
#endelse
print("Time_start", Time_start)

True_anomaly_transit = ((3.0*pi)/2.0) - (omega_arg_periastron*(pi/180.0))

if True_anomaly_transit >= pi:
    ecc_anomaly_transit = 2.0*math.atan(math.tan(True_anomaly_transit/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) + 2.0*pi
#endif
elif True_anomaly_transit <= -pi:
    ecc_anomaly_transit = 2.0*math.atan(math.tan(True_anomaly_transit/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) - 2.0*pi
#endif
else:
    ecc_anomaly_transit = 2.0*math.atan(math.tan(True_anomaly_transit/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc)))
#endelse

if Ecc == 0:
    Time_transit = ((ecc_anomaly_transit*Orbital_period)/(2.0*pi)) + time_peri_passage
#endif
else:
    Time_transit = (((ecc_anomaly_transit - (Ecc*math.sin(ecc_anomaly_transit)))*Orbital_period)/(2.0*pi)) \
                   + time_peri_passage
#endelse
print("Time_transit", Time_transit)




if use_Rp_Rs_ratio=='Y':
    radius_ratio = Rp_Rs_ratio_prior
else:
    radius_ratio = Rp / Rs
#endelse

if use_Impact == 'N' or use_Rorb_Rs_ratio == 'N':
    impact_prior = ((Rorb * math.cos(Inc * (pi / 180.0))) / Rs) * ((1.0 - Ecc ** 2.0)
                    / (1.0 + (Ecc * math.sin(omega_arg_periastron * (pi / 180.0)))))

    Transit_length = ((Orbital_period) / pi) * math.asin((Rs / Rorb) * (math.sqrt(abs((1.0 + radius_ratio) ** 2.0
                        - impact_prior ** 2.0)) / math.sin(Inc * (pi / 180.0)))) * (math.sqrt(1.0 - Ecc ** 2.0) /
                        (1.0 + (Ecc * math.sin(omega_arg_periastron * (pi / 180.0)))))
#END IF

time_avoid = (Transit_length/2.0) + Time_compare_vel
print("Length of time to exclude from comparison (i.e. during planetary transit): ", time_avoid)

RV_offset_datasets_interval_size = abs(RV_offset_datasets_end - RV_offset_datasets_begin)

if other_RV_files == 'Y' and RV_offset_datasets_interval_size > 0:
    #Determine the number of iterations that will be used in this model.
    if RV_offset_datasets_interval > 0 and RV_zero_offset_interval > 0:

        Number_iterations_total = int(((RV_offset_datasets_end - RV_offset_datasets_begin)/RV_offset_datasets_interval)
                                  +1)*int(((RV_zero_offset_end - RV_zero_offset_begin)/RV_zero_offset_interval) + 1)

        Number_RV_offset_iterations = int(((RV_offset_datasets_end - RV_offset_datasets_begin) /
                                           RV_offset_datasets_interval)+ 1)

        Number_RV_zero_offset_iterations = int(((RV_zero_offset_end - RV_zero_offset_begin) /
                                                RV_zero_offset_interval) + 1)

    #endif
    elif RV_offset_datasets_interval > 0 and RV_zero_offset_interval == 0:

        Number_iterations_total = int(((RV_offset_datasets_end - RV_offset_datasets_begin)/RV_offset_datasets_interval)
                                  +1)*int(((RV_zero_offset_end - RV_zero_offset_begin)/1) + 1)

        Number_RV_offset_iterations = int(((RV_offset_datasets_end - RV_offset_datasets_begin) /
                                           RV_offset_datasets_interval) + 1)

        Number_RV_zero_offset_iterations = int(((RV_zero_offset_end - RV_zero_offset_begin)/1) + 1)
    #endif
    elif RV_offset_datasets_interval == 0 and RV_zero_offset_interval > 0:

        Number_iterations_total = int(((RV_offset_datasets_end - RV_offset_datasets_begin)/1)
                                  +1) * int(((RV_zero_offset_end - RV_zero_offset_begin)/RV_zero_offset_interval) + 1)

        Number_RV_offset_iterations = int(((RV_offset_datasets_end - RV_offset_datasets_begin)/1) + 1)

        Number_RV_zero_offset_iterations = int(((RV_zero_offset_end - RV_zero_offset_begin)/RV_zero_offset_interval)+1)

    # endif
    else:
        Number_iterations_total = int(((RV_offset_datasets_end - RV_offset_datasets_begin)/1) + 1) * int(
                                  ((RV_zero_offset_end - RV_zero_offset_begin) /1) + 1)

        Number_RV_offset_iterations = int(((RV_offset_datasets_end - RV_offset_datasets_begin)/1) + 1)

        Number_RV_zero_offset_iterations = int(((RV_zero_offset_end - RV_zero_offset_begin)/1) + 1)
    #endelse




    RV_offset_datasets_new_prior, RV_offset_datasets_new_end, RV_offset_datasets_new_begin, \
    RV_offset_datasets_new_interval, RV_offset_datasets_new_1sigerr = \
    RV_offset(Index, template_RV, Index_other, template_other1, Time_transit, Number_iterations_total,
              Number_RV_offset_iterations, Number_RV_zero_offset_iterations, my_datafilelength, RV_zero_offset_end,
              RV_zero_offset_interval, RV_zero_offset_prior, RV_offset_datasets_end, RV_offset_datasets_interval,
              RV_offset_datasets_begin, datafilelength_other1, Ecc, Orbital_period, time_peri_passage,
              Bessel_function_exit, Rorb, Rorb_star, omega_arg_periastron, Inc, time_avoid, Rs, Rp, RVamplitude,
              num_my_rv, num_rv1, directory_location, Planet_name, my_data_plotting_symbols_array,
              other_data_plotting_symbols_array, RM_data_output, other_data_output)

#endif




if other_RV_files != 'Y' or RV_offset_datasets_interval_size == 0:
    #This part outputs all the values found in the model in a text file.
    offset_between_datasets = open(directory_location + 'offset_between_datasets.txt', 'w')
    offset_between_datasets.write('No other datasets to determine appropriate RV offsets.')
    offset_between_datasets.close()

    #change the IDL plotting symbol markers to python symbol markers.
    for i in range(num_my_rv):
        if my_data_plotting_symbols_array.loc[i, 'Symbols'] == 16:
            #The symbol is a cirlce.
            my_data_plotting_symbols_array.loc[i, 'Symbols'] = 'o'
        #endif
        if my_data_plotting_symbols_array.loc[i, 'Symbols'] == 17:
            #The symbol is a triangle.
            my_data_plotting_symbols_array.loc[i, 'Symbols'] = '^'
        #endif
    #endfor

    # Output array data for Fortran to read in.
    RM_data_output_file = open(RM_data_output, 'w')
    for i in range(num_my_rv):
        RM_data_output_file.write('{:<50.10e}     {:<50.5f}     {:<50.5f}\n'.format(template_RV.loc[i, 'Time'],
                                                                                    template_RV.loc[i, 'RV'],
                                                                                    template_RV.loc[i, 'RV_error']))
    # endfor
    RM_data_output_file.close()

    if other_RV_files == 'Y':

        for i in range(num_rv1):
            if other_data_plotting_symbols_array.loc[i, 'Symbols'] == 36:
                # The symbol is a cirlce.
                other_data_plotting_symbols_array.loc[i, 'Symbols'] = 'x'
            # endif
        # endfor

        # Output array data for Fortran to read in.
        RM_data_output_file2 = open(other_data_output, 'w')
        for i in range(num_rv1):
            RM_data_output_file2.write('{:<50.10e}     {:<50.5f}     {:<50.5f}\n'.format(template_other1.loc[i, 'Time'],
                                                                                         template_other1.loc[i, 'RV'],
                                                                                         template_other1.loc[
                                                                                             i, 'RV_error']))
        # endfor
        RM_data_output_file2.close()

    #endif

    RV_offset_datasets_new_prior = RV_offset_datasets_prior
    RV_offset_datasets_new_end = RV_offset_datasets_end
    RV_offset_datasets_new_begin = RV_offset_datasets_begin
    RV_offset_datasets_new_interval = RV_offset_datasets_interval
    RV_offset_datasets_new_1sigerr = RV_offset_datasets_1sigerr
#endif




#--------------------------------------Output any remaining files for Fortran------------------------------------------#
#plot just your data.
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
for i in range(num_my_rv):
    plt.errorbar(template_RV.loc[i, "Time"], template_RV.loc[i, "RV"],
                 yerr=template_RV.loc[i, "RV_error"], fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'],
                 ecolor='r', color='k', elinewidth=2, markersize=10, capsize=7)
#endfor
plt.xlabel('Time (HJD)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
plt.savefig(directory_location + 'My_data_RVs_' + Planet_name + '.pdf',dpi=150)
plt.close(fig)

if other_RV_files == 'Y' and RV_offset_datasets_interval_size == 0:
    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    for i in range(num_my_rv):
        plt.errorbar(template_RV.loc[i, "Time"], template_RV.loc[i, "RV"], yerr=template_RV.loc[i, "RV_error"],
                     fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='r', mfc='k', elinewidth=2,
                     markersize=6, capsize=7)
    # endfor
    for i in range(num_rv1):
        plt.errorbar(template_other1.loc[i, "Time"], template_other1.loc[i, "RV"],
                     yerr=template_other1.loc[i, "RV_error"], fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='b', color='k', elinewidth=2, markersize=10, capsize=7)
    # endfor
    plt.xlabel('Time (HJD)', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
    plt.savefig(directory_location + 'All_data_RVs_' + Planet_name + '.pdf', dpi=150)
    plt.close(fig)
#endif

#Output your data array for fortran to read in.
RM_data_output_file1 = open(my_data_output, 'w')
for i in range(num_my_rv):
    RM_data_output_file1.write('{:<50.10e}     {:<50.5f}     {:<50.5f}\n'.format(template_RV.loc[i, 'Time'],
                                                                                 template_RV.loc[i, 'RV'],
                                                                                 template_RV.loc[i, 'RV_error']))
#endfor
RM_data_output_file1.close()



#Output file pointing towards parameter.txt file for Fortran to read.
pointer_file_output = open(os.getcwd() + '/pointer.txt', 'w')
pointer_file_output.write('{:>150}\n'.format(directory_input + 'parameters.txt'))
pointer_file_output.close()




#Set test_null_hypothesis to false first.
test_null_hypothesis_text = 'F'
test_null_file_output = open(output_temp_data + 'test_null.txt', 'w')
test_null_file_output.write('{:1}\n'.format(test_null_hypothesis_text))
test_null_file_output.close()

#Output data for Fortran to read in.
param_output = open(directory_input + 'parameters.txt', 'w')
param_output.write('{:>150}\n'.format(RM_data_output))
param_output.write('{:>150}\n'.format(my_data_output))
param_output.write('{:>150}\n'.format(other_data_output))
param_output.write('{:<20d}\n'.format(num_my_rv))
param_output.write('{:<20d}\n'.format(num_rv1))
param_output.write('{:<20d}\n'.format(total_datalength))
param_output.write('{:>150}\n'.format(directory_location))
param_output.write('{:>20}\n'.format(Planet_name))
param_output.write('{:1}\n'.format(other_RV_files))
param_output.write('{:>150}\n'.format(output_temp_data))

param_output.write('{:1}\n'.format(Jupiter_units))
param_output.write('{:<50.10e}\n'.format(Bessel_function_exit))

param_output.write('{:1}\n'.format(phase_shift))
param_output.write('{:1}\n'.format(subtract_JD_mid_time))
param_output.write('{:1}\n'.format(add_RV_zero_offset))
param_output.write('{:<10.5f}\n'.format(Number_orbits))
param_output.write('{:<50.10f}\n'.format(data_plot_model_interval))
param_output.write('{:<50.5f}\n'.format(Phase_angle_start))
param_output.write('{:<50.10f}\n'.format(Time_bef_aft))
param_output.write('{:<50.10f}\n'.format(Time_compare_vel))

param_output.write('{:>150}\n'.format(my_data_plotting_symbols))
param_output.write('{:>150}\n'.format(other_data_plotting_symbols))

param_output.write('{:<50.5f}\n'.format(Ms_solar_prior))
param_output.write('{:<50.5f}\n'.format(Ms_solar_end))
param_output.write('{:<50.5f}\n'.format(Ms_solar_begin))
param_output.write('{:<50.5f}\n'.format(Ms_solar_1sigerr))
param_output.write('{:1}\n'.format(Ms_solar_norm_prior))
param_output.write('{:1}\n'.format(Ms_solar_fixed))

param_output.write('{:<50.5f}\n'.format(Rs_solar_prior))
param_output.write('{:<50.5f}\n'.format(Rs_solar_end))
param_output.write('{:<50.5f}\n'.format(Rs_solar_begin))
param_output.write('{:<50.5f}\n'.format(Rs_solar_1sigerr))
param_output.write('{:1}\n'.format(Rs_solar_norm_prior))
param_output.write('{:1}\n'.format(Rs_solar_fixed))

param_output.write('{:<50.5f}\n'.format(Rp_Rs_ratio_prior))
param_output.write('{:<50.5f}\n'.format(Rp_Rs_ratio_end))
param_output.write('{:<50.5f}\n'.format(Rp_Rs_ratio_begin))
param_output.write('{:<50.5f}\n'.format(Rp_Rs_ratio_1sigerr))
param_output.write('{:1}\n'.format(Rp_Rs_ratio_norm_prior))
param_output.write('{:1}\n'.format(Rp_Rs_ratio_fixed))
param_output.write('{:1}\n'.format(use_Rp_Rs_ratio))

param_output.write('{:<50.5f}\n'.format(vmacro_prior))
param_output.write('{:<50.5f}\n'.format(vmacro_end))
param_output.write('{:<50.5f}\n'.format(vmacro_begin))
param_output.write('{:<50.5f}\n'.format(vmacro_1sigerr))
param_output.write('{:1}\n'.format(vmacro_norm_prior))
param_output.write('{:1}\n'.format(vmacro_fixed))
param_output.write('{:<50.5f}\n'.format(beta_rot))

param_output.write('{:<20d}\n'.format(M_pixels))

param_output.write('{:1}\n'.format(linear_quadratic))
param_output.write('{:<50.5f}\n'.format(u_linear))
param_output.write('{:<50.5f}\n'.format(q_1_prior))
param_output.write('{:<50.5f}\n'.format(q_1_end))
param_output.write('{:<50.5f}\n'.format(q_1_begin))
param_output.write('{:<50.5f}\n'.format(q_1_1sigerr))
param_output.write('{:<50.5f}\n'.format(q_2_prior))
param_output.write('{:<50.5f}\n'.format(q_2_end))
param_output.write('{:<50.5f}\n'.format(q_2_begin))
param_output.write('{:<50.5f}\n'.format(q_2_1sigerr))
param_output.write('{:1}\n'.format(q_norm_prior))
param_output.write('{:1}\n'.format(q_fixed))

param_output.write('{:<50.5f}\n'.format(Inc_prior))
param_output.write('{:<50.5f}\n'.format(Inc_end))
param_output.write('{:<50.5f}\n'.format(Inc_begin))
param_output.write('{:<50.5f}\n'.format(Inc_1sigerr))
param_output.write('{:1}\n'.format(Inc_norm_prior))
param_output.write('{:1}\n'.format(Inc_fixed))

param_output.write('{:<50.5f}\n'.format(omega_arg_periastron_prior))
param_output.write('{:<50.5f}\n'.format(omega_arg_periastron_end))
param_output.write('{:<50.5f}\n'.format(omega_arg_periastron_begin))
param_output.write('{:<50.5f}\n'.format(arg_periastron_1sigerr))
param_output.write('{:1}\n'.format(arg_periastron_norm_prior))
param_output.write('{:1}\n'.format(arg_periastron_fixed))

param_output.write('{:<50.5f}\n'.format(Ecc_prior))
param_output.write('{:<50.5f}\n'.format(Ecc_end))
param_output.write('{:<50.5f}\n'.format(Ecc_begin))
param_output.write('{:<50.5f}\n'.format(Ecc_1sigerr))
param_output.write('{:1}\n'.format(Ecc_norm_prior))
param_output.write('{:1}\n'.format(Ecc_fixed))

param_output.write('{:<50.5f}\n'.format(Mp_prior))
param_output.write('{:<50.5f}\n'.format(Mp_end))
param_output.write('{:<50.5f}\n'.format(Mp_begin))
param_output.write('{:<50.5f}\n'.format(Mp_1sigerr))
param_output.write('{:1}\n'.format(Mp_norm_prior))
param_output.write('{:1}\n'.format(Mp_fixed))

param_output.write('{:<50.5f}\n'.format(Rp_prior))
param_output.write('{:<50.5f}\n'.format(Rp_end))
param_output.write('{:<50.5f}\n'.format(Rp_begin))
param_output.write('{:<50.5f}\n'.format(Rp_1sigerr))
param_output.write('{:1}\n'.format(Rp_norm_prior))
param_output.write('{:1}\n'.format(Rp_fixed))

param_output.write('{:<50.10f}\n'.format(Orbital_period_prior))
param_output.write('{:<50.10f}\n'.format(Orbital_period_end))
param_output.write('{:<50.10f}\n'.format(Orbital_period_begin))
param_output.write('{:<50.10f}\n'.format(Orbital_period_1sigerr))
param_output.write('{:1}\n'.format(Orbital_period_norm_prior))
param_output.write('{:1}\n'.format(Orbital_period_fixed))

param_output.write('{:<50.10f}\n'.format(JD_time_mid_transit_prior))
param_output.write('{:<50.10f}\n'.format(JD_time_mid_transit_end))
param_output.write('{:<50.10f}\n'.format(JD_time_mid_transit_begin))
param_output.write('{:<50.10f}\n'.format(JD_time_mid_transit_1sigerr))
param_output.write('{:1}\n'.format(JD_time_mid_transit_norm_prior))
param_output.write('{:1}\n'.format(JD_time_mid_transit_fixed))

param_output.write('{:<50.5f}\n'.format(Albedo))

param_output.write('{:<50.5f}\n'.format(RV_zero_offset_prior))
param_output.write('{:<50.5f}\n'.format(RV_zero_offset_end))
param_output.write('{:<50.5f}\n'.format(RV_zero_offset_begin))
param_output.write('{:<50.5f}\n'.format(RV_zero_offset_interval))
param_output.write('{:<50.5f}\n'.format(RV_zero_offset_1sigerr))
param_output.write('{:1}\n'.format(RV_zero_offset_norm_prior))
param_output.write('{:1}\n'.format(RV_zero_offset_fixed))

param_output.write('{:<50.5f}\n'.format(RV_offset_datasets_new_prior))
param_output.write('{:<50.5f}\n'.format(RV_offset_datasets_new_end))
param_output.write('{:<50.5f}\n'.format(RV_offset_datasets_new_begin))
param_output.write('{:<50.5f}\n'.format(RV_offset_datasets_new_interval))
param_output.write('{:<50.5f}\n'.format(RV_offset_datasets_new_1sigerr))
param_output.write('{:1}\n'.format(RV_offset_datasets_norm_prior))
param_output.write('{:1}\n'.format(RV_offset_datasets_fixed))

param_output.write('{:<50.5f}\n'.format(K_amp_prior))
param_output.write('{:<50.5f}\n'.format(K_amp_end))
param_output.write('{:<50.5f}\n'.format(K_amp_begin))
param_output.write('{:<50.5f}\n'.format(K_amp_1sigerr))
param_output.write('{:1}\n'.format(K_amp_norm_prior))
param_output.write('{:1}\n'.format(K_amp_fixed))
param_output.write('{:1}\n'.format(use_K_amp))

param_output.write('{:<50.5f}\n'.format(Impact_prior))
param_output.write('{:<50.5f}\n'.format(Impact_end))
param_output.write('{:<50.5f}\n'.format(Impact_begin))
param_output.write('{:<50.5f}\n'.format(Impact_1sigerr))
param_output.write('{:1}\n'.format(Impact_norm_prior))
param_output.write('{:1}\n'.format(Impact_fixed))
param_output.write('{:1}\n'.format(use_Impact))

param_output.write('{:<50.5f}\n'.format(Rorb_Rs_prior))
param_output.write('{:<50.5f}\n'.format(Rorb_Rs_end))
param_output.write('{:<50.5f}\n'.format(Rorb_Rs_begin))
param_output.write('{:<50.5f}\n'.format(Rorb_Rs_1sigerr))
param_output.write('{:1}\n'.format(Rorb_Rs_norm_prior))
param_output.write('{:1}\n'.format(Rorb_Rs_fixed))
param_output.write('{:1}\n'.format(use_Rorb_Rs_ratio))

param_output.write('{:<50.5f}\n'.format(vsini_prior))
param_output.write('{:<50.5f}\n'.format(vsini_end))
param_output.write('{:<50.5f}\n'.format(vsini_begin))
param_output.write('{:<50.5f}\n'.format(vsini_1sigerr))
param_output.write('{:1}\n'.format(vsini_norm_prior))
param_output.write('{:1}\n'.format(impose_prior_vsini))
param_output.write('{:1}\n'.format(vsini_fixed))

param_output.write('{:<50.5f}\n'.format(stellar_rotation_angle_prior))
param_output.write('{:<50.5f}\n'.format(stellar_rotation_angle_end))
param_output.write('{:<50.5f}\n'.format(stellar_rotation_angle_begin))
param_output.write('{:<50.5f}\n'.format(stellar_rotation_angle_1sigerr))
param_output.write('{:1}\n'.format(stellar_rotation_angle_norm_prior))
param_output.write('{:1}\n'.format(impose_prior_stellar_rotation_angle))
param_output.write('{:1}\n'.format(stellar_rotation_angle_fixed))

param_output.write('{:<20d}\n'.format(mcmc_accepted_iteration_size))
param_output.write('{:<20d}\n'.format(number_mcmc_walkers))
param_output.write('{:<50.5f}\n'.format(scale_factor))
param_output.write('{:<50.5f}\n'.format(chi_squared_change))
param_output.write('{:<50.5f}\n'.format(chi_squared_change_fit))
param_output.write('{:1}\n'.format(use_out_transit_rv_for_fit))

param_output.close()




#------------------------------------------Compile and run MCMC Fortran program----------------------------------------#
cmd1 = "gfortran -c -ffree-line-length-none -std=f2008 -O3 ExOSAM_RM_effect_V03.f08 init_seed.f08 numz.f08 ran_mod.f08"
process1 = subprocess.call(cmd1, shell=True, stdout=subprocess.PIPE)

cmd2 = "gfortran ExOSAM_RM_effect_V03.o init_seed.o numz.o ran_mod.o"
process2 = subprocess.call(cmd2, shell=True, stdout=subprocess.PIPE)

cmd3 = "./a.out"
process3 = subprocess.Popen(cmd3, shell=True, stderr=subprocess.PIPE)

while True:
    out = process3.stderr.read(1)
    if out == b'' and process3.poll() != None:
        break
    #endif
    if out != b'':
        sys.stdout.write(out.decode(sys.stdout.encoding))
        sys.stdout.flush()
    #endif
#EndWhileTrue

process3.wait()




#-----------------------------------------Read in the results from the MCMC--------------------------------------------#
input_temp_data = output_temp_data

mask_array_input = input_temp_data + 'mask_array.txt'
best_parameters_mcmc_input = input_temp_data + 'best_parameters_mcmc.txt'
best_parameters_mean_input = input_temp_data + 'best_param_mean.txt'
stand_dev_input = input_temp_data + 'stand_dev.txt'

vsini_mcmc_array_input = input_temp_data + 'vsini_array_mcmc.txt'
spin_orbit_mcmc_array_input = input_temp_data + 'stellar_rotation_angle_array_mcmc.txt'
omega_arg_periastron_mcmc_array_input = input_temp_data + 'omega_arg_periastron_MCMC_array.txt'
Inc_mcmc_array_input = input_temp_data + 'Inc_MCMC_array.txt'
Ecc_mcmc_array_input = input_temp_data + 'Ecc_MCMC_array.txt'
Rp_Rs_ratio_mcmc_array_input = input_temp_data + 'Rp_Rs_ratio_MCMC_array.txt'
Rs_solar_mcmc_array_input = input_temp_data + 'Rs_solar_MCMC_array.txt'
Ms_solar_mcmc_array_input = input_temp_data + 'Ms_solar_MCMC_array.txt'
Rp_mcmc_array_input = input_temp_data + 'Rp_MCMC_array.txt'
Mp_mcmc_array_input = input_temp_data + 'Mp_MCMC_array.txt'
JD_time_mid_transit_mcmc_array_input = input_temp_data + 'JD_time_mid_transit_MCMC_array.txt'
Orbital_period_mcmc_array_input = input_temp_data + 'Orbital_period_MCMC_array.txt'
RV_zero_offset_mcmc_array_input = input_temp_data + 'RV_zero_offset_MCMC_array.txt'
RV_offset_datasets_mcmc_array_input = input_temp_data + 'RV_offset_datasets_MCMC_array.txt'
vmacro_mcmc_array_input = input_temp_data + 'vmacro_MCMC_array.txt'
q_1_mcmc_array_input = input_temp_data + 'q_1_MCMC_array.txt'
q_2_mcmc_array_input = input_temp_data + 'q_2_MCMC_array.txt'
K_amp_mcmc_array_input = input_temp_data + 'K_amp_MCMC_array.txt'
Impact_mcmc_array_input = input_temp_data + 'Impact_MCMC_array.txt'
Rorb_Rs_mcmc_array_input = input_temp_data + 'Rorb_Rs_MCMC_array.txt'

Chi_squared_name_array = input_temp_data + 'chi_squared_array_mcmc.txt'
red_Chi_squared_name_array = input_temp_data + 'reduced_chi_squared_array_mcmc.txt'
likelihood_array_name = input_temp_data + 'likelihood_array_mcmc.txt'

MCMC_properties_input = input_temp_data + 'MCMC_properties.txt'




index_mask_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
columns_mask_array = np.empty(number_mcmc_walkers, dtype='U8')
for i in range(number_mcmc_walkers):
    columns_mask_array[i] = 'walker' + str(i)
#endfor
mask_array = pd.DataFrame(index=index_mask_array, columns=columns_mask_array, dtype='bool')
mask_array_data = np.genfromtxt(mask_array_input, dtype='bool')
length_mask_array_data = len(mask_array_data)
for col in range(number_mcmc_walkers):
    mask_array['walker' + str(col)] = mask_array_data[:,col]
#endfor

#Create a vsini array that combines all the vsini values along a single axis using the mask array.
nonmasked_elements = np.sum(mask_array)
total_nonmasked_elements = np.sum(np.sum(mask_array))
index_mcmc_nomask_array = np.arange(int(total_nonmasked_elements))
#The vsini values that are not masked. Masked values are those beyond which the mcmc iteration reached.
vsini_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['vsini'], dtype='float')
index_vsini_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
columns_for_mcmc_arrays = columns_mask_array
vsini_mcmc_array = pd.DataFrame(index=index_vsini_mcmc_array, columns=columns_for_mcmc_arrays, dtype='float')
vsini_mcmc_data = np.genfromtxt(vsini_mcmc_array_input, dtype='float')
length_vsini_array_data = vsini_mcmc_data.size

spin_orbit_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['angle'], dtype='float')
index_spin_orbit_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
spin_orbit_mcmc_array = pd.DataFrame(index=index_spin_orbit_mcmc_array, columns=columns_for_mcmc_arrays, dtype='float')
spin_orbit_mcmc_data = np.genfromtxt(spin_orbit_mcmc_array_input, dtype='float')
length_spin_orbit_array_data = spin_orbit_mcmc_data.size

omega_arg_periastron_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['angle'], dtype='float')
index_omega_arg_periastron_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
omega_arg_periastron_mcmc_array = pd.DataFrame(index=index_omega_arg_periastron_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
omega_arg_periastron_mcmc_data = np.genfromtxt(omega_arg_periastron_mcmc_array_input, dtype='float')
length_omega_arg_periastron_array_data = len(omega_arg_periastron_mcmc_data)

Inc_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['angle'], dtype='float')
index_Inc_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Inc_mcmc_array = pd.DataFrame(index=index_Inc_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Inc_mcmc_data = np.genfromtxt(Inc_mcmc_array_input, dtype='float')
length_Inc_array_data = len(Inc_mcmc_data)

Ecc_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['value'], dtype='float')
index_Ecc_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Ecc_mcmc_array = pd.DataFrame(index=index_Ecc_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Ecc_mcmc_data = np.genfromtxt(Ecc_mcmc_array_input, dtype='float')
length_Ecc_array_data = len(Ecc_mcmc_data)

Rp_Rs_ratio_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['ratio'], dtype='float')
index_Rp_Rs_ratio_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Rp_Rs_ratio_mcmc_array = pd.DataFrame(index=index_Rp_Rs_ratio_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Rp_Rs_ratio_mcmc_data = np.genfromtxt(Rp_Rs_ratio_mcmc_array_input, dtype='float')
length_Rp_Rs_ratio_array_data = len(Rp_Rs_ratio_mcmc_data)

Rs_solar_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['Rs'], dtype='float')
index_Rs_solar_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Rs_solar_mcmc_array = pd.DataFrame(index=index_Rs_solar_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Rs_solar_mcmc_data = np.genfromtxt(Rs_solar_mcmc_array_input, dtype='float')
length_Rs_solar_array_data = len(Rs_solar_mcmc_data)

Ms_solar_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['Ms'], dtype='float')
index_Ms_solar_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Ms_solar_mcmc_array = pd.DataFrame(index=index_Ms_solar_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Ms_solar_mcmc_data = np.genfromtxt(Ms_solar_mcmc_array_input, dtype='float')
length_Ms_solar_array_data = len(Ms_solar_mcmc_data)

Rp_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['Rp'], dtype='float')
index_Rp_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Rp_mcmc_array = pd.DataFrame(index=index_Rp_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Rp_mcmc_data = np.genfromtxt(Rp_mcmc_array_input, dtype='float')
length_Rp_array_data = len(Rp_mcmc_data)

Mp_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['Mp'], dtype='float')
index_Mp_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Mp_mcmc_array = pd.DataFrame(index=index_Mp_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Mp_mcmc_data = np.genfromtxt(Mp_mcmc_array_input, dtype='float')
length_Mp_array_data = len(Mp_mcmc_data)

JD_time_mid_transit_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['Time'], dtype='float')
index_JD_time_mid_transit_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
JD_time_mid_transit_mcmc_array = pd.DataFrame(index=index_JD_time_mid_transit_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
JD_time_mid_transit_mcmc_data = np.genfromtxt(JD_time_mid_transit_mcmc_array_input, dtype='float')
length_JD_time_mid_transit_array_data = len(JD_time_mid_transit_mcmc_data)

Orbital_period_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['Time'], dtype='float')
index_Orbital_period_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Orbital_period_mcmc_array = pd.DataFrame(index=index_Orbital_period_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Orbital_period_mcmc_data = np.genfromtxt(Orbital_period_mcmc_array_input, dtype='float')
length_Orbital_period_array_data = len(Orbital_period_mcmc_data)

RV_zero_offset_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['RV'], dtype='float')
index_RV_zero_offset_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
RV_zero_offset_mcmc_array = pd.DataFrame(index=index_RV_zero_offset_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
RV_zero_offset_mcmc_data = np.genfromtxt(RV_zero_offset_mcmc_array_input, dtype='float')
length_RV_zero_offset_array_data = len(RV_zero_offset_mcmc_data)

RV_offset_datasets_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['RV'], dtype='float')
index_RV_offset_datasets_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
RV_offset_datasets_mcmc_array = pd.DataFrame(index=index_RV_offset_datasets_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
RV_offset_datasets_mcmc_data = np.genfromtxt(RV_offset_datasets_mcmc_array_input, dtype='float')
length_RV_offset_datasets_array_data = len(RV_offset_datasets_mcmc_data)

vmacro_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['RV'], dtype='float')
index_vmacro_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
vmacro_mcmc_array = pd.DataFrame(index=index_vmacro_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
vmacro_mcmc_data = np.genfromtxt(vmacro_mcmc_array_input, dtype='float')
length_vmacro_array_data = len(vmacro_mcmc_data)

q_1_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['value'], dtype='float')
index_q_1_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
q_1_mcmc_array = pd.DataFrame(index=index_q_1_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
q_1_mcmc_data = np.genfromtxt(q_1_mcmc_array_input, dtype='float')
length_q_1_array_data = len(q_1_mcmc_data)

q_2_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['value'], dtype='float')
index_q_2_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
q_2_mcmc_array = pd.DataFrame(index=index_q_2_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
q_2_mcmc_data = np.genfromtxt(q_2_mcmc_array_input, dtype='float')
length_q_2_array_data = len(q_2_mcmc_data)

K_amp_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['RV'], dtype='float')
index_K_amp_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
K_amp_mcmc_array = pd.DataFrame(index=index_K_amp_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
K_amp_mcmc_data = np.genfromtxt(K_amp_mcmc_array_input, dtype='float')
length_K_amp_array_data = len(K_amp_mcmc_data)

Impact_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['value'], dtype='float')
index_Impact_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Impact_mcmc_array = pd.DataFrame(index=index_Impact_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Impact_mcmc_data = np.genfromtxt(Impact_mcmc_array_input, dtype='float')
length_Impact_array_data = len(Impact_mcmc_data)

Rorb_Rs_mcmc_nomask_array = pd.DataFrame(index=index_mcmc_nomask_array, columns=['value'], dtype='float')
index_Rorb_Rs_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
Rorb_Rs_mcmc_array = pd.DataFrame(index=index_Rorb_Rs_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
Rorb_Rs_mcmc_data = np.genfromtxt(Rorb_Rs_mcmc_array_input, dtype='float')
length_Rorb_Rs_array_data = len(Rorb_Rs_mcmc_data)





begin_element = 0
for col in range(number_mcmc_walkers):
    vsini_mcmc_array['walker' + str(col)] = vsini_mcmc_data[:,col]
    spin_orbit_mcmc_array['walker' + str(col)] = spin_orbit_mcmc_data[:, col]
    omega_arg_periastron_mcmc_array['walker' + str(col)] = omega_arg_periastron_mcmc_data[:,col]
    Inc_mcmc_array['walker' + str(col)] = Inc_mcmc_data[:,col]
    Ecc_mcmc_array['walker' + str(col)] = Ecc_mcmc_data[:,col]
    Rp_Rs_ratio_mcmc_array['walker' + str(col)] = Rp_Rs_ratio_mcmc_data[:,col]
    Rs_solar_mcmc_array['walker' + str(col)] = Rs_solar_mcmc_data[:,col]
    Ms_solar_mcmc_array['walker' + str(col)] = Ms_solar_mcmc_data[:,col]
    Rp_mcmc_array['walker' + str(col)] = Rp_mcmc_data[:,col]
    Mp_mcmc_array['walker' + str(col)] = Mp_mcmc_data[:,col]
    JD_time_mid_transit_mcmc_array['walker' + str(col)] = JD_time_mid_transit_mcmc_data[:,col]
    Orbital_period_mcmc_array['walker' + str(col)] = Orbital_period_mcmc_data[:,col]
    RV_zero_offset_mcmc_array['walker' + str(col)] = RV_zero_offset_mcmc_data[:,col]
    RV_offset_datasets_mcmc_array['walker' + str(col)] = RV_offset_datasets_mcmc_data[:,col]
    vmacro_mcmc_array['walker' + str(col)] = vmacro_mcmc_data[:, col]
    q_1_mcmc_array['walker' + str(col)] = q_1_mcmc_data[:, col]
    q_2_mcmc_array['walker' + str(col)] = q_2_mcmc_data[:, col]
    K_amp_mcmc_array['walker' + str(col)] = K_amp_mcmc_data[:, col]
    Impact_mcmc_array['walker' + str(col)] = Impact_mcmc_data[:, col]
    Rorb_Rs_mcmc_array['walker' + str(col)] = Rorb_Rs_mcmc_data[:, col]

    if col == 0:
        vsini_mcmc_nomask_array.loc[0:nonmasked_elements[col]-1, 'vsini'] = \
            vsini_mcmc_data[0:nonmasked_elements[col], col]
        spin_orbit_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'angle'] = \
            spin_orbit_mcmc_data[0:nonmasked_elements[col], col]
        omega_arg_periastron_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'angle'] = \
            omega_arg_periastron_mcmc_data[0:nonmasked_elements[col], col]
        Inc_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'angle'] = \
            Inc_mcmc_data[0:nonmasked_elements[col], col]
        Ecc_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'value'] = \
            Ecc_mcmc_data[0:nonmasked_elements[col], col]
        Rp_Rs_ratio_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'ratio'] = \
            Rp_Rs_ratio_mcmc_data[0:nonmasked_elements[col], col]
        Rs_solar_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'Rs'] = \
            Rs_solar_mcmc_data[0:nonmasked_elements[col], col]
        Ms_solar_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'Ms'] = \
            Ms_solar_mcmc_data[0:nonmasked_elements[col], col]
        Rp_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'Rp'] = \
            Rp_mcmc_data[0:nonmasked_elements[col], col]
        Mp_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'Mp'] = \
            Rp_mcmc_data[0:nonmasked_elements[col], col]
        JD_time_mid_transit_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'Time'] = \
            JD_time_mid_transit_mcmc_data[0:nonmasked_elements[col], col]
        Orbital_period_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'Time'] = \
            Orbital_period_mcmc_data[0:nonmasked_elements[col], col]
        RV_zero_offset_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'RV'] = \
            RV_zero_offset_mcmc_data[0:nonmasked_elements[col], col]
        RV_offset_datasets_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'RV'] = \
            RV_offset_datasets_mcmc_data[0:nonmasked_elements[col], col]
        vmacro_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'RV'] = \
            vmacro_mcmc_data[0:nonmasked_elements[col], col]
        q_1_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'value'] = \
            q_1_mcmc_data[0:nonmasked_elements[col], col]
        q_2_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'value'] = \
            q_2_mcmc_data[0:nonmasked_elements[col], col]
        K_amp_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'RV'] = \
            K_amp_mcmc_data[0:nonmasked_elements[col], col]
        Impact_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'value'] = \
            Impact_mcmc_data[0:nonmasked_elements[col], col]
        Rorb_Rs_mcmc_nomask_array.loc[0:nonmasked_elements[col] - 1, 'value'] = \
            Rorb_Rs_mcmc_data[0:nonmasked_elements[col], col]
    else:
        vsini_mcmc_nomask_array.loc[begin_element:begin_element+nonmasked_elements[col]-1, 'vsini'] = \
            vsini_mcmc_data[0:nonmasked_elements[col], col]
        spin_orbit_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'angle'] = \
            spin_orbit_mcmc_data[0:nonmasked_elements[col], col]
        omega_arg_periastron_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'angle'] \
            = omega_arg_periastron_mcmc_data[0:nonmasked_elements[col], col]
        Inc_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'angle'] = \
            Inc_mcmc_data[0:nonmasked_elements[col], col]
        Ecc_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'value'] = \
            Ecc_mcmc_data[0:nonmasked_elements[col], col]
        Rp_Rs_ratio_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'ratio'] = \
            Rp_Rs_ratio_mcmc_data[0:nonmasked_elements[col], col]
        Rs_solar_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'Rs'] = \
            Rs_solar_mcmc_data[0:nonmasked_elements[col], col]
        Ms_solar_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'Ms'] = \
            Ms_solar_mcmc_data[0:nonmasked_elements[col], col]
        Rp_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'Rp'] = \
            Rp_mcmc_data[0:nonmasked_elements[col], col]
        Mp_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'Mp'] = \
            Mp_mcmc_data[0:nonmasked_elements[col], col]
        JD_time_mid_transit_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'Time'] = \
            JD_time_mid_transit_mcmc_data[0:nonmasked_elements[col], col]
        Orbital_period_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'Time'] = \
            Orbital_period_mcmc_data[0:nonmasked_elements[col], col]
        RV_zero_offset_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'RV'] = \
            RV_zero_offset_mcmc_data[0:nonmasked_elements[col], col]
        RV_offset_datasets_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'RV'] = \
            RV_offset_datasets_mcmc_data[0:nonmasked_elements[col], col]
        vmacro_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'RV'] = \
            vmacro_mcmc_data[0:nonmasked_elements[col], col]
        q_1_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'value'] = \
            q_1_mcmc_data[0:nonmasked_elements[col], col]
        q_2_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'value'] = \
            q_2_mcmc_data[0:nonmasked_elements[col], col]
        K_amp_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'RV'] = \
            K_amp_mcmc_data[0:nonmasked_elements[col], col]
        Impact_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'value'] = \
            Impact_mcmc_data[0:nonmasked_elements[col], col]
        Rorb_Rs_mcmc_nomask_array.loc[begin_element:begin_element + nonmasked_elements[col] - 1, 'value'] = \
            Rorb_Rs_mcmc_data[0:nonmasked_elements[col], col]
    #endelse
    begin_element = begin_element + nonmasked_elements[col]
#endfor




#Read in the Chi squared array, best fit parameters and their respected standard deviations for each MCMC walker, and
#the overall best fit parameters and standard deviations.
index_chi_squared_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
chi_squared_mcmc_array = pd.DataFrame(index=index_chi_squared_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
chi_squared_mcmc_data = np.genfromtxt(Chi_squared_name_array, dtype='float')
length_chi_squared_array_data = len(chi_squared_mcmc_data)
for col in range(number_mcmc_walkers):
    chi_squared_mcmc_array['walker' + str(col)] = chi_squared_mcmc_data[:,col]
#endfor

index_likelihood_mcmc_array = np.arange(int(mcmc_accepted_iteration_size/0.05))
likelihood_mcmc_array = pd.DataFrame(index=index_likelihood_mcmc_array,
                                               columns=columns_for_mcmc_arrays, dtype='float')
likelihood_mcmc_data = np.genfromtxt(likelihood_array_name, dtype='float')
length_likelihood_array_data = len(likelihood_mcmc_data)
for col in range(number_mcmc_walkers):
    likelihood_mcmc_array['walker' + str(col)] = likelihood_mcmc_data[:,col]
#endfor

index_MCMC_properties_array = np.arange(int(number_mcmc_walkers))
reject_counter_array = pd.DataFrame(index=index_MCMC_properties_array, columns=["Counts"], dtype='int')
total_MCMC_iterations_array = pd.DataFrame(index=index_MCMC_properties_array, columns=["Counts"], dtype='int')
total_accepted_proposals_array = pd.DataFrame(index=index_MCMC_properties_array, columns=["Counts"], dtype='int')
acceptance_rate_array = pd.DataFrame(index=index_MCMC_properties_array, columns=["Rate"], dtype='float')
scale_factor_array = pd.DataFrame(index=index_MCMC_properties_array, columns=["Factor"], dtype='float')
MCMC_properties_data = np.genfromtxt(MCMC_properties_input, dtype='float')
reject_counter_array["Counts"] = MCMC_properties_data[0]
total_MCMC_iterations_array["Counts"] = MCMC_properties_data[1]
total_accepted_proposals_array["Counts"] = MCMC_properties_data[2]
acceptance_rate_array["Rate"] = MCMC_properties_data[3]
scale_factor_array["Factor"] = MCMC_properties_data[4]




#Read in the best over parameters values.
best_overall_parameters = open(best_parameters_mcmc_input, 'r')

format_read1 = ff.FortranRecordReader('(92X, F15.6)')
line = best_overall_parameters.readline()
min_chi_squared_total = format_read1.read(line)[0]
format_read2 = ff.FortranRecordReader('(92X, I10, 1X, I10)')
format_read3 = ff.FortranRecordReader('(92X, I10)')
line = best_overall_parameters.readline()
loc_min_chi_squared_total = format_read2.read(line)
line = best_overall_parameters.readline()
min_r_chi_squared_total = format_read1.read(line)[0]
line = best_overall_parameters.readline()
loc_min_r_chi_squared_total = format_read2.read(line)
line = best_overall_parameters.readline()
max_likelihood_total = format_read1.read(line)[0]
line = best_overall_parameters.readline()
loc_max_likelihood_total = format_read2.read(line)

line = best_overall_parameters.readline()
best_vsini_min_chi_squared = format_read1.read(line)[0]
line = best_overall_parameters.readline()
best_vsini_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
vsini_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_spin_orbit_min_chi_squared = format_read1.read(line)[0]
line = best_overall_parameters.readline()
best_spin_orbit_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
spin_orbit_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_RV_offset_datasets_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
RV_offset_datasets_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_RV_zero_offset_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
RV_zero_offset_stand_dev_mcmc_mean = format_read1.read(line)[0]

format_read3 = ff.FortranRecordReader('(92X, F20.10)')
line = best_overall_parameters.readline()
best_orbital_period_mcmc_mean = format_read3.read(line)[0]
line = best_overall_parameters.readline()
orbital_period_stand_dev_mcmc_mean = format_read3.read(line)[0]

line = best_overall_parameters.readline()
best_JD_time_mid_transit_mcmc_mean = format_read3.read(line)[0]
line = best_overall_parameters.readline()
JD_time_mid_transit_stand_dev_mcmc_mean = format_read3.read(line)[0]

line = best_overall_parameters.readline()
best_Mp_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
Mp_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_Rp_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
Rp_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_Ms_solar_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
Ms_solar_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_Rs_solar_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
Rs_solar_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_Rp_Rs_ratio_mcmc_mean = format_read3.read(line)[0]
line = best_overall_parameters.readline()
Rp_Rs_ratio_stand_dev_mcmc_mean = format_read3.read(line)[0]

line = best_overall_parameters.readline()
best_Ecc_mcmc_mean = format_read3.read(line)[0]
line = best_overall_parameters.readline()
Ecc_stand_dev_mcmc_mean = format_read3.read(line)[0]

line = best_overall_parameters.readline()
best_Inc_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
Inc_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_omega_arg_periastron_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
omega_arg_periastron_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_vmacro_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
vmacro_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_q_1_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
q_1_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_q_2_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
q_2_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_K_amp_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
K_amp_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_Impact_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
Impact_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_Rorb_Rs_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
Rorb_Rs_stand_dev_mcmc_mean = format_read1.read(line)[0]

format_read3 = ff.FortranRecordReader('(92X, I10)')
line = best_overall_parameters.readline()
Number_fit = format_read3.read(line)[0]
line = best_overall_parameters.readline()
Number_RV_points = format_read3.read(line)[0]

best_overall_parameters.close()




#Read in the best parameters for each MCMC walker.
best_walker_parameters = open(best_parameters_mean_input, 'r')

if number_mcmc_walkers == 1:
    format_read_walkers = ff.FortranRecordReader('(F16.8)')
else:
    format_read_walkers = ff.FortranRecordReader(str(number_mcmc_walkers) + '(F16.8, 1X)')
#endelse
line = best_walker_parameters.readline()
RV_offset_datasets_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
RV_zero_offset_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
orbital_period_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
JD_time_mid_transit_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Mp_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Rp_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Ms_solar_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Rs_solar_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Rp_Rs_ratio_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Ecc_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Inc_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
omega_arg_periastron_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
vsini_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
spin_orbit_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
vmacro_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
q_1_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
q_2_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
K_amp_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Impact_mean_walkers = format_read_walkers.read(line)

line = best_walker_parameters.readline()
Rorb_Rs_mean_walkers = format_read_walkers.read(line)

best_walker_parameters.close()




#Read in the standard deviations of each parameter determined for each MCMC walker.
stand_dev_walker_parameters = open(stand_dev_input, 'r')

line = stand_dev_walker_parameters.readline()
RV_offset_datasets_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
RV_zero_offset_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
orbital_period_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
JD_time_mid_transit_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Mp_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Rp_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Ms_solar_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Rs_solar_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Rp_Rs_ratio_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Ecc_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Inc_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
omega_arg_periastron_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
vsini_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
spin_orbit_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
vmacro_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
q_1_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
q_2_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
K_amp_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Impact_stand_dev_walkers = format_read_walkers.read(line)

line = stand_dev_walker_parameters.readline()
Rorb_Rs_stand_dev_walkers = format_read_walkers.read(line)

stand_dev_walker_parameters.close()




#Determine the median of the posterior distributions. Will output median best fit values as well as the mean values.
spin_orbit_median_all_walkers = np.median(spin_orbit_mcmc_nomask_array)
vsini_median_all_walkers = np.median(vsini_mcmc_nomask_array)




test_null_hypothesis = False
#Create a model RV, RM, and transit array with the best fit parameters.
Transit_LC_array, RM_effect_array, RV_array, Timestep, Time_transit_actual, Time_occultation_actual, \
Time_transit_90inc, min_LC_value, amp_RM_effect, max_frac_flux_decrease, RV_full_amp, max_flux_planet, \
Transit_start_no_inc, Time_mid_transit, Transit_end_no_inc, Time_occultation_start, Time_occultation_end = \
    RV_model_array(data_plot_model_interval, Number_orbits, Bessel_function_exit, Jupiter_units, use_Rp_Rs_ratio,
                   best_vsini_mcmc_mean, best_spin_orbit_mcmc_mean, Orbital_period_prior, JD_time_mid_transit_prior,
                   Mp_prior, Mp_fixed, Rp_prior, Ms_solar_prior, Rs_solar_prior, Rp_Rs_ratio_prior, Ecc_prior, Inc_prior,
                   omega_arg_periastron_prior, beta_rot, vmacro_prior, linear_quadratic, q_1_prior, q_2_prior, u_linear,
                   K_amp_prior, use_K_amp, Impact_prior, use_Impact, Rorb_Rs_prior, use_Rorb_Rs_ratio,
                   test_null_hypothesis, Albedo, M_pixels)




#-------------------------------------Compile and run RV residuals Fortran program-------------------------------------#
cmd1 = "gfortran -c -ffree-line-length-none -std=f2008 -O3 RM_effect_residuals_V03.f08"
process1 = subprocess.call(cmd1, shell=True, stdout=subprocess.PIPE)

cmd2 = "gfortran RM_effect_residuals_V03.o"
process2 = subprocess.call(cmd2, shell=True, stdout=subprocess.PIPE)

cmd3 = "./a.out"
process3 = subprocess.Popen(cmd3, shell=True, stderr=subprocess.PIPE)

while True:
    out = process3.stderr.read(1)
    if out == b'' and process3.poll() != None:
        break
    #endif
    if out != b'':
        sys.stdout.write(out.decode(sys.stdout.encoding))
        sys.stdout.flush()
    #endif
#EndWhileTrue

process3.wait()




#------------------------------------------Read in the RV residual arrays----------------------------------------------#
residual_array_input = input_temp_data + 'residual_array.txt'
data_fit_array_input = input_temp_data + 'data_fit_array.txt'
RV_theory_fit_array_input = input_temp_data + 'RV_Theory_fit_array.txt'
residual_all_array_input = input_temp_data + 'residual_all_array.txt'
RV_theory_all_array_input = input_temp_data + 'RV_Theory_all_array.txt'

if other_RV_files == 'Y':
    total_data_points = datafilelength_other1 + my_datafilelength
else:
    total_data_points = my_datafilelength
#endelse

#Read in residual array and differentiate which ones are from my RV's and which ones are from other RV's.
columns_for_residual_arrays = ["Time", "RV", "RV_error"]
residual_data = np.genfromtxt(residual_array_input, dtype='float')
length_residual_data = len(residual_data)
index_residual = np.arange(length_residual_data)
residual_array = pd.DataFrame(index=index_residual, columns=columns_for_residual_arrays, dtype='float')
residual_array['Time'] = residual_data[:,0]
residual_array['RV'] = residual_data[:,1]
residual_array['RV_error'] = residual_data[:,2]

residual_all_data = np.genfromtxt(residual_all_array_input, dtype='float')
length_residual_all_data = len(residual_all_data)
index_residual_all = np.arange(length_residual_all_data)
residuals_all_array = pd.DataFrame(index=index_residual_all, columns=columns_for_residual_arrays, dtype='float')
residuals_all_array['Time'] = residual_all_data[:,0]
residuals_all_array['RV'] = residual_all_data[:,1]
residuals_all_array['RV_error'] = residual_all_data[:,2]

chi_squared_fit_input = input_temp_data + 'chi_squared_value.txt'
#Read in the chi squared and model fit parameters from residuals.
chi_squared_fit = open(chi_squared_fit_input, 'r')
format_read = ff.FortranRecordReader('(92X, F15.6)')
line = chi_squared_fit.readline()
model_data_difference_residuals = format_read.read(line)[0]
line = chi_squared_fit.readline()
chi_2_residuals = format_read.read(line)[0]
line = chi_squared_fit.readline()
r_Chi_2_residuals = format_read.read(line)[0]
line = chi_squared_fit.readline()
likelihood_residuals = format_read.read(line)[0]
chi_squared_fit.close()

#BIC_residuals = chi_2_residuals + (Number_fit)*math.log(Number_RV_points)
BIC_residuals = chi_2_residuals + 3.0*math.log(Number_RV_points)




#---------------------------Test the null hypothesis and determine residuals--------------------------------------------
if test_null == 'Y':
    test_null_hypothesis = True
    # Create a model RV, RM, and transit array with the best fit parameters.
    Transit_LC_array_null, RM_effect_array_null, RV_array_null, Timestep_null, Time_transit_actual_null, \
    Time_occultation_actual_null, Time_transit_90inc_null, min_LC_value_null, amp_RM_effect_null, \
    max_frac_flux_decrease_null, RV_full_amp_null, max_flux_planet_null, Transit_start_no_inc_null, \
    Time_mid_transit_null, Transit_end_no_inc_null, Time_occultation_start_null, Time_occultation_end_null = \
        RV_model_array(data_plot_model_interval, Number_orbits, Bessel_function_exit, Jupiter_units, use_Rp_Rs_ratio,
                       best_vsini_mcmc_mean, best_spin_orbit_mcmc_mean, Orbital_period_prior, JD_time_mid_transit_prior,
                       Mp_prior, Mp_fixed, Rp_prior, Ms_solar_prior, Rs_solar_prior, Rp_Rs_ratio_prior, Ecc_prior,
                       Inc_prior, omega_arg_periastron_prior, beta_rot, vmacro_prior, linear_quadratic, q_1_prior,
                       q_2_prior, u_linear, K_amp_prior, use_K_amp, Impact_prior, use_Impact, Rorb_Rs_prior,
                       use_Rorb_Rs_ratio, test_null_hypothesis, Albedo, M_pixels)

    # Set test_null_hypothesis to false first.
    test_null_hypothesis_text = 'T'
    test_null_file_output = open(output_temp_data + 'test_null.txt', 'w')
    test_null_file_output.write('{:1}\n'.format(test_null_hypothesis_text))
    test_null_file_output.close()

    #-------------------------------------Compile and run RV residuals Fortran program---------------------------------#
    cmd1 = "gfortran -c -ffree-line-length-none -std=f2008 -O3 RM_effect_residuals_V03.f08"
    process1 = subprocess.call(cmd1, shell=True, stdout=subprocess.PIPE)

    cmd2 = "gfortran RM_effect_residuals_V03.o"
    process2 = subprocess.call(cmd2, shell=True, stdout=subprocess.PIPE)

    cmd3 = "./a.out"
    process3 = subprocess.Popen(cmd3, shell=True, stderr=subprocess.PIPE)

    while True:
        out = process3.stderr.read(1)
        if out == b'' and process3.poll() != None:
            break
        #endif
        if out != b'':
            sys.stdout.write(out.decode(sys.stdout.encoding))
            sys.stdout.flush()
        #endif
    #EndWhileTrue

    process3.wait()

    #------------------------------------------Read in the RV residual arrays------------------------------------------#

    # Read in residual array and differentiate which ones are from my RV's and which ones are from other RV's.
    residual_data_null = np.genfromtxt(residual_array_input, dtype='float')
    length_residual_data_null = len(residual_data_null)
    index_residual_null = np.arange(length_residual_data_null)
    residual_array_null = pd.DataFrame(index=index_residual_null, columns=columns_for_residual_arrays, dtype='float')
    residual_array_null['Time'] = residual_data_null[:, 0]
    residual_array_null['RV'] = residual_data_null[:, 1]
    residual_array_null['RV_error'] = residual_data_null[:, 2]

    residual_all_data_null = np.genfromtxt(residual_all_array_input, dtype='float')
    length_residual_all_data_null = len(residual_all_data_null)
    index_residual_all_null = np.arange(length_residual_all_data_null)
    residuals_all_array_null = pd.DataFrame(index=index_residual_all_null, columns=columns_for_residual_arrays,
                                            dtype='float')
    residuals_all_array_null['Time'] = residual_all_data_null[:, 0]
    residuals_all_array_null['RV'] = residual_all_data_null[:, 1]
    residuals_all_array_null['RV_error'] = residual_all_data_null[:, 2]

    chi_squared_fit_input = input_temp_data + 'chi_squared_value.txt'
    #Read in the chi squared and model fit parameters from residuals.
    chi_squared_fit = open(chi_squared_fit_input, 'r')
    format_read = ff.FortranRecordReader('(92X, F15.6)')
    line = chi_squared_fit.readline()
    model_data_difference_null = format_read.read(line)[0]
    line = chi_squared_fit.readline()
    chi_2_null = format_read.read(line)[0]
    line = chi_squared_fit.readline()
    r_Chi_2_null = format_read.read(line)[0]
    line = chi_squared_fit.readline()
    likelihood_null = format_read.read(line)[0]
    chi_squared_fit.close()

    #BIC_null = chi_2_null + (Number_fit - 1.0)*math.log(Number_RV_points)
    BIC_null = chi_2_null + 1.0 * math.log(Number_RV_points)
    BIC_difference = BIC_null - BIC_residuals




#Make sure RV and residual times matches those in the model.
#Rearrange my RV data points so that they are in chronological order with time.
total_interval = round((Orbital_period*Number_orbits)/data_plot_model_interval)
for t in range(my_datafilelength):
    if template_RV.loc[t, "Time"] * (24.0 * 3600.0) < RV_array.loc[0, 'Time']:
        template_RV.loc[t, "Time"] = (RV_array.loc[total_interval - 1, 'Time']/(24.0*3600.0))\
                                     - abs(template_RV.loc[t, 'Time'] - RV_array.loc[0, 'Time']/(24.0*3600.0))
    else:
        template_RV.loc[t, 'Time'] = template_RV.loc[t, 'Time']
    #endelse
#endfor

#Now change the order (if necessary).
template_RV = template_RV.sort_values('Time', ascending=True)
template_RV = template_RV.reset_index(drop=True)

#rearrange other the RV data points so that they are in chronological order with time.
#First change time to match one orbital phase.
if other_RV_files == 'Y':
    for t in range(datafilelength_other1):
        if template_other1.loc[t, "Time"] * (24.0 * 3600.0) < RV_array.loc[0, 'Time']:
            template_other1.loc[t, "Time"] = (RV_array.loc[total_interval - 1, 'Time']/(24.0*3600.0)) - \
                                             abs(template_other1.loc[t, 'Time'] - RV_array.loc[0, 'Time']/(24.0*3600.0))
        else:
            template_other1.loc[t, 'Time'] = template_other1.loc[t, 'Time']
        #endelse
    #endfor

    template_other1 = template_other1.sort_values('Time', ascending=True)
    template_other1 = template_other1.reset_index(drop=True)
#endif

#Rearrange residual RV data points so that they are in chronological order with time.
#First change time to match one orbital phase.
for t in range(total_data_points):
    if residuals_all_array.loc[t, "Time"] * (24.0 * 3600.0) < RV_array.loc[0, 'Time']:
        residuals_all_array.loc[t, "Time"] = (RV_array.loc[total_interval - 1, 'Time']/(24.0*3600.0)) - \
                                              abs(residuals_all_array.loc[t, 'Time'] -
                                                  RV_array.loc[0, 'Time']/(24.0*3600.0))
        if test_null == 'Y':
            residuals_all_array_null.loc[t, "Time"] = (RV_array.loc[total_interval - 1, 'Time'] / (24.0 * 3600.0)) - \
                                                 abs(residuals_all_array_null.loc[t, 'Time'] -
                                                     RV_array.loc[0, 'Time'] / (24.0 * 3600.0))

    else:
        residuals_all_array.loc[t, 'Time'] = residuals_all_array.loc[t, 'Time']
        if test_null == 'Y':
            residuals_all_array_null.loc[t, 'Time'] = residuals_all_array_null.loc[t, 'Time']
    #endelse
#endfor

residuals_all_array = residuals_all_array.sort_values('Time', ascending=True)
residuals_all_array = residuals_all_array.reset_index(drop=True)

if test_null == 'Y':
    residuals_all_array_null = residuals_all_array_null.sort_values('Time', ascending=True)
    residuals_all_array_null = residuals_all_array_null.reset_index(drop=True)
#endif




index_my_RV_residuals = np.arange(my_datafilelength)
my_RV_residuals_array = pd.DataFrame(index=index_my_RV_residuals, columns=columns_for_residual_arrays, dtype='float')
if test_null == 'Y':
    my_RV_residuals_null_array = pd.DataFrame(index=index_my_RV_residuals, columns=columns_for_residual_arrays,
                                         dtype='float')
residual_symbol_mydata = pd.DataFrame(index=index_my_RV_residuals, columns=["Symbols"], dtype='float')

my_rv_num_insert = 0
#Now differentiate the residuals as either from my data or from other data.
for i in range(total_data_points):
    for j in range(my_datafilelength):
        if abs(template_RV.loc[j, 'Time'] - residuals_all_array.loc[i, 'Time']) <= 0.00001:
            my_RV_residuals_array.loc[my_rv_num_insert, 'Time'] = residuals_all_array.loc[i, 'Time']
            my_RV_residuals_array.loc[my_rv_num_insert, 'RV'] = residuals_all_array.loc[i, 'RV']
            my_RV_residuals_array.loc[my_rv_num_insert, 'RV_error'] = residuals_all_array.loc[i, 'RV_error']

            if test_null == 'Y':
                my_RV_residuals_null_array.loc[my_rv_num_insert, 'Time'] = residuals_all_array_null.loc[i, 'Time']
                my_RV_residuals_null_array.loc[my_rv_num_insert, 'RV'] = residuals_all_array_null.loc[i, 'RV']
                my_RV_residuals_null_array.loc[my_rv_num_insert, 'RV_error'] = residuals_all_array_null.loc[i, 'RV_error']

            residual_symbol_mydata.loc[my_rv_num_insert, 'Symbols'] = my_data_plotting_symbols_array.loc[j, 'Symbols']

            #print("my_rv_num_insert: ", my_rv_num_insert)
            #print("my_RV_residuals_array.loc[my_rv_num_insert, 'Time']: ",
            #      my_RV_residuals_array.loc[my_rv_num_insert, 'Time'])
            #print("Time difference: ", abs(template_RV.loc[j, 'Time'] - residuals_all_array.loc[i, 'Time']))
            my_rv_num_insert = my_rv_num_insert + 1

        #endif
    #endfor
#endfor

if other_RV_files == 'Y':
    other_rv_num_insert = 0
    index_other_RV_residuals = np.arange(datafilelength_other1)
    other_RV_residuals_array = pd.DataFrame(index=index_other_RV_residuals, columns=columns_for_residual_arrays,
                                         dtype='float')
    if test_null == 'Y':
        other_RV_residuals_null_array = pd.DataFrame(index=index_other_RV_residuals, columns=columns_for_residual_arrays,
                                                dtype='float')

    residual_symbol_otherdata = pd.DataFrame(index=index_other_RV_residuals, columns=["Symbols"], dtype='float')

    #Now differentiate the residuals as either from my data or from other dataset.
    for i in range(total_data_points):
        for j in range(datafilelength_other1):

            if abs(template_other1.loc[j, 'Time'] - residuals_all_array.loc[i, 'Time']) <= 0.0001:
                other_RV_residuals_array.loc[other_rv_num_insert, 'Time'] = residuals_all_array.loc[i, 'Time']
                other_RV_residuals_array.loc[other_rv_num_insert, 'RV'] = residuals_all_array.loc[i, 'RV']
                other_RV_residuals_array.loc[other_rv_num_insert, 'RV_error'] = residuals_all_array.loc[i, 'RV_error']

                if test_null == 'Y':
                    other_RV_residuals_null_array.loc[other_rv_num_insert, 'Time'] = residuals_all_array_null.loc[i,
                                                                                                                'Time']
                    other_RV_residuals_null_array.loc[other_rv_num_insert, 'RV'] = residuals_all_array_null.loc[i, 'RV']
                    other_RV_residuals_null_array.loc[other_rv_num_insert, 'RV_error'] = residuals_all_array_null.loc[
                        i, 'RV_error']

                residual_symbol_otherdata.loc[other_rv_num_insert, 'Symbols'] = other_data_plotting_symbols_array.loc[
                    j, 'Symbols']

                print("other_rv_num_insert: ", other_rv_num_insert)
                print("other_RV_residuals_array.loc[other_rv_num_insert, 'Time']: ",
                      other_RV_residuals_array.loc[other_rv_num_insert, 'Time'])
                print("Time difference: ", abs(template_other1.loc[j, 'Time'] - residuals_all_array.loc[i, 'Time']))

                other_rv_num_insert = other_rv_num_insert + 1
            #endif
        #endfor
    #endfor
#endif



#Apply velocity offsets to RV data.
RV_my_array_column = ["Time", "RV", "RV_error"]
RV_my_array_index = np.arange(my_datafilelength)
RV_my_array = pd.DataFrame(index=RV_my_array_index, columns=RV_my_array_column)

if other_RV_files == 'Y':
    RV_other_array_column = ["Time", "RV", "RV_error"]
    RV_other_array_index = np.arange(datafilelength_other1)
    RV_other_array = pd.DataFrame(index=RV_other_array_index, columns=RV_other_array_column)
#endif

RV_my_array.loc[:, 'Time'] = template_RV.loc[:, 'Time']
RV_my_array.loc[:, 'RV'] = template_RV.loc[:, 'RV'] + (best_RV_zero_offset_mcmc_mean - RV_zero_offset_prior) \
                           + RV_offset_datasets_new_prior
RV_my_array.loc[:, 'RV_error'] = template_RV.loc[:, 'RV_error']

if other_RV_files == 'Y':
    RV_other_array.loc[:, 'Time'] = template_other1.loc[:, 'Time']
    RV_other_array.loc[:, 'RV'] = template_other1.loc[:, 'RV'] + (best_RV_zero_offset_mcmc_mean - RV_zero_offset_prior)
    RV_other_array.loc[:, 'RV_error'] = template_other1.loc[:, 'RV_error']
#endif




#Create a residual time array with the time in seconds before/after mid transit instead of days.
total_interval = round((Orbital_period*Number_orbits)/data_plot_model_interval)
my_res_time = pd.DataFrame(index=index_my_RV_residuals, columns=["Time"], dtype='float')
for t in range(my_datafilelength):
    if my_RV_residuals_array.loc[t, "Time"]*(24.0*3600.0) < RV_array.loc[0, 'Time']:
        my_res_time.loc[t, 'Time'] = RV_array.loc[total_interval - 1, 'Time'] - abs(my_RV_residuals_array.loc[t, 'Time']
                                                                            *(24.0*3600.0) - RV_array.loc[0, 'Time'])
    else:
        my_res_time.loc[t, 'Time'] = my_RV_residuals_array.loc[t, 'Time']*(24.0*3600.0)
    #endelse
#endfor

index_all_residuals = np.arange(total_data_points)
res_time_all = pd.DataFrame(index=index_all_residuals, columns=["Time"], dtype='float')
for t in range(total_data_points):
    if residuals_all_array.loc[t, "Time"]*(24.0*3600.0) < RV_array.loc[0, 'Time']:
        res_time_all.loc[t, 'Time'] = RV_array.loc[total_interval - 1, 'Time'] - abs(residuals_all_array.loc[t, 'Time']
                                                                            *(24.0*3600.0) - RV_array.loc[0, 'Time'])
    else:
        res_time_all.loc[t, 'Time'] = residuals_all_array.loc[t, 'Time']*(24.0*3600.0)
    #endelse
#endfor

if other_RV_files == 'Y':
    other_res_time = pd.DataFrame(index=index_other_RV_residuals, columns=["Time"], dtype='float')
    for t in range(datafilelength_other1):
        if other_RV_residuals_array.loc[t, "Time"] * (24.0 * 3600.0) < RV_array.loc[0, 'Time']:
            other_res_time.loc[t, 'Time'] = RV_array.loc[total_interval - 1, 'Time'] - abs(
                other_RV_residuals_array.loc[t, 'Time']*(24.0 * 3600.0) - RV_array.loc[0, 'Time'])
        else:
            other_res_time.loc[t, 'Time'] = other_RV_residuals_array.loc[t, 'Time'] * (24.0 * 3600.0)
        #endelse
    #endfor
#endif




#Plot a histogram of the residuals to verify if Gaussian.
RV_residuals_fit_histo = save_data_directory + 'RV_residuals_fit_histo_' + Planet_name + '.pdf'
RV_residuals_abs_fit_histo = save_data_directory + 'RV_residuals_abs_fit_histo_' + Planet_name + '.pdf'
residuals_offset_array = pd.DataFrame(index=index_my_RV_residuals, columns=["RV"], dtype='float')
residuals_offset_abs_array = pd.DataFrame(index=index_my_RV_residuals, columns=["RV"], dtype='float')

for residual_loop in range(my_datafilelength):
    residuals_offset_array.loc[residual_loop, 'RV'] = my_RV_residuals_array.loc[residual_loop, 'RV'] \
                                                          / my_RV_residuals_array.loc[residual_loop, 'RV_error']
    residuals_offset_abs_array.loc[residual_loop, 'RV'] = abs(my_RV_residuals_array.loc[residual_loop, 'RV']) \
                                                      /my_RV_residuals_array.loc[residual_loop, 'RV_error']
#endfor

#Best Gaussian fit to the residuals.
(mu_resids, sigma_resids) = norm.fit(residuals_offset_array)

dx_resids = knuth_bin_width(residuals_offset_array['RV'], return_bins=False, quiet=True)
num_bins_resids = int(np.ceil((max(residuals_offset_array['RV']) - min(residuals_offset_array['RV']))/dx_resids))

#Histogram of the residuals
#n_resids, bins_resids, patches_resids = plt.hist(residuals_offset_array["RV"], bins=num_bins_resids, normed=1,
#                                                 facecolor='green', alpha=0.75)
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
n_resids, bins_resids, patches_resids = plt.hist(residuals_offset_array["RV"],
                                                 int(np.ceil(len(residuals_offset_array)/2.0) + 1), normed=1,
                                                 facecolor='k', alpha=0.75)

# add a 'best fit' line
y_resids = mlab.normpdf(bins_resids, mu_resids, sigma_resids)
l_resids = plt.plot(bins_resids, y_resids, 'r--', linewidth=4)

#plot
plt.xlabel('RV Residuals', fontsize=16)
plt.ylabel('Probability', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu_resids, sigma_resids), fontsize=16)
plt.savefig(RV_residuals_fit_histo,dpi=150)
plt.close(fig)





(mu_resids_abs, sigma_resids_abs) = norm.fit(residuals_offset_abs_array)

dx_resids_abs = knuth_bin_width(residuals_offset_abs_array['RV'], return_bins=False, quiet=True)
num_bins_resids_abs = int(np.ceil((max(residuals_offset_abs_array['RV']) - min(residuals_offset_abs_array['RV']))/
                                  dx_resids_abs))

#Histogram of the residuals
#n_resids, bins_resids, patches_resids = plt.hist(residuals_offset_array["RV"], bins=num_bins_resids, normed=1,
#                                                 facecolor='green', alpha=0.75)
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
n_resids_abs, bins_resids_abs, patches_resids_abs = plt.hist(residuals_offset_abs_array["RV"],
                                                 int(np.ceil(len(residuals_offset_abs_array)/2.0) + 1), normed=1,
                                                 facecolor='k', alpha=0.75)

# add a 'best fit' line
y_resids_abs = mlab.normpdf(bins_resids_abs, mu_resids_abs, sigma_resids_abs)
l_resids_abs = plt.plot(bins_resids_abs, y_resids_abs, 'r--', linewidth=4)

#plot
plt.xlabel('RV Residuals', fontsize=16)
plt.ylabel('Probability', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu_resids_abs, sigma_resids_abs), fontsize=16)
plt.savefig(RV_residuals_abs_fit_histo,dpi=150)
plt.close(fig)




#Plot the best fit model transit lightcurves, RM-effect, and RV curve.
Transit_LC_whole_orbit = save_data_directory + 'Model_LC_whole_orbit_' + Planet_name + '.pdf'
Transit_LC_whole_orbit_zoomed = save_data_directory + 'Model_LC_whole_orbit_zoomed_' + Planet_name + '.pdf'
Transit_centered = save_data_directory + 'Model_LC_centered_transit_' + Planet_name + '.pdf'
occultation_centered = save_data_directory + 'Model_LC_centered_occultation_' + Planet_name + '.pdf'
RM_plot_out = save_data_directory + 'Model_RM_centered_transit_' + Planet_name + '.pdf'
RV_plot_out = save_data_directory + 'Model_RV_curve_' + Planet_name + '.pdf'
all_residuals_plot_out = save_data_directory + 'Residuals_all_data_' + Planet_name + '.pdf'
RVcurve_op_all_data_transit_out = save_data_directory + 'RVcurve_op_all_data_transit_' + Planet_name + '.pdf'
RVcurve_op_all_data_orbit_out = save_data_directory + 'RVcurve_op_all_data_orbit_' + Planet_name + '.pdf'
RVcurve_op_my_data_transit_out = save_data_directory + 'RVcurve_op_my_data_transit_' + Planet_name + '.pdf'
RVcurve_op_my_data_orbit_out = save_data_directory + 'RVcurve_op_my_data_orbit_' + Planet_name + '.pdf'

# RVcurve_op_all_data_transit_null_out = save_data_directory + 'RVcurve_op_all_data_transit_null_' + Planet_name + '.pdf'
# RVcurve_op_all_data_orbit_null_out = save_data_directory + 'RVcurve_op_all_data_orbit_null_' + Planet_name + '.pdf'
# RVcurve_op_my_data_transit_null_out = save_data_directory + 'RVcurve_op_my_data_transit_null_' + Planet_name + '.pdf'
# RVcurve_op_my_data_orbit_null_out = save_data_directory + 'RVcurve_op_my_data_orbit_null_' + Planet_name + '.pdf'

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.plot(Transit_LC_array["Time"]/60.0, Transit_LC_array["Flux"], linestyle='solid', linewidth=3, color="k")
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel('Relative Flux', fontsize=16)
plt.savefig(Transit_LC_whole_orbit,dpi=150)
plt.close(fig)

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.ylim(0.999 - (max(Transit_LC_array["Flux"] - 1.0)),1.0001 + (max(Transit_LC_array["Flux"] - 1.0)))
plt.plot(Transit_LC_array["Time"]/60.0, Transit_LC_array["Flux"], linestyle='solid', linewidth=3, color="k")
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel('Relative Flux', fontsize=16)
plt.savefig(Transit_LC_whole_orbit_zoomed,dpi=150)
plt.close(fig)

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.xlim((Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60.0,
         (Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60.0)
plt.plot(Transit_LC_array["Time"]/60.0, Transit_LC_array["Flux"], linestyle='solid', linewidth=3, color="k")
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel('Relative Flux', fontsize=16)
plt.savefig(Transit_centered,dpi=150)
plt.close(fig)

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.xlim((Time_occultation_start - Time_bef_aft - Time_mid_transit)/60.0,
         (Time_occultation_end + Time_bef_aft - Time_mid_transit)/60.0)
plt.ylim(0.999 - (max(Transit_LC_array["Flux"] - 1.0)),1.0001 + (max(Transit_LC_array["Flux"] - 1.0)))
plt.plot(Transit_LC_array["Time"]/60.0, Transit_LC_array["Flux"], linestyle='solid', linewidth=3, color="k")
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel('Relative Flux', fontsize=16)
plt.savefig(occultation_centered,dpi=150)
plt.close(fig)

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.xlim((Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60.0,
         (Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60.0)
plt.plot(RM_effect_array["Time"]/60.0, RM_effect_array["RV"], linestyle='solid', linewidth=3, color="k")
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
plt.savefig(RM_plot_out,dpi=150)
plt.close(fig)

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.plot(RV_array["Time"]/60.0, RV_array["RV"], linestyle='solid', linewidth=3, color="k")
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
plt.savefig(RV_plot_out, dpi=150)
plt.close(fig)

#Plot the residuals of the RVs to the best model fit.
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(res_time_all["Time"]/60.0, residuals_all_array["RV"], yerr=residuals_all_array["RV_error"], fmt='o',
             ecolor='b', color='k', elinewidth=2, markersize=10, capsize=7)
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
plt.savefig(all_residuals_plot_out, dpi=150)
plt.close(fig)




#change the IDL plotting symbol markers to python symbol markers.
for i in range(num_my_rv):
    if my_data_plotting_symbols_array.loc[i, 'Symbols'] == 16:
        #The symbol is a circle.
        my_data_plotting_symbols_array.loc[i, 'Symbols'] = 'o'
    #endif
    if my_data_plotting_symbols_array.loc[i, 'Symbols'] == 17:
        #The symbol is a triangle.
        my_data_plotting_symbols_array.loc[i, 'Symbols'] = '^'
    #endif
#endfor
#print(my_data_plotting_symbols_array)
if other_RV_files == 'Y':
    for i in range(num_rv1):
        if other_data_plotting_symbols_array.loc[i, 'Symbols'] == 36:
            #The symbol is an x.
            other_data_plotting_symbols_array.loc[i, 'Symbols'] = 'x'
        #endif
    #endfor
#endif




#Now plot the best fit model with the radial velocities.
if other_RV_files == 'Y':

    lower_xlim = (Transit_start_no_inc - Time_bef_aft - Time_mid_transit) / 60.0
    upper_xlim = (Transit_end_no_inc + Time_bef_aft - Time_mid_transit) / 60.0

    position_range_mydata = RV_my_array.where((RV_my_array["Time"]*(24.0*60.0) >= lower_xlim) &
                                              (RV_my_array["Time"]*(24.0*60.0) <= upper_xlim))
    start_index_mydata = position_range_mydata["Time"].first_valid_index()
    end_index_mydata = position_range_mydata["Time"].last_valid_index()

    position_range_otherdata = RV_other_array.where((RV_other_array["Time"] * (24.0 * 60.0) >= lower_xlim) &
                                              (RV_other_array["Time"] * (24.0 * 60.0) <= upper_xlim))
    start_index_otherdata = position_range_otherdata["Time"].first_valid_index()
    end_index_otherdata = position_range_otherdata["Time"].last_valid_index()

    ymax_value1 = np.max(position_range_mydata["RV"]+position_range_mydata["RV_error"]+1.0)
    ymax_value2 = np.max(position_range_otherdata["RV"]+position_range_otherdata["RV_error"]+1.0)
    ymax_value_RV = max(ymax_value1, ymax_value2)

    ymin_value1 = np.min(position_range_mydata["RV"]-position_range_mydata["RV_error"]-1.0)
    ymin_value2 = np.min(position_range_otherdata["RV"]-position_range_otherdata["RV_error"]-1.0)
    ymin_value_RV = min(ymin_value1, ymin_value2)




    position_range_myresid = my_RV_residuals_array.where((my_RV_residuals_array["Time"] * (24.0 * 60.0) >= lower_xlim) &
                                              (my_RV_residuals_array["Time"] * (24.0 * 60.0) <= upper_xlim))
    if test_null == 'Y':
        position_range_myresid_null = my_RV_residuals_null_array.where(
            (my_RV_residuals_null_array["Time"] * (24.0 * 60.0) >= lower_xlim) &
            (my_RV_residuals_null_array["Time"] * (24.0 * 60.0) <= upper_xlim))
    start_index_myresid = position_range_myresid["Time"].first_valid_index()
    end_index_myresid = position_range_myresid["Time"].last_valid_index()

    position_range_otherresid = other_RV_residuals_array.where((other_RV_residuals_array["Time"] * (24.0 * 60.0) >=
                                                                lower_xlim) &
                                                    (other_RV_residuals_array["Time"] * (24.0 * 60.0) <= upper_xlim))
    if test_null == 'Y':
        position_range_otherresid_null = other_RV_residuals_null_array.where((other_RV_residuals_null_array["Time"] *
                                                                              (24.0 * 60.0) >=
                                                                    lower_xlim) &
                                                                   (other_RV_residuals_null_array["Time"] * (
                                                                   24.0 * 60.0) <= upper_xlim))
    start_index_otherresid = position_range_otherresid["Time"].first_valid_index()
    end_index_otherresid = position_range_otherresid["Time"].last_valid_index()

    ymax_value1 = np.max(position_range_myresid["RV"] + position_range_myresid["RV_error"] + 1.0)
    ymax_value2 = np.max(position_range_otherresid["RV"] + position_range_otherresid["RV_error"] + 1.0)
    if test_null == 'Y':
        ymax_value3 = np.max(position_range_myresid_null["RV"] + position_range_myresid_null["RV_error"] + 1.0)
        ymax_value4 = np.max(position_range_otherresid_null["RV"] + position_range_otherresid_null["RV_error"] + 1.0)
    ymax_value_resid = max(ymax_value1, ymax_value2)
    if test_null == 'Y':
        ymax_value_resid = max(ymax_value1, ymax_value2, ymax_value3, ymax_value4)

    ymin_value1 = np.min(position_range_myresid["RV"] - position_range_myresid["RV_error"] - 1.0)
    ymin_value2 = np.min(position_range_otherresid["RV"] - position_range_otherresid["RV_error"] - 1.0)
    if test_null == 'Y':
        ymin_value3 = np.min(position_range_myresid_null["RV"] - position_range_myresid_null["RV_error"] - 1.0)
        ymin_value4 = np.min(position_range_otherresid_null["RV"] - position_range_otherresid_null["RV_error"] - 1.0)
    ymin_value_resid = min(ymin_value1, ymin_value2)
    if test_null == 'Y':
        ymin_value_resid = min(ymin_value1, ymin_value2, ymin_value3, ymin_value4)

    fig = plt.figure()
    gs = gridspec.GridSpec(100, 100)
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    ax1 = plt.subplot(gs[:65, 1:100])
    plt.xlim(lower_xlim, upper_xlim)
    plt.ylim(ymin_value_RV, ymax_value_RV + 50.0)
    plt.setp(plt.gca(), 'xticklabels', [])
    #plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14, horizontalalignment='right')
    plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
    for i in range(num_my_rv):
        ax1.errorbar(RV_my_array.loc[i, "Time"]*(24.0*60.0), RV_my_array.loc[i, "RV"],
                     yerr=RV_my_array.loc[i, "RV_error"], fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='k', color='k', mfc='k', elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0,
                     fillstyle='full', alpha=0.9, markeredgecolor='k', zorder=2)
    # endfor
    for i in range(num_rv1):
        ax1.errorbar(RV_other_array.loc[i, "Time"]*(24.0*60.0), RV_other_array.loc[i, "RV"],
                     yerr=RV_other_array.loc[i, "RV_error"], fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='b', color='k', mfc='k', elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0,
                     fillstyle='full', alpha=0.9, markeredgecolor='k', zorder=1)
    # endfor
    if test_null == 'Y':
        ax1.plot(RV_array_null["Time"] / 60.0, RV_array_null["RV"], linestyle='-.', linewidth=1.25, color="b",
                 zorder=3, alpha=0.6)
    ax1.plot(RV_array["Time"] / 60.0, RV_array["RV"], linestyle='dashed', linewidth=1.25, color="r", zorder=3)
    ax1.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on")
    ax1.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
    ax1.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_major_locator(plt.MaxNLocator('auto'))
    ax1.xaxis.set_major_locator(plt.MaxNLocator('auto'))

    ax2 = plt.subplot(gs[65:100, 1:100])
    plt.xlim(lower_xlim, upper_xlim)
    plt.ylim(ymin_value_resid, ymax_value_resid)
    plt.ylabel(r'O -- C (ms$^{-1}$)', fontsize=14)
    for i in range(num_my_rv):
        ax2.errorbar(my_RV_residuals_array.loc[i, "Time"] * (24.0 * 60.0), my_RV_residuals_array.loc[i, "RV"],
                     yerr=my_RV_residuals_array.loc[i, "RV_error"],
                     fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', mfc='k', color='k',
                     elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9,
                     markeredgecolor='k', zorder=2)
        if test_null == 'Y':
            ax2.errorbar(my_RV_residuals_null_array.loc[i, "Time"] * (24.0 * 60.0),
                         my_RV_residuals_null_array.loc[i, "RV"], yerr=my_RV_residuals_null_array.loc[i, "RV_error"],
                         fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', mfc='k', color='k',
                         elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.5,
                         markeredgecolor='k', zorder=2)
    # endfor
    for i in range(num_rv1):
        ax2.errorbar(other_RV_residuals_array.loc[i, "Time"] * (24.0 * 60.0), other_RV_residuals_array.loc[i, "RV"],
                     yerr=other_RV_residuals_array.loc[i, "RV_error"],
                     fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='b', color='k', mfc='k',
                     elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9,
                     markeredgecolor='k', zorder=1)
        if test_null == 'Y':
            ax2.errorbar(other_RV_residuals_null_array.loc[i, "Time"] * (24.0 * 60.0),
                         other_RV_residuals_null_array.loc[i, "RV"],
                         yerr=other_RV_residuals_null_array.loc[i, "RV_error"],
                         fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='b', color='k', mfc='k',
                         elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.5,
                         markeredgecolor='k', zorder=1)
    # endfor
    ax2.axhline(y=0, linestyle='dashed', linewidth=1.25, color="r", zorder=3)
    plt.xlabel('Time (Minutes From Mid Transit)', fontsize=14)
    ax2.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on")
    ax2.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
    ax2.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_major_locator(plt.MaxNLocator('auto'))
    ax2.xaxis.set_major_locator(plt.MaxNLocator('auto'))

    plt.savefig(RVcurve_op_all_data_transit_out, dpi=150)
    plt.close(fig)




    ymax_value1 = np.max(RV_my_array["RV"] + RV_my_array["RV_error"] + 10.0)
    ymax_value2 = np.max(RV_other_array["RV"] + RV_other_array["RV_error"] + 10.0)
    ymax_value_RV = max(ymax_value1, ymax_value2)

    ymin_value1 = np.min(RV_my_array["RV"] - RV_my_array["RV_error"] - 10.0)
    ymin_value2 = np.min(RV_other_array["RV"] - RV_other_array["RV_error"] - 10.0)
    ymin_value_RV = min(ymin_value1, ymin_value2)

    ymax_value1 = np.max(my_RV_residuals_array["RV"] + my_RV_residuals_array["RV_error"] + 10.0)
    ymax_value2 = np.max(other_RV_residuals_array["RV"] + other_RV_residuals_array["RV_error"] + 10.0)
    if test_null == 'Y':
        ymax_value3 = np.max(my_RV_residuals_null_array["RV"] + my_RV_residuals_null_array["RV_error"] + 10.0)
        ymax_value4 = np.max(other_RV_residuals_null_array["RV"] + other_RV_residuals_null_array["RV_error"] + 10.0)
    ymax_value_resid = max(ymax_value1, ymax_value2)
    if test_null == 'Y':
        ymax_value_resid = max(ymax_value1, ymax_value2, ymax_value3, ymax_value4)

    ymin_value1 = np.min(my_RV_residuals_array["RV"] - my_RV_residuals_array["RV_error"] - 10.0)
    ymin_value2 = np.min(other_RV_residuals_array["RV"] - other_RV_residuals_array["RV_error"] - 10.0)
    if test_null == 'Y':
        ymin_value3 = np.min(my_RV_residuals_null_array["RV"] - my_RV_residuals_null_array["RV_error"] - 10.0)
        ymin_value4 = np.min(other_RV_residuals_null_array["RV"] - other_RV_residuals_null_array["RV_error"] - 10.0)
    ymin_value_resid = min(ymin_value1, ymin_value2)
    if test_null == 'Y':
        ymin_value_resid = min(ymin_value1, ymin_value2, ymin_value3, ymin_value4)

    lower_xlim = RV_array.loc[0, "Time"]/60.0
    upper_xlim = RV_array.loc[total_interval - 1, 'Time']/60.0

    fig = plt.figure()
    gs = gridspec.GridSpec(100, 100)
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    ax1 = plt.subplot(gs[:65, 3:100])
    plt.xlim(lower_xlim, upper_xlim)
    plt.ylim(ymin_value_RV, ymax_value_RV)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
    #plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14, horizontalalignment='right')
    for i in range(num_my_rv):
        ax1.errorbar(RV_my_array.loc[i, "Time"] * (24.0 * 60.0), RV_my_array.loc[i, "RV"],
                     yerr=RV_my_array.loc[i, "RV_error"], fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='k', color='k', mfc='k', elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0,
                     fillstyle='full', alpha=0.9, markeredgecolor='k', zorder=2)
    # endfor
    for i in range(num_rv1):
        ax1.errorbar(RV_other_array.loc[i, "Time"] * (24.0 * 60.0), RV_other_array.loc[i, "RV"],
                     yerr=RV_other_array.loc[i, "RV_error"], fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='b', color='k', mfc='k', elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0,
                     fillstyle='full', alpha=0.9, markeredgecolor='k', zorder=1)
    # endfor
    if test_null == 'Y':
        ax1.plot(RV_array_null["Time"] / 60.0, RV_array_null["RV"], linestyle='-.', linewidth=1.25, color="b",
                 zorder=3, alpha=0.6)
    ax1.plot(RV_array["Time"] / 60.0, RV_array["RV"], linestyle='dashed', linewidth=1.25, color="r", zorder=3)
    ax1.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on")
    ax1.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
    ax1.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.yaxis.set_major_locator(plt.MaxNLocator('auto'))
    ax1.xaxis.set_major_locator(plt.MaxNLocator('auto'))

    ax2 = plt.subplot(gs[65:100, 3:100])
    plt.xlim(lower_xlim, upper_xlim)
    plt.ylim(ymin_value_resid, ymax_value_resid)
    plt.ylabel(r'O -- C (ms$^{-1}$)', fontsize=14)
    for i in range(num_my_rv):
        ax2.errorbar(my_RV_residuals_array.loc[i, "Time"] * (24.0 * 60.0), my_RV_residuals_array.loc[i, "RV"],
                     yerr=my_RV_residuals_array.loc[i, "RV_error"],
                     fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k',
                     elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9,
                     markeredgecolor='k', zorder=2)
        if test_null == 'Y':
            ax2.errorbar(my_RV_residuals_null_array.loc[i, "Time"] * (24.0 * 60.0),
                         my_RV_residuals_null_array.loc[i, "RV"], yerr=my_RV_residuals_null_array.loc[i, "RV_error"],
                         fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k',
                         elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.5,
                         markeredgecolor='k', zorder=2)
    # endfor
    for i in range(num_rv1):
        ax2.errorbar(other_RV_residuals_array.loc[i, "Time"] * (24.0 * 60.0), other_RV_residuals_array.loc[i, "RV"],
                     yerr=other_RV_residuals_array.loc[i, "RV_error"],
                     fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='b', color='k', mfc='k',
                     elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9,
                     markeredgecolor='k', zorder=1)
        if test_null == 'Y':
            ax2.errorbar(other_RV_residuals_null_array.loc[i, "Time"] * (24.0 * 60.0),
                         other_RV_residuals_null_array.loc[i, "RV"],
                         yerr=other_RV_residuals_null_array.loc[i, "RV_error"],
                         fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='b', color='k', mfc='k',
                         elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.5,
                         markeredgecolor='k', zorder=1)
    # endfor
    ax2.axhline(y=0, linestyle='dashed', linewidth=1.25, color="r", zorder=3)
    plt.xlabel('Time (Minutes From Mid Transit)', fontsize=14)
    ax2.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on")
    ax2.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
    ax2.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
    ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax2.yaxis.set_major_locator(plt.MaxNLocator('auto'))
    ax2.xaxis.set_major_locator(plt.MaxNLocator('auto'))

    plt.savefig(RVcurve_op_all_data_orbit_out, dpi=150)
    plt.close(fig)

#endif




#Now just plot the fit with just your RV data.
lower_xlim = (Transit_start_no_inc - Time_bef_aft - Time_mid_transit) / 60.0
upper_xlim = (Transit_end_no_inc + Time_bef_aft - Time_mid_transit) / 60.0

position_range_mydata = RV_my_array.where((RV_my_array["Time"]*(24.0*60.0) >= lower_xlim) &
                                          (RV_my_array["Time"]*(24.0*60.0) <= upper_xlim))
start_index_mydata = position_range_mydata["Time"].first_valid_index()
end_index_mydata = position_range_mydata["Time"].last_valid_index()

ymax_value_RV = np.max(position_range_mydata["RV"]+position_range_mydata["RV_error"]+1.0)

ymin_value_RV = np.min(position_range_mydata["RV"]-position_range_mydata["RV_error"]-1.0)

position_range_myresid = my_RV_residuals_array.where((my_RV_residuals_array["Time"] * (24.0 * 60.0) >= lower_xlim) &

                                                     (my_RV_residuals_array["Time"] * (24.0 * 60.0) <= upper_xlim))

start_index_myresid = position_range_myresid["Time"].first_valid_index()
end_index_myresid = position_range_myresid["Time"].last_valid_index()

ymax_value_resid = np.max(position_range_myresid["RV"] + position_range_myresid["RV_error"] + 1.0)
if test_null == 'Y':
    position_range_myresid_null = my_RV_residuals_null_array.where(
        (my_RV_residuals_null_array["Time"] * (24.0 * 60.0) >= lower_xlim) &
        (my_RV_residuals_null_array["Time"] * (24.0 * 60.0) <= upper_xlim))
    ymax_value1 = np.max(position_range_myresid["RV"] + position_range_myresid["RV_error"] + 1.0)
    ymax_value2 = np.max(position_range_myresid_null["RV"] + position_range_myresid_null["RV_error"] + 1.0)
    ymax_value_resid = max(ymax_value1, ymax_value2)

ymin_value_resid = np.min(position_range_myresid["RV"] - position_range_myresid["RV_error"] - 1.0)
if test_null == 'Y':
    ymin_value1 = np.min(position_range_myresid["RV"] - position_range_myresid["RV_error"] - 1.0)
    ymin_value2 = np.min(position_range_myresid_null["RV"] - position_range_myresid_null["RV_error"] - 1.0)
    ymin_value_resid = min(ymin_value1, ymin_value2)

fig = plt.figure()
gs = gridspec.GridSpec(100, 100)
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
ax1 = plt.subplot(gs[:65, 1:100])
plt.xlim(lower_xlim, upper_xlim)
plt.ylim(ymin_value_RV, ymax_value_RV)
plt.setp(plt.gca(), 'xticklabels', [])
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
for i in range(num_my_rv):
    ax1.errorbar(RV_my_array.loc[i, "Time"]*(24.0*60.0), RV_my_array.loc[i, "RV"], yerr=RV_my_array.loc[i, "RV_error"],
                 fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k', elinewidth=1.5,
                 markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9, markeredgecolor='k',
                 zorder=1)
#endfor
if test_null == 'Y':
    ax1.plot(RV_array_null["Time"] / 60.0, RV_array_null["RV"], linestyle='-.', linewidth=1.25, color="b",
             zorder=3, alpha=0.6)
ax1.plot(RV_array["Time"]/60.0, RV_array["RV"], linestyle='dashed', linewidth=1.25, color="r", zorder=2)
ax1.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on", length=6)
ax1.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
ax1.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_major_locator(plt.MaxNLocator('auto'))
ax1.xaxis.set_major_locator(plt.MaxNLocator('auto'))

ax2 = plt.subplot(gs[65:100, 1:100])
plt.xlim(lower_xlim, upper_xlim)
plt.ylim(ymin_value_resid, ymax_value_resid)
for i in range(num_my_rv):
    ax2.errorbar(my_RV_residuals_array.loc[i, "Time"] * (24.0 * 60.0), my_RV_residuals_array.loc[i, "RV"],
                 yerr=my_RV_residuals_array.loc[i, "RV_error"],
                 fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k', elinewidth=1.5,
                 markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9, markeredgecolor='k',
                 zorder=1)

    if test_null == 'Y':
        ax2.errorbar(my_RV_residuals_null_array.loc[i, "Time"] * (24.0 * 60.0),
                     my_RV_residuals_null_array.loc[i, "RV"], yerr=my_RV_residuals_null_array.loc[i, "RV_error"],
                     fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k',
                     elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.5,
                     markeredgecolor='k', zorder=2)
#endfor
ax2.axhline(y=0, linestyle='dashed', linewidth=1.25, color="r", zorder=2)
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=14)
ax2.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on", length=6)
ax2.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
ax2.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
plt.ylabel(r'O -- C (ms$^{-1}$)', fontsize=14)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_major_locator(plt.MaxNLocator('auto'))
ax2.xaxis.set_major_locator(plt.MaxNLocator('auto'))

plt.savefig(RVcurve_op_my_data_transit_out, dpi=150)
plt.close(fig)




ymax_value_RV = np.max(RV_my_array["RV"] + RV_my_array["RV_error"] + 1.0)

ymin_value_RV = np.min(RV_my_array["RV"] - RV_my_array["RV_error"] - 1.0)

ymax_value_resid = np.max(my_RV_residuals_array["RV"] + my_RV_residuals_array["RV_error"] + 1.0)
if test_null == 'Y':
    ymax_value1 = np.max(my_RV_residuals_array["RV"] + my_RV_residuals_array["RV_error"] + 1.0)
    ymax_value2 = np.max(my_RV_residuals_null_array["RV"] + my_RV_residuals_null_array["RV_error"] + 10.0)
    ymax_value_resid = max(ymax_value1, ymax_value2)

ymin_value_resid = np.min(my_RV_residuals_array["RV"] - my_RV_residuals_array["RV_error"] - 1.0)
if test_null == 'Y':
    ymin_value1 = np.min(my_RV_residuals_array["RV"] - my_RV_residuals_array["RV_error"] - 1.0)
    ymin_value2 = np.min(my_RV_residuals_null_array["RV"] - my_RV_residuals_null_array["RV_error"] - 10.0)
    ymin_value_resid = min(ymin_value1, ymin_value2)

lower_xlim = RV_array.loc[0, "Time"]/60.0
upper_xlim = RV_array.loc[total_interval - 1, 'Time']/60.0

fig = plt.figure()
gs = gridspec.GridSpec(100, 100)
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
ax1 = plt.subplot(gs[:65, 1:100])
plt.xlim(lower_xlim, upper_xlim)
plt.ylim(ymin_value_RV, ymax_value_RV)
plt.setp(plt.gca(), 'xticklabels', [])
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
for i in range(num_my_rv):
    ax1.errorbar(RV_my_array.loc[i, "Time"] * (24.0 * 60.0), RV_my_array.loc[i, "RV"],
                 yerr=RV_my_array.loc[i, "RV_error"],
                 fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k', elinewidth=1.5,
                 markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9, markeredgecolor='k',
                 zorder=1)
#endfor
if test_null == 'Y':
    ax1.plot(RV_array_null["Time"] / 60.0, RV_array_null["RV"], linestyle='-.', linewidth=1.25, color="b",
             zorder=3, alpha=0.6)
ax1.plot(RV_array["Time"] / 60.0, RV_array["RV"], linestyle='dashed', linewidth=1.25, color="r", zorder=2)
ax1.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on", length=6)
ax1.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
ax1.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_major_locator(plt.MaxNLocator('auto'))
ax1.xaxis.set_major_locator(plt.MaxNLocator('auto'))

ax2 = plt.subplot(gs[65:100, 1:100])
plt.xlim(lower_xlim, upper_xlim)
plt.ylim(ymin_value_resid, ymax_value_resid)
for i in range(num_my_rv):
    ax2.errorbar(my_RV_residuals_array.loc[i, "Time"] * (24.0 * 60.0), my_RV_residuals_array.loc[i, "RV"],
                 yerr=my_RV_residuals_array.loc[i, "RV_error"],
                 fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k', elinewidth=1.5,
                 markersize=5, capsize=5, markeredgewidth=1.0, fillstyle='full', alpha=0.9, markeredgecolor='k',
                 zorder=2)

    if test_null == 'Y':
        ax2.errorbar(my_RV_residuals_null_array.loc[i, "Time"] * (24.0 * 60.0),
                     my_RV_residuals_null_array.loc[i, "RV"], yerr=my_RV_residuals_null_array.loc[i, "RV_error"],
                     fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k',
                     elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.5,
                     markeredgecolor='k', zorder=2)
#endfor
ax2.axhline(y=0, linestyle='dashed', linewidth=1.25, color="r", zorder=2)
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=14)
ax2.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on", length=6)
ax2.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
ax2.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
plt.ylabel(r'O -- C (ms$^{-1}$)', fontsize=14)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_major_locator(plt.MaxNLocator('auto'))
ax2.xaxis.set_major_locator(plt.MaxNLocator('auto'))

plt.savefig(RVcurve_op_my_data_orbit_out, dpi=150)
plt.close(fig)




#Plot the posterior distribution.
posterior_out = save_data_directory + 'posterior_distro_' + Planet_name + '.pdf'

fig = plt.figure()
gs = gridspec.GridSpec(1000, 1000)
plt.rc('text', usetex=True)
plt.rc('font', family='Times New Roman')
ax1 = plt.subplot(gs[350:940, 5:650])

# Remove the last ytick label
if total_nonmasked_elements > 10000:
    idx=random.sample(range(total_nonmasked_elements),10000)
else:
    idx = np.arange(total_nonmasked_elements)
#endelse
ax1.plot(spin_orbit_mcmc_nomask_array.loc[idx,'angle'], vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0, ".", color="k",
         alpha=0.20, zorder=1, markersize=2)
#fig.show()
fig.canvas.draw()
ax1.set_xticklabels(ax1.get_xticks())
ax1.set_yticklabels(ax1.get_yticks())
labels_y = [tick.get_text() for tick in ax1.get_yticklabels()]
labels_x = [tick.get_text() for tick in ax1.get_xticklabels()]
ax1.set_yticklabels(labels_y[:-1])
ax1.set_xticklabels(labels_x[:-1])
fig.canvas.draw()
plt.xlabel(r'$\lambda$ ($^{\circ}$)', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12, direction='in', top="on", right="on")
plt.ylabel(r'$v$sin$i_{\star}$ (kms$^{-1}$)', fontsize=14)

#Best Gaussian to the spin-orbit distribution.
(mu_lambda, sigma_lambda) = norm.fit(spin_orbit_mcmc_nomask_array.loc[idx,'angle'])

dx_lambda = knuth_bin_width(spin_orbit_mcmc_nomask_array.loc[idx,'angle'], return_bins=False, quiet=True)
num_bins_lambda = int(np.ceil((max(spin_orbit_mcmc_nomask_array.loc[idx,'angle']) -
                               min(spin_orbit_mcmc_nomask_array.loc[idx,'angle']))/dx_lambda))

#Best Gaussian to the vsini distribution.
(mu_vsini, sigma_vsini) = norm.fit(vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0)

dx_vsini = knuth_bin_width(vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0, return_bins=False, quiet=True)
num_bins_vsini = int(np.ceil((max(vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0) -
                               min(vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0))/dx_vsini))

#Histogram of the spin-orbit distribution.
ax2 = plt.subplot(gs[0:349, 5:650])
plt.setp(plt.gca(), 'xticklabels', [])
n_lambda, bins_lambda, patches_lambda = ax2.hist(spin_orbit_mcmc_nomask_array.loc[idx,'angle'], bins=num_bins_lambda,
                                                 density=True, facecolor='k', alpha=0.75, cumulative=False,
                                                 stacked=False)

# add a 'best fit' line
y_lambda = mlab.normpdf(bins_lambda, mu_lambda, sigma_lambda)
l_lambda = ax2.plot(bins_lambda, y_lambda, 'r--', linewidth=2)

plt.ylabel('Probability', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=12, direction='in', top="on", right="on")

#Histogram of the vsini distribution.
ax3 = plt.subplot(gs[350:940, 651:1000])
plt.setp(plt.gca(), 'yticklabels', [])
n_vsini, bins_vsini, patches_vsini = ax3.hist(vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0, bins=num_bins_vsini,
                                              density=True, facecolor='k', alpha=0.75, cumulative=False,
                                              orientation='horizontal', stacked=False)

# add a 'best fit' line
y_vsini = mlab.normpdf(bins_vsini, mu_vsini, sigma_vsini)
l_vsini = ax3.plot(y_vsini, bins_vsini, 'r--', linewidth=2)

#plot
plt.tick_params(axis='both', which='major', labelsize=12, direction='in', top=True, right=True)
plt.setp(ax3.get_xticklabels(), rotation=45, horizontalalignment='right')
plt.xlabel('Probability', fontsize=14)

#Contour the scatter plot
nbins_x = min(len(bins_lambda), len(bins_vsini))
nbins_y = min(len(bins_vsini), len(bins_lambda))

sigma = 0.1 # this depends on how noisy your data is, play with it!

spin_orbit_mcmc_nomask_array_filt = gaussian_filter(spin_orbit_mcmc_nomask_array.loc[idx,'angle'], sigma)
vsini_mcmc_nomask_array_filt = gaussian_filter(vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0, sigma)

H, xedges, yedges = np.histogram2d(spin_orbit_mcmc_nomask_array_filt,
                                   vsini_mcmc_nomask_array_filt, bins=(nbins_x,nbins_y),
                                   density=True)

#H, xedges, yedges = np.histogram2d(spin_orbit_mcmc_nomask_array.loc[idx,'angle'],
#                                   vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0, bins=(nbins_x,nbins_y),
#                                   normed=True)
x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
pdf = (H*(x_bin_sizes*y_bin_sizes))
one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
levels = sorted([one_sigma, two_sigma, three_sigma])
X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
Z = pdf.T
contour = ax1.contour(X, Y, Z, levels=levels, origin="lower", zorder=2, colors=['r', 'y', 'blue'])

plt.savefig(posterior_out,dpi=150)
plt.close(fig)




posterior_corr_out = save_data_directory + 'posterior_corr_' + Planet_name + '.pdf'

#data = np.vstack([spin_orbit_mcmc_nomask_array.loc[idx,'angle'], vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0,
#                  Inc_mcmc_nomask_array.loc[idx, 'angle'], Rp_Rs_ratio_mcmc_nomask_array.loc[idx, 'ratio'],
#                  Mp_mcmc_nomask_array.loc[idx, 'Mp'], Orbital_period_mcmc_nomask_array.loc[idx, 'Time'],
#                  JD_time_mid_transit_mcmc_nomask_array.loc[idx, 'Time']-2450000.00,
#                  Ecc_mcmc_nomask_array.loc[idx, 'value'], omega_arg_periastron_mcmc_nomask_array.loc[idx, 'angle']])

data = np.vstack([spin_orbit_mcmc_nomask_array.loc[idx,'angle'], vsini_mcmc_nomask_array.loc[idx,'vsini']/1000.0])
labels = np.array([r"$\lambda$ ($^{\circ}$)", r"$v$sin$i_{\star}$ (kms$^{-1}$)"])
if Inc_fixed == 'N':
    data = np.vstack([data, Inc_mcmc_nomask_array.loc[idx, 'angle']])
    labels1 = np.array([r"$i$ ($^{\circ}$)"])
    labels = np.append(labels, labels1)
#endif

if Impact_fixed == 'N':
    data = np.vstack([data, Impact_mcmc_nomask_array.loc[idx, 'value']])
    labels1 = np.array([r"$b$"])
    labels = np.append(labels, labels1)
#endif

if Rorb_Rs_fixed == 'N':
    data = np.vstack([data, Rorb_Rs_mcmc_nomask_array.loc[idx, 'value']])
    labels1 = np.array([r"$a/R_{\star}$"])
    labels = np.append(labels, labels1)
#endif

if Rp_Rs_ratio_fixed == 'N':
    data = np.vstack([data, Rp_Rs_ratio_mcmc_nomask_array.loc[idx, 'ratio']])
    labels1 = np.array([r"$R_{p}/R_{\star}$"])
    labels = np.append(labels, labels1)
#endif

if Mp_fixed == 'N':
    data = np.vstack([data, Mp_mcmc_nomask_array.loc[idx, 'Mp']])
    labels1 = np.array([r"$M_{p}$ $(M_{J})$"])
    labels = np.append(labels, labels1)
#endif

if Orbital_period_fixed == 'N':
    data = np.vstack([data, Orbital_period_mcmc_nomask_array.loc[idx, 'Time']])
    labels1 = np.array([r"$Period$ (days)"])
    labels = np.append(labels, labels1)
#endif

if JD_time_mid_transit_fixed == 'N':
    data = np.vstack([data, JD_time_mid_transit_mcmc_nomask_array.loc[idx, 'Time']-2450000.00])
    labels1 = np.array([r"$T_{\circ}$ (Days-2450000.00)"])
    labels = np.append(labels, labels1)
#endif

if Ecc_fixed == 'N':
    data = np.vstack([data, Ecc_mcmc_nomask_array.loc[idx, 'value']])
    labels1 = np.array([r"$e$"])
    labels = np.append(labels, labels1)
#endif

if arg_periastron_fixed == 'N':
    data = np.vstack([data, omega_arg_periastron_mcmc_nomask_array.loc[idx, 'angle']])
    labels1 = np.array([r"$\varpi$"])
    labels = np.append(labels, labels1)
#endif

if K_amp_fixed == 'N':
    data = np.vstack([data, K_amp_mcmc_nomask_array.loc[idx, 'RV']])
    labels1 = np.array([r"$K$ km\,s$^{-1}$"])
    labels = np.append(labels, labels1)
#endif

if RV_offset_datasets_fixed == 'N':
    data = np.vstack([data, RV_offset_datasets_mcmc_nomask_array.loc[idx, 'RV']])
    labels1 = np.array([r"$V_{d}$ m\,s$^{-1}$"])
    labels = np.append(labels, labels1)
#endif

if RV_zero_offset_fixed == 'N':
    data = np.vstack([data, RV_zero_offset_mcmc_nomask_array.loc[idx, 'RV']])
    labels1 = np.array([r"$V_{0}$ m\,s$^{-1}$"])
    labels = np.append(labels, labels1)
#endif

data2 = data.T

# Plot it.
figure = corner.corner(data2, labels=labels, quantiles=[0.05, 0.32, 0.50, 0.68, 0.95], show_titles=True,
                       title_kwargs={"fontsize": 12})
figure.savefig(posterior_corr_out,dpi=150)
close(figure)




#Output the final results in a nicely formatted text file.
results_output = open(save_data_directory + 'results_out.txt', 'w')
header = "***********************************************Best Fit Model Properties*************************************\
***********"
results_output.write('{:>120}\n'.format(header))
results_output.write('\n')
results_output.write('{:>80}'.format('Full Amplitude of the Rossiter-McLaughlin Effect Velocity Anomaly (m/s): '))
results_output.write('{:<20.1f}\n'.format(amp_RM_effect))
results_output.write('{:>80}'.format('Full Amplitude of the Keplerian Orbit (m/s): '))
results_output.write('{:<20.1f}\n'.format(RV_full_amp))
results_output.write('{:>80}'.format('Fractional Transit Depth: '))
results_output.write('{:<20.6f}\n'.format(1.0 - min_LC_value))
results_output.write('{:>80}'.format('Fractional Occultation Depth: '))
results_output.write('{:<20.6f}\n'.format(max_flux_planet))
results_output.write('{:>80}'.format('Transit Length (s): '))
results_output.write('{:<20.0f}\n'.format(Time_transit_actual))
results_output.write('{:>80}'.format('Occultation Length (s): '))
results_output.write('{:<20.0f}\n'.format(Time_occultation_actual))
results_output.write('\n')
results_output.write('\n')
header = "****************************************************MCMC Properties******************************************\
***********"
results_output.write('{:>120}\n'.format(header))
results_output.write('\n')
results_output.write('{:>80}'.format('Minimum Chi Squared Value: '))
results_output.write('{:<20.2f}\n'.format(min_chi_squared_total))
results_output.write('{:>80}'.format('Minimum Reduced Chi Squared Value: '))
results_output.write('{:<20.2f}\n'.format(min_r_chi_squared_total))
results_output.write('{:>80}'.format('Maximum Likelihood: '))
results_output.write('{:<20.4f}\n'.format(max_likelihood_total))

if test_null == 'Y':
    results_output.write('{:>80}'.format('Chi Squared Value of Null Hypothesis: '))
    results_output.write('{:<20.2f}\n'.format(chi_2_null))
    results_output.write('{:>80}'.format('Reduced Chi Squared Value of Null Hypothesis: '))
    results_output.write('{:<20.2f}\n'.format(r_Chi_2_null))
    results_output.write('{:>80}'.format('Likelihood of Null Hypothesis: '))
    results_output.write('{:<20.4f}\n'.format(likelihood_null))
    results_output.write('{:>80}'.format('BIC of Residuals: '))
    results_output.write('{:<20.4f}\n'.format(BIC_residuals))
    results_output.write('{:>80}'.format('BIC of Residuals for Null Hypothesis: '))
    results_output.write('{:<20.4f}\n'.format(BIC_null))
    results_output.write('{:>80}'.format('Difference in BIC between RM and Null models: '))
    results_output.write('{:<20.4f}\n'.format(BIC_difference))

results_output.write('{:>80}'.format('Total Number of Accepted MCMC Iterations: '))
results_output.write('{:<20d}\n'.format(int(np.sum(total_accepted_proposals_array)[0])))
results_output.write('{:>80}'.format('Total MCMC Iterations: '))
results_output.write('{:<20d}\n'.format(int(np.sum(total_MCMC_iterations_array)[0])))
results_output.write('{:>80}'.format('Total Rejected MCMC Proposals: '))
results_output.write('{:<20d}\n'.format(int(np.sum(reject_counter_array)[0])))
results_output.write('{:>80}'.format('Average Acceptance Rate of MCMC Proposals: '))
results_output.write('{:<20.6f}\n'.format(np.mean(acceptance_rate_array)[0]))
results_output.write('{:>80}'.format('Average Scale Factor MCMC Proposals: '))
results_output.write('{:<20.6f}\n'.format(np.mean(scale_factor_array)[0]))
results_output.write('\n')
results_output.write('\n')
header = "****************************************************Best Fit Results*****************************************\
***********"
results_output.write('{:>120}\n'.format(header))
results_output.write('\n')
results_output.write('{:>80}'.format('Best Spin-Orbit Angle From Minimum Chi^2 (deg): '))
results_output.write('{:<20.1f}\n'.format(best_spin_orbit_min_chi_squared))
results_output.write('{:>80}'.format('Best Spin-Orbit Angle From MCMC mean (deg): '))
results_output.write('{:<20.1f}\n'.format(best_spin_orbit_mcmc_mean))
results_output.write('{:>80}'.format('Best Spin-Orbit Angle From MCMC median (deg): '))
results_output.write('{:<20.1f}\n'.format(spin_orbit_median_all_walkers))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Spin-Orbit Angle From MCMC Variance (deg): '))
results_output.write('{:<20.1f}\n'.format(spin_orbit_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Best Vsini Value From Minimum Chi^2 (m/s): '))
results_output.write('{:<20.0f}\n'.format(best_vsini_min_chi_squared))
results_output.write('{:>80}'.format('Best Vsini From MCMC mean (m/s): '))
results_output.write('{:<20.0f}\n'.format(best_vsini_mcmc_mean))
results_output.write('{:>80}'.format('Best Vsini From MCMC median (m/s): '))
results_output.write('{:<20.0f}\n'.format(vsini_median_all_walkers))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Vsini From MCMC Variance (m/s): '))
results_output.write('{:<20.0f}\n'.format(vsini_stand_dev_mcmc_mean))
results_output.write('\n')
results_output.write('\n')
header = "****************************************************Best Fit Priors******************************************\
***********"
results_output.write('{:>120}\n'.format(header))
results_output.write('\n')
results_output.write('{:>80}'.format('Orbital Period Prior (days): '))
results_output.write('{:<20.6f}\n'.format(Orbital_period_prior))
results_output.write('{:>80}'.format('Best Orbital Period From MCMC mean (days): '))
results_output.write('{:<20.6f}\n'.format(best_orbital_period_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Orbital Period From MCMC Variance (days): '))
results_output.write('{:<20.6f}\n'.format(orbital_period_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Mid Transit Time Prior (days): '))
results_output.write('{:<20.6f}\n'.format(JD_time_mid_transit_prior))
results_output.write('{:>80}'.format('Best Mid Transit Time From MCMC mean (days): '))
results_output.write('{:<20.6f}\n'.format(best_JD_time_mid_transit_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Mid Transit Time From MCMC Variance (days): '))
results_output.write('{:<20.6f}\n'.format(JD_time_mid_transit_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('RV Zero Offset Prior (m/s): '))
results_output.write('{:<20.1f}\n'.format(RV_zero_offset_prior))
results_output.write('{:>80}'.format('Best RV Zero Offset From MCMC mean (m/s): '))
results_output.write('{:<20.1f}\n'.format(best_RV_zero_offset_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty RV Zero Offset From MCMC Variance (m/s): '))
results_output.write('{:<20.1f}\n'.format(RV_zero_offset_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('RV Offset Datasets Prior (m/s): '))
results_output.write('{:<20.1f}\n'.format(RV_offset_datasets_prior))
results_output.write('{:>80}'.format('Best RV Offset Datasets From MCMC mean (m/s): '))
results_output.write('{:<20.1f}\n'.format(best_RV_offset_datasets_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty RV Offset Datasets From MCMC Variance (m/s): '))
results_output.write('{:<20.1f}\n'.format(RV_offset_datasets_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Mp Prior (M_J or M_E): '))
results_output.write('{:<20.3f}\n'.format(Mp_prior))
results_output.write('{:>80}'.format('Best Mp From MCMC mean (M_J or M_E): '))
results_output.write('{:<20.3f}\n'.format(best_Mp_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Mp From MCMC Variance (M_J or M_E): '))
results_output.write('{:<20.3f}\n'.format(Mp_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Rp Prior (R_J or R_E): '))
results_output.write('{:<20.3f}\n'.format(Rp_prior))
results_output.write('{:>80}'.format('Best Rp From MCMC mean (R_J or R_E): '))
results_output.write('{:<20.3f}\n'.format(best_Rp_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Rp From MCMC Variance (R_J or R_E): '))
results_output.write('{:<20.3f}\n'.format(Rp_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Rp to Rs Ratio Prior: '))
results_output.write('{:<20.6f}\n'.format(Rp_Rs_ratio_prior))
results_output.write('{:>80}'.format('Best Rp to Rs Ratio From MCMC mean: '))
results_output.write('{:<20.6f}\n'.format(best_Rp_Rs_ratio_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Rp to Rs Ratio From MCMC Variance: '))
results_output.write('{:<20.6f}\n'.format(Rp_Rs_ratio_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('M* Prior (M_s): '))
results_output.write('{:<20.3f}\n'.format(Ms_solar_prior))
results_output.write('{:>80}'.format('Best M* From MCMC mean (M_s): '))
results_output.write('{:<20.3f}\n'.format(best_Ms_solar_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty M* From MCMC Variance (M_s): '))
results_output.write('{:<20.3f}\n'.format(Ms_solar_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('R* Prior (R_s): '))
results_output.write('{:<20.3f}\n'.format(Rs_solar_prior))
results_output.write('{:>80}'.format('Best R* From MCMC mean (R_s): '))
results_output.write('{:<20.3f}\n'.format(best_Rs_solar_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty R* From MCMC Variance (R_s): '))
results_output.write('{:<20.3f}\n'.format(Rs_solar_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Eccentricity Prior: '))
results_output.write('{:<20.4f}\n'.format(Ecc_prior))
results_output.write('{:>80}'.format('Best Eccentricity From MCMC mean: '))
results_output.write('{:<20.4f}\n'.format(best_Ecc_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Eccentricity From MCMC Variance: '))
results_output.write('{:<20.4f}\n'.format(Ecc_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Inclination Angle Prior (deg): '))
results_output.write('{:<20.2f}\n'.format(Inc_prior))
results_output.write('{:>80}'.format('Best Inclination Angle From MCMC mean (deg): '))
results_output.write('{:<20.2f}\n'.format(best_Inc_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Inclination Angle From MCMC Variance (deg): '))
results_output.write('{:<20.2f}\n'.format(Inc_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Argument Periastron Prior (deg): '))
results_output.write('{:<20.1f}\n'.format(omega_arg_periastron_prior + 180.0))
results_output.write('{:>80}'.format('Best Argument Periastron From MCMC mean (deg): '))
results_output.write('{:<20.1f}\n'.format(best_omega_arg_periastron_mcmc_mean + 180.0))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Argument Periastron From MCMC Variance (deg): '))
results_output.write('{:<20.1f}\n'.format(omega_arg_periastron_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Vmacro Prior (km/s): '))
results_output.write('{:<20.2f}\n'.format(vmacro_prior))
results_output.write('{:>80}'.format('Best Vmacro From MCMC mean (km/s): '))
results_output.write('{:<20.2f}\n'.format(best_vmacro_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Vmacro From MCMC Variance (km/s): '))
results_output.write('{:<20.2f}\n'.format(vmacro_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Limb Darkening q_1 Prior: '))
results_output.write('{:<20.4f}\n'.format(q_1_prior))
results_output.write('{:>80}'.format('Best q_1 From MCMC mean: '))
results_output.write('{:<20.4f}\n'.format(best_q_1_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty q_1 From MCMC Variance: '))
results_output.write('{:<20.4f}\n'.format(q_1_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Limb Darkening q_2 Prior: '))
results_output.write('{:<20.4f}\n'.format(q_2_prior))
results_output.write('{:>80}'.format('Best q_2 From MCMC mean: '))
results_output.write('{:<20.4f}\n'.format(best_q_2_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty q_2 From MCMC Variance: '))
results_output.write('{:<20.4f}\n'.format(q_2_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Impact Prior: '))
results_output.write('{:<20.4f}\n'.format(Impact_prior))
results_output.write('{:>80}'.format('Best Impact From MCMC mean: '))
results_output.write('{:<20.4f}\n'.format(best_Impact_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Impact From MCMC Variance: '))
results_output.write('{:<20.4f}\n'.format(Impact_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('Rorb over Rs ratio Prior: '))
results_output.write('{:<20.4f}\n'.format(Rorb_Rs_prior))
results_output.write('{:>80}'.format('Best Rorb over Rs ratio From MCMC mean: '))
results_output.write('{:<20.4f}\n'.format(best_Rorb_Rs_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty Rorb over Rs ratio From MCMC Variance: '))
results_output.write('{:<20.4f}\n'.format(Rorb_Rs_stand_dev_mcmc_mean))
results_output.write('\n')

results_output.write('{:>80}'.format('K amplitude Prior (km/s): '))
results_output.write('{:<20.2f}\n'.format(K_amp_prior))
results_output.write('{:>80}'.format('Best K amplitude From MCMC mean (km/s): '))
results_output.write('{:<20.4f}\n'.format(best_K_amp_mcmc_mean))
results_output.write('{:>80}'.format('1 Sigma Uncertainty K amplitude From MCMC Variance (km/s): '))
results_output.write('{:<20.4f}\n'.format(K_amp_stand_dev_mcmc_mean))
results_output.close()



results_output_theory = open(save_data_directory + 'theoretical_RV.txt', 'w')
time_string = "Time"
RV_string = "RV"
results_output_theory.write('{:50}'.format(time_string))
results_output_theory.write('     {:50}\n'.format(RV_string))
#results_output_theory.write('{:<50}     {:<50}\n'.format(time_string, RV_string)
for i in range(len(RV_array)):
    results_output_theory.write('{:<50.10e}     {:<50.5f}\n'.format(RV_array.loc[i,"Time"]/60.0, RV_array.loc[i, "RV"]))
#endfor
results_output_theory.close()

if test_null == 'Y':
    results_output_theory_null = open(save_data_directory + 'theoretical_RV_null.txt', 'w')
    time_string = "Time"
    RV_string = "RV"
    results_output_theory_null.write('{:50}'.format(time_string))
    results_output_theory_null.write('     {:50}\n'.format(RV_string))
    # results_output_theory.write('{:<50}     {:<50}\n'.format(time_string, RV_string)
    for i in range(len(RV_array_null)):
        results_output_theory_null.write(
            '{:<50.10e}     {:<50.5f}\n'.format(RV_array_null.loc[i, "Time"] / 60.0, RV_array_null.loc[i, "RV"]))
    # endfor
    results_output_theory_null.close()
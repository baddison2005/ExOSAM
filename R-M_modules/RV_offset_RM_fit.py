#RV zero offset function for RMLEAPWE.

#This python function is used to determine the best fit RV zero offset and the offset between datasets by simulating
#the out-of-transit and in-transit RVs. The fitting is usually done out of transit given the uncertianty of the RM
#effect.



#---------------------------------------Import required Python packages.-----------------------------------------------#
import scipy.special as sci
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd




#--------------------------------Grab parameters from the RM analysis program.-----------------------------------------#
def RV_offset(Index, template_RV, Index_other, template_other1, Time_transit, Number_iterations_total,
              Number_RV_offset_iterations, Number_RV_zero_offset_iterations, my_datafilelength, RV_zero_offset_end,
              RV_zero_offset_interval, RV_zero_offset_prior, RV_offset_datasets_end, RV_offset_datasets_interval,
              RV_offset_datasets_begin, datafilelength_other1, Ecc, Orbital_period, time_peri_passage,
              Bessel_function_exit, Rorb, Rorb_star, omega_arg_periastron, Inc, time_avoid, Rs, Rp, RVamplitude,
              num_my_rv, num_rv1, directory_location, Planet_name, my_data_plotting_symbols_array,
              other_data_plotting_symbols_array, RM_data_output, other_data_output):

#------------------------------------------------------constants-------------------------------------------------------#
    Rss = 6.9634E8  # Radius of our Sun (in meters).
    AU = 1.4960E11  # One astronomical unit in meters.
    pi = math.pi
    G = 6.6738E-11  # Gravitation constant.
    Mss = 1.9891E30  # Mass of the sun in kilograms.
    Time = 0  # The time for the FOR loop.
    Me = 5.9722E24  # Mass of the Earth in kg.
    Mj = 1.8986E27  # Mass of Jupiter in kg.
    Rj = 7.1492E7  # Radius of Jupiter (in meters).
    RE = 6.3781E6  # Radius of Earth (in meters if planet is given in Earth radii).
    Long_ascend_node = 90  # Set to 90. Only measuring the relative or projected inclination of orbit.

    #Setup dataframes for the analysis.
    time_data_array_columns = ["Time"]
    time_data_array = pd.DataFrame(index=Index, columns=time_data_array_columns)
    time_data_array['Time'] = template_RV['Time']
    print("time_data_array before: ", time_data_array)
    time_data_array.loc[:, 'Time'] = time_data_array.loc[:, 'Time'] * 24.0 * 3600.0 + Time_transit
    print("time_data_array after: ", time_data_array)

    time_data_array_other_columns = ["Time"]
    time_data_array_other = pd.DataFrame(index=Index_other, columns=time_data_array_other_columns)
    time_data_array_other['Time'] = template_other1['Time']
    time_data_array_other.loc[:, 'Time'] = time_data_array_other.loc[:, 'Time'] * 24.0 * 3600.0 + Time_transit

    Number_iterations_total_index = np.arange(Number_iterations_total)
    Number_RV_offset_iterations_index = np.arange(Number_RV_offset_iterations)
    Number_array_RV_zero_offset_index = np.arange(Number_RV_zero_offset_iterations)
    iteration_column = ["Int_number"]
    Number_array_total = pd.DataFrame(index=Number_iterations_total_index, columns=iteration_column)
    Number_array_RV_offset = pd.DataFrame(index=Number_RV_offset_iterations_index, columns=iteration_column)
    Number_array_RV_zero_offset = pd.DataFrame(index=Number_array_RV_zero_offset_index, columns=iteration_column)

    RV_offset_column = ["RV"]
    RV_offset_datasets_array = pd.DataFrame(index=Number_RV_offset_iterations_index, columns=RV_offset_column)
    RV_zero_offset_array = pd.DataFrame(index=Number_array_RV_zero_offset_index, columns=RV_offset_column)

    Number_RV_offset_datasets = 0
    total_interval = my_datafilelength
    Number_RV_zero_offset = 0
    Number = 0
    Iteration_number = 0
    fit_count = 0

    parameter_columns = ["RV_zero_offset", "RV_offset_datasets"]
    Parameter_array_fit = pd.DataFrame(index=Number_RV_offset_iterations_index, columns=parameter_columns)
    Parameter_array_total = pd.DataFrame(index=Number_iterations_total_index, columns=parameter_columns)

    fit_index_column = ["fit_indeces"]
    fit_array_indeces = pd.DataFrame(index=Number_RV_offset_iterations_index, columns=fit_index_column)

    chi2_column = ["chi2"]
    Rchi2_column = ["Rchi2"]
    diff_column = ["diff"]
    Chi_squared_array_total = pd.DataFrame(index=Number_iterations_total_index, columns=chi2_column)
    Reduced_Chi_Squared_array_total = pd.DataFrame(index=Number_iterations_total_index, columns=Rchi2_column)
    model_data_difference_total = pd.DataFrame(index=Number_iterations_total_index, columns=diff_column)

    Chi_squared_array_fit = pd.DataFrame(index=Number_RV_offset_iterations_index, columns=chi2_column)
    Reduced_Chi_Squared_array_fit = pd.DataFrame(index=Number_RV_offset_iterations_index, columns=Rchi2_column)
    model_data_difference_fit_total = pd.DataFrame(index=Number_RV_offset_iterations_index, columns=diff_column)

    Number_iterations_left_total = Number_iterations_total

    for RV_zero_offset_loop in range(Number_RV_zero_offset_iterations):
        k = 0
        RV_zero_offset = RV_zero_offset_end - (RV_zero_offset_interval * (RV_zero_offset_loop))
        # Put RV_zero_offset values into the RV_zero_offset array at the position Number_RV_zero_offset.
        RV_zero_offset_array.loc[Number_RV_zero_offset, "RV"] = RV_zero_offset
        # Increase Number_RV_zero_offset by 1 each time this do loop runs.
        Number_RV_zero_offset = Number_RV_zero_offset + 1

        if RV_zero_offset < RV_zero_offset_prior + (RV_zero_offset_interval/10.0) and \
           RV_zero_offset > RV_zero_offset_prior - (RV_zero_offset_interval/10.0):
            fit_flag = True
        #endif
        else:
            fit_flag = False
        #endelse

        for RV_offset_datasets_loop in range(k, Number_RV_offset_iterations):
            RV_offset_datasets = RV_offset_datasets_end - (RV_offset_datasets_interval * (RV_offset_datasets_loop))

            if Number_RV_offset_datasets < math.floor(((RV_offset_datasets_end - RV_offset_datasets_begin)
                                                      /RV_offset_datasets_interval) + 1.0):
                #Put RV_offset_datasets values into RV_offset_datasets array at the position Number_RV_offset_datasets.
                RV_offset_datasets_array.loc[Number_RV_offset_datasets, "RV"] = RV_offset_datasets
            #endif

            # Number_RV_offset_datasets by 1 each time this do loop runs.
            Number_RV_offset_datasets = Number_RV_offset_datasets + 1

            Number_array_total.loc[Number, "Int_number"] = Number

            RV_my_array_column = ["Time", "RV", "RV_error"]
            RV_my_array_index = np.arange(my_datafilelength)
            RV_my_array = pd.DataFrame(index=RV_my_array_index, columns=RV_my_array_column)

            RV_other_array_column = ["Time", "RV", "RV_error"]
            RV_other_array_index = np.arange(datafilelength_other1)
            RV_other_array = pd.DataFrame(index=RV_other_array_index, columns=RV_other_array_column)

            RV_theory_column = ["Time", "RV"]
            RV_theory = pd.DataFrame(index=RV_my_array_index, columns=RV_theory_column)
            Time_column = ['Time']
            Timestep = pd.DataFrame(index=RV_my_array_index, columns=Time_column)
            Time_of_loop = pd.DataFrame(index=RV_my_array_index, columns=Time_column)
            data_points_compare = pd.DataFrame(index=RV_my_array_index, columns=['Points'])
            Planet_star_distance_columns = ['Time', 'Distance']
            Planet_star_distance_array = pd.DataFrame(index=RV_my_array_index, columns=Planet_star_distance_columns)
            Phase_array = pd.DataFrame(index=RV_my_array_index, columns=['Phase'])
            Phase_act_array = pd.DataFrame(index=RV_my_array_index, columns=['Phase'])
            Phase_array_n = pd.DataFrame(index=RV_my_array_index, columns=['Phase'])
            True_anomaly_array = pd.DataFrame(index=RV_my_array_index, columns=['Phase'])
            ecc_anomaly_array = pd.DataFrame(index=RV_my_array_index, columns=['Phase'])
            Xpos_array = pd.DataFrame(index=RV_my_array_index, columns=['Position'])
            Ypos_array = pd.DataFrame(index=RV_my_array_index, columns=['Position'])
            Zpos_array = pd.DataFrame(index=RV_my_array_index, columns=['Position'])

            Parameter_array_total.loc[Number, "RV_zero_offset"] = RV_zero_offset
            Parameter_array_total.loc[Number, "RV_offset_datasets"] = RV_offset_datasets

            if fit_flag == True:
                Parameter_array_fit.loc[fit_count, "RV_zero_offset"] = RV_zero_offset
                Parameter_array_fit.loc[fit_count, "RV_offset_datasets"] = RV_offset_datasets
                fit_array_indeces.loc[fit_count, "fit_indeces"] = Number
                fit_count = fit_count + 1
            #endif

            #print("fit_flag: ", fit_flag)

            #Now apply the velocity offset to your data.
            RV_my_array.loc[:, 'Time'] = template_RV.loc[:, 'Time']
            RV_my_array.loc[:, 'RV'] = template_RV.loc[:, 'RV'] + (RV_zero_offset - RV_zero_offset_prior) \
                                       + RV_offset_datasets
            RV_my_array.loc[:, 'RV_error'] = template_RV.loc[:, 'RV_error']

            RV_other_array.loc[:, 'Time'] = template_other1.loc[:, 'Time']
            RV_other_array.loc[:, 'RV'] = template_other1.loc[:, 'RV'] + (RV_zero_offset - RV_zero_offset_prior)
            RV_other_array.loc[:, 'RV_error'] = template_other1.loc[:, 'RV_error']

            num_compare = 0

            data_points_compare.loc[:, 'Points'] = False

#***********************************************************************************************************************

            for Time_loop in range(my_datafilelength):

                #The time will come from the first data point which must be normalized as a fraction of a day before
                # mid transit.
                time = time_data_array.loc[Time_loop, "Time"]

                #print("time: ", time)

                #Set transit flag to False when outside of transit event.
                transit_flag = False
                #Set occultation flag to False when outside of secondary transit event.
                time_avoid_flag = False

                if Ecc == 0:
                    ecc_anomaly = ((2.0*pi)/Orbital_period)*(time - time_peri_passage)
                #endif
                else:
                    #Calculate the eccentric anomaly using the functions found near the end of the program.
                    sum_ecc_anomaly = 0
                    ecc_anomaly_before = 1000000
                    for order in range(1, 21):
                        #Calculate the value of the Bessel function which is used to find the eccentric anomaly.
                        Bessel_value = sci.jn(order, order*Ecc)
                        #Bessel_value = BESELJ(order*Ecc, order, /DOUBLE)
                        #PRINT *, "Bessel_value = ", Bessel_value
                        if order > 1:
                            ecc_anomaly_before = sum_ecc_anomaly
                        #endif
                        sum_ecc_anomaly = sum_ecc_anomaly + ((2.0/order)*Bessel_value*math.sin(order*(((2.0*pi)
                                          /Orbital_period)*(time - time_peri_passage))))
                        #PRINT *, "sum_ecc_anomaly = ", sum_ecc_anomaly
                        ecc_anomaly_after = sum_ecc_anomaly
                        if order > 1 and math.fabs(ecc_anomaly_after - ecc_anomaly_before) <= Bessel_function_exit:
                            break
                        #endif
                    #endfor
                    ecc_anomaly = (((2.0*pi)/Orbital_period)*(time - time_peri_passage)) + sum_ecc_anomaly
                #endelse

                #Now use the ecc_anomaly to determine the True_anomaly and Phase_angle for a given data point.
                True_anomaly = 2.0*(math.atan(math.tan(ecc_anomaly/2.0)*(math.sqrt((1.0 + Ecc)/(1.0 - Ecc)))))

                if Ecc == 0:
                    #The distance between the center of the planet to the center of the star in a circular orbit.
                    Planet_star_distance = Rorb + Rorb_star
                #endif
                else:
                    #The distance between the center of the planet to the center of the star in an eccentric orbit.
                    Planet_star_distance = ((Rorb + Rorb_star)*(1.0 - Ecc**2.0))/(1.0 + (Ecc*math.cos(True_anomaly)))
                #endelse
                #Time in reference to mid transit time

                #print('Planet_star_distance: ', Planet_star_distance)
                Time_ref = Time - Time_transit

                #The position of the planet on the x-axis.
                Xpos = Planet_star_distance*((-math.sin(True_anomaly)*math.sin((omega_arg_periastron)*(pi/180.0)))
                       + (math.cos(True_anomaly)*math.cos((omega_arg_periastron)*(pi/180.0))))
                #The position of the planet on the y-axis.
                Ypos = Planet_star_distance*((math.cos(True_anomaly + pi)*math.cos((Inc)*(pi/180.0))*
                       math.sin((omega_arg_periastron)*(pi/180.0))) + (math.sin(True_anomaly + pi)*
                       math.cos((Inc)*(pi/180.0))*math.cos((omega_arg_periastron)*(pi/180.0))))
                #The position of the planet on the z-axis.
                Zpos = Planet_star_distance*((-math.cos((omega_arg_periastron)*(pi/180.0))*
                       math.sin((Inc)*(pi/180.0))*math.sin(True_anomaly)) - (math.cos(True_anomaly)*
                       math.sin((Inc)*(pi/180.0))*math.sin((omega_arg_periastron)*(pi/180.0))))
                Dist2 = Xpos**2.0 + Ypos**2.0                   #Square of the planet-star apparent seperation.
                Distance_center = math.sqrt(Dist2)              #Apparent seperation between the planet and the star.

                # print('Distance_center: ', Distance_center)
                # print('Rs + Rp: ', Rs + Rp)
                # print('Distance_center - (Rs + Rp): ', Distance_center - (Rs + Rp))

                # if Xpos <= 0 and Zpos >= 0:
                #     #Planet is currently in quadrant three so add pi.
                #     phase_angle = math.atan(Xpos/Zpos) + pi
                #     #Taking into account orbital inclination.
                #     #Phase_angle_observed = math.atan(-(math.sqrt(Xpos**2.0 + Ypos**2.0))/Zpos) + pi
                # elif Xpos >= 0 and Zpos >= 0:
                #     #Planet is currently in quadrant four so add pi.
                #     phase_angle = math.atan(Xpos/Zpos) + pi
                #     #Taking into account orbital inclination.
                #     #Phase_angle_observed = math.atan((math.sqrt(Xpos**2.0 + Ypos**2.0))/Zpos) + pi
                # elif Xpos >= 0 and Zpos <= 0:
                #     #Planet is currently in quadrant one so add 2pi.
                #     phase_angle = 2.0*pi + math.atan(Xpos/Zpos)
                #     #Taking into account orbital inclination.
                #     #Phase_angle_observed = 2.0*pi + math.atan((math.sqrt(Xpos**2.0 + Ypos**2.0))/Zpos)
                # elif Xpos <= 0 and Zpos <= 0:
                #     #Planet is currently in quadrant two so add 2pi.
                #     phase_angle = math.atan(Xpos/Zpos)
                #     #Taking into account orbital inclination.
                #     #Phase_angle_observed = math.atan(-(math.sqrt(Xpos**2.0 + Ypos**2.0))/Zpos)
                # #endif

                True_phase = math.acos(math.sin(True_anomaly + (omega_arg_periastron*(pi/180.0)))*
                             math.sin(Inc*(pi/180.0)))
                # Phase_orbit_n = math.acos(math.sin(True_anomaly + (omega_arg_periastron*(pi/180.0))))

                if Time_ref >= -time_avoid and Time_ref <= time_avoid:
                    time_avoid_flag = True
                    #PRINT, "Time_ref", Time_ref
                    #PRINT, "length_time_compare", length_time_compare
                #endif

                if Distance_center <= (Rs + Rp) and Zpos > 0:
                    #If the seperation between the disk of the star and planet are less then Rs + Rp and
                    #Zpos is positive then the transit begins.
                    transit_flag = True
                    #Planet is transiting. Will be used to determine when to do the chi squared analysis.
                #endif

                #If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is negative
                #then the secondary transit (occulation) begins.
                if Distance_center <= (Rs + Rp) and Zpos < 0:
                    occultation_flag = True                  #Planet is occulting.
                #endif

                if Ecc == 0:
                    #The radial velocity of the star which includes adding the RM effect for a circular orbit.
                    RV = RVamplitude*(math.cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0))))
                else:
                    #The radial velocity of the star which includes adding the RM effect for an eccentric orbit.
                    RV = RVamplitude*(math.cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0)))
                         + (Ecc*math.cos((omega_arg_periastron)*(pi/180.0) + pi)))
                #endelse

                # print('RV: ', RV)
                # print('transit_flag: ', transit_flag)
                # print('time_avoid_flag: ', time_avoid_flag)

                #If transit flag is Y, then indicate the data point to compare.
                if transit_flag == False and time_avoid_flag == False:
                    num_compare = num_compare + 1
                    data_points_compare.loc[Time_loop, 'Points'] = True
                #endif

                Timestep.loc[Time_loop, 'Time'] = Time_ref
                Time_of_loop.loc[Time_loop, 'Time'] = Time

                RV_theory.loc[Time_loop, 'Time'] = Time_ref
                RV_theory.loc[Time_loop, 'RV'] = RV

                Planet_star_distance_array.loc[Time_loop, 'Distance'] = Planet_star_distance
                Planet_star_distance_array.loc[Time_loop, 'Time'] = Time_ref

                # Phase_array.loc[Time_loop, 'Phase'] = phase_angle
                Phase_act_array.loc[Time_loop, 'Phase'] = True_phase
                # Phase_array_n.loc[Time_loop, 'Phase'] = Phase_orbit_n
                True_anomaly_array.loc[Time_loop, 'Phase'] = True_anomaly
                ecc_anomaly_array.loc[Time_loop, 'Phase'] = ecc_anomaly
                Xpos_array.loc[Time_loop, 'Position'] = Xpos
                Ypos_array.loc[Time_loop, 'Position'] = Ypos
                Zpos_array.loc[Time_loop, 'Position'] = Zpos

            #endfor

            Number_points1 = 0
            chi_2 = 0          #Set chi squared equal to zero.
            model_data_difference = 0

            diffs_column = ['Num']
            num_compare_index = np.arange(num_compare)
            Model_theory_difference_array = pd.DataFrame(index=num_compare_index, columns=diffs_column)

            #print('num_compare: ', num_compare)

            for s in range(total_interval):
                if data_points_compare.loc[s, 'Points'] == True:

                    chi_2 = chi_2 + ((RV_my_array.loc[s, 'RV'] - RV_theory.loc[s, 'RV'])
                            /RV_my_array.loc[s, 'RV_error'])**2.0

                    Model_theory_difference_array.loc[Number_points1, 'Num'] = \
                        math.fabs(RV_my_array.loc[s, 'RV'] - RV_theory.loc[s, 'RV'])
                    model_data_difference = model_data_difference \
                                            + math.fabs(RV_my_array.loc[s, 'RV'] - RV_theory.loc[s, 'RV'])

                    #print('chi_2: ', chi_2)

                    #print('model_data_difference: ', model_data_difference)

                    # Counter to determine the number of elements used to calculate sigma squared.
                    Number_points1 = Number_points1 + 1
                #endif

            #endfor

            print("Chi squared not reduced ", chi_2)
            r_Chi_2 = chi_2/(Number_points1 - 2 - 1)
            print("Chi squared reduced ", r_Chi_2)

            Chi_squared_array_total.loc[Number, "chi2"] = chi_2              #Puts the chi squared values into an array.
            # Puts the reduced chi squared values into an array.
            Reduced_Chi_Squared_array_total.loc[Number, "Rchi2"] = r_Chi_2
            model_data_difference_total.loc[Number, "diff"] = model_data_difference

            if fit_flag == True:
                #Puts the chi squared values into an array.
                Chi_squared_array_fit.loc[fit_count - 1, "chi2"] = chi_2
                #Puts the reduced chi squared values into an array.
                Reduced_Chi_Squared_array_fit.loc[fit_count - 1, "Rchi2"] = r_Chi_2
                model_data_difference_fit_total.loc[fit_count - 1, "diff"] = model_data_difference
            #endif

            #Let's the user know the progress of the program.
            Number_iterations_left_total = Number_iterations_left_total - 1
            print("There are ",Number_iterations_left_total," interations left out of ",Number_iterations_total)
            Iteration_number = Iteration_number + 1                    #Let's the user know the progress of the program.
            print("Iteration number is ", Iteration_number)
            print("Iteration fit number is ", fit_count)

            percent_processed = (100.0 * Iteration_number)/Number_iterations_total
            string_bar = ''

            if percent_processed >= 4.0:
                for percent_processed_loop in range(math.floor(percent_processed/4.0)):
                    string_bar = '|' + string_bar                           #At every 4% add a bar to the progress bar.
                #endfor
            #endif

            print("Processing data: [", '{:25}'.format(string_bar), "]"," ", '{:d}'.format(round(percent_processed)),
                  "% complete")

            Number = Number + 1                                         #Counter variable.

            print("RV_zero_offset: ", RV_zero_offset)
            print("RV_offset_datasets: ", RV_offset_datasets)
        #endfor
    #endfor

    print("number of total iterations: ", Number)
    print("number of iterations there are suppose to be in total: ", Number_iterations_total)
    print("number of fit iterations: ", fit_count)
    print("number of iterations there are suppose to be in fit count: ", Number_RV_offset_iterations)

    #Determine the lowest value for Chi square, reduced Chi square, and where these occur
    #for the full parameter space and for just the parameters that will be used for the fit.
    #Also find the parameters which match these lowest values.
    min_chi_squared_total = Chi_squared_array_total["chi2"].min()
    loc_min_chi_squared_total = Chi_squared_array_total["chi2"].values.argmin()
    print("mininum chi squared value full parameter space: ", min_chi_squared_total)
    print("mininum chi squared location full parameter space: ", loc_min_chi_squared_total)
    min_reduced_chi_squared_total = Reduced_Chi_Squared_array_total["Rchi2"].min()
    loc_min_reduced_chi_squared_total = Reduced_Chi_Squared_array_total["Rchi2"].values.argmin()
    print("mininum reduced chi squared value of full parameter space: ", min_reduced_chi_squared_total)
    print("mininum reduced chi squared location full parameter space: ", loc_min_reduced_chi_squared_total)

    #Determine the lowest value for Chi square, reduced Chi square, and where these occur
    #only in the parameter space used for the fit. Also find the parameters which match these lowest values.
    min_chi_squared_fit = Chi_squared_array_fit["chi2"].min()
    loc_min_chi_squared_fit = Chi_squared_array_fit["chi2"].values.argmin()
    print("mininum chi squared value restrictive fitting parameter space: ", min_chi_squared_fit)
    print("mininum chi squared location restrictive parameter space: ", loc_min_chi_squared_fit)
    min_reduced_chi_squared_fit = Reduced_Chi_Squared_array_fit["Rchi2"].min()
    loc_min_reduced_chi_squared_fit = Reduced_Chi_Squared_array_fit["Rchi2"].values.argmin()
    print("mininum reduced chi squared value restrictive fitting parameter space: ", min_reduced_chi_squared_fit)
    print("mininum reduced chi squared location restrictive parameter space: ", loc_min_reduced_chi_squared_fit)

    best_RV_zero_offset_total = Parameter_array_total.loc[loc_min_chi_squared_total, "RV_zero_offset"]
    print("The best RV zero offset corresponding to the lowest value of chi squared in full parameter space: ",
          best_RV_zero_offset_total, "m/s")
    best_RV_offset_datasets_total = Parameter_array_total.loc[loc_min_chi_squared_total, "RV_offset_datasets"]
    print("The best RV offset between the two data sets corresponding to the lowest value of chi squared in full "
          "parameter space: ", best_RV_offset_datasets_total, "m/s")

    best_RV_zero_offset_fit = Parameter_array_fit.loc[loc_min_chi_squared_fit, "RV_zero_offset"]
    print("The best RV zero offset corresponding to the lowest value of chi squared in restricted parameter space: ",
          best_RV_zero_offset_fit, "m/s")
    best_RV_offset_datasets_fit = Parameter_array_fit.loc[loc_min_chi_squared_fit, "RV_offset_datasets"]
    print("The best RV offset between the two data sets corresponding to the lowest value of chi squared in restricted"
          " parameter space: ", best_RV_offset_datasets_fit, "m/s")

    count_chi_total = 0

    #print("Chi_squared_array_total: ", Chi_squared_array_total)

    #print("Number of Nans: ", Chi_squared_array_total["chi2"].isnull().sum().sum())


    Chi_1sigma_array_indices_total = Chi_squared_array_total[Chi_squared_array_total.loc[:,"chi2"]
                                                             <= min_chi_squared_total + 2.296].index.tolist()

    print("Chi_1sigma_array_indices_total: ", Chi_1sigma_array_indices_total)
    count_chi_total = len(Chi_1sigma_array_indices_total)
    print("count_chi_total: ", count_chi_total)


    Parameter_chi_total_array_index = np.arange(count_chi_total)
    Parameter_chi_total_columns = ["RV_zero_offset", "RV_offset_datasets", "chi2"]
    Parameter_chi_total_array = pd.DataFrame(index=Parameter_chi_total_array_index, columns=Parameter_chi_total_columns)

    for k in range(count_chi_total):
        Parameter_chi_total_array.loc[k, "RV_zero_offset"] = \
            Parameter_array_total.loc[Chi_1sigma_array_indices_total[k], "RV_zero_offset"]

        Parameter_chi_total_array.loc[k, "RV_offset_datasets"] = \
            Parameter_array_total.loc[Chi_1sigma_array_indices_total[k], "RV_offset_datasets"]

        Parameter_chi_total_array.loc[k, "chi2"] = \
            Chi_squared_array_total.loc[Chi_1sigma_array_indices_total[k], "chi2"]

    #endfor

    # Parameter_chi_total_array.loc[0:count_chi_total-1, "RV_zero_offset"] = \
    #     Parameter_array_total.loc[Chi_1sigma_array_indices_total[0:count_chi_total-1], "RV_zero_offset"]
    #
    # Parameter_chi_total_array.loc[0:count_chi_total-1, "RV_offset_datasets"] = \
    #     Parameter_array_total.loc[Chi_1sigma_array_indices_total[0:count_chi_total-1], "RV_offset_datasets"]
    #
    # Parameter_chi_total_array.loc[0:count_chi_total-1, "chi2"] = \
    #     Chi_squared_array_total.loc[Chi_1sigma_array_indices_total[0:count_chi_total-1], "chi2"]


    print("Parameter_chi_total_array: ", Parameter_chi_total_array)


    RV_offset_datasets_plus_error_total = Parameter_chi_total_array.loc[0:count_chi_total-1, "RV_offset_datasets"].max()
    RV_offset_datasets_subscript_max = Parameter_chi_total_array.loc[0:count_chi_total-1, "RV_offset_datasets"]\
        .values.argmax()
    print("RV offset between data sets plus error full parameter space: ",
          RV_offset_datasets_plus_error_total - best_RV_offset_datasets_total)
    chi_square_RV_offset_datasets_plus_total = Parameter_chi_total_array.loc[RV_offset_datasets_subscript_max, "chi2"]
    print("Chi squared value at RV offset between data sets plus error full parameter space: ",
          chi_square_RV_offset_datasets_plus_total)

    RV_offset_datasets_minus_error_total = Parameter_chi_total_array.loc[0:count_chi_total-1,
                                           "RV_offset_datasets"].min()
    RV_offset_datasets_subscript_min = Parameter_chi_total_array.loc[0:count_chi_total-1, "RV_offset_datasets"]\
        .values.argmin()
    print("RV offset between data sets minus error full parameter space: ",
          RV_offset_datasets_minus_error_total - best_RV_offset_datasets_total)
    chi_square_RV_offset_datasets_minus_total = Parameter_chi_total_array.loc[RV_offset_datasets_subscript_min, "chi2"]
    print("Chi squared value at RV offset between data sets minus error full parameter space: ",
          chi_square_RV_offset_datasets_minus_total)

    counter_error_loop = 0
    for i in range(loc_min_chi_squared_fit - 1, -1, -1):
        flag_exit = False
        if Parameter_array_fit.loc[i, 'RV_zero_offset'] == best_RV_zero_offset_fit and \
            Parameter_array_fit.loc[i, 'RV_offset_datasets'] > best_RV_offset_datasets_fit and \
            Chi_squared_array_fit.loc[i, 'chi2'] >= min_chi_squared_fit + 1.0:

            RV_offset_datasets_plus_error_fit = Parameter_array_fit.loc[i, 'RV_offset_datasets']
            print("RV offset between data sets plus error restricted parameter space: ",
                  RV_offset_datasets_plus_error_fit - best_RV_offset_datasets_fit)
            chi_square_RV_offset_datasets_plus_fit = Chi_squared_array_fit.loc[i, 'chi2']
            print("Chi squared value at RV offset between data sets plus error restricted parameter space: ",
                  chi_square_RV_offset_datasets_plus_fit)
            flag_exit = True
        #endif
        if flag_exit == True:
            break       #Exit do loop. Found plus error.
        #endif
        counter_error_loop = counter_error_loop + 1
    #endfor


    counter_error_loop = 0
    for i in range(loc_min_chi_squared_fit - 1, Number_RV_offset_iterations):
        flag_exit = False
        if Parameter_array_fit.loc[i, 'RV_zero_offset'] == best_RV_zero_offset_fit and \
                    Parameter_array_fit.loc[i, 'RV_offset_datasets'] < best_RV_offset_datasets_fit and \
                    Chi_squared_array_fit.loc[i, 'chi2'] >= min_chi_squared_fit + 1.0:
            RV_offset_datasets_minus_error_fit = Parameter_array_fit.loc[i, 'RV_offset_datasets']
            print("RV offset between data sets minus error restricted parameter space: ",
                  RV_offset_datasets_minus_error_fit - best_RV_offset_datasets_fit)
            chi_square_RV_offset_datasets_minus_fit = Chi_squared_array_fit.loc[i, 'chi2']
            print("Chi squared value at RV offset between data sets minus error restricted parameter space: ",
                chi_square_RV_offset_datasets_minus_fit)
            flag_exit = True
        #endif
        if flag_exit == True:
            break   #Exit do loop. Found minus error.
        #endif
        counter_error_loop = counter_error_loop + 1
    #endfor

    #Apply the best RV offset between datasets.
    z = num_my_rv

    RV_my_array1_column = ["Time", "RV", "RV_error"]
    RV_my_array1_index = np.arange(z)
    RV_my_array1 = pd.DataFrame(index=RV_my_array1_index, columns=RV_my_array1_column)
    RV_my_array1.loc[:, "Time"] = template_RV.loc[:, 'Time']
    RV_my_array1.loc[:, 'RV'] = template_RV.loc[:, 'RV'] + best_RV_offset_datasets_fit
    RV_my_array1.loc[:, 'RV_error'] = template_RV.loc[:, 'RV_error']

    total_datalength = num_my_rv + num_rv1
    z = total_datalength

    RV_my_array_time = pd.DataFrame(index=RV_my_array1_index, columns=['Time'])
    RV_my_array_RV = pd.DataFrame(index=RV_my_array1_index, columns=['RV'])
    RV_my_array_RV_error = pd.DataFrame(index=RV_my_array1_index, columns=['RV_error'])
    RV_my_array_time.loc[:, "Time"] = template_RV.loc[:, 'Time']
    RV_my_array_RV.loc[:, "RV"] = template_RV.loc[:, 'RV'] + best_RV_offset_datasets_fit
    RV_my_array_RV_error.loc[:, "RV_error"] = template_RV.loc[:, 'RV_error']

    #Now combine all the RV from different arrays into one array.
    RV_column = ["Time", "RV", "RV_error"]
    RV_index = np.arange(z)
    RV_array = pd.DataFrame(index=RV_index, columns=RV_column)
    RV_array.loc[0:z-1, "Time"] = pd.concat([template_other1.loc[:, "Time"], RV_my_array_time.loc[:, "Time"]],
                                            ignore_index=True)
    RV_array.loc[0:z-1, "RV"] = pd.concat([template_other1.loc[:, "RV"], RV_my_array_RV.loc[:, "RV"]],
                                            ignore_index=True)
    RV_array.loc[0:z-1, "RV_error"] = pd.concat([template_other1.loc[:, "RV_error"],
                                                 RV_my_array_RV_error.loc[:, "RV_error"]], ignore_index=True)

    #Sort the RV array in time to feed into Fortran routine.
    RV_array = RV_array.sort_values('Time', ascending=True)
    RV_array = RV_array.reset_index(drop=True)

    #This part outputs all the values found in the model in a text file.
    offset_between_datasets = open(directory_location + 'offset_between_datasets.txt', 'w')
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format('Best RV offset datasets whole parameter space:',
                                                              best_RV_offset_datasets_total))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format(
        'Positive uncertainty on RV offset datasets whole parameter space:',
        RV_offset_datasets_plus_error_total - best_RV_offset_datasets_total))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format(
        'Negative uncertainty on RV offset datasets whole parameter space:',
        RV_offset_datasets_minus_error_total - best_RV_offset_datasets_total))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format('Best RV offset datasets fitted parameter space:',
                                                            best_RV_offset_datasets_fit))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format(
        'Positive uncertainty on RV offset datasets fitted parameter space:',
        RV_offset_datasets_plus_error_fit - best_RV_offset_datasets_fit))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format(
        'Negative uncertainty on RV offset datasets fitted parameter space:',
        RV_offset_datasets_minus_error_fit - best_RV_offset_datasets_fit))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format('Best RV zero offset whole parameter space:',
                                                              best_RV_zero_offset_total))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format('Min. value of Chi^2 whole parameter space:',
                                                              min_chi_squared_total))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format('Min. value of Reduced Chi^2 whole parameter space:',
                                                              min_reduced_chi_squared_total))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format('Min. value of Chi^2 restricted parameter space:',
                                                              min_chi_squared_fit))
    offset_between_datasets.write('{:>75} {:<20.8f}\n'.format('Min. value of Reduced Chi^2 restricted parameter space:',
                                                              min_reduced_chi_squared_fit))
    offset_between_datasets.write('{:>75} {:<20d}\n'.format('location of chi^2 whole parameter space:',
                                                              loc_min_chi_squared_total))
    offset_between_datasets.write('{:>75} {:<20d}\n'.format('location of Reduced chi^2 whole parameter space:',
                                                            loc_min_reduced_chi_squared_total))
    offset_between_datasets.write('{:>75} {:<20d}\n'.format('location of chi^2 restricted parameter space:',
                                                            loc_min_chi_squared_fit))
    offset_between_datasets.write('{:>75} {:<20d}\n'.format('location of Reduced chi^2 restricted parameter space:',
                                                            loc_min_reduced_chi_squared_fit))
    offset_between_datasets.close()

    #Plot the RV data.
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
    #print(my_data_plotting_symbols_array)
    for i in range(num_rv1):
        if other_data_plotting_symbols_array.loc[i, 'Symbols'] == 36:
            #The symbol is a cirlce.
            other_data_plotting_symbols_array.loc[i, 'Symbols'] = 'x'
        #endif
    #endfor

    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    for i in range(num_my_rv):
        plt.errorbar(RV_my_array1.loc[i, "Time"], RV_my_array1.loc[i, "RV"], yerr=RV_my_array1.loc[i, "RV_error"],
                     fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='r', mfc='k', elinewidth=2,
                     markersize=6, capsize=7)
    #endfor
    for i in range(num_rv1):
        plt.errorbar(template_other1.loc[i, "Time"], template_other1.loc[i, "RV"],
                     yerr=template_other1.loc[i, "RV_error"], fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='b', color='k', elinewidth=2, markersize=10, capsize=7)
    #endfor
    plt.xlabel('Time (HJD)', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
    plt.savefig(directory_location + 'RV offset all RV data' + Planet_name + '.pdf',dpi=150)
    plt.close(fig)

    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    for i in range(num_rv1):
        plt.errorbar(template_other1.loc[i, "Time"], template_other1.loc[i, "RV"],
                     yerr=template_other1.loc[i, "RV_error"], fmt=other_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='b', color='k', elinewidth=2, markersize=10, capsize=7)
    #endfor
    plt.xlabel('Time (HJD)', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
    plt.savefig(directory_location + 'RV offset other RV data' + Planet_name + '.pdf',dpi=150)
    plt.close(fig)

    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    for i in range(num_my_rv):
        plt.errorbar(RV_my_array1.loc[i, "Time"], RV_my_array1.loc[i, "RV"],
                     yerr=RV_my_array1.loc[i, "RV_error"], fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'],
                     ecolor='r', color='k', elinewidth=2, markersize=10, capsize=7)
    #endfor
    plt.xlabel('Time (HJD)', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=16)
    plt.savefig(directory_location + 'RV offset my RV data' + Planet_name + '.pdf',dpi=150)
    plt.close(fig)

    #Output array data for Fortran to read in.
    RM_data_output_file = open(RM_data_output, 'w')
    for i in range(total_datalength):
        RM_data_output_file.write('{:<50.10e}     {:<50.5f}     {:<50.5f}\n'.format(RV_array.loc[i, 'Time'],
                                                                                    RV_array.loc[i, 'RV'],
                                                                                    RV_array.loc[i, 'RV_error']))
    #endfor
    RM_data_output_file.close()

    #Output array data for Fortran to read in.
    RM_data_output_file2 = open(other_data_output, 'w')
    for i in range(num_rv1):
        RM_data_output_file2.write('{:<50.10e}     {:<50.5f}     {:<50.5f}\n'.format(template_other1.loc[i,'Time'],
                                                                                     template_other1.loc[i,'RV'],
                                                                                     template_other1.loc[i,'RV_error']))
    #endfor
    RM_data_output_file2.close()

    RV_offset_datasets_plus_error = RV_offset_datasets_plus_error_total - best_RV_offset_datasets_total
    RV_offset_datasets_minus_error = RV_offset_datasets_minus_error_total - best_RV_offset_datasets_total

    if abs(RV_offset_datasets_plus_error) >= abs(RV_offset_datasets_minus_error):
        max_RV_offset_datasets_error = RV_offset_datasets_plus_error*2.0
    #endif
    else:
        max_RV_offset_datasets_error = abs(RV_offset_datasets_minus_error)*2.0
    #endelse
    RV_offset_datasets_new_prior = best_RV_offset_datasets_fit
    RV_offset_datasets_new_end = best_RV_offset_datasets_fit + (max_RV_offset_datasets_error)
    RV_offset_datasets_new_begin = best_RV_offset_datasets_fit - (max_RV_offset_datasets_error)
    RV_offset_datasets_new_interval = max_RV_offset_datasets_error/2.0
    RV_offset_datasets_new_1sigerr = (abs(RV_offset_datasets_plus_error)+abs(RV_offset_datasets_minus_error))/2.0

    return RV_offset_datasets_new_prior, RV_offset_datasets_new_end, RV_offset_datasets_new_begin, \
           RV_offset_datasets_new_interval, RV_offset_datasets_new_1sigerr


#enddef
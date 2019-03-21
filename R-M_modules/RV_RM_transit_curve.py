#RV, RM, Transit lightcurve model array function for RMLEAPWE.

#This python function is used to generate the best fit Rossiter-McLaughlin effect anomaly curve, RV Keplerian orbit,
#transit light curve, and Occultation light curve. It uses the best fit spin-orbit angle and vsini values determined in
#the MCMC analysis and the priors given for the other parameters to generate the model.



#---------------------------------------Import required Python packages.-----------------------------------------------#
import scipy.special as sci
import math
import numpy as np
import pandas as pd




#--------------------------------Grab parameters from the RM analysis program.-----------------------------------------#
def RV_model_array(data_plot_model_interval, Number_orbits, Bessel_function_exit, Jupiter_units, use_Rp_Rs_ratio, vsini,
                   stellar_rotation_angle, orbital_period1, JD_time_mid_transit, Mp_best, Mp_fixed, Rp_best, Ms_solar,
                   Rs_solar, Rp_Rs_ratio, Ecc, Inc, omega_arg_periastron, beta_rot, vmacro, linear_quadratic, q_1, q_2,
                   u, K_amp, use_K_amp, Impact, use_Impact, Rorb_Rs, use_Rorb_Rs_ratio, test_null_hypothesis, Albedo,
                   M_pixels):




# ------------------------------------------------------constants-------------------------------------------------------#
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




#---------------------------------------Planetary orbital calculations-------------------------------------------------#
    Ms = Mss * Ms_solar  # Mass of the star (in kg).
    Orbital_period_day = orbital_period1
    Orbital_period = orbital_period1 * 3600.0 * 24.0
    if use_Impact == 'Y' and use_Rorb_Rs_ratio == 'Y':
        Inc = math.acos((Impact/(Rorb_Rs))*((1.0 +
                                            (Ecc*math.sin(omega_arg_periastron*(pi/180.0))))/(1.0 - Ecc**2.0)))*\
                                            (180.0/pi)
    #END IF

    if use_K_amp == 'Y' and Mp_fixed == 'Y' and Ecc == 0:
        Mp = ((Orbital_period)/(2.0*pi*G))**(1.0/3.0)*((Ms**(2.0/3.0)*K_amp)/math.sin(Inc*(pi/180.0)))
    else:
        if Jupiter_units == 'Y':
            Mp = Mp_best*Mj
        #END IF

        if Jupiter_units == 'N':
            Mp = Mp_best*Me
        #END IF
    #END IF

    Rorb = ((Orbital_period**2.0*G*(Ms + Mp))/(4.0*pi**2.0))**(1.0/3.0)     #Semi-major axis.
    #Semi-major of stellar axis.
    Rorb_star = ((Orbital_period**2.0*G*((Mp**(3.0))/(Ms + Mp)**(2.0)))/(4.0*pi**2.0))**(1.0/3.0)

    if use_Rorb_Rs_ratio == 'Y':
        Rs = Rorb/(Rorb_Rs)
    else:
        Rs = Rss * Rs_solar
    #END IF

    #Based on given priors, determine the approximate time in the simulation to start and finish comparing
    #RV's with model.
    if use_Rp_Rs_ratio == 'Y':
        Rp = Rp_Rs_ratio * Rs
    else:
        if Jupiter_units == 'Y':
            Rp = Rp_best*Rj
        #END IF

        if Jupiter_units == 'N':
            Rp = Rp_best*RE
        #END IF
    #END IF

    RV = 0.0                                        #Set the RV to zero.
    Rs2 = Rs**2.0                                   #Variable to speed up calculations.
    Rp2 = Rp**2.0                                   #Variable to speed up calculations.
    L_total_star = 1.0                              #Normalize total luminosity of star to 1.
    Total_L = 1.0                                   #Set the total amount of flux from the star to 1.
    Total_RM = 0.0                                  #Set the Rossiter-McLaughlin anomaly to zero.

    if linear_quadratic == 'q':
        Io = 6.0 / (pi*Rs2*(6.0 - (2.0*q_1) - q_2)) #Initial light intensity equation with limb darkening (quadratic)
    else:
        Io = 1.0 / (pi*Rs2*(1.0-(u/3.0)))           #The initial light intensity equation with limb darkening (linear)
                                                    #(normalize Io such that total star luminosity is 1).
    #endelse

    Aplan = pi*Rp2                                  #Surface area of the planet.
    vturb = math.sqrt(beta_rot**2.0 + vmacro**2.0)    #Velocity width of spectral line due to mechanisms other than
                                                    #rotation (i.e. micro and macro turbulence).


    if test_null_hypothesis:
        vsini = 0.0
        vturb = 0.0

    True_anomaly_start = pi - omega_arg_periastron*(pi/180.0)
    print("True_anomaly_start: ", True_anomaly_start*(180.0/pi))

    if True_anomaly_start >= pi:
        ecc_anomaly_start = 2.0*math.atan(math.tan(True_anomaly_start/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) + 2.0*pi
    elif True_anomaly_start <= -pi:
        ecc_anomaly_start = 2.0*math.atan(math.tan(True_anomaly_start/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) - 2.0*pi
    else:
        ecc_anomaly_start = 2.0*math.atan(math.tan(True_anomaly_start/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc)))
    #endelse

    #Calculate the amount of time (in seconds) after the passage of periastron has occurred when the planet is at mid
    #transit. This is determined from the transit mid time and the argument of periastron.
    if omega_arg_periastron > 180.0:
        Mean_anomaly_transit = (2.0*pi) - (omega_arg_periastron*(pi/180.0)) + pi
    else:
        Mean_anomaly_transit = pi - (omega_arg_periastron*(pi/180.0))
    #endelse

    JD_time_peri = (JD_time_mid_transit*24.0*3600.0) - ((Mean_anomaly_transit*Orbital_period)/(2.0*pi))
    print("JD_time_peri: ", JD_time_peri)

    time_peri_passage = (JD_time_mid_transit*24.0*3600.0) - JD_time_peri
    print("time_peri_passage: ", time_peri_passage)

    if use_K_amp == 'Y':
        RVamplitude = K_amp
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

    if Ecc == 0:
        Time_start = ((ecc_anomaly_start * Orbital_period) / (2.0 * pi)) + time_peri_passage
    else:
        Time_start = (((ecc_anomaly_start-(Ecc*math.sin(ecc_anomaly_start)))*Orbital_period)/(2.0*pi))+time_peri_passage
    #endelse

    print("ecc_anomaly_start: ", ecc_anomaly_start*(180.0/pi))
    print("Time_start: ", Time_start)

    True_anomaly_transit = ((3.0*pi)/2.0) - (omega_arg_periastron*(pi/180.0))
    print("True_anomaly_transit: ", True_anomaly_transit*(180.0/pi))

    if True_anomaly_transit >= pi:
        ecc_anomaly_transit = 2.0*math.atan(math.tan(True_anomaly_transit/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) + 2.0*pi
    elif True_anomaly_transit <= -pi:
        ecc_anomaly_transit = 2.0*math.atan(math.tan(True_anomaly_transit/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc))) - 2.0*pi
    else:
        ecc_anomaly_transit = 2.0*math.atan(math.tan(True_anomaly_transit/2.0)*math.sqrt((1.0-Ecc)/(1.0+Ecc)))
    #endelse

    if Ecc == 0:
        Time_transit = ((ecc_anomaly_transit*Orbital_period)/(2.0*pi)) + time_peri_passage
    else:
        Time_transit = (((ecc_anomaly_transit - (Ecc*math.sin(ecc_anomaly_transit)))*Orbital_period)/(2.0*pi)) \
                       + time_peri_passage
    #endelse







#***********************************************************************************************************************
    total_interval = round((Orbital_period*Number_orbits)/data_plot_model_interval)

    Time_transit_start = 0                #Set the transit start time to zero.
    Transit_start_no_inc = 0              #Set the transit start time for 90 inc to zero.
    Time_occultation_start = 0            #Set the occultation start time to zero.
    Time_transit_end = 0                  #Set the transit end time to zero.
    Transit_end_no_inc = 0                #Set the transit end time for 90 inc to zero.
    Time_occultation_end = 0              #Set the occultation end time to zero.
    Time_mid_transit = 0                  #Set the mid transit time to zero.
    Time_mid_occultation = 0              #Set the mid occultation time to zero.
    Planet_star_distance = 0              #Set the distance between center of the star to the center of the planet
                                          #(orbital radius) to zero.
    Transit_start_no_inc_position = 0
    transit_start_position = 0
    Occultation_start_position = 0
    transit_End_position = 0
    Occultation_End_position = 0
    transit_mid_position = 0
    occultation_mid_position = 0
    transit_End_position_no_inc = 0
    occultation_End_position = 0

    LC_model_columns = ["Time", "Flux"]
    RV_columns = ["Time", "RV"]
    Index = np.arange(total_interval)
    Transit_LC_array = pd.DataFrame(index=Index, columns=LC_model_columns)
    RM_effect_array = pd.DataFrame(index=Index, columns=RV_columns)
    RV_array = pd.DataFrame(index=Index, columns=RV_columns)
    Timestep = pd.DataFrame(index=Index, columns=["Time"])




#***********************************************************************************************************************

    for Time_loop in range(total_interval):

        Time = Time_start + (Time_loop*data_plot_model_interval)

        #Set transit flag to False when outside of transit event.
        transit_flag = False
        #Set occultation flag to False when outside of secondary transit event.
        occultation_flag = False
        time_avoid_flag = False

        if Ecc == 0:
            ecc_anomaly = ((2.0*pi)/Orbital_period)*(Time - time_peri_passage)
        else:
            #Calculate the eccentric anomaly using the functions found near the end of the program.
            sum_ecc_anomaly = 0
            ecc_anomaly_before = 1000000
            for order in range(1, 21):
                #Calculate the value of the Bessel function which is used to find the eccentric anomaly.
                Bessel_value = sci.jn(order, order*Ecc)
                if order > 1:
                    ecc_anomaly_before = sum_ecc_anomaly
                #endif
                sum_ecc_anomaly = sum_ecc_anomaly + ((2.0/order)*Bessel_value*math.sin(order*(((2.0*pi)
                                  /Orbital_period)*(Time - time_peri_passage))))
                ecc_anomaly_after = sum_ecc_anomaly
                if order > 1 and math.fabs(ecc_anomaly_after - ecc_anomaly_before) <= Bessel_function_exit:
                    break
                #endif
            #endfor
            ecc_anomaly = (((2.0*pi)/Orbital_period)*(Time - time_peri_passage)) + sum_ecc_anomaly
        #endelse

        #Now use the ecc_anomaly to determine the True_anomaly and Phase_angle for a given data point.
        True_anomaly = 2.0*(math.atan(math.tan(ecc_anomaly/2.0)*(math.sqrt((1.0 + Ecc)/(1.0 - Ecc)))))

        if Ecc == 0:
            #The time of the similation for a specific true_anomaly with the mid transit time equal to 0.
            #Time_check = ((ecc_anomaly*Orbital_period)/(2.0*pi)) + time_peri_passage
            #The distance between the center of the planet to the center of the star in a circular orbit.
            Planet_star_distance = Rorb + Rorb_star
        else:
            #The time of the similation for a specific true_anomaly.
            #Time_check = (((ecc_anomaly - (Ecc*math.sin(ecc_anomaly)))*Orbital_period)/(2.0*pi)) + time_peri_passage
            #The distance between the center of the planet to the center of the star in an eccentric orbit.
            Planet_star_distance = ((Rorb + Rorb_star)*(1.0 - Ecc**2.0))/(1.0 + (Ecc*math.cos(True_anomaly)))
        #endelse

        #Time in reference to mid transit time
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
        Dist2 = Xpos**2.0 + Ypos**2.0                   #Square of the planet-star apparent separation.
        Distance_center = math.sqrt(Dist2)              #Apparent separation between the planet and the star.
        Lblocked = 0.0
        Lblocked2 = 0.0                                 #A variable for the Anomalous velocity equation.
        v_rm = 0.0                                      #Anomalous velocity of each pixel set to zero.
        Total_RM = 0.0

        True_phase = math.acos(math.sin(True_anomaly + (omega_arg_periastron*(pi/180.0)))*math.sin(Inc*(pi/180.0)))

        #If the planet is neither in front of the star (transit) or behind the star (occultation), then calculate the
        #flux being reflected off the surface of the exoplanet based on its bond albedo, radius, phase angle, etc.
        if Distance_center > (Rs + Rp):
            Aplan = pi*Rp2                                               #Surface area of the planet.
            Total_L = 1.0 + ((Aplan/(4.0*pi*Planet_star_distance**2.0))*Albedo*(0.5*(1.0+math.cos(True_phase))))
        #endif

        if math.fabs(Xpos) <= (Rs + Rp) and Zpos > 0 and Transit_start_no_inc == 0:
            #If the separation between the disk of the star and planet @ an inclination of 90 are less then
            #Rs + Rp and Zpos is positive then the transit begins.
            Transit_start_no_inc = Time_ref              #Indicates when the transit starts at 90 inc. If statement
                                                         #prevents this variable from being overwritten every time the
                                                         #loop runs.
        #endif

        if Distance_center <= (Rs + Rp) and Zpos > 0:
            #If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is positive
            #then the transit begins.

            if Time_transit_start == 0:
                Time_transit_start = Time_ref          #Indicates when the transit starts. If statement prevents this
                                                       #variable from being overwritten every time the loop runs.
            #endif

            if (Rp2/Rs2) >= 0.030:
                Radius_planet_array = M_pixels / 2.0      #Radius of planet in the pixel array.
                Pixel = Rp/Radius_planet_array            #The number of meters per pixel.
                Area_pixel = Pixel**2.0                   #The area of each pixel.
                Io_Pixel = Io * Area_pixel                #Variable to speed up calculations (also physically
                                                          #represents the luminosity of the brightest pixel).

                for i in range(M_pixels):
                    X_pixel = (Pixel * i) + (Xpos - Rp)   #Calculates the location of the pixel on the x-axis.
                    X_pixel2 = X_pixel**2.0
                    XpXp = X_pixel * Xpos                 #temporary var for speed calculation

                    for j in range(M_pixels):
                        Y_pixel = (Pixel * j) + math.fabs(Ypos) - Rp     #Calculates the location of the pixel on the
                                                                         # y-axis.
                        # Calculates the location of the pixel along the x-axis of the rotation axis of the star.
                        X_pixel_prime = (X_pixel*math.cos((stellar_rotation_angle*pi)/180.0)) + \
                                        (Y_pixel*math.sin((stellar_rotation_angle*pi)/180.0))

                        Dist_center_pixel = X_pixel2 + Y_pixel**2.0    #squared distance of pixel from star.
                        #Calculates the values of the pixels according to how far away they are from the center of the
                        #star and the plane.
                        #squared distance of pixel from planet using a limb darkening equation.
                        Dist_planet_pixel = Dist_center_pixel  - 2.0*(XpXp + Y_pixel*Ypos) + Dist2
                        Sub_planet_velocity = vsini*(X_pixel_prime/Rs)
                        Sub_vturb = vturb                                    #The turbulence velocity blocked by planet.

                        if Dist_center_pixel <= Rs2 and Dist_planet_pixel <= Rp2:
                            if linear_quadratic == 'q':
                                #Quadratic limb darkening equation.
                                Lblocked2 = Io_Pixel*(1.0-q_1*(1.0-math.sqrt(math.fabs(1.0-(Dist_center_pixel/Rs2)))) -
                                            q_2*(1.0 - math.sqrt(math.fabs(1.0-(Dist_center_pixel/Rs2))))**2.0)
                            else:
                                #First order linear limb darkening equation.
                                Lblocked2 = Io_Pixel*(1.0-u*(1.0-math.sqrt(math.fabs(1.0-(Dist_center_pixel/Rs2)))))
                            #endelse

                            Lblocked = Lblocked + Lblocked2                     #First order limb darkening equation.
                            # Anomalous velocity of each pixel.
                            if test_null_hypothesis:
                                v_rm = 0.0
                            else:
                                v_rm = v_rm - ((Lblocked2*Sub_planet_velocity)*((((2.0*vturb**2.0)+(2.0*vsini**2.0))/
                                       ((2.0*vturb**2.0) + vsini**2.0))**(3.0/2.0))*(1.0-((Sub_planet_velocity**2.0)/
                                       ((2.0*vturb**2.0)+(vsini**2.0)))))

                        #endif
                    #endfor
                #endfor

                Total_L = 1.0 - Lblocked
                Total_RM = 0.0 + v_rm                       #Total anomalous velocity for all the pixels.

            elif Rp2/Rs2 <= 0.030:

                #Calculates the location of the center of the planet along the x-axis of the rotation axis of the star.
                X_prime = (Xpos*math.cos((stellar_rotation_angle*pi)/180.0)) + \
                          (Ypos*math.sin((stellar_rotation_angle*pi)/180.0))
                set_distance_center = Distance_center     #The limb darkening equation will use this distance as long as
                                                          #the center of the planet is inside the radius of the star.
                Sub_planet_velocity = vsini*(X_prime/Rs)      #Calculate the subplanetary velocity (the stellar velocity
                                                              #blocked by the planetary disc) from vsini times the
                                                              #from the center of the planet to the stellar rotation
                                                              #axis divided by the radius of the star.
                Sub_vturb = vturb                             #The turbulence velocity blocked by planet.
                Io_planet = 0.0

                #Start the planet on the x-axis so that the planet is just far enough away not to touch the disk of the star.
                if Distance_center <= (Rs + Rp) and Distance_center >= (Rs - Rp):

                    dist_cent1_int = (Distance_center**2.0 + Rs2 - Rp2)/(2.0*Distance_center)
                    #Location on the x-axis for the first intersection point.
                    X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) \
                                  + ((Ypos/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))
                    #Location on the y-axis for the first intersection point.
                    Y_int_1 = ((Ypos*dist_cent1_int)/Distance_center) \
                                  - ((Xpos/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))
                    #Location on the x-axis for the second intersection point.
                    X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) \
                                  - ((Ypos/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))
                    #Location on the y-axis for the second intersection point.
                    Y_int_2 = ((Ypos*dist_cent1_int)/Distance_center) \
                                  + ((Xpos/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))

                    #The limb darkening equation will use this distance if any part of the disc of the planet is outside
                    #the radius of the star.
                    #This is the distance between the center of the star to the center of the area inside the star
                    #blocked by the planet.
                    set_distance_center = ((math.fabs(Distance_center) - Rp) + Rs)/2.0
                    #Next the program calculates the length between the two intersection point.
                    Length = math.sqrt((X_int_1 - X_int_2)**2.0 + (Y_int_1 - Y_int_2)**2.0)
                    #Calculate the angle between the Y position of the center of planet and the x position of the center
                    #of planet.
                    #This is used to determine the center of the area inside the star blocked by the planet in the
                    #stellar rotational axis coordinate system.
                    Theta_inside = math.atan(Ypos/Xpos)
                    #Calculates the distance from the center of the star (along x axis) to the center of the area inside
                    #the star blocked by the planet.
                    #the pi factor added in this equation guarantees that Xpos_in has the correct sign.
                    Xpos_in = set_distance_center*math.cos(Theta_inside + pi)
                    Xpos_inside = Xpos_in

                    #This makes sure Xpos_inside has the correct sign.
                    if Xpos >= 0:
                            Xpos_inside = -Xpos_in
                    #endif

                    #Calculates the distance from the center of the star (along y axis) to the center of the area inside
                    #the star blocked by the planet.
                    Ypos_inside = math.fabs(set_distance_center*math.sin(Theta_inside))
                    #Changes the x-coordinate to the stellar rotation axis by an angle formed between the orbital plane
                    #of the planet and the stellar roatation plane of the star.
                    x_prime_distance = (Xpos_inside*math.cos((stellar_rotation_angle*pi)/180.0)) \
                                           + (Ypos_inside*math.sin((stellar_rotation_angle*pi)/180.0))
                    Sub_planet_velocity = vsini*(x_prime_distance/Rs)    #Calculate the subplanetary velocity (the
                                                                         #stellar velocity blocked by the planetary
                    #disk) from vsini times the distance from the center of the planet to the stellar rotation axis
                    #divided by the radius of the star.
                    Beta1 = 2.0*math.asin((0.5*Length)/Rp)         #Angle used to calculate the partial area of the
                                                                   #planet in of the star.
                    Alpha1 = 2.0*math.asin((0.5*Length)/Rs)        #Angle used to calculate the partial area of the
                                                                   #planet in of the star.

                    if Distance_center >= math.sqrt(Rs2 - Rp2):
                        #The surface area of the planet when the center is outside the disk of the star.
                        Aplan = ((0.5 * Rs2 * (Alpha1 - math.sin(Alpha1))) + (0.5 * Rp2 * (Beta1 - math.sin(Beta1))))
                        #Normalized area of the planet with the given limb darkening function
                        Io_planet = Io * Aplan
                    #endif

                    if Distance_center < math.sqrt(Rs2 - Rp2):
                        Aplan = ((pi*Rp2 + (0.5 * Rs2 * (Alpha1 - math.sin(Alpha1)))) -
                                    (0.5 * Rp2 * (Beta1 - math.sin(Beta1))))
                        #The surface area of the planet when the center is inside the disk of the star.
                        Io_planet = Io * Aplan
                    #endif
                #endif

                if Distance_center < (Rs - Rp):
                    Aplan = (pi*Rp2)                        #The surface area of the planet.
                    Io_planet = Io * Aplan
                #endif

                #The ratio of the area of the planet blocking the star to the area of the star utilizing the first order
                #limb darkening equation.
                if linear_quadratic == 'q':
                    #Quadratic limb darkening equation.
                    Lblocked = Io_planet*(1.0-q_1*(1.0-math.sqrt(math.fabs(1.0-(set_distance_center**2.0/Rs2)))) -
                                   q_2*(1.0 - math.sqrt(math.fabs(1.0-(set_distance_center**2.0/Rs2))))**2.0)
                else:
                    #linear limb darkening equation.
                    Lblocked = (Io_planet)*(1.0-(u*(1.0-math.sqrt(math.fabs(1.0-(set_distance_center**2.0/Rs2))))))
                #endelse
                if test_null_hypothesis:
                    v_rm = 0.0
                else:
                    v_rm = - ((Lblocked*Sub_planet_velocity)*((((2.0*vturb**2.0)+(2.0*vsini**2.0))/((2.0*vturb**2.0)
                           + vsini**2.0))**(3.0/2.0))*(1.0-((Sub_planet_velocity**2.0)/((2.0*vturb**2.0)+(vsini**2.0)))))
                Total_L = 1.0 - Lblocked                      #Total amount of light blocked by the planet.
                Total_RM = 0.0 + v_rm                         #Total anomalous velocity.

            #endif

        #endif

        #If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is negative then the
        #secondary transit (occulation) begins.
        if (Distance_center <= (Rs + Rp)) and (Zpos < 0):

            if Time_occultation_start == 0:
                Time_occultation_start = Time_ref           #Indicates when the secondary transit or occultation starts.
                                                            #If statement prevents this variable from being overwritten
                                                            #every time the loop runs.
                Occultation_start_position = Time_loop
            #endif

            if Distance_center >= (Rs - Rp):           #Partial secondary transit.
                dist_cent1_int = (Distance_center**2.0 + Rs2 - Rp2)/(2.0*Distance_center)
                Ypos2 = math.fabs(Ypos)   #Set Ypos to a positive value to avoid discontinuity.
                #Location on the x-axis for the first intersection point.
                X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) + \
                          ((Ypos2/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))
                #Location on the y-axis for the first intersection point.
                Y_int_1 = ((Ypos2*dist_cent1_int)/Distance_center) - \
                          ((Xpos/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))
                #Location on the x-axis for the second intersection point.
                X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) - \
                          ((Ypos2/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))
                #Location on the y-axis for the second intersection point.
                Y_int_2 = ((Ypos2*dist_cent1_int)/Distance_center) + \
                          ((Xpos/Distance_center)*math.sqrt(Rs2 - dist_cent1_int**2.0))
                #The limb darkening equation will use this distance if any part of the disc of the planet is outside the
                #radius of the star. This is the distance between the center of the star to the center of the area
                #inside the planet blocked by the star.
                set_distance_center = ((math.fabs(Distance_center) - Rp) + Rs)/2.0
                #Next the program calculates the length between the two intersection point.
                Length = math.sqrt((X_int_1 - X_int_2)**2.0 + (Y_int_1 - Y_int_2)**2.0)
                #Calculate the angle between the Y position of the center of planet and the x position of the center of
                #planet.
                Theta_inside = math.atan(Ypos2/Xpos)
                #Calculates the distance from the center of the star (along x axis) to the center of the area inside the
                #star blocked by the planet.
                #The pi factor added in this equation guarneetes that Xpos_in has the correct sign.
                Xpos_in = set_distance_center*math.cos(Theta_inside + pi)
                Xpos_inside = Xpos_in

                #This makes sure Xpos_inside has the correct sign.
                if Xpos >= 0:
                    Xpos_inside = -Xpos_in
                #endif

                #Calculates the distance from the center of the star (along y axis) to the center of the area inside the
                #star blocked by the planet.
                Ypos_inside = set_distance_center*math.sin(Theta_inside)
                #Angle used to calculate the partial area of the planet in the star.
                Beta1 = 2.0*math.asin((0.5*Length)/Rp)
                #Angle used to calculate the partial area of the planet in the star.
                Alpha1 = 2.0*math.asin((0.5*Length)/Rs)

                if Distance_center >= math.sqrt(Rs2 - Rp2):
                    #The surface area of the planet that is not behind the star when the center of the planets disk is
                    #visible.
                    Aplan = (pi*Rp2) - ((0.5*Rp2*(Beta1 - math.sin(Beta1))) + (0.5*Rs2*(Alpha1 - math.sin(Alpha1))))
                #endif

                if Distance_center < math.sqrt(Rs2 - Rp2):
                    #The surface area of the planet that is not behind the star when the center of the planets disk is
                    #behind the star.
                    Aplan = 0.5 * Rp2 * (Beta1 - math.sin(Beta1))
                #endif
            #endif

            if Distance_center < (Rs - Rp):
                Aplan = 0.0                      #No flux coming from planet since it's behind the star.
            #endif

            #Total amount of light blocked by the planet.
            Total_L = 1.0 + ((Aplan/(4.0*pi*Planet_star_distance**2.0))*Albedo*(0.5*(1.0 + math.cos(True_phase))))
        #endif

        if Total_L >= 1.001:
            print('Total_L: ', Total_L)
            print('Aplan: ', Aplan)
            print('pi*Rp2: ', pi*Rp2)
            print('Length: ', Length)
            print('Alpha1: ', Alpha1)
            print('Beta1: ', Beta1)
            print('Beta1 - math.sin(Beta1): ', Beta1 - math.sin(Beta1))
            print('((0.5*Rp2*(Beta1 - math.sin(Beta1))) + (0.5*Rs2*(Alpha1 - math.sin(Alpha1)))): ',
                  ((0.5*Rp2*(Beta1 - math.sin(Beta1))) + (0.5*Rs2*(Alpha1 - math.sin(Alpha1)))))
            print('Planet_star_distance: ', Planet_star_distance)
            print('True_phase: ', True_phase)
            print('math.cos(True_phase): ', math.cos(True_phase))
            print('Time_ref: ', Time_ref)
            print('            ')
        #endif

        #print("Test print")

        if Time_transit_start != 0 and Time_transit_end == 0 and Distance_center >= (Rs + Rp):
            Time_transit_end = Time_ref               #Indicates when the transit Ends. If statement prevents this
                                                      #variable from being overwritten every time the loop runs.
            transit_End_position = Time_loop
        #endif

        if Time_mid_transit == 0 and Xpos >= 0:
            Time_mid_transit = Time_ref         #Indicates when the mid transit point occurs. If statement prevents
                                                #this variable from being overwritten every time the loop runs.
            transit_mid_position = Time_loop
        #endif

        if Time_mid_occultation == 0 and Xpos <= 0 and Time_mid_transit != 0:
            Time_mid_occultation = Time_ref      #Indicates when the mid occultation point occurs. If statement prevents
                                                 #this variable from being overwritten every time the loop runs.
            occultation_mid_position = Time_loop
        #endif

        if Transit_start_no_inc != 0 and Transit_end_no_inc == 0 and math.fabs(Xpos) >= (Rs + Rp):
            Transit_end_no_inc = Time_ref        #Indicates when the transit Ends @ 90 inc. If statement prevents this
                                                 #variable from being overwritten every time the loop runs.
            transit_End_position_no_inc = Time_loop
        #endif

        if Time_occultation_start != 0 and Time_occultation_end == 0 and Distance_center >= (Rs + Rp):
            Time_occultation_end = Time_ref        #Indicates when the transit Ends. If statement prevents this variable
                                                   #from being overwritten every time the loop runs.
            occultation_End_position = Time_loop
        #endif

        if Ecc == 0:
            #The radial velocity of the star which includes adding the RM effect for a circular orbit.
            RV = RVamplitude*(math.cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0)))) + Total_RM
        else:
            #The radial velocity of the star which includes adding the RM effect for an eccentric orbit.
            RV = RVamplitude*(math.cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0)))
                     + (Ecc*math.cos((omega_arg_periastron)*(pi/180.0) + pi))) + Total_RM
        #endif

        #Put data into the Lightcurve array from the total luminosity.
        Transit_LC_array.loc[Time_loop, "Flux"] = Total_L
        Transit_LC_array.loc[Time_loop, "Time"] = Time_ref
        #Put all the anomalous velocity data into the anomalous velocity array.
        RM_effect_array.loc[Time_loop, "RV"] = Total_RM
        RM_effect_array.loc[Time_loop, "Time"] = Time_ref
        #Put all the RV data into the RV array.
        RV_array.loc[Time_loop, "RV"] = RV
        RV_array.loc[Time_loop, "Time"] = Time_ref
        #Set timestep to start at negative time.  The mid point of transit is at zero and everything after that is
        #positive time in seconds.
        Timestep.loc[Time_loop, "Time"] = Time_ref

    #endfor




    Time_transit_actual = Time_transit_end - Time_transit_start
    print("Time_transit_actual: ", Time_transit_actual)
    #The length of the occultation which cannot be determined ahead of time if the planet has an eccentric orbit.
    Time_occultation_actual = Time_occultation_end - Time_occultation_start
    print("Time_occultation_actual: ", Time_occultation_actual)
    Time_transit_90inc = Transit_end_no_inc - Transit_start_no_inc
    print("Time_transit_90inc: ", Time_transit_90inc)

    min_LC_value = Transit_LC_array["Flux"].min()
    print("min_LC_value: ", min_LC_value)
    amp_RM_effect = RM_effect_array["RV"].max() + math.fabs(RM_effect_array["RV"].min())
    print("amp_RM_effect: ", amp_RM_effect)
    max_frac_flux_decrease = 1.0 - Transit_LC_array["Flux"].min()
    print("max_frac_flux_decrease: ", max_frac_flux_decrease)
    RV_full_amp = RV_array["RV"].max() + math.fabs(RV_array["RV"].min())
    print("RV_full_amp: ", RV_full_amp)
    # os.system("pause")
    max_flux_planet = Transit_LC_array["Flux"].max() - 1.0
    print("max_flux_planet: ", max_flux_planet)




    return Transit_LC_array, RM_effect_array, RV_array, Timestep, Time_transit_actual, Time_occultation_actual, \
           Time_transit_90inc, min_LC_value, amp_RM_effect, max_frac_flux_decrease, RV_full_amp, max_flux_planet, \
           Transit_start_no_inc, Time_mid_transit, Transit_end_no_inc, Time_occultation_start, Time_occultation_end

#enddef
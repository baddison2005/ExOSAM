PROGRAM RM_effect_residuals
    
!This program determines the residuals between the model from the best fit spin-orbit angle and vsini (as determined
!from the MCMC routine) and the observed RVs. The priors on the other parameters are used to fit the model and 
!calculate the residuals. This is because the best fit values for the other parameters that were determined from
!fitting the RM effect are less accurate than being determined from other types of data (e.g., transit light curves,
!RVs taken over the planets orbit, etc.).

implicit none

	!   Constant and variable decleration

double precision, parameter :: Rss = 6.9634D8                !Radius of our Sun (in meters).
double precision, parameter :: Rj = 7.1492D7                 !Radius of Jupiter (in meters).
double precision, parameter :: RE = 6.3781D6                 !Radius of Earth (in meters if planet is given in Earth radii).
double precision, parameter :: AU = 1.4960D11                !One astronomical unit in meters.
double precision, parameter :: pi = 3.1415926535           !Constant value of pi.
double precision, parameter :: G = 6.6738D-11                !Gravitation constant.
double precision, parameter :: Mss = 1.9891D30               !Mass of the sun in kilograms.
double precision, parameter :: Me = 5.9722D24                !Mass of the Earth in kg.
double precision, parameter :: Mj = 1.8986D27                !Mass of Jupiter in kg.
double precision, parameter :: day_sec = 86400.0D0         !Length of a day in seconds.

double precision :: a, b, c         !variables used to read in data
double precision :: Albedo
double precision :: Alpha1
double precision :: Aplan
double precision :: Area_pixel
double precision :: Bessel_function_exit
double precision :: Bessel_value
double precision :: best_spin_orbit_angle
double precision :: best_vsini
double precision :: Beta1
double precision :: Center_of_planet_x
double precision :: Center_of_planet_y
double precision :: chi_squared_change          !Check whether this variable will remain. May change file read in.
double precision :: chi_squared_change_fit
double precision :: data_plot_model_interval
double precision :: Dist2
double precision :: Distance_center
double precision :: Dist_center_pixel
double precision :: dist_cent1_int
double precision :: Dist_planet_pixel
double precision :: Ecc
double precision :: ecc_anomaly
double precision :: ecc_anomaly_after
double precision :: ecc_anomaly_before
double precision :: ecc_anomaly_start
double precision :: ecc_anomaly_transit
double precision :: Ecc_begin
double precision :: Ecc_end
double precision :: Ecc_prior
double precision :: Ecc_1sigerr
double precision :: impact_prior
double precision :: Inc
double precision :: Inc_begin
double precision :: Inc_end
double precision :: Inc_prior
double precision :: Inc_1sigerr
double precision :: Io
double precision :: Io_Pixel
double precision :: Io_planet
double precision :: JD_time_peri
double precision :: JD_time_mid_transit
double precision :: JD_time_mid_transit_begin
double precision :: JD_time_mid_transit_end
double precision :: JD_time_mid_transit_prior
double precision :: JD_time_mid_transit_1sigerr
double precision :: Lblocked
double precision :: Lblocked2
double precision :: Length
double precision :: length_time_compare
double precision :: Mean_anomaly_transit
double precision :: model_data_difference
double precision :: Mp
double precision :: Mp_begin
double precision :: Mp_end
double precision :: Mp_prior
double precision :: Mp_1sigerr
double precision :: Ms
double precision :: Ms_solar_prior                                !MCMC variables.
double precision :: Ms_solar_end                                  !MCMC variables.
double precision :: Ms_solar_begin                                !MCMC variables.
double precision :: Ms_solar_1sigerr                              !MCMC variables.
double precision :: Number_orbits
double precision :: omega_arg_periastron
double precision :: omega_arg_periastron_begin
double precision :: omega_arg_periastron_end
double precision :: omega_arg_periastron_prior
double precision :: omega_arg_periastron_1sigerr
double precision :: Orbital_period
double precision :: Orbital_period_1
double precision :: Orbital_period_begin
double precision :: Orbital_period_end
double precision :: Orbital_period_prior
double precision :: Orbital_period_1sigerr
double precision :: Phase_angle
double precision :: Phase_angle_observed
double precision :: Phase_angle_start
double precision :: Phase_orbit_n
double precision :: Pixel
double precision :: Planet_star_distance
double precision :: q_1
double precision :: q_2
double precision :: Radius_planet_array
double precision :: radius_ratio_prior
double precision :: Rorb
double precision :: Rorb_prior
double precision :: Rorb_star
double precision :: Rorb_star_prior
double precision :: Rp
double precision :: Rp_begin
double precision :: Rp_end
double precision :: Rp_prior
double precision :: Rp2
double precision :: Rp_1sigerr
double precision :: Rp_Rs_ratio_prior                             !MCMC variables.
double precision :: Rp_Rs_ratio_end                               !MCMC variables.
double precision :: Rp_Rs_ratio_begin                             !MCMC variables.
double precision :: Rp_Rs_ratio_1sigerr                           !MCMC variables.
double precision :: Rs
double precision :: Rs2
double precision :: Rs_solar_prior                                !MCMC variables.
double precision :: Rs_solar_end                                  !MCMC variables.
double precision :: Rs_solar_begin                                !MCMC variables.
double precision :: Rs_solar_1sigerr                              !MCMC variables.
double precision :: RV
double precision :: RVamplitude
double precision :: RV_offset_datasets_begin
double precision :: RV_offset_datasets_end
double precision :: RV_offset_datasets_interval
double precision :: RV_offset_datasets_prior
double precision :: RV_offset_datasets_1sigerr
double precision :: RV_zero_offset
double precision :: RV_zero_offset_mean_mcmc
double precision :: RV_zero_offset_begin
double precision :: RV_zero_offset_end
double precision :: RV_zero_offset_interval
double precision :: RV_zero_offset_prior
double precision :: RV_zero_offset_1sigerr
double precision :: scale_factor_in
double precision :: set_distance_center
double precision :: stellar_rotation_angle_begin
double precision :: stellar_rotation_angle_end
double precision :: stellar_rotation_angle_prior
double precision :: stellar_rotation_angle_mean_mcmc
double precision :: stellar_rotation_angle_1sigerr
double precision :: Sub_planet_velocity
double precision :: sum_ecc_anomaly
double precision :: time
double precision :: Time_bef_aft
double precision :: Time_check
double precision :: Time_compare_vel
double precision :: time_peri_passage
double precision :: Time_ref
double precision :: Time_start
double precision :: Time_transit
double precision :: Theta_inside
double precision :: Transit_length_prior
double precision :: True_anomaly
double precision :: True_anomaly_start
double precision :: True_anomaly_transit
double precision :: True_phase
double precision :: Total_L
double precision :: Total_RM
double precision :: u
double precision :: v_rm
double precision :: vmacro
double precision :: beta_rot
double precision :: vsini_mean_mcmc
double precision :: vsini_begin
double precision :: vsini_end
double precision :: vsini_prior
double precision :: vsini_1sigerr
double precision :: vturb
double precision :: X_int_1
double precision :: X_int_2
double precision :: X_pixel
double precision :: X_pixel2
double precision :: X_pixel_prime
double precision :: X_prime
double precision :: x_prime_distance
double precision :: Xpos
double precision :: Xpos_in
double precision :: Xpos_inside
double precision :: XpXp
double precision :: Y_int_1
double precision :: Y_int_2
double precision :: Y_pixel
double precision :: Ypos
double precision :: Ypos2
double precision :: Ypos_inside
double precision :: Zpos




integer(8) :: num_compare = 0
integer(8) :: num_compare_out
integer(8) :: Number = 1               !Counter.
integer(8) :: Number_points1 = 1                           !The number of points that are used from the predicted light curve.
integer(8) :: mcmc_accepted_iteration_size
integer(8) :: number_mcmc_walkers
integer(8) :: M_pixels
integer(8) :: datafilelength
integer(8) :: total_interval
integer(8) :: num_my_rv
integer(8) :: num_other_rv

!Loop variables
integer(8) :: order, Time_loop
integer(8) :: i, l, j, s, f, k, zz




character(LEN=1) :: add_RV_zero_offset
character(LEN=1) :: Jupiter_Earth_units
character(LEN=1) :: Phase_shift
character(LEN=1) :: time_compare_flag = 'N'
character(LEN=1) :: transit_flag
character(LEN=1) :: occultation_flag
character(LEN=1) :: other_RV_files
character(LEN=1) :: subtract_JD_mid_time
character(LEN=60) :: current_directory
character(LEN=150) :: input_RM_filename
character(LEN=150) :: input_my_data_filename
character(LEN=150) :: input_other_data_filename
character(LEN=1) :: Ms_solar_norm_prior
character(LEN=1) :: Ms_solar_fixed
character(LEN=1) :: Rs_solar_norm_prior
character(LEN=1) :: Rs_solar_fixed
character(LEN=1) :: Rp_Rs_ratio_norm_prior
character(LEN=1) :: Rp_Rs_ratio_fixed
character(LEN=1) :: use_Rp_Rs_ratio
character(LEN=1) :: Inc_norm_prior
character(LEN=1) :: Inc_fixed
character(LEN=1) :: omega_arg_periastron_norm_prior
character(LEN=1) :: omega_arg_periastron_fixed
character(LEN=1) :: Ecc_norm_prior
character(LEN=1) :: Ecc_fixed
character(LEN=1) :: Mp_norm_prior
character(LEN=1) :: Mp_fixed
character(LEN=1) :: Rp_norm_prior
character(LEN=1) :: Rp_fixed
character(LEN=1) :: Orbital_period_norm_prior
character(LEN=1) :: Orbital_period_fixed
character(LEN=1) :: JD_time_mid_transit_norm_prior
character(LEN=1) :: JD_time_mid_transit_fixed
character(LEN=1) :: RV_zero_offset_norm_prior
character(LEN=1) :: RV_zero_offset_fixed
character(LEN=1) :: RV_offset_datasets_norm_prior
character(LEN=1) :: RV_offset_datasets_fixed
character(LEN=1) :: vsini_norm_prior
character(LEN=1) :: vsini_fixed
character(LEN=1) :: stellar_rotation_angle_norm_prior
character(LEN=1) :: impose_prior_stellar_rotation_angle
character(LEN=1) :: stellar_rotation_angle_fixed
character(LEN=20) :: Planet_name
character(LEN=150) :: output_temp_data
character(LEN=150) :: directory_location
character(LEN=150) :: output_data_fit_array_filename
character(LEN=150) :: output_RV_theory_fit_array_filename
character(LEN=150) :: output_RV_theory_all_array_filename
character(LEN=150) :: output_residual_all_array_filename
character(LEN=150) :: output_residual_array_filename
character(LEN=150) :: input_best_parameters_mcmc_filename
character(LEN=60) :: input_data_filename
character(LEN=150) :: my_data_plotting_symbols
character(LEN=150) :: other_data_plotting_symbols
character(LEN=1) :: linear_quadratic
character(LEN=1) :: impose_prior_vsini
character(LEN=1) :: use_out_transit_rv_for_fit




DOUBLE PRECISION, Allocatable :: Data(:,:), Data_my_rv(:,:), Data_adjusted(:,:), RV_theory(:,:), sorted_data_array(:,:)
DOUBLE PRECISION, Allocatable :: Data_my_rv_offset(:,:), time_data_array(:), Residuals_array(:,:), data_in_fit(:,:)
DOUBLE PRECISION, Allocatable :: RV_Theory_fit(:,:), Residuals_all_array(:,:), RV_Theory_all(:,:), RV_offset_data_array(:,:)
DOUBLE PRECISION, Allocatable :: Data_other_rv(:,:)
INTEGER(8), DIMENSION(:), ALLOCATABLE :: index_array
CHARACTER(len=1), Allocatable :: data_points_compare(:)




!Get the locations of files.
CALL get_environment_variable('PWD', current_directory)

current_directory = TRIM(current_directory)

OPEN(unit=99, FILE=TRIM(TRIM(current_directory)  // '/pointer.txt'), status='old', action='read')
READ(99,46,end=47) input_data_filename
46 FORMAT(A150)
47 CLOSE(99)

input_data_filename = TRIM(ADJUSTL(input_data_filename))
PRINT *, 'input_data_filename = ', input_data_filename




!Read in output data from IDL and assaign values to parameters and variables.
OPEN(unit=99, FILE=TRIM(ADJUSTL(input_data_filename)), status='old', action='read')
!Now read line by line of data file and assaign values to parameters and variable.
READ(99,49,end=50) input_RM_filename, input_my_data_filename, input_other_data_filename, num_my_rv, num_other_rv, &
datafilelength, directory_location, Planet_name, other_RV_files, output_temp_data, &
Jupiter_Earth_units, Bessel_function_exit, Phase_shift, subtract_JD_mid_time, add_RV_zero_offset, &
Number_orbits, data_plot_model_interval, Phase_angle_start, Time_bef_aft, Time_compare_vel, &
my_data_plotting_symbols, other_data_plotting_symbols, Ms_solar_prior, Ms_solar_end, Ms_solar_begin, &
Ms_solar_1sigerr, Ms_solar_norm_prior, Ms_solar_fixed, Rs_solar_prior, Rs_solar_end, &
Rs_solar_begin, Rs_solar_1sigerr, Rs_solar_norm_prior, Rs_solar_fixed, Rp_Rs_ratio_prior, &
Rp_Rs_ratio_end, Rp_Rs_ratio_begin, Rp_Rs_ratio_1sigerr, Rp_Rs_ratio_norm_prior, Rp_Rs_ratio_fixed, &
use_Rp_Rs_ratio, vmacro, beta_rot, M_pixels, linear_quadratic, &
u, q_1, q_2, Inc_prior, Inc_end, &
Inc_begin, Inc_1sigerr, Inc_norm_prior, Inc_fixed, omega_arg_periastron_prior, &
omega_arg_periastron_end, omega_arg_periastron_begin, omega_arg_periastron_1sigerr, omega_arg_periastron_norm_prior, omega_arg_periastron_fixed, &
Ecc_prior, Ecc_end, Ecc_begin, Ecc_1sigerr, Ecc_norm_prior, &
Ecc_fixed, Mp_prior, Mp_end, Mp_begin, Mp_1sigerr, &
Mp_norm_prior, Mp_fixed, Rp_prior, Rp_end, Rp_begin, &
Rp_1sigerr, Rp_norm_prior, Rp_fixed, Orbital_period_prior, Orbital_period_end, &
Orbital_period_begin, Orbital_period_1sigerr, Orbital_period_norm_prior, Orbital_period_fixed, JD_time_mid_transit_prior, &
JD_time_mid_transit_end, JD_time_mid_transit_begin, JD_time_mid_transit_1sigerr, JD_time_mid_transit_norm_prior, JD_time_mid_transit_fixed, &
Albedo, RV_zero_offset_prior, RV_zero_offset_end, RV_zero_offset_begin, RV_zero_offset_interval, &
RV_zero_offset_1sigerr, RV_zero_offset_norm_prior, RV_zero_offset_fixed, RV_offset_datasets_prior, RV_offset_datasets_end, &
RV_offset_datasets_begin, RV_offset_datasets_interval, RV_offset_datasets_1sigerr, RV_offset_datasets_norm_prior, RV_offset_datasets_fixed, &
vsini_prior, vsini_end, vsini_begin, vsini_1sigerr, vsini_norm_prior, &
impose_prior_vsini, vsini_fixed, stellar_rotation_angle_prior, stellar_rotation_angle_end, stellar_rotation_angle_begin, &
stellar_rotation_angle_1sigerr, stellar_rotation_angle_norm_prior, impose_prior_stellar_rotation_angle, stellar_rotation_angle_fixed, mcmc_accepted_iteration_size, &
number_mcmc_walkers, scale_factor_in, chi_squared_change, chi_squared_change_fit, use_out_transit_rv_for_fit
49 FORMAT(A150, /, A150, /, A150, /, I20, /, I20, /, &
          I20, /, A150, /, A20, /, A1, /, A150, /, &
          A1, /, E50.10, /, A1, /, A1, /, A1, /, &
          F10.5, /, F50.10, /, F50.5, /, F50.10, /, F50.10, /, &
          A150, /, A150, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, A1, /, A1, /, F50.5, /, F50.5, /, &
          F50.5, /, F50.5, /, A1, /, A1, /, F50.5, /, &
          F50.5, /, F50.5, /, F50.5, /, A1, /, A1, /, &
          A1, /, F50.5, /, F50.5, /, I20, /, A1, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, F50.5, /, A1, /, A1, /, F50.5, /, &
          F50.5, /, F50.5, /, F50.5, /, A1, /, A1, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, A1, /, &
          A1, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          A1, /, A1, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, A1, /, A1, /, F50.10, /, F50.10, /, &
          F50.10, /, F50.10, /, A1, /, A1, /, F50.10, /, &
          F50.10, /, F50.10, /, F50.10, /, A1, /, A1, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, A1, /, A1, /, F50.5, /, F50.5, /, &
          F50.5, /, F50.5, /, F50.5, /, A1, /, A1, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, A1, /, &
          A1, /, A1, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, A1, /, A1, /, A1, /, I20, /, &
          I20, /, F50.5, /, F50.5, /, F50.5, /, A1)
50 CLOSE(99)




input_RM_filename = ADJUSTL(input_RM_filename)
input_RM_filename = TRIM(input_RM_filename)
PRINT *, 'input_RM_filename = ', input_RM_filename

input_my_data_filename = ADJUSTL(input_my_data_filename)
input_my_data_filename = TRIM(input_my_data_filename)
PRINT *, 'input_my_data_filename = ', input_my_data_filename

input_other_data_filename = ADJUSTL(input_other_data_filename)
input_other_data_filename = TRIM(input_other_data_filename)
PRINT *, 'input_other_data_filename = ', input_other_data_filename

directory_location = ADJUSTL(directory_location)
directory_location = TRIM(directory_location)
PRINT *, 'directory_location = ', directory_location
Planet_name = ADJUSTL(Planet_name)
Planet_name = TRIM(Planet_name)
PRINT *, 'Planet_name = ', Planet_name
output_temp_data = ADJUSTL(output_temp_data)
output_temp_data = TRIM(output_temp_data)
PRINT *, 'output_temp_data = ', output_temp_data

!MCMC file names. ****************************************************************************
output_residual_array_filename = TRIM(TRIM(output_temp_data) // 'residual_array.txt')
output_data_fit_array_filename = TRIM(TRIM(output_temp_data) // 'data_fit_array.txt')
output_RV_theory_fit_array_filename = TRIM(TRIM(output_temp_data) // 'RV_Theory_fit_array.txt')
output_residual_all_array_filename = TRIM(TRIM(output_temp_data) // 'residual_all_array.txt')
output_RV_theory_all_array_filename = TRIM(TRIM(output_temp_data) // 'RV_Theory_all_array.txt')
input_best_parameters_mcmc_filename = TRIM(TRIM(output_temp_data) // 'best_parameters_mcmc.txt')




Allocate(Data_other_rv(num_other_rv,3))
Allocate(Data_my_rv(num_my_rv,3), time_data_array(datafilelength))
Allocate(RV_offset_data_array(datafilelength,3))
Allocate(Data_adjusted(datafilelength,3), sorted_data_array(datafilelength,3))




OPEN(unit=99, FILE=input_best_parameters_mcmc_filename, status='old', action='read')
READ(99,245) vsini_mean_mcmc
READ(99,246) stellar_rotation_angle_mean_mcmc
READ(99,247) RV_zero_offset_mean_mcmc
245  FORMAT(/, /, /, /, /, /, /, 92X, F15.6)
246  FORMAT(/, /, 92X, F15.6)
247  FORMAT(/, /, /, 92X, F15.6)
CLOSE(99)

PRINT *,  'vsini_mean_mcmc: ', vsini_mean_mcmc
PRINT *,  'stellar_rotation_angle_mean_mcmc: ', stellar_rotation_angle_mean_mcmc




!Read in output data from IDL and assaign values to the data array.
OPEN(unit=99, FILE=input_my_data_filename, status='old', action='read')
!Now read line by line of data file and dump values into the data array.

DO l = 1, num_my_rv
    READ (99, 53) a, b, c
    53 FORMAT(E50.10, 5X, F50.5, 5X, F50.5)
    Data_my_rv(l,1) = a
    Data_my_rv(l,2) = b
    Data_my_rv(l,3) = c
END DO 

CLOSE(99)

IF (other_RV_files == 'Y') THEN
    !Read in output data from IDL and assaign values to the data array.
    OPEN(unit=99, FILE=input_other_data_filename, status='old', action='read')
    !Now read line by line of data file and dump values into the data array.

    DO l = 1, num_other_rv
        READ (99, 54) a, b, c
        54 FORMAT(E50.10, 5X, F50.5, 5X, F50.5)
        Data_other_rv(l,1) = a
        Data_other_rv(l,2) = b
        Data_other_rv(l,3) = c
    END DO 

    CLOSE(99)
END IF




!Now apply the RV offset to my data.
Allocate(Data_my_rv_offset(num_my_rv,3), data(datafilelength,3))

DO l = 1, num_my_rv
    Data_my_rv_offset(l,1) = Data_my_rv(l,1)
    Data_my_rv_offset(l,2) = Data_my_rv(l,2) + (RV_zero_offset_mean_mcmc - RV_zero_offset_prior) + RV_offset_datasets_prior
    Data_my_rv_offset(l,3) = Data_my_rv(l,3)
END DO

PRINT *, "RV_zero_offset_mean_mcmc - RV_zero_offset_prior: ", RV_zero_offset_mean_mcmc - RV_zero_offset_prior
PRINT *, "RV_zero_offset_mean_mcmc: ", RV_zero_offset_mean_mcmc
PRINT *, "RV_zero_offset_prior: ", RV_zero_offset_prior

IF (other_RV_files == 'Y') THEN
    DO l = 1, datafilelength
        IF (l <= num_other_rv) THEN
            Data(l,1) = Data_other_rv(l,1)
            Data(l,2) = Data_other_rv(l,2)
            Data(l,3) = Data_other_rv(l,3)
        ELSE
            zz = l - num_other_rv
            Data(l,1) = Data_my_rv_offset(zz,1)
            Data(l,2) = Data_my_rv_offset(zz,2) + (RV_zero_offset_mean_mcmc - RV_zero_offset_prior)
            Data(l,3) = Data_my_rv_offset(zz,3)
        END IF
    END DO
ELSE
    DO l = 1, datafilelength
        Data(l,1) = Data_my_rv_offset(l,1)
        Data(l,2) = Data_my_rv_offset(l,2)
        Data(l,3) = Data_my_rv_offset(l,3)
    END DO
END IF

!**********************************************************************************************************************************************
!Based on given priors, determine the approximate time in the simulation to start and finish comparing RV's with model.
Ms = Ms_solar_prior * Mss
Rs = Rss * Rs_solar_prior

IF (use_Rp_Rs_ratio == 'Y') THEN
    Rp = Rp_Rs_ratio_prior * (Rs_solar_prior * Rss)
ELSE
    IF (Jupiter_Earth_units == 'Y') THEN
        Rp = Rp_prior*Rj
    END IF

    IF (Jupiter_Earth_units == 'N') THEN
        Rp = Rp_prior*RE
    END IF
END IF

IF (Jupiter_Earth_units == 'Y') THEN
    Mp = Mp_prior*Mj
END IF

IF (Jupiter_Earth_units == 'N') THEN
    Mp = Mp_prior*Me
END IF

Rorb_prior = (((Orbital_period_prior*day_sec)**2.0d0*G*(Ms + Mp))/(4.0d0*pi**2.0d0))**(1.0d0/3.0d0)     !Semi-major axis.
Rorb_star_prior = (((Orbital_period_prior*day_sec)**2.0D0*G*((Mp**(3.0D0))/(Ms + Mp)**(2.0D0)))/(4.0d0*pi**2.0D0))**(1.0d0/3.0d0)     !Semi-major of stellar axis.
radius_ratio_prior = Rp_Rs_ratio_prior
impact_prior = ((Rorb_prior*cos(Inc_prior*(pi/180.0D0)))/Rs)*((1.0D0 - Ecc_prior**2.0D0)/(1.0D0 + (Ecc_prior*sin(omega_arg_periastron_prior*(pi/180.0D0)))))
Transit_length_prior = ((Orbital_period_prior*day_sec)/pi) * asin((Rs/Rorb_prior)*(sqrt(abs((1.0D0 + radius_ratio_prior)**2.0D0 - impact_prior**2.0D0))/sin(Inc_prior*(pi/180.0D0)))) * &
                       (sqrt(1.0D0 - Ecc_prior**2.0D0)/(1.0D0 + (Ecc_prior*sin(omega_arg_periastron_prior*(pi/180.0D0)))))
                       
PRINT *, "Length of transit with given priors: ", Transit_length_prior
length_time_compare = (Transit_length_prior/2.0D0) + Time_compare_vel
PRINT *, "Length of time to compare: ", length_time_compare




!**********************************************************************************************************************************************
best_vsini = vsini_mean_mcmc
best_spin_orbit_angle = stellar_rotation_angle_mean_mcmc

Orbital_period = Orbital_period_prior*day_sec               !Convert orbital period in days to seconds.
Orbital_period_1 = Orbital_period_prior
PRINT *, "Orbital_period seconds: ", Orbital_period
PRINT *, "Orbital_period days: ", Orbital_period_1

Rorb = ((Orbital_period**2.0d0*G*(Ms + Mp))/(4.0d0*pi**2.0d0))**(1.0d0/3.0d0)     !Semi-major axis.
Rorb_star = ((Orbital_period**2.0D0*G*((Mp**(3.0D0))/(Ms + Mp)**(2.0D0)))/(4.0d0*pi**2.0D0))**(1.0d0/3.0d0)     !Semi-major of stellar axis.

PRINT *, "Semi-major axis: ", Rorb/AU
PRINT *, "Stellar Semi-major axis: ", Rorb_star/AU

!Run model again but with the best parameters.
num_compare = 0
Aplan = pi*Rp**2.0D0                         !Surface area of the planet.

Rs2 = Rs**2.0D0                              !Square the radius of star to speed up calculations.
Rp2 = Rp**2.0D0                              !Square the radius of planet to speed up calculations.

IF (linear_quadratic == 'q') THEN
    Io = 6.0D0 / (pi*Rs2*(6.0D0 - (2.0D0*q_1) - q_2))     !The initial light intensity equation with limb darkening (quadratic)
ELSE
    Io = 1.0D0 / (pi*Rs2*(1.0D0-(u/3.0D0)))      !The initial light intensity equation with limb darkening (linear)
                                                !(normalize Io such that total star luminosity is 1). 
END IF

Ecc = Ecc_prior

PRINT *, "Ecc: ", Ecc

Inc = Inc_prior

PRINT *, "Inc: ", Inc

IF (Ecc == 0) THEN 
    !Maximum amplitude caused by the exoplanet in a circular orbit.
    RVamplitude = SQRT((G*(Mp**(3.0D0))*(sin(Inc*(pi/180.0D0)))**(3.0D0))/(Rorb_star*sin(Inc*(pi/180.0D0))*(Ms + Mp)**(2.0D0)))
ELSE
    !Maximum amplitude caused by the exoplanet in an eccentric orbit.
    RVamplitude = SQRT((G*(Mp**(3.0D0))*(sin(Inc*(pi/180.0D0)))**(3.0D0))/((1.0D0 - Ecc**2.0D0)*Rorb_star*sin(Inc*(pi/180.0D0))*(Ms + Mp)**(2.0D0)))
END IF

PRINT *, "RVamplitude: ", RVamplitude

vturb = SQRT(beta_rot**2.0D0 + vmacro**2.0D0)

omega_arg_periastron = omega_arg_periastron_prior

True_anomaly_start = pi - (omega_arg_periastron*(pi/180.0D0))
PRINT *, "True_anomaly_start", True_anomaly_start*(180.0D0/pi)

IF (True_anomaly_start >= pi) THEN
    ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) + 2.0D0*pi 
ELSE IF (True_anomaly_start <= -pi) THEN
    ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) - 2.0D0*pi
ELSE
    ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc)))
END IF

PRINT *, "ecc_anomaly_start: ", ecc_anomaly_start*(180.0D0/pi)

!Calculate the amount of time (in seconds) after the passage of periastron has occurred when the planet is at mid transit.
!This is determined from the transit mid time and the argument of periastron.
IF (omega_arg_periastron > 180.0D0) THEN
    Mean_anomaly_transit = (2.0D0*pi) - (omega_arg_periastron*(pi/180.0D0)) + pi
ELSE
    Mean_anomaly_transit = pi - (omega_arg_periastron*(pi/180.0D0))
END IF

PRINT *, "Mean_anomaly_transit: ", Mean_anomaly_transit*(180.0D0/pi)

JD_time_mid_transit = JD_time_mid_transit_prior

JD_time_peri = (JD_time_mid_transit*day_sec) - ((Mean_anomaly_transit*Orbital_period)/(2.0D0*pi))

PRINT *, "JD_time_peri: ", JD_time_peri

time_peri_passage = (JD_time_mid_transit*day_sec) - JD_time_peri

PRINT *, "time_peri_passage: ", time_peri_passage

IF (Ecc == 0) THEN
    Time_start = ((ecc_anomaly_start*Orbital_period)/(2.0D0*pi)) + time_peri_passage 
ELSE
    Time_start = (((ecc_anomaly_start - (Ecc*sin(ecc_anomaly_start)))*Orbital_period)/(2.0D0*pi)) + time_peri_passage
END IF 

PRINT *, "Time_start: ", Time_start

!      PRINT *, "ecc_anomaly_start", ecc_anomaly_start*(180.0D0/pi)
!      PRINT *, "Time_start", Time_start

True_anomaly_transit = ((3.0D0*pi)/2.0D0) - (omega_arg_periastron*(pi/180.0D0))
!      PRINT *, "True_anomaly_transit", True_anomaly_transit*(180.0D0/pi)

PRINT *, "True_anomaly_transit: ", True_anomaly_transit

IF (True_anomaly_transit >= pi) THEN
    ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) + 2.0D0*pi 
ELSE IF (True_anomaly_transit <= -pi) THEN
    ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) - 2.0D0*pi
ELSE
    ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc)))
END IF

PRINT *, "ecc_anomaly_transit: ", ecc_anomaly_transit

IF (Ecc == 0) THEN
    Time_transit = ((ecc_anomaly_transit*Orbital_period)/(2.0D0*pi)) + time_peri_passage
ELSE
    Time_transit = (((ecc_anomaly_transit - (Ecc*sin(ecc_anomaly_transit)))*Orbital_period)/(2.0D0*pi)) + time_peri_passage
END IF

PRINT *, "Time_transit: ", Time_transit

!Adjust data time relative to mid transit time based on the new orbital period and JD time mid transit.
DO l = 1, datafilelength
    Data_adjusted(l,1) = Data(l,1)
    Data_adjusted(l,2) = Data(l,2)
    Data_adjusted(l,3) = Data(l,3)
END DO

!Verify that data array is ascending order as a function of time
CALL Selection_sort(Data_adjusted(:,1),index_array)
  
DO l = 1, datafilelength
    sorted_data_array(l,1) = Data_adjusted(index_array(l), 1)
    sorted_data_array(l,2) = Data_adjusted(index_array(l), 2)
    sorted_data_array(l,3) = Data_adjusted(index_array(l), 3)
END DO

!Create a new time array for the data which has been created for the JD time of mid transit and ecc anomalies.
DO l = 1, datafilelength
   time_data_array(l) = (sorted_data_array(l,1)*day_sec) + Time_transit
END DO

RV_zero_offset = RV_zero_offset_prior
DO l = 1, datafilelength
    RV_offset_data_array(l,1) = sorted_data_array(l,1)
    RV_offset_data_array(l,2) = sorted_data_array(l,2)
    RV_offset_data_array(l,3) = sorted_data_array(l,3)
END DO

total_interval = datafilelength

Allocate(RV_theory(total_interval, 2), data_points_compare(total_interval))

DO k = 1, total_interval
   data_points_compare(k) = 'N'
END DO

DEALLOCATE(Data, Data_my_rv, Data_adjusted, sorted_data_array, Data_my_rv_offset, Data_other_rv)




!*******************************************************************************************************************************************       
DO Time_loop = 1, total_interval
    
   !time = Time_start + ((Time_loop - 1.0D0)*Time_plot_interval)
   !The time will come from the first data point which must be normalized as a fraction of a day before mid transit.
   time = time_data_array(Time_loop)

   !Set transit flag to N when outside of transit event.
   transit_flag = 'N'
   !Set occultation flag to N when outside of secondary transit event.
   occultation_flag = 'N'  

   time_compare_flag = 'N' 
   
   IF (Ecc == 0) THEN
      ecc_anomaly = ((2.0D0*pi)/Orbital_period)*(time - time_peri_passage)
   ELSE
      !Calculate the eccentric anomaly using the functions found near the end of the program.
      sum_ecc_anomaly = 0
      ecc_anomaly_before = 1000000
      DO order = 1, 20
         !Calculate the value of the Bessel function which is used to find the eccentric anomaly.
         !Bessel_value = BESSJ(order, order*Ecc)
         !Fortran 2008 Bessel function.
         Bessel_value = bessel_jn(order, order*Ecc)
         !PRINT *, "Bessel_value = ", Bessel_value
         IF (order > 1) THEN
            ecc_anomaly_before = sum_ecc_anomaly
         END IF
         sum_ecc_anomaly = sum_ecc_anomaly + ((2.0D0/order)*Bessel_value*sin(order*(((2.0D0*pi) &
                           /Orbital_period)*(time - time_peri_passage))))
         !PRINT *, "sum_ecc_anomaly = ", sum_ecc_anomaly
         ecc_anomaly_after = sum_ecc_anomaly
         IF ((order > 1) .AND. (ABS(ecc_anomaly_after - ecc_anomaly_before) <= Bessel_function_exit)) EXIT 
      END DO
      ecc_anomaly = (((2.0D0*pi)/Orbital_period)*(time - time_peri_passage)) + sum_ecc_anomaly
   END IF

   !Now use the ecc_anomaly to determine the True_anomaly and Phase_angle for a given data point.
   True_anomaly = 2.0D0*(atan(tan(ecc_anomaly/2.0D0)*(sqrt((1.0D0 + Ecc)/(1.0D0 - Ecc)))))
   
   IF (Ecc == 0) THEN
      !The time of the similation for a specific true_anomaly with the mid transit time equal to 0.
      Time_check = ((ecc_anomaly*Orbital_period)/(2.0D0*pi)) + time_peri_passage
      !The distance between the center of the planet to the center of the star in a circular orbit.
      !Planet_star_distance = Rorb
      Planet_star_distance = Rorb + Rorb_star
   ELSE
      !The time of the similation for a specific true_anomaly.
      Time_check = (((ecc_anomaly - (Ecc*sin(ecc_anomaly)))*Orbital_period)/(2.0D0*pi)) + time_peri_passage         
      !The distance between the center of the planet to the center of the star in an eccentric orbit.
      !Planet_star_distance = (Rorb*(1.0D0 - Ecc**2.0D0))/(1.0D0 + (Ecc*cos(True_anomaly)))
      Planet_star_distance = ((Rorb + Rorb_star)*(1.0D0 - Ecc**2.0D0))/(1.0D0 + (Ecc*cos(True_anomaly)))
   END IF
   !Time in reference to mid transit time
   Time_ref = Time - Time_transit

   !The position of the planet on the x-axis.
   Xpos = Planet_star_distance*((-sin(True_anomaly)*sin((omega_arg_periastron)*(pi/180.0D0))) &
          + (cos(True_anomaly)*cos((omega_arg_periastron)*(pi/180.0D0))))
   !The position of the planet on the y-axis.
   Ypos = Planet_star_distance*((cos(True_anomaly + pi)*cos((Inc)*(pi/180.0D0))*sin((omega_arg_periastron)*(pi/180.0D0))) &
          + (sin(True_anomaly + pi)*cos((Inc)*(pi/180.0D0))*cos((omega_arg_periastron)*(pi/180.0D0))))
   !The position of the planet on the z-axis.
   Zpos = Planet_star_distance*((-cos((omega_arg_periastron)*(pi/180.0D0))*sin((Inc)*(pi/180.0D0))*sin(True_anomaly)) &
          - (cos(True_anomaly)*sin((Inc)*(pi/180.0D0))*sin((omega_arg_periastron)*(pi/180.0D0))))
   Dist2 = Xpos**2.0D0 + Ypos**2.0D0                   !Square of the planet-star apparent seperation.
   Distance_center = sqrt(Dist2)               !Apparent seperation between the planet and the star.
   Lblocked = 0.0D0
   Lblocked2 = 0.0D0                           !A variable for the Anomalous velocity equation.
   v_rm = 0.0D0                                !Anomalous velocity of each pixel set to zero.
   Total_RM = 0.0D0

   IF (Xpos <= 0 .AND. Zpos >= 0) THEN
      !Planet is currently in quadrant three so add pi.
      IF (Zpos == 0) THEN
         phase_angle = pi/2.0D0
         Phase_angle_observed = pi/2.0D0
      ELSE
         phase_angle = atan(Xpos/Zpos) + pi
         !Taking into account orbital inclination.
         Phase_angle_observed = atan(-(sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos) + pi
      END IF
   ELSE IF (Xpos >= 0 .AND. Zpos >= 0) THEN
      !Planet is currently in quadrant four so add pi.
      IF (Zpos == 0) THEN
         phase_angle = pi/2.0D0 + pi
         Phase_angle_observed = pi/2.0D0 + pi
      ELSE
         phase_angle = atan(Xpos/Zpos) + pi
         !Taking into account orbital inclination.
         Phase_angle_observed = atan((sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos) + pi
      END IF
   ELSE IF (Xpos >= 0 .AND. Zpos <= 0) THEN
      !Planet is currently in quadrant one so add 2pi.
      IF (Zpos == 0) THEN
         phase_angle = pi/2.0D0
         Phase_angle_observed = pi/2.0D0
      ELSE
         phase_angle = 2.0D0*pi + atan(Xpos/Zpos)
         !Taking into account orbital inclination.
         Phase_angle_observed = 2.0D0*pi + atan((sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos)
      END IF
   ELSE IF (Xpos <= 0 .AND. Zpos <= 0) THEN
      !Planet is currently in quadrant two so add 2pi.
      IF (Zpos == 0) THEN
         phase_angle = pi/2.0D0 + pi
         Phase_angle_observed = pi/2.0D0 + pi
      ELSE
         phase_angle = atan(Xpos/Zpos)
         !Taking into account orbital inclination.
         Phase_angle_observed = atan(-(sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos)
      END IF
   END IF

   True_phase = acos(sin(True_anomaly + (omega_arg_periastron*(pi/180.0D0)))*sin(Inc*(pi/180.0D0)))
   Phase_orbit_n = acos(sin(True_anomaly + (omega_arg_periastron*(pi/180.0D0))))

   !If the planet is neither infront of the star (transit) or behind the star (occultation), then calculate the flux being reflected 
   !off the surface of the exoplanet based on its bond albedo, radius, phase angle, etc.
   IF (Distance_center > (Rs + Rp)) THEN
      Aplan = pi*Rp2                                               !Surface area of the planet.
      Total_L = 1.0D0 + ((Aplan/(4.0D0*pi*Planet_star_distance**2.0D0))*Albedo*(0.5D0*(1.0D0+cos(True_phase))))
   END IF
   
   !Hopefully this will ensure that the same data points are compared despite the length of the transit changing due to verying the parameters.
   IF ((Time_ref >= -length_time_compare) .AND. (Time_ref <= length_time_compare)) THEN
      time_compare_flag = 'Y'
   END IF 
   
   IF ((Distance_center <= (Rs + Rp)) .AND. (Zpos > 0)) THEN            
      !If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is positive then the transit begins.
      transit_flag = 'Y'                  !Planet is transiting. Will be used to determine when to do the chi squared analysis.

      !PRINT *, 'Time_ref during transit: ', Time_ref

      IF ((Rp2/Rs2) >= 0.030) THEN
         Radius_planet_array = M_pixels / 2.0D0         !Radius of planet in the pixel array. 
         Center_of_planet_x = M_pixels / 2.0D0          !The center of the planet on the X-axis.
         Center_of_planet_y = M_pixels / 2.0D0          !The center of the planet on the Y-axis.
         Pixel = Rp/Radius_planet_array               !The number of meters per pixel.
         Area_pixel = Pixel**2.0D0                        !The area of each pixel.
         Io_Pixel = Io * Area_pixel                   !Variable to speed up calculations (also physically represents the
                                                      !luminosity of the brightest pixel).    

         DO i = 1, M_pixels
            X_pixel = (pixel * i) + (Xpos - Rp)           !Calculates the location of the pixel on the x-axis.
            X_pixel2 = X_pixel**2.0D0
            XpXp = X_pixel * Xpos                         !temporary var for speed calculation
            DO j = 1, M_pixels
               Y_pixel = (pixel * j) + abs(Ypos) - Rp        !Calculates the location of the pixel on the y-axis.
               !Calculates the location of the pixel along the x-axis of the rotation axis of the star.
               X_pixel_prime = (X_pixel*cos((best_spin_orbit_angle*pi)/180.0D0)) + &
                               (Y_pixel*sin((best_spin_orbit_angle*pi)/180.0D0))     
               Dist_center_pixel = X_pixel2 + Y_pixel**2.0D0         !squared distance of pixel from star
               !Calculates the values of the pixels according to how far away they are from the center of the star and the plane.
               !squared distance of pixel from planet using a limb darkening equation.
               Dist_planet_pixel = Dist_center_pixel  - (2.0D0*(XpXp + (Y_pixel*Ypos))) + Dist2
               Sub_planet_velocity = best_vsini*(X_pixel_prime/Rs)
               IF ((Dist_center_pixel <= Rs2) .AND. (Dist_planet_pixel <= Rp2)) THEN                
                  IF (linear_quadratic == 'q') THEN
                     Lblocked2 = Io_Pixel*(1.0D0-q_1*(1.0D0-sqrt(abs(1.0D0-(Dist_center_pixel/Rs2)))) - &
                                 q_2*(1.0D0 - sqrt(abs(1.0D0-(Dist_center_pixel/Rs2))))**2.0D0)          !Quadratic limb darkening equation.
                  ELSE
                     Lblocked2 = Io_Pixel*(1.0D0-u*(1.0D0-sqrt(abs(1.0D0-(Dist_center_pixel/Rs2)))))     !First order limb darkening equation.          
                  END IF
                  Lblocked = Lblocked + Lblocked2                                         !First order limb darkening equation.
                  v_rm = v_rm - ((Lblocked2*Sub_planet_velocity)*((((2.0D0*vturb**2.0D0)+(2.0D0*best_vsini**2.0D0))/((2.0D0*vturb**2.0D0) &
                         + best_vsini**2.0D0))**(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity**2.0D0)/((2.0D0*vturb**2.0D0)+(best_vsini**2.0D0)))))   
                  !Anomalous velocity of each pixel.
               END IF
            END DO
         END DO
         Total_L = 1.0D0 - Lblocked
         Total_RM = 0.0D0 + v_rm                       !Total anomalous velocity for all the pixels.
         
      ELSE IF ((Rp2/Rs2) <= 0.030) THEN
  
         !Calculates the location of the center of the planet along the x-axis of the rotation axis of the star.
         X_prime = (Xpos*cos((best_spin_orbit_angle*pi)/180.0D0)) &
                   + (Ypos*sin((best_spin_orbit_angle*pi)/180.0D0))     
         set_distance_center = Distance_center         !The limb darkening equation will use this distance as long as the center of the 
                                                       !planet is inside the radius of the star.
         Sub_planet_velocity = best_vsini*(X_prime/Rs)      !Calculate the subplanetary velocity (the stellar velocity blocked by the 
                                !planetary disc) from vsini times the distance from the center of the planet to 
                                !the stellar rotation axis divided by the radius of the star.
         Io_planet = 0.0D0
            
         !Start the planet on the x-axis so that the planet is just far enough away not to touch the disk of the star. 
         IF ((Distance_center <= (Rs + Rp)) .AND. (Distance_center >= (Rs - Rp))) THEN    
                 
            dist_cent1_int = (Distance_center**2.0D0 + Rs2 - Rp2)/(2.0D0*Distance_center)
            !Location on the x-axis for the first intersection point.
            X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) &
                      + ((Ypos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
            !Location on the y-axis for the first intersection point.
            Y_int_1 = ((Ypos*dist_cent1_int)/Distance_center) &
                      - ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))            
            !Location on the x-axis for the second intersection point.
            X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) &
                      - ((Ypos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))             
            !Location on the y-axis for the second intersection point.
            Y_int_2 = ((Ypos*dist_cent1_int)/Distance_center) &
                      + ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
            !The limb darkening equation will use this distance if any part of the disc of the planet is outside the radius of the 
            !star.
            !This is the distance between the center of the star to the center of the area inside the star blocked by the planet.           
            set_distance_center = ((abs(Distance_center) - Rp) + Rs)/2.0D0
            !Next the program calculates the length between the two intersection point.
            Length = sqrt((X_int_1 - X_int_2)**2.0D0 + (Y_int_1 - Y_int_2)**2.0D0)
            !Calculate the angle between the Y position of the center of planet and the x position of the center of planet.
            !This is used to determine the center of the area inside the star blocked by the planet in the stellar rotational
            !axis coordinate system.
            Theta_inside = atan(Ypos/Xpos)
            !Calculates the distance from the center of the star (along x axis) to the center of the area inside the star blocked by 
            !the planet.
            !the pi factor added in this equation guarantees that Xpos_in has the correct sign.
            Xpos_in = set_distance_center*cos(Theta_inside + pi)
            Xpos_inside = Xpos_in

            !This makes sure Xpos_inside has the correct sign.
            IF (Xpos >= 0) THEN
               Xpos_inside = -Xpos_in
            END IF

            !Calculates the distance from the center of the star (along y axis) to the center of the area inside the star blocked by 
            !the planet.
            Ypos_inside = abs(set_distance_center*sin(Theta_inside))
            !Changes the x-coordinate to the stellar rotation axis by an angle formed between the orbital plane of the planet and the 
            !stellar roatation plane of the star.
            x_prime_distance = (Xpos_inside*cos((best_spin_orbit_angle*pi)/180.0D0)) &
                               + (Ypos_inside*sin((best_spin_orbit_angle*pi)/180.0D0))
            Sub_planet_velocity = best_vsini*(x_prime_distance/Rs)    !Calculate the subplanetary velocity (the stellar velocity blocked by 
            !the planetary disc) from vsini times the distance from the center of the planet to the stellar rotation axis divided by 
            !the radius of the star.
            Beta1 = 2.0D0*asin((0.5D0*Length)/Rp)         !Angle used to calculate the partial area of the planet in of the star.
            Alpha1 = 2.0D0*asin((0.5D0*Length)/Rs)        !Angle used to calculate the partial area of the planet in of the star.
                    
            IF (Distance_center >= sqrt(Rs2 - Rp2)) THEN   
               !The surface area of the planet when the center is outside the disk of the star.        
               Aplan = (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1))) + (0.5D0 * Rp2 * (Beta1 - sin(Beta1)))
               !Normalized area of the planet with the given limb darkening function
               Io_planet = Io * Aplan
            END IF
     
            IF (Distance_center < sqrt(Rs2 - Rp2)) THEN
               Aplan = ((pi*Rp2 + (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1)))) - (0.5D0 * Rp2 * (Beta1 - sin(Beta1))))
               !The surface area of the planet when the center is inside the disk of the star.
               Io_planet = Io * Aplan
            END IF
         END IF

         IF (Distance_center < (Rs - Rp)) THEN
            Aplan = (pi*Rp2)                        !The surface area of the planet.
            Io_planet = Io * Aplan
         END IF
         !The ratio of the area of the planet blocking the star to the area of the star utilizing the first order limb darkening 
         !equation.
         IF (linear_quadratic == 'q') THEN
            Lblocked = Io_planet*(1.0D0-q_1*(1.0D0-sqrt(abs(1.0D0-(set_distance_center**2.0D0/Rs2)))) - &
                        q_2*(1.0D0 - sqrt(abs(1.0D0-(set_distance_center**2.0D0/Rs2))))**2.0D0)          !Quadratic limb darkening equation.
         ELSE
            Lblocked = (Io_planet)*(1.0D0-(u*(1.0D0-sqrt(abs(1.0D0-(set_distance_center**2.0D0/Rs2))))))  
         END IF                                 
         v_rm = - ((Lblocked*Sub_planet_velocity)*((((2.0D0*vturb**2.0D0)+(2.0D0*best_vsini**2.0D0))/((2.0D0*vturb**2.0D0) &
                + best_vsini**2.0D0))**(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity**2.0D0)/((2.0D0*vturb**2.0D0)+(best_vsini**2.0D0)))))
         Total_L = 1.0D0 - Lblocked                                             !Total amount of light blocked by the planet.
         Total_RM = 0.0D0 + v_rm                                                  !Total anomalous velocity.
      END IF

   END IF 

   !If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is negative then the secondary transit 
   !(occulation) begins.
   IF ((Distance_center <= (Rs + Rp)) .AND. (Zpos < 0)) THEN     

      occultation_flag = 'Y'                  !Planet is occulting.       

      IF (Distance_center >= (Rs - Rp)) THEN           !Partial secondary transit.             
         dist_cent1_int = (Distance_center**2.0D0 + Rs2 - Rp2)/(2.0D0*Distance_center)
         Ypos2 = abs(Ypos)   !Set Ypos to a positive value to avoid discontinuity.  
         !Location on the x-axis for the first intersection point.
         X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) + ((Ypos2/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
         !Location on the y-axis for the first intersection point.
         Y_int_1 = ((Ypos2*dist_cent1_int)/Distance_center) - ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))            
         !Location on the x-axis for the second intersection point.
         X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) - ((Ypos2/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))       
         !Location on the y-axis for the second intersection point.
         Y_int_2 = ((Ypos2*dist_cent1_int)/Distance_center) + ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
         !The limb darkening equation will use this distance if any part of the disc of the planet is outside the radius of the 
         !star. This is the distance between the center of the star to the center of the area inside the planet blocked by the star.     

         set_distance_center = ((abs(Distance_center) - Rp) + Rs)/2.0D0
         !Next the program calculates the length between the two intersection point.
         Length = sqrt((X_int_1 - X_int_2)**2.0D0 + (Y_int_1 - Y_int_2)**2.0D0)
         !Calculate the angle between the Y position of the center of planet and the x position of the center of planet.
         Theta_inside = atan(Ypos2/Xpos)
         !Calculates the distance from the center of the star (along x axis) to the center of the area inside the star blocked by 
         !the planet.
         !the pi factor added in this equation guarneetes that Xpos_in has the correct sign.
         Xpos_in = set_distance_center*cos(Theta_inside + pi)
         Xpos_inside = Xpos_in
   
         !This makes sure Xpos_inside has the correct sign.
         IF (Xpos >= 0) THEN
            Xpos_inside = -Xpos_in
         END IF

         !Calculates the distance from the center of the star (along y axis) to the center of the area inside the star blocked by 
         !the planet.
         Ypos_inside = set_distance_center*sin(Theta_inside)
         Beta1 = 2.0D0*asin((0.5D0*Length)/Rp)         !Angle used to calculate the partial area of the planet in the star.
         Alpha1 = 2.0D0*asin((0.5D0*Length)/Rs)        !Angle used to calculate the partial area of the planet in the star.
                    
         IF (Distance_center >= sqrt(Rs2 - Rp2)) THEN  
            !The surface area of the planet that is not behind the star when the center of the planets disk is visible.        
            Aplan = (pi*Rp2) - ((0.5D0 * Rp2 * (Beta1 - sin(Beta1))) + (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1))))
         END IF
     
         IF (Distance_center < sqrt(Rs2 - Rp2)) THEN
            !The surface area of the planet that is not behind the star when the center of the planets disk is behind the star.
            Aplan = (0.5D0 * Rp2 * (Beta1 - sin(Beta1)))
         END IF
      END IF
   
      IF (Distance_center < (Rs - Rp)) THEN
         Aplan = 0.0D0                      !No flux coming from planet since it's behind the star.
      END IF
  
      !Total amount of light blocked by the planet.
      Total_L = 1.0D0 + ((Aplan/(4.0D0*pi*Planet_star_distance**2.0D0))*Albedo*(0.5D0*(1.0D0 + cos(True_phase))))
   END IF 

   IF (Ecc == 0) THEN
      !The radial velocity of the star which includes adding the RM effect for a circular orbit.
      RV = RVamplitude*(cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0D0)))) + Total_RM
   ELSE
      !The radial velocity of the star which includes adding the RM effect for an eccentric orbit.
      RV = RVamplitude*(cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0D0))) &
           + (Ecc*cos((omega_arg_periastron)*(pi/180.0D0) + pi))) + Total_RM
   END IF

   !If transit flag is Y, then indicate the data point to compare.
   IF ((transit_flag == 'Y')  .OR. (time_compare_flag == 'Y')) THEN
      num_compare = num_compare + 1
      data_points_compare(Time_loop) = 'Y'
   END IF
   
   !Put data into the Lightcurve array from the total luminosity.
   RV_theory(Time_loop,1) = Time_ref
   RV_theory(Time_loop,2) = RV

END DO
                           
!*******************************************************************************************************************************

PRINT *, "num_compare: ", num_compare

Number_points1 = 1
model_data_difference = 0

Allocate(Residuals_array(num_compare,3))
Allocate(Data_in_fit(num_compare,3))
Allocate(RV_Theory_fit(num_compare,2))
Allocate(Residuals_all_array(total_interval,3))
Allocate(RV_Theory_all(total_interval,2))

num_compare_out = num_compare

DO s = 1, total_interval
   IF (data_points_compare(s) == 'Y') THEN
      Residuals_array(Number_points1,1) = RV_offset_data_array(s,1)
      Residuals_array(Number_points1,2) = RV_offset_data_array(s,2) - RV_Theory(s,2)
      Residuals_array(Number_points1,3) = RV_offset_data_array(s,3)
      Data_in_fit(Number_points1,1) = RV_offset_data_array(s,1)
      Data_in_fit(Number_points1,2) = RV_offset_data_array(s,2)
      Data_in_fit(Number_points1,3) = RV_offset_data_array(s,3)
      RV_Theory_fit(Number_points1,1) = RV_Theory(s,1)/day_sec
      RV_Theory_fit(Number_points1,2) = RV_Theory(s,2)
      model_data_difference = model_data_difference + ABS(RV_offset_data_array(s,2) - RV_Theory(s,2)) 
      Number_points1 = Number_points1 + 1              !Counter to determine the number of elements used to calculate sigma squared.
   END IF
END DO

PRINT *, "model_data_difference: ", model_data_difference
PRINT *, "Number_points1: ", Number_points1 - 1

DO s = 1, total_interval
   Residuals_all_array(s,1) = RV_offset_data_array(s,1)
   Residuals_all_array(s,2) = RV_offset_data_array(s,2) - RV_Theory(s,2)
   Residuals_all_array(s,3) = RV_offset_data_array(s,3)
   RV_Theory_all(s,1) = RV_Theory(s,1)/day_sec
   RV_Theory_all(s,2) = RV_Theory(s,2)
END DO   
      
DEALLOCATE(RV_theory, time_data_array, RV_offset_data_array)
DEALLOCATE(data_points_compare) 
            
!Write parameter output array so that IDL can read the array in.
OPEN(unit=99, FILE=output_residual_array_filename, status='replace', action='write')
!Now write array line by line into parameter array file.
DO l = 1, num_compare_out
   WRITE (99,200) Residuals_array(l,1), Residuals_array(l,2), Residuals_array(l,3)
   200 FORMAT(F15.6, 1X, F15.6, 1X, F15.6)
   END DO 
CLOSE(99)

!Write parameter output array so that IDL can read the array in.
OPEN(unit=99, FILE=output_data_fit_array_filename, status='replace', action='write')
!Now write array line by line into parameter array file.
DO l = 1, num_compare_out
   WRITE (99,201) Data_in_fit(l,1), Data_in_fit(l,2), Data_in_fit(l,3)
   201 FORMAT(F15.6, 1X, F15.6, 1X, F15.6)
   END DO 
CLOSE(99)

!Write parameter output array so that IDL can read the array in.
OPEN(unit=99, FILE=output_RV_theory_fit_array_filename, status='replace', action='write')
!Now write array line by line into parameter array file.
DO l = 1, num_compare_out
   WRITE (99,202) RV_Theory_fit(l,1), RV_Theory_fit(l,2)
   202 FORMAT(F15.6, 1X, F15.6)
   END DO 
CLOSE(99)

!Write parameter output array so that IDL can read the array in.
OPEN(unit=99, FILE=output_residual_all_array_filename, status='replace', action='write')
!Now write array line by line into parameter array file.
DO l = 1, total_interval
   WRITE (99,203) Residuals_all_array(l,1), Residuals_all_array(l,2), Residuals_all_array(l,3)
   203 FORMAT(F15.6, 1X, F15.6, 1X, F15.6)
   END DO 
CLOSE(99)

!Write parameter output array so that IDL can read the array in.
OPEN(unit=99, FILE=output_RV_theory_all_array_filename, status='replace', action='write')
!Now write array line by line into parameter array file.
DO l = 1, total_interval
   WRITE (99,204) RV_Theory_all(l,1), RV_Theory_all(l,2)
   204 FORMAT(F15.6, 1X, F15.6)
   END DO 
CLOSE(99)         

DEALLOCATE(Residuals_array)
DEALLOCATE(Data_in_fit)
DEALLOCATE(RV_Theory_fit)
DEALLOCATE(Residuals_all_array)
DEALLOCATE(RV_Theory_all)
DEALLOCATE(index_array)




CONTAINS

!Subroutine for ranking an array in ascending order. 
SUBROUTINE Selection_sort(a,index_array)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: a(:)
INTEGER(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: index_array
DOUBLE PRECISION :: temp_array(size(a))
INTEGER(8) :: i, minIndex, temp2
DOUBLE PRECISION :: temp
Allocate(index_array(SIZE(a)))
 
DO i = 1, SIZE(a)
   index_array(i) = i
   temp_array(i) = a(i)
END DO
 
DO i = 1, SIZE(a) - 1
   minIndex = MINLOC(temp_array(i:), 1) + i - 1
   IF (temp_array(i) > temp_array(minIndex)) THEN
      temp = temp_array(i)
      temp2 = index_array(i)
      temp_array(i) = temp_array(minIndex)
      index_array(i) = index_array(minIndex)
      temp_array(minIndex) = temp
      index_array(minIndex) = temp2
   END IF
END DO
END SUBROUTINE Selection_sort

END PROGRAM RM_effect_residuals
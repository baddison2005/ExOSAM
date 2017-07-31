PROGRAM ExOSAM_RM_effect_V02

!   Purpose: To determine the best fit values for the spin-orbit alignment and rotational velocity of the host star from a 
!        given data set of the radial velocities during the RM effect as well as out-of-transit RVs and input priors. This 
!        program accepts input data from the Python wrapper RMLEAPWE (see input data), does all the theoretical model calculations,
!        and outputs the best fit model and uncertainties back to RMLEAPWE for additional analysis and plotting (see output data).
!        This code is written in modern FORTRAN to increase speed and computational efficiency.

use ran_mod
use init_seed
implicit none

!Constants and parameter declerations.
double precision, parameter :: Rss = 6.9634D8                !Radius of our Sun (in meters).
double precision, parameter :: Rj = 7.1492D7                 !Radius of Jupiter (in meters).
double precision, parameter :: RE = 6.3781D6                 !Radius of Earth (in meters if planet is given in Earth radii).
double precision, parameter :: AU = 1.4960D11                !One astronomical unit in meters.
double precision, parameter :: pi = 3.1415926535D0           !Constant value of pi.
double precision, parameter :: G = 6.67384D-11               !Gravitation constant.
double precision, parameter :: Mss = 1.9891D30               !Mass of the sun in kilograms.
double precision, parameter :: Me = 5.9722D24                !Mass of the Earth in kg.
double precision, parameter :: Mj = 1.8986D27                !Mass of Jupiter in kg.
double precision, parameter :: day_sec = 86400.0D0           !Length of a day in seconds.

!Double precision variable declerations.
double precision :: a, b, c         !variables used to read in data
double precision :: Albedo
double precision :: Alpha1
double precision :: Aplan
double precision :: Area_pixel
double precision :: Bessel_function_exit
double precision :: Bessel_value
double precision :: Beta1
double precision :: Center_of_planet_x
double precision :: Center_of_planet_y
double precision :: chi_squared_change          !Check whether this variable will remain. May change file read in.
double precision :: chi_squared_change_fit
double precision :: chi_2
double precision :: clock_rate
double precision :: data_plot_model_interval
double precision :: Dist2
double precision :: Distance_center
double precision :: Dist_center_pixel
double precision :: dist_cent1_int
double precision :: Dist_planet_pixel
double precision :: Ecc
double precision :: Ecc_prior
double precision :: Ecc_end
double precision :: Ecc_begin
double precision :: Ecc_1sigerr
double precision :: ecc_anomaly
double precision :: ecc_anomaly_after
double precision :: ecc_anomaly_before
double precision :: ecc_anomaly_start
double precision :: ecc_anomaly_transit
double precision :: impact_prior
double precision :: Inc_prior
double precision :: Inc_end
double precision :: Inc_begin
double precision :: Inc_1sigerr
double precision :: Inc
double precision :: Io
double precision :: Io_Pixel
double precision :: Io_planet
double precision :: JD_time_peri
double precision :: JD_time_mid_transit
double precision :: JD_time_mid_transit_prior
double precision :: JD_time_mid_transit_end
double precision :: JD_time_mid_transit_begin
double precision :: JD_time_mid_transit_1sigerr
double precision :: Lblocked
double precision :: Lblocked2
double precision :: Length
double precision :: length_time_compare
double precision :: Mean_anomaly_transit
double precision :: model_data_difference
double precision :: Mp
double precision :: Mpp
double precision :: Mp_prior
double precision :: Mp_end
double precision :: Mp_begin
double precision :: Mp_1sigerr
double precision :: Ms
double precision :: Ms_solar
double precision :: Number_orbits
double precision :: omega_arg_periastron
double precision :: omega_arg_periastron_prior
double precision :: omega_arg_periastron_end
double precision :: omega_arg_periastron_begin
double precision :: omega_arg_periastron_1sigerr
double precision :: Orbital_period
double precision :: Orbital_period_1
double precision :: Orbital_period_prior
double precision :: Orbital_period_end
double precision :: Orbital_period_begin
double precision :: Orbital_period_1sigerr
double precision :: percent_processed
double precision :: phase_angle
double precision :: Phase_angle_observed
double precision :: Phase_angle_start
double precision :: Phase_orbit_n
double precision :: phase_time
double precision :: Pixel
double precision :: Planet_star_distance
double precision :: q_1
double precision :: q_2
double precision :: r_Chi_2
double precision :: Radius_planet_array
double precision :: radius_ratio_prior
double precision :: Rorb
double precision :: Rorb_prior
double precision :: Rorb_star
double precision :: Rorb_star_prior
double precision :: Rp
double precision :: Rpp
double precision :: Rp_prior
double precision :: Rp_end
double precision :: Rp_begin
double precision :: Rp_1sigerr
double precision :: Rp2
double precision :: Rs
double precision :: Rs_solar
double precision :: Rs2
double precision :: RV
double precision :: RVamplitude
double precision :: RV_offset_datasets
double precision :: RV_offset_datasets_begin
double precision :: RV_offset_datasets_end
double precision :: RV_offset_datasets_interval
double precision :: RV_offset_datasets_prior
double precision :: RV_offset_datasets_1sigerr
double precision :: RV_zero_offset
double precision :: RV_zero_offset_begin
double precision :: RV_zero_offset_end
double precision :: RV_zero_offset_interval
double precision :: RV_zero_offset_prior
double precision :: RV_zero_offset_1sigerr
double precision :: min_chi_squared_all_walkers
double precision :: min_reduced_chi_squared_all_walkers
double precision :: max_likelihood_all_walkers
double precision :: best_RV_offset_datasets_chi2_all_walkers
double precision :: best_RV_zero_offset_chi2_all_walkers
double precision :: best_orbital_period_chi2_all_walkers
double precision :: best_JD_time_mid_transit_chi2_all_walkers
double precision :: best_Mp_chi2_all_walkers
double precision :: best_Rp_chi2_all_walkers
double precision :: best_Ms_solar_chi2_all_walkers
double precision :: best_Rs_solar_chi2_all_walkers
double precision :: best_Rp_Rs_ratio_chi2_all_walkers
double precision :: best_Ecc_chi2_all_walkers
double precision :: best_Inc_chi2_all_walkers
double precision :: best_omega_arg_periastron_chi2_all_walkers
double precision :: best_vsini_chi2_all_walkers
double precision :: best_spin_orbit_chi2_all_walkers




!Variables created for MCMC.
double precision :: accept_rate = 0.0D0                           !MCMC variables.
double precision :: best_RV_offset_datasets_chi2                  !MCMC variables.
double precision :: best_orbital_period_chi2                      !MCMC variables.
double precision :: best_JD_time_mid_transit_chi2                 !MCMC variables.
double precision :: best_Mp_chi2                                  !MCMC variables.
double precision :: best_Rp_chi2                                  !MCMC variables.
double precision :: best_Ecc_chi2                                 !MCMC variables.
double precision :: best_Inc_chi2                                 !MCMC variables.
double precision :: best_omega_arg_periastron_chi2                !MCMC variables.
double precision :: best_RV_zero_offset_chi2                      !MCMC variables.
double precision :: best_vsini_chi2                               !MCMC variables.
double precision :: best_spin_orbit_chi2                          !MCMC variables.
double precision :: best_Ms_solar_chi2                            !MCMC variables.
double precision :: best_Rs_solar_chi2                            !MCMC variables.
double precision :: best_Rp_Rs_ratio_chi2                         !MCMC variables.
! ********************************** change scale_factor to user input ************************************************
double precision :: scale_factor                                  !Scale factor that controls the step-size of the MCMC.
double precision :: scale_factor_in
! ********************************** change scale_factor to user input ************************************************
double precision :: vsini_variance                                !MCMC variables.
double precision :: vsini_variance_MCMC                           !MCMC variables.
double precision :: vsini_sum                                     !MCMC variables.
double precision :: vsini_mean                                    !MCMC variables.
double precision :: vsini_diff_mean_2                             !MCMC variables.
double precision :: stellar_rotation_angle_variance               !MCMC variables.
double precision :: stellar_rotation_angle_variance_MCMC          !MCMC variables.
double precision :: likelihood                                    !MCMC variables.
double precision :: culumative_likelihood                         !MCMC variables.
double precision :: culumative_all_likelihood                     !MCMC variables.
double precision :: acceptance_param_previous                     !MCMC variables.
double precision :: acceptance_param_current                      !MCMC variables.
double precision :: stellar_rotation_angle_stand_dev              !MCMC variables.
double precision :: stellar_rotation_angle_stand_dev_MCMC         !MCMC variables.
double precision :: vsini_stand_dev                               !MCMC variables.
double precision :: vsini_stand_dev_MCMC                          !MCMC variables.
double precision :: accept_prob                                   !MCMC variables.
double precision :: random_draw                                   !MCMC variables.
double precision :: stellar_rotation_angle_diff_mean_2            !MCMC variables.
double precision :: stellar_rotation_angle_mean                   !MCMC variables.
double precision :: stellar_rotation_angle_sum                    !MCMC variables.
double precision :: RV_offset_datasets_variance                   !MCMC variables.
double precision :: RV_offset_datasets_variance_MCMC              !MCMC variables.
double precision :: RV_offset_datasets_stand_dev                  !MCMC variables.
double precision :: RV_offset_datasets_stand_dev_MCMC             !MCMC variables.
double precision :: RV_offset_datasets_sum                        !MCMC variables.
double precision :: RV_offset_datasets_mean                       !MCMC variables.
double precision :: RV_offset_datasets_diff_mean_2                !MCMC variables.
double precision :: Orbital_period_variance                       !MCMC variables.
double precision :: Orbital_period_variance_MCMC                  !MCMC variables.
double precision :: Orbital_period_stand_dev                      !MCMC variables.
double precision :: Orbital_period_stand_dev_MCMC                 !MCMC variables.
double precision :: Orbital_period_sum                            !MCMC variables.
double precision :: Orbital_period_mean                           !MCMC variables.
double precision :: Orbital_period_diff_mean_2                    !MCMC variables.
double precision :: JD_time_mid_transit_variance                  !MCMC variables.
double precision :: JD_time_mid_transit_variance_MCMC             !MCMC variables.
double precision :: JD_time_mid_transit_stand_dev                 !MCMC variables.
double precision :: JD_time_mid_transit_stand_dev_MCMC            !MCMC variables.
double precision :: JD_time_mid_transit_sum                       !MCMC variables.
double precision :: JD_time_mid_transit_mean                      !MCMC variables.
double precision :: JD_time_mid_transit_diff_mean_2               !MCMC variables.
double precision :: Mp_variance                                   !MCMC variables.
double precision :: Mp_variance_MCMC                              !MCMC variables.
double precision :: Mp_stand_dev                                  !MCMC variables.
double precision :: Mp_stand_dev_MCMC                             !MCMC variables.
double precision :: Mp_sum                                        !MCMC variables.
double precision :: Mp_mean                                       !MCMC variables.
double precision :: Mp_diff_mean_2                                !MCMC variables.
double precision :: Rp_variance                                   !MCMC variables.
double precision :: Rp_variance_MCMC                              !MCMC variables.
double precision :: Rp_stand_dev                                  !MCMC variables.
double precision :: Rp_stand_dev_MCMC                             !MCMC variables.
double precision :: Rp_sum                                        !MCMC variables.
double precision :: Rp_mean                                       !MCMC variables.
double precision :: Rp_diff_mean_2                                !MCMC variables.
double precision :: Ecc_variance                                  !MCMC variables.
double precision :: Ecc_variance_MCMC                             !MCMC variables.
double precision :: Ecc_stand_dev                                 !MCMC variables.
double precision :: Ecc_stand_dev_MCMC                            !MCMC variables.
double precision :: Ecc_sum                                       !MCMC variables.
double precision :: Ecc_mean                                      !MCMC variables.
double precision :: Ecc_diff_mean_2                               !MCMC variables.
double precision :: Inc_variance                                  !MCMC variables.
double precision :: Inc_variance_MCMC                             !MCMC variables.
double precision :: Inc_stand_dev                                 !MCMC variables.
double precision :: Inc_stand_dev_MCMC                            !MCMC variables.
double precision :: Inc_sum                                       !MCMC variables.
double precision :: Inc_mean                                      !MCMC variables.
double precision :: Inc_diff_mean_2                               !MCMC variables.
double precision :: omega_arg_periastron_variance                 !MCMC variables.
double precision :: omega_arg_periastron_variance_MCMC            !MCMC variables.
double precision :: omega_arg_periastron_stand_dev                !MCMC variables.
double precision :: omega_arg_periastron_stand_dev_MCMC           !MCMC variables.
double precision :: omega_arg_periastron_sum                      !MCMC variables.
double precision :: omega_arg_periastron_mean                     !MCMC variables.
double precision :: omega_arg_periastron_diff_mean_2              !MCMC variables.
double precision :: RV_zero_offset_variance                       !MCMC variables.
double precision :: RV_zero_offset_variance_MCMC                  !MCMC variables.
double precision :: RV_zero_offset_stand_dev                      !MCMC variables.
double precision :: RV_zero_offset_stand_dev_MCMC                 !MCMC variables.
double precision :: RV_zero_offset_sum                            !MCMC variables.
double precision :: RV_zero_offset_mean                           !MCMC variables.
double precision :: RV_zero_offset_diff_mean_2                    !MCMC variables.
double precision :: Ms_solar_prior                                !MCMC variables.
double precision :: Ms_solar_end                                  !MCMC variables.
double precision :: Ms_solar_begin                                !MCMC variables.
double precision :: Ms_solar_1sigerr                              !MCMC variables.
double precision :: Rs_solar_prior                                !MCMC variables.
double precision :: Rs_solar_end                                  !MCMC variables.
double precision :: Rs_solar_begin                                !MCMC variables.
double precision :: Rs_solar_1sigerr                              !MCMC variables.
double precision :: Rp_Rs_ratio_prior                             !MCMC variables.
double precision :: Rp_Rs_ratio_end                               !MCMC variables.
double precision :: Rp_Rs_ratio_begin                             !MCMC variables.
double precision :: Rp_Rs_ratio_1sigerr                           !MCMC variables.
double precision :: Rp_Rs_ratio_variance                          !MCMC variables.
double precision :: Rp_Rs_ratio_variance_MCMC                     !MCMC variables.
double precision :: Rp_Rs_ratio_stand_dev                         !MCMC variables.
double precision :: Rp_Rs_ratio_stand_dev_MCMC                    !MCMC variables.
double precision :: Rp_Rs_ratio                                   !MCMC variables.
double precision :: Rp_Rs_ratio_sum                               !MCMC variables.
double precision :: Rp_Rs_ratio_mean                              !MCMC variables.
double precision :: Rp_Rs_ratio_diff_mean_2                       !MCMC variables.
double precision :: Rs_solar_variance                             !MCMC variables.
double precision :: Rs_solar_variance_MCMC                        !MCMC variables.
double precision :: Rs_solar_stand_dev                            !MCMC variables.
double precision :: Rs_solar_stand_dev_MCMC                       !MCMC variables.
double precision :: Rs_solar_sum                                  !MCMC variables.
double precision :: Rs_solar_mean                                 !MCMC variables.
double precision :: Rs_solar_diff_mean_2                          !MCMC variables.
double precision :: Ms_solar_variance                             !MCMC variables.
double precision :: Ms_solar_variance_MCMC                        !MCMC variables.
double precision :: Ms_solar_stand_dev                            !MCMC variables.
double precision :: Ms_solar_stand_dev_MCMC                       !MCMC variables.
double precision :: Ms_solar_sum                                  !MCMC variables.
double precision :: Ms_solar_mean                                 !MCMC variables.
double precision :: Ms_solar_diff_mean_2                          !MCMC variables.




double precision :: stellar_rotation_angle_1sigerr_temp           !MCMC variables.
double precision :: vsini_1sigerr_temp                            !MCMC variables.
double precision :: RV_offset_datasets_1sigerr_temp               !MCMC variables.
double precision :: RV_zero_offset_1sigerr_temp                   !MCMC variables.
double precision :: Orbital_period_1sigerr_temp                   !MCMC variables.
double precision :: JD_time_mid_transit_1sigerr_temp              !MCMC variables.
double precision :: Mp_1sigerr_temp                               !MCMC variables.
double precision :: Rp_1sigerr_temp                               !MCMC variables.
double precision :: Ecc_1sigerr_temp                              !MCMC variables.
double precision :: Inc_1sigerr_temp                              !MCMC variables.
double precision :: Ms_solar_1sigerr_temp                         !MCMC variables.
double precision :: Rs_solar_1sigerr_temp                         !MCMC variables.
double precision :: Rp_Rs_ratio_1sigerr_temp                      !MCMC variables.




double precision :: set_distance_center
double precision :: stellar_rotation_angle_begin
double precision :: stellar_rotation_angle_end
double precision :: stellar_rotation_angle_prior
double precision :: stellar_rotation_angle_1sigerr
double precision :: stellar_rotation_angle
double precision :: Sub_planet_velocity
double precision :: sum_ecc_anomaly
double precision :: Transit_length_prior
double precision :: Time
double precision :: Time_bef_aft
double precision :: Time_check
double precision :: Time_compare_vel
double precision :: time_peri_passage
double precision :: Time_ref
double precision :: time_remaining = 0.0D0
double precision :: Time_transit
double precision :: Theta_inside
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
double precision :: vsini
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

double precision :: RV_offset_datasets_sum3
double precision :: RV_offset_datasets_mean2
double precision :: RV_zero_offset_sum3
double precision :: RV_zero_offset_mean2
double precision :: Orbital_period_sum3
double precision :: Orbital_period_mean2
double precision :: JD_time_mid_transit_sum3
double precision :: JD_time_mid_transit_mean2
double precision :: Mp_sum3
double precision :: Mp_mean2
double precision :: Rp_sum3
double precision :: Rp_mean2
double precision :: Ms_solar_sum3
double precision :: Ms_solar_mean2
double precision :: Rs_solar_sum3
double precision :: Rs_solar_mean2
double precision :: Rp_Rs_ratio_sum3
double precision :: Rp_Rs_ratio_mean2
double precision :: Ecc_sum3
double precision :: Ecc_mean2
double precision :: Inc_sum3
double precision :: Inc_mean2
double precision :: omega_arg_periastron_sum3
double precision :: omega_arg_periastron_mean2
double precision :: vsini_sum3
double precision :: vsini_mean2
double precision :: stellar_rotation_angle_sum3
double precision :: stellar_rotation_angle_mean2




!Integer variable declerations.
integer(8) :: mcmc_accepted_iteration_size      
integer(8) :: number_mcmc_walkers                 
integer(8) :: num_compare = 0
integer(8) :: Number_iterations_left_total
integer(8) :: Number_points1 = 1                           !The number of points that are used from the predicted light curve.
integer(8) :: num_my_rv
integer(8) :: num_other_rv
integer(8) :: Number_fit
integer(8) :: order
integer(8) :: percent_count=1
integer(4) :: exitstat
integer(4) :: cmdstat
integer(8) :: Time_loop
integer(8) :: M_pixels
integer(8) :: i, l, j, s, e, k, zz            !Loop variables
integer(8) :: percent_processed_loop
integer(8) :: datafilelength
integer(8) :: string_bar_num = 1 
integer(8) :: total_interval
integer(8) :: reject_counter = 0                                  !MCMC variables.
integer(8) :: MCMC_loop = 1                                           !MCMC variables.
integer(8) :: all_loop = 1                                        !MCMC variables.
integer(8) :: total_loop = 1
integer(8) :: MCMC_print_updates = 0                              !MCMC variables.
integer(8) :: MCMC_updates = 0                                    !MCMC variables.
integer(8) :: count_tested_proposals = 0                          !MCMC variables.
integer(8) :: walker
integer(8) :: Ecc_draw_loop
integer(8) :: vsini_draw_loop
integer(8) :: RV_offset_datasets_draw_loop
integer(8) :: Orbital_period_draw_loop
integer(8) :: JD_time_mid_transit_draw_loop
integer(8) :: Mp_draw_loop
integer(8) :: Rp_draw_loop
integer(8) :: Rp_Rs_ratio_draw_loop
integer(8) :: Rs_solar_draw_loop
integer(8) :: Ms_solar_draw_loop
integer(8) :: Inc_draw_loop
integer(8) :: omega_arg_periastron_draw_loop
integer(8) :: RV_zero_offset_draw_loop
integer(8) :: stellar_rotation_angle_draw_loop
integer(8) :: cr, cm, system_time_1, system_time_it1, system_time_it2
integer(8) :: system_time_4

!character variable declerations.
character(LEN=100) :: cmd1
character(len=1), Allocatable :: cmdmsg(:)
character(LEN=1) :: add_RV_zero_offset
character(LEN=1) :: Phase_shift                          !This parameter needs to be removed in the input file.
character(LEN=1) :: subtract_JD_mid_time                 !This parameter needs to be removed in the input file.

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

character(LEN=25) :: string_bar
character(LEN=1) :: Jupiter_Earth_units
character(LEN=1) :: blank = ' '
character(LEN=1) :: transit_flag
character(LEN=1) :: time_compare_flag = 'N'
character(LEN=1) :: other_RV_files
character(LEN=1) :: array_resorted = 'N'
character(LEN=60) :: input_data_filename
character(LEN=60) :: current_directory
character(LEN=150) :: input_RM_filename
character(LEN=150) :: input_my_data_filename
character(LEN=150) :: input_other_data_filename
character(LEN=20) :: Planet_name
character(LEN=150) :: output_temp_data
character(LEN=150) :: directory_location
character(LEN=150) :: my_data_plotting_symbols
character(LEN=150) :: other_data_plotting_symbols
character(LEN=1) :: linear_quadratic
character(LEN=1) :: Reject_flag='F'             !MCMC parameter.
character(LEN=1) :: impose_prior_vsini          !MCMC parameter.
character(LEN=1) :: use_out_transit_rv_for_fit



!Character strings created for MCMC.
character(LEN=150) :: output_chi_squared_array_all_filename
character(LEN=150) :: output_chi_squared_array_MCMC_filename
character(LEN=150) :: output_reduced_chi_array_all_filename
character(LEN=150) :: output_reduced_chi_array_MCMC_filename
character(LEN=150) :: output_likelihood_array_all_filename
character(LEN=150) :: output_likelihood_array_MCMC_filename
character(LEN=150) :: output_vsini_all_array_filename
character(LEN=150) :: output_vsini_MCMC_array_filename
character(LEN=150) :: output_spin_orbit_all_array_filename
character(LEN=150) :: output_spin_orbit_MCMC_array_filename
character(LEN=150) :: output_spin_orbit_norm_array_filename
character(LEN=150) :: output_parameter_array_all_filename
character(LEN=150) :: output_parameter_array_MCMC_filename
character(LEN=150) :: output_MCMC_number_filename
character(LEN=150) :: output_all_number_filename
character(LEN=150) :: output_parameter_cutoff_filename
character(LEN=150) :: output_best_parameters_mcmc_filename
character(LEN=150) :: output_unweight_variance_vsini_lambda_filename
character(LEN=150) :: output_weight_variance_vsini_lambda_filename
character(LEN=150) :: output_delta_chi_uncert_vsini_lambda_filename
character(LEN=150) :: output_sigma_uncert_vsini_lambda_filename
character(LEN=150) :: output_mask_array_filename
character(LEN=150) :: output_RV_offset_datasets_MCMC_array_filename
character(LEN=150) :: output_RV_offset_datasets_all_array_filename
character(LEN=150) :: output_RV_zero_offset_MCMC_array_filename
character(LEN=150) :: output_RV_zero_offset_all_array_filename
character(LEN=150) :: output_Orbital_period_MCMC_array_filename
character(LEN=150) :: output_Orbital_period_all_array_filename
character(LEN=150) :: output_JD_time_mid_transit_MCMC_array_filename
character(LEN=150) :: output_JD_time_mid_transit_all_array_filename
character(LEN=150) :: output_Mp_MCMC_array_filename
character(LEN=150) :: output_Mp_all_array_filename
character(LEN=150) :: output_Rp_MCMC_array_filename
character(LEN=150) :: output_Rp_all_array_filename
character(LEN=150) :: output_Ms_solar_MCMC_array_filename
character(LEN=150) :: output_Ms_solar_all_array_filename
character(LEN=150) :: output_Rs_solar_MCMC_array_filename
character(LEN=150) :: output_Rs_solar_all_array_filename
character(LEN=150) :: output_Rp_Rs_ratio_MCMC_array_filename
character(LEN=150) :: output_Rp_Rs_ratio_all_array_filename
character(LEN=150) :: output_Ecc_MCMC_array_filename
character(LEN=150) :: output_Ecc_all_array_filename
character(LEN=150) :: output_Inc_MCMC_array_filename
character(LEN=150) :: output_Inc_all_array_filename
character(LEN=150) :: output_omega_arg_periastron_MCMC_array_filename
character(LEN=150) :: output_omega_arg_periastron_all_array_filename
character(LEN=150) :: output_MCMC_properties_filename
character(LEN=150) :: output_best_param_chi_squared_filename
character(LEN=150) :: output_best_param_mean_filename
character(LEN=150) :: output_stand_dev_filename
character(LEN=150) :: output_scale_factor_array_filename




!Data arrays.
DOUBLE PRECISION, Allocatable :: Data(:,:), Data_my_rv_offset(:,:), Data_my_rv(:,:), Data_other_rv(:,:), Data_adjusted(:,:)
DOUBLE PRECISION, Allocatable :: RV_theory(:,:), time_data_array(:), sorted_data_array(:,:), RV_offset_data_array(:,:)
INTEGER(8), ALLOCATABLE :: index_array(:)
CHARACTER(len=1), Allocatable :: data_points_compare(:)



!Arrays created for MCMC
DOUBLE PRECISION, Allocatable :: RV_offset_datasets_accept_mcmc_array(:,:), RV_offset_datasets_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: RV_offset_datasets_sum2(:)
DOUBLE PRECISION, Allocatable :: Orbital_period_accept_mcmc_array(:,:), Orbital_period_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Orbital_period_sum2(:)
DOUBLE PRECISION, Allocatable :: JD_time_mid_transit_accept_mcmc_array(:,:), JD_time_mid_transit_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: JD_time_mid_transit_sum2(:)
DOUBLE PRECISION, Allocatable :: Mp_accept_mcmc_array(:,:), Mp_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Mp_sum2(:)
DOUBLE PRECISION, Allocatable :: Rp_accept_mcmc_array(:,:), Rp_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Rp_sum2(:)
DOUBLE PRECISION, Allocatable :: Rp_Rs_ratio_accept_mcmc_array(:,:), Rp_Rs_ratio_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Rp_Rs_ratio_sum2(:)
DOUBLE PRECISION, Allocatable :: Ecc_accept_mcmc_array(:,:), Ecc_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Ecc_sum2(:)
DOUBLE PRECISION, Allocatable :: Inc_accept_mcmc_array(:,:), Inc_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Inc_sum2(:)
DOUBLE PRECISION, Allocatable :: omega_arg_periastron_accept_mcmc_array(:,:), omega_arg_periastron_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: omega_arg_periastron_sum2(:)
DOUBLE PRECISION, Allocatable :: RV_zero_offset_accept_mcmc_array(:,:), RV_zero_offset_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: RV_zero_offset_sum2(:)
DOUBLE PRECISION, Allocatable :: Rs_solar_accept_mcmc_array(:,:), Rs_solar_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Rs_solar_sum2(:)
DOUBLE PRECISION, Allocatable :: Ms_solar_accept_mcmc_array(:,:), Ms_solar_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: Ms_solar_sum2(:)
DOUBLE PRECISION, Allocatable :: vsini_accept_mcmc_array(:,:), vsini_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: vsini_sum2(:)
DOUBLE PRECISION, Allocatable :: stellar_rotation_angle_accept_mcmc_array(:,:), stellar_rotation_angle_all_mcmc_array(:,:)
DOUBLE PRECISION, Allocatable :: stellar_rotation_angle_sum2(:)
DOUBLE PRECISION, Allocatable :: vsini_variance_current_array(:,:), vsini_variance_MCMC_array(:,:)
DOUBLE PRECISION, Allocatable :: stellar_rotation_angle_variance_current_array(:,:), stellar_rotation_angle_variance_MCMC_array(:,:)
DOUBLE PRECISION, Allocatable :: stellar_rotation_angle_stand_dev_current_array(:,:)
DOUBLE PRECISION, Allocatable :: stellar_rotation_angle_stand_dev_MCMC_array(:,:)
DOUBLE PRECISION, Allocatable :: vsini_stand_dev_current_array(:,:)
DOUBLE PRECISION, Allocatable :: vsini_stand_dev_MCMC_array(:,:)
DOUBLE PRECISION, Allocatable :: Chi_squared_array_MCMC(:,:), Reduced_Chi_Squared_array_MCMC(:,:)
DOUBLE PRECISION, Allocatable :: likelihood_array_MCMC(:,:), culumative_likelihood_array_MCMC(:,:)
integer(8), Allocatable :: loc_min_chi_squared_MCMC(:), loc_min_reduced_chi_squared_MCMC(:), loc_max_likelihood_MCMC(:)
DOUBLE PRECISION, Allocatable :: walker_parameter_array(:,:), walker_parameter_array_chi2(:,:), walker_parameter_array_mean(:,:)
DOUBLE PRECISION, Allocatable :: walker_parameter_array_variance(:,:), walker_parameter_array_stand_dev(:,:)
DOUBLE PRECISION, Allocatable :: min_chi_squared_MCMC(:), min_reduced_chi_squared_MCMC(:), max_likelihood_MCMC(:)
DOUBLE PRECISION, Allocatable :: scale_factor_array(:,:)
LOGICAL, Allocatable :: mask_array(:,:)

integer(8), DIMENSION(2) :: loc_min_chi_squared_all_walkers, loc_min_reduced_chi_squared_all_walkers
integer(8), DIMENSION(2) :: loc_max_likelihood_all_walkers




DOUBLE PRECISION, Allocatable :: normal1(:,:), normal2(:,:), normal3(:,:), normal4(:,:), normal5(:,:), normal6(:,:)
DOUBLE PRECISION, Allocatable :: normal7(:,:), normal8(:,:), normal9(:,:), normal10(:,:), normal11(:,:), normal12(:,:)
DOUBLE PRECISION, Allocatable :: normal13(:,:), normal14(:,:), random_draw_array(:,:)

DOUBLE PRECISION, Allocatable :: spread_rand_RV_offset(:,:), spread_rand_Orbital_period(:,:), spread_rand_JD(:,:)
DOUBLE PRECISION, Allocatable :: spread_rand_Mp(:,:), spread_rand_Rp(:,:), spread_rand_Rp_Rs_ratio(:,:)
DOUBLE PRECISION, Allocatable :: spread_rand_Rs(:,:), spread_rand_Ms(:,:), spread_rand_Ecc(:,:)
DOUBLE PRECISION, Allocatable :: spread_rand_Inc(:,:), spread_rand_omega(:,:), spread_rand_RV_zero(:,:)
DOUBLE PRECISION, Allocatable :: spread_rand_vsini(:,:), spread_rand_spin_orbit(:,:)




!Get the locations of files.
!CALL getcwd(current_directory)
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
output_chi_squared_array_all_filename = TRIM(TRIM(output_temp_data) // 'chi_squared_array_all.txt')
output_chi_squared_array_MCMC_filename = TRIM(TRIM(output_temp_data) // 'chi_squared_array_mcmc.txt')
output_reduced_chi_array_all_filename = TRIM(TRIM(output_temp_data) // 'reduced_chi_squared_array_all.txt')
output_reduced_chi_array_MCMC_filename = TRIM(TRIM(output_temp_data) // 'reduced_chi_squared_array_mcmc.txt')
output_likelihood_array_all_filename = TRIM(TRIM(output_temp_data) // 'likelihood_array_all.txt')
output_likelihood_array_MCMC_filename = TRIM(TRIM(output_temp_data) // 'likelihood_array_mcmc.txt')
output_vsini_all_array_filename = TRIM(TRIM(output_temp_data) // 'vsini_array_all.txt')
output_vsini_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'vsini_array_mcmc.txt')
output_spin_orbit_all_array_filename = TRIM(TRIM(output_temp_data) // 'stellar_rotation_angle_array_all.txt')
output_spin_orbit_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'stellar_rotation_angle_array_mcmc.txt')
output_parameter_array_all_filename = TRIM(TRIM(output_temp_data) // 'parameter_array_all.txt')
output_spin_orbit_norm_array_filename = TRIM(TRIM(output_temp_data) // 'spin_orbit_norm_array.txt')
output_parameter_array_MCMC_filename = TRIM(TRIM(output_temp_data) // 'parameter_array_mcmc.txt')
output_MCMC_number_filename = TRIM(TRIM(output_temp_data) // 'number_mcmc.txt')
output_all_number_filename = TRIM(TRIM(output_temp_data) // 'number_all.txt')
output_parameter_cutoff_filename = TRIM(TRIM(output_temp_data) // 'parameter_cutoff_mcmc.txt')
output_best_parameters_mcmc_filename = TRIM(TRIM(output_temp_data) // 'best_parameters_mcmc.txt')
output_unweight_variance_vsini_lambda_filename = TRIM(TRIM(output_temp_data) // 'unweighted_uncertainty_vsini_lambda.txt')
output_weight_variance_vsini_lambda_filename = TRIM(TRIM(output_temp_data) // 'weighted_uncertainty_vsini_lambda.txt')
output_delta_chi_uncert_vsini_lambda_filename = TRIM(TRIM(output_temp_data) // 'delta_chi_uncertainty_vsini_lambda.txt')
output_sigma_uncert_vsini_lambda_filename = TRIM(TRIM(output_temp_data) // 'sigma_uncertainty_vsini_lambda.txt')
output_mask_array_filename = TRIM(TRIM(output_temp_data) // 'mask_array.txt')
output_RV_offset_datasets_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'RV_offset_datasets_MCMC_array.txt')
output_RV_offset_datasets_all_array_filename = TRIM(TRIM(output_temp_data) // 'RV_offset_datasets_all_array.txt')
output_RV_zero_offset_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'RV_zero_offset_MCMC_array.txt')
output_RV_zero_offset_all_array_filename = TRIM(TRIM(output_temp_data) // 'RV_zero_offset_all_array.txt')
output_Orbital_period_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Orbital_period_MCMC_array.txt')
output_Orbital_period_all_array_filename = TRIM(TRIM(output_temp_data) // 'Orbital_period_all_array.txt')
output_JD_time_mid_transit_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'JD_time_mid_transit_MCMC_array.txt')
output_JD_time_mid_transit_all_array_filename = TRIM(TRIM(output_temp_data) // 'JD_time_mid_transit_all_array.txt')
output_Mp_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Mp_MCMC_array.txt')
output_Mp_all_array_filename = TRIM(TRIM(output_temp_data) // 'Mp_all_array.txt')
output_Rp_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Rp_MCMC_array.txt')
output_Rp_all_array_filename = TRIM(TRIM(output_temp_data) // 'Rp_all_array.txt')
output_Ms_solar_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Ms_solar_MCMC_array.txt')
output_Ms_solar_all_array_filename = TRIM(TRIM(output_temp_data) // 'Ms_solar_all_array.txt')
output_Rs_solar_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Rs_solar_MCMC_array.txt')
output_Rs_solar_all_array_filename = TRIM(TRIM(output_temp_data) // 'Rs_solar_all_array.txt')
output_Rp_Rs_ratio_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Rp_Rs_ratio_MCMC_array.txt')
output_Rp_Rs_ratio_all_array_filename = TRIM(TRIM(output_temp_data) // 'Rp_Rs_ratio_all_array.txt')
output_Ecc_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Ecc_MCMC_array.txt')
output_Ecc_all_array_filename = TRIM(TRIM(output_temp_data) // 'Ecc_all_array.txt')
output_Inc_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'Inc_MCMC_array.txt')
output_Inc_all_array_filename = TRIM(TRIM(output_temp_data) // 'Inc_all_array.txt')
output_omega_arg_periastron_MCMC_array_filename = TRIM(TRIM(output_temp_data) // 'omega_arg_periastron_MCMC_array.txt')
output_omega_arg_periastron_all_array_filename = TRIM(TRIM(output_temp_data) // 'omega_arg_periastron_all_array.txt')
output_MCMC_properties_filename = TRIM(TRIM(output_temp_data) // 'MCMC_properties.txt')
output_best_param_chi_squared_filename = TRIM(TRIM(output_temp_data) // 'best_param_chi_squared.txt')
output_best_param_mean_filename = TRIM(TRIM(output_temp_data) // 'best_param_mean.txt')
output_stand_dev_filename = TRIM(TRIM(output_temp_data) // 'stand_dev.txt')
output_scale_factor_array_filename = TRIM(TRIM(output_temp_data) // 'scale_factor_array.txt')
!*************************************************************************************************

Allocate(Data_my_rv(num_my_rv,3))
Allocate(Data_other_rv(num_other_rv,3))

!Now let's read in the data.
!Read in output data from IDL and assaign values to the data array.
OPEN(unit=99, FILE=input_my_data_filename, status='old', action='read')
!Now read line by line of data file and dump values into the data array.

DO l = 1, num_my_rv
    READ (99, 51) a, b, c
    51 FORMAT(E50.10, 5X, F50.5, 5X, F50.5)
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
        READ (99, 52) a, b, c
        52 FORMAT(E50.10, 5X, F50.5, 5X, F50.5)
        Data_other_rv(l,1) = a
        Data_other_rv(l,2) = b
        Data_other_rv(l,3) = c
    END DO 

    CLOSE(99)
END IF

e = 1
k = 1

vturb = SQRT(beta_rot**2.0D0 + vmacro**2.0D0)            !Velocity width of spectral line due to mechanisms other than rotation (i.e. micro and macro turbulence).

total_interval = datafilelength

!MCMC array allocation. Create arrays that accept all rejected iterations large enough to hold them. Assuming 1% acceptance rate (hopefully it will
!be much higher than this).
Allocate(RV_offset_datasets_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(RV_offset_datasets_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Orbital_period_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Orbital_period_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(JD_time_mid_transit_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(JD_time_mid_transit_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Mp_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Mp_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Rp_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Rp_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Rp_Rs_ratio_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Rp_Rs_ratio_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Rs_solar_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Rs_solar_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Ms_solar_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Ms_solar_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))

Allocate(Ecc_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Ecc_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Inc_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Inc_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(omega_arg_periastron_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(omega_arg_periastron_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(RV_zero_offset_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(RV_zero_offset_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))

Allocate(vsini_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(vsini_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(vsini_variance_current_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(vsini_variance_MCMC_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(vsini_stand_dev_current_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(vsini_stand_dev_MCMC_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))

Allocate(stellar_rotation_angle_accept_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(stellar_rotation_angle_all_mcmc_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(stellar_rotation_angle_variance_current_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(stellar_rotation_angle_variance_MCMC_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(stellar_rotation_angle_stand_dev_current_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(stellar_rotation_angle_stand_dev_MCMC_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(walker_parameter_array(5, number_mcmc_walkers))
Allocate(walker_parameter_array_chi2(14, number_mcmc_walkers), walker_parameter_array_mean(14, number_mcmc_walkers))
Allocate(walker_parameter_array_variance(14, number_mcmc_walkers), walker_parameter_array_stand_dev(14, number_mcmc_walkers))
Allocate(mask_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))

Allocate(RV_offset_datasets_sum2(number_mcmc_walkers))
Allocate(RV_zero_offset_sum2(number_mcmc_walkers))
Allocate(Orbital_period_sum2(number_mcmc_walkers))
Allocate(JD_time_mid_transit_sum2(number_mcmc_walkers))
Allocate(Mp_sum2(number_mcmc_walkers))
Allocate(Rp_sum2(number_mcmc_walkers))
Allocate(Ms_solar_sum2(number_mcmc_walkers))
Allocate(Rs_solar_sum2(number_mcmc_walkers))
Allocate(Rp_Rs_ratio_sum2(number_mcmc_walkers))
Allocate(Ecc_sum2(number_mcmc_walkers))
Allocate(Inc_sum2(number_mcmc_walkers))
Allocate(omega_arg_periastron_sum2(number_mcmc_walkers))
Allocate(vsini_sum2(number_mcmc_walkers))
Allocate(stellar_rotation_angle_sum2(number_mcmc_walkers))

Allocate(Chi_squared_array_MCMC(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(Reduced_Chi_Squared_array_MCMC(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(likelihood_array_MCMC(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(culumative_likelihood_array_MCMC(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))

Allocate(loc_min_chi_squared_MCMC(number_mcmc_walkers), loc_min_reduced_chi_squared_MCMC(number_mcmc_walkers))
Allocate(loc_max_likelihood_MCMC(number_mcmc_walkers), min_chi_squared_MCMC(number_mcmc_walkers))
Allocate(min_reduced_chi_squared_MCMC(number_mcmc_walkers), max_likelihood_MCMC(number_mcmc_walkers))

Allocate(normal1(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal2(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal3(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal4(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal5(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal6(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal7(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal8(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal9(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal10(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal11(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal12(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal13(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(normal14(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(random_draw_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))

Allocate(spread_rand_RV_offset(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Orbital_period(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_JD(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Mp(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Rp(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Rp_Rs_ratio(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Rs(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Ms(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Ecc(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_Inc(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_omega(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_RV_zero(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_vsini(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))
Allocate(spread_rand_spin_orbit(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))

Allocate(scale_factor_array(NINT(mcmc_accepted_iteration_size/0.05D0), number_mcmc_walkers))




call system_clock(count_rate=cr)
call system_clock(count_max=cm)
clock_rate = DBLE(cr)

call system_clock(system_time_1)
call system_clock(system_time_it1)
call system_clock(system_time_it2)

mask_array(1:NINT(mcmc_accepted_iteration_size/0.05D0), 1:number_mcmc_walkers) = .FALSE.

Rs = Rss * Rs_solar_prior
Ms = Mss * Ms_solar_prior 

Rs2 = Rs**2.0D0                              !Square the radius of star to speed up calculations.

IF (linear_quadratic == 'q') THEN
   Io = 6.0D0 / (pi*Rs2*(6.0D0 - (2.0D0*q_1) - q_2))     !The initial light intensity equation with limb darkening (quadratic)
ELSE
   Io = 1.0D0 / (pi*Rs2*(1.0D0-(u/3.0D0)))      !The initial light intensity equation with limb darkening (linear)
                                                !(normalize Io such that total star luminosity is 1). 
END IF

!Based on given priors, determine the approximate time in the simulation to start and finish comparing RV's with model.
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
Transit_length_prior = ((Orbital_period_prior*day_sec)/pi) * asin((Rs/Rorb_prior)*(sqrt((1.0D0 + radius_ratio_prior)**2.0D0 - impact_prior**2.0D0)/sin(Inc_prior*(pi/180.0D0)))) * &
                       (sqrt(1.0D0 - Ecc_prior**2.0D0)/(1.0D0 + (Ecc_prior*sin(omega_arg_periastron_prior*(pi/180.0D0)))))
PRINT *, "Orbital_period_prior: ", Orbital_period_prior
PRINT *, "Orbital_period_prior * day_sec: ", Orbital_period_prior*day_sec
PRINT *, "Rorb_prior: ", Rorb_prior
PRINT *, "Rorb_star_prior: ", Rorb_star_prior
PRINT *, "Inc_prior: ", Inc_prior
PRINT *, "Rp_prior: ", Rp_prior
PRINT *, "Mp_prior: ", Mp_prior
PRINT *, "Ecc_prior: ", Ecc_prior
PRINT *, "omega_arg_periastron_prior: ", omega_arg_periastron_prior
PRINT *, "vsini_prior: ", vsini_prior
PRINT *, "stellar_rotation_angle_prior: ", stellar_rotation_angle_prior
PRINT *, "Rs: ", Rs
PRINT *, "Impact parameter prior: ", impact_prior
PRINT *, "Radius ratio prior: ", radius_ratio_prior
PRINT *, "Length of transit with given priors: ", Transit_length_prior
length_time_compare = (Transit_length_prior/2.0D0) + Time_compare_vel
PRINT *, "Length of time to compare: ", length_time_compare

PRINT *, "number_mcmc_walkers: ", number_mcmc_walkers
PRINT *, "Total potential number of iterations per walker: ", NINT(mcmc_accepted_iteration_size/0.05D0)
PRINT *, "Total potential number of all iterations: ", NINT(mcmc_accepted_iteration_size/0.05D0) * number_mcmc_walkers
PRINT *, "Total number of accepted iterations: ", mcmc_accepted_iteration_size * number_mcmc_walkers

Number_fit = 0
!Determine the number of free parameters.
IF (stellar_rotation_angle_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (vsini_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (RV_offset_datasets_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (RV_zero_offset_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Orbital_period_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (JD_time_mid_transit_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Mp_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Rp_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Ecc_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Inc_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Ms_solar_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Rs_solar_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF

IF (Rp_Rs_ratio_fixed == 'N') THEN
    Number_fit = Number_fit + 1
END IF




!Allocate random numbers outside the do concurrent loop.
DO i = 1, number_mcmc_walkers
    !Generate a random seed for each MCMC walker.
    CALL init_random_seed()
    
    DO j = 1, NINT(mcmc_accepted_iteration_size/0.05D0)
    
        IF ((RV_offset_datasets_fixed == 'N') .AND. (RV_offset_datasets_norm_prior == 'Y')) THEN
            normal1(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((RV_offset_datasets_fixed == 'N') .AND. (RV_offset_datasets_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_RV_offset(j,i) = spread_rand(RV_offset_datasets_begin,RV_offset_datasets_end)
            ELSE
                spread_rand_RV_offset(j,i) = (ran1() - 0.5D0) * abs(RV_offset_datasets_begin - RV_offset_datasets_end)
            END IF
        END IF
        
        IF ((Orbital_period_fixed == 'N') .AND. (Orbital_period_norm_prior == 'Y')) THEN
            normal2(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Orbital_period_fixed == 'N') .AND. (Orbital_period_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Orbital_period(j,i) = spread_rand(Orbital_period_begin,Orbital_period_end)
            ELSE
                spread_rand_Orbital_period(j,i) = (ran1() - 0.5D0) * abs(Orbital_period_begin - Orbital_period_end)
            END IF
        END IF
        
        IF ((JD_time_mid_transit_fixed == 'N') .AND. (JD_time_mid_transit_norm_prior == 'Y')) THEN
            normal3(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((JD_time_mid_transit_fixed == 'N') .AND. (JD_time_mid_transit_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_JD(j,i) = spread_rand(JD_time_mid_transit_begin,JD_time_mid_transit_end)
            ELSE
                spread_rand_JD(j,i) = (ran1() - 0.5D0) * abs(JD_time_mid_transit_begin - JD_time_mid_transit_end)
            END IF
        END IF
        
        IF ((Mp_fixed == 'N') .AND. (Mp_norm_prior == 'Y')) THEN
            normal4(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Mp_fixed == 'N') .AND. (Mp_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Mp(j,i) = spread_rand(Mp_begin,Mp_end)
            ELSE
                spread_rand_Mp(j,i) = (ran1() - 0.5D0) * abs(Mp_begin - Mp_end)
            END IF
        END IF
        
        IF ((Rp_fixed == 'N') .AND. (Rp_norm_prior == 'Y')) THEN
            normal5(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Rp_fixed == 'N') .AND. (Rp_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Rp(j,i) = spread_rand(Rp_begin,Rp_end)
            ELSE
                spread_rand_Rp(j,i) = (ran1() - 0.5D0) * abs(Rp_begin - Rp_end)
            END IF
        END IF
        
        IF ((Rp_Rs_ratio_fixed == 'N') .AND. (Rp_Rs_ratio_norm_prior == 'Y')) THEN
            normal6(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Rp_Rs_ratio_fixed == 'N') .AND. (Rp_Rs_ratio_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Rp_Rs_ratio(j,i) = spread_rand(Rp_Rs_ratio_begin,Rp_Rs_ratio_end)
            ELSE
                spread_rand_Rp_Rs_ratio(j,i) = (ran1() - 0.5D0) * abs(Rp_Rs_ratio_begin - Rp_Rs_ratio_end)
            END IF
        END IF
        
        IF ((Rs_solar_fixed == 'N') .AND. (Rs_solar_norm_prior == 'Y')) THEN
            normal7(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Rs_solar_fixed == 'N') .AND. (Rs_solar_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Rs(j,i) = spread_rand(Rs_solar_begin,Rs_solar_end)
            ELSE
                spread_rand_Rs(j,i) = (ran1() - 0.5D0) * abs(Rs_solar_begin - Rs_solar_end)
            END IF
        END IF
        
        IF ((Ms_solar_fixed == 'N') .AND. (Ms_solar_norm_prior == 'Y')) THEN
            normal8(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Ms_solar_fixed == 'N') .AND. (Ms_solar_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Ms(j,i) = spread_rand(Ms_solar_begin,Ms_solar_end)
            ELSE
                spread_rand_Ms(j,i) = (ran1() - 0.5D0) * abs(Ms_solar_begin - Ms_solar_end)
            END IF
        END IF
        
        IF ((Ecc_fixed == 'N') .AND. (Ecc_norm_prior == 'Y')) THEN
            normal9(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Ecc_fixed == 'N') .AND. (Ecc_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Ecc(j,i) = spread_rand(Ecc_begin,Ecc_end)
            ELSE
                spread_rand_Ecc(j,i) = (ran1() - 0.5D0) * abs(Ecc_begin - Ecc_end)
            END IF
        END IF
        
        IF ((Inc_fixed == 'N') .AND. (Inc_norm_prior == 'Y')) THEN
            normal10(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((Inc_fixed == 'N') .AND. (Inc_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_Inc(j,i) = spread_rand(Inc_begin,Inc_end)
            ELSE
                spread_rand_Inc(j,i) = (ran1() - 0.5D0) * abs(Inc_begin - Inc_end)
            END IF
        END IF
        
        IF ((omega_arg_periastron_fixed == 'N') .AND. (omega_arg_periastron_norm_prior == 'Y')) THEN
            normal11(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((omega_arg_periastron_fixed == 'N') .AND. (omega_arg_periastron_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_omega(j,i) = spread_rand(omega_arg_periastron_begin,omega_arg_periastron_end)
            ELSE
                spread_rand_omega(j,i) = (ran1() - 0.5D0) * abs(omega_arg_periastron_begin - omega_arg_periastron_end)
            END IF
        END IF
        
        IF ((RV_zero_offset_fixed == 'N') .AND. (RV_zero_offset_norm_prior == 'Y')) THEN
            normal12(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((RV_zero_offset_fixed == 'N') .AND. (RV_zero_offset_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_RV_zero(j,i) = spread_rand(RV_zero_offset_begin,RV_zero_offset_end)
            ELSE
                spread_rand_RV_zero(j,i) = (ran1() - 0.5D0) * abs(RV_zero_offset_begin - RV_zero_offset_end)
            END IF
        END IF
        
        IF ((vsini_fixed == 'N') .AND. (vsini_norm_prior == 'Y')) THEN
            normal13(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((vsini_fixed == 'N') .AND. (vsini_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_vsini(j,i) = spread_rand(vsini_begin,vsini_end)
            ELSE
                spread_rand_vsini(j,i) = (ran1() - 0.5D0) * abs(vsini_begin - vsini_end)
            END IF
        END IF
        
        IF ((stellar_rotation_angle_fixed == 'N') .AND. (stellar_rotation_angle_norm_prior == 'Y')) THEN
            normal14(j,i) = normal(0.0D0, 1.0D0)
        ELSE IF ((stellar_rotation_angle_fixed == 'N') .AND. (stellar_rotation_angle_norm_prior == 'N')) THEN
            IF (j ==1) THEN
                spread_rand_spin_orbit(j,i) = spread_rand(stellar_rotation_angle_begin,stellar_rotation_angle_end)
            ELSE
                spread_rand_spin_orbit(j,i) = (ran1() - 0.5D0) * abs(stellar_rotation_angle_begin - stellar_rotation_angle_end)
            END IF
        END IF
        
        random_draw_array(j,i) = ran1()
        
        
        

        percent_processed = (100.0D0 * (NINT(mcmc_accepted_iteration_size/0.05D0) * (i - 1.0D0) + j))/(NINT(mcmc_accepted_iteration_size/0.05D0) * number_mcmc_walkers)
        string_bar = ''
        string_bar_num = 1
        IF ((percent_processed < 4.0D0) .AND. (percent_processed > 0.05D0)) THEN
            IF ((percent_processed - FLOOR(percent_processed) < 0.2499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.0000001D0)) THEN
                string_bar(string_bar_num:) = ':'
            END IF
            IF ((percent_processed - FLOOR(percent_processed) < 0.4999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.2499999D0)) THEN
                string_bar(string_bar_num:) = '/'
            END IF
            IF ((percent_processed - FLOOR(percent_processed) < 0.7499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.4999999D0)) THEN
                string_bar(string_bar_num:) = '-'
            END IF
            IF ((percent_processed - FLOOR(percent_processed) < 0.9999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.7499999D0)) THEN
                string_bar(string_bar_num:) = '\'
            END IF
        END IF   
            
        IF (percent_processed >= 4.0D0) THEN
            DO percent_processed_loop = 1, FLOOR(percent_processed/4.0D0)
                string_bar(string_bar_num:) = '|'
                string_bar_num = string_bar_num + 1
                IF ((percent_processed < 99.9D0) .AND. (percent_processed_loop == FLOOR(percent_processed/4.0D0))) THEN
                    IF ((percent_processed - FLOOR(percent_processed) < 0.2499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.0000001D0)) THEN
                        string_bar(string_bar_num:) = ':'
                    END IF
                    IF ((percent_processed - FLOOR(percent_processed) < 0.4999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.2499999D0)) THEN
                        string_bar(string_bar_num:) = '/'
                    END IF
                    IF ((percent_processed - FLOOR(percent_processed) < 0.7499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.4999999D0)) THEN
                        string_bar(string_bar_num:) = '-'
                    END IF
                    IF ((percent_processed - FLOOR(percent_processed) < 0.9999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.7499999D0)) THEN
                        string_bar(string_bar_num:) = '\'
                    END IF
                END IF
            END DO
        END IF 
            
        PRINT 500, string_bar, FLOOR(percent_processed)
        500 FORMAT("Drawing random numbers: [", A25, "]", 1X, I3, "% complete")

        IF (percent_processed < 4.0) THEN
            PRINT *, 'Estimated time remaining: ', 'calculating....'
        ENDIF
        IF (percent_processed >= 4.0) THEN
            IF (FLOOR(percent_processed) - (percent_count*4) >= 0) THEN
                call system_clock(system_time_it2)
                time_remaining = (((system_time_it2 - system_time_it1)/clock_rate)/percent_processed)*(100.0D0 - percent_processed)
                percent_count = percent_count + 1
            END IF
            WRITE (*,501) 'Estimated time remaining: ', time_remaining/60.0D0, ' minutes'
            501 FORMAT(A30,F8.2,A8)
        END IF
        
        WRITE (*,502) 'Iteration number: ', (i-1)*NINT(mcmc_accepted_iteration_size/0.05D0) + j, &
        ' out of ', NINT(mcmc_accepted_iteration_size/0.05) * number_mcmc_walkers, ' interations'
        502 FORMAT(A20,I8,A8,I8,A12)

        PRINT *," "
        PRINT *,"***********************************************************************************************"
        PRINT *," "
    
    END DO
END DO




!MCMC walkers. Parallelize each MCMC walker.
!************************************************************************************************************************************************
DO CONCURRENT (walker = 1:number_mcmc_walkers)

    !MCMC routines
    !************************************************************************************************************************************************
    !This MCMC routine first draws random values for the priors based on their reported literature values. Then values for the proposal parameters
    !are randomly drawn based on previously accepted proposal values. The MCMC performs a random walk on the proposal parameters. It does not do a
    !random walk for the priors.

    !Rejection flag for the proposal parameters. Initialize as false.
    Reject_flag = 'F'
    MCMC_loop = 1     !Initialize MCMC_loop to 1
    Number_iterations_left_total = mcmc_accepted_iteration_size
    string_bar_num = 1 
    reject_counter = 0                                 
    all_loop = 1                                        
    MCMC_print_updates = 0                             
    MCMC_updates = 0                                              
    count_tested_proposals = 0
    accept_rate = 0.0D0
    scale_factor = scale_factor_in
    Ecc_draw_loop = 1
    vsini_draw_loop = 1
    RV_offset_datasets_draw_loop = 1
    Orbital_period_draw_loop = 1
    JD_time_mid_transit_draw_loop = 1
    Mp_draw_loop = 1
    Rp_draw_loop = 1
    Rp_Rs_ratio_draw_loop = 1
    Rs_solar_draw_loop = 1
    Ms_solar_draw_loop = 1
    Inc_draw_loop = 1
    omega_arg_periastron_draw_loop = 1
    RV_zero_offset_draw_loop = 1
    stellar_rotation_angle_draw_loop = 1
                         
    
    
   
    !mcmc_accepted_iteration_size is the sample size of the accepted proposal runs (i.e., accepted main parameter values). This loop will run
    !until the number of accepted proposals equals the mcmc_accepted_iteration_size. Note, this does not guarentee convergence or well mixing
    !of the mcmc.
    !************************************************************** MCMC start ***********************************************************************************!
    DO WHILE (MCMC_loop <= mcmc_accepted_iteration_size)
   
        !*********************************************** Draw values for proposal parameters start *******************************************************!
        !First check whether this is the first iteration. If it is, then randomaly select proposal values from the guess at the mean
        !and variance. Then check whether the rejection flag is true or false. For the second iteration, the rejection flag will be
        !false. If rejection flag is false, use previous proposal values and variance to generate new proposal values. If rejection
        !flag is true, the proposal values from the previous iteration were rejected. Generate new random proposal values from the
        !previously accepted proposal values.
   
        IF (all_loop == 1) THEN
            !Check if this is the first iteration.  
        
            RV_offset_datasets_variance = RV_offset_datasets_1sigerr**2.0D0
            RV_offset_datasets_variance_MCMC = RV_offset_datasets_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            RV_offset_datasets_stand_dev = SQRT(RV_offset_datasets_variance)
            RV_offset_datasets_stand_dev_MCMC = SQRT(RV_offset_datasets_variance_MCMC)  
            
            IF (RV_offset_datasets_fixed == 'Y') THEN
                RV_offset_datasets = RV_offset_datasets_prior
            ELSE IF ((RV_offset_datasets_fixed == 'N') .AND. (RV_offset_datasets_norm_prior == 'Y')) THEN
                RV_offset_datasets = RV_offset_datasets_prior + (RV_offset_datasets_stand_dev_MCMC*scale_factor*normal1(all_loop, walker))
            ELSE IF ((RV_offset_datasets_fixed == 'N') .AND. (RV_offset_datasets_norm_prior == 'N')) THEN
                !Draw randomaly between min and max values.
                RV_offset_datasets = spread_rand_RV_offset(all_loop, walker)
            END IF
            
            RV_offset_datasets_sum = RV_offset_datasets
            RV_offset_datasets_mean = RV_offset_datasets
            RV_offset_datasets_diff_mean_2 = 0.0D0
            
            RV_offset_datasets_accept_mcmc_array(all_loop, walker) = RV_offset_datasets
            RV_offset_datasets_all_mcmc_array(all_loop, walker) = RV_offset_datasets 
            
            
            
            
            Orbital_period_variance = Orbital_period_1sigerr**2.0D0
            Orbital_period_variance_MCMC = Orbital_period_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Orbital_period_stand_dev = SQRT(Orbital_period_variance)
            Orbital_period_stand_dev_MCMC = SQRT(Orbital_period_variance_MCMC)
            
            IF (Orbital_period_fixed == 'Y') THEN
                Orbital_period_1 = Orbital_period_prior
            ELSE IF ((Orbital_period_fixed == 'N') .AND. (Orbital_period_norm_prior == 'Y')) THEN
                Orbital_period_1 = Orbital_period_prior + (Orbital_period_stand_dev_MCMC*scale_factor*normal2(all_loop, walker))
            ELSE IF ((Orbital_period_fixed == 'N') .AND. (Orbital_period_norm_prior == 'N')) THEN
                !Draw randomaly between min and max values.
                Orbital_period_1 = spread_rand_Orbital_period(all_loop, walker)
            END IF
            
            Orbital_period_sum = Orbital_period_1
            Orbital_period_mean = Orbital_period_1
            Orbital_period_diff_mean_2 = 0.0D0
            Orbital_period_accept_mcmc_array(all_loop, walker) = Orbital_period_1
            Orbital_period_all_mcmc_array(all_loop, walker) = Orbital_period_1 
            Orbital_period = Orbital_period_1*day_sec               !Convert orbital period in days to seconds.  
            
            
            
            
            JD_time_mid_transit_variance = JD_time_mid_transit_1sigerr**2.0D0
            JD_time_mid_transit_variance_MCMC = JD_time_mid_transit_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            JD_time_mid_transit_stand_dev = SQRT(JD_time_mid_transit_variance)
            JD_time_mid_transit_stand_dev_MCMC = SQRT(JD_time_mid_transit_variance_MCMC)

            IF (JD_time_mid_transit_fixed == 'Y') THEN
                JD_time_mid_transit = JD_time_mid_transit_prior
            ELSE IF ((JD_time_mid_transit_fixed == 'N') .AND. (JD_time_mid_transit_norm_prior == 'Y')) THEN
                JD_time_mid_transit = JD_time_mid_transit_prior + (JD_time_mid_transit_stand_dev_MCMC*scale_factor*normal3(all_loop, walker))
            ELSE IF ((JD_time_mid_transit_fixed == 'N') .AND. (JD_time_mid_transit_norm_prior == 'N')) THEN
                JD_time_mid_transit = spread_rand_JD(all_loop, walker)
            END IF
            
            JD_time_mid_transit_sum = JD_time_mid_transit
            JD_time_mid_transit_mean = JD_time_mid_transit
            JD_time_mid_transit_diff_mean_2 = 0.0D0
            JD_time_mid_transit_accept_mcmc_array(all_loop, walker) = JD_time_mid_transit
            JD_time_mid_transit_all_mcmc_array(all_loop, walker) = JD_time_mid_transit
            
            
            
            
            Mp_variance = Mp_1sigerr**2.0D0
            Mp_variance_MCMC = Mp_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Mp_stand_dev = SQRT(Mp_variance)
            Mp_stand_dev_MCMC = SQRT(Mp_variance_MCMC)

            IF (Mp_fixed == 'Y') THEN
                Mpp = Mp_prior
            ELSE IF ((Mp_fixed == 'N') .AND. (Mp_norm_prior == 'Y')) THEN
                Mpp = Mp_prior + (Mp_stand_dev_MCMC*scale_factor*normal4(all_loop, walker))
            ELSE IF ((Mp_fixed == 'N') .AND. (Mp_norm_prior == 'N')) THEN
                Mpp = spread_rand_Mp(all_loop, walker)
            END IF
            Mp_sum = Mpp
            Mp_mean = Mpp
            Mp_diff_mean_2 = 0.0D0
            Mp_accept_mcmc_array(all_loop, walker) = Mpp
            Mp_all_mcmc_array(all_loop, walker) = Mpp
            
            
            
            
            Rp_variance = Rp_1sigerr**2.0D0
            Rp_variance_MCMC = Rp_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Rp_stand_dev = SQRT(Rp_variance)
            Rp_stand_dev_MCMC = SQRT(Rp_variance_MCMC)

            IF (Rp_fixed == 'Y') THEN
                Rpp = Rp_prior
            ELSE IF ((Rp_fixed == 'N') .AND. (Rp_norm_prior == 'Y')) THEN
                Rpp = Rp_prior + (Rp_stand_dev_MCMC*scale_factor*normal5(all_loop, walker))
            ELSE IF ((Rp_fixed == 'N') .AND. (Rp_norm_prior == 'N')) THEN
                Rpp = spread_rand_Rp(all_loop, walker)
            END IF
            Rp_sum = Rpp
            Rp_mean = Rpp
            Rp_diff_mean_2 = 0.0D0
            Rp_accept_mcmc_array(all_loop, walker) = Rpp
            Rp_all_mcmc_array(all_loop, walker) = Rpp
            
            
            
            
            Rp_Rs_ratio_variance = Rp_Rs_ratio_1sigerr**2.0D0
            Rp_Rs_ratio_variance_MCMC = Rp_Rs_ratio_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Rp_Rs_ratio_stand_dev = SQRT(Rp_Rs_ratio_variance)
            Rp_Rs_ratio_stand_dev_MCMC = SQRT(Rp_Rs_ratio_variance_MCMC)

            IF (Rp_Rs_ratio_fixed == 'Y') THEN
                Rp_Rs_ratio = Rp_Rs_ratio_prior
            ELSE IF ((Rp_Rs_ratio_fixed == 'N') .AND. (Rp_Rs_ratio_norm_prior == 'Y')) THEN
                Rp_Rs_ratio = Rp_Rs_ratio_prior + (Rp_Rs_ratio_stand_dev_MCMC*scale_factor*normal6(all_loop, walker))
            ELSE IF ((Rp_Rs_ratio_fixed == 'N') .AND. (Rp_Rs_ratio_norm_prior == 'N')) THEN
                Rp_Rs_ratio = spread_rand_Rp_Rs_ratio(all_loop, walker)
            END IF
            Rp_Rs_ratio_sum = Rp_Rs_ratio
            Rp_Rs_ratio_mean = Rp_Rs_ratio
            Rp_Rs_ratio_diff_mean_2 = 0.0D0
            Rp_Rs_ratio_accept_mcmc_array(all_loop, walker) = Rp_Rs_ratio
            Rp_Rs_ratio_all_mcmc_array(all_loop, walker) = Rp_Rs_ratio
            
            
            
            
            Rs_solar_variance = Rs_solar_1sigerr**2.0D0
            Rs_solar_variance_MCMC = Rs_solar_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Rs_solar_stand_dev = SQRT(Rs_solar_variance)
            Rs_solar_stand_dev_MCMC = SQRT(Rs_solar_variance_MCMC)

            IF (Rs_solar_fixed == 'Y') THEN
                Rs_solar = Rs_solar_prior
            ELSE IF ((Rs_solar_fixed == 'N') .AND. (Rs_solar_norm_prior == 'Y')) THEN
                Rs_solar = Rs_solar_prior + (Rs_solar_stand_dev_MCMC*scale_factor*normal7(all_loop, walker))
            ELSE IF ((Rs_solar_fixed == 'N') .AND. (Rs_solar_norm_prior == 'N')) THEN
                Rs_solar = spread_rand_Rs(all_loop, walker)
            END IF
            Rs_solar_sum = Rs_solar
            Rs_solar_mean = Rs_solar
            Rs_solar_diff_mean_2 = 0.0D0
            Rs_solar_accept_mcmc_array(all_loop, walker) = Rs_solar
            Rs_solar_all_mcmc_array(all_loop, walker) = Rs_solar
            
            
            
            
            Ms_solar_variance = Ms_solar_1sigerr**2.0D0
            Ms_solar_variance_MCMC = Ms_solar_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Ms_solar_stand_dev = SQRT(Ms_solar_variance)
            Ms_solar_stand_dev_MCMC = SQRT(Ms_solar_variance_MCMC)

            IF (Ms_solar_fixed == 'Y') THEN
                Ms_solar = Ms_solar_prior
            ELSE IF ((Ms_solar_fixed == 'N') .AND. (Ms_solar_norm_prior == 'Y')) THEN
                Ms_solar = Ms_solar_prior + (Ms_solar_stand_dev_MCMC*scale_factor*normal8(all_loop, walker))
            ELSE IF ((Ms_solar_fixed == 'N') .AND. (Ms_solar_norm_prior == 'N')) THEN
                Ms_solar = spread_rand_Ms(all_loop, walker)
            END IF
            Ms_solar_sum = Ms_solar
            Ms_solar_mean = Ms_solar
            Ms_solar_diff_mean_2 = 0.0D0
            Ms_solar_accept_mcmc_array(all_loop, walker) = Ms_solar
            Ms_solar_all_mcmc_array(all_loop, walker) = Ms_solar

            
            
            
            Ecc_variance = Ecc_1sigerr**2.0D0
            Ecc_variance_MCMC = Ecc_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Ecc_stand_dev = SQRT(Ecc_variance)
            Ecc_stand_dev_MCMC = SQRT(Ecc_variance_MCMC)

            IF (Ecc_fixed == 'Y') THEN
                Ecc = Ecc_prior
            ELSE IF ((Ecc_fixed == 'N') .AND. (Ecc_norm_prior == 'Y')) THEN
                Ecc = Ecc_prior + (Ecc_stand_dev_MCMC*scale_factor*normal9(ecc_draw_loop, walker))
                !If a negative Ecc is drawn, keep drawing an Ecc value until it is positive. Cannot have a negative Ecc.
                IF (Ecc < 0) THEN
                    DO WHILE (Ecc < 0.0D0)
                        ecc_draw_loop = ecc_draw_loop + 1
                        Ecc = Ecc_prior + (Ecc_stand_dev_MCMC*scale_factor*normal9(ecc_draw_loop, walker))
                    END DO
                END IF
            ELSE IF ((Ecc_fixed == 'N') .AND. (Ecc_norm_prior == 'N')) THEN
                Ecc = spread_rand_Ecc(ecc_draw_loop, walker)
            END IF
            Ecc_sum = Ecc
            Ecc_mean = Ecc
            Ecc_diff_mean_2 = 0.0D0
            Ecc_accept_mcmc_array(all_loop, walker) = Ecc
            Ecc_all_mcmc_array(all_loop, walker) = Ecc
            
            
            
            
            Inc_variance = Inc_1sigerr**2.0D0
            Inc_variance_MCMC = Inc_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            Inc_stand_dev = SQRT(Inc_variance)
            Inc_stand_dev_MCMC = SQRT(Inc_variance_MCMC)

            IF (Inc_fixed == 'Y') THEN
                Inc = Inc_prior
            ELSE IF ((Inc_fixed == 'N') .AND. (Inc_norm_prior == 'Y')) THEN
                Inc = Inc_prior + (Inc_stand_dev_MCMC*scale_factor*normal10(all_loop, walker))
                IF (Inc > 90.0D0) THEN
                    Inc = 90.0D0 - abs(Inc - 90.0D0)
                END IF
            ELSE IF ((Inc_fixed == 'N') .AND. (Inc_norm_prior == 'N')) THEN
                Inc = spread_rand_Inc(all_loop, walker)
            END IF
            Inc_sum = Inc
            Inc_mean = Inc
            Inc_diff_mean_2 = 0.0D0
            Inc_accept_mcmc_array(all_loop, walker) = Inc
            Inc_all_mcmc_array(all_loop, walker) = Inc
            
            
            
            
            omega_arg_periastron_variance = omega_arg_periastron_1sigerr**2.0D0
            omega_arg_periastron_variance_MCMC = omega_arg_periastron_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            omega_arg_periastron_stand_dev = SQRT(omega_arg_periastron_variance)
            omega_arg_periastron_stand_dev_MCMC = SQRT(omega_arg_periastron_variance_MCMC)

            IF (omega_arg_periastron_fixed == 'Y') THEN
                omega_arg_periastron = omega_arg_periastron_prior
            ELSE IF ((omega_arg_periastron_fixed == 'N') .AND. (omega_arg_periastron_norm_prior == 'Y')) THEN
                omega_arg_periastron = omega_arg_periastron_prior + (omega_arg_periastron_stand_dev_MCMC*scale_factor*normal11(all_loop, walker))
            ELSE IF ((omega_arg_periastron_fixed == 'N') .AND. (omega_arg_periastron_norm_prior == 'N')) THEN
                omega_arg_periastron = spread_rand_omega(all_loop, walker)
            END IF
            omega_arg_periastron_sum = omega_arg_periastron
            omega_arg_periastron_mean = omega_arg_periastron
            omega_arg_periastron_diff_mean_2 = 0.0D0
            omega_arg_periastron_accept_mcmc_array(all_loop, walker) = omega_arg_periastron
            omega_arg_periastron_all_mcmc_array(all_loop, walker) = omega_arg_periastron
            
            
            
            
            RV_zero_offset_variance = RV_zero_offset_1sigerr**2.0D0
            RV_zero_offset_variance_MCMC = RV_zero_offset_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            RV_zero_offset_stand_dev = SQRT(RV_zero_offset_variance)
            RV_zero_offset_stand_dev_MCMC = SQRT(RV_zero_offset_variance_MCMC)

            IF (RV_zero_offset_fixed == 'Y') THEN
                RV_zero_offset = RV_zero_offset_prior
            ELSE IF ((RV_zero_offset_fixed == 'N') .AND. (RV_zero_offset_norm_prior == 'Y')) THEN
                RV_zero_offset = RV_zero_offset_prior + (RV_zero_offset_stand_dev_MCMC*scale_factor*normal12(all_loop, walker))
            ELSE IF ((RV_zero_offset_fixed == 'N') .AND. (RV_zero_offset_norm_prior == 'Y')) THEN
                RV_zero_offset = spread_rand_RV_zero(all_loop, walker)
            END IF
            RV_zero_offset_sum = RV_zero_offset
            RV_zero_offset_mean = RV_zero_offset
            RV_zero_offset_diff_mean_2 = 0.0D0
            RV_zero_offset_accept_mcmc_array(all_loop, walker) = RV_zero_offset
            RV_zero_offset_all_mcmc_array(all_loop, walker) = RV_zero_offset
            
            
            
            
            !Set the vsini_variance_MCMC equal to vsini_1sigerr (the initial guess at the vsini uncertainty) for the first 100 iterations.
            !Set the stellar_rotation_angle_variance_MCMC equal to stellar_rotation_angle_1sigerr (the initial guess at the vsini uncertainty) for the first 100 iterations.
            !scale_factor is used to ensure acceptance rate of ~25%. The initial value for the scale_factor is set by the user.
            !The scale factor will be updated every 100 successful proposal draws (i.e., 100 accepted proposal draws).
            vsini_variance = vsini_1sigerr**2.0D0
            vsini_variance_MCMC = vsini_1sigerr**2.0D0       !Update after every 100 accepted proposal draws.
            vsini_stand_dev = SQRT(vsini_variance)
            vsini_stand_dev_MCMC = SQRT(vsini_variance_MCMC)
            IF (vsini_fixed == 'Y') THEN
                vsini = vsini_prior
            ELSE IF ((vsini_fixed == 'N') .AND. (vsini_norm_prior == 'Y')) THEN
                vsini = vsini_prior + (vsini_stand_dev_MCMC*scale_factor*normal13(vsini_draw_loop, walker))
                !If a negative vsini is drawn, keep drawing a vsini value until it is positive. Cannot have a negative vsini.
                IF (vsini < 0) THEN
                    DO WHILE (vsini < 0.0D0)
                        vsini_draw_loop = vsini_draw_loop + 1
                        vsini = vsini_prior + (vsini_stand_dev_MCMC*scale_factor*normal13(vsini_draw_loop, walker))
                    END DO
                END IF
            ELSE IF ((vsini_fixed == 'N') .AND. (vsini_norm_prior == 'Y')) THEN
                vsini = spread_rand_vsini(vsini_draw_loop, walker)
            END IF
            vsini_sum = vsini
            vsini_mean = vsini
            vsini_diff_mean_2 = 0.0D0
            vsini_accept_mcmc_array(all_loop, walker) = vsini
            vsini_all_mcmc_array(all_loop, walker) = vsini
            vsini_variance_current_array(all_loop, walker) = vsini_variance             !The variance updated at every iteration.
            vsini_variance_MCMC_array(all_loop, walker) = vsini_variance_MCMC           !The variance updated every 100 iterations.
            vsini_stand_dev_current_array(all_loop, walker) = vsini_stand_dev           !The standard deviation updated at every iteration.
            vsini_stand_dev_MCMC_array(all_loop, walker) = vsini_stand_dev_MCMC         !The standard deviation updated every 100 iterations.            
            
            
            
            
            stellar_rotation_angle_variance = stellar_rotation_angle_1sigerr**2.0D0
            stellar_rotation_angle_variance_MCMC = stellar_rotation_angle_1sigerr**2.0D0      !Update after every 100 accepted proposal draws.
            stellar_rotation_angle_stand_dev = SQRT(stellar_rotation_angle_variance)
            stellar_rotation_angle_stand_dev_MCMC = SQRT(stellar_rotation_angle_variance_MCMC)
            IF (stellar_rotation_angle_fixed == 'Y') THEN
                stellar_rotation_angle = stellar_rotation_angle_prior
            ELSE IF ((stellar_rotation_angle_fixed == 'N') .AND. (stellar_rotation_angle_norm_prior == 'Y')) THEN
                stellar_rotation_angle = stellar_rotation_angle_prior + (stellar_rotation_angle_stand_dev_MCMC*scale_factor*normal14(all_loop, walker)) 
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
            ELSE IF ((stellar_rotation_angle_fixed == 'N') .AND. (stellar_rotation_angle_norm_prior == 'N')) THEN
                stellar_rotation_angle = spread_rand_spin_orbit(all_loop, walker)
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
            END IF
            stellar_rotation_angle_sum = stellar_rotation_angle
            stellar_rotation_angle_mean = stellar_rotation_angle
            stellar_rotation_angle_diff_mean_2 = 0.0D0
            stellar_rotation_angle_accept_mcmc_array(all_loop, walker) = stellar_rotation_angle
            stellar_rotation_angle_all_mcmc_array(all_loop, walker) = stellar_rotation_angle
            stellar_rotation_angle_variance_current_array(all_loop, walker) = stellar_rotation_angle_variance         !The variance updated at every iteration.
            stellar_rotation_angle_variance_MCMC_array(all_loop, walker) = stellar_rotation_angle_variance_MCMC       !The variance updated every 100 iterations.
            stellar_rotation_angle_stand_dev_current_array(all_loop, walker) = stellar_rotation_angle_stand_dev           !The standard deviation updated at every iteration.
            stellar_rotation_angle_stand_dev_MCMC_array(all_loop, walker) = stellar_rotation_angle_stand_dev_MCMC         !The standard deviation updated every 100 iterations.
            
            
            

        ELSE IF (all_loop >= 2) THEN
            !Use the previously accepted proposal values of the priors as the input to determine the proposal values for this iteration.      
            !RV_offset_datasets = RV_offset_datasets_accept_mcmc_array(MCMC_loop - 1) + (RV_offset_datasets_stand_dev_MCMC*scale_factor*normal(0.0D0, 1.0D0))           
            
            IF (RV_offset_datasets_fixed == 'Y') THEN
                RV_offset_datasets = RV_offset_datasets_prior
            ELSE IF ((RV_offset_datasets_fixed == 'N') .AND. (RV_offset_datasets_norm_prior == 'Y')) THEN
                RV_offset_datasets = RV_offset_datasets_accept_mcmc_array(all_loop - 1, walker) + (RV_offset_datasets_stand_dev_MCMC*scale_factor*normal1(all_loop, walker))
            ELSE IF ((RV_offset_datasets_fixed == 'N') .AND. (RV_offset_datasets_norm_prior == 'N')) THEN
                RV_offset_datasets = RV_offset_datasets_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_RV_offset(RV_offset_datasets_draw_loop, walker))
                IF ((RV_offset_datasets < RV_offset_datasets_begin) .OR. (RV_offset_datasets > RV_offset_datasets_end)) THEN
                    DO WHILE ((RV_offset_datasets < RV_offset_datasets_begin) .OR. (RV_offset_datasets > RV_offset_datasets_end))
                        RV_offset_datasets_draw_loop = RV_offset_datasets_draw_loop + 1
                        RV_offset_datasets = RV_offset_datasets_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_RV_offset(RV_offset_datasets_draw_loop, walker))
                    END DO
                END IF
            END IF
            RV_offset_datasets_all_mcmc_array(all_loop, walker) = RV_offset_datasets 
      
          
      

            IF (Orbital_period_fixed == 'Y') THEN
                Orbital_period_1 = Orbital_period_prior
            ELSE IF ((Orbital_period_fixed == 'N') .AND. (Orbital_period_norm_prior == 'Y')) THEN
                Orbital_period_1 = Orbital_period_accept_mcmc_array(all_loop - 1, walker) + (Orbital_period_stand_dev_MCMC*scale_factor*normal2(all_loop, walker))
            ELSE IF ((Orbital_period_fixed == 'N') .AND. (Orbital_period_norm_prior == 'N')) THEN
                Orbital_period_1 = Orbital_period_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Orbital_period(Orbital_period_draw_loop, walker))
                IF ((Orbital_period_1 < Orbital_period_begin) .OR. (Orbital_period_1 > Orbital_period_end)) THEN
                    DO WHILE ((Orbital_period_1 < Orbital_period_begin) .OR. (Orbital_period_1 > Orbital_period_end))
                        Orbital_period_draw_loop = Orbital_period_draw_loop + 1
                        Orbital_period_1 = Orbital_period_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Orbital_period(Orbital_period_draw_loop, walker))
                    END DO
                END IF
            END IF
            Orbital_period_all_mcmc_array(all_loop, walker) = Orbital_period_1 
            Orbital_period = Orbital_period_1*day_sec               !Convert orbital period in days to seconds.
            
            
            
            
            IF (JD_time_mid_transit_fixed == 'Y') THEN
                JD_time_mid_transit = JD_time_mid_transit_prior
            ELSE IF ((JD_time_mid_transit_fixed == 'N') .AND. (JD_time_mid_transit_norm_prior == 'Y')) THEN
                JD_time_mid_transit = JD_time_mid_transit_accept_mcmc_array(all_loop - 1, walker) + (JD_time_mid_transit_stand_dev_MCMC*scale_factor*normal3(all_loop, walker))
            ELSE IF ((JD_time_mid_transit_fixed == 'N') .AND. (JD_time_mid_transit_norm_prior == 'N')) THEN
                JD_time_mid_transit = JD_time_mid_transit_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_JD(JD_time_mid_transit_draw_loop, walker))
                IF ((JD_time_mid_transit < JD_time_mid_transit_begin) .OR. (JD_time_mid_transit > JD_time_mid_transit_end)) THEN
                    DO WHILE ((JD_time_mid_transit < JD_time_mid_transit_begin) .OR. (JD_time_mid_transit > JD_time_mid_transit_end))
                        JD_time_mid_transit_draw_loop = JD_time_mid_transit_draw_loop + 1
                        JD_time_mid_transit = JD_time_mid_transit_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_JD(JD_time_mid_transit_draw_loop, walker))
                    END DO
                END IF
            END IF
            JD_time_mid_transit_all_mcmc_array(all_loop, walker) = JD_time_mid_transit      
      

      

            IF (Mp_fixed == 'Y') THEN
                Mpp = Mp_prior
            ELSE IF ((Mp_fixed == 'N') .AND. (Mp_norm_prior == 'Y')) THEN
                Mpp = Mp_accept_mcmc_array(all_loop - 1, walker) + (Mp_stand_dev_MCMC*scale_factor*normal4(all_loop, walker))
            ELSE IF ((Mp_fixed == 'N') .AND. (Mp_norm_prior == 'N')) THEN
                Mpp = Mp_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Mp(Mp_draw_loop, walker))
                IF ((Mpp < Mp_begin) .OR. (Mpp > Mp_end)) THEN
                    DO WHILE ((Mpp < Mp_begin) .OR. (Mpp > Mp_end))
                        Mp_draw_loop = Mp_draw_loop + 1
                        Mpp = Mp_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Mp(Mp_draw_loop, walker))
                    END DO
                END IF
            END IF
            Mp_all_mcmc_array(all_loop, walker) = Mpp

      
      

            IF (Rp_fixed == 'Y') THEN
                Rpp = Rp_prior
            ELSE IF ((Rp_fixed == 'N') .AND. (Rp_norm_prior == 'Y')) THEN
                Rpp = Rp_accept_mcmc_array(all_loop - 1, walker) + (Rp_stand_dev_MCMC*scale_factor*normal5(all_loop, walker))
            ELSE IF ((Rp_fixed == 'N') .AND. (Rp_norm_prior == 'N')) THEN
                Rpp = Rp_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Rp(Rp_draw_loop, walker))
                IF ((Rpp < Rp_begin) .OR. (Rpp > Rp_end)) THEN
                    DO WHILE ((Rpp < Rp_begin) .OR. (Rpp > Rp_end))
                        Rp_draw_loop = Rp_draw_loop + 1
                        Rpp = Rp_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Rp(Rp_draw_loop, walker))
                    END DO
                END IF
            END IF
            Rp_all_mcmc_array(all_loop, walker) = Rpp
      
      
      

            IF (Rp_Rs_ratio_fixed == 'Y') THEN
                Rp_Rs_ratio = Rp_Rs_ratio_prior
            ELSE IF ((Rp_Rs_ratio_fixed == 'N') .AND. (Rp_Rs_ratio_norm_prior == 'Y')) THEN
                Rp_Rs_ratio = Rp_Rs_ratio_accept_mcmc_array(all_loop - 1, walker) + (Rp_Rs_ratio_stand_dev_MCMC*scale_factor*normal6(all_loop, walker))
            ELSE IF ((Rp_Rs_ratio_fixed == 'N') .AND. (Rp_Rs_ratio_norm_prior == 'N')) THEN
                Rp_Rs_ratio = Rp_Rs_ratio_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Rp_Rs_ratio(Rp_Rs_ratio_draw_loop, walker))
                IF ((Rp_Rs_ratio < Rp_Rs_ratio_begin) .OR. (Rp_Rs_ratio > Rp_Rs_ratio_end)) THEN
                    DO WHILE ((Rp_Rs_ratio < Rp_Rs_ratio_begin) .OR. (Rp_Rs_ratio > Rp_Rs_ratio_end))
                        Rp_Rs_ratio_draw_loop = Rp_Rs_ratio_draw_loop + 1
                        Rp_Rs_ratio = Rp_Rs_ratio_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Rp_Rs_ratio(Rp_Rs_ratio_draw_loop, walker))
                    END DO
                END IF
            END IF
            Rp_Rs_ratio_all_mcmc_array(all_loop, walker) = Rp_Rs_ratio
            
            


            IF (Rs_solar_fixed == 'Y') THEN
                Rs_solar = Rs_solar_prior
            ELSE IF ((Rs_solar_fixed == 'N') .AND. (Rs_solar_norm_prior == 'Y')) THEN
                Rs_solar = Rs_solar_accept_mcmc_array(all_loop - 1, walker) + (Rs_solar_stand_dev_MCMC*scale_factor*normal7(all_loop, walker))
            ELSE IF ((Rs_solar_fixed == 'N') .AND. (Rs_solar_norm_prior == 'N')) THEN
                Rs_solar = Rs_solar_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Rs(Rs_solar_draw_loop, walker))
                IF ((Rs_solar < Rs_solar_begin) .OR. (Rs_solar > Rs_solar_end)) THEN
                    DO WHILE ((Rs_solar < Rs_solar_begin) .OR. (Rs_solar > Rs_solar_end))
                        Rs_solar_draw_loop = Rs_solar_draw_loop + 1
                        Rs_solar = Rs_solar_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Rs(Rs_solar_draw_loop, walker))
                    END DO
                END IF
            END IF
            Rs_solar_all_mcmc_array(all_loop, walker) = Rs_solar
            
            


            IF (Ms_solar_fixed == 'Y') THEN
                Ms_solar = Ms_solar_prior
            ELSE IF ((Ms_solar_fixed == 'N') .AND. (Ms_solar_norm_prior == 'Y')) THEN
                Ms_solar = Ms_solar_accept_mcmc_array(all_loop - 1, walker) + (Ms_solar_stand_dev_MCMC*scale_factor*normal8(all_loop, walker))
            ELSE IF ((Ms_solar_fixed == 'N') .AND. (Ms_solar_norm_prior == 'N')) THEN
                Ms_solar = Ms_solar_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Ms(Ms_solar_draw_loop, walker))
                IF ((Ms_solar < Ms_solar_begin) .OR. (Ms_solar > Ms_solar_end)) THEN
                    DO WHILE ((Ms_solar < Ms_solar_begin) .OR. (Ms_solar > Ms_solar_end))
                        Ms_solar_draw_loop = Ms_solar_draw_loop + 1
                        Ms_solar = Ms_solar_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Ms(Ms_solar_draw_loop, walker))
                    END DO
                END IF
            END IF
            Ms_solar_all_mcmc_array(all_loop, walker) = Ms_solar
            
            

      
            IF (Ecc_fixed == 'Y') THEN
                Ecc = Ecc_prior
            ELSE IF ((Ecc_fixed == 'N') .AND. (Ecc_norm_prior == 'Y')) THEN
                Ecc = Ecc_accept_mcmc_array(all_loop - 1, walker) + (Ecc_stand_dev_MCMC*scale_factor*normal9(Ecc_draw_loop, walker))
                !If a negative Ecc is drawn, keep drawing an Ecc value until it is positive. Cannot have a negative Ecc.
                IF (Ecc < 0) THEN
                    DO WHILE (Ecc < 0.0D0)
                        Ecc_draw_loop = Ecc_draw_loop + 1
                        Ecc = Ecc_accept_mcmc_array(all_loop - 1, walker) + (Ecc_stand_dev_MCMC*scale_factor*normal9(Ecc_draw_loop, walker))
                    END DO
                END IF
            ELSE IF ((Ecc_fixed == 'N') .AND. (Ecc_norm_prior == 'N')) THEN
                Ecc = Ecc_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Ecc(Ecc_draw_loop, walker))
                IF ((Ecc < Ecc_begin) .OR. (Ecc > Ecc_end)) THEN
                    DO WHILE ((Ecc < Ecc_begin) .OR. (Ecc > Ecc_end))
                        Ecc_draw_loop = Ecc_draw_loop + 1
                        Ecc = Ecc_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Ecc(Ecc_draw_loop, walker))
                    END DO
                END IF
            END IF
            Ecc_all_mcmc_array(all_loop, walker) = Ecc
            
            

            
            IF (Inc_fixed == 'Y') THEN
                Inc = Inc_prior
            ELSE IF ((Inc_fixed == 'N') .AND. (Inc_norm_prior == 'Y')) THEN
                Inc = Inc_accept_mcmc_array(all_loop - 1, walker) + (Inc_stand_dev_MCMC*scale_factor*normal10(all_loop, walker))
                IF (Inc > 90.0D0) THEN
                    Inc = 90.0D0 - abs(Inc - 90.0D0)
                END IF
            ELSE IF ((Inc_fixed == 'N') .AND. (Inc_norm_prior == 'N')) THEN
                Inc = Inc_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Inc(Inc_draw_loop, walker))
                IF ((Inc < Inc_begin) .OR. (Inc > Inc_end)) THEN
                    DO WHILE ((Inc < Inc_begin) .OR. (Inc > Inc_end))
                        Inc_draw_loop = Inc_draw_loop + 1
                        Inc = Inc_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_Inc(Inc_draw_loop, walker))
                    END DO
                END IF
            END IF
            Inc_all_mcmc_array(all_loop, walker) = Inc
            
            

      
            IF (omega_arg_periastron_fixed == 'Y') THEN
                omega_arg_periastron = omega_arg_periastron_prior
            ELSE IF ((omega_arg_periastron_fixed == 'N') .AND. (omega_arg_periastron_norm_prior == 'Y')) THEN
                omega_arg_periastron = omega_arg_periastron_accept_mcmc_array(all_loop - 1, walker) + (omega_arg_periastron_stand_dev_MCMC*scale_factor*normal11(all_loop, walker))
            ELSE IF ((omega_arg_periastron_fixed == 'N') .AND. (omega_arg_periastron_norm_prior == 'N')) THEN
                omega_arg_periastron = omega_arg_periastron_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_omega(omega_arg_periastron_draw_loop, walker))
                IF ((omega_arg_periastron < omega_arg_periastron_begin) .OR. (omega_arg_periastron > omega_arg_periastron_end)) THEN
                    DO WHILE ((omega_arg_periastron < omega_arg_periastron_begin) .OR. (omega_arg_periastron > omega_arg_periastron_end))
                        omega_arg_periastron_draw_loop = omega_arg_periastron_draw_loop + 1
                        omega_arg_periastron = omega_arg_periastron_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_omega(omega_arg_periastron_draw_loop, walker))
                    END DO
                END IF
            END IF
            omega_arg_periastron_all_mcmc_array(all_loop, walker) = omega_arg_periastron
      
      
      
      
            IF (RV_zero_offset_fixed == 'Y') THEN
                RV_zero_offset = RV_zero_offset_prior
            ELSE IF ((RV_zero_offset_fixed == 'N') .AND. (RV_zero_offset_norm_prior == 'Y')) THEN
                RV_zero_offset = RV_zero_offset_accept_mcmc_array(all_loop - 1, walker) + (RV_zero_offset_stand_dev_MCMC*scale_factor*normal12(all_loop, walker))
            ELSE IF ((RV_zero_offset_fixed == 'N') .AND. (RV_zero_offset_norm_prior == 'Y')) THEN
                RV_zero_offset = RV_zero_offset_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_RV_zero(RV_zero_offset_draw_loop, walker))
                IF ((RV_zero_offset < RV_zero_offset_begin) .OR. (RV_zero_offset > RV_zero_offset_end)) THEN
                    DO WHILE ((RV_zero_offset < RV_zero_offset_begin) .OR. (RV_zero_offset > RV_zero_offset_end))
                        RV_zero_offset_draw_loop = RV_zero_offset_draw_loop + 1
                        RV_zero_offset = RV_zero_offset_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_RV_zero(RV_zero_offset_draw_loop, walker))
                    END DO
                END IF
            END IF
            RV_zero_offset_all_mcmc_array(all_loop, walker) = RV_zero_offset
      
      
      
      
            !Set the vsini_variance and stellar_rotation_angle_variance equal to the proposal variance of the previously accepted proposal.
            !Use the previously accepted proposal values of vsini and stellar_rotation_angle as the input to determine the proposal values for this iteration.
            IF (vsini_fixed == 'Y') THEN
                vsini = vsini_prior
            ELSE IF ((vsini_fixed == 'N') .AND. (vsini_norm_prior == 'Y')) THEN
                vsini = vsini_accept_mcmc_array(all_loop - 1, walker) + (vsini_stand_dev_MCMC*scale_factor*normal13(vsini_draw_loop, walker))      
                !If a negative vsini is drawn, keep drawing a vsini value until it is positive. Cannot have a negative vsini.
                IF (vsini < 0) THEN
                    DO WHILE (vsini < 0.0D0)
                        vsini_draw_loop = vsini_draw_loop + 1
                        vsini = vsini_accept_mcmc_array(all_loop - 1, walker) + (vsini_stand_dev_MCMC*scale_factor*normal13(vsini_draw_loop, walker))
                    END DO
                END IF
            ELSE IF ((vsini_fixed == 'N') .AND. (vsini_norm_prior == 'Y')) THEN
                vsini = vsini_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_vsini(vsini_draw_loop, walker))
                IF ((vsini < vsini_begin) .OR. (vsini > vsini_end)) THEN
                    DO WHILE ((vsini < vsini_begin) .OR. (vsini > vsini_end))
                        vsini_draw_loop = vsini_draw_loop + 1
                        vsini = vsini_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_vsini(vsini_draw_loop, walker))
                    END DO
                END IF
            END IF
            vsini_all_mcmc_array(all_loop, walker) = vsini    !Put the vsini value in the array holding all proposal values. Wait on putting the vsini
                                                      !value in the accepted proposal array until the likihood has been calculated to determine
                                                      !if the proposal is accepted or not.
            
            
            
            
            IF (stellar_rotation_angle_fixed == 'Y') THEN
                stellar_rotation_angle = stellar_rotation_angle_prior
            ELSE IF ((stellar_rotation_angle_fixed == 'N') .AND. (stellar_rotation_angle_norm_prior == 'Y')) THEN
                stellar_rotation_angle = stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) + (stellar_rotation_angle_stand_dev_MCMC*scale_factor*normal14(all_loop, walker))
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF
            ELSE IF ((stellar_rotation_angle_fixed == 'N') .AND. (stellar_rotation_angle_norm_prior == 'N')) THEN
                stellar_rotation_angle = stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_spin_orbit(stellar_rotation_angle_draw_loop, walker))
                IF ((stellar_rotation_angle < stellar_rotation_angle_begin) .OR. (stellar_rotation_angle > stellar_rotation_angle_end)) THEN
                    DO WHILE ((stellar_rotation_angle < stellar_rotation_angle_begin) .OR. (stellar_rotation_angle > stellar_rotation_angle_end))
                        stellar_rotation_angle_draw_loop = stellar_rotation_angle_draw_loop + 1
                        stellar_rotation_angle = stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) + (scale_factor*spread_rand_spin_orbit(stellar_rotation_angle_draw_loop, walker))
                    END DO
                END IF
                IF (stellar_rotation_angle <= -180.0D0) THEN
                    stellar_rotation_angle = 360.0D0 - stellar_rotation_angle
                END IF
                IF (stellar_rotation_angle > 180.0D0) THEN
                    stellar_rotation_angle = stellar_rotation_angle - 360.0D0
                END IF 
            END IF
            stellar_rotation_angle_all_mcmc_array(all_loop, walker) = stellar_rotation_angle      !Put the stellar_rotation_angle value in the array holding all proposal values. 
                                                                                                  !Wait on putting the stellar_rotation_angle value in the accepted proposal array 
                                                                                                  !until the likihood has been calculated to determine if the proposal is accepted 
                                                                                                  !or not.
      
              
              
            
        END IF

        count_tested_proposals = count_tested_proposals + 1

        !Now apply the RV offset to my data.
        Allocate(Data_my_rv_offset(num_my_rv,3), data(datafilelength,3))

        DO l = 1, num_my_rv
            Data_my_rv_offset(l,1) = Data_my_rv(l,1)
            Data_my_rv_offset(l,2) = Data_my_rv(l,2) + RV_offset_datasets
            Data_my_rv_offset(l,3) = Data_my_rv(l,3)
        END DO
   
        IF (other_RV_files == 'Y') THEN
            DO l = 1, datafilelength
                IF (l <= num_other_rv) THEN
                    Data(l,1) = Data_other_rv(l,1)
                    Data(l,2) = Data_other_rv(l,2)
                    Data(l,3) = Data_other_rv(l,3)
                ELSE
                    zz = l - num_other_rv
                    Data(l,1) = Data_my_rv_offset(zz,1)
                    Data(l,2) = Data_my_rv_offset(zz,2)
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
   
        Allocate(Data_adjusted(datafilelength,3), sorted_data_array(datafilelength,3))
        !Adjust data time relative to mid transit time based on the new orbital period and JD time mid transit.
        DO l = 1, datafilelength
            IF (((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1) - (FLOOR((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1)) < 0.5) THEN

                phase_time = ((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1) - &
                              (FLOOR((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1)) + 1.0D0

            ELSE
                phase_time = ((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1) - &
                              (FLOOR((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1))
            END IF

            Data_adjusted(l,1) = (Orbital_period_1*phase_time) - Orbital_period_1
            Data_adjusted(l,2) = Data(l,2)
            Data_adjusted(l,3) = Data(l,3)
        END DO

   
   
   
        !Verify that data array is ascending order as a function of time
        CALL Selection_sort(Data_adjusted(:,1),index_array)
  
        DO l = 1, datafilelength
            sorted_data_array(l,1) = Data_adjusted(index_array(l), 1)
            sorted_data_array(l,2) = Data_adjusted(index_array(l), 2)
            sorted_data_array(l,3) = Data_adjusted(index_array(l), 3)
            IF (index_array(l) /= l) THEN
                array_resorted = 'Y'
            ELSE
                array_resorted = 'N'
            END IF
        END DO
   
   
   
   
        IF (Jupiter_Earth_units == 'Y') THEN
            Mp = Mpp*Mj
        END IF

        IF (Jupiter_Earth_units == 'N') THEN
            Mp = Mpp*Me
        END IF
        
        Ms = Ms_solar*Mss
        Rs = Rs_solar*Rss
        
        
        Rorb = ((Orbital_period**2.0d0*G*(Ms + Mp))/(4.0d0*pi**2.0d0))**(1.0d0/3.0d0)     !Semi-major axis.
        Rorb_star = ((Orbital_period**2.0D0*G*((Mp**(3.0D0))/(Ms + Mp)**(2.0D0)))/(4.0d0*pi**2.0D0))**(1.0d0/3.0d0)     !Semi-major of stellar axis.

        IF (use_Rp_Rs_ratio == 'Y') THEN
             Rp = Rp_Rs_ratio * (Rs_solar * Rss)
        ELSE
            IF (Jupiter_Earth_units == 'Y') THEN
                Rp = Rpp*Rj
            END IF

            IF (Jupiter_Earth_units == 'N') THEN
                Rp = Rpp*RE
            END IF
        END IF
        
        Aplan = pi*Rp**2.0D0                         !Surface area of the planet.

        Rp2 = Rp**2.0D0                              !Square the radius of planet to speed up calculations.
        
        Rs2 = Rs**2.0D0                              !Square the radius of star to speed up calculations.
        
        IF (Ecc == 0) THEN 
            !Maximum amplitude caused by the exoplanet in a circular orbit.
            RVamplitude = SQRT((G*(Mp**(3.0D0))*(sin(Inc*(pi/180.0D0)))**(3.0D0))/(Rorb_star*sin(Inc*(pi/180.0D0))*(Ms + Mp)**(2.0D0)))
        ELSE
            !Maximum amplitude caused by the exoplanet in an eccentric orbit.
            RVamplitude = SQRT((G*(Mp**(3.0D0))*(sin(Inc*(pi/180.0D0)))**(3.0D0))/((1.0D0 - Ecc**2.0D0)*Rorb_star*sin(Inc*(pi/180.0D0))*(Ms + Mp)**(2.0D0)))
        END IF
        
        True_anomaly_start = pi - (omega_arg_periastron*(pi/180.0D0))

        IF (True_anomaly_start >= pi) THEN
            ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) + 2.0D0*pi 
        ELSE IF (True_anomaly_start <= -pi) THEN
            ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) - 2.0D0*pi
        ELSE
            ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc)))
        END IF
        
        !Calculate the amount of time (in seconds) after the passage of periastron has occurred when the planet is at mid transit.
        !This is determined from the transit mid time and the argument of periastron.
        IF (omega_arg_periastron > 180.0D0) THEN
            Mean_anomaly_transit = (2.0D0*pi) - (omega_arg_periastron*(pi/180.0D0)) + pi
        ELSE
            Mean_anomaly_transit = pi - (omega_arg_periastron*(pi/180.0D0))
        END IF
        
        JD_time_peri = (JD_time_mid_transit*day_sec) - ((Mean_anomaly_transit*Orbital_period)/(2.0D0*pi))
        time_peri_passage = (JD_time_mid_transit*day_sec) - JD_time_peri
        True_anomaly_transit = ((3.0D0*pi)/2.0D0) - (omega_arg_periastron*(pi/180.0D0))
        
        IF (True_anomaly_transit >= pi) THEN
            ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) + 2.0D0*pi 
        ELSE IF (True_anomaly_transit <= -pi) THEN
            ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) - 2.0D0*pi
        ELSE
            ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc)))
        END IF

        IF (Ecc == 0) THEN
            Time_transit = ((ecc_anomaly_transit*Orbital_period)/(2.0D0*pi)) + time_peri_passage
        ELSE
            Time_transit = (((ecc_anomaly_transit - (Ecc*sin(ecc_anomaly_transit)))*Orbital_period)/(2.0D0*pi)) + time_peri_passage
        END IF
        
        !Create a new time array for the data which has been created for the JD time of mid transit and ecc anomalies.
        Allocate(time_data_array(datafilelength))
        DO l = 1, datafilelength
            time_data_array(l) = (sorted_data_array(l,1)*day_sec) + Time_transit
        END DO

        !Now apply the change in the RV zero offset from the prior offset to the data.
        Allocate(RV_offset_data_array(datafilelength,3))
        DO l = 1, datafilelength
            RV_offset_data_array(l,1) = sorted_data_array(l,1)
            RV_offset_data_array(l,2) = sorted_data_array(l,2) + (RV_zero_offset - RV_zero_offset_prior)
            RV_offset_data_array(l,3) = sorted_data_array(l,3)
        END DO
   
        Allocate(RV_theory(total_interval, 2), data_points_compare(total_interval))

        num_compare = 0
        Planet_star_distance = 0                       !Set the distance between center of the star to the center of the planet (orbital radius) to zero.
        
        !******************************************************************************************************************************************* 

        DO Time_loop = 1, total_interval

            !The time will come from the first data point which must be normalized as a fraction of a day before mid transit.
            time = time_data_array(Time_loop)

            !Set transit flag to N when outside of transit event.
            transit_flag = 'N'
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
                    IF (order > 1) THEN
                        ecc_anomaly_before = sum_ecc_anomaly
                    END IF
                    sum_ecc_anomaly = sum_ecc_anomaly + ((2.0D0/order)*Bessel_value*sin(order*(((2.0D0*pi) &
                                      /Orbital_period)*(time - time_peri_passage))))
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
                Planet_star_distance = Rorb + Rorb_star
            ELSE
                !The time of the similation for a specific true_anomaly.
                Time_check = (((ecc_anomaly - (Ecc*sin(ecc_anomaly)))*Orbital_period)/(2.0D0*pi)) + time_peri_passage         
                !The distance between the center of the planet to the center of the star in an eccentric orbit.
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
            
            IF ((Xpos <= 0) .AND. (Zpos >= 0)) THEN
                !Planet is currently in quadrant three so add pi.
                IF (Zpos == 0) THEN
                    phase_angle = pi/2.0D0
                    Phase_angle_observed = pi/2.0D0
                ELSE
                    phase_angle = atan(Xpos/Zpos) + pi
                    !Taking into account orbital inclination.
                    Phase_angle_observed = atan(-(sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos) + pi
                END IF
            ELSE IF ((Xpos >= 0) .AND. (Zpos >= 0)) THEN
                !Planet is currently in quadrant four so add pi.
                IF (Zpos == 0) THEN
                    phase_angle = pi/2.0D0 + pi
                    Phase_angle_observed = pi/2.0D0 + pi
                ELSE
                    phase_angle = atan(Xpos/Zpos) + pi
                    !Taking into account orbital inclination.
                    Phase_angle_observed = atan((sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos) + pi
                END IF
            ELSE IF ((Xpos >= 0) .AND. (Zpos <= 0)) THEN
                !Planet is currently in quadrant one so add 2pi.
                IF (Zpos == 0) THEN
                    phase_angle = pi/2.0D0
                    Phase_angle_observed = pi/2.0D0
                ELSE
                    phase_angle = 2.0D0*pi + atan(Xpos/Zpos)
                    !Taking into account orbital inclination.
                    Phase_angle_observed = 2.0D0*pi + atan((sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos)
                END IF
            ELSE IF ((Xpos <= 0) .AND. (Zpos <= 0)) THEN
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
                
                IF ((Rp2/Rs2) >= 0.030) THEN
                    Radius_planet_array = M_pixels / 2.0D0       !Radius of planet in the pixel array. 
                    Center_of_planet_x = M_pixels / 2.0D0        !The center of the planet on the X-axis.
                    Center_of_planet_y = M_pixels / 2.0D0        !The center of the planet on the Y-axis.
                    Pixel = Rp/Radius_planet_array               !The number of meters per pixel.
                    Area_pixel = Pixel**2.0D0                    !The area of each pixel.
                    Io_Pixel = Io * Area_pixel                   !Variable to speed up calculations (also physically represents the
                                                                 !luminosity of the brightest pixel).
                                                                 
                    DO i = 1, M_pixels
                        X_pixel = (pixel * i) + (Xpos - Rp)           !Calculates the location of the pixel on the x-axis.
                        X_pixel2 = X_pixel**2.0D0
                        XpXp = X_pixel * Xpos                         !temporary var for speed calculation
                        
                        DO j = 1, M_pixels
                            Y_pixel = (pixel * j) + abs(Ypos) - Rp        !Calculates the location of the pixel on the y-axis.
                            !Calculates the location of the pixel along the x-axis of the rotation axis of the star.
                            X_pixel_prime = (X_pixel*cos((stellar_rotation_angle*pi)/180.0D0)) + &
                                            (Y_pixel*sin((stellar_rotation_angle*pi)/180.0D0))     
                            Dist_center_pixel = X_pixel2 + Y_pixel**2.0D0         !squared distance of pixel from star
                            !Calculates the values of the pixels according to how far away they are from the center of the star and the plane.
                            !squared distance of pixel from planet using a limb darkening equation.
                            Dist_planet_pixel = Dist_center_pixel  - (2.0D0*(XpXp + (Y_pixel*Ypos))) + Dist2
                            Sub_planet_velocity = vsini*(X_pixel_prime/Rs)
                            
                            IF ((Dist_center_pixel <= Rs2) .AND. (Dist_planet_pixel <= Rp2)) THEN
                                IF (linear_quadratic == 'q') THEN
                                    Lblocked2 = Io_Pixel*(1.0D0-q_1*(1.0D0-sqrt(abs(1.0D0-(Dist_center_pixel/Rs2)))) - &
                                                q_2*(1.0D0 - sqrt(abs(1.0D0-(Dist_center_pixel/Rs2))))**2.0D0)          !Quadratic limb darkening equation.
                                ELSE
                                    Lblocked2 = Io_Pixel*(1.0D0-u*(1.0D0-sqrt(abs(1.0D0-(Dist_center_pixel/Rs2)))))     !First order limb darkening equation.          
                                END IF

                                Lblocked = Lblocked + Lblocked2                                         
                                v_rm = v_rm - ((Lblocked2*Sub_planet_velocity)*((((2.0D0*vturb**2.0D0)+(2.0D0*vsini**2.0D0))/((2.0D0*vturb**2.0D0) &
                                       + vsini**2.0D0))**(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity**2.0D0)/((2.0D0*vturb**2.0D0)+(vsini**2.0D0)))))   
                                !Anomalous velocity of each pixel.
                            END IF
                        END DO
                    END DO
                    
                    Total_L = 1.0D0 - Lblocked
                    Total_RM = 0.0D0 + v_rm                       !Total anomalous velocity for all the pixels.
                    
                ELSE IF ((Rp2/Rs2) <= 0.030) THEN
  
                    !Calculates the location of the center of the planet along the x-axis of the rotation axis of the star.
                    X_prime = (Xpos*cos((stellar_rotation_angle*pi)/180.0D0)) &
                              + (Ypos*sin((stellar_rotation_angle*pi)/180.0D0))     
                    set_distance_center = Distance_center         !The limb darkening equation will use this distance as long as the center of the 
                                                                  !planet is inside the radius of the star.
                    Sub_planet_velocity = vsini*(X_prime/Rs)      !Calculate the subplanetary velocity (the stellar velocity blocked by the 
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
                                  
                        !The limb darkening equation will use this distance if any part of the disc of the planet is outside the radius of the star.
                        !This is the distance between the center of the star to the center of the area inside the star blocked by the planet.           
                        set_distance_center = ((abs(Distance_center) - Rp) + Rs)/2.0D0
                        !Next the program calculates the length between the two intersection point.
                        Length = sqrt((X_int_1 - X_int_2)**2.0D0 + (Y_int_1 - Y_int_2)**2.0D0)
                        !Calculate the angle between the Y position of the center of planet and the x position of the center of planet.
                        !This is used to determine the center of the area inside the star blocked by the planet in the stellar rotational
                        !axis coordinate system.
                        Theta_inside = atan(Ypos/Xpos)
                        !Calculates the distance from the center of the star (along x axis) to the center of the area inside the star blocked by the planet.
                        !the pi factor added in this equation guarantees that Xpos_in has the correct sign.
                        Xpos_in = set_distance_center*cos(Theta_inside + pi)
                        Xpos_inside = Xpos_in
                        
                        !This makes sure Xpos_inside has the correct sign.
                        IF (Xpos >= 0) THEN
                            Xpos_inside = -Xpos_in
                        END IF
               
                        !Calculates the distance from the center of the star (along y axis) to the center of the area inside the star blocked by the planet.
                        Ypos_inside = abs(set_distance_center*sin(Theta_inside))
                        !Changes the x-coordinate to the stellar rotation axis by an angle formed between the orbital plane of the planet and the 
                        !stellar roatation plane of the star.
                        x_prime_distance = (Xpos_inside*cos((stellar_rotation_angle*pi)/180.0D0)) &
                                           + (Ypos_inside*sin((stellar_rotation_angle*pi)/180.0D0))
                        Sub_planet_velocity = vsini*(x_prime_distance/Rs)    !Calculate the subplanetary velocity (the stellar velocity blocked by 
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
                    v_rm = - ((Lblocked*Sub_planet_velocity)*((((2.0D0*vturb**2.0D0)+(2.0D0*vsini**2.0D0))/((2.0D0*vturb**2.0D0) &
                           + vsini**2.0D0))**(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity**2.0D0)/((2.0D0*vturb**2.0D0)+(vsini**2.0D0)))))
                    Total_L = 1.0D0 - Lblocked                                             !Total amount of light blocked by the planet.
                    Total_RM = 0.0D0 + v_rm                                                !Total anomalous velocity.
                END IF

            END IF
            
            !If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is negative then the secondary transit 
            !(occulation) begins.
            IF ((Distance_center <= (Rs + Rp)) .AND. (Zpos < 0)) THEN           

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
            IF ((transit_flag == 'Y') .OR. (time_compare_flag == 'Y') .OR. (use_out_transit_rv_for_fit == 'Y')) THEN
                num_compare = num_compare + 1
                data_points_compare(Time_loop) = 'Y'
            END IF               

            RV_theory(Time_loop,1) = Time_ref
            RV_theory(Time_loop,2) = RV

        END DO
         
        !*******************************************************************************************************************************
   
        Number_points1 = 1
        chi_2 = 0.0D0          !Set chi squared equal to zero.
        model_data_difference = 0.0D0        

        DO s = 1, total_interval
            IF (data_points_compare(s) == 'Y') THEN
                !Now do chi squared analysis but only for data points during transit.
                !Now calculate chi squared from the adjusted RV model.

                chi_2 = chi_2 + ((RV_offset_data_array(s,2) - RV_Theory(s,2))/RV_offset_data_array(s,3))**2.0D0    
                 
                model_data_difference = model_data_difference + ABS(RV_offset_data_array(s,2) - RV_Theory(s,2))

                Number_points1 = Number_points1 + 1              !Counter to determine the number of elements used to calculate sigma squared.

            END IF
        END DO
        
        IF (stellar_rotation_angle_1sigerr <= 0.0D0) THEN
            stellar_rotation_angle_1sigerr_temp = 1.0D0
        ELSE
            stellar_rotation_angle_1sigerr_temp = stellar_rotation_angle_1sigerr
        END IF
        
        IF (vsini_1sigerr <= 0.0D0) THEN
            vsini_1sigerr_temp = 1.0D0
        ELSE
            vsini_1sigerr_temp = vsini_1sigerr
        END IF
        
        IF (RV_offset_datasets_1sigerr <= 0.0D0) THEN
            RV_offset_datasets_1sigerr_temp = 1.0D0
        ELSE
            RV_offset_datasets_1sigerr_temp = RV_offset_datasets_1sigerr
        END IF
        
        IF (RV_zero_offset_1sigerr <= 0.0D0) THEN
            RV_zero_offset_1sigerr_temp = 1.0D0
        ELSE
            RV_zero_offset_1sigerr_temp = RV_zero_offset_1sigerr
        END IF
        
        IF (Orbital_period_1sigerr <= 0.0D0) THEN
            Orbital_period_1sigerr_temp = 1.0D0
        ELSE
            Orbital_period_1sigerr_temp = Orbital_period_1sigerr
        END IF
        
        IF (JD_time_mid_transit_1sigerr <= 0.0D0) THEN
            JD_time_mid_transit_1sigerr_temp = 1.0D0
        ELSE
            JD_time_mid_transit_1sigerr_temp = JD_time_mid_transit_1sigerr
        END IF
        
        IF (Mp_1sigerr <= 0.0D0) THEN
            Mp_1sigerr_temp = 1.0D0
        ELSE
            Mp_1sigerr_temp = Mp_1sigerr
        END IF
        
        IF (Rp_1sigerr <= 0.0D0) THEN
            Rp_1sigerr_temp = 1.0D0
        ELSE
            Rp_1sigerr_temp = Rp_1sigerr
        END IF
        
        IF (Ecc_1sigerr <= 0.0D0) THEN
            Ecc_1sigerr_temp = 1.0D0
        ELSE
            Ecc_1sigerr_temp = Ecc_1sigerr
        END IF
        
        IF (Inc_1sigerr <= 0.0D0) THEN
            Inc_1sigerr_temp = 1.0D0
        ELSE
            Inc_1sigerr_temp = Inc_1sigerr
        END IF
        
        IF (Ms_solar_1sigerr <= 0.0D0) THEN
            Ms_solar_1sigerr_temp = 1.0D0
        ELSE
            Ms_solar_1sigerr_temp = Ms_solar_1sigerr
        END IF
        
        IF (Rs_solar_1sigerr <= 0.0D0) THEN
            Rs_solar_1sigerr_temp = 1.0D0
        ELSE
            Rs_solar_1sigerr_temp = Rs_solar_1sigerr
        END IF
        
        IF (Rp_Rs_ratio_1sigerr <= 0.0D0) THEN
            Rp_Rs_ratio_1sigerr_temp = 1.0D0
        ELSE
            Rp_Rs_ratio_1sigerr_temp = Rp_Rs_ratio_1sigerr
        END IF
        
        IF ((impose_prior_vsini == 'Y') .AND. (impose_prior_stellar_rotation_angle == 'Y')) THEN
                
            chi_2 = chi_2 + ((stellar_rotation_angle - stellar_rotation_angle_prior)/stellar_rotation_angle_1sigerr_temp)**2.0D0 &
                    + ((vsini - vsini_prior)/vsini_1sigerr_temp)**2.0D0 &
                    + ((RV_offset_datasets - RV_offset_datasets_prior)/RV_offset_datasets_1sigerr_temp)**2.0D0 &
                    + ((RV_zero_offset - RV_zero_offset_prior)/RV_zero_offset_1sigerr_temp)**2.0D0 &
                    + ((Orbital_period_1 - Orbital_period_prior)/Orbital_period_1sigerr_temp)**2.0D0 &
                    + ((JD_time_mid_transit - JD_time_mid_transit_prior)/JD_time_mid_transit_1sigerr_temp)**2.0D0 &
                    + ((Mpp - Mp_prior)/Mp_1sigerr_temp)**2.0D0 + ((Rpp - Rp_prior)/Rp_1sigerr_temp)**2.0D0 &
                    + ((Ecc - Ecc_prior)/Ecc_1sigerr_temp)**2.0D0 + ((Inc - Inc_prior)/Inc_1sigerr_temp)**2.0D0 &
                    + ((Ms_solar - Ms_solar_prior)/Ms_solar_1sigerr_temp)**2.0D0 + ((Rs_solar - Rs_solar_prior)/Rs_solar_1sigerr_temp)**2.0D0 &
                    + ((Rp_Rs_ratio - Rp_Rs_ratio_prior)/Rp_Rs_ratio_1sigerr_temp)**2.0D0
                        
        ELSE IF ((impose_prior_vsini == 'Y') .AND. (impose_prior_stellar_rotation_angle == 'N')) THEN
                
            chi_2 = chi_2 + ((vsini - vsini_prior)/vsini_1sigerr_temp)**2.0D0 & 
                    + ((RV_offset_datasets - RV_offset_datasets_prior)/RV_offset_datasets_1sigerr_temp)**2.0D0 &
                    + ((Orbital_period_1 - Orbital_period_prior)/Orbital_period_1sigerr_temp)**2.0D0 &
                    + ((RV_zero_offset - RV_zero_offset_prior)/RV_zero_offset_1sigerr_temp)**2.0D0 &
                    + ((JD_time_mid_transit - JD_time_mid_transit_prior)/JD_time_mid_transit_1sigerr_temp)**2.0D0 &
                    + ((Mpp - Mp_prior)/Mp_1sigerr_temp)**2.0D0 + ((Rpp - Rp_prior)/Rp_1sigerr_temp)**2.0D0 &
                    + ((Ecc - Ecc_prior)/Ecc_1sigerr_temp)**2.0D0 + ((Inc - Inc_prior)/Inc_1sigerr_temp)**2.0D0 &
                    + ((Ms_solar - Ms_solar_prior)/Ms_solar_1sigerr_temp)**2.0D0 + ((Rs_solar - Rs_solar_prior)/Rs_solar_1sigerr_temp)**2.0D0 &
                    + ((Rp_Rs_ratio - Rp_Rs_ratio_prior)/Rp_Rs_ratio_1sigerr_temp)**2.0D0
                        
        ELSE IF ((impose_prior_vsini == 'N') .AND. (impose_prior_stellar_rotation_angle == 'Y')) THEN
                
            chi_2 = chi_2 + ((stellar_rotation_angle - stellar_rotation_angle_prior)/stellar_rotation_angle_1sigerr_temp)**2.0D0 &
                    + ((RV_offset_datasets - RV_offset_datasets_prior)/RV_offset_datasets_1sigerr_temp)**2.0D0 &
                    + ((RV_zero_offset - RV_zero_offset_prior)/RV_zero_offset_1sigerr_temp)**2.0D0 &
                    + ((Orbital_period_1 - Orbital_period_prior)/Orbital_period_1sigerr_temp)**2.0D0 &
                    + ((JD_time_mid_transit - JD_time_mid_transit_prior)/JD_time_mid_transit_1sigerr_temp)**2.0D0 &
                    + ((Mpp - Mp_prior)/Mp_1sigerr_temp)**2.0D0 + ((Rpp - Rp_prior)/Rp_1sigerr_temp)**2.0D0 &
                    + ((Ecc - Ecc_prior)/Ecc_1sigerr_temp)**2.0D0 + ((Inc - Inc_prior)/Inc_1sigerr_temp)**2.0D0 &
                    + ((Ms_solar - Ms_solar_prior)/Ms_solar_1sigerr_temp)**2.0D0 + ((Rs_solar - Rs_solar_prior)/Rs_solar_1sigerr_temp)**2.0D0 &
                    + ((Rp_Rs_ratio - Rp_Rs_ratio_prior)/Rp_Rs_ratio_1sigerr_temp)**2.0D0
                        
        ELSE
                
            chi_2 = chi_2 + ((RV_offset_datasets - RV_offset_datasets_prior)/RV_offset_datasets_1sigerr_temp)**2.0D0 &
                    + ((RV_zero_offset - RV_zero_offset_prior)/RV_zero_offset_1sigerr_temp)**2.0D0 &
                    + ((Orbital_period_1 - Orbital_period_prior)/Orbital_period_1sigerr_temp)**2.0D0 &
                    + ((JD_time_mid_transit - JD_time_mid_transit_prior)/JD_time_mid_transit_1sigerr_temp)**2.0D0 &
                    + ((Mpp - Mp_prior)/Mp_1sigerr_temp)**2.0D0 + ((Rpp - Rp_prior)/Rp_1sigerr_temp)**2.0D0 &
                    + ((Ecc - Ecc_prior)/Ecc_1sigerr_temp)**2.0D0 + ((Inc - Inc_prior)/Inc_1sigerr_temp)**2.0D0 &
                    + ((Ms_solar - Ms_solar_prior)/Ms_solar_1sigerr_temp)**2.0D0 + ((Rs_solar - Rs_solar_prior)/Rs_solar_1sigerr_temp)**2.0D0 &
                    + ((Rp_Rs_ratio - Rp_Rs_ratio_prior)/Rp_Rs_ratio_1sigerr_temp)**2.0D0
              
        END IF
        
        
        
        !The reduced chi squared.
        r_Chi_2 = chi_2/(Number_points1 - Number_fit - 1)
   
        !The single likelihood.
        likelihood = dexp(-chi_2/2.0D0)
        
        !Now check if this is the first iteration of MCMC. If it is the values are accepted by default.
        IF (all_loop == 1) THEN
            culumative_likelihood = likelihood
            culumative_all_likelihood = culumative_likelihood
      
            likelihood_array_MCMC(all_loop, walker) = likelihood
            culumative_likelihood_array_MCMC(all_loop, walker) = culumative_likelihood
            Chi_squared_array_MCMC(all_loop, walker) = chi_2
            Reduced_Chi_Squared_array_MCMC(all_loop, walker) = r_Chi_2
                  
            !Rejection flag for the proposal parameters is set to false.
            Reject_flag = 'F'
            stellar_rotation_angle_mean = stellar_rotation_angle
           
        ELSE IF (all_loop > 1) THEN
            !If the MCMC loop is greater than 1 (i.e., second iteration and beyond), we need to compare the likelihood
            !of the current proposals with the likelihood computed from the previous proposal.
            acceptance_param_previous = Chi_squared_array_MCMC(all_loop - 1, walker)
            acceptance_param_current = chi_2
            
            accept_prob = 1.0D0
            random_draw = 1.0D0
            
            IF (acceptance_param_current > acceptance_param_previous) THEN   
                !Calculate the probability of acceptance. Then draw a random number between 0 and 1.
                !If the acceptance probability is greater than or equal to the randomly drawn number
                !than the proposal is accepted. Otherwise, the proposal is rejected.
                accept_prob = dexp(-(acceptance_param_current - acceptance_param_previous)/2.0D0)
                random_draw = random_draw_array(all_loop, walker)
                !random_draw = ran1()
            END IF
            
            !If acceptance_param_current is less than or equal to acceptance_param_previous, the new proposal values are automatically accepted.
            !If acceptance_param_current is greater than acceptance_param_previous, the new proposal values are accepted
            !with probability exp(-(acceptance_param_current - acceptance_param_previous)/2).
            
            IF ((acceptance_param_current <= acceptance_param_previous) .OR. (accept_prob >= random_draw)) THEN
                
                culumative_likelihood = culumative_likelihood + likelihood/(all_loop + culumative_likelihood)
                !culumative_all_likelihood = culumative_likelihood
   
                likelihood_array_MCMC(all_loop, walker) = likelihood
                culumative_likelihood_array_MCMC(all_loop, walker) = culumative_likelihood
                Chi_squared_array_MCMC(all_loop, walker) = chi_2
                Reduced_Chi_Squared_array_MCMC(all_loop, walker) = r_Chi_2
            
                vsini_sum = vsini_sum + vsini                
                vsini_mean = vsini_sum/all_loop
                vsini_diff_mean_2 = vsini_diff_mean_2 + abs(vsini - vsini_mean)**2.0D0
                vsini_variance = vsini_diff_mean_2/(all_loop - 1.0D0)  !sample variance.
                vsini_stand_dev = SQRT(vsini_variance)
            
                IF ((stellar_rotation_angle_mean > 90) .AND. (stellar_rotation_angle < -90)) THEN
                    stellar_rotation_angle_sum = stellar_rotation_angle_sum + stellar_rotation_angle + 360.0D0
                ELSE IF ((stellar_rotation_angle_mean < -90) .AND. (stellar_rotation_angle > 90)) THEN
                    stellar_rotation_angle_sum = stellar_rotation_angle_sum + stellar_rotation_angle - 360.0D0
                ELSE  
                    stellar_rotation_angle_sum = stellar_rotation_angle_sum + stellar_rotation_angle
                END IF
                stellar_rotation_angle_mean = stellar_rotation_angle_sum/all_loop
         
                stellar_rotation_angle_diff_mean_2 = stellar_rotation_angle_diff_mean_2 + abs(stellar_rotation_angle - stellar_rotation_angle_mean)**2.0D0
         
                IF (stellar_rotation_angle_mean <= -180) THEN
                    stellar_rotation_angle_mean = 360.0D0 - stellar_rotation_angle_mean
                ELSE IF (stellar_rotation_angle_mean > 180) THEN
                    stellar_rotation_angle_mean = stellar_rotation_angle_mean - 360.0D0
                END IF
         
                stellar_rotation_angle_variance = stellar_rotation_angle_diff_mean_2/(all_loop - 1.0D0)  !sample variance.
                stellar_rotation_angle_stand_dev = SQRT(stellar_rotation_angle_variance)
            
                RV_offset_datasets_sum = RV_offset_datasets_sum + RV_offset_datasets
                RV_offset_datasets_mean = RV_offset_datasets_sum/all_loop
                RV_offset_datasets_diff_mean_2 = RV_offset_datasets_diff_mean_2 + abs(RV_offset_datasets - RV_offset_datasets_mean)**2.0D0
                RV_offset_datasets_variance = RV_offset_datasets_diff_mean_2/(all_loop - 1.0D0)  !sample variance.
                RV_offset_datasets_stand_dev = SQRT(RV_offset_datasets_variance)
            
                RV_zero_offset_sum = RV_zero_offset_sum + RV_zero_offset
                RV_zero_offset_mean = RV_zero_offset_sum/all_loop
                RV_zero_offset_diff_mean_2 = RV_zero_offset_diff_mean_2 + abs(RV_zero_offset - RV_zero_offset_mean)**2.0D0
                RV_zero_offset_variance = RV_zero_offset_diff_mean_2/(all_loop - 1.0D0) 
                RV_zero_offset_stand_dev = SQRT(RV_zero_offset_variance)
            
                Orbital_period_sum = Orbital_period_sum + Orbital_period_1
                Orbital_period_mean = Orbital_period_sum/all_loop
                Orbital_period_diff_mean_2 = Orbital_period_diff_mean_2 + abs(Orbital_period_1 - Orbital_period_mean)**2.0D0
                Orbital_period_variance = Orbital_period_diff_mean_2/(all_loop - 1.0D0)
                Orbital_period_stand_dev = SQRT(Orbital_period_variance)  

                JD_time_mid_transit_sum = JD_time_mid_transit_sum + JD_time_mid_transit
                JD_time_mid_transit_mean = JD_time_mid_transit_sum/all_loop
                JD_time_mid_transit_diff_mean_2 = JD_time_mid_transit_diff_mean_2 + abs(JD_time_mid_transit - JD_time_mid_transit_mean)**2.0D0
                JD_time_mid_transit_variance = JD_time_mid_transit_diff_mean_2/(all_loop - 1.0D0)    
                JD_time_mid_transit_stand_dev = SQRT(JD_time_mid_transit_variance)

                Mp_sum = Mp_sum + Mpp
                Mp_mean = Mp_sum/all_loop
                Mp_diff_mean_2 = Mp_diff_mean_2 + abs(Mpp - Mp_mean)**2.0D0
                Mp_variance = Mp_diff_mean_2/(all_loop - 1.0D0)
                Mp_stand_dev = SQRT(Mp_variance)

                Rp_sum = Rp_sum + Rpp
                Rp_mean = Rp_sum/all_loop
                Rp_diff_mean_2 = Rp_diff_mean_2 + abs(Rpp - Rp_mean)**2.0D0
                Rp_variance = Rp_diff_mean_2/(all_loop - 1.0D0)
                Rp_stand_dev = SQRT(Rp_variance)

                Ecc_sum = Ecc_sum + Ecc
                Ecc_mean = Ecc_sum/all_loop
                Ecc_diff_mean_2 = Ecc_diff_mean_2 + abs(Ecc - Ecc_mean)**2.0D0
                Ecc_variance = Ecc_diff_mean_2/(all_loop - 1.0D0)
                Ecc_stand_dev = SQRT(Ecc_variance)

                Inc_sum = Inc_sum + Inc
                Inc_mean = Inc_sum/all_loop
                Inc_diff_mean_2 = Inc_diff_mean_2 + abs(Inc - Inc_mean)**2.0D0
                Inc_variance = Inc_diff_mean_2/(all_loop - 1.0D0)
                Inc_stand_dev = SQRT(Inc_variance)

                omega_arg_periastron_sum = omega_arg_periastron_sum + omega_arg_periastron
                omega_arg_periastron_mean = omega_arg_periastron_sum/all_loop
                omega_arg_periastron_diff_mean_2 = omega_arg_periastron_diff_mean_2 + abs(omega_arg_periastron - omega_arg_periastron_mean)**2.0D0
                omega_arg_periastron_variance = omega_arg_periastron_diff_mean_2/(all_loop - 1.0D0)     
                omega_arg_periastron_stand_dev = SQRT(omega_arg_periastron_variance)    
            
                Ms_solar_sum = Ms_solar_sum + Ms_solar
                Ms_solar_mean = Ms_solar_sum/all_loop
                Ms_solar_diff_mean_2 = Ms_solar_diff_mean_2 + abs(Ms_solar - Ms_solar_mean)**2.0D0
                Ms_solar_variance = Ms_solar_diff_mean_2/(all_loop - 1.0D0)
                Ms_solar_stand_dev = SQRT(Ms_solar_variance)
            
                Rs_solar_sum = Rs_solar_sum + Rs_solar
                Rs_solar_mean = Rs_solar_sum/all_loop
                Rs_solar_diff_mean_2 = Rs_solar_diff_mean_2 + abs(Rs_solar - Rs_solar_mean)**2.0D0
                Rs_solar_variance = Rs_solar_diff_mean_2/(all_loop - 1.0D0)
                Rs_solar_stand_dev = SQRT(Rs_solar_variance)
            
                Rp_Rs_ratio_sum = Rp_Rs_ratio_sum + Rp_Rs_ratio
                Rp_Rs_ratio_mean = Rp_Rs_ratio_sum/all_loop
                Rp_Rs_ratio_diff_mean_2 = Rp_Rs_ratio_diff_mean_2 + abs(Rp_Rs_ratio - Rp_Rs_ratio_mean)**2.0D0
                Rp_Rs_ratio_variance = Rp_Rs_ratio_diff_mean_2/(all_loop - 1.0D0)
                Rp_Rs_ratio_stand_dev = SQRT(Rp_Rs_ratio_variance)
                
                vsini_accept_mcmc_array(all_loop, walker) = vsini
                vsini_variance_current_array(all_loop, walker) = vsini_variance             !The variance updated at every iteration.
                vsini_stand_dev_current_array(all_loop, walker) = vsini_stand_dev           !The standard deviation updated at every iteration.
         
                stellar_rotation_angle_accept_mcmc_array(all_loop, walker) = stellar_rotation_angle
                stellar_rotation_angle_variance_current_array(all_loop, walker) = stellar_rotation_angle_variance         !The variance updated at every iteration.
                stellar_rotation_angle_stand_dev_current_array(all_loop, walker) = stellar_rotation_angle_stand_dev       !The standard deviation updated at every iteration. 
         
                RV_offset_datasets_accept_mcmc_array(all_loop, walker) = RV_offset_datasets
                RV_zero_offset_accept_mcmc_array(all_loop, walker) = RV_zero_offset
                Orbital_period_accept_mcmc_array(all_loop, walker) = Orbital_period_1
                JD_time_mid_transit_accept_mcmc_array(all_loop, walker) = JD_time_mid_transit
                Mp_accept_mcmc_array(all_loop, walker) = Mpp
                Rp_accept_mcmc_array(all_loop, walker) = Rpp
                Ms_solar_accept_mcmc_array(all_loop, walker) = Ms_solar
                Rs_solar_accept_mcmc_array(all_loop, walker) = Rs_solar
                Rp_Rs_ratio_accept_mcmc_array(all_loop, walker) = Rp_Rs_ratio
                Ecc_accept_mcmc_array(all_loop, walker) = Ecc
                Inc_accept_mcmc_array(all_loop, walker) = Inc
                omega_arg_periastron_accept_mcmc_array(all_loop, walker) = omega_arg_periastron

                !Rejection flag for the proposal parameters is set to false.
                Reject_flag = 'F'
                
                IF (FLOOR(DBLE(MCMC_loop/100.0D0)) - MCMC_updates >= 1) THEN
                    !Update the MCMC variance, standard deviation, and scale factor at every 100th iteration.
                    vsini_variance_MCMC = vsini_variance
                    vsini_stand_dev_MCMC = vsini_stand_dev
            
                    stellar_rotation_angle_variance_MCMC = stellar_rotation_angle_variance
                    stellar_rotation_angle_stand_dev_MCMC = stellar_rotation_angle_stand_dev
            
                    RV_offset_datasets_variance_MCMC = RV_offset_datasets_variance
                    RV_offset_datasets_stand_dev_MCMC = RV_offset_datasets_stand_dev
               
                    RV_zero_offset_variance_MCMC = RV_zero_offset_variance
                    RV_zero_offset_stand_dev_MCMC = RV_zero_offset_stand_dev
            
                    Orbital_period_variance_MCMC = Orbital_period_variance
                    Orbital_period_stand_dev_MCMC = Orbital_period_stand_dev
            
                    JD_time_mid_transit_variance_MCMC = JD_time_mid_transit_variance
                    JD_time_mid_transit_stand_dev_MCMC = JD_time_mid_transit_stand_dev
            
                    Mp_variance_MCMC = Mp_variance
                    Mp_stand_dev_MCMC = Mp_stand_dev
            
                    Rp_variance_MCMC = Rp_variance
                    Rp_stand_dev_MCMC = Rp_stand_dev
               
                    Ms_solar_variance_MCMC = Ms_solar_variance
                    Ms_solar_stand_dev_MCMC = Ms_solar_stand_dev
            
                    Rs_solar_variance_MCMC = Rs_solar_variance
                    Rs_solar_stand_dev_MCMC = Rs_solar_stand_dev
               
                    Rp_Rs_ratio_variance_MCMC = Rp_Rs_ratio_variance
                    Rp_Rs_ratio_stand_dev_MCMC = Rp_Rs_ratio_stand_dev
            
                    Ecc_variance_MCMC = Ecc_variance
                    Ecc_stand_dev_MCMC = Ecc_stand_dev
            
                    Inc_variance_MCMC = Inc_variance
                    Inc_stand_dev_MCMC = Inc_stand_dev
            
                    omega_arg_periastron_variance_MCMC = omega_arg_periastron_variance
                    omega_arg_periastron_stand_dev_MCMC = omega_arg_periastron_stand_dev
            
                    !Update the scale factor for drawing new proposals. This is based on Cameron et al. 2007.
                    scale_factor = (400.0D0*scale_factor)/DBLE(count_tested_proposals)
            
                    count_tested_proposals = 0
                    MCMC_updates = MCMC_updates + 1         
                END IF
                
                vsini_variance_MCMC_array(all_loop, walker) = vsini_variance_MCMC           !The variance updated every 100 iterations.
                vsini_stand_dev_MCMC_array(all_loop, walker) = vsini_stand_dev_MCMC         !The standard deviation updated every 100 iterations.
                
                stellar_rotation_angle_variance_MCMC_array(all_loop, walker) = stellar_rotation_angle_variance_MCMC       !The variance updated every 100 iterations.
                stellar_rotation_angle_stand_dev_MCMC_array(all_loop, walker) = stellar_rotation_angle_stand_dev_MCMC     !The standard deviation updated every 100 iterations.
                
                
                
                
            ELSE IF (accept_prob < random_draw) THEN
                !The proposal is rejected.
                reject_counter = reject_counter + 1
                !Rejection flag for the proposal parameters is set to true.
                Reject_flag = 'T'
                
                !Grab previous likelihood value.
                culumative_likelihood = culumative_likelihood + likelihood_array_MCMC(all_loop - 1, walker)/(all_loop + culumative_likelihood)
                !culumative_all_likelihood = culumative_likelihood
   
                likelihood_array_MCMC(all_loop, walker) = likelihood_array_MCMC(all_loop - 1, walker)
                culumative_likelihood_array_MCMC(all_loop, walker) = culumative_likelihood
                Chi_squared_array_MCMC(all_loop, walker) = Chi_squared_array_MCMC(all_loop - 1, walker)
                Reduced_Chi_Squared_array_MCMC(all_loop, walker) = Reduced_Chi_Squared_array_MCMC(all_loop - 1, walker)
                
                vsini_sum = vsini_sum + vsini_accept_mcmc_array(all_loop - 1, walker)             
                vsini_mean = vsini_sum/all_loop
                vsini_diff_mean_2 = vsini_diff_mean_2 + abs(vsini_accept_mcmc_array(all_loop - 1, walker) - vsini_mean)**2.0D0
                vsini_variance = vsini_diff_mean_2/(all_loop - 1.0D0)  !sample variance.
                vsini_stand_dev = SQRT(vsini_variance)
            
                IF ((stellar_rotation_angle_mean > 90) .AND. (stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) < -90)) THEN
                    stellar_rotation_angle_sum = stellar_rotation_angle_sum + stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) + 360.0D0
                ELSE IF ((stellar_rotation_angle_mean < -90) .AND. (stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) > 90)) THEN
                    stellar_rotation_angle_sum = stellar_rotation_angle_sum + stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) - 360.0D0
                ELSE  
                    stellar_rotation_angle_sum = stellar_rotation_angle_sum + stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker)
                END IF
                stellar_rotation_angle_mean = stellar_rotation_angle_sum/all_loop
         
                stellar_rotation_angle_diff_mean_2 = stellar_rotation_angle_diff_mean_2 &
                                                     + abs(stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker) - stellar_rotation_angle_mean)**2.0D0
         
                IF (stellar_rotation_angle_mean <= -180) THEN
                    stellar_rotation_angle_mean = 360.0D0 - stellar_rotation_angle_mean
                ELSE IF (stellar_rotation_angle_mean > 180) THEN
                    stellar_rotation_angle_mean = stellar_rotation_angle_mean - 360.0D0
                END IF
         
                stellar_rotation_angle_variance = stellar_rotation_angle_diff_mean_2/(all_loop - 1.0D0)  !sample variance.
                stellar_rotation_angle_stand_dev = SQRT(stellar_rotation_angle_variance)
            
                RV_offset_datasets_sum = RV_offset_datasets_sum + RV_offset_datasets_accept_mcmc_array(all_loop - 1, walker)
                RV_offset_datasets_mean = RV_offset_datasets_sum/all_loop
                RV_offset_datasets_diff_mean_2 = RV_offset_datasets_diff_mean_2 + abs(RV_offset_datasets_accept_mcmc_array(all_loop - 1, walker) - RV_offset_datasets_mean)**2.0D0
                RV_offset_datasets_variance = RV_offset_datasets_diff_mean_2/(all_loop - 1.0D0)  !sample variance.
                RV_offset_datasets_stand_dev = SQRT(RV_offset_datasets_variance)
            
                RV_zero_offset_sum = RV_zero_offset_sum + RV_zero_offset_accept_mcmc_array(all_loop - 1, walker)
                RV_zero_offset_mean = RV_zero_offset_sum/all_loop
                RV_zero_offset_diff_mean_2 = RV_zero_offset_diff_mean_2 + abs(RV_zero_offset_accept_mcmc_array(all_loop - 1, walker) - RV_zero_offset_mean)**2.0D0
                RV_zero_offset_variance = RV_zero_offset_diff_mean_2/(all_loop - 1.0D0) 
                RV_zero_offset_stand_dev = SQRT(RV_zero_offset_variance)
            
                Orbital_period_sum = Orbital_period_sum + Orbital_period_accept_mcmc_array(all_loop - 1, walker)
                Orbital_period_mean = Orbital_period_sum/all_loop
                Orbital_period_diff_mean_2 = Orbital_period_diff_mean_2 + abs(Orbital_period_accept_mcmc_array(all_loop - 1, walker) - Orbital_period_mean)**2.0D0
                Orbital_period_variance = Orbital_period_diff_mean_2/(all_loop - 1.0D0)
                Orbital_period_stand_dev = SQRT(Orbital_period_variance)  

                JD_time_mid_transit_sum = JD_time_mid_transit_sum + JD_time_mid_transit_accept_mcmc_array(all_loop - 1, walker)
                JD_time_mid_transit_mean = JD_time_mid_transit_sum/all_loop
                JD_time_mid_transit_diff_mean_2 = JD_time_mid_transit_diff_mean_2 &
                                                  + abs(JD_time_mid_transit_accept_mcmc_array(all_loop - 1, walker) - JD_time_mid_transit_mean)**2.0D0
                JD_time_mid_transit_variance = JD_time_mid_transit_diff_mean_2/(all_loop - 1.0D0)    
                JD_time_mid_transit_stand_dev = SQRT(JD_time_mid_transit_variance)

                Mp_sum = Mp_sum + Mp_accept_mcmc_array(all_loop - 1, walker)
                Mp_mean = Mp_sum/all_loop
                Mp_diff_mean_2 = Mp_diff_mean_2 + abs(Mp_accept_mcmc_array(all_loop - 1, walker) - Mp_mean)**2.0D0
                Mp_variance = Mp_diff_mean_2/(all_loop - 1.0D0)
                Mp_stand_dev = SQRT(Mp_variance)

                Rp_sum = Rp_sum + Rp_accept_mcmc_array(all_loop - 1, walker)
                Rp_mean = Rp_sum/all_loop
                Rp_diff_mean_2 = Rp_diff_mean_2 + abs(Rp_accept_mcmc_array(all_loop - 1, walker) - Rp_mean)**2.0D0
                Rp_variance = Rp_diff_mean_2/(all_loop - 1.0D0)
                Rp_stand_dev = SQRT(Rp_variance)

                Ecc_sum = Ecc_sum + Ecc_accept_mcmc_array(all_loop - 1, walker)
                Ecc_mean = Ecc_sum/all_loop
                Ecc_diff_mean_2 = Ecc_diff_mean_2 + abs(Ecc_accept_mcmc_array(all_loop - 1, walker) - Ecc_mean)**2.0D0
                Ecc_variance = Ecc_diff_mean_2/(all_loop - 1.0D0)
                Ecc_stand_dev = SQRT(Ecc_variance)

                Inc_sum = Inc_sum + Inc_accept_mcmc_array(all_loop - 1, walker)
                Inc_mean = Inc_sum/all_loop
                Inc_diff_mean_2 = Inc_diff_mean_2 + abs(Inc_accept_mcmc_array(all_loop - 1, walker) - Inc_mean)**2.0D0
                Inc_variance = Inc_diff_mean_2/(all_loop - 1.0D0)
                Inc_stand_dev = SQRT(Inc_variance)

                omega_arg_periastron_sum = omega_arg_periastron_sum + omega_arg_periastron_accept_mcmc_array(all_loop - 1, walker)
                omega_arg_periastron_mean = omega_arg_periastron_sum/all_loop
                omega_arg_periastron_diff_mean_2 = omega_arg_periastron_diff_mean_2 &
                                                   + abs(omega_arg_periastron_accept_mcmc_array(all_loop - 1, walker) - omega_arg_periastron_mean)**2.0D0
                omega_arg_periastron_variance = omega_arg_periastron_diff_mean_2/(all_loop - 1.0D0)     
                omega_arg_periastron_stand_dev = SQRT(omega_arg_periastron_variance)    
            
                Ms_solar_sum = Ms_solar_sum + Ms_solar_accept_mcmc_array(all_loop - 1, walker)
                Ms_solar_mean = Ms_solar_sum/all_loop
                Ms_solar_diff_mean_2 = Ms_solar_diff_mean_2 + abs(Ms_solar_accept_mcmc_array(all_loop - 1, walker) - Ms_solar_mean)**2.0D0
                Ms_solar_variance = Ms_solar_diff_mean_2/(all_loop - 1.0D0)
                Ms_solar_stand_dev = SQRT(Ms_solar_variance)
            
                Rs_solar_sum = Rs_solar_sum + Rs_solar_accept_mcmc_array(all_loop - 1, walker)
                Rs_solar_mean = Rs_solar_sum/all_loop
                Rs_solar_diff_mean_2 = Rs_solar_diff_mean_2 + abs(Rs_solar_accept_mcmc_array(all_loop - 1, walker) - Rs_solar_mean)**2.0D0
                Rs_solar_variance = Rs_solar_diff_mean_2/(all_loop - 1.0D0)
                Rs_solar_stand_dev = SQRT(Rs_solar_variance)
            
                Rp_Rs_ratio_sum = Rp_Rs_ratio_sum + Rp_Rs_ratio_accept_mcmc_array(all_loop - 1, walker)
                Rp_Rs_ratio_mean = Rp_Rs_ratio_sum/all_loop
                Rp_Rs_ratio_diff_mean_2 = Rp_Rs_ratio_diff_mean_2 + abs(Rp_Rs_ratio_accept_mcmc_array(all_loop - 1, walker) - Rp_Rs_ratio_mean)**2.0D0
                Rp_Rs_ratio_variance = Rp_Rs_ratio_diff_mean_2/(all_loop - 1.0D0)
                Rp_Rs_ratio_stand_dev = SQRT(Rp_Rs_ratio_variance)

                vsini_accept_mcmc_array(all_loop, walker) = vsini_accept_mcmc_array(all_loop - 1, walker)
                vsini_variance_current_array(all_loop, walker) = vsini_variance             !The variance updated at every iteration.
                vsini_stand_dev_current_array(all_loop, walker) = vsini_stand_dev           !The standard deviation updated at every iteration.
         
                stellar_rotation_angle_accept_mcmc_array(all_loop, walker) = stellar_rotation_angle_accept_mcmc_array(all_loop - 1, walker)
                stellar_rotation_angle_variance_current_array(all_loop, walker) = stellar_rotation_angle_variance         !The variance updated at every iteration.
                stellar_rotation_angle_stand_dev_current_array(all_loop, walker) = stellar_rotation_angle_stand_dev       !The standard deviation updated at every iteration. 
         
                RV_offset_datasets_accept_mcmc_array(all_loop, walker) = RV_offset_datasets_accept_mcmc_array(all_loop - 1, walker)
                RV_zero_offset_accept_mcmc_array(all_loop, walker) = RV_zero_offset_accept_mcmc_array(all_loop - 1, walker)
                Orbital_period_accept_mcmc_array(all_loop, walker) = Orbital_period_accept_mcmc_array(all_loop - 1, walker)
                JD_time_mid_transit_accept_mcmc_array(all_loop, walker) = JD_time_mid_transit_accept_mcmc_array(all_loop - 1, walker)
                Mp_accept_mcmc_array(all_loop, walker) = Mp_accept_mcmc_array(all_loop - 1, walker)
                Rp_accept_mcmc_array(all_loop, walker) = Rp_accept_mcmc_array(all_loop - 1, walker)
                Ms_solar_accept_mcmc_array(all_loop, walker) = Ms_solar_accept_mcmc_array(all_loop - 1, walker)
                Rs_solar_accept_mcmc_array(all_loop, walker) = Rs_solar_accept_mcmc_array(all_loop - 1, walker)
                Rp_Rs_ratio_accept_mcmc_array(all_loop, walker) = Rp_Rs_ratio_accept_mcmc_array(all_loop - 1, walker)
                Ecc_accept_mcmc_array(all_loop, walker) = Ecc_accept_mcmc_array(all_loop - 1, walker)
                Inc_accept_mcmc_array(all_loop, walker) = Inc_accept_mcmc_array(all_loop - 1, walker)
                omega_arg_periastron_accept_mcmc_array(all_loop, walker) = omega_arg_periastron_accept_mcmc_array(all_loop - 1, walker)
                
                vsini_variance_MCMC_array(all_loop, walker) = vsini_variance_MCMC           !The variance updated every 100 iterations.
                vsini_stand_dev_MCMC_array(all_loop, walker) = vsini_stand_dev_MCMC         !The standard deviation updated every 100 iterations.
                
                stellar_rotation_angle_variance_MCMC_array(all_loop, walker) = stellar_rotation_angle_variance_MCMC       !The variance updated every 100 iterations.
                stellar_rotation_angle_stand_dev_MCMC_array(all_loop, walker) = stellar_rotation_angle_stand_dev_MCMC     !The standard deviation updated every 100 iterations.
            
            END IF
            
        END IF
        
        scale_factor_array(all_loop, walker) = scale_factor
        
        DEALLOCATE(RV_theory)
        DEALLOCATE(data_points_compare)   
        DEALLOCATE(RV_offset_data_array)
        DEALLOCATE(Data_my_rv_offset, data)
        DEALLOCATE(Data_adjusted, sorted_data_array) 
        DEALLOCATE(time_data_array)
        
        accept_rate = DBLE(MCMC_loop)/DBLE(all_loop)
        
        !Give the user some essential information on the progress of the MCMC. Do this for every 25 accepted proposal iterations.
        IF (FLOOR(DBLE(MCMC_loop/25)) - MCMC_print_updates >= 1) THEN
            PRINT *, "Orbital_period_1: ", Orbital_period_1
            PRINT *, "Orbital_period_mean: ", Orbital_period_mean
            PRINT *, "JD_time_mid_transit: ", JD_time_mid_transit
            PRINT *, "JD_time_mid_transit_mean: ", JD_time_mid_transit_mean
            PRINT *, "Semi-major axis: ", Rorb/AU
            PRINT *, "Stellar Semi-major axis: ", Rorb_star/AU
            PRINT *, "Mass of planet: ", Mpp
            PRINT *, "Mp_mean: ", Mp_mean
            PRINT *, "Radius of planet: ", Rpp
            PRINT *, "Rp_mean: ", Rp_mean
            PRINT *, "Mass of star: ", Ms_solar
            PRINT *, "Mp_mean: ", Ms_solar_mean
            PRINT *, "Radius of star: ", Rs_solar
            PRINT *, "Rp_mean: ", Rs_solar_mean
            PRINT *, "Radius of planet to star ratio: ", Rp_Rs_ratio
            PRINT *, "Radius of planet to star ratio mean: ", Rp_Rs_ratio_mean
            PRINT *, "Eccentricity of planet: ", Ecc
            PRINT *, "Ecc_mean: ", Ecc_mean
            PRINT *, "Inclination of planet: ", Inc
            PRINT *, "Inc_mean: ", Inc_mean
            PRINT *, "omega_arg_periastron: ", omega_arg_periastron
            PRINT *, "omega_arg_periastron_mean: ", omega_arg_periastron_mean
            PRINT *, "RV_zero_offset: ", RV_zero_offset
            PRINT *, "RV_zero_offset_mean: ", RV_zero_offset_mean
            PRINT *, "RV_offset_datasets: ", RV_offset_datasets
            PRINT *, "RV_offset_datasets_mean: ", RV_offset_datasets_mean
            PRINT *, "vsini: ", vsini
            PRINT *, "vsini_mean: ", vsini_mean
            PRINT *, "vsini standard deviation MCMC: ", vsini_stand_dev_MCMC
            PRINT *, "vsini standard deviation current: ", vsini_stand_dev
            PRINT *, "spin-orbit angle: ", stellar_rotation_angle
            PRINT *, "stellar_rotation_angle_mean: ", stellar_rotation_angle_mean
            PRINT *, "spin-orbit angle standard deviation MCMC: ", stellar_rotation_angle_stand_dev_MCMC
            PRINT *, "spin-orbit angle standard deviation current: ", stellar_rotation_angle_stand_dev
   
            PRINT *,"There are ", Number_iterations_left_total," interations left out of ", mcmc_accepted_iteration_size
            
            PRINT *,"The MCMC walker number is: ", walker
            PRINT *,"The number of accepted proposals are: ", MCMC_loop
            PRINT *,"The number of tested proposals are: ", all_loop
            PRINT *,"Current proposal acceptance rate: ", accept_rate
            PRINT *,"The number of rejected proposals are: ", reject_counter
            PRINT *,"The cumulative likelihood: ", culumative_likelihood
            PRINT *,"The current likelihood: ", likelihood
            PRINT *,"Previous chi squared value: ", acceptance_param_previous
            PRINT *,"Current chi squared value: ", acceptance_param_current
            PRINT *,"Scale_factor: ", scale_factor
            
            
                
                
            percent_processed = (100.0D0 * (mcmc_accepted_iteration_size * (walker - 1.0D0) + MCMC_loop))/(mcmc_accepted_iteration_size * number_mcmc_walkers)
            string_bar = ''
            string_bar_num = 1
            IF ((percent_processed < 4.0D0) .AND. (percent_processed > 0.05D0)) THEN
                IF ((percent_processed - FLOOR(percent_processed) < 0.2499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.0000001D0)) THEN
                    string_bar(string_bar_num:) = ':'
                END IF
                IF ((percent_processed - FLOOR(percent_processed) < 0.4999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.2499999D0)) THEN
                    string_bar(string_bar_num:) = '/'
                END IF
                IF ((percent_processed - FLOOR(percent_processed) < 0.7499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.4999999D0)) THEN
                    string_bar(string_bar_num:) = '-'
                END IF
                IF ((percent_processed - FLOOR(percent_processed) < 0.9999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.7499999D0)) THEN
                    string_bar(string_bar_num:) = '\'
                END IF
            END IF   
            
            IF (percent_processed >= 4.0D0) THEN
                DO percent_processed_loop = 1, FLOOR(percent_processed/4.0D0)
                    string_bar(string_bar_num:) = '|'
                    string_bar_num = string_bar_num + 1
                    IF ((percent_processed < 99.9D0) .AND. (percent_processed_loop == FLOOR(percent_processed/4.0D0))) THEN
                        IF ((percent_processed - FLOOR(percent_processed) < 0.2499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.0000001D0)) THEN
                            string_bar(string_bar_num:) = ':'
                        END IF
                        IF ((percent_processed - FLOOR(percent_processed) < 0.4999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.2499999D0)) THEN
                            string_bar(string_bar_num:) = '/'
                        END IF
                        IF ((percent_processed - FLOOR(percent_processed) < 0.7499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.4999999D0)) THEN
                            string_bar(string_bar_num:) = '-'
                        END IF
                        IF ((percent_processed - FLOOR(percent_processed) < 0.9999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.7499999D0)) THEN
                            string_bar(string_bar_num:) = '\'
                        END IF
                    END IF
                END DO
            END IF 
            
            PRINT 1000, string_bar, FLOOR(percent_processed)
            1000 FORMAT("Processing data: [", A25, "]", 1X, I3, "% complete")

            PRINT *," "
            PRINT *,"***********************************************************************************************"
            PRINT *," "
      
            MCMC_print_updates = MCMC_print_updates + 1  
        END IF    
        
        mask_array(all_loop, walker) = .TRUE.
         
        IF (Reject_flag == 'F') THEN
            mcmc_loop = mcmc_loop + 1        !Proposal was accepted so increase mcmc_loop variable by 1.
            Number_iterations_left_total = Number_iterations_left_total - 1         !Let's the user know the progress of the program.
        END IF
        all_loop = all_loop + 1    !Increase all_loop counter by one for every iteration.
        total_loop = total_loop + 1
        Ecc_draw_loop = Ecc_draw_loop + 1
        vsini_draw_loop = vsini_draw_loop + 1
        RV_offset_datasets_draw_loop = RV_offset_datasets_draw_loop + 1
        Orbital_period_draw_loop = Orbital_period_draw_loop + 1
        JD_time_mid_transit_draw_loop = JD_time_mid_transit_draw_loop + 1
        Mp_draw_loop = Mp_draw_loop + 1
        Rp_draw_loop = Rp_draw_loop + 1
        Rp_Rs_ratio_draw_loop = Rp_Rs_ratio_draw_loop + 1
        Rs_solar_draw_loop = Rs_solar_draw_loop + 1
        Ms_solar_draw_loop = Ms_solar_draw_loop + 1
        Inc_draw_loop = Inc_draw_loop + 1
        omega_arg_periastron_draw_loop = omega_arg_periastron_draw_loop + 1
        RV_zero_offset_draw_loop = RV_zero_offset_draw_loop + 1
        stellar_rotation_angle_draw_loop = stellar_rotation_angle_draw_loop + 1
        
    END DO
    
!**************************************************************End of MCMC Walker iteration***********************************************************************************!    
    
    walker_parameter_array(1, walker) = reject_counter
    walker_parameter_array(2, walker) = all_loop - 1
    walker_parameter_array(3, walker) = MCMC_loop - 1
    walker_parameter_array(4, walker) = accept_rate
    walker_parameter_array(5, walker) = scale_factor
    
    PRINT *, "number of accepted MCMC iterations: ", MCMC_loop - 1
    PRINT *, "number of iterations in total: ", all_loop - 1
    PRINT *, "total_loop: ", total_loop - 1
    PRINT *, "number of rejected proposals: ", reject_counter
    PRINT *, "overall acceptance rate: ", accept_rate

    PRINT *," "
    PRINT *,"***********************************************************************************************"
    PRINT *," "
 
    !Determine the lowest value for Chi square, reduced Chi square, and where these occur. Also determine the greatest value for the likelihood
    !and where this occurs. Then determine the best proposal parameter values based on the mean value from all the accepted proposal runs
    !and from the minimum value of Chi square. The uncertainties in the proposal parameters will be determined from the standard deviation
    !and from using the covariance matrix method. Report the values of the priors that correspond to the lowest value of Chi square and
    !maximum value of the likelihood.
    
    min_chi_squared_MCMC =  minval(Chi_squared_array_MCMC, dim=1, mask=mask_array)
    print *, "Mininum chi squared value for walker number ", walker, "of the MCMC is: ", min_chi_squared_MCMC(walker)
    loc_min_chi_squared_MCMC = minloc(Chi_squared_array_MCMC, dim=1, mask=mask_array)
    print *, "Mininum chi squared location for walker number ", walker, "of the MCMC is: ", loc_min_chi_squared_MCMC(walker)
    min_reduced_chi_squared_MCMC = minval(Reduced_Chi_Squared_array_MCMC, dim=1, mask=mask_array)
    print *, "Mininum reduced chi squared value for walker number ", walker, "of the MCMC is: ", min_reduced_chi_squared_MCMC(walker)
    loc_min_reduced_chi_squared_MCMC = minloc(Reduced_Chi_Squared_array_MCMC, dim=1, mask=mask_array)
    print *, "Mininum reduced chi squared location for walker number ", walker, "of the MCMC is: ", loc_min_reduced_chi_squared_MCMC(walker)

    max_likelihood_MCMC =  maxval(likelihood_array_MCMC, dim=1, mask=mask_array)
    print *, "Maximum likelihood value for walker number ", walker, "of the MCMC is: ", max_likelihood_MCMC(walker)
    loc_max_likelihood_MCMC = maxloc(likelihood_array_MCMC, dim=1, mask=mask_array)
    print *, "Maximum likelihood location for walker number ", walker, "of the MCMC is: ", loc_max_likelihood_MCMC(walker)
    
    !Best fit priors and proposal parameter values based on the location of the minimum value of Chi squared.
    best_RV_offset_datasets_chi2 = RV_offset_datasets_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(1, walker) = best_RV_offset_datasets_chi2
    PRINT *, "The best RV offset between datasets corresponding to the lowest value of chi squared: ", best_RV_offset_datasets_chi2, "m/s", " for walker number ", walker
    
    best_RV_zero_offset_chi2 = RV_zero_offset_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(2, walker) = best_RV_zero_offset_chi2
    PRINT *, "The best RV zero offset velocity corresponding to the lowest value of chi squared: ", best_RV_zero_offset_chi2, "m/s", " for walker number ", walker
    
    best_orbital_period_chi2 = Orbital_period_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(3, walker) = best_orbital_period_chi2
    PRINT *, "The best orbital period corresponding to the lowest value of chi squared: ", best_orbital_period_chi2, "days", " for walker number ", walker
    
    best_JD_time_mid_transit_chi2 = JD_time_mid_transit_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(4, walker) = best_JD_time_mid_transit_chi2
    PRINT *, "The best mid transit time corresponding to the lowest value of chi squared: ", best_JD_time_mid_transit_chi2, "BJD", " for walker number ", walker

    best_Mp_chi2 = Mp_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(5, walker) = best_Mp_chi2
    PRINT *, "The best planetary mass corresponding to the lowest value of chi squared: ", best_Mp_chi2, "M_J", " for walker number ", walker
    
    best_Rp_chi2 = Rp_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(6, walker) = best_Rp_chi2
    PRINT *, "The best planetary radius corresponding to the lowest value of chi squared: ", best_Rp_chi2, "R_J", " for walker number ", walker
    
    best_Ms_solar_chi2 = Ms_solar_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(7, walker) = best_Ms_solar_chi2
    PRINT *, "The best stellar mass corresponding to the lowest value of chi squared: ", best_Ms_solar_chi2, "M_Sol", " for walker number ", walker
    
    best_Rs_solar_chi2 = Rs_solar_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(8, walker) = best_Rs_solar_chi2
    PRINT *, "The best stellar radius corresponding to the lowest value of chi squared: ", best_Rs_solar_chi2, "R_Sol", " for walker number ", walker
    
    best_Rp_Rs_ratio_chi2 = Rp_Rs_ratio_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(9, walker) = best_Rp_Rs_ratio_chi2
    PRINT *, "The best planet to stellar radius ratio corresponding to the lowest value of chi squared: ", best_Rp_Rs_ratio_chi2, " for walker number ", walker
    
    best_Ecc_chi2 = Ecc_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(10, walker) = best_Ecc_chi2
    PRINT *, "The best eccentricity value corresponding to the lowest value of chi squared: ", best_Ecc_chi2, " for walker number ", walker
    
    best_Inc_chi2 = Inc_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(11, walker) = best_Inc_chi2
    PRINT *, "The best inclination angle corresponding to the lowest value of chi squared: ", best_Inc_chi2, "deg", " for walker number ", walker
    
    best_omega_arg_periastron_chi2 = omega_arg_periastron_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(12, walker) = best_omega_arg_periastron_chi2
    PRINT *, "The best argument of periastron angle corresponding to the lowest value of chi squared: ", best_omega_arg_periastron_chi2, "deg", " for walker number ", walker
    
    best_vsini_chi2 = vsini_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(13, walker) = best_vsini_chi2
    PRINT *, "The best vsini value corresponding to the lowest value of chi squared: ", best_vsini_chi2, "m/s", " for walker number ", walker
    
    best_spin_orbit_chi2 = stellar_rotation_angle_accept_mcmc_array(loc_min_chi_squared_MCMC(walker), walker)
    walker_parameter_array_chi2(14, walker) = best_spin_orbit_chi2
    PRINT *, "The best spin-orbit angle corresponding to the lowest value of chi squared: ", best_spin_orbit_chi2, "deg", " for walker number ", walker
    
    
    
    
    !Determine the best vsini and lambda from the mean of the MCMC runs.
    !vsini_sum and stellar_rotation_angle_mean are calculated in the MCMC loop.
    RV_offset_datasets_sum2 = sum(RV_offset_datasets_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "RV_offset_datasets_sum2: ", RV_offset_datasets_sum2
    RV_offset_datasets_mean = RV_offset_datasets_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(1, walker) = RV_offset_datasets_mean
    
    RV_zero_offset_sum2 = sum(RV_zero_offset_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "RV_zero_offset_sum2: ", RV_zero_offset_sum2
    RV_zero_offset_mean = RV_zero_offset_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(2, walker) = RV_zero_offset_mean
    
    Orbital_period_sum2 = sum(Orbital_period_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Orbital_period_sum2: ", Orbital_period_sum2
    Orbital_period_mean = Orbital_period_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(3, walker) = Orbital_period_mean
    
    JD_time_mid_transit_sum2 = sum(JD_time_mid_transit_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "JD_time_mid_transit_sum2: ", JD_time_mid_transit_sum2
    JD_time_mid_transit_mean = JD_time_mid_transit_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(4, walker) = JD_time_mid_transit_mean
    
    Mp_sum2 = sum(Mp_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Mp_sum2: ", Mp_sum2
    Mp_mean = Mp_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(5, walker) = Mp_mean
    
    Rp_sum2 = sum(Rp_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Rp_sum2: ", Rp_sum2
    Rp_mean = Rp_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(6, walker) = Rp_mean
    
    Ms_solar_sum2 = sum(Ms_solar_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Ms_solar_sum2: ", Ms_solar_sum2
    Ms_solar_mean = Ms_solar_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(7, walker) = Ms_solar_mean
    
    Rs_solar_sum2 = sum(Rs_solar_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Rs_solar_sum2: ", Rs_solar_sum2
    Rs_solar_mean = Rs_solar_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(8, walker) = Rs_solar_mean
    
    Rp_Rs_ratio_sum2 = sum(Rp_Rs_ratio_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Rp_Rs_ratio_sum2: ", Rp_Rs_ratio_sum2
    Rp_Rs_ratio_mean = Rp_Rs_ratio_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(9, walker) = Rp_Rs_ratio_mean
    
    Ecc_sum2 = sum(Ecc_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Ecc_sum2: ", Ecc_sum2
    Ecc_mean = Ecc_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(10, walker) = Ecc_mean
    
    Inc_sum2 = sum(Inc_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "Inc_sum2: ", Inc_sum2
    Inc_mean = Inc_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(11, walker) = Inc_mean
    
    omega_arg_periastron_sum2 = sum(omega_arg_periastron_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "omega_arg_periastron_sum2: ", omega_arg_periastron_sum2
    omega_arg_periastron_mean = omega_arg_periastron_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(12, walker) = omega_arg_periastron_mean
    
    vsini_sum2 = sum(vsini_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "vsini_sum2: ", vsini_sum2
    vsini_mean = vsini_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(13, walker) = vsini_mean
    
    stellar_rotation_angle_sum2 = sum(stellar_rotation_angle_accept_mcmc_array, dim=1, mask=mask_array)
    PRINT *, "stellar_rotation_angle_sum2: ", stellar_rotation_angle_sum2
    stellar_rotation_angle_mean = stellar_rotation_angle_sum2(walker)/(all_loop - 1.0D0)
    walker_parameter_array_mean(14, walker) = stellar_rotation_angle_mean
    
    
    
    
    !Calculate the standard deviation from the unweighted means.
    vsini_variance = 0.0D0  !reset to zero.
    stellar_rotation_angle_variance = 0.0D0
    RV_offset_datasets_variance = 0.0D0
    RV_zero_offset_variance = 0.0D0
    Orbital_period_variance = 0.0D0
    JD_time_mid_transit_variance = 0.0D0
    Mp_variance = 0.0D0
    Rp_variance = 0.0D0
    Ms_solar_variance = 0.0D0
    Rs_solar_variance = 0.0D0
    Rp_Rs_ratio_variance = 0.0D0
    Ecc_variance = 0.0D0
    Inc_variance = 0.0D0
    omega_arg_periastron_variance = 0.0D0
    DO l = 1, all_loop - 1
        vsini_variance = vsini_variance + abs(vsini_accept_mcmc_array(l, walker) - vsini_mean)**2.0D0
        stellar_rotation_angle_variance = stellar_rotation_angle_variance + abs(stellar_rotation_angle_accept_mcmc_array(l, walker) &
                                          - stellar_rotation_angle_mean)**2.0D0
        RV_offset_datasets_variance = RV_offset_datasets_variance + abs(RV_offset_datasets_accept_mcmc_array(l, walker) &
                                      - RV_offset_datasets_mean)**2.0D0
        RV_zero_offset_variance = RV_zero_offset_variance + abs(RV_zero_offset_accept_mcmc_array(l, walker) - RV_zero_offset_mean)**2.0D0
        Orbital_period_variance = Orbital_period_variance + abs(Orbital_period_accept_mcmc_array(l, walker) - Orbital_period_mean)**2.0D0
        JD_time_mid_transit_variance = JD_time_mid_transit_variance + abs(JD_time_mid_transit_accept_mcmc_array(l, walker) &
                                       - JD_time_mid_transit_mean)**2.0D0
        Mp_variance = Mp_variance + abs(Mp_accept_mcmc_array(l, walker) - Mp_mean)**2.0D0
        Rp_variance = Rp_variance + abs(Rp_accept_mcmc_array(l, walker) - Rp_mean)**2.0D0
        Ms_solar_variance = Ms_solar_variance + abs(Ms_solar_accept_mcmc_array(l, walker) - Ms_solar_mean)**2.0D0
        Rs_solar_variance = Rs_solar_variance + abs(Rs_solar_accept_mcmc_array(l, walker) - Rs_solar_mean)**2.0D0
        Rp_Rs_ratio_variance = Rp_Rs_ratio_variance + abs(Rp_Rs_ratio_accept_mcmc_array(l, walker) - Rp_Rs_ratio_mean)**2.0D0
        Ecc_variance = Ecc_variance + abs(Ecc_accept_mcmc_array(l, walker) - Ecc_mean)**2.0D0
        Inc_variance = Inc_variance + abs(Inc_accept_mcmc_array(l, walker) - Inc_mean)**2.0D0
        omega_arg_periastron_variance = omega_arg_periastron_variance + abs(omega_arg_periastron_accept_mcmc_array(l, walker) - omega_arg_periastron_mean)**2.0D0
    END DO
    
    walker_parameter_array_variance(1, walker) = RV_offset_datasets_variance
    walker_parameter_array_variance(2, walker) = RV_zero_offset_variance
    walker_parameter_array_variance(3, walker) = Orbital_period_variance
    walker_parameter_array_variance(4, walker) = JD_time_mid_transit_variance
    walker_parameter_array_variance(5, walker) = Mp_variance
    walker_parameter_array_variance(7, walker) = Ms_solar_variance
    walker_parameter_array_variance(8, walker) = Rs_solar_variance
    walker_parameter_array_variance(9, walker) = Rp_Rs_ratio_variance
    walker_parameter_array_variance(10, walker) = Ecc_variance
    walker_parameter_array_variance(11, walker) = Inc_variance
    walker_parameter_array_variance(12, walker) = omega_arg_periastron_variance
    walker_parameter_array_variance(13, walker) = vsini_variance
    walker_parameter_array_variance(14, walker) = stellar_rotation_angle_variance
    
    
    
    
    vsini_stand_dev = SQRT(vsini_variance/(all_loop - 1.0D0))  !sample variance.
    stellar_rotation_angle_stand_dev = SQRT(stellar_rotation_angle_variance/(all_loop - 1.0D0))  !sample variance.
    RV_offset_datasets_stand_dev = SQRT(RV_offset_datasets_variance/(all_loop - 1.0D0))
    RV_zero_offset_stand_dev = SQRT(RV_zero_offset_variance/(all_loop - 1.0D0))
    Orbital_period_stand_dev = SQRT(Orbital_period_variance/(all_loop - 1.0D0))
    JD_time_mid_transit_stand_dev = SQRT(JD_time_mid_transit_variance/(all_loop - 1.0D0))
    Mp_stand_dev = SQRT(Mp_variance/(all_loop - 1.0D0))
    Rp_stand_dev = SQRT(Rp_variance/(all_loop - 1.0D0))
    Ms_solar_stand_dev = SQRT(Ms_solar_variance/(all_loop - 1.0D0))
    Rs_solar_stand_dev = SQRT(Rs_solar_variance/(all_loop - 1.0D0))
    Rp_Rs_ratio_stand_dev = SQRT(Rp_Rs_ratio_variance/(all_loop - 1.0D0))
    Ecc_stand_dev = SQRT(Ecc_variance/(all_loop - 1.0D0))
    Inc_stand_dev = SQRT(Inc_variance/(all_loop - 1.0D0))
    omega_arg_periastron_stand_dev = SQRT(omega_arg_periastron_variance/(all_loop - 1.0D0))
    
    walker_parameter_array_stand_dev(1, walker) = RV_offset_datasets_stand_dev
    walker_parameter_array_stand_dev(2, walker) = RV_zero_offset_stand_dev
    walker_parameter_array_stand_dev(3, walker) = Orbital_period_stand_dev
    walker_parameter_array_stand_dev(4, walker) = JD_time_mid_transit_stand_dev
    walker_parameter_array_stand_dev(5, walker) = Mp_stand_dev
    walker_parameter_array_stand_dev(7, walker) = Ms_solar_stand_dev
    walker_parameter_array_stand_dev(8, walker) = Rs_solar_stand_dev
    walker_parameter_array_stand_dev(9, walker) = Rp_Rs_ratio_stand_dev
    walker_parameter_array_stand_dev(10, walker) = Ecc_stand_dev
    walker_parameter_array_stand_dev(11, walker) = Inc_stand_dev
    walker_parameter_array_stand_dev(12, walker) = omega_arg_periastron_stand_dev
    walker_parameter_array_stand_dev(13, walker) = vsini_stand_dev
    walker_parameter_array_stand_dev(14, walker) = stellar_rotation_angle_stand_dev

    PRINT *, 'vsini_mean: ', vsini_mean
    PRINT *, 'vsini_stand_dev: ', vsini_stand_dev
    
    PRINT *, 'stellar_rotation_angle_mean: ', stellar_rotation_angle_mean
    PRINT *, 'stellar_rotation_angle_stand_dev: ', stellar_rotation_angle_stand_dev

    PRINT *, 'RV_offset_datasets_mean: ', RV_offset_datasets_mean
    PRINT *, 'RV_offset_datasets_stand_dev: ', RV_offset_datasets_stand_dev
    
    PRINT *, 'RV_zero_offset_mean: ', RV_zero_offset_mean
    PRINT *, 'RV_zero_offset_stand_dev: ', RV_zero_offset_stand_dev

    PRINT *, 'Orbital_period_mean: ', Orbital_period_mean
    PRINT *, 'Orbital_period_stand_dev: ', Orbital_period_stand_dev

    PRINT *, 'JD_time_mid_transit_mean: ', JD_time_mid_transit_mean
    PRINT *, 'JD_time_mid_transit_stand_dev: ', JD_time_mid_transit_stand_dev

    PRINT *, 'Mp_mean: ', Mp_mean
    PRINT *, 'Mp_stand_dev: ', Mp_stand_dev

    PRINT *, 'Rp_mean: ', Rp_mean
    PRINT *, 'Rp_stand_dev: ', Rp_stand_dev
    
    PRINT *, 'Ms_solar_mean: ', Ms_solar_mean
    PRINT *, 'Ms_solar_stand_dev: ', Ms_solar_stand_dev

    PRINT *, 'Rs_solar_mean: ', Rs_solar_mean
    PRINT *, 'Rs_solar_stand_dev: ', Rs_solar_stand_dev
    
    PRINT *, 'Rp_Rs_ratio_mean: ', Rp_Rs_ratio_mean
    PRINT *, 'Rp_Rs_ratio_stand_dev: ', Rp_Rs_ratio_stand_dev

    PRINT *, 'Ecc_mean: ', Ecc_mean
    PRINT *, 'Ecc_stand_dev: ', Ecc_stand_dev

    PRINT *, 'Inc_mean: ', Inc_mean
    PRINT *, 'Inc_stand_dev: ', Inc_stand_dev

    PRINT *, 'omega_arg_periastron_mean: ', omega_arg_periastron_mean
    PRINT *, 'omega_arg_periastron_stand_dev: ', omega_arg_periastron_stand_dev
    


    
    PRINT *, ''

    PRINT *," "
    PRINT *,"***********************************************************************************************"
    PRINT *," "
    
    !stop

END DO

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

cmd1 = "say 'Finished Markov Chain Monte Carlo iterations.'"
cmd1 = TRIM(cmd1)
call execute_command_line(cmd1, wait=.false., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)

!************************************************************** MCMC end ***********************************************************************************!
!Determine the mean and uncertainties from all the MCMC walkers. Then write out results to text files.

!Find the minimum chi squared value of the from the whole parameter space sampled by all the MCMC walkers.
min_chi_squared_all_walkers =  minval(Chi_squared_array_MCMC, mask=mask_array)
print *, "Mininum chi squared value for whole MCMC parameter space is: ", min_chi_squared_all_walkers
loc_min_chi_squared_all_walkers = minloc(Chi_squared_array_MCMC, mask=mask_array)
print *, "Mininum chi squared location for whole MCMC parameter space is: (", loc_min_chi_squared_all_walkers(1), ',', loc_min_chi_squared_all_walkers(2), ')'
min_reduced_chi_squared_all_walkers = minval(Reduced_Chi_Squared_array_MCMC, mask=mask_array)
print *, "Mininum reduced chi squared value for whole MCMC parameter space is: ", min_reduced_chi_squared_all_walkers
loc_min_reduced_chi_squared_all_walkers = minloc(Reduced_Chi_Squared_array_MCMC, mask=mask_array)
print *, "Mininum reduced chi squared location for whole MCMC parameter space is: (", loc_min_reduced_chi_squared_all_walkers(1), ',', loc_min_reduced_chi_squared_all_walkers(2), ')'

max_likelihood_all_walkers =  maxval(likelihood_array_MCMC, mask=mask_array)
print *, "Maximum likelihood value for whole MCMC parameter space is: ", max_likelihood_all_walkers
loc_max_likelihood_all_walkers = maxloc(likelihood_array_MCMC, mask=mask_array)
print *, "Maximum likelihood location for whole MCMC parameter space is: (", loc_max_likelihood_all_walkers(1), ',', loc_max_likelihood_all_walkers(2), ')'

!Best fit priors and proposal parameter values based on the location of the minimum value of Chi squared.
best_RV_offset_datasets_chi2_all_walkers = RV_offset_datasets_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best RV offset between datasets corresponding to the lowest value of chi squared: ", best_RV_offset_datasets_chi2_all_walkers, "m/s"

best_RV_zero_offset_chi2_all_walkers = RV_zero_offset_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best RV zero offset velocity corresponding to the lowest value of chi squared: ", best_RV_zero_offset_chi2_all_walkers, "m/s"

best_orbital_period_chi2_all_walkers = Orbital_period_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best orbital period corresponding to the lowest value of chi squared: ", best_orbital_period_chi2_all_walkers, "days"

best_JD_time_mid_transit_chi2_all_walkers = JD_time_mid_transit_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best mid transit time corresponding to the lowest value of chi squared: ", best_JD_time_mid_transit_chi2_all_walkers, "BJD"

best_Mp_chi2_all_walkers = Mp_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best planetary mass corresponding to the lowest value of chi squared: ", best_Mp_chi2_all_walkers, "M_J"

best_Rp_chi2_all_walkers = Rp_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best planetary radius corresponding to the lowest value of chi squared: ", best_Rp_chi2_all_walkers, "R_J"

best_Ms_solar_chi2_all_walkers = Ms_solar_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best stellar mass corresponding to the lowest value of chi squared: ", best_Ms_solar_chi2_all_walkers, "M_Sol"

best_Rs_solar_chi2_all_walkers = Rs_solar_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best stellar radius corresponding to the lowest value of chi squared: ", best_Rs_solar_chi2_all_walkers, "R_Sol"

best_Rp_Rs_ratio_chi2_all_walkers = Rp_Rs_ratio_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best planet to stellar radius ratio corresponding to the lowest value of chi squared: ", best_Rp_Rs_ratio_chi2_all_walkers

best_Ecc_chi2_all_walkers = Ecc_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best eccentricity value corresponding to the lowest value of chi squared: ", best_Ecc_chi2_all_walkers

best_Inc_chi2_all_walkers = Inc_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best inclination angle corresponding to the lowest value of chi squared: ", best_Inc_chi2_all_walkers, "deg"

best_omega_arg_periastron_chi2_all_walkers = omega_arg_periastron_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best argument of periastron angle corresponding to the lowest value of chi squared: ", best_omega_arg_periastron_chi2_all_walkers, "deg"

best_vsini_chi2_all_walkers = vsini_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best vsini value corresponding to the lowest value of chi squared: ", best_vsini_chi2_all_walkers, "m/s"

best_spin_orbit_chi2_all_walkers = stellar_rotation_angle_accept_mcmc_array(loc_min_chi_squared_all_walkers(1), loc_min_chi_squared_all_walkers(2))
PRINT *, "The best spin-orbit angle corresponding to the lowest value of chi squared: ", best_spin_orbit_chi2_all_walkers, "deg"




!Determine the best vsini and lambda from the mean of the all MCMC walker chains.
!vsini_sum and stellar_rotation_angle_mean are calculated in the MCMC loop.
RV_offset_datasets_sum3 = sum(RV_offset_datasets_accept_mcmc_array, mask=mask_array)
PRINT *, "RV_offset_datasets_sum3: ", RV_offset_datasets_sum3
RV_offset_datasets_mean2 = RV_offset_datasets_sum3/(total_loop - 1.0D0)

RV_zero_offset_sum3 = sum(RV_zero_offset_accept_mcmc_array, mask=mask_array)
PRINT *, "RV_zero_offset_sum3: ", RV_zero_offset_sum3
RV_zero_offset_mean2 = RV_zero_offset_sum3/(total_loop - 1.0D0)

Orbital_period_sum3 = sum(Orbital_period_accept_mcmc_array, mask=mask_array)
PRINT *, "Orbital_period_sum3: ", Orbital_period_sum3
Orbital_period_mean2 = Orbital_period_sum3/(total_loop - 1.0D0)

JD_time_mid_transit_sum3 = sum(JD_time_mid_transit_accept_mcmc_array, mask=mask_array)
PRINT *, "JD_time_mid_transit_sum3: ", JD_time_mid_transit_sum3
JD_time_mid_transit_mean2 = JD_time_mid_transit_sum3/(total_loop - 1.0D0)

Mp_sum3 = sum(Mp_accept_mcmc_array, mask=mask_array)
PRINT *, "Mp_sum3: ", Mp_sum3
Mp_mean2 = Mp_sum3/(total_loop - 1.0D0)

Rp_sum3 = sum(Rp_accept_mcmc_array, mask=mask_array)
PRINT *, "Rp_sum3: ", Rp_sum3
Rp_mean2 = Rp_sum3/(total_loop - 1.0D0)

Ms_solar_sum3 = sum(Ms_solar_accept_mcmc_array, mask=mask_array)
PRINT *, "Ms_solar_sum3: ", Ms_solar_sum3
Ms_solar_mean2 = Ms_solar_sum3/(total_loop - 1.0D0)

Rs_solar_sum3 = sum(Rs_solar_accept_mcmc_array, mask=mask_array)
PRINT *, "Rs_solar_sum3: ", Rs_solar_sum3
Rs_solar_mean2 = Rs_solar_sum3/(total_loop - 1.0D0)

Rp_Rs_ratio_sum3 = sum(Rp_Rs_ratio_accept_mcmc_array, mask=mask_array)
PRINT *, "Rp_Rs_ratio_sum3: ", Rp_Rs_ratio_sum3
Rp_Rs_ratio_mean2 = Rp_Rs_ratio_sum3/(total_loop - 1.0D0)

Ecc_sum3 = sum(Ecc_accept_mcmc_array, mask=mask_array)
PRINT *, "Ecc_sum3: ", Ecc_sum3
Ecc_mean2 = Ecc_sum3/(total_loop - 1.0D0)

Inc_sum3 = sum(Inc_accept_mcmc_array, mask=mask_array)
PRINT *, "Inc_sum3: ", Inc_sum3
Inc_mean2 = Inc_sum3/(total_loop - 1.0D0)

omega_arg_periastron_sum3 = sum(omega_arg_periastron_accept_mcmc_array, mask=mask_array)
PRINT *, "omega_arg_periastron_sum3: ", omega_arg_periastron_sum3
omega_arg_periastron_mean2 = omega_arg_periastron_sum3/(total_loop - 1.0D0)

vsini_sum3 = sum(vsini_accept_mcmc_array, mask=mask_array)
PRINT *, "vsini_sum3: ", vsini_sum3
vsini_mean2 = vsini_sum3/(total_loop - 1.0D0)

stellar_rotation_angle_sum3 = sum(stellar_rotation_angle_accept_mcmc_array, mask=mask_array)
PRINT *, "stellar_rotation_angle_sum3: ", stellar_rotation_angle_sum3
stellar_rotation_angle_mean2 = stellar_rotation_angle_sum3/(total_loop - 1.0D0)




!Calculate the standard deviation from the unweighted means.
vsini_variance = 0.0D0  !reset to zero.
stellar_rotation_angle_variance = 0.0D0
RV_offset_datasets_variance = 0.0D0
RV_zero_offset_variance = 0.0D0
Orbital_period_variance = 0.0D0
JD_time_mid_transit_variance = 0.0D0
Mp_variance = 0.0D0
Rp_variance = 0.0D0
Ms_solar_variance = 0.0D0
Rs_solar_variance = 0.0D0
Rp_Rs_ratio_variance = 0.0D0
Ecc_variance = 0.0D0
Inc_variance = 0.0D0
omega_arg_periastron_variance = 0.0D0
DO i = 1, number_mcmc_walkers
    DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
        IF (mask_array(l, i) .EQV. .TRUE.) THEN
            vsini_variance = vsini_variance + abs(vsini_accept_mcmc_array(l, i) - vsini_mean2)**2.0D0
            stellar_rotation_angle_variance = stellar_rotation_angle_variance + abs(stellar_rotation_angle_accept_mcmc_array(l, i) &
                                              - stellar_rotation_angle_mean2)**2.0D0
            RV_offset_datasets_variance = RV_offset_datasets_variance + abs(RV_offset_datasets_accept_mcmc_array(l, i) &
                                          - RV_offset_datasets_mean2)**2.0D0
            RV_zero_offset_variance = RV_zero_offset_variance + abs(RV_zero_offset_accept_mcmc_array(l, i) - RV_zero_offset_mean2)**2.0D0
            Orbital_period_variance = Orbital_period_variance + abs(Orbital_period_accept_mcmc_array(l, i) - Orbital_period_mean2)**2.0D0
            JD_time_mid_transit_variance = JD_time_mid_transit_variance + abs(JD_time_mid_transit_accept_mcmc_array(l, i) &
                                           - JD_time_mid_transit_mean2)**2.0D0
            Mp_variance = Mp_variance + abs(Mp_accept_mcmc_array(l, i) - Mp_mean2)**2.0D0
            Rp_variance = Rp_variance + abs(Rp_accept_mcmc_array(l, i) - Rp_mean2)**2.0D0
            Ms_solar_variance = Ms_solar_variance + abs(Ms_solar_accept_mcmc_array(l, i) - Ms_solar_mean2)**2.0D0
            Rs_solar_variance = Rs_solar_variance + abs(Rs_solar_accept_mcmc_array(l, i) - Rs_solar_mean2)**2.0D0
            Rp_Rs_ratio_variance = Rp_Rs_ratio_variance + abs(Rp_Rs_ratio_accept_mcmc_array(l, i) - Rp_Rs_ratio_mean2)**2.0D0
            Ecc_variance = Ecc_variance + abs(Ecc_accept_mcmc_array(l, i) - Ecc_mean2)**2.0D0
            Inc_variance = Inc_variance + abs(Inc_accept_mcmc_array(l, i) - Inc_mean2)**2.0D0
            omega_arg_periastron_variance = omega_arg_periastron_variance + abs(omega_arg_periastron_accept_mcmc_array(l, i) - omega_arg_periastron_mean2)**2.0D0
        END IF
    END DO
END DO




vsini_stand_dev = SQRT(vsini_variance/(total_loop - 1.0D0))  !sample variance.
stellar_rotation_angle_stand_dev = SQRT(stellar_rotation_angle_variance/(total_loop - 1.0D0))  !sample variance.
RV_offset_datasets_stand_dev = SQRT(RV_offset_datasets_variance/(total_loop - 1.0D0))
RV_zero_offset_stand_dev = SQRT(RV_zero_offset_variance/(total_loop - 1.0D0))
Orbital_period_stand_dev = SQRT(Orbital_period_variance/(total_loop - 1.0D0))
JD_time_mid_transit_stand_dev = SQRT(JD_time_mid_transit_variance/(total_loop - 1.0D0))
Mp_stand_dev = SQRT(Mp_variance/(total_loop - 1.0D0))
Rp_stand_dev = SQRT(Rp_variance/(total_loop - 1.0D0))
Ms_solar_stand_dev = SQRT(Ms_solar_variance/(total_loop - 1.0D0))
Rs_solar_stand_dev = SQRT(Rs_solar_variance/(total_loop - 1.0D0))
Rp_Rs_ratio_stand_dev = SQRT(Rp_Rs_ratio_variance/(total_loop - 1.0D0))
Ecc_stand_dev = SQRT(Ecc_variance/(total_loop - 1.0D0))
Inc_stand_dev = SQRT(Inc_variance/(total_loop - 1.0D0))
omega_arg_periastron_stand_dev = SQRT(omega_arg_periastron_variance/(total_loop - 1.0D0))

PRINT *, 'vsini_mean all walkers: ', vsini_mean2
PRINT *, 'vsini_stand_dev all walkers: ', vsini_stand_dev

PRINT *, 'stellar_rotation_angle_mean all walkers: ', stellar_rotation_angle_mean2
PRINT *, 'stellar_rotation_angle_stand_dev all walkers: ', stellar_rotation_angle_stand_dev

PRINT *, 'RV_offset_datasets_mean all walkers: ', RV_offset_datasets_mean2
PRINT *, 'RV_offset_datasets_stand_dev all walkers: ', RV_offset_datasets_stand_dev

PRINT *, 'RV_zero_offset_mean all walkers: ', RV_zero_offset_mean2
PRINT *, 'RV_zero_offset_stand_dev all walkers: ', RV_zero_offset_stand_dev

PRINT *, 'Orbital_period_mean all walkers: ', Orbital_period_mean2
PRINT *, 'Orbital_period_stand_dev all walkers: ', Orbital_period_stand_dev

PRINT *, 'JD_time_mid_transit_mean all walkers: ', JD_time_mid_transit_mean2
PRINT *, 'JD_time_mid_transit_stand_dev all walkers: ', JD_time_mid_transit_stand_dev

PRINT *, 'Mp_mean all walkers: ', Mp_mean2
PRINT *, 'Mp_stand_dev all walkers: ', Mp_stand_dev

PRINT *, 'Rp_mean all walkers: ', Rp_mean2
PRINT *, 'Rp_stand_dev all walkers: ', Rp_stand_dev

PRINT *, 'Ms_solar_mean all walkers: ', Ms_solar_mean2
PRINT *, 'Ms_solar_stand_dev all walkers: ', Ms_solar_stand_dev

PRINT *, 'Rs_solar_mean all walkers: ', Rs_solar_mean2
PRINT *, 'Rs_solar_stand_dev all walkers: ', Rs_solar_stand_dev

PRINT *, 'Rp_Rs_ratio_mean all walkers: ', Rp_Rs_ratio_mean2
PRINT *, 'Rp_Rs_ratio_stand_dev all walkers: ', Rp_Rs_ratio_stand_dev

PRINT *, 'Ecc_mean all walkers: ', Ecc_mean2
PRINT *, 'Ecc_stand_dev all walkers: ', Ecc_stand_dev

PRINT *, 'Inc_mean all walkers: ', Inc_mean2
PRINT *, 'Inc_stand_dev all walkers: ', Inc_stand_dev

PRINT *, 'omega_arg_periastron_mean all walkers: ', omega_arg_periastron_mean2
PRINT *, 'omega_arg_periastron_stand_dev all walkers: ', omega_arg_periastron_stand_dev
    
    
    
                                        
!STOP

!*****************************************************************************************************************************************

!Write the mask array indicating which values in the arrays to use.
OPEN(unit=99, FILE=output_mask_array_filename, status='replace', action='write')
!Now write array line by line into chi squared array file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i < number_mcmc_walkers) THEN
            IF (mask_array(l, i) .EQV. .TRUE.) THEN
                WRITE (99, "(A4, 2X)", advance='no') 'True'
            ELSE IF (mask_array(l, i) .EQV. .FALSE.) THEN
                WRITE (99, "(A5, 1X)", advance='no') 'False'
            END IF
        ELSE IF (i == number_mcmc_walkers) THEN
            IF (mask_array(l, i) .EQV. .TRUE.) THEN
                WRITE (99, "(A4)", advance='no') 'True'
            ELSE IF (mask_array(l, i) .EQV. .FALSE.) THEN
                WRITE (99, "(A5)", advance='no') 'False'
            END IF
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write the randomly drawn array used for the spin-orbit angle.
OPEN(unit=99, FILE=output_spin_orbit_norm_array_filename, status='replace', action='write')
!Now write array line by line into chi squared array file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F30.15)", advance='no') normal14(l, 1)
        ELSE
            WRITE (99, "(1X, F30.15)", advance='no') normal14(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write the scale factor array.
OPEN(unit=99, FILE=output_scale_factor_array_filename, status='replace', action='write')
!Now write array line by line into chi squared array file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F30.15)", advance='no') scale_factor_array(l, 1)
        ELSE
            WRITE (99, "(1X, F30.15)", advance='no') scale_factor_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write output chi squared array so that IDL can read the array in. This is the array with only the chi squared values from the accepted proposals.
OPEN(unit=99, FILE=output_chi_squared_array_MCMC_filename, status='replace', action='write')
!Now write array line by line into chi squared array file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F30.15)", advance='no') Chi_squared_array_MCMC(l, 1)
        ELSE
            WRITE (99, "(1X, F30.15)", advance='no') Chi_squared_array_MCMC(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write output reduced chi squared array so that IDL can read the array in. This is the array with only the reduced chi squared values
!from the accepted proposals.
OPEN(unit=99, FILE=output_reduced_chi_array_MCMC_filename, status='replace', action='write')
!Now write array line by line into reduced chi squared array file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F20.10)", advance='no') Reduced_Chi_Squared_array_MCMC(l, 1)
        ELSE
            WRITE (99, "(1X, F20.10)", advance='no') Reduced_Chi_Squared_array_MCMC(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write output likelihood array so that IDL can read the array in. This is the array with only the likelihood values
!from the accepted proposals.
OPEN(unit=99, FILE=output_likelihood_array_MCMC_filename, status='replace', action='write')
!Now write array line by line into likelihood array file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F20.10)", advance='no') likelihood_array_MCMC(l, 1)
        ELSE
            WRITE (99, "(1X, F20.10)", advance='no') likelihood_array_MCMC(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write output file of the vsini array from accepted proposals.
OPEN(unit=99, FILE=output_vsini_MCMC_array_filename, status='replace', action='write')
!Now write vsini array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') vsini_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') vsini_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write output file of the vsini array from both accepted and rejected proposals.
OPEN(unit=99, FILE=output_vsini_all_array_filename, status='replace', action='write')
!Now write vsini array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') vsini_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') vsini_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)
 
!Write output file of the spin-orbit array from accepted proposals.
OPEN(unit=99, FILE=output_spin_orbit_MCMC_array_filename, status='replace', action='write')
!Now write the spin-orbit array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') stellar_rotation_angle_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') stellar_rotation_angle_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the spin-orbit array from accepted proposals.
OPEN(unit=99, FILE=output_spin_orbit_all_array_filename, status='replace', action='write')
!Now write the spin-orbit array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') stellar_rotation_angle_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') stellar_rotation_angle_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)      

!Write output file of the RV offset datasets array from accepted proposals.
OPEN(unit=99, FILE=output_RV_offset_datasets_MCMC_array_filename, status='replace', action='write')
!Now write RV offset datasets array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') RV_offset_datasets_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') RV_offset_datasets_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write output file of the RV offset datasets array from both accepted and rejected proposals.
OPEN(unit=99, FILE=output_RV_offset_datasets_all_array_filename, status='replace', action='write')
!Now write RV offset datasets array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') RV_offset_datasets_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') RV_offset_datasets_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)
 
!Write output file of the RV zero offset array from accepted proposals.
OPEN(unit=99, FILE=output_RV_zero_offset_MCMC_array_filename, status='replace', action='write')
!Now write the RV zero offset array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') RV_zero_offset_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') RV_zero_offset_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the RV zero offset array from accepted proposals.
OPEN(unit=99, FILE=output_RV_zero_offset_all_array_filename, status='replace', action='write')
!Now write the RV zero offset array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') RV_zero_offset_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') RV_zero_offset_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)        

!Write output file of the orbital period array from accepted proposals.
OPEN(unit=99, FILE=output_Orbital_period_MCMC_array_filename, status='replace', action='write')
!Now write the orbital period array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F20.10)", advance='no') Orbital_period_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F20.10)", advance='no') Orbital_period_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the orbital period array from accepted proposals.
OPEN(unit=99, FILE=output_Orbital_period_all_array_filename, status='replace', action='write')
!Now write the orbital period array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F20.10)", advance='no') Orbital_period_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F20.10)", advance='no') Orbital_period_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)   

!Write output file of the JD mid transit array from accepted proposals.
OPEN(unit=99, FILE=output_JD_time_mid_transit_MCMC_array_filename, status='replace', action='write')
!Now write the JD mid transit array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F20.10)", advance='no') JD_time_mid_transit_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F20.10)", advance='no') JD_time_mid_transit_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the JD mid transit array from accepted proposals.
OPEN(unit=99, FILE=output_JD_time_mid_transit_all_array_filename, status='replace', action='write')
!Now write the JD mid transit array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F20.10)", advance='no') JD_time_mid_transit_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F20.10)", advance='no') JD_time_mid_transit_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)           

!Write output file of the Mp array from accepted proposals.
OPEN(unit=99, FILE=output_Mp_MCMC_array_filename, status='replace', action='write')
!Now write the Mp array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Mp_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Mp_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the Mp array from accepted proposals.
OPEN(unit=99, FILE=output_Mp_all_array_filename, status='replace', action='write')
!Now write the Mp array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Mp_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Mp_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)      
   
!Write output file of the Rp array from accepted proposals.
OPEN(unit=99, FILE=output_Rp_MCMC_array_filename, status='replace', action='write')
!Now write the Rp array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Rp_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Rp_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the Rp array from accepted proposals.
OPEN(unit=99, FILE=output_Rp_all_array_filename, status='replace', action='write')
!Now write the Rp array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Rp_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Rp_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)  
   
!Write output file of the Ms_solar array from accepted proposals.
OPEN(unit=99, FILE=output_Ms_solar_MCMC_array_filename, status='replace', action='write')
!Now write the Ms_solar array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Ms_solar_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Ms_solar_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the Ms_solar array from accepted proposals.
OPEN(unit=99, FILE=output_Ms_solar_all_array_filename, status='replace', action='write')
!Now write the Ms_solar array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Ms_solar_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Ms_solar_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)
     
!Write output file of the Rs_solar array from accepted proposals.
OPEN(unit=99, FILE=output_Rs_solar_MCMC_array_filename, status='replace', action='write')
!Now write the Rs_solar array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Rs_solar_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Rs_solar_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the Rs_solar array from accepted proposals.
OPEN(unit=99, FILE=output_Rs_solar_all_array_filename, status='replace', action='write')
!Now write the Rs_solar array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Rs_solar_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Rs_solar_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)       
   
!Write output file of the Rp_Rs_ratio array from accepted proposals.
OPEN(unit=99, FILE=output_Rp_Rs_ratio_MCMC_array_filename, status='replace', action='write')
!Now write the Rp_Rs_ratio array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Rp_Rs_ratio_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Rp_Rs_ratio_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the Rp_Rs_ratio array from accepted proposals.
OPEN(unit=99, FILE=output_Rp_Rs_ratio_all_array_filename, status='replace', action='write')
!Now write the Rp_Rs_ratio array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Rp_Rs_ratio_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Rp_Rs_ratio_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)

!Write output file of the Ecc array from accepted proposals.
OPEN(unit=99, FILE=output_Ecc_MCMC_array_filename, status='replace', action='write')
!Now write the Ecc array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Ecc_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Ecc_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the Ecc array from accepted proposals.
OPEN(unit=99, FILE=output_Ecc_all_array_filename, status='replace', action='write')
!Now write the Ecc array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Ecc_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Ecc_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)
      
!Write output file of the Inc array from accepted proposals.
OPEN(unit=99, FILE=output_Inc_MCMC_array_filename, status='replace', action='write')
!Now write the Inc array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Inc_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Inc_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the Inc array from accepted proposals.
OPEN(unit=99, FILE=output_Inc_all_array_filename, status='replace', action='write')
!Now write the Inc array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') Inc_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') Inc_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)
   
!Write output file of the omega_arg_periastron array from accepted proposals.
OPEN(unit=99, FILE=output_omega_arg_periastron_MCMC_array_filename, status='replace', action='write')
!Now write the omega_arg_periastron array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') omega_arg_periastron_accept_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') omega_arg_periastron_accept_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99) 
            
!Write output file of the omega_arg_periastron array from accepted proposals.
OPEN(unit=99, FILE=output_omega_arg_periastron_all_array_filename, status='replace', action='write')
!Now write the omega_arg_periastron array line by line into a file.
DO l = 1, NINT(mcmc_accepted_iteration_size/0.05)
    DO i = 1, number_mcmc_walkers
        IF (i == 1) THEN
            WRITE (99, "(F15.6)", advance='no') omega_arg_periastron_all_mcmc_array(l, 1)
        ELSE
            WRITE (99, "(1X, F15.6)", advance='no') omega_arg_periastron_all_mcmc_array(l, i)
        END IF        
    END DO
    WRITE (99, "(A1)", advance='yes') blank
END DO 
CLOSE(99)




!Output the MCMC run properties such as the number proposals rejected, the total number of MCMC iterations,
!number of accepted iterations, acceptance rate, and scale factor.
OPEN(unit=99, FILE=output_MCMC_properties_filename, status='replace', action='write')
DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F15.6)", advance='no') walker_parameter_array(1, 1)
    ELSE
        WRITE (99, "(1X, F15.6)", advance='no') walker_parameter_array(1, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F15.6)", advance='no') walker_parameter_array(2, 1)
    ELSE
        WRITE (99, "(1X, F15.6)", advance='no') walker_parameter_array(2, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F15.6)", advance='no') walker_parameter_array(3, 1)
    ELSE
        WRITE (99, "(1X, F15.6)", advance='no') walker_parameter_array(3, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F15.6)", advance='no') walker_parameter_array(4, 1)
    ELSE
        WRITE (99, "(1X, F15.6)", advance='no') walker_parameter_array(4, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F15.6)", advance='no') walker_parameter_array(5, 1)
    ELSE
        WRITE (99, "(1X, F15.6)", advance='no') walker_parameter_array(5, i)
    END IF        
END DO
CLOSE(99)

 


!Output the best fit parameters for each walker based on the minimum chi squared.
OPEN(unit=99, FILE=output_best_param_chi_squared_filename, status='replace', action='write')
DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(1, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(1, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(2, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(2, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(3, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(3, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(4, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(4, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(5, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(5, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(6, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(6, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(7, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(7, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(8, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(8, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(9, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(9, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(10, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(10, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(11, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(11, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(12, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(12, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(13, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(13, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_chi2(14, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_chi2(14, i)
    END IF        
END DO
CLOSE(99)


 

!Output the best fit parameters for each walker based on the mean of all MCMC iterations.
OPEN(unit=99, FILE=output_best_param_mean_filename, status='replace', action='write')
DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(1, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(1, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(2, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(2, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(3, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(3, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(4, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(4, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(5, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(5, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(6, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(6, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(7, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(7, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(8, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(8, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(9, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(9, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(10, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(10, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(11, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(11, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(12, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(12, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(13, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(13, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_mean(14, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_mean(14, i)
    END IF        
END DO
CLOSE(99)




!Output the standard deviation of each parameters for each walker based on all of MCMC iterations.
OPEN(unit=99, FILE=output_stand_dev_filename, status='replace', action='write')
DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(1, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(1, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(2, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(2, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(3, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(3, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(4, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(4, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(5, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(5, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(6, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(6, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(7, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(7, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(8, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(8, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(9, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(9, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(10, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(10, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(11, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(11, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(12, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(12, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(13, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(13, i)
    END IF        
END DO
WRITE (99, "(A1)", advance='yes') blank

DO i = 1, number_mcmc_walkers
    IF (i == 1) THEN
        WRITE (99, "(F16.8)", advance='no') walker_parameter_array_stand_dev(14, 1)
    ELSE
        WRITE (99, "(1X, F16.8)", advance='no') walker_parameter_array_stand_dev(14, i)
    END IF        
END DO
CLOSE(99)




!Write out best vsini, best spin_orbit, minimum chi squared, minimum reduced chi squared, 
!location of minimum chi squared, uncertainty for vsini and spin_orbit, and locations for 
!these uncertainties.
OPEN(unit=99, FILE=output_best_parameters_mcmc_filename, status='replace', action='write')
WRITE(99,133) 'Minimum Chi Squared all walkers:', min_chi_squared_all_walkers
WRITE(99,134) 'Location of Minimum Chi Squared all walkers:', loc_min_chi_squared_all_walkers(1), ',', loc_min_chi_squared_all_walkers(2)
WRITE(99,135) 'Minimum Reduced Chi Squared all walkers:', min_reduced_chi_squared_all_walkers
WRITE(99,136) 'Location of Minimum Reduced Chi Squared all walkers:', loc_min_reduced_chi_squared_all_walkers(1), ',', loc_min_reduced_chi_squared_all_walkers(2)
WRITE(99,137) 'Maximum likelihood all walkers:', max_likelihood_all_walkers
WRITE(99,138) 'Location of maximum likelihood all walkers:', loc_max_likelihood_all_walkers(1), ',', loc_max_likelihood_all_walkers(2)

WRITE(99,139) 'Best vsin(i) value from minimum chi square all walkers:', best_vsini_chi2_all_walkers, 'km/s'
WRITE(99,140) 'Best vsin(i) value from the mean of MCMC runs all walkers:', vsini_mean2, 'km/s'
WRITE(99,141) 'Uncertainty on vsin(i) value from variance of mean all walkers:', vsini_stand_dev, 'km/s'

WRITE(99,143) 'Best spin-orbit alignment angle value from minimum chi square all walkers:', best_spin_orbit_chi2_all_walkers, 'deg'
WRITE(99,144) 'Best spin-orbit alignment angle value from mean of MCMC runs all walkers:', stellar_rotation_angle_mean2, 'deg'
WRITE(99,145) 'Uncertainty on spin-orbit angle value from variance of mean all walkers:', stellar_rotation_angle_stand_dev, 'deg'

WRITE(99,147) 'Best RV offset between datasets:', RV_offset_datasets_mean2, 'm/s'
WRITE(99,148) 'Variance of RV offset between datasets:', RV_offset_datasets_stand_dev, 'm/s'

WRITE(99,149) 'Best RV zero offset:', RV_zero_offset_mean2, 'm/s'
WRITE(99,150) 'Variance of RV zero offset:', RV_zero_offset_stand_dev, 'm/s'

WRITE(99,151) 'Best orbital period:', Orbital_period_mean2, 'days'
WRITE(99,152) 'Variance of orbital period:', Orbital_period_stand_dev, 'days'

WRITE(99,153) 'Best JD time mid transit:', JD_time_mid_transit_mean2, 'days'
WRITE(99,154) 'Variance of JD time mid transit:', JD_time_mid_transit_stand_dev, 'days'

WRITE(99,155) 'Best Mp:', Mp_mean2, 'Mj'
WRITE(99,156) 'Variance of Mp:', Mp_stand_dev, 'Mj'

WRITE(99,157) 'Best Rp:', Rp_mean2, 'Rj'
WRITE(99,158) 'Variance of Rp:', Rp_stand_dev, 'Rj'

WRITE(99,159) 'Best Ms solar:', Ms_solar_mean2, 'Ms'
WRITE(99,160) 'Variance of Ms solar:', Ms_solar_stand_dev, 'Ms'

WRITE(99,161) 'Best Rs solar:', Rs_solar_mean2, 'Rs'
WRITE(99,162) 'Variance of Rs solar:', Rs_solar_stand_dev, 'Rs'

WRITE(99,163) 'Best Rp to Rs ratio:', Rp_Rs_ratio_mean2
WRITE(99,164) 'Variance of Rp to Rs ratio:', Rp_Rs_ratio_stand_dev

WRITE(99,165) 'Best Ecc:', Ecc_mean2
WRITE(99,166) 'Variance of Ecc:', Ecc_stand_dev

WRITE(99,167) 'Best Inc:', Inc_mean2, 'deg'
WRITE(99,168) 'Variance of Inc:', Inc_stand_dev, 'deg'

WRITE(99,169) 'Best omega arg periastron:', omega_arg_periastron_mean2, 'deg'
WRITE(99,170) 'Variance of omega arg periastron:', omega_arg_periastron_stand_dev, 'deg'

133  FORMAT(A90, 2X, F15.6)
134  FORMAT(A90, 2X, I10, A1, I10)
135  FORMAT(A90, 2X, F15.6)
136  FORMAT(A90, 2X, I10, A1, I10)
137  FORMAT(A90, 2X, F15.6)
138  FORMAT(A90, 2X, I10, A1, I10)

139  FORMAT(A90, 2X, F15.6, 5X, A4)
140  FORMAT(A90, 2X, F15.6, 5X, A4)
141  FORMAT(A90, 2X, F15.6, 5X, A4)

143  FORMAT(A90, 2X, F15.6, 5X, A3)
144  FORMAT(A90, 2X, F15.6, 5X, A3)
145  FORMAT(A90, 2X, F15.6, 5X, A3)

147  FORMAT(A90, 2X, F15.6, 5X, A3)
148  FORMAT(A90, 2X, F15.6, 5X, A3)

149  FORMAT(A90, 2X, F15.6, 5X, A3)
150  FORMAT(A90, 2X, F15.6, 5X, A3)

151  FORMAT(A90, 2X, F20.10, 5X, A4)
152  FORMAT(A90, 2X, F20.10, 5X, A4)

153  FORMAT(A90, 2X, F20.10, 5X, A4)
154  FORMAT(A90, 2X, F20.10, 5X, A4)

155  FORMAT(A90, 2X, F15.6, 5X, A2)
156  FORMAT(A90, 2X, F15.6, 5X, A2)

157  FORMAT(A90, 2X, F15.6, 5X, A2)
158  FORMAT(A90, 2X, F15.6, 5X, A2)

159  FORMAT(A90, 2X, F15.6, 5X, A2)
160  FORMAT(A90, 2X, F15.6, 5X, A2)

161  FORMAT(A90, 2X, F15.6, 5X, A2)
162  FORMAT(A90, 2X, F15.6, 5X, A2)

163  FORMAT(A90, 2X, F20.10)
164  FORMAT(A90, 2X, F20.10)

165  FORMAT(A90, 2X, F20.10)
166  FORMAT(A90, 2X, F20.10)

167  FORMAT(A90, 2X, F15.6, 5X, A3)
168  FORMAT(A90, 2X, F15.6, 5X, A3)

169  FORMAT(A90, 2X, F15.6, 5X, A3)
170  FORMAT(A90, 2X, F15.6, 5X, A3)
CLOSE(99)


   

call system_clock(system_time_4)
PRINT *, 'The time taken to run is: ', (system_time_4 - system_time_1)/clock_rate, ' seconds'

DEALLOCATE(Data_my_rv)
DEALLOCATE(Data_other_rv)
DEALLOCATE(RV_offset_datasets_accept_mcmc_array)
DEALLOCATE(RV_offset_datasets_all_mcmc_array)
DEALLOCATE(RV_offset_datasets_sum2)
DEALLOCATE(Orbital_period_accept_mcmc_array)
DEALLOCATE(Orbital_period_all_mcmc_array)
DEALLOCATE(Orbital_period_sum2)
DEALLOCATE(JD_time_mid_transit_accept_mcmc_array)
DEALLOCATE(JD_time_mid_transit_all_mcmc_array)
DEALLOCATE(JD_time_mid_transit_sum2)
DEALLOCATE(Mp_accept_mcmc_array)
DEALLOCATE(Mp_all_mcmc_array)
DEALLOCATE(Mp_sum2)
DEALLOCATE(Rp_accept_mcmc_array)
DEALLOCATE(Rp_all_mcmc_array)
DEALLOCATE(Rp_sum2)
DEALLOCATE(Rp_Rs_ratio_accept_mcmc_array)
DEALLOCATE(Rp_Rs_ratio_all_mcmc_array)
DEALLOCATE(Rp_Rs_ratio_sum2)
DEALLOCATE(Ecc_accept_mcmc_array)
DEALLOCATE(Ecc_all_mcmc_array)
DEALLOCATE(Ecc_sum2)
DEALLOCATE(Inc_accept_mcmc_array)
DEALLOCATE(Inc_all_mcmc_array)
DEALLOCATE(Inc_sum2)
DEALLOCATE(omega_arg_periastron_accept_mcmc_array)
DEALLOCATE(omega_arg_periastron_all_mcmc_array)
DEALLOCATE(omega_arg_periastron_sum2)
DEALLOCATE(RV_zero_offset_accept_mcmc_array)
DEALLOCATE(RV_zero_offset_all_mcmc_array)
DEALLOCATE(RV_zero_offset_sum2)
DEALLOCATE(Rs_solar_accept_mcmc_array)
DEALLOCATE(Rs_solar_all_mcmc_array)
DEALLOCATE(Rs_solar_sum2)
DEALLOCATE(Ms_solar_accept_mcmc_array)
DEALLOCATE(Ms_solar_all_mcmc_array)
DEALLOCATE(Ms_solar_sum2)
DEALLOCATE(vsini_accept_mcmc_array)
DEALLOCATE(vsini_all_mcmc_array)
DEALLOCATE(vsini_sum2)
DEALLOCATE(vsini_variance_current_array)
DEALLOCATE(vsini_variance_MCMC_array)
DEALLOCATE(vsini_stand_dev_current_array)
DEALLOCATE(vsini_stand_dev_MCMC_array)
DEALLOCATE(stellar_rotation_angle_accept_mcmc_array)
DEALLOCATE(stellar_rotation_angle_all_mcmc_array)
DEALLOCATE(stellar_rotation_angle_sum2)
DEALLOCATE(stellar_rotation_angle_variance_current_array)
DEALLOCATE(stellar_rotation_angle_variance_MCMC_array)
DEALLOCATE(stellar_rotation_angle_stand_dev_current_array)
DEALLOCATE(stellar_rotation_angle_stand_dev_MCMC_array)
DEALLOCATE(Chi_squared_array_MCMC)
DEALLOCATE(Reduced_Chi_Squared_array_MCMC)
DEALLOCATE(likelihood_array_MCMC)
DEALLOCATE(culumative_likelihood_array_MCMC)
DEALLOCATE(index_array)
DEALLOCATE(walker_parameter_array)
DEALLOCATE(loc_min_chi_squared_MCMC)
DEALLOCATE(loc_min_reduced_chi_squared_MCMC)
DEALLOCATE(loc_max_likelihood_MCMC)
DEALLOCATE(walker_parameter_array_chi2)
DEALLOCATE(walker_parameter_array_mean)
DEALLOCATE(walker_parameter_array_variance)
DEALLOCATE(walker_parameter_array_stand_dev)
DEALLOCATE(min_chi_squared_MCMC)
DEALLOCATE(min_reduced_chi_squared_MCMC)
DEALLOCATE(max_likelihood_MCMC)
DEALLOCATE(mask_array)
DEALLOCATE(normal1)
DEALLOCATE(normal2)
DEALLOCATE(normal3)
DEALLOCATE(normal4)
DEALLOCATE(normal5)
DEALLOCATE(normal6)
DEALLOCATE(normal7)
DEALLOCATE(normal8)
DEALLOCATE(normal9)
DEALLOCATE(normal10)
DEALLOCATE(normal11)
DEALLOCATE(normal12)
DEALLOCATE(normal13)
DEALLOCATE(normal14)
DEALLOCATE(random_draw_array)
DEALLOCATE(spread_rand_RV_offset)
DEALLOCATE(spread_rand_Orbital_period)
DEALLOCATE(spread_rand_JD)
DEALLOCATE(spread_rand_Mp)
DEALLOCATE(spread_rand_Rp)
DEALLOCATE(spread_rand_Rp_Rs_ratio)
DEALLOCATE(spread_rand_Rs)
DEALLOCATE(spread_rand_Ms)
DEALLOCATE(spread_rand_Ecc)
DEALLOCATE(spread_rand_Inc)
DEALLOCATE(spread_rand_omega)
DEALLOCATE(spread_rand_RV_zero)
DEALLOCATE(spread_rand_vsini)
DEALLOCATE(spread_rand_spin_orbit)
DEALLOCATE(scale_factor_array)




cmd1 = "say 'Finished processing data.'"
cmd1 = TRIM(cmd1)
call execute_command_line(cmd1, wait=.false., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)

CONTAINS

    !Subroutine for ranking an array in ascending order. 
    PURE SUBROUTINE Selection_sort(a,index_array)
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

END PROGRAM ExOSAM_RM_effect_V02
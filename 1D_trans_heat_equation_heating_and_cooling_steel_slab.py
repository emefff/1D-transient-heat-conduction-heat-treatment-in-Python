"""
The below script solves the transient 1D heat equation with nonlinear HTC, 
thermal conductivity and heat capacity. 

It is based on https://github.com/alexlib/trans_heat_cond by Alex Liberzon
but modified for nonlinear coefficients which are preferred for heat treatment
 of steel etc. The reason being that temperature differences are quite high. 
 Only the 'slab' option from the original library is used.
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS IS STILL WORK IN PROGRESS, USE AT YOUR OWN RISK.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

"""

import numpy as np
import matplotlib.pyplot as plt
from funcTheta import theta
import tqdm


###############################################################################
def interpolate(x_list, y_list, x_eval):
    """
    Returns a linearly interpolated value of some x_eval. X- and Y- coordinates
    of a function to interpolate are given by x_list and y_list.

    Parameters
    ----------
    x_list : TYPE list of floats
        DESCRIPTION. List of x-coords of funtion to be interpolated
    y_list : TYPE list of floats
        DESCRIPTION. List of y-coords of funtion to be interpolated
    x_eval : TYPE float
        DESCRIPTION. x-coord to be evaluated. Must be within the boundaries
                     of x_list, that is it cannot be smaller than min(x_list)
                     and larger than max(x_list). THIS IS NOT CHECKED!!

    Returns 
    -------
    y_val : TYPE float
    DESCRIPTION. y-coord evaluated via linear interpolation.

    """
    # printing = True
    printing = False

    x_numpy = np.array(x_list)
    y_numpy = np.array(y_list)
    x_nearest = min(x_list, key=lambda x:abs(x-x_eval))
    
    # we need the nearest larger element
    x_larger = x_numpy[x_numpy > x_eval].min()  
    x_larger_idx = list(x_numpy).index(x_larger)

    # we also need the nearest lower element
    x_smaller = x_numpy[x_numpy < x_eval].max()
    x_smaller_idx = list(x_numpy).index(x_smaller)
    # what if x_eval is exactly on one point?? then x_smaller is wrong!
    if (x_eval == x_nearest):
        x_smaller = x_eval
        x_smaller_idx = list(x_numpy).index(x_smaller)
        
    # we need the y-values of our found x_larger and x_smaller
    y_smaller = y_numpy[x_smaller_idx]
    y_larger = y_numpy[x_larger_idx]
    
    # now we calculate the y_eval at our x_eval via simple interpolation
    k = (y_larger - y_smaller) / (x_larger - x_smaller)
    y_eval = k * (x_eval - x_smaller) + (y_smaller)
    
    if printing == True:
        print("###########################")
        print(f"{x_eval = }, {type(x_eval) = }")
        print(f"{x_nearest = }")
        print(f"{x_smaller = }, {x_smaller_idx = }, {y_smaller = }")
        print(f"{x_larger = }, {x_larger_idx = }, {y_larger = }")    
        print(f"{y_eval = }\n")
        print("***************************\n")
    return y_eval

def flatten(xss):
    return [x for xs in xss for x in xs]
###############################################################################


# plt.close('all')

# PARAMETERS AND INPUT DATA
rhow = 7850     # density of steel, 7850 kg/m3
d = 0.25E0      # sphere/slab/cylinder diameter, m
# cpw = 466       # specific heat capacity steel, J/(kg K)
# kw = 31.8       # thermal conductivity steel, W/(m K)
Ti = 25+273.15  # uniform initial temp of slab, K
Tinf = 870+273  # surrounding fluid or gas temp, K
t_heating = 8*3600   # heating time, s
t_cooling = 8*3600   # cooling time in s
tmax = t_heating + t_cooling # total time
# h = 375         # heat transfer coefficent, W/(m2 K)

# SOME INITIAL CALCULATIONS
ro = (d/2)                          # radius of sphere (a.k.a outer radius), m
rs = ro/ro                          # dimensionless surface radius, (-)
rc = 1e-12/ro                       # dimensionless center radius, (-)

num_steps_heating = 5000
step_size_heating = t_heating/num_steps_heating
# alpha = kw/(rhow*cpw)               # thermal diffusivity steel, m^2/s
t = np.arange(0, t_heating, step_size_heating) # time range for simulation, s
z = np.arange(0, 1250, 0.1)         # range to evaluate the zeta, Bi equation
z[0] = 1e-12                        # prevent divide by zero warning

# Bi = (h*ro)/kw                      # Biot number, (-)
# Fo = (alpha * t) / (ro**2)          # Fourier number, (-)

# SLAB TEMPERATURE PROFILE
b = 0   # shape factor where 2 sphere, 1 cylinder, 0 slab

h_list_temperature = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900]
h_list_values_heating = [102, 111, 126, 143, 163, 186, 212, 243, 277, 297]
h_list_values_cooling = [40, 40, 60, 80, 100, 120, 140, 160, 180, 200]
# h_list_values_cooling = h_list_values_heating
# h_list_values_cooling = [300, 300, 300, 300, 300, 300, 300, 300, 300, 300]
# h_list_values_heating = [500, 500, 500, 500, 500, 500, 500, 500, 500, 500]
# h_list_values_heating = h_list_values_radiation

print("HEAT TRANSFER COEFFICIENT vs TEMPERATURE")
print(f"{h_list_temperature = }")
print(f"{h_list_values_heating = }")
print(f"{h_list_values_cooling = }")
# plt.figure(figsize=(15,9))
# plt.plot(h_list_temperature, h_list_values_heating, 'r-', label="heating")
# plt.plot(h_list_temperature, h_list_values_cooling, 'b-', label="cooling")
# plt.title("Heat transfer coefficients")
# plt.ylabel('HTC / W/(m²*K)')
# plt.xlabel('Temperature / °C')
# plt.show()

# we also vary the thermal conductivity of the steel
# values taken from https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783433601570.app1 
# for carbon steel

kw_list_temperature = [20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900]
kw_list_values = []
for i, temp in enumerate(kw_list_temperature): # having these formulae to calc the data makes the following interpolation
    if 20 <= temp < 800:                       # obsolete, but we may not always have such good data
        kw = 54 -3.33E-2 * temp                # mostly we will have lists.
    elif 800 <= temp < 1200:
        kw = 27.3
    else:
        kw = None
        print("For this temperature {temp} we do not have kw-values!")
    kw_list_values.append(kw)
print("\nTHERMAL CONDUCTIVITY vs TEMPERATURE")
print(f"{kw_list_temperature = }")
print(f"{kw_list_values = }")
# plt.figure(figsize=(15,9))
# plt.plot(kw_list_temperature, kw_list_values, 'k-')
# plt.title("Thermal conductivity")
# plt.ylabel('Thermal conductivity / W/(m*K)')
# plt.xlabel('Temperature / °C')
# plt.show()

# we vary the heat capacity of the steel
# values taken from https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783433601570.app1 
# for carbon steel
cpw_list_temperature = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900] # in these data, it is better to have less temps
cpw_list_values = []                                                     # the curve has a sharp spike, if we hit it we'll
for i, temp in enumerate(cpw_list_temperature):                          # get a lot of oscillations in the result
    if 0<temp<600:
        cpw = 425 + 7.73E-1*temp - 1.69E-3*temp**2 + 2.22E-6*temp**3
    elif 600<=temp<735:
        cpw = 666 + 13002/(738-temp)
    elif 735<=temp<900:
        cpw = 545 + 17820/(temp-731)
    elif 900<=temp<=1200:
        cpw = 650
    else:
        cpw = None
        print("For this temperature {temp} we do not have cpw-values!")
    cpw_list_values.append(cpw)
print(f"\nHEAT CAPACITY vs TEMPERATURE")
print(f"{cpw_list_temperature = }")
print(f"{cpw_list_values = }")
# plt.figure(figsize=(15,9))
# plt.plot(cpw_list_temperature, cpw_list_values, 'k-')
# plt.title(f"Heat capacity")
# plt.ylabel('Heat capacity / J/(kg*K)')
# plt.xlabel('Temperature / °C')
# plt.show()

# we want to calculate with nonlinear h, and alpha-values thus we split the calculation
# time interval according to the number of intervals. We cannot split the temperature
# intervals, because it is unknown before the calculation. BUT: Splitting the time and taking
# the h from the h_list values is better than using on single, constant h.

num_intervals = 100
num_time_intervals_heating = num_intervals # the larger the smoother the curves but can lead to bifurcating oscillations if too small
# must be a divisor of step_size_cooling
time_interval_heating = t_heating / (num_time_intervals_heating)

# we need some lists for the split results
t_list = []
thetaRo_list = []
To_slab_list = []
thetaR_list = []
Tr_slab_list = []

# just for evalution and debugging
debugging = False
if debugging == True:
    h_list1 = []
    h_list2 = []
    kw_list1 = []
    kw_list2 = []
    cpw_list1 = []
    cpw_list2 = []
    alpha_list1 = []
    alpha_list2 = []
    T_cent_init_list = []

t_init = 0
t_end = time_interval_heating # for the first iteration

T_surf_init = Ti - 273
T_cent_init = T_surf_init
T_avg = T_surf_init

# HEATING
print("\nCalculating heating step...")
for i in tqdm.tqdm(range(num_time_intervals_heating)):
    # rs, rc, b, z are constant, the other parameters t, Bi , Fo must be adapted
    # per timestep!!
    
    # this is our new split time interval
    t = np.arange(t_init, t_end, step_size_heating) # time range for simulation, s there was a +0.002 in t_end
    
    # we must compute a new h for Bi and for the surface, we have to separate the temperatures also
    h_new = interpolate(h_list_temperature, h_list_values_heating, T_surf_init)
    # h_list1.append(h_new)
    Bi = (h_new*ro)/kw                      # Biot number, (-)
    
    # we also adapt alpha in the same way as h!
    kw_new = interpolate(kw_list_temperature, kw_list_values, T_surf_init)
    cpw_new = interpolate(cpw_list_temperature, cpw_list_values, T_surf_init)
    alpha_new = kw_new / (rhow*cpw_new) # if you want to keep it constant just do kw/( rhwo * cpw) here!
    
    # for debugging
    # kw_list1.append(kw_new)
    # cpw_list1.append(cpw_new)
    # alpha_list1.append(alpha_new)
    # alpha_new = kw / (rhow * cpw)
    
    Fo = (alpha_new * t) / (ro**2)   
    
    # surface temperature where ro for outer surface
    thetaRo = theta(rs, b, z, Bi, Fo)   # dimensionless temperature profile
    To_slab = Tinf + thetaRo*(Ti-Tinf)  # convert theta to temperature in Kelvin, K
    
    # we must compute a new h for Bi and for the center, we have to separate the temperatures also
    h_new = interpolate(h_list_temperature, h_list_values_heating, T_cent_init)
    # h_list2.append(h_new)
    # T_cent_init_list.append(T_cent_init)
    Bi = (h_new*ro)/kw                      # Biot number, (-)
    
    # we also adapt alpha in the same way as h!
    kw_new = interpolate(kw_list_temperature, kw_list_values, T_cent_init)
    cpw_new = interpolate(cpw_list_temperature, cpw_list_values, T_cent_init)
    alpha_new = kw_new / (rhow*cpw_new) # if you want to keep it constant just do kw/( rhwo * cpw) here!
    
    # for debugging
    # kw_list2.append(kw_new)
    # cpw_list2.append(cpw_new)
    # alpha_list2.append(alpha_new)
    # alpha_new = kw / (rhow * cpw)
    
    Fo = (alpha_new * t) / (ro**2)   
        
    # center temperature where r for center
    thetaR = theta(rc, b, z, Bi, Fo)    # dimensionless temperature profile
    Tr_slab = Tinf + thetaR*(Ti-Tinf)   # convert theta to temperature in Kelvin, K
    
    # new values for next iteration
    t_init = t_end
    t_end = t_init + time_interval_heating
    T_surf_init = To_slab[-1] - 273.15 # we need celcius for the list, not pretty :(
    T_cent_init = Tr_slab[-1] - 273.15 # we need celcius for the list, not pretty :(
    T_surf_avg = np.average(To_slab) - 273.15
    T_cent_avg = np.average(Tr_slab) - 273.15
    
    T_surf_init = T_surf_avg
    T_cent_init = T_cent_avg
    
    # we need to store all results in some lists
    t_list.append(list(t))
    thetaRo_list.append(list(thetaRo))
    To_slab_list.append(list(To_slab))
    thetaR_list.append(list(thetaR))
    Tr_slab_list.append(list(Tr_slab))
  

# COOLING
print("\nCalculating cooling step...")
num_time_intervals_cooling = num_intervals # must be a divisor of step_size_cooling
time_interval_cooling = t_cooling / num_time_intervals_cooling

Ti_new = To_slab_list[-1][-1] # 850+273.15  # uniform initial temp of slab, K
Tinf_new = 25+273  # surrounding fluid or gas temp, K

t_init = 0
t_end = t_init + time_interval_cooling # for the first iteration

t_init_new = t_heating
t_end_new = t_init_new + time_interval_cooling # for the first iteration

T_surf_init_new = Ti_new - 273
T_cent_init_new = Tr_slab_list[-1][-1] - 273.15

num_steps_cooling = 4000
step_size_cooling = t_cooling/num_steps_cooling

for i in tqdm.tqdm(range(num_time_intervals_cooling)):
    # rs, rc, b, z are constant, the other parameters t, Bi , Fo must be adapted
    # per timestep!!
    
    # this is our new split time interval
    t = np.arange(t_init, t_end, step_size_cooling) # time range for simulation, s there was a +0.002 in t_end
    t_new = np.arange(t_init_new, t_end_new, step_size_cooling) # only for continuous time in final t-list, otherwise Fo is completely wrong!
    
    # we must compute a new h for Bi for the surface
    h_new = interpolate(h_list_temperature, h_list_values_cooling, T_surf_init_new)
    Bi = (h_new*ro)/kw                      # Biot number, (-)
    
    # we also adapt alpha in the same way as h!
    kw_new = interpolate(kw_list_temperature, kw_list_values, T_surf_init_new)
    cpw_new = interpolate(cpw_list_temperature, cpw_list_values, T_surf_init_new)
    alpha_new = kw_new / (rhow*cpw_new) # if you want to keep it constant just do kw/( rhwo * cpw) here!
       
    # alpha_new = kw / (rhow * cpw)
    # print(f"{alpha_new = }")
    Fo = (alpha_new * t) / (ro**2)   
    
    # surface temperature where ro for outer surface
    thetaRo = theta(rs, b, z, Bi, Fo)   # dimensionless temperature profile
    To_slab = Tinf_new + thetaRo*(Ti_new-Tinf_new)  # convert theta to temperature in Kelvin, K
    
    # we must compute a new h for Bi for the center
    h_new = interpolate(h_list_temperature, h_list_values_cooling, T_cent_init_new)
    Bi = (h_new*ro)/kw                      # Biot number, (-)
    
    # we also adapt alpha in the same way as h!
    kw_new = interpolate(kw_list_temperature, kw_list_values, T_cent_init_new)
    cpw_new = interpolate(cpw_list_temperature, cpw_list_values, T_cent_init_new)
    alpha_new = kw_new / (rhow*cpw_new) # if you want to keep it constant just do kw/( rhwo * cpw) here!
    # alpha_new = kw / (rhow * cpw)
    # print(f"{alpha_new = }")
    Fo = (alpha_new * t) / (ro**2)   

    # center temperature where r for center
    thetaR = theta(rc, b, z, Bi, Fo)    # dimensionless temperature profile
    Tr_slab = Tinf_new + thetaR*(Ti_new-Tinf_new)   # convert theta to temperature in Kelvin, K
    
    # new values for next iteration
    t_init = t_end
    t_end = t_init + time_interval_cooling
    
    # new values for continuous time
    t_init_new = t_end_new
    t_end_new = t_init_new + time_interval_cooling

    T_surf_init_new = To_slab[-1] - 273.15 # we need celcius for the list, not pretty :(
    T_cent_init_new = Tr_slab[-1] - 273.15 # we need celcius for the list, not pretty :(
    T_surf_avg = np.average(To_slab) - 273.15
    T_cent_avg = np.average(Tr_slab) - 273.15
    
    T_surf_init = T_surf_avg
    T_cent_init = T_cent_avg

    # we need to store all results in some lists
    t_list.append(list(t_new))
    thetaRo_list.append(list(thetaRo))
    To_slab_list.append(list(To_slab))
    thetaR_list.append(list(thetaR))
    Tr_slab_list.append(list(Tr_slab))
    

# PLOT RESULTS
# we need to flatten the lists of lists
t_list = flatten(t_list)
thetaRo_list = flatten(thetaR_list)
To_slab_list = flatten(To_slab_list)
thetaR_list = flatten(thetaR_list)
Tr_slab_list = flatten(Tr_slab_list)
  
# we convert the temperatures to °C, as engineers we do not have much use for °K
for i, temp in enumerate(To_slab_list):
    temp_new = temp - 273.15
    To_slab_list[i] = temp_new
for i, temp in enumerate(Tr_slab_list):
    temp_new = temp - 273.15
    Tr_slab_list[i] = temp_new
Th = Tinf - 273.15

# filtering oscillating results is not very elegant, but....
from scipy.signal import savgol_filter
window_size = num_time_intervals_cooling * 3
polyorder = 3
To_slab_list_filter = savgol_filter(To_slab_list, window_size, polyorder)
Tr_slab_list_filter = savgol_filter(Tr_slab_list, window_size, polyorder)
 
# in heat-treatment processes, heating and cooling rates are also interesting 
# for practical reasons
def heat_cool_rate(x_list, y_list, num_points):
    """
    Calculates heating/cooling rates from a given temperature curve (°C and seconds)
    and converts the heating/cooling rate to °C/minute.

    Parameters
    ----------
    x_list : TYPE list of floats
        DESCRIPTION. x-coordinate time min seconds
    y_list : TYPE list of floats
        DESCRIPTION.y-coordinate temperature
    num_points : TYPE int
        DESCRIPTION. number of datapoints

    Returns
    -------
    deriv_list : TYPE
        DESCRIPTION.

    """
    deriv_list = []
    x_list_len = len(x_list)
    y_list_len = len(y_list)
    if x_list_len != y_list_len:
        print("x_list and y_list are not equal in length!")
    
    for i, (x_coord, y_coord) in enumerate(zip(x_list, y_list)):
        if i < num_points-1:
            deriv = (y_list[i+1] - y_list[i]) / (x_list[i+1] - x_list[i]) * 60
            deriv_list.append(deriv)
    return deriv_list

heating_rate_surface = heat_cool_rate(t_list, To_slab_list_filter, num_steps_heating+num_steps_cooling)
heating_rate_center = heat_cool_rate(t_list, Tr_slab_list_filter, num_steps_heating+num_steps_cooling)
heating_rate_surface_filtered = savgol_filter(heating_rate_surface, window_size*2, polyorder+1)
heating_rate_center_filtered = savgol_filter(heating_rate_center, window_size*2, polyorder+1)

# we plot the results for the slab
plt.figure(figsize=(15,9))
plt.plot(t_list, To_slab_list_filter, '-g', lw=2, label='surface filtered')
plt.plot(t_list, Tr_slab_list_filter, '-r', lw=2, label='center filtered')
plt.plot(t_list, To_slab_list, '--g', lw=0.5, label='surface')
plt.plot(t_list, Tr_slab_list, '--r', lw=0.5, label='center')
plt.title(f"Slab with diameter {d} meters")
plt.ylabel('Temperature / °C')
plt.xlabel('Time / s')
plt.xlim([0, tmax])
plt.axhline(y=Th, color='k', linestyle='--', label=r'T$_\infty$')
plt.rcParams['xtick.major.pad'] = 6
plt.rcParams['ytick.major.pad'] = 6
plt.legend(loc='best', numpoints=1)
plt.grid()
plt.tight_layout()
plt.show()

# we plot the heating and cooling rates
plt.figure(figsize=(15,9))
plt.plot(t_list[:-1], heating_rate_surface_filtered, '-g', lw=2, label='heating/cooling rate filt. surface')
plt.plot(t_list[:-1], heating_rate_center_filtered, '-r', lw=2, label='heating/cooling rate filt. center')
plt.plot(t_list[:-1], heating_rate_surface, '--g', lw=0.5, label='heating/cooling rate surface')
plt.plot(t_list[:-1], heating_rate_center, '--r', lw=0.5, label='heating/cooling rate center')
plt.title(f"Slab with diameter {d} meters, heating/cooling rates")
plt.ylabel('Heating/Cooling rates / °C/minute')
plt.xlabel('Time /s')
# plt.yscale('log')
plt.legend(loc='best', numpoints=1)
plt.legend(loc='best', numpoints=1)
plt.tight_layout()
plt.show()

# for debugging 
if debugging == True:
    plt.figure(figsize=(15,9)) # looks ok
    plt.plot(h_list1, 'ro')
    plt.plot(h_list2, 'bo')
    plt.title("h_list")
    plt.show()
    plt.figure(figsize=(15,9)) # looks ok
    plt.plot(kw_list1, 'gx')
    plt.plot(kw_list2, 'yx')
    plt.title("kw_list")
    plt.show()
    plt.figure(figsize=(15,9)) # spiky!
    plt.plot(cpw_list1, 'b^')
    plt.plot(cpw_list2, 'r^')
    plt.title("cpw_list")
    plt.show()
    plt.figure(figsize=(15,9)) # looks ~ok
    plt.plot(alpha_list1, 'y.')
    plt.plot(alpha_list2, 'g.')
    plt.title("alpha_list")
    plt.show()
    plt.figure(figsize=(15,9)) # 
    plt.plot(T_cent_init_list, 'rv')
    # plt.plot(alpha_list2, 'g.')
    plt.title("T_cent_init_list")
    plt.show()

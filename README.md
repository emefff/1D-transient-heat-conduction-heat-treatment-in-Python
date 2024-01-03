# 1D-transient-heat-conduction-heat-treatment-in-Python

**This work is based on https://github.com/alexlib/trans_heat_cond/blob/master/LICENSE, thanks a lot. 
The attached code is still work in progress, it shows some weird behaviour under certain conditions (feedback on temperature by oscillating data that feeds back to data....).
To-Dos: functionalizing, cleaning up, find reason for weird oscillations etc.
Use at your own risk, any help is greatly appreciated! Thank you!**

The 1D transient heat equation is solved for conduction and heat transfer to the surface of a steel slab. This is useful for heat treatment, especially for estimating the soaking time. It is is mostly estimated if FEA is not available (rule of thumb for example for austenitizing: 1 hour per inch of steel thickness). 
While these rules of thumb are for sure enough for a first estimate, it is not economically viable and in certain steels it can even lead to detrimental results (for example: grain coarsening in tool steels that need very high Ta).

In the shared script, a slab with a thicknes of 250mm is heated in a surrounding medium of 870°C with a soaking time of 5 hours (half the time that above rule of thumb recommends). The needed 850°C are reached for about 5 minutes in the core of the slab, this can be enough for some steels. Compared to the rule of thumb, only half the time is needed. 
Subsequently, the slab is cooled with different heat transfer curve until is reaches room temperature. Some of the input data is still not physically correct, only assumed. 

As a result we get the following temperature curves for core and surface temperature (filtered results due to oscillations in the original data):

![temperature_curves_heating_and cooling](https://github.com/emefff/1D-transient-heat-conduction-heat-treatment-in-Python/assets/89903493/8a452010-0967-4baa-8180-616ff71776c7)

From that we can easily derive the heating and cooling rates (unit: °C / minute) which we can use for comparison with TTT-diagrams etc.:

![heating_cooling_rates](https://github.com/emefff/1D-transient-heat-conduction-heat-treatment-in-Python/assets/89903493/868ea857-99b0-4da4-8a55-46c10aa83724)

emefff@gmx.at



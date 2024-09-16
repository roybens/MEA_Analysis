TITLE calculates Na+/K+ charge overlap and excess Na+ influx COMMENT	Hallermann, de Kock, Stuart and Kole, Nature Neuroscience, 2012
	doi:10.1038/nn.3132ENDCOMMENT
NEURON {
        SUFFIX charge_     : changed "charge" to "charge_" because of conflicts with NEURON's "charge"
	USEION na READ ina
	USEION k READ ik
        RANGE vmax, vmin, tmax, tmin
        RANGE na_ch, na_ch_overl, overl
	RANGE na_ch_before_peak
	RANGE na_ch_after_peak
	RANGE na_ch_excess_ratio
        RANGE peak_reached
        RANGE peak_time
}

PARAMETER {	tStart (ms)
	tEnd (ms)
	peak_tolerance (mV)
	peak_lowest (mV)
}

ASSIGNED {
        v (millivolt)
        vmin (millivolt)
        tmin (ms)
        vmax (millivolt)
        tmax (ms)
        na_ch (milliamp/cm2)
        na_ch_overl (milliamp/cm2)
        na_ch_overl_tmp (milliamp/cm2)
	overl

	na_ch_excess_ratio
	na_ch_before_peak
	na_ch_after_peak
        peak_reached
        peak_time
	
	ina  (milliamp/cm2)
	ik  (milliamp/cm2)
}


INITIAL {
        vmin = 1e6
        tmin = 0
        vmax = -1e6
	tmax = 0
        peak_reached = 0
        peak_time = 0
	na_ch = 0
	na_ch_before_peak = 0
	na_ch_after_peak = 0
	na_ch_excess_ratio = 0
	na_ch_overl = 0
	overl = 0
:	tStart = 500
:	tEnd = 1000	
:	peak_tolerance = 0.1	(millivolt)
:	peak_lowest = -60	(millivolt)
}


BREAKPOINT {
VERBATIM
      if (t > tStart) {
		if (t < tEnd) {
			if (v < vmin) {
        		        vmin = v;
                		tmin = t;
		        }
		        if (v > vmax) {
        		        vmax = v;
	                	tmax = t;
		        }
 			na_ch = na_ch + ina;
			na_ch_overl_tmp = ina;				
			if (-ik > ina) {
                		na_ch_overl_tmp = -ik;
		        }
			na_ch_overl = na_ch_overl + na_ch_overl_tmp;
			if (na_ch !=  0) {	//na_ch is negative
                		overl = (na_ch - na_ch_overl) / na_ch;
			}
			if ( (v < vmax - peak_tolerance) && (v > peak_lowest) && (peak_reached == 0) ) {
				peak_reached = 1;
				peak_time = t;
			}
			if (peak_reached == 0) {
				na_ch_before_peak = na_ch_before_peak + ina;			
			} else {
				na_ch_after_peak = na_ch_after_peak + ina;
			}
			if (na_ch_before_peak != 0) {
				na_ch_excess_ratio = (na_ch_before_peak + na_ch_after_peak) / na_ch_before_peak;
			}
		}
	}
ENDVERBATIM
}


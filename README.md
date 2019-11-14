# regionalCO2

This code is based on modifications of CESM1.2.2/CAM4 and has been tested using the 2 degree finite volume grid.

No guarantees are given for its functionality. Please test thoroughly and double-check that your output looks as expected.

If you use this code, please cite the following two papers. If you only have space to cite one, please cite the second one (note that volume number, page numbers, and DOI needs to be added as the manuscript was under review at time of code publication):

1) Stuecker, M. F., C. M. Bitz, K. C. Armour, C. Proistosescu, S. M. Kang, S.-P. Xie, D. Kim, S. McGregor, W. Zhang, S. Zhao, W. Cai, Y. Dong, and F.-F. Jin (2018): Polar amplification dominated by local forcing and feedbacks, Nature Climate Change, 8, 1076-1081, doi:10.1038/s41558-018-0339-y

2) Stuecker, M. F., A. Timmermann, F.-F. Jin, C. Proistosescu, S. M. Kang, D. Kim, K.-S. Yun, E.-S. Chung, J.-E. Chu, C. M. Bitz, K. C. Armour, and M. Hayashi (2019): Strong remote control of future equatorial warming by off-equatorial forcing, Nature Climate Change, under review

The code reads in a dummy variable from a Netcdf file, which is used in the radiation code instead of the CO2 variable.

The netcdf forcing file should have the following format:

dimensions:
	time = 12 ;
	lat = 96 ;
	lon = 144 ;
variables:
	float lat(lat) ;
		lat:units = "degrees_north" ;
	float lon(lon) ;
		lon:units = "degrees_east" ;
	float DUMMY(time, lat, lon) ;
		DUMMY:units = "mol/mol" ;
		DUMMY:cell_method = "time: mean" ;
	float time(time) ;
		time:calendar = "noleap" ;
		time:units = "days since 1850-01-01 00:00" ;
	int date(time) ;
	int datesec(time) ;
		datesec:units = "seconds" ;
}

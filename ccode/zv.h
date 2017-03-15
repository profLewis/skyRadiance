#define FLOAT	double
#define LAM_TOL 0.00001
#define Sigma1MinCosDiff 0.00001

typedef struct{
	FLOAT	lambda_start;
	FLOAT	lambda_end;
	FLOAT	lambda_gap;
	FLOAT	angres[2];
	int	no_of_wavebands;
}Lambda;

/*
**	TTHG
*/

typedef struct{
	FLOAT	g1,g2;		/* asymmetry parameters */
	FLOAT	alpha;
}TTHG;

typedef struct{
	char	*filename;
	FLOAT	*albedo;
}Albedo;

/*
**	Aerosol Inputs:
*/

typedef struct{
	FLOAT	omega_A;	/* Aerosol single-scattering albedo: 	**
				** 0.6 -> urban / industrial		**
				** 0.9 -> continental			**
				** 1.0 -> maritime (pure scattering)	*/

	FLOAT	alpha;		/* Angstrom wavelength exponent:	**
				** generally in the range 0.5-1.6	**
				** though can be less for polluted 	**
				** atmosphere. 				**
				** av. conditions: alpha = 1.3 +/- 0.2	*/

	FLOAT	beta; 		/* Angstrom turbidity coefficient:	**
				** [Leckner, 1978]:			**
				** lat:	| 60N	| 45N	| 30N	| 0	**
				** -------------------------------------**
				** low 	| 0.01	| 0.047	| 0.047	| 0.047	**
				** med	| 0.057	| 0.077	| 0.097	| 0.117	**
				** high	| 0.093	| 0.187	| 0.375	| 0.375	**
				**					**
				** generally peaks in summer.		**
				** average mean values given by 	**
				** [Angstrom, 1961]:			**
				**					**
				** beta=[0.025+0.1cos^2(lat)]exp(-0.7h)	**
				**					**
				** h 	- height in km			**
				** lat 	- latitude			*/
	FLOAT	H_Penndorf;	/* Penndorf scale height parameter */
}Aerosol_Inputs;

/*
**	Ozone Inputs
*/

typedef struct{
	char	*filename;	/* filename of Ozone data */
	FLOAT	Co3;		/* ozone concentration (cm NTP)		**
				** typically, 0.4cm			**
				** usually in the range [0.2,0.6]	*/

	FLOAT	*Ko3;		/* Spectral Absorption coefficient	**
				** of Ozone (cm^-1)			*/
	FLOAT	b,c;		/* Ozone scale height parameters 	    **
				** b: altitude at which ozone conc. is max. **
				** c: ratio of maximum ozone conc. to	    **
				**    total ozone amount		    */
}Ozone_Inputs;

/*
**	Water Vapour Inputs
*/

typedef struct{
	char	*filename;	/* filename of Kw data */
	FLOAT	Cw;		/* Precipitable water content of atm.	**
				** (cm)					**
				** [Leckner, 1978]:			**
				**					**
				**	lat.		Cw		**
				** -------------------------------------**
				**	15N		4.1		**
				**	45N (summer)	2.9		**
				**	45N (winter)	0.9		**
				**	60N (summer)	2.1		**
				**	60N (winter)	0.4		**
				**	US62		1.4		**
				**					**
				** or [Leckner, 1978]:			**
				**					**
				** Cw = 181 density_W(0) / 0.795	**
				**					**
				** for standard atmosphere		**
				**					**
				** density_w = phi * Ps / (R * T(0))	**
				**					**
				** phi - relative humidity (fractional) **
				** R - Gas constant (461.51)		**
				** T(0) - absolute temperature (K)	**
				** Ps - water vapour saturation pressure**
				** Over water,				**
				** Ps = exp(-5416/T(0) + 26.23)	(Pa)	**
				**					**
				** other sources: e.g.			**
				** [Science data book, 1979]		**
				** or other tables of thermal data	*/

	FLOAT	*Kw;		/* spectral absorption coefficient of	**
				** water (cm-1)				**
				** tabular data [Leckner,1978]		*/

	FLOAT	beta_0;		/* water vapour scale height parameter	**
				**					**
				** beta_0 = exp(-(ln(Cw)-0.061Td)+0.715)**
				**	          _____________________ **
				**		 (       0.983	       )**
				**					*/

	FLOAT	Td;		/* dew point temperature 		*/
 
}WaterVapour_Inputs;

/*
**	Mixed Gasses
*/

typedef struct{
	char	*filename;	/* Kg filemane */
	FLOAT	*Kg;		/* effective absorption coefficient of	**
				** uniformly-mixed gasses		**
				** tabular data [Leckner,1978]		*/
}MixedGas_Inputs;

typedef struct{
	char	*filename;	/* EOD filemane */
	FLOAT	*EOD;
	
}EOD_Inputs;

typedef struct{

	FLOAT	E;		/* Extraterrestrial Irradiance */
	FLOAT	tau_R;		/* Rayleigh optical depth */
	FLOAT 	tau_A;		/* Aerosol optical depth */
	FLOAT	tau_As;		/* Aerosol scattering optical depth */
	FLOAT	tau_Aa;		/* Aerosol absorption optical depth */
	FLOAT	tau_O3;		/* Ozone optical depth */
	FLOAT	tau_W;		/* Water Vapour optical depth */
	FLOAT	tau_G;		/* Mixed Gasses optical depth */

	FLOAT	tau_a;		/* absorption optical depth */
	FLOAT	tau_s;		/* scattering optical depth */
}OpticalDepth;

typedef struct{
	FLOAT	A,B,C,D;	/* see [Frohlich & Shaw, 1980, p. 1775] */
}Rayleigh;

typedef struct{
	Rayleigh		rayleigh;
	EOD_Inputs		Eod;
	WaterVapour_Inputs	water_vapour;
	Ozone_Inputs		ozone;
	MixedGas_Inputs		mixed_gas;
	Aerosol_Inputs		aerosol;
	Albedo			albedo;
	TTHG			tthg;
	FLOAT			Legendre_coefficients[3];
	FLOAT			altitude[2];	/* altitude[0]: ground **
						** altitude[1]: sensor */
	int			altitude_flag;
}Atmospheric_Inputs;


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <hipl_format.h>
#include <time.h>
#include "useful3.h"
#include "wavefront.h"
#include "allocate.h"
#include "zv.h"
#include "error.h"

FLOAT	attenuation(tau)
FLOAT	tau;
{
	return(exp(-tau));
}

FLOAT	E_Julian(E0,D)
FLOAT	E0;
int	D;
{
	FLOAT	d,t;

	d=2.0*M_PI*(D-3)/365.0;
	t=0.01672*cos(d) + 1.0;
	t *= t;
	return(E0*t);
}

FLOAT	linear_interpolation(current_wavelength,wavelength,value)
FLOAT	current_wavelength,*wavelength,*value;
{
	FLOAT	gap,out;
	if(wavelength[1]==wavelength[0])return(value[1]);

	gap=1.0 - (wavelength[1]-current_wavelength)/(wavelength[1]-wavelength[0]);
	out= gap*(value[1]-value[0]) + value[0];
	return(out);
}

void	read_spectral_data_file(fp,data,lambda)
FLOAT	*data;
Lambda	*lambda;
FILE 	*fp;
{
	FLOAT	wavelength[2],value[2],current_wavelength,current_value;
	int	count=0;
	char	buffer[500];

	wavelength[0]=0.0;
	value[0]=0.0;
	
/*
**	lambda data expected in nm
*/
	current_wavelength=lambda->lambda_start;

	while(fgets(buffer,500,fp)){
		if(buffer[0]!='#'){
			if(sscanf(buffer,"%lf %lf",&wavelength[1],&value[1])!=2)error1("read_spectral_data_inputs:\terror reading data");
			while(wavelength[1]>=current_wavelength){
				current_value=linear_interpolation(current_wavelength,wavelength,value);
				data[count]=current_value;
				current_wavelength += lambda->lambda_gap;
				if(current_wavelength - LAM_TOL>lambda->lambda_end){
					fclose(fp);
					return;
				}
				count++;
			}
			wavelength[0]=wavelength[1];
			value[0]=value[1];
		}
				
	}
	fclose(fp);
	return;
}

void	read_height_file(N_h,fp,data,height_lut)
FLOAT	*data;
FLOAT	*height_lut;
FILE 	*fp;
int	N_h;
{
	FLOAT	height[2],value[2],current_height,current_value;
	int	count=0;
	char	buffer[500];

	height[0]=0.0;
	value[0]=0.0;
	current_height=height_lut[0];
	
/*
**	height data expected in km
*/

	while(fgets(buffer,500,fp)){
		if(buffer[0]!='#'){
			if(sscanf(buffer,"%lf %lf",&height[1],&value[1])!=2)error1("read_spectral_data_inputs:\terror reading data");
			while(height[1]>=current_height){
				current_value=linear_interpolation(current_height,height,value);
				data[count]=current_value;
				current_height = height_lut[++count];
				if(count>=N_h){
					fclose(fp);
					return;
				}
			}
			height[0]=height[1];
			value[0]=value[1];
		}
				
	}
	fclose(fp);
	return;
}

FILE	*open_spectral_data_file(ipfilename,str)
char	*ipfilename, *str;
{
	char	filename[1000];
	FILE	*fp;

	if(ipfilename==NULL){

		/* use default	*/
		strcpy(filename,getenv("ATMOS"));
		if(strlen(filename)==0)error1("open_spectral_data_file:\terror reading environmental variable ATMOS - has it been set?");
		strcat(filename,str);
	}else
		strcpy(filename,ipfilename);

	if(strcmp(filename,"-")==0)
		fp=stdin;
	else
		if(!(fp=fopen(filename,"r")))error2("open_spectral_data_file:\terror opening file",filename);

	return(fp);
}

FLOAT	*read_atmospheric_inputs(lambda,filename,str)
char	*filename,*str;
Lambda	*lambda;
{
	FILE	*fp;
	FLOAT	*out;
	
	fp=open_spectral_data_file(filename,str);

	out=d_allocate(lambda->no_of_wavebands);

	read_spectral_data_file(fp,out,lambda);
	return(out);
}

void	read_atmospheric_parameters_file(atmos_file,atmospheric_inputs)
Atmospheric_Inputs	*atmospheric_inputs;
char			*atmos_file;
{
	FILE	*fp;
	char	buffer[1000],str[100],filename[1000];
	int	start=0;
/*
**	load default values
*/
	atmospheric_inputs->ozone.filename=NULL;
	atmospheric_inputs->water_vapour.filename=NULL;
	atmospheric_inputs->mixed_gas.filename=NULL;
	atmospheric_inputs->albedo.filename=NULL;
	
	atmospheric_inputs->water_vapour.Cw=1.4;	/* [Leckner, 1978] */
	atmospheric_inputs->water_vapour.Td=10.02;	/* degrees C */

	atmospheric_inputs->ozone.Co3=0.3;		/* [Leckner, 1978] */
	atmospheric_inputs->ozone.b=20;			/* [Guzzi et al., 1987] */
	atmospheric_inputs->ozone.c=5;			/* [Guzzi et al., 1987] */
	atmospheric_inputs->aerosol.H_Penndorf=0.97;	/* [Guzzi et al., 1987] */
	atmospheric_inputs->aerosol.omega_A=0.9;	/* [Deepak & Gerber, 1983] */
	atmospheric_inputs->aerosol.alpha=1.3;		/* [Angstrom, 1961] */
	atmospheric_inputs->aerosol.beta=0.1;		/* [Angstrom, 1961] */
	atmospheric_inputs->rayleigh.A=0.00838;		/* [Frolich & Shaw, 1980] */
	atmospheric_inputs->rayleigh.B=3.916;		/* [Frolich & Shaw, 1980] */
	atmospheric_inputs->rayleigh.C=0.074;		/* [Frolich & Shaw, 1980] */
	atmospheric_inputs->rayleigh.D=0.050;		/* [Frolich & Shaw, 1980] */

	atmospheric_inputs->tthg.alpha=0.9618;		/* [Kattawar,] */
	atmospheric_inputs->tthg.g1=0.713;
	atmospheric_inputs->tthg.g2= -0.7598;
/*
**	read atmospheric parameter file
**	either from the file 'atmos_file'
**	or (if atmos_file == NULL) from '$ATMOS/default_irradiance_parameters'
*/
	fp=open_spectral_data_file(atmos_file,"/default_irradiance_parameters");
	while(fgets(buffer,1000,fp)){
		buffer[strlen(buffer)-1]=0;
		if(!start && strlen(buffer)!=0 && buffer[0] != '#' && (strcmp(buffer,"Atmosphere {")==0 || strcmp(buffer,"Atmosphere{")==0))start=1;
		if(start && strlen(buffer)!=0 && buffer[0] != '#'){
			if(sscanf(buffer,"%s",str)==1){
				if(strcmp(buffer,"Atmosphere {")==0 || strcmp(buffer,"Atmosphere{")==0){
					;
				}else if(strcmp(str,"rayleigh")==0 ||strcmp(str,"Rayleigh")==0){
/*
**	expect:
**		rayleigh A B C D
*/
					if(sscanf(buffer,"%s %lf %lf %lf %lf",str,&(atmospheric_inputs->rayleigh.A),&(atmospheric_inputs->rayleigh.B),&(atmospheric_inputs->rayleigh.C),&(atmospheric_inputs->rayleigh.D))!=5)error2("read_atmospheric_parameters_file:\terror reading rayleigh parameters in file",atmos_file);
				}else if(strcmp(str,"TTHG")==0 ||strcmp(str,"tthg")==0 ){
/*
**	expect:
**		TTHG alpha g1 g2
*/
					if(sscanf(buffer,"%s %lf %lf %lf",str,&(atmospheric_inputs->tthg.alpha),&(atmospheric_inputs->tthg.g1),&(atmospheric_inputs->tthg.g2))!=4)error2("read_atmospheric_parameters_file:\terror reading TTHG parameters in file",atmos_file);
				}else if(strcmp(str,"E0D")==0 ||strcmp(str,"Solar_Irradiance")==0 ||strcmp(str,"Solar_irradiance")==0 ||strcmp(str,"solar_irradiance")==0){
					if(sscanf(buffer,"%s %s",str,filename)==2){
						atmospheric_inputs->Eod.filename=c_allocate(strlen(filename)+1);
						strcpy(atmospheric_inputs->Eod.filename,filename);
					}else error2("read_atmospheric_parameters_file:\terror reading solar irradiance parameters in file",atmos_file);

				}else if(strcmp(str,"aerosol")==0 ||strcmp(str,"aerosols")==0 ||strcmp(str,"Aerosol")==0 ||strcmp(str,"Aerosols")==0){
/*
**	expect:
**		aerosol alpha beta single_scattering_albedo
*/
					if(sscanf(buffer,"%s %lf %lf %lf %lf",str,&(atmospheric_inputs->aerosol.alpha),&(atmospheric_inputs->aerosol.beta),&(atmospheric_inputs->aerosol.omega_A),&(atmospheric_inputs->aerosol.H_Penndorf))!=5){
						;
					}else if(sscanf(buffer,"%s %lf %lf %lf",str,&(atmospheric_inputs->aerosol.alpha),&(atmospheric_inputs->aerosol.beta),&(atmospheric_inputs->aerosol.omega_A))!=4)
						error2("read_atmospheric_parameters_file:\terror reading aerosol parameters in file",atmos_file);
				}else if(strcmp(str,"ozone")==0||strcmp(str,"Ozone")==0){
					if(sscanf(buffer,"%s %lf %lf %lf %s",str,&(atmospheric_inputs->ozone.Co3),&(atmospheric_inputs->ozone.b),&(atmospheric_inputs->ozone.c),filename)==5){
						atmospheric_inputs->ozone.filename=c_allocate(strlen(filename)+1);
						strcpy(atmospheric_inputs->ozone.filename,filename);

					}else if(sscanf(buffer,"%s %lf %lf %lf",str,&(atmospheric_inputs->ozone.Co3),&(atmospheric_inputs->ozone.b),&(atmospheric_inputs->ozone.c))!=4){
						;
					}else if(sscanf(buffer,"%s %lf %s",str,&(atmospheric_inputs->ozone.Co3),filename)==3){
						atmospheric_inputs->ozone.filename=c_allocate(strlen(filename)+1);
						strcpy(atmospheric_inputs->ozone.filename,filename);
					}else 	if(sscanf(buffer,"%s %lf",str,&(atmospheric_inputs->ozone.Co3))!=2)
						error2("read_atmospheric_parameters_file:\terror reading ozone parameters in file",atmos_file);
				}else if(strcmp(str,"albedo")==0||strcmp(str,"Albedo")==0){	/* surface albedo */
					if(sscanf(buffer,"%s %s",str,filename)==2){
						atmospheric_inputs->albedo.filename=c_allocate(strlen(filename)+1);
						strcpy(atmospheric_inputs->albedo.filename,filename);
					}else error2("read_atmospheric_parameters_file:\terror reading albedo parameters in file",atmos_file);
				}else if(strcmp(str,"water_vapour")==0||strcmp(str,"Water_vapour")==0||strcmp(str,"water_vapor")==0||strcmp(str,"Water_vapor")==0){
					if(sscanf(buffer,"%s %lf %lf %s",str,&(atmospheric_inputs->water_vapour.Cw),&(atmospheric_inputs->water_vapour.Td),filename)==4){
						atmospheric_inputs->water_vapour.filename=c_allocate(strlen(filename)+1);
						strcpy(atmospheric_inputs->water_vapour.filename,filename);
					} else if(sscanf(buffer,"%s %lf %lf",str,&(atmospheric_inputs->water_vapour.Cw),&(atmospheric_inputs->water_vapour.Td))==3){
						;
					}else 	if(sscanf(buffer,"%s %lf",str,&(atmospheric_inputs->water_vapour.Cw))!=2)error2("read_atmospheric_parameters_file:\terror reading water_vapour parameters in file",atmos_file);
				}else if(strcmp(str,"mixed_gas")==0||strcmp(str,"Mixed_gas")==0){
					if(sscanf(buffer,"%s %s",str,filename)!=2)error2("read_atmospheric_parameters_file:\terror reading mixed_gas parameters in file",atmos_file);
					atmospheric_inputs->mixed_gas.filename=c_allocate(strlen(filename)+1);
					strcpy(atmospheric_inputs->mixed_gas.filename,filename);
				}else if(strcmp(str,"}")==0)start=0;
				else error3("read_atmospheric_parameters_file:\tunrecognised parameter in file",atmos_file,str);
			}
		}
	}
	if(atmospheric_inputs->water_vapour.Cw!=0.0)atmospheric_inputs->water_vapour.beta_0 = exp(-(log(atmospheric_inputs->water_vapour.Cw) - 0.061*atmospheric_inputs->water_vapour.Td + 0.715)/0.983);else atmospheric_inputs->water_vapour.beta_0=0.0;
	return;
}

void	read_tau_file(JulianDay,tau_file,lambda,atmospheric_inputs,tau)
char			*tau_file;
Atmospheric_Inputs	*atmospheric_inputs;
Lambda			*lambda;
OpticalDepth		*tau;
int			JulianDay;
{
	FILE	*fp;
	char	buffer[1000];
	FLOAT	wavelength_nm,wavelength_um,H_R=1.0,H_A=1.0,H_O3=1.0,H_W=1.0,Aerosol_scale_height(),Rayleigh_scale_height(),Ozone_scale_height(),Water_vapour_scale_height();
	int	count=0;

	if(atmospheric_inputs->altitude_flag==1){
		H_R=Rayleigh_scale_height(atmospheric_inputs->altitude);
		H_A=Aerosol_scale_height(atmospheric_inputs->altitude,atmospheric_inputs->aerosol.H_Penndorf);
		H_O3=Ozone_scale_height(atmospheric_inputs->altitude,atmospheric_inputs->ozone.b,atmospheric_inputs->ozone.c);
		H_W=Water_vapour_scale_height(atmospheric_inputs->altitude,atmospheric_inputs->water_vapour.beta_0);
	}

	fp=open_spectral_data_file(tau_file,"");
	atmospheric_inputs->ozone.Ko3=d_allocate(lambda->no_of_wavebands);
	atmospheric_inputs->water_vapour.Kw=d_allocate(lambda->no_of_wavebands);
	atmospheric_inputs->mixed_gas.Kg=d_allocate(lambda->no_of_wavebands);
/*
**	albedo is not part of the default tau ip data -> read it from file
*/
	atmospheric_inputs->albedo.albedo=read_atmospheric_inputs(lambda,atmospheric_inputs->albedo.filename,"/ALBEDO");

	while(fgets(buffer,1000,fp)){
		buffer[strlen(buffer)-1]=0;
		if(strlen(buffer)!=0 && buffer[0] != '#'){
			if(sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&wavelength_nm,&wavelength_um,&(tau[count].tau_R),&(tau[count].tau_A),&(atmospheric_inputs->ozone.Ko3[count]),&(tau[count].tau_O3),&(atmospheric_inputs->water_vapour.Kw[count]),&(tau[count].tau_W),&(atmospheric_inputs->mixed_gas.Kg[count]),&(tau[count].tau_G),&(tau[count].E))==11){
				tau[count].E=E_Julian(tau[count].E,JulianDay);
				tau[count].tau_As=atmospheric_inputs->aerosol.omega_A * tau[count].tau_A;
				tau[count].tau_Aa = (1.0 - atmospheric_inputs->aerosol.omega_A) * tau[count].tau_A;
				tau[count].tau_a = tau[count].tau_O3 + tau[count].tau_W + tau[count].tau_G + tau[count].tau_Aa;
				tau[count].tau_s = tau[count].tau_R + tau[count].tau_As;

				count++;
				if(count>lambda->no_of_wavebands)error2("read_tau_file:\ttoo many lines in tau file",tau_file);
			}
		}
	}	
	
	return;
}

void	read_tau_data(JulianDay,tau_file,atmos_file,lambda,atmospheric_inputs,tau)
char			*tau_file;
Atmospheric_Inputs	*atmospheric_inputs;
Lambda			*lambda;
char			*atmos_file;
OpticalDepth		*tau;
int			JulianDay;
{

	read_atmospheric_parameters_file(atmos_file,atmospheric_inputs);
/*
**	read optical depth data
*/
	read_tau_file(JulianDay,tau_file,lambda,atmospheric_inputs,tau);
	
	return;
}

void	read_atmospheric_input_data(atmos_file,lambda,atmospheric_inputs)
Atmospheric_Inputs	*atmospheric_inputs;
Lambda			*lambda;
char			*atmos_file;
{

	read_atmospheric_parameters_file(atmos_file,atmospheric_inputs);
	
/*
**	read spectral absorption coefficient data
*/
	atmospheric_inputs->water_vapour.Kw=read_atmospheric_inputs(lambda,atmospheric_inputs->water_vapour.filename,"/KW");
	atmospheric_inputs->ozone.Ko3=read_atmospheric_inputs(lambda,atmospheric_inputs->ozone.filename,"/KO3");
	atmospheric_inputs->mixed_gas.Kg=read_atmospheric_inputs(lambda,atmospheric_inputs->mixed_gas.filename,"/KG");
	atmospheric_inputs->Eod.EOD=read_atmospheric_inputs(lambda,atmospheric_inputs->Eod.filename,"/E0D");
	
	atmospheric_inputs->albedo.albedo=read_atmospheric_inputs(lambda,atmospheric_inputs->albedo.filename,"/ALBEDO");
	return;
}

/*
**	Aerosol scale height
**	[Guzzi et al., 1987]
*/
FLOAT	Aerosol_scale_height1(height,H_Penndorf)
FLOAT	height,H_Penndorf;
{
	return(attenuation(height/H_Penndorf));
}

FLOAT	Aerosol_scale_height(height,H_Penndorf)
FLOAT	*height,H_Penndorf;
{
/*
**	height[0]:	lower height
**	height[1]:	upper height
*/
	return(fabs(Aerosol_scale_height1(height[0],H_Penndorf)-Aerosol_scale_height1(height[1],H_Penndorf)));
}

/*
**	Rayleigh scale height
**	[Guzzi et al., 1987]
*/

FLOAT	Rayleigh_scale_height1(height)
FLOAT	height;
{
	return(attenuation(0.1188*height + 0.00116*height*height));
}

FLOAT	Rayleigh_scale_height(height)
FLOAT	*height;
{
/*
**	height[0]:	lower height
**	height[1]:	upper height
*/
	return(fabs(Rayleigh_scale_height1(height[0])-Rayleigh_scale_height1(height[1])));
}

/*
**	Ozone scale height
**	[Guzzi et al., 1987]
**	[Lacis and Hansen, 1974]
*/

FLOAT	Ozone_scale_height1(height,b,c)
FLOAT	height,b,c;
{
	return((1.0 + exp(-b/c))/(1.0+exp((height-b)/c)));
}

FLOAT	Ozone_scale_height(height,b,c)
FLOAT	*height,b,c;
{
/*
**	height[0]:	lower height
**	height[1]:	upper height
*/
	return(fabs(Ozone_scale_height1(height[0],b,c)-Ozone_scale_height1(height[1],b,c)));
}

/*
**	water vapour scale height
**	[Guzzi et al., 1987]
*/

FLOAT	Water_vapour_scale_height1(height,beta_0)
FLOAT	height,beta_0;
{
	return(attenuation(-beta_0*height));
}

FLOAT	Water_vapour_scale_height(height,beta_0)
FLOAT	*height,beta_0;
{
/*
**	height[0]:	lower height
**	height[1]:	upper height
*/
	return(fabs(Water_vapour_scale_height1(height[0],beta_0)-Water_vapour_scale_height1(height[1],beta_0)));
}

void	calculate_tau(JulianDay,i,tau,lambda,atmospheric_inputs,scale_height_lut)
OpticalDepth		*tau;
Atmospheric_Inputs	*atmospheric_inputs;
FLOAT			lambda;
int			i,JulianDay;
Image_characteristics	*scale_height_lut;
{

	FLOAT	conc,H_R=1.0,H_A=1.0,H_O3=1.0,H_W=1.0;
/*
**	lambda in um
*/
	if(atmospheric_inputs->altitude_flag==1){
		if(scale_height_lut->data.fdata==NULL){
/*
**	model scale heights (normal incidence)
*/
			H_R=Rayleigh_scale_height(atmospheric_inputs->altitude);
			H_A=Aerosol_scale_height(atmospheric_inputs->altitude,atmospheric_inputs->aerosol.H_Penndorf);
			H_O3=Ozone_scale_height(atmospheric_inputs->altitude,atmospheric_inputs->ozone.b,atmospheric_inputs->ozone.c);
			H_W=Water_vapour_scale_height(atmospheric_inputs->altitude,atmospheric_inputs->water_vapour.beta_0);
		}else{
/*
**	LUT-generated scale heights accounted for later -> set H* = 1.0
*/
			H_R=H_A=H_O3=H_W=1.0;
		}
	}
		
	tau->tau_R = H_R*atmospheric_inputs->rayleigh.A*pow(lambda,(-(atmospheric_inputs->rayleigh.B + atmospheric_inputs->rayleigh.C*lambda + atmospheric_inputs->rayleigh.D/lambda)));
	tau->tau_A = H_A*atmospheric_inputs->aerosol.beta * pow(lambda,-atmospheric_inputs->aerosol.alpha);
	tau->tau_As = atmospheric_inputs->aerosol.omega_A * tau->tau_A;
	tau->tau_Aa = (1.0 - atmospheric_inputs->aerosol.omega_A) * tau->tau_A;
	
	tau->tau_O3 = H_O3*atmospheric_inputs->ozone.Co3 * atmospheric_inputs->ozone.Ko3[i];
	
	conc = atmospheric_inputs->water_vapour.Cw * atmospheric_inputs->water_vapour.Kw[i];
	tau->tau_W = 0.2385 * H_W* conc / pow(1.0 + 20.07*H_W*conc,0.45);
	/*tau->tau_W = 0.2385 * H_W* conc / pow(1.0 + 20.07*conc,0.45);*/
	
	conc = atmospheric_inputs->mixed_gas.Kg[i];
	tau->tau_G = 1.41 * H_R*conc / pow(1.0 + 118.93*H_R*conc,0.45);
	/*tau->tau_G = 1.41 * H_R*conc / pow(1.0 + 118.93*conc,0.45);*/

	tau->tau_a = tau->tau_O3 + tau->tau_W + tau->tau_G + tau->tau_Aa;
	tau->tau_s = tau->tau_R + tau->tau_As;

	tau->E = E_Julian(atmospheric_inputs->Eod.EOD[i],JulianDay);

	return;

}

void	calculate_optical_depth(JulianDay,tau,lambda,atmospheric_inputs,scale_height_lut)
OpticalDepth		*tau;
Atmospheric_Inputs	*atmospheric_inputs;
Lambda			*lambda;
int			JulianDay;
Image_characteristics	*scale_height_lut;
{
	int	i;
	FLOAT	wavelength;


	for(i=0;i<lambda->no_of_wavebands;i++){
		/* utilise wavelength in um */
		wavelength=(lambda->lambda_start + i*lambda->lambda_gap)/1000.0;
		calculate_tau(JulianDay,i,&(tau[i]),wavelength,atmospheric_inputs,scale_height_lut);
	}

	return;
		
}

void	output_tau_data(tau,lambda,atmospheric_inputs)
OpticalDepth		*tau;
Lambda			*lambda;
Atmospheric_Inputs	*atmospheric_inputs;
{
	int	i;
	FLOAT	wavelength;

	fprintf(stdout,"# Lewis version of tau: Version 1.0 (Mar 15 1993)\n");
	fprintf(stdout,"# \n");
	fprintf(stdout,"# LambdaSList: Tabulation of atmosphere parameter only dependent components of Z & V sky model.\n");
	fprintf(stdout,"# Ozone concentration coefficient : C_O3 = %.4lf\n",atmospheric_inputs->ozone.Co3);
	fprintf(stdout,"# Precipital Water vapour concentration coefficient : C_W =  %.4lf\n",atmospheric_inputs->water_vapour.Cw);
	fprintf(stdout,"# Angstrom formula Total Aerosol optical thickness coefficients : alpha = %.4lf, beta = %l.4f, omega_A = %.4lf\n",atmospheric_inputs->aerosol.alpha,atmospheric_inputs->aerosol.beta,atmospheric_inputs->aerosol.omega_A);
	fprintf(stdout,"#\tnm\tum\ttau_R\t\ttau_A\t\tK_O3\t\ttau_O3\t\tK_W\t\ttau_W\t\tK_G\t\ttau_G\t\tE\n");
	fprintf(stdout,"# \n");

	for(i=0;i<lambda->no_of_wavebands;i++){
		wavelength=(lambda->lambda_start + i*lambda->lambda_gap);
		fprintf(stdout,"\t%.1lf\t%.4lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.2lf\n",wavelength,wavelength/1000.0,tau[i].tau_R,tau[i].tau_A,atmospheric_inputs->ozone.Ko3[i],tau[i].tau_O3,atmospheric_inputs->water_vapour.Kw[i],tau[i].tau_W,atmospheric_inputs->mixed_gas.Kg[i],tau[i].tau_G,tau[i].E);
		fflush(stdout);
	}

}

FLOAT	Henyey_Greenstein(g,cost,sign)
FLOAT	g,cost,sign;
{
	FLOAT	out,g2;

	g2=g*g;

	out=(1.0-g2)/pow(1.0 + g2 + sign*(2.0*g*cost),1.5);
	
	return(out);
}



/*
**	two-term Henyey-Greenstein phase function
**	for aerosols
*/

FLOAT	TwoTermHenyeyGreenstein(tthg,cost)		/* cost -> cos(phase) */
FLOAT	cost;
TTHG	*tthg;
{
/*
**	note that we expect g1 to be +ve (forward scattering)
**	and g2 -ve (back scattering)
*/
	FLOAT	out,g12,g22,numer1,numer2,denom1,denom2;
	 
	g12=tthg->g1*tthg->g1;
	g22=tthg->g2*tthg->g2;

	numer1=tthg->alpha*(1.0-g12);
	numer2=(1.0-tthg->alpha)*(1.0-g22);
	denom1=pow((1.0+g12-(2.0*tthg->g1*cost)),1.5);
	denom2=pow((1.0+g22-(2.0*tthg->g2*cost)),1.5);
	out=(numer1/denom1)+(numer2/denom2);	
	return(out);
}

/*
**	calculate Rayleigh (molecular) phase function
*/

FLOAT	Phase_Rayleigh(cos2t)	/* cos2t -> cos^2(phase) */
FLOAT	cos2t;
{
	FLOAT	out;
	out = 0.75*(1.0+cos2t);
	return(out);
}

/*
**	calculate total phase function
*/

FLOAT	PF_Total(tau_R, p_R, tau_As, p_A)
FLOAT	tau_R, p_R, tau_As, p_A;
{
	FLOAT	out;

	out = (tau_R*p_R + tau_As*p_A)/(tau_R+tau_As);

	return(out);
}

/*
** Calculate the Zeroth, First & Second Legendre moments of the total scattering function 
** As well as requiring the integrated value of p1 the components of p, 
** the total scattering phase function,
** are required. p0 is required to normalise the values of p.
**
** Integrate from 0 to Pi:
**              p0 = 0.5 * p(theta)*sin(theta)*dtheta
**              p1 = 1.5 * p(theta)*cos(theta)*sin(theta)*dtheta.
**              p2 = 0.5 * p(theta)*sin(theta)*(3*cos^2(theta)-1) *dtheta
**
** The normalisation condition for p is p0 = 1.
**
** As the tau's vary with wavelength this is done for each wavelength ..
**
**	**--> this function taken (almost verbatim) from A.Newton's sky (g++) code <--**
**
*/

void	Legendre(Legendre_coeffs,tthg,tau_R,tau_As,dn)
FLOAT	tau_R;
FLOAT	tau_As;
int	dn;
TTHG	*tthg;
FLOAT	*Legendre_coeffs;
{
	int	i=0;
	FLOAT	dtheta,theta=0.0,cost=0.0,sint=0.0,p_R=0.0,p_A=0.0,p_T=0.0; 
	FLOAT	pt=0.0;
	FLOAT	cos2t,sin2t;

	dtheta=M_PI/(FLOAT)dn;
	for(i=0;i<3;i++)Legendre_coeffs[i]=0.0;

	for( i=0, theta=dtheta/2.0; i<dn; i++, theta += dtheta ) {
		cost=cos(theta);
		sint=sin(theta);
		cos2t=cost*cost;
		sin2t=sint*sint;

/*
** Rayleigh & Aerosol phase function components (independent of lambda)
*/
		p_R=Phase_Rayleigh(cos2t);		/* Rayleigh scattering phase fn */
		p_A=TwoTermHenyeyGreenstein(tthg,cost);	/* Aerosol scattering phase function */
		p_T=PF_Total(tau_R, p_R, tau_As, p_A);
		pt=p_T*sint;
		Legendre_coeffs[0]+=pt; Legendre_coeffs[1]+=(pt*cost); Legendre_coeffs[2]+=(pt*0.5*((3*cos2t)-1.0));
	}
	Legendre_coeffs[0]*=0.5*dtheta; Legendre_coeffs[1]*=1.5*dtheta; Legendre_coeffs[2]*=dtheta;

	return;
}

/*
**	take account of Earth orbit about the Sun
**	to find E0 for Julian day D
*/
	

FLOAT	R_aux(edtheta,cos_theta)
FLOAT	edtheta,cos_theta;
{
	FLOAT	out;

	out = 1.0 + 1.5*cos_theta + (1.0 - 1.5*cos_theta)*edtheta;

	return(out);
}

FLOAT	m_Kasten(theta,a,b,c)	/* theta in radians */
FLOAT	theta,a,b,c;
{
/*
**	see [Kasten, 1965]:
**	Kasten, F., "A new table and approximation formula
**	for the relative optical air mass", 
*/
	FLOAT	gamma;	/* solar elevation */
	FLOAT	m;	/* relative optical air mass */

	gamma = (M_PI/2.0) - theta;

	m = 1.0/(sin(gamma) + a*pow((gamma+b),-c));

	return(m);

}

FLOAT	m_Iqbal(theta)	/* in degrees */
FLOAT	theta;
{
	FLOAT	out,a,b;

	a=93.885-theta;

	b=0.15*pow(a,-1.253);

	out=1.0/(cos(DTOR(theta))+b);

	return(out);
}
/* single scattering transmission */
FLOAT	sigma1_ZV(edtheta,costheta)
FLOAT	*edtheta,*costheta;
{
	FLOAT	numer,denom,sigma1;

	numer=edtheta[1]-edtheta[0];
	denom=costheta[1]-costheta[0];	/* should already be checked so as not to == 0 */

	sigma1=0.25*numer/denom;
	return(sigma1);
}

FLOAT	Sigma1(tau_s,theta,edtheta,costheta)
FLOAT	tau_s,theta,*edtheta,*costheta;
{
	FLOAT	sigma1,sp,sm,ctp[2],ctm[2],ep[2],em[2],cosdiff,dtheta=0.01;

	cosdiff=costheta[1]-costheta[0];
/*
**	lifted from A. Newton's code in
**	atmosphere.cc
*/
	if(fabs(cosdiff)<Sigma1MinCosDiff){
		/* linearly interpolate sigma1 from interval +/- dtheta on either side of theta0 */
		ctp[1]=cos(theta+dtheta);
		ctm[1]=cos(theta-dtheta);
		ep[1]=attenuation(tau_s/ctp[1]);
		em[1]=attenuation(tau_s/ctm[1]);
		ep[0]=em[0]=edtheta[0];
		ctp[0]=ctm[0]=costheta[0];
		sp=sigma1_ZV(ep,ctp);
		sm=sigma1_ZV(em,ctm);
		sigma1=(sm+sp)/2.0;
	}else{
		sigma1=sigma1_ZV(edtheta,costheta);
	}

	return(sigma1);
}

FLOAT	E_d(E_0D,tau,m)
FLOAT	E_0D,tau,m;
{
	return(E_0D*exp(-tau*m));
}

FLOAT	E0Dcostheta0byPi(E_0D,cos_theta0)
FLOAT	E_0D,cos_theta0;
{
	return(E_0D*cos_theta0/M_PI);
}

FLOAT	exptau_abycostheta0(tau_a,cos_theta0)
FLOAT	tau_a,cos_theta0;
{
	return(attenuation(tau_a/cos_theta0));
}

FLOAT	Sigma(i,tau,atmospheric_inputs,p_R,p_A,cos_theta,edtheta,Rtheta,sigma1)
int			i;
OpticalDepth		*tau;
Atmospheric_Inputs	*atmospheric_inputs;
FLOAT			*cos_theta,*edtheta,*Rtheta,p_R,p_A,sigma1;
{
	FLOAT	t1,t2,t3,omr;

	t2=0.5*(edtheta[1]+edtheta[0]);
	t3=(PF_Total(tau->tau_R,p_R,tau->tau_As,p_A)-((3.0+atmospheric_inputs->Legendre_coefficients[1])*cos_theta[1]*cos_theta[0]));
	omr=1.0-atmospheric_inputs->albedo.albedo[i];
	t1=((omr*Rtheta[1] + (2.0*atmospheric_inputs->albedo.albedo[i]))*Rtheta[0])/(4.0+((3.0-atmospheric_inputs->Legendre_coefficients[1])*omr*tau->tau_s));
	return(t1-t2+(t3*sigma1));
}

FLOAT	Lsky(energy,sigma_scat,absorption)
FLOAT	energy,sigma_scat,absorption;
{
	return(energy*sigma_scat*absorption);
}

FLOAT	CosPsi(theta,phi,theta0,phi0)
FLOAT	theta,phi,theta0,phi0;
{
	FLOAT	st,st0,ct,ct0,cos_of_angle;

	st=sin(theta);
	ct=cos(theta);
	st0=sin(theta0);
	ct0=cos(theta0);

	cos_of_angle=(ct*ct0)-(st*st0*cos(phi-phi0));

	cos_of_angle=MIN(MAX(-1.0,cos_of_angle),1.0);

	return(cos_of_angle);
}

void	write_direct_header()
{
	fprintf(stdout,"# Sky Radiance Data from Lewis' Zibbordi & Voss model :\n");
	fprintf(stdout,"# Columns: 1. Lambda (nm), 2. Direct, 3. Diffuse (total), 4. Diffuse (cosine), 5. Ratio Direct to Diffuse, 6. Percent diffuse(1), 7. Total Optical Depth, 8. downward Direct path Transmittance, 9. upward path Transmittance (vertical), 10. Percent diffuse(2) on a horizontal plane, 11. Hemispherical integral cf. 2 Pi.\n");
	return;
}

FLOAT	lut_scale_height(h,frame_no,scale_height_lut,theta)
int	frame_no;
FLOAT	theta;	/* radians */
FLOAT	h;
Image_characteristics	*scale_height_lut;
{
	FLOAT	out;
	int	itheta,iheight;

	itheta=(int)(theta/(0.5*M_PI)+0.5);

	if(h<=19.6)
		iheight=(int)((h/0.1)+0.5);
        else if(h<=50.0)
		iheight=(int)(((h-19.6)/0.2)+0.5);
        else iheight=(int)(((h-50.0)/0.5)+0.5);

	out= *(scale_height_lut->data.fdata + scale_height_lut->hd.rows*scale_height_lut->hd.cols*frame_no + scale_height_lut->hd.cols*iheight + itheta);

	return(out);
}

void	create_skymap(direct_op,sun_position,lambda,skydata,atmospheric_inputs,tau,JulianDay,scale_height_lut)
Atmospheric_Inputs	*atmospheric_inputs;
Image_characteristics	*skydata,*scale_height_lut;
OpticalDepth		*tau;
Lambda			*lambda;
FLOAT			*sun_position;		/* theta_0, phi_0 */
int			direct_op,JulianDay;
{
	OpticalDepth	Tau;
	int	i,j,k,dn=100;
	FLOAT	p_R,p_A,cos_theta[2],sin_theta[2],theta,phi,m,m0,cos_phase_angle,cos2_phase_angle;
	FLOAT	sigma1,out,phi0,theta0;
	FLOAT	trans_down,trans_up,tau_T,E_direct,E_direct_cosine;
	FLOAT	Rtheta[2],edtheta[2],absorption,energy,scattering,E_diffuse,E_diffuse_cosine,E_diffuse_total,halfHemiInt;
	FLOAT	dtheta,dphi,domega,domega1,domega2;
	FLOAT	wavelength,direct2diffuse,diffusepercent1,diffusepercent2;

	theta0=DTOR(sun_position[0]);
	cos_theta[0]=cos(theta0);
	sin_theta[0]=sin(theta0);
	phi0=DTOR(sun_position[1]);

	if(scale_height_lut->data.fdata==NULL){
		m0=1.0/cos_theta[0];		/* default plane parallel atmosphere */

		m=m_Iqbal(sun_position[0]);
	}

	dtheta=(1.0/(FLOAT)(skydata->hd.rows-1))*M_PI*0.5;
	dphi=(1.0/(FLOAT)(skydata->hd.cols-1))*2.0*M_PI;
	domega1=dtheta*dphi;
	domega2=4.0*sin(dtheta/2.0)*sin(dphi/2.0);

	if(direct_op){
		write_direct_header();
	}

	for(i=0;i<lambda->no_of_wavebands;i++){
		wavelength=(lambda->lambda_start + i*lambda->lambda_gap);

		E_diffuse_cosine=0.0;E_diffuse_total=0.0;halfHemiInt=0.0;

		if(scale_height_lut->data.fdata==NULL){
			tau_T = tau[i].tau_a + tau[i].tau_s;
			E_direct=E_d(tau[i].E,tau_T,m);
			E_direct_cosine=E_direct*cos_theta[0];

			trans_down=attenuation(tau_T*m);
			trans_up=attenuation(tau_T);

			energy=E0Dcostheta0byPi(tau[i].E,cos_theta[0]);

			absorption=exptau_abycostheta0(tau[i].tau_a,cos_theta[0]);

			edtheta[0]=attenuation(tau[i].tau_s/cos_theta[0]);

			Rtheta[0]=R_aux(edtheta[0],cos_theta[0]);
/*
**	numerically evaluate first 3 terms of Legendre polymonial expansion
*/
			Legendre(&(atmospheric_inputs->Legendre_coefficients[0]),&(atmospheric_inputs->tthg),tau[i].tau_R,tau[i].tau_As,dn);
		}

		for(j=0;j<skydata->hd.rows;j++){

			theta= j*dtheta;

			if(scale_height_lut->data.fdata!=NULL){

				Tau.tau_R=lut_scale_height(0,scale_height_lut,theta)*tau[i].tau_R;
				Tau.tau_G=lut_scale_height(0,scale_height_lut,theta)*tau[i].tau_G;
				Tau.tau_W=lut_scale_height(1,scale_height_lut,theta)*tau[i].tau_W;
				Tau.tau_A=lut_scale_height(2,scale_height_lut,theta)*tau[i].tau_A;
				Tau.tau_As=lut_scale_height(2,scale_height_lut,theta)*tau[i].tau_As;
				Tau.tau_Aa=lut_scale_height(2,scale_height_lut,theta)*tau[i].tau_Aa;
				Tau.tau_O3=lut_scale_height(3,scale_height_lut,theta)*tau[i].tau_O3;

				tau[i].tau_a = Tau.tau_O3 + Tau.tau_W + Tau.tau_G + Tau.tau_Aa;
				tau[i].tau_s = Tau.tau_R + Tau.tau_As;
				tau_T = tau[i].tau_a + tau[i].tau_s;
				E_direct=E_d(tau[i].E,tau_T,m);
				E_direct_cosine=E_direct*cos_theta[0];

				trans_down=attenuation(tau_T);
				trans_up=attenuation(tau_T/m);

				energy=E0Dcostheta0byPi(tau[i].E,cos_theta[0]);

				absorption=exptau_abycostheta0(tau[i].tau_a,cos_theta[0]);

				edtheta[0]=attenuation(tau[i].tau_s/cos_theta[0]);

				Rtheta[0]=R_aux(edtheta[0],cos_theta[0]);
/*
**	numerically evaluate first 3 terms of Legendre polymonial expansion
*/
				Legendre(&(atmospheric_inputs->Legendre_coefficients[0]),&(atmospheric_inputs->tthg),tau[i].tau_R,tau[i].tau_As,dn);
			}
			sin_theta[1]=sin(theta);
			cos_theta[1]=cos(theta);

			if(cos_theta[1]==0.0)cos_theta[1]=1e-10;

			domega=sin_theta[1]*domega2;

			edtheta[1]=attenuation(tau[i].tau_s/cos_theta[1]);

			Rtheta[1]=R_aux(edtheta[1],cos_theta[1]);

			sigma1=Sigma1(tau[i].tau_s,theta,edtheta,cos_theta);

			for(k=0;k<skydata->hd.cols;k++){
				phi=k*dphi + M_PI;
				cos_phase_angle=CosPsi(theta,phi,theta0,phi0);
				cos2_phase_angle = cos_phase_angle*cos_phase_angle;
				p_R=Phase_Rayleigh(cos2_phase_angle);	/* Rayleigh scattering phase fn */
				p_A=TwoTermHenyeyGreenstein(&(atmospheric_inputs->tthg),cos_phase_angle); /* Aerosol scattering phase function */

/* scattering transmission coefficient of the atmosphere */
				/*scattering=Sigma(tau[i].tau_s,atmospheric_inputs->albedo.albedo[i],atmospheric_inputs->Legendre_coefficients[1],tau[i].tau_R,tau[i].tau_As,p_R,p_A,cos_theta0,cos_theta,edtheta0,edtheta,Rtheta0,Rtheta,sigma1);*/
				scattering=Sigma(i,&(tau[i]),atmospheric_inputs,p_R,p_A,cos_theta,edtheta,Rtheta,sigma1);

				out=Lsky(energy,scattering,absorption);

				*(skydata->data.fdata + i*skydata->hd.cols*skydata->hd.rows + j*skydata->hd.cols + k) = out;

				if(k==0 || k==skydata->hd.cols -1){
					E_diffuse=out*domega/2.0;
					halfHemiInt+=sin_theta[1]*domega2/2.0;
				}else{
					E_diffuse=out*domega;
					halfHemiInt+=sin_theta[1]*domega2;
				}
				E_diffuse_total += E_diffuse;
				E_diffuse_cosine += E_diffuse*cos_theta[1];
				
			}
		}
		if(direct_op){
			direct2diffuse=E_direct/E_diffuse_total;
			diffusepercent1=(100.0*E_diffuse_total/(E_direct+E_diffuse_total));
			diffusepercent2=(100.0*E_diffuse_cosine/(E_direct_cosine+E_diffuse_cosine));
			fprintf(stdout,"%.2lf\t%.2lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n",wavelength,E_direct,E_diffuse_total,E_diffuse_cosine,direct2diffuse,diffusepercent1,tau_T,trans_down,trans_up,diffusepercent2,halfHemiInt/(2.0*M_PI));
			fflush(stdout);
		}

	}

			

	return;
}		


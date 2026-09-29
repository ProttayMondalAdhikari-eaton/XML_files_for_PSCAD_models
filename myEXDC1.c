VERSION:
3.001

#include "myEXDC1.h"

STATIC:

	double dt,washoutG,wBase;
	double Vmax, Vref, tmp;
	double efldi,Verr,v_init;
	double mon1,mon2,mon3,mon4,mon5,mon6,mon7,mon8,mon9,mon10;
	int mon_array[10] = {1,2,3,4,5,6,7,8,9,10};

#include <builtin_gcc_generatorStructures.h>
	
	struct realpolex realpole1;
	struct realpolex realpole2;
	struct integratorx integrator1;
	struct washoutx washout1;
	struct leadlagx leadlag1;
	struct satKEx satKE1;
	
RAM_FUNCTIONS:
	static void initIntegrator(double T,double Xic,double Yic,struct integratorx* intG)
	{
		double dt = getTimeStep();
		intG->time_const = dt / (2.0 * T);
		intG->x_old = Xic;
		intG->y_old = Yic;
	}

	static void initWashout(double G, double T, double Xic, struct washoutx* wshout)
	{
		double K1;
		double dt = getTimeStep();

		if (T <= 0.0)
		{
			printf("\nError: Washout Time Constant must be greater than 0.0.");
			reportError("initWashout ",2);
		}

		wshout->x_old 	= Xic * G;
		wshout->y_old = 0.0;
		K1 	= 1.0/(1.0+dt/(2.0*T));	/* K1= 1/(1+dt/2T)*/
		wshout->K1 = K1;
		wshout->K2 = K1*(1.0-dt/(2.0*T)); /* K2= (1-dt/2T) / (1+dt/2T)*/
	}
	
	static void initRealPole(double G,double T,double Xic,double Yic,struct realpolex* rlpole)
	{
		double K;
		double dt = getTimeStep();

		rlpole->x_old	= Xic * G;
		rlpole->y_old = Yic;

		// Approximate cases with small time-constants as pure gains.
		if(T < dt)
			rlpole->pureGain = 1;
		else
			rlpole->pureGain = 0;

		K = (2.0 * T) / dt;

		rlpole->K1 = 1.0/(1.0+K);
		rlpole->K2 = (K-1.0)/(K+1.0);
	}

	static void initLeadLag(double Tb,double Tc,double Xic,double Yic,struct leadlagx* ldlag)
	{
		double x1,x2,x3,x4;
		double dt = getTimeStep();

			x1 = 1.0/(2.0*Tb/dt + 1.0);
			x2 = 2.0*Tc/dt + 1.0;
			x3 = 2.0*Tc/dt - 1.0;
			x4 = 2.0*Tb/dt - 1.0;
		
			ldlag->K1= x2*x1;
			ldlag->K2= x3*x1;
			ldlag->K3= x4*x1;

			ldlag->x_old = Xic;
			ldlag->y_old = Yic;
	}

	static double initSatKE (double E1,double SE1,double E2, double SE2,double Efd, double Ke, int Cal, struct satKEx* satKEx)
	{
	   double A,B,X,Se;

	/* Compute constants A and B from E1,Se1,E2,Se2 */
		X= sqrt( (SE1*E1)/(SE2*E2)  );
		if (Cal == 1)
		{
			A= ((X*E2 - E1) / (X - 1.0));
		}
		else
		{
			A= fabs ((X*E2 - E1) / (X - 1.0));
		}
		B= (SE2*E2)/((E2-A)*(E2-A));

		satKEx->A = A;
		satKEx->B = B;
		satKEx->Ke = Ke;

		
	/* Compute Se Function: Se= Efd* [B*(Efd-A)**2/Efd] */
		Efd= fabs(Efd);
		Se= 0.0;
		if (Cal == 1)
		{
			if (Efd > A && Efd != 0.0) Se= B*(Efd-A)*(Efd-A)/Efd;
		}
		else
		{
			if (Efd > fabs(A)) Se= B*(Efd-A)*(Efd-A)/Efd;
		}
		Se= Efd*Se + Efd*Ke;

		return Se;
	}

	static double initSatKE0 (double E1,double SE1,double E2, double SE2,double Efd, double Vrmax, int Cal, struct satKEx* satKEx)
	{

		double A,B,X,Se,Sat0,Ke_,Ke;
		double SXTMSX,FA,FB,FC;

	/* Compute constants A and B from E1,Se1,E2,Se2 */
		X= sqrt( (SE1*E1)/(SE2*E2)  );
		if (Cal == 1)
		{
			A= ((X*E2 - E1) / (X - 1.0));
		}
		else
		{
			A= fabs((X*E2 - E1) / (X - 1.0));
		}     
		B= (SE2*E2)/((E2-A)*(E2-A));

		satKEx->A = A;
		satKEx->B = B;

		if (SE1 == 0.0)
		{
			Sat0= 0.0;
		}
		else
		{
			SXTMSX = SE2 * E2;
			FC = SXTMSX/(SE1*E1);
			FA = FC - 1.0;
			FB = 2.0 * ( FC*E1 - E2 );
			FC = FC * E1*E1 - E2*E2;
			FA = (FB-sqrt(FB*FB-4.0*FA*FC))/(2.0*FA);
				
			if (Efd <= FA)
				Sat0= 0.0;
			else
				Sat0 =(((Efd-FA)*(Efd-FA))*SXTMSX/(Efd*((E2-FA)*(E2-FA))));
		}

		Ke_ = Vrmax/(10.0*Efd) - Sat0;

		satKEx->Ke = Ke_;	

		
	/* Compute Se Function: Se= Efd* [B*(Efd-A)**2/Efd] */
		Efd= fabs(Efd);
		Se= 0.0;
		
		if (Cal == 1)
		{
			if (Efd > A && Efd != 0.0) Se= B*(Efd-A)*(Efd-A)/Efd;
		}
		else
		{
			if (Efd > fabs(A)) Se= B*(Efd-A)*(Efd-A)/Efd;
		}
		Se= Efd*Se + Efd*Ke_;

	   return Se;
	}

	static double computeInitVcomp (char *Gen)
	{
		double Rpu,Xpu,Pinit,Qinit,genMVA,Ppu,Qpu,Vi;
		double Ireal,Iimag, Vreal, Vimag, Vre, Vim, Vcinit;
		double type;
		Vcinit = 1.0;
			type = get_data_compType(Gen,"CompType");
			if (type==1.0)
			{
			Vi = get_data(Gen,"Vpu");
			Rpu = get_data(Gen,"VcompRpu");
			Xpu = get_data(Gen,"VcompXpu");
			Pinit = get_data(Gen,"P0");
			Qinit = get_data(Gen,"Q0");
			genMVA = get_data(Gen,"MVA");
			Ppu = Pinit/genMVA;
			Qpu = Qinit/genMVA;
			Ireal = Ppu/Vi;
			Iimag = -Qpu/Vi;  
			Vreal = Vi + (Rpu*Ireal) - (Xpu*Iimag); 
			Vimag = Rpu*Iimag + Xpu*Ireal;
			Vre = Vreal*Vreal;
			Vim = Vimag*Vimag;
			Vcinit = sqrt(Vre + Vim);
			}
			else if(type==2.0)
			{
				Vcinit = get_data(Gen,"rBusVoltage");
			}
			else
			{
				Vcinit = 1.0;
			}
			return Vcinit;
	}

RAM:

	dt= getTimeStep();
	wBase = 2.0*PI*HTZ;

// Calculation of vrmax value.
	if(Vrmax <= 0.0)
	{
		if(E2< E1)
		{
			printf("\nWarning: Vrmax is less than zero.");
			printf("\n         E2 must be greater than E1");
			printf("\n         Subsystem %d, Machine %s, exciter type ESDC1A.\n", getSubsystem(),Gen);
			reportError("EXDC1",1);
		}
		Vrmax = (Se2 + Ke)*E2;	
		Vrmin = -1*Vrmax;
	}

	if (Vrmax < Vrmin)
	{
	tmp = Vrmin;
	Vrmin = Vrmax;
	Vrmax = tmp;
	printf("\nWarning: Vrmax is less than Vrmin.");
	printf("\n                Vrmin set to Vrmax and Vrmax set to Vrmin");
	printf("\n                Subsystem %d, Machine %s, exciter type EXDC1.\n", getSubsystem(),Gen);
	reportError("EXDC1",1);
	}

/* read initial field voltage */
	efldi = get_data(Gen,"Efldi");

/* Compute initial Vref value */
	if (Vi < 0.0 || LFInit==1)
	{
		Vi= get_data(Gen,"Vpu");
		if (LDComp == 1)
		{
			Vi = computeInitVcomp(Gen);
		}
	}
/* Compute initial Se + KE */
	if (Ke <= 0.0)
	{
		v_init = initSatKE0(E1,Se1,E2,Se2,efldi,Vrmax,Cal,&satKE1);
	}
	else
	{
		v_init = initSatKE(E1,Se1,E2,Se2,efldi,Ke,Cal,&satKE1);
	}
     
      Verr = v_init/Ka;	
      Vref = Verr + Vi;

	Vmax= 1.2;
	if (Vref  > Vmax) Vmax= 1.2*Vref;

/* the create_slider function is in the .h file in the INPUTS section */

/* Initialize variables for real pole #1 */

	initRealPole(1.0,Tr,Vi,Vi,&realpole1);

/* Initialize variables for lead - lag */

	initLeadLag(Tb,Tc,Verr,Verr,&leadlag1);
 
/* Initialize variables for real pole #2 */

	initRealPole(Ka,Ta,Verr,v_init,&realpole2);

/* Initialize variables for wash out */
	washoutG = Kf/Tf1;
	initWashout(washoutG,Tf1,efldi,&washout1);

/* Initialize variables for integrator */

	initIntegrator(Te,0.0,efldi,&integrator1);
	
	
/*  Check initial output values are within Vrmax and Vrmin limits */
	if (v_init > Vrmax || v_init < Vrmin)
	{
		printf("\nWarning: Initial value of mon6 is not within the entered limits.");
		printf("\n VRmax = %f, VRmin = %f, mon6 initial value = %f.",Vrmax,Vrmin,v_init);
		printf("\n Subsystem %d, Machine %s, exciter type EXDC1x.\n", getSubsystem(),Gen);
		reportError("EXDC1x",1);
	}	
	
/* Initialize all monitoring variables that are going to be referenced in the code section before they are assigned a value */
	mon9 = v_init;
	mon10 = 0.0;

CODE_FUNCTIONS:
#include <builtin_gcc.h>

CODE:

	double puSpeed;

/* Include load compensation */
	if (LDComp==1)
	{
		mon1 = VC;
	}
	else
	{
		mon1 = VPU;
	}
	
/* realpole #1*/
	mon2 = realPole(mon1, 1.0, 0, 0.0, 0.0,0,0,0.0, &realpole1); 

/* Summing junction */
	mon3 = Vref1 - mon2;

/* Optionally include stabilizer input */
	if (PSS == 1)
	{
		mon3 = mon3 + VS;
	}
	mon4 = mon3 - mon10; 

/* lead-lag */
	mon5 = leadLag(mon4,1.0,0,0.0,0.0,0,0,0.0,&leadlag1);

/* realpole #2 */
	mon6 = realPole(mon5,Ka,1,Vrmin,Vrmax,0,0,0.0,&realpole2);

/* Summing Junction */
	mon7 = mon6 - mon9; 

/* Integrator */
	mon8 = integrator(mon7,0,0.0,0.0,0,0,0.0,&integrator1); 

/* Saturation Function */
	mon9 = satKE(mon8, &satKE1); 

/* Washout */
	mon10 = washOut(mon8,washoutG,0,0.0,0.0,0,0,0.0,&washout1); 

/* Efd */
	OUT = mon8;
    if (spdMult == 1)
    {
		puSpeed = Speed/wBase;
		OUT = mon8 * puSpeed; 
	}

/* check if monitoring is enabled, if so assign the output based on the dial position */

	if (Mon ==1)
	{
		if (dial1 == 2)
		{
			internalV = mon2;
		}
		else if (dial1 == 3)
		{
			internalV = mon3;
		}
		else if (dial1 == 4)
		{
			internalV = mon4;
		}
		else if (dial1 == 5)
		{
			internalV = mon5;
		}
		else if (dial1 == 6)
		{
			internalV = mon6;
		}
		else if (dial1 == 7)
		{
			internalV = mon7;
		}
		else if (dial1 == 8)
		{
			internalV = mon8;
		}
		else if (dial1==9)
		{
			internalV = mon9;
		}
		else if (dial1==10)
		{
			internalV = mon10;
		}
		else
		{
			internalV = mon1;
		}
	} 






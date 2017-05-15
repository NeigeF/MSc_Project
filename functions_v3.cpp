#include "includeall_v3.h"

// Valid on 01/11/2016
// Modified on 24/01/2017 (add Kerr horizon & Slim disk)

// Functions to be used in the main program

// Schwarzchild radius
double RSCH(parameters* p) {
	double m1 = p->m1;
	double spin = p->spin;
	double G = 6.67e-11;		// gravitational constant (SI)
	double c = 299792458;		// speed of light (SI)
	double Ms = 1.989e30; 		// solar mass in kg
//	return 2*G*m1*Ms*1e-8 /(c*c);	 // in 10^10 cm
	return ( G*m1*Ms/(c*c) + sqrt((G*m1*Ms/(c*c))*(G*m1*Ms/(c*c))-spin*spin) )*1e-10*1e2; // in 10^10cm
	}

// ISCO Inner most stable circular orbit
double ISCO(parameters* p){
	double spin = p->spin;
	double G = 6.67e-11;		// gravitational constant (SI)
	double c = 299792458;		// speed of light (SI)
	double Ms = 1.989e30; 		// solar mass in kg
	double m = p->m1*G*Ms/(c*c);	// mass in GM/c^2 units
	double A1 = 1 + pow(1-spin*spin/(m*m),0.33333)*(pow(1+spin/m,0.33333) + pow(1-spin/m,0.3333));
	double A2 = sqrt(3*spin*spin/(m*m) + A1*A1);
	return m*(3+A2-sqrt((3 - A1)*(3 + A1 + 2*A2)))*1e-10*1e2;
}

// Eddington accretion rate
double dmed(double m1)	{			// m1 given in solar mass
	double G = 6.67e-11;			// gravitational constant (SI)
	double sigma = 6.7e-29;			// Thomson cross section (SI)
	double eta = 0.1;			// accretion efficiency
	double mp = 1.627e-27;			 // mass proton (kg)
	double c = 299792458;			// speed of light (SI)
	double pi = M_PI;
	double Ms = 1.989e30; 			// solar mass in kg
	double dm = 4*pi*G*m1*Ms*mp/(sigma*eta*c);//(kg/s)
	dm = dm * 1e-13 ;			 // (1e16g/s)
	return dm;
	}


/*****************************************************************************************************************************/
//                         SHAKURA SUNYAEV SOLUTION THIN DISK
/****************************************************************************************************************************/

// factor less than 1 closer to the center
double compute_f(double R, double R_star){
	double ff;
	ff = 1-sqrt(R_star/R);
	ff = pow(ff,.25);
	return ff;
	}
// surface density (g/cm^2)
double compute_density(double r, parameters* p, double f) {
	double alpha = p->alpha;
	double m1 = p-> m1;
	double dm = p->dm; //*dmed(m1);		
	return ( 5.2*pow(alpha, -0.8)*pow(dm, .7)*pow(m1, .25)*pow(r, -.75)*pow(f, 2.8) );
	};
// H = thickness (cm)
double compute_H (double r, parameters*p, double f) {
	double alpha = p->alpha;		
	double m1 = p-> m1;
	double dm = p->dm;//*dmed(m1);
	return ( 1.7e8 * pow(alpha,-.1)*pow(dm, .15)*pow(m1, -0.375)*pow(r, 1.125) * pow(f,0.6) );
	};
// density (g/cm3)
double compute_rho (double r, parameters* p, double f){
	double alpha = p->alpha;
	double m1 = p-> m1;
	double dm = p->dm;//*dmed(m1);	
	return (3.1e-8 * pow(alpha, -.7)*pow(dm, .55)*pow(m1, 0.625)*pow(r, -1.875) * pow(f,2.2));
	};
// temperature (K)
double compute_T (double r, parameters* p, double f){
	double alpha = p->alpha;
	double m1 = p-> m1;
	double dm = p->dm;//*dmed(m1);
	return ( 1.4e4 * pow(alpha,-0.2)*pow(dm, .3)*pow(m1, 0.25)*pow(r, -0.75) * pow(f,1.2) ) ;
	};
// optical depth
double compute_tau (parameters* p, double f){
	double alpha = p->alpha;
	double dm = p->dm;	
	double m1 = p-> m1;
	return(  190 * pow(alpha,-0.8)*pow(dm, .2)*pow(m1, 0.25)*pow(f,.8));
	};
// viscosity cm^2/s
double compute_nu (double r, parameters* p, double f){
	double alpha = p->alpha;
	double dm = p->dm;	
	double m1 = p-> m1;
	return ( 1.8e14 * pow(alpha,0.8)*pow(dm, .3)*pow(m1, -0.25)*pow(r, 0.75) * pow(f,1.2));
	};
// radial velocity (cm/s)
double compute_vr (double r, parameters* p, double f){
	double alpha = p->alpha;
	double dm = p->dm;	
	double m1 = p-> m1;
	return(- 2.7e4 * pow(alpha,0.8)*pow(dm, .3)*pow(m1, -0.25)*pow(r, -0.25) * pow(f,-2.8));
	};

// Radius (10^10cm) in function of temperature (K)
double RofT(parameters* p, double T){
	double alpha = p->alpha;
	double dm = p->dm;
	double m1 = p->m1;
	double R =  (1/1.4)*1e-4*pow(alpha,0.2)*pow(dm,-0.3)*pow(m1,-0.25)*T;
	R = pow(R,-1.33333);
	return R;
	}

// initialise state vector. Need a state vector with everything 0 (or else) exept the initial position R0.
void init_state(states*s, parameters*p) {
	double C = 1477; // = GMs/c2
	double alpha = p->alpha;
	double dm = p->dm;	
	double m1 = p-> m1;
	double R = s->R;
	double R_star = RSCH(p); // schwarzschild radius in 10^10cm
	double f = compute_f(R,R_star);
	double* mass_density = &s->mass_density;
	double* H = &s->H;
	double* rho = &s->rho;
	double* T = &s->T;
	double* tau = &s->tau;
	double* nu = &s->nu;
	*mass_density = compute_density(R,p,f);
	*H = compute_H(R,p,f);
	*rho = compute_rho(R,p,f);
	*T = compute_T(R,p,f);
	*tau = compute_tau(p,f);
	*nu = compute_nu(R,p,f);
//cout << f << endl;
	};

/***********************************************************************************************************/
// Yuan 2012 Hot accretion flow (scaling law)
/***********************************************************************************************************/

// Pressure(gcm-1s-2)
double compute_pressure_Yuan(double r, parameters* p){	 // r needs to be in scharwschild units
	double alpha = p->alpha;			 // alpha viscosity SS prescription
	double dm = p->dm; 				// need to be accretion onto the BH in eddington
	double m1 = p->m1; 				// mass of the BH in Solar masses
	double wind = p->wind; 				// wind parameter (0 for no wind)
	return(1.7e16*(1/alpha)*(1/m1)*dm*pow(r,wind-2.5));
	};

// sound speed squared (cm2s-2)
double compute_soundspeed2_Yuan(double r, parameters* p){ // r in Rsch units
	double alpha = p->alpha; 			// alpha viscosity SS prescription
	double dm = p->dm;				 // need to be accretion onto the BH in eddington
	double m1 = p->m1; 				// mass of the BH in Solar masses
	return 1.4e20*(1/r);
	};

// density (cs2 = P/rho) (g/cm3)
double compute_density_Yuan(double r, parameters* p){ 	// r in Rsch units
	double alpha = p->alpha; 			// alpha viscosity SS prescription
	double dm = p->dm; 				// need to be accretion onto the BH in eddington
	double m1 = p->m1; 				// mass of the BH in Solar masses
	return compute_pressure_Yuan(r,p)/compute_soundspeed2_Yuan(r,p);
	}

// temperature (K)
double compute_temperature_Yuan(double r, parameters* p){ // r in Rsch units
	double alpha = p->alpha; 			// alpha viscosity SS prescription
	double dm = p->dm;				// need to be accretion onto the BH in eddington
	double m1 = p->m1;				// mass of the BH in Solar masses
	double mp = 1.627e-27;				// mass proton (kg)
	double k = 1.38064852e-23;			// Boltzmann constant (J/K)
	double c = 299792458;				// light speed (m/s)	
	return mp*c*c/(6*k*r);
}
// radial velocity (cm/s)


/***********************************************************************************************************/
// Slim accretion disk Watarai & Fukue 1999
/***********************************************************************************************************/

double eps_prime(parameters* p){
	double gamma = p->gamma;
	double f_slim = p->f_slim;
	return (5./3. - gamma)/((gamma-1)*f_slim);
}

double g_slim(parameters *p){
	double alpha = p->alpha;
	double gamma = p->gamma;
	double f_slim = p->f_slim;
	double eps = eps_prime(p);
	return sqrt(1 + 18*alpha*alpha/((5+2*eps)*(5+2*eps))) - 1;
	}

// Free fall velocity (r is input in 10^10cm units but should be put back in m)
double v_kepl(parameters *p, double r){
	double mBH = p->m1;
	double G = 6.67e-11; 	// Gravitational constant in USI
	double Ms = 2e30; 	// Mass Sun in kg
	r = r*1e10*1e-2; 	// r from 10^10cm to m
	return sqrt(G*mBH*Ms/r); // keplerian velocity in m/s
}

// coefficients:
double c1(parameters *p){
	double alpha = p->alpha;
	double gamma = p->gamma;
	double f_slim = p->f_slim;
	double eps = eps_prime(p);
	return (5+2*eps)/(3*alpha*alpha)*g_slim(p);
//	return 3*alpha/(5+2*eps);
}

double c2(parameters *p){
	double alpha = p->alpha;
	double gamma = p->gamma;
	double f_slim = p->f_slim;
	double eps = eps_prime(p);
	return sqrt(2*eps*(5+2*eps)/(9*alpha*alpha)*g_slim(p));
}

double c3(parameters *p){
	double alpha = p->alpha;
	double gamma = p->gamma;
	double f_slim = p->f_slim;
	double eps = eps_prime(p);
	return 2*(5+2*eps)/(9*alpha*alpha)*g_slim(p);
//	return 2./(5+2*eps);
}

// Drift velocity
double vr_slim(parameters *p, double r){
	return -c1(p)*p->alpha*v_kepl(p,r); // in m/s
}

// Optical depth
double tau_slim(parameters *p, double r){
	double dm = p->dm/dmed(p->m1); // dm in units of Eddington limit
	return 1./(sqrt(2)*c1(p)*p->alpha)*dm*sqrt(RSCH(p)/r); // RSCH in 10^10cm and r in 10^10cm 
}

// Effective temperature
double Teff_slim(parameters* p, double r){
	double alpha = p->alpha;
	double gamma = p->gamma;
	double f_slim = p->f_slim;
	double m_BH = p->m1;
	r = r*1e10*1e-2; 			// r from 10^10cm to m
	double G = 6.67e-11;			// gravitational constant (SI)
	double sigma = 6.7e-29;			// Thomson cross section (SI)
	double sigma_boltzmann = 5.670367e-8;	// Stefan-Bolzmann constant (SI)
	double eta = 0.1;			// accretion efficiency
	double mp = 1.627e-27;			 // mass proton (kg)
	double c = 299792458;			// speed of light (SI)
	double pi = M_PI;
	double Ms = 1.989e30; 			// solar mass in kg
	double L_eddington = 4*pi*G*m_BH*Ms*mp*c/(sigma*eta);
	double T = pow((3*sqrt(c3(p))*L_eddington/(16*pi*sigma_boltzmann)),0.25)*pow(r,-.5);
//	double T = T_slim(p,r)/pow(tau_slim(p,r), 0.25);
//	cout << T << endl;
	return T;
}

// Temperature
double T_slim(parameters* p, double r){
	//return Teff_slim(p,r)*pow(tau_slim(p,r),0.25);
	double T;
	double alpha = p->alpha;
	double m_BH = p->m1;
	double dm = p->dm*1e16*1e-3;		// in kg/s
	r = r*1e10*1e-2; 			// r from 10^10cm to m
	double G = 6.67e-11;			// gravitational constant (SI)
	double sigma_boltzmann = 5.670367e-8;	// Stefan-Bolzmann constant (SI)
	double mp = 1.627e-27;			 // mass proton (kg)
	double c = 299792458;			// speed of light (SI)
	double pi = M_PI;
	double Ms = 1.989e30; 			// solar mass in kg
//cout << "alpha = " << alpha << "  m = " << m_BH << " r = " << r << "  dm = " << dm << endl;
//	T = pow((3*dm*c*sqrt(c3(p))*sqrt(G*m_BH*Ms/(r*r*r*r*r))/(16*pi*sigma_boltzmann*c1(p)*alpha)),0.25);
	T = G*m_BH*Ms*mp/(6*1.38e-23*r);
	return T;
}

// Scale height
double H_slim(parameters *p, double r){
	r = r*1e10*1e-2;					// r in m instead of 10^10 cm
	return sqrt(c3(p))*r;					// (m)
}

// Surface density
double surface_density_slim(parameters *p, double r){
	r = r*1e10*1e-2;					// r in m instead of 10^10 cm
	double G = 6.67e-11;					// gravitational constant (SI) 
	double Ms = 1.989e30; 					// solar mass in kg
	double alpha = p->alpha;
	double dm = p->dm*1e16*1e-3;				// dm is given in 10^16g/s, we convert it in kg/s 
	double m = p->m1*Ms;					// BH mass in kg 
	return dm/(2*M_PI)*1/(c1(p)*alpha)*1/sqrt(G*m*r); 	// (kg/m^2)
}		

// Density
double density_slim(parameters *p, double r){
	return surface_density_slim(p,r)/H_slim(p,r);
}

// Radiation pressure
double pressure_slim(parameters *p, double r){
	r = r*1e10*1e-2;
	double sigma_boltzmann = 5.670367e-8;			// Stefan-Bolzmann constant (SI)
	double Ms = 1.989e30; 					// solar mass in kg
	double c = 299792458;					// speed of light (SI)
	double m = p->m1*Ms;	
	return 4*sigma_boltzmann*2*sqrt(c3(p))*r/(c*3);
}

// initialise state vector. Need a state vector with everything 0 (or else) exept the initial position R0.
void init_state_slim(states*s, parameters*p) {
	double C = 1477; // = GMs/c2
	double alpha = p->alpha;
	double dm = p->dm;	
	double m1 = p-> m1;
	double R = s->R;
	double R_star = RSCH(p); // schwarzschild radius in 10^10cm
	double f = compute_f(R,R_star);
	double* mass_density = &s->mass_density;
	double* H = &s->H;
	double* rho = &s->rho;
	double* T = &s->T;
	double* tau = &s->tau;
	double* nu = &s->nu;
	*mass_density = surface_density_slim(p,R); // SURFACE DENSITY
	*H = H_slim(p,R);
	*rho = density_slim(p,R);
	*T = T_slim(p,R);			
//	*T = compute_T(R, p, 1);
//	*tau = compute_tau(p,f);
//	*nu = compute_nu(R,p,f);
//cout << f << endl;
	};



/***********************************************************************************************************/
// Integration functions
/***********************************************************************************************************/


// RK4 integrator for R
void RK4 (states* s, parameters* p, double f, double timestep) {
	double k1, k2, k3, k4;
	double* r = &s->R;
	double m1 = p->m1;
	double vr;
	double r1, f1;
	// first coef
	vr = compute_vr(*r,p,f);
	k1 = timestep*vr*1e-10;
//cout << vr <<" " << *r << endl;
	// second coef
	r1 = *r + k1/2;
	f1 = compute_f(r1,RSCH(p));
	vr = compute_vr(r1,p,f);
	k2 = timestep * vr*1e-10;
//cout <<  vr << " " << *r << endl;
	// third coef
	r1 = *r +k2/2;
	f = compute_f(r1,RSCH(p));
	vr = compute_vr(r1,p,f);
	k3 = timestep * vr*1e-10;
//cout <<vr << " " << *r << endl;
	// fourth coef
	r1 = *r + k3;
	f = compute_f(r1,RSCH(p));
	vr = compute_vr(r1,p,f);
	k4 = timestep * vr*1e-10;
//cout << vr << " " <<*r << endl;
	// result:
	*r = *r + (k1 + 2*k2 + 2*k3 + k4)/6;//e-10 because R is in unit R10/cm and vR is in cm/s
	}

void RK4_slim(states* s, parameters* p, double f, double timestep) {
	double k1, k2, k3, k4;
	double* r = &s->R;
	double m1 = p->m1;
	double vr;
	double r1, f1;
	// first coef
	vr = vr_slim(p,r1);
	k1 = timestep*vr*1e-10*1e2;
//cout << vr <<" " << *r << endl;
	// second coef
	r1 = *r + k1/2;
	vr = vr_slim(p,r1);
	k2 = timestep * vr*1e-10*1e2;
//cout <<  vr << " " << *r << endl;
	// third coef
	r1 = *r +k2/2;
	vr = vr_slim(p,r1);
	k3 = timestep * vr*1e-10*1e2;
//cout <<vr << " " << *r << endl;
	// fourth coef
	r1 = *r + k3;
	vr = vr_slim(p,r1);
	k4 = timestep * vr*1e-10*1e2;
//cout << vr << " " <<*r << endl;
	// result:
	*r = *r + (k1 + 2*k2 + 2*k3 + k4)/6;//1e10*1e-2 because R is in unit R10/cm and vR is in m/s
	}


/****************************************************************************************************************************/
// Nuclear physics (differentiation and integration of nuclear reactions)
/****************************************************************************************************************************/

// Compute the nuclear reaction rates (alpha network, incomplete 15/09/16)
void compute_alpha_rates(states* s, rates_alpha* rate){
	double T = s->T;
	T = T*1e-9;
	rate->r6 = rate6(T);
	rate->rv6 = rev6(T);
	rate->r7 = rate7(T);
	rate->rv7 = rev7(T);
	rate->r8 = rate8(T);
	rate->rv8 = rev8(T);
	rate->r9 = rate9(T);
	rate->rv9 = rev9(T);
	rate->r10 = rate10(T);
	rate->rv10 = rev10(T);
	rate->r11 = rate11(T);
	rate->rv11 = rev11(T);
	rate->r12 = rate12(T);
	rate->rv12 = rev12(T);
	rate->r13 = rate13(T);
	rate->rv13 = rev13(T);
	rate->r14 = rate14(T);
	rate->rv14 = rev14(T);
	rate->r15 = rate15(T);
	rate->rv15 = rev15(T);
	rate->r16 = rate16(T);
	rate->rv16 = rev16(T);
	rate->r17 = rate17(T);
	rate->rv17 = rev17(T);
	rate->rC12C12 = rateC12C12(T);
	rate->rC12C12a = rateC12C12(T);
	rate->rC12O16 = rateC12O16(T);
	rate->rC12O16a = rateC12O16(T);
	rate->rO16O16 = rateO16O16(T);
	rate->rO16O16a = rateO16O16(T);
}

// energy per unit time per unit mass (J/g/s)

double get_energy_He(double rho, LaVectorDouble compo, rates_alpha *r){
	// energy involved in each reaction in MeV
	double Q6 = 7.274;
	double Qv6 = -7.274;
	double Q7 = 1.943;
	double Qv7 = -1.943;
	double Q8 = 4.730;
	double Qv8 = -4.730;
	double Q9 = 9.316;
	double Qv9 = -9.316;
	double Q10 = 9.986;
	double Qv10 = -9.986;
	double Q11 = 6.949;
	double Qv11 = -6.949;
	double Q12 = 6.642;
	double Qv12 = -6.642;
	double Q13 = 7.042;
	double Qv13 = -7.041;
	double Q14 = 5.128;
	double Qv14 = -5.128;
	double Q15 = 7.694;
	double Qv15 = -7.694;
	double Q16 = 7.943;
	double Qv16 = -7.943;
	double Q17 = 8.001;
	double Qv17 = -8.001;
	// Avogadro constant
	double NA = 6.022140857e23;
	// 1MeV in J
	double unitfactor = 1.602176565e-19*1e6; 
	// Total energy, contribution from each forward reaction in alpha network 
	double Q  = 1./6.*rho*compo(0)*compo(0)*compo(0)*Q6*r->r6 + compo(0)*compo(1)*Q7*r->r7 + compo(0)*compo(2)*Q8*r->r8 + compo(0)*compo(3)*Q9*r->r9 
		+ compo(0)*compo(4)*Q10*r->r10 + compo(0)*compo(5)*Q11*r->r11 + compo(0)*compo(6)*Q12*r->r12 + compo(0)*compo(7)*Q13*r->r13
		+ compo(0)*compo(8)*Q14*r->r14 + compo(0)*compo(9)*Q15*r->r15 + compo(0)*compo(10)*Q16*r->r16 + compo(0)*compo(11)*Q17*r->r17;
	Q = Q*rho*NA;
	Q = Q*unitfactor;
	return Q;
}

double get_energy_CC(double rho, LaVectorDouble compo, rates_alpha *r){
	double QC12C12 = 13.931;
	double NA = 6.022140857e23;
	double unitfactor = 1.602176565e-19*1e6; // 1MeV in J
	double Q = compo(1)*compo(1)*QC12C12*r->rC12C12;	
	Q = Q*rho*NA;
	Q = Q*unitfactor;
	return Q;
}

double get_energy_CO(double rho, LaVectorDouble compo, rates_alpha *r){
	double QC12O16 = 16.754;
	double NA = 6.022140857e23;
	double unitfactor = 1.602176565e-19*1e6; // 1MeV in J
	double Q = compo(1)*compo(2)*QC12O16*r->rC12O16;	
	Q = Q*rho*NA;
	Q = Q*unitfactor;
	return Q;
}

double get_energy_OO(double rho, LaVectorDouble compo, rates_alpha *r){
	double QO16O16 = 16.541;
	double NA = 6.022140857e23;
	double unitfactor = 1.602176565e-19*1e6; // 1MeV in J
	double Q = compo(2)*compo(2)*QO16O16*r->rO16O16;	
	Q = Q*rho*NA;
	Q = Q*unitfactor;
	return Q;
}

// energy of reaction in J/s
double get_energy_rev(double rho, LaVectorDouble compo, rates_alpha *r){
	// energy involved in each reaction in MeV
	double Q6 = 7.274; 
	double Qv6 = -7.274;
	double Q7 = 1.943;
	double Qv7 = -1.943;
	double Q8 = 4.730;
	double Qv8 = -4.730;
	double Q9 = 9.316;
	double Qv9 = -9.316;
	double Q10 = 9.986;
	double Qv10 = -9.986;
	double Q11 = 6.949;
	double Qv11 = -6.949;
	double Q12 = 6.642;
	double Qv12 = -6.642;
	double Q13 = 7.042;
	double Qv13 = -7.041;
	double Q14 = 5.128;
	double Qv14 = -5.128;
	double Q15 = 7.694;
	double Qv15 = -7.694;
	double Q16 = 7.943;
	double Qv16 = -7.943;
	double Q17 = 8.001;
	double Qv17 = -8.001;
	// Avogadro constant
	double NA = 6.022140857e23;
	// 1MeV in J
	double unitfactor = 1.602176565e-19*1e6;
	// Total energy, contribution from each reverse reaction in alpha network 
	double Q = compo(1)*Qv6*r->rv6 + compo(2)*Qv7*r->rv7 + compo(3)*Qv8*r->rv8 + compo(4)*Qv9*r->rv9 + compo(5)*Qv10*r->rv10 
		+ compo(6)*Qv11*r->rv11 + compo(7)*Qv12*r->rv12 + compo(8)*Qv13*r->rv13 + compo(9)*Qv14*r->rv14 
		+ compo(10)*Qv15*r->rv15 + compo(11)*Qv16*r->rv16 + compo(12)*Qv17*r->rv17;
	Q = Q*NA;
	Q = Q*unitfactor;
	return Q;
}



// Rate of change of reacting nuclear species in the disk with alpha chain
void rate_of_change_alpha(states* s, rates_alpha* r, LaVectorDouble compo, LaVectorDouble& deriv_compo_expl){
	double rho = s->rho;
	deriv_compo_expl(0) = -3./6.*compo(0)*compo(0)*compo(0)*rho*rho*r->r6 + 3*compo(1)*r->rv6 - compo(0)*compo(1)*rho*r->r7 + compo(2)*r->rv7 - compo(0)*compo(2)*rho*r->r8 + compo(3)*r->rv8 - compo(0)*compo(3)*rho*r->r9 + compo(4)*r->rv9 - compo(0)*compo(4)*rho*r->r10 + compo(5)*r->rv10 - compo(0)*compo(5)*rho*r->r11 + compo(6)*r->rv11 - compo(0)*compo(6)*rho*r->r12 + compo(7)*r->rv12 - compo(0)*compo(7)*rho*r->r13 + compo(8)*r->rv13 - compo(0)*compo(8)*rho*r->r14 + compo(9)*r->rv14 - compo(0)*compo(9)*rho*r->r15 + compo(10)*r->rv15 - compo(0)*compo(10)*rho*r->r16 + compo(11)*r->rv16 - compo(0)*compo(11)*rho*r->r17 + compo(12)*r->rv17 + 0.5*compo(1)*compo(1)*rho*r->rC12C12a + compo(1)*compo(2)*rho*r->rC12O16a + 0.5*compo(2)*compo(2)*rho*r->rO16O16a;
	deriv_compo_expl(1) =  compo(0)*compo(0)*compo(0)*rho*rho*r->r6/6. - compo(1)*r->rv6 - compo(0)*compo(1)*rho*r->r7 + compo(2)*r->rv7 - compo(1)*compo(1)*rho*r->rC12C12 - compo(1)*compo(2)*rho*r->rC12O16;
	deriv_compo_expl(2) =  compo(0)*compo(1)*rho*r->r7 - compo(2)*r->rv7 - compo(0)*compo(2)*rho*r->r8 + compo(3)*r->rv8 - compo(1)*compo(2)*rho*r->rC12O16 - compo(2)*compo(2)*rho*r->rO16O16;
	deriv_compo_expl(3) =  compo(0)*compo(2)*rho*r->r8 - compo(3)*r->rv8 - compo(0)*compo(3)*rho*r->r9 + compo(4)*r->rv9 + 0.5*compo(1)*compo(1)*rho*r->rC12C12a;
	deriv_compo_expl(4) =  compo(0)*compo(3)*rho*r->r9 - compo(4)*r->rv9 - compo(0)*compo(4)*rho*r->r10 + compo(5)*r->rv10 + compo(1)*compo(2)*rho*r->rC12O16a;
	deriv_compo_expl(5) =  compo(0)*compo(4)*rho*r->r10 - compo(5)*r->rv10 - compo(0)*compo(5)*rho*r->r11 + compo(6)*r->rv11 + 0.5*compo(2)*compo(2)*rho*r->rO16O16a;
	deriv_compo_expl(6) =  compo(0)*compo(5)*rho*r->r11 - compo(6)*r->rv11 - compo(0)*compo(6)*rho*r->r12 + compo(7)*r->rv12;
	deriv_compo_expl(7) =  compo(0)*compo(6)*rho*r->r12 - compo(7)*r->rv12 - compo(0)*compo(7)*rho*r->r13 + compo(8)*r->rv13;
	deriv_compo_expl(8) =  compo(0)*compo(7)*rho*r->r13 - compo(8)*r->rv13 - compo(0)*compo(8)*rho*r->r14 + compo(9)*r->rv14;
	deriv_compo_expl(9) =  compo(0)*compo(8)*rho*r->r14 - compo(9)*r->rv14 - compo(0)*compo(9)*rho*r->r15 + compo(10)*r->rv15;
	deriv_compo_expl(10) = compo(0)*compo(9)*rho*r->r15 - compo(10)*r->rv15 - compo(0)*compo(10)*rho*r->r16 + compo(11)*r->rv16;
	deriv_compo_expl(11) = compo(0)*compo(10)*rho*r->r16 - compo(11)*r->rv16 - compo(0)*compo(11)*rho*r->r17 + compo(12)*r->rv17;
	deriv_compo_expl(12) = compo(0)*compo(11)*rho*r->r17 - compo(12)*r->rv17;
};

// Rate of change of reacting nuclear species in the disk with alpha chain
void rate_of_change_alpha_simple(states* s, rates_alpha* r, LaVectorDouble compo, LaVectorDouble& deriv_compo_expl){
	double rho = s->rho;
	deriv_compo_expl(0) = -3./6.*compo(0)*compo(0)*compo(0)*rho*rho*r->r6  - compo(0)*compo(1)*rho*r->r7 - compo(0)*compo(2)*rho*r->r8;
	deriv_compo_expl(1) =  compo(0)*compo(0)*compo(0)*rho*rho*r->r6/6. - compo(0)*compo(1)*rho*r->r7;
	deriv_compo_expl(2) =  compo(0)*compo(1)*rho*r->r7;
	deriv_compo_expl(3) =  compo(0)*compo(2)*rho*r->r8 - compo(3)*r->rv8 - compo(0)*compo(3)*rho*r->r9 + compo(4)*r->rv9 + 0.5*compo(1)*compo(1)*rho*r->rC12C12a;
};


// Compute the mass fraction X from the nucleon fraction Y (= compo). X = Y*A where A is the number of nucleons
void compute_mass_fraction(LaVectorDouble Y, LaVectorDouble& X){
	X(0) = Y(0)*4;		X(1) = Y(1)*12;		X(2) = Y(2)*16;		X(3) = Y(3)*20; 
	X(4) = Y(4)*24;		X(5) = Y(5)*28; 	X(6) = Y(6)*32; 	X(7) = Y(7)*36;
	X(8) = Y(8)*40;		X(9) = Y(9)*44; 	X(10) = Y(10)*48;	X(11) = Y(11)*52;
	X(12) = Y(12)*56;
};
	
// choose timestep WRONG DO NOT CONSIDER ELEMENTS WHICH ARE PRODUCED !!!! (division by zero)  TO FIX LATER IF TIME
//void adapt_timestep(double* dt, LaVectorDouble compo, LaVectorDouble deriv_compo_expl){
//	for (double i = 0; i<3; ++i){
//		double time_scale = abs(compo(i)/deriv_compo_expl(i))/10;
//cout << "timescale = " << time_scale << endl;
//		if ( time_scale < *dt) {*dt = time_scale;} 
//	}
//};
	
// Linearisation: Matrix to be inverted 

void linearize_alpha(states* s, LaVectorDouble compo, rates_alpha* r, double timestep, LaGenMatDouble& A){
	double rho = s->rho;
	A(0,0) = 1 + (3*3./6.*compo(0)*compo(0)*rho*rho*r->r6 + compo(1)*rho*r->r7 + compo(2)*rho*r->r8 + compo(3)*rho*r->r9 + compo(4)*rho*r->r10 + compo(5)*rho*r->r11 + compo(6)*rho*r->r12 + compo(7)*rho*r->r13 + compo(8)*rho*r->r14 + compo(9)*rho*r->r15 + compo(10)*rho*r->r16 + compo(11)*rho*r->r17)*timestep;
	A(0,1) = compo(0)*rho*r->r7*timestep - 3*r->rv6*timestep - compo(1)*rho*r->rC12C12a*timestep - compo(2)*rho*r->rC12O16a*timestep;
	A(0,2) = compo(0)*rho*r->r8*timestep - r->rv7*timestep - compo(1)*rho*r->rC12O16a*timestep - compo(2)*rho*r->rO16O16a*timestep;
	A(0,3) = compo(0)*rho*r->r9*timestep - r->rv8*timestep;
	A(0,4) = compo(0)*rho*r->r10*timestep -r->rv9*timestep;
	A(0,5) = compo(0)*rho*r->r11*timestep - r->rv10*timestep;
	A(0,6) = compo(0)*rho*r->r12*timestep - r->rv11*timestep;
	A(0,7) = compo(0)*rho*r->r13*timestep - r->rv12*timestep;
	A(0,8) = compo(0)*rho*r->r14*timestep - r->rv13*timestep;
	A(0,9) = compo(0)*rho*r->r15*timestep - r->rv14*timestep;
	A(0,10) = compo(0)*rho*r->r16*timestep - r->rv15*timestep;
	A(0,11) = compo(0)*rho*r->r17*timestep - r->rv16*timestep;
	A(0,12) = -r->rv17*timestep;
	A(1,0) = (-3./6.*compo(0)*compo(0)*rho*rho*r->r6 + compo(1)*rho*r->r7)*timestep;
	A(1,1) = 1 + compo(0)*rho*r->r7*timestep + r->rv6*timestep + 2*compo(1)*rho*r->rC12C12*timestep + compo(2)*rho*r->rC12O16*timestep;
	A(1,2) = -r->rv7*timestep + compo(1)*rho*r->rC12O16*timestep;
	A(1,3) = 0;
	A(1,4) = 0;
	A(1,5) = 0;
	A(1,6) = 0;
	A(1,7) = 0;
	A(1,8) = 0;
	A(1,9) = 0;
	A(1,10) = 0;
	A(1,11) = 0;
	A(1,12) = 0;
	A(2,0) = (-compo(1)*rho*r->r7 + compo(2)*rho*r->r8)*timestep;
	A(2,1) = -compo(0)*rho*r->r8*timestep + compo(2)*rho*r->rC12O16*timestep;
	A(2,2) = 1+compo(0)*rho*r->r8*timestep + r->rv7*timestep + compo(1)*rho*r->rC12O16*timestep + 2*compo(2)*rho*r->rO16O16*timestep;
	A(2,3) = -r->rv8*timestep;
	A(2,4) = 0;
	A(2,5) = 0;
	A(2,6) = 0;
	A(2,7) = 0;
	A(2,8) = 0;
	A(2,9) = 0;
	A(2,10) = 0;
	A(2,11) = 0;
	A(2,12) = 0;
	A(3,0) = (-compo(2)*rho*r->r8+compo(3)*rho*r->r9)*timestep;
	A(3,1) = -compo(1)*rho*r->rC12C12a*timestep;
	A(3,2) = -compo(0)*rho*r->r8*timestep;
	A(3,3) = 1+compo(0)*rho*r->r9*timestep + r->rv8*timestep;
	A(3,4) = -r->rv9*timestep;
	A(3,5) = 0;
	A(3,6) = 0;
	A(3,7) = 0;
	A(3,8) = 0;
	A(3,9) = 0;
	A(3,10) = 0;
	A(3,11) = 0;
	A(3,12) = 0;
	A(4,0) = (-compo(3)*rho*r->r9+compo(4)*rho*r->r10)*timestep;
	A(4,1) = -compo(2)*rho*r->rC12O16a*timestep;
	A(4,2) = -compo(1)*rho*r->rC12O16a*timestep;
	A(4,3) = -compo(0)*rho*r->r9*timestep;
	A(4,4) = 1+compo(0)*rho*r->r10*timestep + r->rv9*timestep;
	A(4,5) = -r->rv10*timestep;
	A(4,6) = 0;
	A(4,7) = 0;
	A(4,8) = 0;
	A(4,9) = 0;
	A(4,10) = 0;
	A(4,11) = 0;
	A(4,12) = 0;
	A(5,0) = (-compo(4)*rho*r->r10+compo(5)*rho*r->r11)*timestep;
	A(5,1) = 0;
	A(5,2) = -compo(2)*rho*r->rO16O16a*timestep;
	A(5,3) = 0;
	A(5,4) = -compo(0)*rho*r->r10*timestep;
	A(5,5) = 1+compo(0)*rho*r->r11*timestep + r->rv10*timestep;
	A(5,6) = -r->rv11*timestep;
	A(5,7) = 0;
	A(5,8) = 0;
	A(5,9) = 0;
	A(5,10) = 0;
	A(5,11) = 0;
	A(5,12) = 0;
	A(6,0) = (-compo(5)*rho*r->r11+compo(6)*rho*r->r12)*timestep;
	A(6,1) = 0;
	A(6,2) = 0;
	A(6,3) = 0;
	A(6,4) = 0;
	A(6,5) = -compo(0)*rho*r->r11*timestep;
	A(6,6) = 1+compo(0)*rho*r->r12*timestep + r->rv11*timestep;
	A(6,7) = -r->rv12*timestep;
	A(6,8) = 0;
	A(6,9) = 0;
	A(6,10) = 0;
	A(6,11) = 0;
	A(6,12) = 0;
	A(7,0) = (-compo(6)*rho*r->r12+compo(7)*rho*r->r13)*timestep;
	A(7,1) = 0;
	A(7,2) = 0;
	A(7,3) = 0;
	A(7,4) = 0;
	A(7,5) = 0;
	A(7,6) = -compo(0)*rho*r->r12*timestep;
	A(7,7) = 1+compo(0)*rho*r->r13*timestep + r->rv12*timestep;
	A(7,8) = -r->rv13*timestep;
	A(7,9) = 0;
	A(7,10) = 0;
	A(7,11) = 0;
	A(7,12) = 0;
	A(8,0) = (-compo(7)*rho*r->r13+compo(8)*rho*r->r14)*timestep;
	A(8,1) = 0;
	A(8,2) = 0;
	A(8,3) = 0;
	A(8,4) = 0;
	A(8,5) = 0;
	A(8,6) = 0;
	A(8,7) = -compo(0)*rho*r->r13*timestep;
	A(8,8) = 1+compo(0)*rho*r->r14*timestep + r->rv13*timestep;
	A(8,9) = -r->rv14*timestep;
	A(8,10) = 0;
	A(8,11) = 0;
	A(8,12) = 0;
	A(9,0) = (-compo(8)*rho*r->r14+compo(9)*rho*r->r15)*timestep;
	A(9,1) = 0;
	A(9,2) = 0;
	A(9,3) = 0;
	A(9,4) = 0;
	A(9,5) = 0;
	A(9,6) = 0;
	A(9,7) = 0;
	A(9,8) = -compo(0)*rho*r->r14*timestep;
	A(9,9) = 1+compo(0)*rho*r->r15*timestep + r->rv14*timestep;
	A(9,10) = -r->rv15*timestep;
	A(9,11) = 0;
	A(9,12) = 0;
	A(10,0) = (-compo(9)*rho*r->r15+compo(10)*rho*r->r16)*timestep;
	A(10,1) = 0;
	A(10,2) = 0;
	A(10,3) = 0;
	A(10,4) = 0;
	A(10,5) = 0;
	A(10,6) = 0;
	A(10,7) = 0;
	A(10,8) = 0;
	A(10,9) = -compo(0)*rho*r->r15*timestep;
	A(10,10) = 1+compo(0)*rho*r->r16*timestep + r->rv15*timestep;
	A(10,11) = -r->rv16*timestep;
	A(10,12) = 0;
	A(11,0) = (-compo(10)*rho*r->r16+compo(11)*rho*r->r17)*timestep;
	A(11,1) = 0;
	A(11,2) = 0;
	A(11,3) = 0;
	A(11,4) = 0;
	A(11,5) = 0;
	A(11,6) = 0;
	A(11,7) = 0;
	A(11,8) = 0;
	A(11,9) = 0;
	A(11,10) = -compo(0)*rho*r->r16*timestep;
	A(11,11) = 1+compo(0)*rho*r->r17*timestep + r->rv16*timestep;
	A(11,12) = -r->rv17*timestep;
	A(12,0) = (-compo(11)*rho*r->r17)*timestep;
	A(12,1) = 0;
	A(12,2) = 0;
	A(12,3) = 0;
	A(12,4) = 0;
	A(12,5) = 0;
	A(12,6) = 0;
	A(12,7) = 0;
	A(12,8) = 0;
	A(12,9) = 0;
	A(12,10) = 0;
	A(12,11) = -compo(0)*rho*r->r17*timestep;
	A(12,12) = 1 + r->rv17*timestep;

};






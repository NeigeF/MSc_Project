#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <sstream>
#include <laslv.h>
#include <lapackpp.h>
#include "include_rates.h"
using namespace std;

struct states {
	double mass_density;	// surface density (g/cm^2) BAD NAME!!!
	double H;		// scale height
	double rho;		// density (g/cm^3)
	double T;		// temperature (K)
	double tau;		// optical depth
	double nu;		// viscosity
	double R;		//  10^10 cm
};
struct parameters {
	double alpha;		// Shakura-Sunyaev prescription
	double dm;		// (10^16g/s)
	double m1;		// (Msun)
	double wind;		// wind prescription
	double gamma;		// calorific constant ratio
	double f_slim;		// fraction of energy advected
	double spin;		// black hole spin
};

struct rates_alpha {
	double r6;
	double rv6;
	double r7;
	double rv7;
	double r8;
	double rv8;
	double r9;
	double rv9;
	double r10;
	double rv10;
	double r11;
	double rv11;
	double r12;
	double rv12;
	double r13;
	double rv13;
	double r14;
	double rv14;
	double r15;
	double rv15;
	double r16;
	double rv16;
	double r17;
	double rv17;
	double rC12C12;	
	double rC12C12a;
	double rC12C12p;
	double rC12O16;
	double rC12O16a;
	double rO16O16;
	double rO16O16a;
};

double RSCH(parameters*);
double ISCO(parameters*);
double dmed(double);

// Shakura-Sunyaev solution (thin disk)
double compute_f(double, double);
double compute_density(double, parameters*, double);
double compute_H(double, parameters*, double );
double compute_rho(double, parameters*, double);
double compute_T(double, parameters*, double);
double compute_tau(double, parameters*, double);
double compute_nu(double, parameters*, double);
double compute_vr(double, parameters*, double);
double compute_P(double, parameters*, double);
double RofT(parameters*, double);

// Yuan 2014 solution (hot accretion flows)
double compute_pressure_Yuan(double, parameters*);
double compute_soundspeed2_Yuan(double, parameters*);
double compute_temperature_Yuan(double, parameters*);
double compute_density_Yuan(double, parameters*);

// Slim disk solution
double eps_prime(parameters* p);
double g_slim(parameters*);
double v_kepl(parameters*, double);
double c1(parameters*);
double c2(parameters*);
double c3(parameters*);
double vr_slim(parameters*, double);
double tau_slim(parameters*, double);
double Teff_slim(parameters*, double);
double T_slim(parameters*, double);
double H_slim(parameters*, double);
double surface_density_slim(parameters*, double);
double density_slim(parameters*, double);

// Integration
void init_state(states*, parameters*);
void init_state_slim(states*, parameters*);
void RK4(states*, parameters*, double, double);
void RK4_slim(states*, parameters*, double, double);

// Nuclear physics: differentiation and intgration
void compute_alpha_rates(states* , rates_alpha*);
void rate_of_change_alpha(states*, rates_alpha*, LaVectorDouble, LaVectorDouble&);
void rate_of_change_alpha_simple(states*, rates_alpha*, LaVectorDouble, LaVectorDouble&);
void linearize_alpha(states*, LaVectorDouble, rates_alpha*, double, LaGenMatDouble&);
void adapt_timestep(double*, LaVectorDouble, LaVectorDouble);
void compute_mass_fraction(LaVectorDouble, LaVectorDouble&);

// get the nuclear energy (J/g/s)
double get_energy_He(double, LaVectorDouble, rates_alpha *);
double get_energy_CC(double, LaVectorDouble, rates_alpha *);
double get_energy_CO(double, LaVectorDouble, rates_alpha *);
double get_energy_OO(double, LaVectorDouble, rates_alpha *);
double get_energy_rev(double, LaVectorDouble, rates_alpha *);


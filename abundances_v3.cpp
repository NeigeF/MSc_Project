// Steady alpha accretion disk (Shakura-Sunavey)
// 1st attempt to integrate the abundances in the disk with alpha chain network
// 15/09/2016
// 26/09/2016
// 28/11/2016 firsts runs doable non-rotating BHs, alpha chain
// 00/02/2017 add spinning BHs
// 26/02/2017 make input file to avoid recompiling
// 21/03/2017 fixed a bug on BH spin

# include "includeall_v3.h"

int main(int argc, char* argv[]) {

int save = 0;
cout << "Save in output files? yes 1 no 0 " << endl;
cin >> save;

// useful constants
double G = 6.67e-11;		// Gravitational constant (SI)
double c = 3e8;			// Speed of light (SI)
double Ms = 1.989e30; 		// Solar mass (kg)
double Rs = 695700000;		// Solar radius (m)

// Declare variables
double alpha;			// alpha parameter (less than 1)
double m1;			// mass of the central object (solar mass)
double dm;			// *dmed(m1); 	// accretion rate (10^16 g/s);
double spin;			// *G*m1*Ms/(c*c);	// BH spin
double m_WD;			// mass of the binary companion (solar mass)
double R_WD ;			// 10^10cm
int n;				// number of isotopes involved in the network
int i;				// index
double dt_max;			// maximum timestep allowed by user

// open, read, close input file
ifstream input_file;
input_file.open(argv[1]);
input_file >> m1 >> m_WD;
input_file >> dm >> alpha >> spin;
input_file >> n >> dt_max;

   LaVectorDouble compo(n);
   LaVectorDouble mass_fraction_initial(n);
for(i=0;i<n;++i){ input_file >> compo(i);}
    mass_fraction_initial = compo;
// mass fraction X to molar fraction Y:
compo(0) = compo(0)/4; compo(1) = compo(1)/12; compo(2) = compo(2)/16; compo(3) = compo(3)/20;
   
string s;
input_file >> s;
input_file.close();


// White dwarf mass
double M_WD = m_WD*Ms*1e3*1e-16;// mass of the binary companion (white dwarf) (10^16g)
dm = dm*dmed(m1);		
		
// Mass radius relation for WD (Rappaport 1988 ?)    // R_WD calculated with Verbunt&Rappaport(1988) eq(15) 
double Mp = 0.00057; 		// solar masses
double Mch = 1.44; 		// chandrasekkhar mass in solar masses
R_WD = 0.0114*sqrt(pow(m_WD/Mch,-0.666667)-pow(m_WD/Mch,0.66666667))*pow((1+3.5*pow(m_WD/Mp,-0.6666667)+(Mp/m_WD)),0.666667)*Rs*1e3*1e2*1e-10; //10^10cm

// Mass of elements in the disk
    LaVectorDouble mass(n);
    LaVectorDouble mass_fraction(n);

// Integration initialisation
double R_min, R_max;
double r, f, vr;		// radius (10^10 cm, relat correction, radial velocity)	
double t=0;			// time
double M = 0;			// mass processed in the disk
double M_disk; 			// total disk mass
double M_el = 0; 		// sum of the masses of each isotopes in the disk

// Nuclear power
double energy_He, energy_CC, energy_CO, energy_OO, energy_rev;

// viscous disspation rate and gravitational energy release
double energy_visc, energy_grav;


// initialise parameters of the system : computer
    parameters param;
    param.alpha = alpha;
    param.dm = dm;
    param.m1 = m1;
    param.spin = spin*(G*m1*Ms/(c*c));
    param.wind = 0;				// haven't managed to make Watarai (1999 or 2000) disk
    states state;
    rates_alpha rate_a;
    LaGenMatDouble deriv_linear(n,n);
    LaVectorDouble deriv_compo_expl(n);
    LaVectorDouble deriv_compo_impl(n);


// initialise output file of the run
stringstream ss,cc;
ss << s;			// Folder name (must exist already. type mkdir name in terminal)
cc << s;
ss << "spin-" << spin;
ss << "mdot-" << dm/dmed(m1);
ss << "_m-" << m1;
ss << ".txt";
cc << "results_spin-"<< spin << "_mdot-"<< dm/dmed(m1) << "_m-" <<m1 << ".txt";

FILE *output = fopen(ss.str().c_str(),"a");
FILE *results = fopen(cc.str().c_str(),"a");



// ========================================================================================================
// disk properties
// ========================================================================================================
double q = m_WD/m1;		// Mass ratio
// Roche lobe overflow : Eggleton (1983) (10^10cm)
double R_RLOF = R_WD*(pow(q,0.66667) + log(1+pow(q,0.3333)))/(0.48*pow(q,0.666667)); // orbital separation at roche contact
// Circularisation radius = outer disk (10^10cm) (calculation in logbook)
double R_disk = R_RLOF/(1+q);
// Total disk mass
// Integration zone in the disk:
R_min = RofT(&param, 1.5e10);		// in 10^10 cm
R_max = RofT(&param, 7e8);		// in 10^10 cm
if (R_min < ISCO(&param)){R_min = ISCO(&param);}
M_disk = 2*M_PI*3.1*1.7*1e20*1e-16/1.25 *pow(alpha,-0.8)*pow(dm,0.7)*pow(m1,0.25)*(pow(R_disk,1.25)-pow(R_min,1.25));

// check Print out
cout << endl;
cout << endl;
cout << "------------------------------------------------------------------------------"<< endl;
cout << "Disk properties " << endl;
cout << "------------------------------------------------------------------------------"<< endl;
cout << "Tmax= 1.5e10 K, Rmin = " << R_min <<" 10^10cm = " << R_min/RSCH(&param) << " Rsch "  << endl;
cout << "Tmin = 7e8 K, Rmax = " << R_max <<" 10^10cm = " << R_max/RSCH(&param) << " Rsch " << endl;
cout << "Disk radius = " << R_disk << " 10^10cm = "<< R_disk/RSCH(&param)  << " Rsch" << endl;
cout << "Total disk mass = " << M_disk << " 10^16g = " << M_disk/Ms << " Ms" << endl;
cout << "BH spin = " << param.spin << endl;
cout << "ISCO = " << ISCO(&param) << endl;


// =======================================================================================================
// INTEGRATE R(t) (AND STATE(t))
// =======================================================================================================
// TIMESTEP
double dt = dt_max;
int index = 0;

// Radius range around the BH //R_min = 3*RSCH(m1);	// in 10^10 cm //R_max = 5e2*RSCH(m1);	// in 10^10 cm

// Initialise state object
r = R_max; state.R = r; cout<< "r = "<< r << endl;
init_state(&state,&param);

f = compute_f(r,RSCH(&param));
vr = -compute_vr(r,&param,f);

// Calculate nuclear reaction rates
compute_alpha_rates(&state, &rate_a);

// compute variation elements
rate_of_change_alpha(&state, &rate_a, compo, deriv_compo_expl);

// Linearize (compute matrix to invert)
linearize_alpha(&state, compo, &rate_a, dt, deriv_linear); //cout << "deriv linear matrix" << endl; //cout << deriv_linear << endl;

// Solve Linear system
LaLinearSolve(deriv_linear, deriv_compo_impl, deriv_compo_expl*dt);
// New abundances
compo = compo + deriv_compo_impl;
compute_mass_fraction(compo, mass_fraction);

M = M + dm*dt; // in unit of 10^16g
mass = mass + 2*M_PI* mass_fraction * state.rho * r * state.H * vr * dt * 1e-16*1e10;

// Go through the disk
while (state.R > R_min){		
	// Integrate blob composition in the disk

	// Calculate nuclear reaction rates
	compute_alpha_rates(&state, &rate_a);

	// compute variation elements
	rate_of_change_alpha(&state, &rate_a, compo, deriv_compo_expl);

	// Linearize (compute matrix to invert)
	linearize_alpha(&state, compo, &rate_a, dt, deriv_linear);

	// Solve Linear system
	LaLinearSolve(deriv_linear, deriv_compo_impl, deriv_compo_expl*dt);

	// New abundances
	compo = compo + deriv_compo_impl;	// implicit scheme
//	compo = compo + deriv_compo_expl*dt;	// explicit scheme
	compute_mass_fraction(compo, mass_fraction);

	// energy per surface unit released per unit time : burning, viscous, gravitational
	 energy_He = get_energy_He(state.rho,compo,&rate_a)*state.mass_density;
	 energy_CC = get_energy_CC(state.rho,compo,&rate_a)*state.mass_density;
	 energy_CO = get_energy_CO(state.rho,compo,&rate_a)*state.mass_density;
	 energy_OO = get_energy_OO(state.rho,compo,&rate_a)*state.mass_density;
	 energy_rev = get_energy_rev(state.rho,compo,&rate_a)*state.mass_density;
	 energy_visc = (9./4.)*state.mass_density*1e-3*state.nu*1e-4*(G*m1*Ms)/(r*1e10*1e-2*r*1e10*1e-2*r*1e10*1e-2);
	 energy_grav = state.mass_density*1e-3*G*m1*Ms*(vr*1e-2)/((r*1e10*1e-2)*(r*1e10*1e-2));

	// Mass which has been transformed
	M = M + dm*dt; // in unit of 10^16g
	mass = mass + 2*M_PI* mass_fraction * state.rho * r * state.H * vr * dt*1e-16*1e10; // in 10^16g

	// Abundances are positive or zero: check for bugs
	//for (i = 0; i<n; ++i) { if (compo(i) < 0) { cout << i << " " << compo(i) <<"impl " << deriv_compo_impl(i)<< "expl " <<  endl;}}

	// Update disk state
	r = state.R;	
	f = compute_f(r,RSCH(&param));
	RK4(&state,&param,f,dt);
	init_state(&state,&param);
	vr = -compute_vr(r,&param,f);

	// Increment time
	t = t+dt;

	// Save results
//	if(index % 5 ==0){
// test mass fraction sum
M_el = 0;
for (i = 0; i<n; ++i){ M_el = M_el + 2*M_PI* mass_fraction(i) * state.rho * r * state.H * vr * dt*1e-16*1e10;}
//cout << M_el/(2*M_PI* state.rho * r * state.H * vr * dt*1e-16*1e10) << endl;
//	if(index % 10000 == 0) {cout << "t = " << t << "  n = " << index << " T = " << state.T << endl; }
	if (save==1){
	fprintf(output,"%.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \n ",t, r, state.H, state.rho, state.T, compo(0)*4, compo(1)*12, compo(2)*16, compo(3)*20, compo(4)*24, compo(5)*28, compo(6)*32, compo(7)*36, compo(8)*40, compo(9)*44, compo(10)*48, compo(11)*52, compo(12)*56, energy_He, energy_CC, energy_CO, energy_OO, energy_rev, energy_visc,energy_grav);
	}
//	}
	index++;
}

// Summed mass of elements in the disk
M_el = 0;
for (i = 0; i<n; ++i){ M_el = M_el + mass(i); }

// Total mass of elements in the whole disk (not only where integration took place):
// For newly produced elements: same
// For initial elments: different: need to add: initial mass fraction * mass outer disk

// mass outer disk:
double M_outer = 2*M_PI*3.1*1.7*1e20*1e-16/1.25 *pow(alpha,-0.8)*pow(dm,0.7)*pow(m1,0.25)*(pow(R_disk,1.25)-pow(R_max,1.25));
    LaVectorDouble total_mass(n);
    total_mass = mass + M_outer*mass_fraction_initial;

// mass innner disk:
if (ISCO(&param) < R_min){
	double M_inner = 2*M_PI*3.1*1.7*1e20*1e-16/1.25 *pow(alpha,-0.8)*pow(dm,0.7)*pow(m1,0.25)*(pow(R_disk,1.25)-pow(R_max,1.25));
	total_mass = total_mass + M_inner*mass_fraction;
	}

// Total disk mass normally is
double M222 = 0;
for (i = 0; i<n; ++i){ M222 = M222 + total_mass(i); }

cout << endl;
cout << endl;
cout << endl;
cout << endl;
cout << "-------------------------------------------------------------------------------" << endl;
cout << "RESULTS: M = " << m1 <<  " Ms, Mdot = " << dm << " 1e16g/s = " << dm/dmed(m1) << " Eddington "<< "M_WD = " << m_WD << " Ms " << endl;
cout << "-------------------------------------------------------------------------------" << endl;
cout << endl;	
//cout << "Mass processed in the disk: " << M << " 10^16g =  " << M/Ms/1e3*1e16 << " Msun " <<endl;
cout << endl;
cout << "Mass of each isotope in 10^16g in the burning disk fraction:" << endl;
cout << mass << endl;
cout << endl;
cout << "Total Mass of each isotope in 10^16g in the whole disk:" << endl;
cout << total_mass << endl;
cout << endl;
//cout << "Sum M_el = " << M_el << " In the whole disk = " << M222 << endl;
cout << "Final temperature " << state.T << "K, final density " << state.rho << " g/cm3" << endl;
cout << "Time of integration: " << t << " s" << endl;
cout << endl;
cout << endl;
cout << endl;
cout << endl;

// Write final results in a new file
if (save == 1){
fprintf(results, "%.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n %d \t %.6e \n",m1, m_WD, dm/dmed(m1), M_disk, M_disk/M_WD, R_max, R_min, R_disk, t, 4, total_mass(0), 12, total_mass(1), 16, total_mass(2), 20, total_mass(3), 24, total_mass(4), 28, total_mass(5), 32, total_mass(6), 36, total_mass(7), 40, total_mass(8), 44, total_mass(9), 48, total_mass(10), 52, total_mass(11), 56, total_mass(12)); 
}

}


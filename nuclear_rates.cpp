/************************************************
 		Main functions
************************************************/

// Continuously updated nuclear rates -- 17/11/2016
// Source: Angulo 1999 from rateO to rate9
// To do: check if find same result with rates from other sources
// Find other sources for heavier nuclei fusion


// To be copied and pasted in the other folders of work



# include "include_rates.h"

// Nuclear reaction rates -----------------------------------------------------------------------------------------------------------------------------------------

// 0. H+n -> h2 + gamma
double rate0 (double T9) {
	double r;
//	r = 44060*( 1 + 0.106597*sqrt(T9) - 2.75037*T9 + 4.62949*pow(T9, 1.5) - 3.52204 * T9*T9 + 1.34596 * pow(T9,2.5) - 0.209351*T9*T9*T9);
	r = 4.742e4 * ( 1 - 0.85*sqrt(T9) + 0.49*T9 - 0.0962*pow(T9, 1.5) + 8.47e-3 *T9*T9 - 2.8e-4 * pow(T9, 2.5) );
	return r;
	}


// ----------  === PP chain ==  -------------

// 1. H+H -> e+ + nu + H2 (Angulo 1999)
double rate1 (double T9) {
	double r;
	r = 4.08e-15*pow(T9,-0.6667)*exp(-3.381*pow(T9,-0.333))*(1+3.82*T9 + 1.51*T9*T9 + 0.144*pow(T9,3) - 1.14e-2*pow(T9,4));
	return r;
	}

// 2. H2+H2 -> n + He3 (Angulo 1999)
double rate2 (double T9) {
	double r;
	r = 4.67e8*pow(T9,-0.66667)*exp(-4.259*pow(T9, -0.333)) * (1 + 1.079*T9 - 0.1124*T9*T9 + 5.68e-3 * T9*T9*T9);
	return r;
	} 

// 3. H2 + H -> gamma + He3 (Angulo 1999)
double rate3(double T9) {
	double r;
	if (T9 <= 0.11){
	r = 1.81*pow(10, 3)*pow(T9, -0.66667)*exp(-3.721*pow(T9, -0.333))*(1+14.3*T9 - 90.5*pow(T9, 2) + 395 *pow(T9, 3)); }
	else if (T9 > 0.11) {
	r = 2.58*pow(10, 3)*pow(T9, -0.66667)*exp(-3.721*pow(T9, -0.333))*(1+3.96*T9 + 0.116*pow(T9, 2)); }
	return r;
	}

double rev3(double T9){
	double r;
	r = rate9(T9)*1.63e10*pow(T9,1.5)*exp(-63.749/T9);
	return r;
	}

// dd. He3+He3 -> He4 + 2H  (PPI) (Angulo 1999)
double ratedd(double T9){
	return 5.59e10*pow(T9,-0.66667)*exp(-12.277*pow(T9, -0.333))*(1 - 0.135*T9 + 2.54e-2*T9*T9 - 1.29e-3*T9*T9*T9);
	} 

double revdd(double T9){
	return ratedd(T9)*3.392e-10*pow(T9,-1.5)*exp(-149.23/T9);
	}


// dHe4. He3 + He4 -> Be7 + gamma (Angulo 1999)   (PPII)
double ratedHe4(double T9){
	return 5.46e6*pow(T9,-1.5)*exp(-12.827*pow(T9,-0.333))*(1 - 0.307*T9 + 8.81e-2*T9*T9 - 1.06e-2*T9*T9*T9 + 4.46e-4*T9*T9*T9*T9);
	}
double revdHe4(double T9){
	return ratedHe4(T9)*1.113e10*pow(T9, 1.5)*exp(-18.412/T9);
	}

// eBe7. e- + Be7 -> Li7 + nu (Fowler, Caughan, Zimmerman 1975)
double rateeBe7(double T9){
	return 1.34e-10/sqrt(T9)*( 1 - 0.537*pow(T9,0.3333) + 3.86*pow(T9, 0.666667) + 1.20*T9 +0.0027/T9 * exp(2.515e-3/T9) );
	}

// pLi7. H + Li7 -> 2 He4 (Angulo 1999)
double ratepLi7(double T9){
	return 7.20e8*pow(T9, -0.666667)*exp(-8.473*pow(T9,-0.3333) - (T9/6.5)*(T9/6.5))*(1 + 1.05*T9 - 0.653*T9*T9 + 0.185*T9*T9*T9 - 2.12e-2*T9*T9*T9*T9 + 9.30e-4*T9*T9*T9*T9*T9) + 9.85e6*pow(T9,0.576)*exp(-10.415/T9);
	 } 

double revpLi7(double T9){
	return ratepLi7(T9)*4.676*exp(-201.30/T9);
	}

// pBe7. H + Be7 -> B8 + gamma
double ratepBe7(double T9){
	return 2.61e5*pow(T9, -0.666667)*exp(-10.264*pow(T9,-0.333))*(1 - 5.11e-2*T9 + 4.68e-2*T9*T9 - 6.60e-3*T9*T9*T9 + 3.12e-4*T9*T9*T9*T9) + 2.05e3*pow(T9,-1.5)*exp(-7.345/T9);
	}

double revpBe7(double T9){
	return ratepBe7(T9)*1.306e10*pow(T9,1.5)*exp(-1.594/T9);
	}

// B8. B8 -> e+ + nu + 2He4 decay of B8. (Wagoner 1967), in s-1
double rateB8(double T9){
	return 0.89;	
	}


// ----------  === CNO CYCLE ==  -------------

// C12 + H -> N13 + gamma (Angulo 1999)
double ratepC12(double T9){
	return 2.00e7*pow(T9, -0.666667)*exp(-13.692*pow(T9, -0.3333) - (T9/0.46)*(T9/0.46)) * (1 + 9.89*T9 - 59.8*T9*T9 + 266*T9*T9*T9) + 1.00e5*pow(T9, -0.66667)*exp(-4.913/T9) + 4.24e5*pow(T9,-0.66667)*exp(-21.62/T9);
	}

double revpC12(double T9){
	return 8.847e9*pow(T9,1.5)*exp(-22.553/T9)*ratepC12(T9);
	}

// N13 -> C13 + e+  + nu (Angulo 1999)

// C13 + H -> N14 + gamma (Angulo 1999)
double ratepC13(double T9){
	double r = 9.57e7*pow(T9, -0.66667)*(1 + 3.56*T9)*exp(-13.720*pow(T9,-0.3333) - T9*T9) + 1.5e6*pow(T9,-0.3333)*exp(-5.930/T9) + 6.83e5*pow(T9, -0.864)*exp(-12.057/T9);
	return r* (1 - 2.070*exp(-37.838/T9));
	}

double revpC13(double T9){
	return 1.190e10*pow(T9,1.5)*exp(-87.619/T9)*ratepC13(T9);
	}

// N14 + H -> O15 + gamma (Angulo 1999)
double ratepN14(double T9){
	return 4.83e7*pow(T9,-0.6667)*exp(-15.231*pow(T9,-0.3333) - (T9/0.8)*(T9/0.8))*(1 - 2.00*T9 + 3.41*T9*T9 -2.43*T9*T9*T9) + 2.36e3*pow(T9, -1.5)*exp(-3.010/T9) + 6.72e3*pow(T9,0.380)*exp(-9.530/T9);
	}

double revpN14(double T9){
	return 2.699e10*pow(T9,1.5)*exp(-84.677/T9)*ratepN14(T9);
	}

// O15 -> N15 + e+ + nu (Angulo 1999)

// N15 + H -> C12 + He4 (Angulo 1999)
double ratepN15C12(double T9){
	double r;
	if (T9 <= 2.5){ r = 1.12e12*pow(T9,-0.66667)*exp(-15.253*pow(T9,-0.3333) - (T9/0.28)*(T9/0.28))*(1 + 4.95*T9 + 143*T9*T9) + 1.01e8*pow(T9, - 1.5)*exp(-3.643/T9) + 1.19e9*pow(T9,-1.5)*exp(-7.406/T9);}
	else if (T9 > 2.5){ r = 4.17e7*pow(T9,0.917)*exp(-3.292/T9);}
	return r;
	}

double revpN15C12(double T9){
	return 7.060e-1*exp(-57.622/T9)*ratepN15C12(T9);
	}


// N15 + H -> O16 + gamma (Angulo 1999)
double ratepN15O16(double T9){
	double r;
	if (T9 <= 3.5){ r = 1.08e9*pow(T9, -0.66667)*exp(-15.254*pow(T9,-0.3333) - (T9/0.34)*(T9/0.34)) * (1 + 6.15*T9 + 16.4 * T9*T9) + 9.23e3*pow(T9, -1.5)*exp(-3.597/T9) + 3.27e6*pow(T9, -1.5)*exp(-11.024/T9);}
	else if (T9 > 3.5){ r = 3.54e4*pow(T9,0.095)*exp(-2.306/T9);}
	return r;
	}

double revpN15016(double T9){
	return 3.622e10*pow(T9,1.5)*exp(-140.73/T9)*ratepN15O16(T9);
	}

// O16 + H -> F17 + gamma (Angulo 1999)
double ratepO16(double T9){
	double r = 7.37e7*exp(-16.696*pow(T9,-0.3333))*pow(T9,-0.82);
	return r*(1 + 202*exp(-70.348/T9 - 0.161/T9));
	}

double revpO16(double T9){
	return 3.037e9*pow(T9,1.5)*exp(-6.966/T9)*ratepO16(T9);
	}


// F17 -> O17 + e+ + nu (Angulo 1999)
// O17 + H -> N14 + He4 (Angulo 1999)
double ratepO17(double T9){
	double r;
	if (T9 <= 6) {r = 9.20e8*pow(T9,-0.66667)*exp(-16.715*pow(T9,-0.3333) - (T9/0.06)*(T9/0.06))*(1 - 80.31*T9 + 2211*T9*T9) + 9.13e-4*pow(T9,-1.5)*exp(0.7667/T9) + 9.68*pow(T9,-1.5)*exp(-2.083/T9) + 8.13e6*pow(T9,-1.5)*exp(-5.685/T9) + 1.85e6*pow(T9,1.591)*exp(-4.848/T9);}
	else if (T9 > 6) {r = 8.73e6*pow(T9, 0.950)*exp(-7.508/T9);}
	return r*(1 + 1.033*exp(-10.034/T9 - 0.165*T9));
	}

double revpO17(double T9){
	return 6.759e-1*exp(-13.829/T9)*ratepO17(T9);
	}



// 5. C12 + alpha -> O16 + He4
double rate5(double T9){
	double r1 = 6.67e7*pow(T9, -2)*exp(-32.123*pow(T9, -0.333) - (T9/4.6)*(T9/4.6) ) * (1+2.54*T9 + 1.04*T9*T9 - 0.226*T9*T9*T9) + 1.39e3*pow(T9,-1.5)*exp(-28.930/T9);
	double r2 = 6.56e7*pow(T9, -2) * exp(-32.123*pow(T9, -0.333) - (T9/1.3)*(T9/1.3) ) * (1+9.23*T9 - 13.7*T9*T9 + 7.4*T9*T9*T9);
	double r3 = 12.2*T9*T9*exp(-26.9/T9);
	double r = r1 + r2 + r3;
	return r;
	}

// aa. He4 + He4 
double rateaa(double T9){
	double r = 2.43e9*pow(T9, -0.6667)*exp(-13.490*pow(T9, -0.3333) - (T9/0.15)*(T9/0.15) ) * (1+74.5*T9) + 6.09e5*pow(T9, -1.5)*exp(-1.054/T9);
	return r;
	}

// aBe8 He4 + Be8
double rateaBe8(double T9){
	double r;
	r = 2.76e7*pow(T9,-0.6667)*exp(-23.57*pow(T9, -0.333)  - (T9/0.4)*(T9/0.4)) * (1+5.47*T9+326*T9*T9) + 130.7*pow(T9,-1.5)*exp(-3.338/T9)+2.51e4*pow(T9,-1.5)*exp(-20.307/T9);
	return r;
	}



//-----  ==== ALPHA CHAIN ====  -----
	
// 6. He4 + He4 + He4 -> gamma + C12 (Angulo 1999)
double rate6(double T9){
	double r;
	r = rateaa(T9) * rateaBe8(T9);
	if (T9 <= 0.03)	{ r = r * 3.07e-16*(1-29.1*T9 + 1308*T9*T9);}
	else if (T9 > 0.03) {r = r * 3.44e-16 * (1+0.0158*pow(T9,-0.65));}
	return r;
	}

double rev6(double T9){
	return rate6(T9)*2.003e20*T9*T9*T9*exp(-84.415/T9);
	}

//6.bis He4 + He4 + He4 -> gamma + C12 (Fowler, Caughlan & Zimmerman 1975)
double rate6bis(double T9) {
	double a = 0.5; // between 0 and 1
	return 2.49e-8/(T9*T9*T9)*exp(-4.4109/T9) + a*1.35e-7/pow(T9, 1.5)*exp(-24.811/T9);
	} 

double rev6bis(double T9){
	return 2.00e20*T9*T9*T9*exp(-84.424/T9);
	}

// 7. C12 + He4 -> gamma + O16 (Angulo 1999)
// 0.06 < T9 < 10
double rate7(double T9){
	double r, rE1, rE2, rRes;
	rE1 = 6.66e7*pow(T9, -2)*exp(-32.123*pow(T9, -0.333) - ((T9/4.6)*(T9/4.6)) ) * (1 + 2.54*T9 + 1.04*T9*T9 - 0.226*T9*T9*T9) + 1.39 * 1000*pow(T9, -1.5) * exp(-28.930/T9);
	rE2 = 6.56e7*pow(T9, -2)*exp(-32.123*pow(T9, -0.333) - ((T9/1.3)*(T9/1.3)) ) * (1 + 9.23*T9 - 13.7*T9*T9 + 7.4*T9*T9*T9);
	rRes = 19.2*T9*T9*exp(-26.9/T9);
	r = rE1 + rE2 + rRes;
	return r;
	}

double rev7(double T9){
	return rate7(T9)*5.132e10*pow(T9,1.5)*exp(-83.109/T9);
	}

// 8. O16 + He4 -> gamma + Ne20 (Angulo 1999)
// 0.1 < T9 < 10
double rate8(double T9){
	double r = 2.68e10 * pow(T9, -0.666)*exp(-39.76*pow(T9, -0.333) - ((T9/1.6)*(T9/1.6))) + 51.1*pow(T9, -1.5)*exp(-10.32/T9) + 616.1*pow(T9, -1.5)*exp(-12.2/T9) + 0.41*pow(T9, 2.966)*exp(-11.9/T9);
	return r;
	}

double rev8(double T9){
	return rate8(T9)*5.653e10*pow(T9,1.5)*exp(-54.886/T9);
	}

// 9.Ne20 + He4 -> gamma + Mg24 (Angulo 1999)
double rate9(double T9){
	double r;
	if(T9<=1){ r = 8.72*pow(T9, -0.532)*exp(-8.995/T9);}
	else if(T9>1) { r= 3.74e2*pow(T9, 2.229)*exp(-12.681/T9);}
	r = r*(1-7.787*exp(-19.821/T9 - 0.114*T9));
	return r;
	}

double rev9(double T9){
	return rate9(T9)*6.010e10*pow(T9,1.5)*exp(-108.11/T9);
	}

// 10. Mg24 + He4 -> gamma + Si28 (Harris, Fowler, Caughlan, Zimmerman 1983)
double rate10(double T9){
	double rc121 = .1; // factor [0;1], .1 is value in cococubed approx19 (WHY???)
	double r;
	double GT9 = 1. + 5.*exp(-15.882/T9);
	r = 4.78e1/pow(T9, 1.5)*exp(-13.506/T9) + 2.38e3/pow(T9, 1.5)*exp(-15.218/T9) + 2.47e2*pow(T9,1.5)*exp(-15.147/T9) + rc121*1.72e-9/pow(T9,1.5)*exp(-5.028/T9) + 1.25e-3/pow(T9,1.5)*exp(-7.929/T9) + 2.43e1/T9*exp(-11.523/T9);
	r = r/GT9;
	return r;
	}

double rev10(double T9){
	return rate10(T9)*6.27e10*pow(T9,1.5)*exp(-115.881/T9);
	}

// 11. Si28 + He4 -> gamma + S32 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rate11(double T9){
	double tau, A, B, C, D;
	if (T9 < 1) {
	tau = 59.49;
	A = 49.60;
	B = 1.27e-2;
	C = 4.133e-3;
	D = 2.791e-3;
	}
	else if (T9 >= 1) {
	tau = 61.02;
	A = 53.60;
	B = 6.34e-2;
	C = 2.541e-3;
	D = -2.9e-4;
	}

	double r = pow(T9, -0.666)*exp(A- (tau/pow(T9, .333))*(1+B*T9 + C*T9*T9 + D*T9*T9*T9));
	return r;}

double rev11(double T9){
	double GSi28, GS32;
	double ASi28 = -2.002e1;
	double BSi28 = 1.324;
	double CSi28 = 3.104e-2;
	double DSi28 = 0;
	double AS32 = -2.319e1;
	double BS32 = 4.331e-1;
	double CS32 = 1.321e-1;
	double DS32 = 0;
	GSi28 = exp((ASi28 + BSi28*T9 + CSi28*T9*T9 + DSi28*T9*T9*T9)/T9)+1;
	GS32 = exp((AS32 + BS32*T9 + CS32*T9*T9 + DS32*T9*T9*T9)/T9)+1;
	double Q = 6.949;
	double rev = 6.461e10;
	double r = rev*(GSi28/GS32)*pow(T9,1.5)*rate11(T9)*exp(-11.605*Q/T9);
	return r;
	}
	

// 11bis. Si28 + He4 -> gamma + S32 cococubed approx19
double rate11bis(double T9){
	double z = T9;
        double z2 = z*z;
        double z3 = z2*z;
        double aa = 1. + 6.340e-2*z + 2.541e-3*z2 - 2.900e-4*z3;
        double r = 4.82e22 * pow(T9, -0.666) * exp(-61.015 *pow(T9, -0.333) * aa);
	return r;}



// 12. S32 + He4 -> gamma + Ar36 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rate12(double T9){
	double tau, A, B, C, D;
	if (T9 < 1) {
	tau = 65.37;
	A = 52.09;
	B = 1.821e-2;
	C = -5.033e-3;
	D = 5.584e-3;
	}
	else if (T9 >= 1) {
	tau = 66.69;
	A = 55.41;
	B = 4.913e-2;
	C = 4.637e-3;
	D = -4.067e-4;
	}

	double r = pow(T9, -0.666)*exp(A- (tau/pow(T9, .333))*(1+B*T9 + C*T9*T9 + D*T9*T9*T9));
	return r;}

double rev12(double T9){
	double GAr36,GS32;
	double AAr36 = -1.894e1;
	double BAr36 = -5.814e-3;
	double CAr36 = 1.697e-1;
	double DAr36 = 0;
	double AS32 = -2.319e1;
	double BS32 = 4.331e-1;
	double CS32 = 1.321e-1;
	double DS32 = 0;
	double rev = 6.616e10;
	double Q = 6.642;
	GAr36 = exp((AAr36 + BAr36*T9 + CAr36*T9*T9 + DAr36*T9*T9*T9)/T9)+1;
	GS32 = exp((AS32 + BS32*T9 + CS32*T9*T9 + DS32*T9*T9*T9)/T9)+1;
	double r = rev*(GS32/GAr36)*pow(T9,1.5)*rate12(T9)*exp(-11.605*Q/T9);
	return r;
	}

// 13. Ar36 + He4 -> gamma + Ca40 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rate13(double T9){
	double tau,A,B,C,D;
	if (T9 < 1){
	tau = 71.01;
	A = 54.48;
	B = 2.676e-2;
	C = -3.300e-2;
	D = 3.361e-2;
	}
	else if (T9 >= 1){
	tau = 78.27;
	A = 70.11;
	B = 1.458e-1;
	C = -1.069e-2;
	D = 3.790e-4;
	}
	double r = pow(T9, -0.666)*exp(A- (tau/pow(T9, .333))*(1+B*T9 + C*T9*T9 + D*T9*T9*T9));
	return r;
	}

double rev13(double T9){
	double rev = 6.74e10;
	double Q = 7.041;	
	double GAr36, GCa40;
	double AAr36 = -1.894e1;
	double BAr36 = -5.814e-3;
	double CAr36 = 1.697e-1;
	double DAr36 = 0;
	double ACa40 = -4.150e1;
	double BCa40 = 1.636;
	double CCa40 = 1.483e-1;
	double DCa40 = 0;
	GAr36 = exp((AAr36 + BAr36*T9 + CAr36*T9*T9 + DAr36*T9*T9*T9)/T9)+1;
	GCa40 = exp((ACa40 + BCa40*T9 + CCa40*T9*T9 + DCa40*T9*T9*T9)/T9)+1;
	double r = rev*(GAr36/GCa40)*pow(T9,1.5)*rate13(T9)*exp(-11.605*Q/T9);
	return r;
	}


// 14. Ca40 + He4 -> gamma + Ti44 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rate14(double T9){
	double tau,A,B,C,D;
	tau = 76.44;
	A = 56.80;
	B = 1.650e-2;
	C = 5.973e-3;
	D = -3.889e-4;
	double r = pow(T9, -0.666)*exp(A- (tau/pow(T9, .333))*(1+B*T9 + C*T9*T9 + D*T9*T9*T9));
	return r;
	}

double rev14(double T9){
	double rev = 6.843e10;
	double Q = 5.128;
	double GCa40, GTi44;
	double ACa40 = -4.150e1;
	double BCa40 = 1.636;
	double CCa40 = 1.483e-1;
	double DCa40 = 0;
	double ATi44 = -1.111e1;
	double BTi44 = 6.293e-1;
	double CTi44 = 1.732e-1;
	double DTi44 = 0;
	GCa40 = exp((ACa40 + BCa40*T9 + CCa40*T9*T9 + DCa40*T9*T9*T9)/T9)+1;
	GTi44 = exp((ATi44 + BTi44*T9 + CTi44*T9*T9 + DTi44*T9*T9*T9)/T9)+1;
	double r = rev*(GCa40/GTi44)*pow(T9,1.5)*rate14(T9)*exp(-11.605*Q/T9);
	return r;
	}
	
// 15. Ti44 + He4 -> gamma + Cr48 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rate15(double T9){
	double tau, A, B, C, D;
	if (T9<1){
	tau = 81.66;	
	A = 56.98;
	B = -8.364e-2;
	C = 2.085e-1;
	D = -7.477e-2;
	}
	else if (T9>1){
	tau = 81.23;
	A = 60.18;
	B = 1.066e-1;
	C = -1.102e-2;
	D = 5.324e-4;
	}
	double r = pow(T9, -0.666)*exp(A- (tau/pow(T9, .333))*(1+B*T9 + C*T9*T9 + D*T9*T9*T9));
	return r;
	}

double rev15(double T9){
	double GTi44, GCr48;
	double rev = 6.928e10;
	double Q = 7.694;	
	double ATi44 = -1.111e1;
	double BTi44 = 6.293e-1;
	double CTi44 = 1.732e-1;
	double DTi44 = 0;
	double ACr48 = -6.100;
	double BCr48 = -3.738e-2;
	double CCr48 = 2.125e-1;
	double DCr48 = 0;
	GTi44 = exp((ATi44 + BTi44*T9 + CTi44*T9*T9 + DTi44*T9*T9*T9)/T9)+1;
	GCr48 = exp((ACr48 + BCr48*T9 + CCr48*T9*T9 + DCr48*T9*T9*T9)/T9)+1;	
	double r = rev*(GTi44/GCr48)*pow(T9,1.5)*rate15(T9)*exp(-11.605*Q/T9);
	return r;
	}

// 16. Cr48 + He4 -> gamma + Fe52 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rate16(double T9){
	double tau, A, B, C, D;
	tau = 81.42;
	A = 53.00;
	B = 6.325e-2;
	C = -5.671e-3;
	D = 2.848e-4;
	double r = pow(T9, -0.666)*exp(A- (tau/pow(T9, .333))*(1+B*T9 + C*T9*T9 + D*T9*T9*T9));
	return r;
	}

double rev16(double T9){
	double GCr48, GFe52;
	double rev = 7.001e10;
	double Q = 7.943;
	double ACr48 = -6.100;
	double BCr48 = -3.738e-2;
	double CCr48 = 2.125e-1;
	double DCr48 = 0;
	double AFe52 = -9.298;
	double BFe52 = 1.309;
	double CFe52 = 1.125e-1;
	double DFe52 = 0;
	GCr48 = exp((ACr48 + BCr48*T9 + CCr48*T9*T9 + DCr48*T9*T9*T9)/T9)+1;	
	GFe52 = exp((AFe52 + BFe52*T9 + CFe52*T9*T9 + DFe52*T9*T9*T9)/T9)+1;	
	double r = rev*(GCr48/GFe52)*pow(T9,1.5)*rate16(T9)*exp(-11.605*Q/T9);
	return r;	
	}	

// 17. Fe52 + He4 -> gamma + Ni56 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rate17(double T9){
	double tau, A, B, C, D;

	tau = 91.67;
	A = 62.22;
	B = 7.846e-2;
	C = -7.430e-3;
	D = 3.723e-4;
	double r = pow(T9, -0.666)*exp(A- (tau/pow(T9, .333))*(1+B*T9 + C*T9*T9 + D*T9*T9*T9));
	return r;
	}

double rev17(double T9){
	double GFe52, GNi56;
	double rev = 7.064e10;
	double Q = 8.001;
	double AFe52 = -9.298;
	double BFe52 = 1.309;
	double CFe52 = 1.125e-1;
	double DFe52 = 0;
	double ANi56 = -2.539e1;
	double BNi56 = -1.067;
	double CNi56 = 3.467e-1;
	double DNi56 = 0;
	GFe52 = exp((AFe52 + BFe52*T9 + CFe52*T9*T9 + DFe52*T9*T9*T9)/T9)+1;	
	GNi56 = exp((ANi56 + BNi56*T9 + CNi56*T9*T9 + DNi56*T9*T9*T9)/T9)+1;	
	double r = rev*(GFe52/GNi56)*pow(T9,1.5)*rate17(T9)*exp(-11.605*Q/T9);
	return r;	
	}

// Heavier reactions

// C12C12. C12 + C12 -> Mg24 + gamma (Fowler, Caughlan & Zimmerman 1975)
/*
double rateC12C12(double T9){
	double r;
	if (T9 < 1) { r = 0;}
	else if (T9 > 1){
	double T9A = T9/(1. + 0.067*T9);
	r = 1.26e27*pow(T9A, 0.8333)/pow(T9,1.5)*exp(-84.165/pow(T9A,0.3333))/exp(-0.01*T9A*T9A*T9A*T9A) + 5.56e-3*exp(1.685*pow(T9A,0.66666));
	}
	return r;
	}
*/

// C12C12. C12 + C12 -> Mg24 + gamma (Caughlan & Fowler 1988)
double rateC12C12(double T9){	
	double T9A = T9/(1+0.0396*T9);
	return 4.27e26*pow(T9A, 0.83333)/pow(T9, 1.5)*exp(-84.165/pow(T9A, 0.33333) - 2.12e-3*T9*T9*T9);
	}

// C12C12a. C12 + C12 -> He4 + Ne20 (Caughlan & Fowler 1988)
double rateC12C12a(double T9){
	return rateC12C12(T9)*0.65;
	}

// C12C12p. C12 + C12 -> H + Na23 (Caughlan & Fowler1988)
double rateC12C12p(double T9){
	return rateC12C12(T9)*0.3;
	}

// C12O16. C12 + O16 -> Si28 + gamma (Caughlan & Fowler 1988) 
double rateC12O16(double T9){
	double r;
	if (T9 < 4) { r = 0;}
	else if (T9 > 4) {
	double T9A = T9/(1. + 0.067*T9);
	r = 1.72e31*pow(T9A, 0.8333)/pow(T9, 1.5)*exp(-106.594/pow(T9A,0.3333))/exp(-0.180*T9A*T9A) + 1.06e-3*exp(2.562*pow(T9A,0.66666));
	}
	return r; 
	}

// C12O16a. C12 + O16 -> He4 + Mg24 (Caughlan & Fowler 1988) 
double rateC12O16a(double T9){
	double r = rateC12O16(T9);
//	if (T9 < 2.6) {r = r*0.37;}
//	else if (T9 > 2.6) {r = r*0.38;}
	r = r*0.4;
	return r;
	}

// C12O16p. C12 + O16 -> p + Al27 (Caughlan & Fowler 1988) 
double rateC12O16p(double T9){
	double r = rateC12O16(T9);
//	if (T9 < 2.6) {r = r*0.54;}
//	else if (T9 > 2.6) {r = r*0.68;}
	r = r*0.5;
	return r;
	}

// C12O16n. C12 + O16 -> n + Si27 (Caughlan & Fowler 1988) 
double rateC12O16n(double T9){
	double r = rateC12O16(T9);
//	if (T9 < 1.75) {r = r*0.09;}
//	else if (T9 > 2.6) {r = r*0.14;}
	r = r*0.1;
	return r;
	}
		
//O16O16. O16 + O16 -> S32 + gamma (Caughlan & Fowler1988)
double rateO16O16(double T9){
	double r;
	if (T9 < 6) { r = 0;}
	else if (T9 > 6){
	double T9A = T9/(1. + 0.067*T9);
	r = 3.61e37*pow(T9A, 0.8333)/pow(T9,1.5)*exp(-135.93/pow(T9A,0.3333))/exp(-0.032*T9A*T9A*T9A*T9A) + 3.89e-4*exp(2.659*pow(T9A,0.66666));
	}
	return r;
	}

// O16O16a. O16 + O16 -> He4 + Si28 (Caughlan & Fowler1988)
double rateO16O16a(double T9){
	double r = rateO16O16(T9);
//	if (T9 < 2.6) {r = r*0.30;}
//	else if (T9 > 2.6) {r = r*0.32;}
	r = r*0.29;
	return r;
	}

// O16O16p. O16 + O16 -> p + P31 (Caughlan & Fowler1988)
double rateO16O16p(double T9){
	double r = rateO16O16(T9);
//	if (T9 < 2.6) {r = r*0.85;}
//	else if (T9 > 2.6) {r = r*1.09;}
	r = r*0.8;
	return r;
	}

// O16O16n. O16 + O16 -> n + S31 (Caughlan & Fowler1988)
double rateO16O16n(double T9){
	double r = rateO16O16(T9);
//	if (T9 < 2.6) {r = r*0.25;}
//	else if (T9 > 2.6) {r = r*0.16;}
	r = r*0.14;
	return r;
	}




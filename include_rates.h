// include all

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std;

double RG (double (double, double, double) , double , double,  double, double, double);

double rate0 (double);		// 0. H+n -> h2 + gamma
// pp chain
double rate1 (double);		// 1. H+H -> e+ + nu + H2 (Angulo 1999)
double rate2 (double);		// 2. H2+H2 -> n + He3 (Angulo 1999)
double rate3 (double);		// 3. H2 + H -> gamma + He3 (Angulo 1999)
double rev3(double);
double ratedd(double);		// dd. He3+He3 -> He4 + 2H  (PPI) (Angulo 1999)
double revdd(double);
double ratedHe4(double);	// dHe4. He3 + He4 -> Be7 + gamma (Angulo 1999)   (PPII)
double revdHe4(double);
double rateeBe7(double);	// eBe7. e- + Be7 -> Li7 + nu (Fowler, Caughan, Zimmerman 1975)
double ratepLi7(double);	// pLi7. H + Li7 -> 2 He4 (Angulo 1999)
double revpLi7(double);
double ratepBe7(double);	// pBe7. H + Be7 -> B8 + gamma
double revpBe7(double);

double rate5 (double);		// 5. C12 + alpha -> O16 + He4
double rateaa(double);		// aa. He4 + He4
double rateaBe8(double);	// aBe8 He4 + Be8

// CNO cycle
double ratepC12(double);	// C12 + H -> N13 + gamma (Angulo 1999)
double revpC12(double);
double ratepC13(double);	// C13 + H -> N14 + gamma (Angulo 1999)
double revpC13(double);
double ratepN14(double);	// N14 + H -> O15 + gamma (Angulo 1999)
double revpN14(double);	
double ratepN15C12(double);	// N15 + H -> C12 + He4 (Angulo 1999)
double revpN15C12(double);
double ratepN15O16(double);	// N15 + H -> O16 + gamma (Angulo 1999)
double revpN15016(double);
double ratepO16(double);	// O16 + H -> F17 + gamma (Angulo 1999)
double revpO16(double);
double ratepO17(double);	// O17 + H -> N14 + He4 (Angulo 1999)
double revpO17(double);

// alpha chain
double rate6(double);		// 6. He4 + He4 + He4 -> gamma + C12 (Angulo 1999)
double rev6(double);
double rate6bis(double);	// 6bis. He4 + He4 + He4 -> gamma + C12 (Fowler, Caughlan & Zimmerman 1975)
double rev6bis(double);
double rate7(double);		// 7. C12 + He4 -> gamma + O16 (Angulo 1999)
double rev7(double);
double rate8(double);		// 8. O16 + He4 -> gamma + Ne20 (Angulo 1999)
double rev8(double);
double rate9(double);		// 9.Ne20 + He4 -> gamma + Mg24 (Angulo 1999)
double rev9(double);
double rate10(double);		// 10. Mg24 + He4 -> gamma + Si28 (Harris, Fowler, Caughlan, Zimmerman 1983)
double rev10(double);
double rate11(double);		// 11. Si28 + He4 -> gamma + S32 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rev11(double);
double rate11bis(double);	// 11bis. Si28 + He4 -> gamma + S32 cococubed approx19
double rate12(double);		// 12. S32 + He4 -> gamma + Ar36 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rev12(double);
double rate13(double);		// 13. Ar36 + He4 -> gamma + Ca40 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rev13(double);
double rate14(double);		// 14. Ca40 + He4 -> gamma + Ti44 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rev14(double);
double rate15(double);		// 15. Ti44 + He4 -> gamma + Cr48 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rev15(double);
double rate16(double);		// 16. Cr48 + He4 -> gamma + Fe52 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rev16(double);
double rate17(double);		// 17. Fe52 + He4 -> gamma + Ni56 (Woosley, Fowler, Holmes, Zimmerman 1978)
double rev17(double);
double rateC12C12(double);	// C12C12. C12 + C12 -> gamma + Mg24 (Fowler, Caughlan & Zimmerman 1975)
double rateC12C12bis(double);
double rateC12C12a(double);	// C12C12a. C12 + C12 -> He4 + Ne20 (Fowler, Caughlan & Zimmerman 1975)
double rateC12C12p(double);	// C12C12p. C12 + C12 -> H + Na23 (Fowler, Caughlan & Zimmerman 1975)
double rateC12O16(double);	// C12O16. C12 + O16 -> gamma + Si28 (Fowler, Caughlan & Zimmerman 1975) 
double rateC12O16a(double);	// C12O16a. C12 + O16 -> a + Mg24 (Fowler, Caughlan & Zimmerman 1975) 
double rateC12O16p(double);	// C12O16p. C12 + O16 -> p + Al27 (Fowler, Caughlan & Zimmerman 1975) 
double rateC12O16n(double);	// C12O16n. C12 + O16 -> n + Si27 (Fowler, Caughlan & Zimmerman 1975) 
double rateO16O16(double);	// O16O16. O16 + O16 -> gamma + S32  (Fowler, Caughlan & Zimmerman 1975)
double rateO16O16a(double);	// O16O16a. O16 + O16 -> a + Si28  (Fowler, Caughlan & Zimmerman 1975)
double rateO16O16p(double);	// O16O16p. O16 + O16 -> p + P31  (Fowler, Caughlan & Zimmerman 1975)
double rateO16O16n(double);	// O16O16n. O16 + O16 -> n + S31  (Fowler, Caughlan & Zimmerman 1975)
	


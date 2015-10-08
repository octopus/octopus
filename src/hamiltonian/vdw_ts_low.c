/*
** Copyright (C) 2015 Xavier Andrade, Carlos Borca, Alfredo Correa
** 
** This library is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This library is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
*/

#include <stdio.h>
#include <math.h>

#ifndef _TEST
#include "config.h"
#endif

/* Function to retrieve Van der Waals parameters of the free atoms. */
void get_vdw_params (const int zatom, 
		     double * c6, double * alpha, double * r0) {
  
  switch(zatom) {
    
    // Hydrogen (H)
    case 1:
      *alpha = 4.500000;
      *c6 = 6.500000;
      *r0 = 3.100000;
      break;
    
    // Helium (He)
    case 2:
      *alpha = 1.380000;
      *c6 = 1.460000;
      *r0 = 2.650000;
      break;
    
    // Lithium (Li)
    case 3:
      *alpha = 164.200000;
      *c6 = 1387.000000;
      *r0 = 4.160000;
      break;
    
    // Berillium (Be)
    case 4:
      *alpha = 38.000000;
      *c6 = 214.000000;
      *r0 = 4.170000;
      break;
    
    // Boron (B)
    case 5:
      *alpha = 21.000000;
      *c6 = 99.500000;
      *r0 = 3.890000;
      break;
    
    // Carbon (C)
    case 6:
      *alpha = 12.000000;
      *c6 = 46.600000;
      *r0 = 3.590000;
      break;
      
    // Nitrogen (N)
    case 7:
      *alpha = 7.400000;
      *c6 = 24.200000;
      *r0 = 3.340000;
      break;
      
    // Oxygen (O)
    case 8:
      *alpha = 5.400000;
      *c6 = 15.600000;
      *r0 = 3.190000;
      break;
      
    // Fluorine (F)
    case 9:
      *alpha = 3.800000;
      *c6 = 9.520000;
      *r0 = 3.040000;
      break;
      
    // Neon (Ne)
    case 10:
      *alpha = 2.670000;
      *c6 = 6.380000;
      *r0 = 2.910000;
      break;
      
    // Sodium (Na)
    case 11:
      *alpha = 162.700000;
      *c6 = 1556.000000;
      *r0 = 3.730000;
      break;
      
    // Magnesium (Mg)
    case 12:
      *alpha = 71.000000;
      *c6 = 627.000000;
      *r0 = 4.270000;
      break;
      
    // Aluminium (Al)
    case 13:
      *alpha = 60.000000;
      *c6 = 528.000000;
      *r0 = 4.330000;
      break;
      
    // Silicon (Si)
    case 14:
      *alpha = 37.000000;
      *c6 = 305.000000;
      *r0 = 4.200000;
      break;
      
    // Phosphorus (P)
    case 15:
      *alpha = 25.000000;
      *c6 = 185.000000;
      *r0 = 4.010000;
      break;
      
    // Sulphur (S)
    case 16:
      *alpha = 19.600000;
      *c6 = 134.000000;
      *r0 = 3.860000;
      break;
      
    // Clorine (Cl)
    case 17:
      *alpha = 15.000000;
      *c6 = 94.600000;
      *r0 = 3.710000;
      break;
      
    // Argon (Ar)
    case 18:
      *alpha = 11.100000;
      *c6 = 64.300000;
      *r0 = 3.550000;
      break;
      
    // Potassium (K)
    case 19:
      *alpha = 292.900000;
      *c6 = 3897.000000;
      *r0 = 3.710000;
      break;
      
    // Calcium (Ca)
    case 20:
      *alpha = 160.000000;
      *c6 = 2221.000000;
      *r0 = 4.650000;
      break;
      
    // Scandium (Sc)
    case 21:
      *alpha = 120.000000;
      *c6 = 1383.000000;
      *r0 = 4.590000;
      break;
      
    // Titanium (Ti)
    case 22:
      *alpha = 98.000000;
      *c6 = 1044.000000;
      *r0 = 4.510000;
      break;
      
    // Vanadium (V)
    case 23:
      *alpha = 84.000000;
      *c6 = 832.000000;
      *r0 = 4.440000;
      break;
      
    // Chromium (Cr)
    case 24:
      *alpha = 78.000000;
      *c6 = 602.000000;
      *r0 = 3.990000;
      break;
      
    // Manganese (Mn)
    case 25:
      *alpha = 63.000000;
      *c6 = 552.000000;
      *r0 = 3.970000;
      break;
      
    // Iron (Fe)
    case 26:
      *alpha = 56.000000;
      *c6 = 482.000000;
      *r0 = 4.230000;
      break;
      
    // Cobalt (Co)
    case 27:
      *alpha = 50.000000;
      *c6 = 408.000000;
      *r0 = 4.180000;
      break;
      
    // Nickel (Ni)
    case 28:
      *alpha = 48.000000;
      *c6 = 373.000000;
      *r0 = 3.820000;
      break;
      
    // Copper (Cu)
    case 29:
      *alpha = 42.000000;
      *c6 = 253.000000;
      *r0 = 3.760000;
      break;
      
    // Zinc (Zn)
    case 30:
      *alpha = 40.000000;
      *c6 = 284.000000;
      *r0 = 4.020000;
      break;
      
    // Gallium (Ga)
    case 31:
      *alpha = 60.000000;
      *c6 = 498.000000;
      *r0 = 4.190000;
      break;
      
    // Germanium (Ge)
    case 32:
      *alpha = 41.000000;
      *c6 = 354.000000;
      *r0 = 4.200000;
      break;
      
    // Arsenic (As)
    case 33:
      *alpha = 29.000000;
      *c6 = 246.000000;
      *r0 = 4.110000;
      break;
      
    // Selenium (Se)
    case 34:
      *alpha = 25.000000;
      *c6 = 210.000000;
      *r0 = 4.040000;
      break;
      
    // Bromine (Br)
    case 35:
      *alpha = 20.000000;
      *c6 = 162.000000;
      *r0 = 3.930000;
      break;
      
    // Krypton (Kr)
    case 36:
      *alpha = 16.800000;
      *c6 = 129.600000;
      *r0 = 3.820000;
      break;
      
    // Rubidium (Rb)
    case 37:
      *alpha = 319.200000;
      *c6 = 4691.000000;
      *r0 = 3.720000;
      break;
      
    // Strontium (Sr)
    case 38:
      *alpha = 199.000000;
      *c6 = 3170.000000;
      *r0 = 4.540000;
      break;
    
    // Elements from 39 - Yttrium (Y) to 44 Ruthenium (Ru) are not included.
    
    // Rhodium (Rh)
    case 45:
      *alpha = 56.1;
      *c6 = 469.0;
      *r0 = 3.95;
      break;
      
    // Palladium (Pd)
    case 46:
      *alpha = 23.680000;
      *c6 = 157.500000;
      *r0 = 3.66000;
      break;
      
    // Silver (Ag)
    case 47:
      *alpha = 50.600000;
      *c6 = 339.000000;
      *r0 = 3.820000;
      break;
      
    // Cadmium (Cd)
    case 48:
      *alpha = 39.7;
      *c6 = 452.0;
      *r0 = 3.99;
      break;
      
    // Elements from 49 - Indium (In) to 51 - Antimony (Sb) are not included.
      
    // Tellurium (Te)
    case 52:
      *alpha = 37.65;
      *c6 = 396.0;
      *r0 = 4.22;
      break;
      
    // Iodine (I)
    case 53:
      *alpha = 35.000000;
      *c6 = 385.000000;
      *r0 = 4.170000;
      break;
      
    // Xenon (Xe)
    case 54:
      *alpha = 27.300000;
      *c6 = 285.900000;
      *r0 = 4.080000;
      break;
      
    // Element 55 - Caesium (Cs) is not included.
    
    // Barium (Ba)
    case 56:
      *alpha = 275.0;
      *c6 = 5727.0;
      *r0 = 4.77;
      break;
      
    // Elements from 57 - Lanthanum (La) to 76 - Osmium (Os) are not included.
      
    // Iridium (Ir)
    case 77:
      *alpha = 42.51;
      *c6 = 359.1;
      *r0 = 4.00;
      break;
      
    // Platinum (Pt)
    case 78:
      *alpha = 39.68;
      *c6 = 347.1;
      *r0 = 3.92;
      break;
      
    // Gold (Au)
    case 79:
      *alpha = 36.5;
      *c6 = 298.0;
      *r0 = 3.86;
      break;
      
    // Mercury (Hg)
    case 80:
      *alpha = 33.9;
      *c6 = 392.0;
      *r0 = 3.98;
      break;
      
    // Element 81 - Thallium (Tl) is not included.
      
    // Lead (Pb)
    case 82:
      *alpha = 61.8;
      *c6 = 697.0;
      *r0 = 4.31;
      break;
      
    // Bismuth (Bi)
    case 83:
      *alpha = 49.02;
      *c6 = 571.0;
      *r0 = 4.32;
      
    // Elements from 84 - Polonium (Po) to 118 - Ununoctium (Uuo) are not included.
    
  }
}


/* Damping function. */
void fdamp (const double rr, const double r0ab, 
	    double * ff, double * dffdrab, double * dffdr0) {

  const double dd = 20.0;
  const double sr = 0.94; // Value for PBE. Should be 0.96 for PBE0.

  // Calculate the damping coefficient.
  double ee = exp(-dd*((rr/(sr*r0ab)) - 1.0));
  *ff = 1.0/(1.0 + ee);
  double dee = ee*(*ff)*(*ff);
  
  // Calculate the derivative of the damping function with respect to the distance between atoms A and B.
  *dffdrab = (dd/(sr*r0ab))*dee;
  
  // Calculate the derivative of the damping function with respect to the distance between the van der Waals radius.
  *dffdr0 = -dd*rr/(sr*r0ab*r0ab)*dee;

}


/* Calculation of the square of separations. */
void distance (const int iatom, const int jatom, const double coordinates[], 
	       double * rr, double * rr2, double * rr6, double *rr7) {
  
  double x_ij = coordinates[3*iatom + 0] - coordinates[3*jatom + 0];
  double y_ij = coordinates[3*iatom + 1] - coordinates[3*jatom + 1];
  double z_ij = coordinates[3*iatom + 2] - coordinates[3*jatom + 2];
  
  *rr2 = x_ij*x_ij + y_ij*y_ij + z_ij*z_ij;
  *rr6 = rr2[0]*rr2[0]*rr2[0]; // This is the same as: *rr6 = (*rr2)*(*rr2)*(*rr2)
  *rr = sqrt(rr2[0]); // This is the same as: *rr = sqrt(*rr2)
  *rr7 = rr6[0]*rr[0]; // This is the same as: *rr6 = (*rr6)*(*rr)
  
  // Print information controls.
  //printf("R_(%i-%i)= %f\n", iatom+1, jatom+1, *rr);
  //printf("R_(%i-%i)^2 = %f\n", iatom+1, jatom+1, *rr2);
  //printf("R_(%i-%i)^6 = %f\n", iatom+1, jatom+1, *rr6);
  
}


/* Function to calculate the Van der Waals energy... and forces */
void vdw_calculate (const int natoms, const int zatom[], const double coordinates[], const double volume_ratio[], 
		    double * energy, double force[], double derivative_coeff[]) {
  
  int ia;

  *energy = 0.0;

  // Loop to calculate the pair-wise Van der Waals energy correction.
  for (ia = 0; ia < natoms; ia++) {

    double c6_a, alpha_a, r0_a;
    int ib;
    
    force[3*ia + 0] = 0.0;
    force[3*ia + 1] = 0.0;
    force[3*ia + 2] = 0.0;

    derivative_coeff[ia] = 0.0;
    
    get_vdw_params(zatom[ia], &c6_a, &alpha_a, &r0_a);
    
    for (ib = 0; ib < natoms; ib++) {

      double c6_b, alpha_b, r0_b;

      if (ia == ib) continue;

      // Pair-wise calculation of separations.
      double rr, rr2, rr6, rr7;
      distance(ia, ib, coordinates, &rr, &rr2, &rr6, &rr7);
      
      get_vdw_params(zatom[ib], &c6_b, &alpha_b, &r0_b);
      
      // Determination of c6abfree, for isolated atoms a and b.
      double num = 2.0*c6_a*c6_b;
      double den = (alpha_b/alpha_a)*c6_a + (alpha_a/alpha_b)*c6_b;

      double c6abfree = num/den;

      // Determination of c6ab, for bonded atoms a and b.
      double c6ab = volume_ratio[ia]*volume_ratio[ib]*c6abfree;
      
      // Determination of the effective radius of atom a.
      double r0ab = cbrt(volume_ratio[ia])*r0_a + cbrt(volume_ratio[ib])*r0_b;

      // Pair-wise calculation of the damping coefficient
      double ff;
      double dffdrab;
      double dffdr0;
      fdamp(rr, r0ab, &ff, &dffdrab, &dffdr0);
      
      // Pair-wise correction to energy.
      *energy += -0.5*ff*c6ab/rr6;
      
      // Calculation of the pair-wise partial energy derivative with respect to the distance between atoms A and B.
      double deabdrab = -dffdrab*c6ab/rr6 + 6.0*ff*c6ab/rr7;
      
      // Derivative of the AB van der Waals separation with respect to the volume ratio of atom A.
      double dr0dvra = r0_a/(3.0*pow(volume_ratio[ia], 2.0/3.0));
      
      // Derivative of the damping function with respecto to the volume ratio of atom A.
      double dffdvra = dffdr0*dr0dvra;
      
      // Calculation of the pair-wise partial energy derivative with respect to the volume ratio of atom A.
      double deabdvra = -dffdvra*c6ab/rr6 - ff*volume_ratio[ib]*c6abfree/rr6;
      
      force[3*ia + 0] += -deabdrab*(coordinates[3*ia + 0] - coordinates[3*ib + 0])/rr;
      force[3*ia + 1] += -deabdrab*(coordinates[3*ia + 1] - coordinates[3*ib + 1])/rr;
      force[3*ia + 2] += -deabdrab*(coordinates[3*ia + 2] - coordinates[3*ib + 2])/rr;
      
      derivative_coeff[ia] += deabdvra;
      
      // Print information controls.
      //printf("Distance between atoms %i and %i = %f.\n", ia+1, ib+1, rr);
      //printf("Atom %i, c6= %f, alpha= %f, r0= %f.\n", ia+1, c6_a, alpha_a, r0_a);
      //printf("Atom %i, c6= %f, alpha= %f, r0= %f.\n", ib+1, c6_b, alpha_b, r0_b);
      //printf("For atoms %i and %i, c6ab= %f.\n", ia+1, ib+1, c6abfree);
    }
  }
  //printf("THE FINAL VAN DER WAALS ENERGY CORRECTION IS = %f.\n", *energy);
  //printf("THE FINAL VAN DER WAALS FORCE CORRECTION IS = %f.\n", *force);
}

#ifndef _TEST
/* This is a wrapper to be called from Fortran. */
void FC_FUNC_(f90_vdw_calculate, F90_VDW_CALCULATE) (const int * natoms, const int zatom[], const double coordinates[], const double volume_ratio[],
						     double * energy, double force[], double derivative_coeff[]) {
  vdw_calculate(*natoms, zatom, coordinates, volume_ratio, energy, force, derivative_coeff);
}
#endif

/* Main test function. */
#ifdef _TEST
int main () {
  const int natoms = 3;
  const int zatom[] = {23, 29, 31};
  const double volume_ratio[] = {1.0, 1.0, 1.0};
  double energy;
  double force[natoms*3];
  double derivative_coeff[natoms];
  
  double x;
  for(x = 0.1; x < 10; x += 0.1) {
    double coordinates[] = {0.2, -0.3, 0.5,  -0.7, 1.1, -1.3,  x, 0.0, 0.0};
    
    vdw_calculate(natoms, zatom, coordinates, volume_ratio, &energy, force, derivative_coeff);
        
    coordinates[5] += 0.001;
    double energy_2;
    
    vdw_calculate(natoms, zatom, coordinates, volume_ratio, &energy_2, force, derivative_coeff);
    printf("%f %f %f %f\n", x, energy, force[5], -(energy_2-energy)/0.001);
  }
}
#endif

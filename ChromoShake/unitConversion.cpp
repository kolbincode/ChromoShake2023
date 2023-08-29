#include "unitConversion.h"
#include <math.h>

//---------------------------------------------------------------------
// Helper functions and constant

#ifndef	M_PI
// define Pi as 2*arccosine(0) ; arccosine of 0 is Pi/2 radians
const double M_PI = 2*acos(0.0);
#endif

static inline	double	sqr(double a) { return (a)*(a); }

//---------------------------------------------------------------------
// Functions declared in the header

double dna_tensile_spring_constant( double segment_length_meters,       //default:    10e-9 meters; mass separation = 10nm
                                    double dna_youngs_modulus_pa,       //default:    2e9 Pascals; modulus_gigaPascal = 2; set constant
                                    double dna_radius_meters)           //default:    6e-10 meters; dna_radius_nanometers = 0.6
{
    return dna_youngs_modulus_pa
           * (M_PI * sqr(dna_radius_meters))
           / segment_length_meters;
}

double dna_bending_rigidity(    double dna_youngs_modulus_pa,           //will vary with persistence length; 0.2 GPa is 5nm, 2.0 GPa is 50nm, 20 Gpa is 500nm
                                double dna_radius_meters)               // see above
{
    return ( dna_youngs_modulus_pa * M_PI * pow(dna_radius_meters, 4.0) )
                / 4.0;
}

double dna_bending_spring_constant( double segment_length_meters,                   //
                                    double dna_youngs_modulus_pa,                   //
                                    double dna_radius_meters)                       //
{
    return 2 * dna_bending_rigidity( dna_youngs_modulus_pa, dna_radius_meters )     //
        / sqr(segment_length_meters);
}

double Brownian_force_equivalent(   double effective_radius_meters,         //default:    rest_length 10 nanometers (??)
                                    double temperature_kelvin,              //default:    temperature_Celsius = 25
                                    double viscosity_pascal_seconds)        /*default:    viscosity_centiPoise = 1
                                                                                        (10 Poise = 1cPoise = 1 Pa * sec) */
{
    double kBT = 6 * M_PI * viscosity_pascal_seconds * effective_radius_meters * temperature_kelvin * Boltzmann_constant;
    return sqrt(2*kBT);
}

extern double mass_damping_equivalent(  double viscosity_pascal_seconds,
                                        double sphere_radius_meters)
{
    double drag_coefficient = 6 * M_PI * viscosity_pascal_seconds * sphere_radius_meters;    // where is sphere_radius_meters?
    return drag_coefficient;
}


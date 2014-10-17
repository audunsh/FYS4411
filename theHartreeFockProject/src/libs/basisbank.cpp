#include <basisbank.h>
#include <armadillo>
#include <basis.h>
#include <primitive.h>
using namespace std;
using namespace arma;
basisbank::basisbank(basis BS){
    bs = BS;}

//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:51 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:he STO-6G
void basisbank::add_STO6G_He(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,65.984568,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,12.098198,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,3.384640,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,1.162715,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,0.451516,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,0.185959,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
}


//# STO-6G EMSL Basis Set Exchange Library 10/16/14 3:08 PM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:h STO-6G
void basisbank::add_STO6G_H(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,35.523221,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,6.513144,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,1.822143,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,0.625955,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,0.243077,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,0.100112,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:55 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:f STO-6G
void basisbank::add_STO6G_F(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,1728.626574,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,316.941790,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,88.668891,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,30.460157,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,11.828570,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,4.871659,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,67.032281,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,13.267438,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,4.123510,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,1.586463,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.689002,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.315820,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,67.032281,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,67.032281,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,67.032281,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,13.267438,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,13.267438,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,13.267438,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,4.123510,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,4.123510,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,4.123510,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,1.586463,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,1.586463,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,1.586463,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.689002,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.689002,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.689002,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.315820,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.315820,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.315820,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


//# 6-311G EMSL Basis Set Exchange Library 10/17/14 4:49 AM
//# Elements References
//# -------- ----------
//# H, Li - Ne: R. Krishnan, J.S. Binkley, R. Seeger and J.A. Pople,
//# J. Chem. Phys. 72, 650 (1980)
//# Na - Ar: A.D. McLean and G.S. Chandler J. Chem. Phys. 72, 5639, (1980).
//# K - Ca: J-P. Blaudeau, M. P. McGrath, L.A. Curtiss and L. Radom,
//# J. Chem. Phys. 107, 5016 (1997).
//# Ga - Kr: L. A. Curtiss, M. P. McGrath, J-P. Blandeau, N. E. Davis,
//# R. C. Binning, Jr. L. Radom, J. Chem. Phys. 103, 6104 (1995).
//# I : M.N. Glukhovstev, A. pross, M.P. McGrath, L. Radom, J. Chem. Phys.
//# 103, 1878 (1995)
//#
//Ignored the following information from file:h 6-311G
void basisbank::add_6311G_H(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.025494,33.865000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.190373,5.094790,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.852161,1.158790,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(1.000000,0.325840,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    bs.add_state();
    Primitive S0A2 = bs.turbomolePrimitive(1.000000,0.102741,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:53 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:b STO-6G
void basisbank::add_STO6G_B(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,506.011147,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,92.776814,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,25.955658,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,8.916445,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,3.462507,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,1.426054,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,23.194575,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,4.590810,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,1.426819,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,0.548948,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.238410,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.109280,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,23.194575,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,23.194575,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,23.194575,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,4.590810,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,4.590810,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,4.590810,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,1.426819,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,1.426819,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,1.426819,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,0.548948,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,0.548948,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,0.548948,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.238410,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.238410,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.238410,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.109280,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.109280,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.109280,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:57 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:mg STO-6G
void basisbank::add_STO6G_Mg(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,2600.756771,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,476.845907,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,133.404301,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,45.827978,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,17.796345,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,7.329518,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,124.842404,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,24.709570,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,7.679716,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,2.954664,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,1.283212,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.588190,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive S0A2 = bs.turbomolePrimitive(-0.007943,9.433006,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A2);
    Primitive S1A2 = bs.turbomolePrimitive(-0.071003,2.526244,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A2);
    Primitive S2A2 = bs.turbomolePrimitive(-0.178503,0.947368,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A2);
    Primitive S3A2 = bs.turbomolePrimitive(0.151064,0.424059,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A2);
    Primitive S4A2 = bs.turbomolePrimitive(0.735491,0.209845,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A2);
    Primitive S5A2 = bs.turbomolePrimitive(0.276059,0.108147,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A2);
    bs.add_state();
    Primitive P0A3 = bs.turbomolePrimitive(0.003760,124.842404,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A3);
    Primitive P0B3 = bs.turbomolePrimitive(0.003760,124.842404,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B3);
    Primitive P0C3 = bs.turbomolePrimitive(0.003760,124.842404,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C3);
    Primitive P1A3 = bs.turbomolePrimitive(0.037679,24.709570,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A3);
    Primitive P1B3 = bs.turbomolePrimitive(0.037679,24.709570,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B3);
    Primitive P1C3 = bs.turbomolePrimitive(0.037679,24.709570,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C3);
    Primitive P2A3 = bs.turbomolePrimitive(0.173897,7.679716,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A3);
    Primitive P2B3 = bs.turbomolePrimitive(0.173897,7.679716,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B3);
    Primitive P2C3 = bs.turbomolePrimitive(0.173897,7.679716,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C3);
    Primitive P3A3 = bs.turbomolePrimitive(0.418036,2.954664,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A3);
    Primitive P3B3 = bs.turbomolePrimitive(0.418036,2.954664,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B3);
    Primitive P3C3 = bs.turbomolePrimitive(0.418036,2.954664,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C3);
    Primitive P4A3 = bs.turbomolePrimitive(0.425860,1.283212,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A3);
    Primitive P4B3 = bs.turbomolePrimitive(0.425860,1.283212,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B3);
    Primitive P4C3 = bs.turbomolePrimitive(0.425860,1.283212,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C3);
    Primitive P5A3 = bs.turbomolePrimitive(0.101708,0.588190,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A3);
    Primitive P5B3 = bs.turbomolePrimitive(0.101708,0.588190,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B3);
    Primitive P5C3 = bs.turbomolePrimitive(0.101708,0.588190,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C3);
    bs.add_state();
    Primitive P0A4 = bs.turbomolePrimitive(-0.007139,9.433006,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A4);
    Primitive P0B4 = bs.turbomolePrimitive(-0.007139,9.433006,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B4);
    Primitive P0C4 = bs.turbomolePrimitive(-0.007139,9.433006,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C4);
    Primitive P1A4 = bs.turbomolePrimitive(-0.018293,2.526244,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A4);
    Primitive P1B4 = bs.turbomolePrimitive(-0.018293,2.526244,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B4);
    Primitive P1C4 = bs.turbomolePrimitive(-0.018293,2.526244,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C4);
    Primitive P2A4 = bs.turbomolePrimitive(0.076216,0.947368,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A4);
    Primitive P2B4 = bs.turbomolePrimitive(0.076216,0.947368,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B4);
    Primitive P2C4 = bs.turbomolePrimitive(0.076216,0.947368,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C4);
    Primitive P3A4 = bs.turbomolePrimitive(0.414510,0.424059,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A4);
    Primitive P3B4 = bs.turbomolePrimitive(0.414510,0.424059,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B4);
    Primitive P3C4 = bs.turbomolePrimitive(0.414510,0.424059,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C4);
    Primitive P4A4 = bs.turbomolePrimitive(0.488962,0.209845,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A4);
    Primitive P4B4 = bs.turbomolePrimitive(0.488962,0.209845,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B4);
    Primitive P4C4 = bs.turbomolePrimitive(0.488962,0.209845,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C4);
    Primitive P5A4 = bs.turbomolePrimitive(0.105882,0.108147,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A4);
    Primitive P5B4 = bs.turbomolePrimitive(0.105882,0.108147,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B4);
    Primitive P5C4 = bs.turbomolePrimitive(0.105882,0.108147,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C4);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:52 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:be STO-6G
void basisbank::add_STO6G_Be(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,312.870494,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,57.364463,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,16.048509,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,5.513096,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,2.140897,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,0.881739,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,13.633247,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,2.698375,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,0.838653,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,0.322660,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.140131,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.064233,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,13.633247,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,13.633247,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,13.633247,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,2.698375,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,2.698375,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,2.698375,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,0.838653,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,0.838653,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,0.838653,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,0.322660,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,0.322660,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,0.322660,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.140131,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.140131,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.140131,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.064233,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.064233,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.064233,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:53 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:c STO-6G
void basisbank::add_STO6G_C(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,742.737049,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,136.180025,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,38.098264,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,13.087782,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,5.082369,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,2.093200,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,30.497239,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,6.036200,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,1.876046,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,0.721783,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.313471,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.143687,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,30.497239,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,30.497239,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,30.497239,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,6.036200,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,6.036200,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,6.036200,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,1.876046,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,1.876046,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,1.876046,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,0.721783,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,0.721783,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,0.721783,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.313471,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.313471,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.313471,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.143687,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.143687,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.143687,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:56 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:na STO-6G
void basisbank::add_STO6G_Na(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,2600.756771,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,476.845907,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,133.404301,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,45.827978,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,17.796345,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,7.329518,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,124.842404,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,24.709570,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,7.679716,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,2.954664,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,1.283212,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.588190,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive S0A2 = bs.turbomolePrimitive(-0.007943,9.433006,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A2);
    Primitive S1A2 = bs.turbomolePrimitive(-0.071003,2.526244,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A2);
    Primitive S2A2 = bs.turbomolePrimitive(-0.178503,0.947368,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A2);
    Primitive S3A2 = bs.turbomolePrimitive(0.151064,0.424059,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A2);
    Primitive S4A2 = bs.turbomolePrimitive(0.735491,0.209845,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A2);
    Primitive S5A2 = bs.turbomolePrimitive(0.276059,0.108147,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A2);
    bs.add_state();
    Primitive P0A3 = bs.turbomolePrimitive(0.003760,124.842404,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A3);
    Primitive P0B3 = bs.turbomolePrimitive(0.003760,124.842404,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B3);
    Primitive P0C3 = bs.turbomolePrimitive(0.003760,124.842404,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C3);
    Primitive P1A3 = bs.turbomolePrimitive(0.037679,24.709570,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A3);
    Primitive P1B3 = bs.turbomolePrimitive(0.037679,24.709570,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B3);
    Primitive P1C3 = bs.turbomolePrimitive(0.037679,24.709570,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C3);
    Primitive P2A3 = bs.turbomolePrimitive(0.173897,7.679716,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A3);
    Primitive P2B3 = bs.turbomolePrimitive(0.173897,7.679716,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B3);
    Primitive P2C3 = bs.turbomolePrimitive(0.173897,7.679716,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C3);
    Primitive P3A3 = bs.turbomolePrimitive(0.418036,2.954664,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A3);
    Primitive P3B3 = bs.turbomolePrimitive(0.418036,2.954664,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B3);
    Primitive P3C3 = bs.turbomolePrimitive(0.418036,2.954664,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C3);
    Primitive P4A3 = bs.turbomolePrimitive(0.425860,1.283212,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A3);
    Primitive P4B3 = bs.turbomolePrimitive(0.425860,1.283212,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B3);
    Primitive P4C3 = bs.turbomolePrimitive(0.425860,1.283212,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C3);
    Primitive P5A3 = bs.turbomolePrimitive(0.101708,0.588190,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A3);
    Primitive P5B3 = bs.turbomolePrimitive(0.101708,0.588190,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B3);
    Primitive P5C3 = bs.turbomolePrimitive(0.101708,0.588190,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C3);
    bs.add_state();
    Primitive P0A4 = bs.turbomolePrimitive(-0.007139,9.433006,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A4);
    Primitive P0B4 = bs.turbomolePrimitive(-0.007139,9.433006,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B4);
    Primitive P0C4 = bs.turbomolePrimitive(-0.007139,9.433006,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C4);
    Primitive P1A4 = bs.turbomolePrimitive(-0.018293,2.526244,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A4);
    Primitive P1B4 = bs.turbomolePrimitive(-0.018293,2.526244,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B4);
    Primitive P1C4 = bs.turbomolePrimitive(-0.018293,2.526244,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C4);
    Primitive P2A4 = bs.turbomolePrimitive(0.076216,0.947368,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A4);
    Primitive P2B4 = bs.turbomolePrimitive(0.076216,0.947368,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B4);
    Primitive P2C4 = bs.turbomolePrimitive(0.076216,0.947368,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C4);
    Primitive P3A4 = bs.turbomolePrimitive(0.414510,0.424059,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A4);
    Primitive P3B4 = bs.turbomolePrimitive(0.414510,0.424059,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B4);
    Primitive P3C4 = bs.turbomolePrimitive(0.414510,0.424059,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C4);
    Primitive P4A4 = bs.turbomolePrimitive(0.488962,0.209845,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A4);
    Primitive P4B4 = bs.turbomolePrimitive(0.488962,0.209845,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B4);
    Primitive P4C4 = bs.turbomolePrimitive(0.488962,0.209845,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C4);
    Primitive P5A4 = bs.turbomolePrimitive(0.105882,0.108147,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A4);
    Primitive P5B4 = bs.turbomolePrimitive(0.105882,0.108147,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B4);
    Primitive P5C4 = bs.turbomolePrimitive(0.105882,0.108147,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C4);
}


//# 6-31++G EMSL Basis Set Exchange Library 10/17/14 4:47 AM
//# Elements References
//# -------- ----------
//# H - He: W.J. Hehre, R. Ditchfield and J.A. Pople, J. Chem. Phys. 56,
//# Li - Ne: 2257 (1972). Note: Li and B come from J.D. Dill and J.A.
//# Pople, J. Chem. Phys. 62, 2921 (1975).
//# Na - Ar: M.M. Francl, W.J. Petro, W.J. Hehre, J.S. Binkley, M.S. Gordon,
//# D.J. DeFrees and J.A. Pople, J. Chem. Phys. 77, 3654 (1982)
//# K - Zn: V. Rassolov, J.A. Pople, M. Ratner and T.L. Windus, J. Chem. Phys.
//# 109, 1223 (1998)
//# Note: He and Ne are unpublished basis sets taken from the Gaussian
//# program
//#
//# Elements Reference
//# -------- ----------
//# H, Li-Cl: T. Clark, J. Chandrasekhar, G.W. Spitznagel, P.V.R. Schleyer,
//# J. Comp. Chem. 4, 294 (1983).
//#
//Ignored the following information from file:h 6-31++G
void basisbank::add_631ppG_H(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.033495,18.731137,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.234727,2.825394,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.813757,0.640122,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(1.000000,0.161278,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    bs.add_state();
    Primitive S0A2 = bs.turbomolePrimitive(1.000000,0.036000,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:54 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:o STO-6G
void basisbank::add_STO6G_O(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,1355.584234,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,248.544885,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,69.533902,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,23.886772,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,9.275933,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,3.820341,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,52.187762,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,10.329320,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,3.210345,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,1.235135,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.536420,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.245881,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,52.187762,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,52.187762,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,52.187762,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,10.329320,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,10.329320,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,10.329320,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,3.210345,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,3.210345,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,3.210345,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,1.235135,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,1.235135,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,1.235135,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.536420,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.536420,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.536420,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.245881,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.245881,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.245881,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:54 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:n STO-6G
void basisbank::add_STO6G_N(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,1027.828458,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,188.451223,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,52.721861,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,18.111382,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,7.033180,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,2.896652,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,39.198808,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,7.758467,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,2.411326,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,0.927724,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.402911,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.184684,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,39.198808,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,39.198808,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,39.198808,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,7.758467,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,7.758467,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,7.758467,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,2.411326,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,2.411326,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,2.411326,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,0.927724,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,0.927724,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,0.927724,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.402911,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.402911,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.402911,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.184684,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.184684,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.184684,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:56 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:ne STO-6G
void basisbank::add_STO6G_Ne(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,2146.955475,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,393.641936,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,110.126828,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,37.831538,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,14.691093,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,6.050603,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,85.504429,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,16.923558,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,5.259829,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,2.023646,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.878871,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.402851,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,85.504429,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,85.504429,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,85.504429,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,16.923558,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,16.923558,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,16.923558,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,5.259829,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,5.259829,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,5.259829,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,2.023646,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,2.023646,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,2.023646,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.878871,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.878871,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.878871,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.402851,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.402851,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.402851,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


//# STO-6G EMSL Basis Set Exchange Library 10/17/14 4:51 AM
//# Elements References
//# -------- ----------
//# H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 51, 2657
//# (1969).
//# Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart and J.A. Pople,
//# J. Chem. Phys. 52, 2769 (1970).
//#
//Ignored the following information from file:li STO-6G
void basisbank::add_STO6G_Li(vec3 corePos){

    bs.add_state();
    Primitive S0A0 = bs.turbomolePrimitive(0.009164,167.175846,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A0);
    Primitive S1A0 = bs.turbomolePrimitive(0.049361,30.651508,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A0);
    Primitive S2A0 = bs.turbomolePrimitive(0.168538,8.575187,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A0);
    Primitive S3A0 = bs.turbomolePrimitive(0.370563,2.945808,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A0);
    Primitive S4A0 = bs.turbomolePrimitive(0.416492,1.143944,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A0);
    Primitive S5A0 = bs.turbomolePrimitive(0.130334,0.471139,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A0);
    bs.add_state();
    Primitive S0A1 = bs.turbomolePrimitive(-0.013253,6.597564,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S0A1);
    Primitive S1A1 = bs.turbomolePrimitive(-0.046992,1.305830,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S1A1);
    Primitive S2A1 = bs.turbomolePrimitive(-0.033785,0.405851,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S2A1);
    Primitive S3A1 = bs.turbomolePrimitive(0.250242,0.156146,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S3A1);
    Primitive S4A1 = bs.turbomolePrimitive(0.595117,0.067814,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S4A1);
    Primitive S5A1 = bs.turbomolePrimitive(0.240706,0.031084,0,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, S5A1);
    bs.add_state();
    Primitive P0A2 = bs.turbomolePrimitive(0.003760,6.597564,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0A2);
    Primitive P0B2 = bs.turbomolePrimitive(0.003760,6.597564,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0B2);
    Primitive P0C2 = bs.turbomolePrimitive(0.003760,6.597564,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P0C2);
    Primitive P1A2 = bs.turbomolePrimitive(0.037679,1.305830,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1A2);
    Primitive P1B2 = bs.turbomolePrimitive(0.037679,1.305830,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1B2);
    Primitive P1C2 = bs.turbomolePrimitive(0.037679,1.305830,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P1C2);
    Primitive P2A2 = bs.turbomolePrimitive(0.173897,0.405851,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2A2);
    Primitive P2B2 = bs.turbomolePrimitive(0.173897,0.405851,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2B2);
    Primitive P2C2 = bs.turbomolePrimitive(0.173897,0.405851,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P2C2);
    Primitive P3A2 = bs.turbomolePrimitive(0.418036,0.156146,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3A2);
    Primitive P3B2 = bs.turbomolePrimitive(0.418036,0.156146,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3B2);
    Primitive P3C2 = bs.turbomolePrimitive(0.418036,0.156146,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P3C2);
    Primitive P4A2 = bs.turbomolePrimitive(0.425860,0.067814,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4A2);
    Primitive P4B2 = bs.turbomolePrimitive(0.425860,0.067814,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4B2);
    Primitive P4C2 = bs.turbomolePrimitive(0.425860,0.067814,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P4C2);
    Primitive P5A2 = bs.turbomolePrimitive(0.101708,0.031084,1,0,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5A2);
    Primitive P5B2 = bs.turbomolePrimitive(0.101708,0.031084,0,1,0,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5B2);
    Primitive P5C2 = bs.turbomolePrimitive(0.101708,0.031084,0,0,1,corePos);
    bs.add_primitive_to_state(bs.Nstates-1, P5C2);
}


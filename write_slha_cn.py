#!/usr/bin/env python

import sys, os
from math import *

mN1 = sys.argv[1]
mN2 = sys.argv[2]
mN3 = sys.argv[3]
mN4 = sys.argv[4]
mC1 = sys.argv[5]
mC2 = sys.argv[6]
msf = sys.argv[7]

U11, U12, U21, U22 = sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11]
V11, V12, V21, V22 = sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15]

N11, N12, N13, N14 = sys.argv[16], sys.argv[17], sys.argv[18], sys.argv[19]
N21, N22, N23, N24 = sys.argv[20], sys.argv[21], sys.argv[22], sys.argv[23]
N31, N32, N33, N34 = sys.argv[24], sys.argv[25], sys.argv[26], sys.argv[27]
N41, N42, N43, N44 = sys.argv[28], sys.argv[29], sys.argv[30], sys.argv[31]

output = """

# SOFTSUSY3.6.1 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.6.1       # version number
Block MODSEL  # Select model
     1    0   # nonUniversal
Block SMINPUTS             # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # mb(mb)
     6    1.73300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    1.00000000e+01   # tanb, DRbar, Feynman gauge
Block EXTPAR               # non-universal SUSY breaking parameters
     0    -1.00000000e+00  # Set MX=MSUSY
     1     2.80000000e+02  # M_1(MX)
     2     2.90000000e+02  # M_2(MX)
     3     5.00000000e+03  # M_3(MX)
     11    0.00000000e+00  # At(MX)
     12    0.00000000e+00  # Ab(MX)
     13    0.00000000e+00  # Atau(MX)
     23    3.00000000e+02  # mu(MX)
     26    5.00000000e+03  # mA(pole)
     31    5.00000000e+03  # meL(MX)
     32    5.00000000e+03  # mmuL(MX)
     33    5.00000000e+03  # mtauL(MX)
     34    5.00000000e+03  # meR(MX)
     35    5.00000000e+03  # mmuR(MX)
     36    5.00000000e+03  # mtauR(MX)
     41    5.00000000e+03  # mqL1(MX)
     42    5.00000000e+03  # mqL2(MX)
     43    5.00000000e+03  # mqL3(MX)
     44    5.00000000e+03  # muR(MX)
     45    5.00000000e+03  # mcR(MX)
     46    5.00000000e+03  # mtR(MX)
     47    5.00000000e+03  # mdR(MX)
     48    5.00000000e+03  # msR(MX)
     49    5.00000000e+03  # mbR(MX)
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=6.17365905e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03828269e+01   # MW
        25     1.23721699e+02   # h0
        35     5.00001450e+03   # H0
        36     4.99999960e+03   # A0
        37     5.00106632e+03   # H+
   1000021     5.27583016e+03   # ~g
   1000022     {mN1}   # ~neutralino(1)
   1000023     {mN2}   # ~neutralino(2)
   1000025     {mN3}   # ~neutralino(3)
   1000035     {mN4}   # ~neutralino(4)
   1000024     {mC1}   # ~chargino(1)   
   1000037     {mC2}   # ~chargino(2)
   1000001     {msf}   # ~d_L
   1000002     {msf}   # ~u_L
   1000003     {msf}   # ~s_L
   1000004     {msf}   # ~c_L
   1000005     {msf}   # ~b_1
   1000006     {msf}   # ~t_1
   1000011     {msf}   # ~e_L
   1000012     {msf}   # ~nue_L
   1000013     {msf}   # ~mu_L
   1000014     {msf}   # ~numu_L
   1000015     {msf}   # ~stau_1
   1000016     {msf}   # ~nu_tau_L
   2000001     {msf}   # ~d_R
   2000002     {msf}   # ~u_R
   2000003     {msf}   # ~s_R
   2000004     {msf}   # ~c_R
   2000005     {msf}   # ~b_2
   2000006     {msf}   # ~t_2
   2000011     {msf}   # ~e_R
   2000013     {msf}   # ~mu_R
   2000015     {msf}   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05326799e-01       # alpha - evaluated at p^2=0
Block nmix                  # neutralino mixing matrix
  1  1     {N11}   # N_11
  1  2     {N12}   # N_12
  1  3     {N13}   # N_13
  1  4     {N14}   # N_14
  2  1     {N21}   # N_21
  2  2     {N22}   # N_22
  2  3     {N23}   # N_23
  2  4     {N24}   # N_24
  3  1     {N31}   # N_31
  3  2     {N32}   # N_32
  3  3     {N33}   # N_33
  3  4     {N34}   # N_34
  4  1     {N41}   # N_41
  4  2     {N42}   # N_42
  4  3     {N43}   # N_43
  4  4     {N44}   # N_44
Block Umix                  # chargino U mixing matrix 
  1  1    {U11}    # U_11
  1  2    {U12}    # U_12
  2  1    {U21}    # U_21
  2  2    {U22}    # U_22
Block Vmix                  # chargino V mixing matrix 
  1  1    {V11}    # V_11
  1  2    {V12}    # V_12
  2  1    {V21}    # V_21
  2  2    {V22}    # V_22
Block stopmix               # stop mixing matrix
  1  1     1.15156258e-01   # F_11
  1  2     9.93347390e-01   # F_12
  2  1     9.93347390e-01   # F_21
  2  2    -1.15156258e-01   # F_22
Block sbotmix               # sbottom mixing matrix
  1  1     9.99696809e-01   # F_11
  1  2     2.46229472e-02   # F_12
  2  1    -2.46229472e-02   # F_21
  2  2     9.99696809e-01   # F_22
Block staumix               # stau mixing matrix
  1  1     1.36416555e-01   # F_11
  1  2     9.90651565e-01   # F_12
  2  1     9.90651565e-01   # F_21
  2  2    -1.36416555e-01   # F_22
Block gauge Q= 5.00173615e+03  # SM gauge couplings
     1     3.66174627e-01   # g'(Q)MSSM DRbar
     2     6.41230998e-01   # g(Q)MSSM DRbar
     3     9.83379571e-01   # g3(Q)MSSM DRbar
Block yu Q= 5.00173615e+03  
  3  3     8.22475532e-01   # Yt(Q)MSSM DRbar
Block yd Q= 5.00173615e+03  
  3  3     1.28440097e-01   # Yb(Q)MSSM DRbar
Block ye Q= 5.00173615e+03  
  3  3     9.80771750e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 5.00173615e+03 # Higgs mixing parameters
     1     3.00000000e+02    # mu(Q)MSSM DRbar
     2     9.46875110e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.40119137e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     2.48604720e+07    # mA^2(Q)MSSM DRbar
Block msoft Q= 5.00173615e+03  # MSSM DRbar SUSY breaking parameters
     1     2.80000000e+02      # M_1(Q)
     2     2.90000000e+02      # M_2(Q)
     3     5.00000000e+03      # M_3(Q)
    21     2.45919951e+07      # mH1^2(Q)
    22     8.37809441e+05      # mH2^2(Q)
    31     5.00000000e+03      # meL(Q)
    32     5.00000000e+03      # mmuL(Q)
    33     5.00000000e+03      # mtauL(Q)
    34     5.00000000e+03      # meR(Q)
    35     5.00000000e+03      # mmuR(Q)
    36     5.00000000e+03      # mtauR(Q)
    41     4.99999997e+03      # mqL1(Q)
    42     4.99999997e+03      # mqL2(Q)
    43     4.99999996e+03      # mqL3(Q)
    44     4.99999997e+03      # muR(Q)
    45     4.99999997e+03      # mcR(Q)
    46     4.99999994e+03      # mtR(Q)
    47     4.99999997e+03      # mdR(Q)
    48     4.99999997e+03      # msR(Q)
    49     4.99999997e+03      # mbR(Q)
Block au Q= 5.00173615e+03  
  1  1     2.91282396e-05      # Au(Q)MSSM DRbar
  2  2     2.91286552e-05      # Ac(Q)MSSM DRbar
  3  3     5.09381165e-05      # At(Q)MSSM DRbar
Block ad Q= 5.00173615e+03  
  1  1     7.57957895e-06      # Ad(Q)MSSM DRbar
  2  2     7.58001385e-06      # As(Q)MSSM DRbar
  3  3     1.49253965e-05      # Ab(Q)MSSM DRbar
Block ae Q= 5.00173615e+03  
  1  1     7.12153835e-07      # Ae(Q)MSSM DRbar
  2  2     7.12162891e-07      # Amu(Q)MSSM DRbar
  3  3     7.14829898e-07      # Atau(Q)MSSM DRbar""".format( 
    mN1=mN1, mN2=mN2, mN3=mN3, mN4=mN4, 
    mC1=mC1, mC2=mC2, msf=msf,
    U11=U11,  U12=U12,  U21=U21,  U22=U22, 
    V11=V11,  V12=V12,  V21=V21,  V22=V22,     
    N11=N11,  N12=N12,  N13=N13,  N14=N14, 
    N21=N21,  N22=N22,  N23=N23,  N24=N24, 
    N31=N31,  N32=N32,  N33=N33,  N34=N34, 
    N41=N41,  N42=N42,  N43=N43,  N44=N44, 
  )                

outFile = open('slha.out','w')
outFile.write(output)




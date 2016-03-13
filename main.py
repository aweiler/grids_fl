#! /usr/bin/env python

import os, sys, random
import numpy as np
from numpy import linalg as la
from numpy import random as rn
import iminuit as im
from math import *
import subprocess as sbproc

FNULL = open(os.devnull, 'w')

def get_line(List):
    line = str(List[0])
    for i in range(1, len(List)):
        line += '  ' + str(List[i])
    return line

def get_xsec_from_prospino(mNneu, mCha, msf, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch):

    ep = 10**(-15)
    ep = 0
    # if U[0,0] == 0: U[0,0] = ep
    # if V[0,0] == 0: V[0,0] = ep
    # if N[0,0] == 0: N[0,0] = ep
    # if N[0,1] == 0: N[0,1] = N[1,0] = ep

    for i in xrange(4):
        for j in xrange(4):
            #x = random.uniform(-1, 1)
            if N[i,j] == 0: N[i,j] = ep

    for i in xrange(2):
        for j in xrange(2):
            #x = random.uniform(-1, 1)            
            if U[i,j] == 0: U[i,j] = ep
            #x = random.uniform(-1, 1)
            if V[i,j] == 0: V[i,j] = ep

    arguments = [
                str(mNeu[0]), str(mNeu[1]), str(mNeu[2]), str(mNeu[3]),
                str(mCha[0]), str(mCha[1]), str(msf),
                str(U[0,0]), str(U[0,1]), str(U[1,0]), str(U[1,1]),
                str(V[0,0]), str(V[0,1]), str(V[1,0]), str(V[1,1]),                     
                str(N[0,0]), str(N[0,1]), str(N[0,2]), str(N[0,3]),
                str(N[1,0]), str(N[1,1]), str(N[1,2]), str(N[1,3]),
                str(N[2,0]), str(N[2,1]), str(N[2,2]), str(N[2,3]),
                str(N[3,0]), str(N[3,1]), str(N[3,2]), str(N[3,3]),
                ]

    input_array = ['python', 'write_slha_cn.py'] + arguments
    sbproc.call( input_array )

    input_argument = ['slha.out', str(i1_pros), str(i2_pros)]
    input_argument += [sqrtS, nlo_string, error_switch]    
    sbproc.call( ['sh', 'run_prospino.sh'] + input_argument, stdout=FNULL )
    #sbproc.call( ['sh', 'run_prospino.sh'] + input_argument )

    #sbproc.call( ['sh', 'run_prospino.sh', 'slha.out', str(i1_pros), str(i2_pros)] )

    #sbproc.call( ['cp', 'slha.out', 'slha_tmp/' + str(i) + '.out'], stdout=FNULL )

    f = open('prospino.out','r')
    if nlo_string == 'LO': column = 14
    if nlo_string == 'NLO': column = 15
    xval = []
    if error_switch == '0':
        xval.append(float( f.readline().split()[column]))
    elif error_switch == '1':
        xval.append(float( f.readline().split()[column]))
        xval.append(float( f.readline().split()[column]))
        xval.append(float( f.readline().split()[column]))        
    else:
        print 'ERROR!! error_switch has to be 0 or 1!'
        exit()

    return xval 

def get_zero_matrices():

    U = np.array([[0.,0.],[0.,0.]])
    V = np.array([[0.,0.],[0.,0.]])

    N = np.array([
        [0., 0., 0., 0.], 
        [0., 0., 0., 0.], 
        [0., 0., 0., 0.], 
        [0., 0., 0., 0.]
        ])

    return U, V, N

def get_N_matrix():

    a0, a1 = 0., 2.*pi
    nsub = []
    ang = []
    for i in range(6): 
        nsub.append(np.identity(4))
        a = a0 + a1 * rn.rand(1)[0]
        ang.append( a )

    ic = 0
    for i in range(3):
        for j in range(i+1, 4):
            #print i, j, ic
            a = ang[ic]
            nmat = nsub[ic]
            nmat[i,i] = cos(a)
            nmat[i,j] = sin(a)
            nmat[j,i] = -sin(a)
            nmat[j,j] = cos(a)
            nsub[ic] = nmat            
            ic += 1

    Nmat = np.identity(4)
    for i in range(6):
        Nmat = np.dot( Nmat, nsub[i] )

    return Nmat

def get_2by2( ang ):
    o11 = cos(ang)
    o12 = sin(ang)
    o22 = o11
    o21 = -o12
    return np.array( [[o11, o12], 
                      [o21, o22]] )

def get_VU():
    a1, a2 = 2. * pi * rn.rand(2)
    V = get_2by2(a1)
    U = get_2by2(a2)
    return V, U


sw = sqrt(1. - 8.02463984**2 / 9.11876000**2)
cw = sqrt(1 - sw**2)
r2 = sqrt(2.)

deno = sw * cw
Lu = ( 1./2. - 2./3. * sw**2 ) / deno
Ld = (-1./2 - (-1./3.) * sw**2 ) / deno
Ru = (-2./3. * sw**2 ) / deno
Rd = ( - (-1./3.) * sw**2 ) / deno


def Lut(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (2./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Lu

def Ldt(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (-1./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Ld

def Rut(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (2./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Ru

def Rdt(i, n): 
    return (n[i, 0]*cw + n[i, 1]*sw) * (-1./3.) + (-n[i, 0]*sw + n[i, 1]*cw) * Rd

def gNC( i, j, n, v, u ):
    elem1 = n[i, 1] * v[j, 0] - n[i, 3] * v[j, 1]/sqrt(2.) 
    elem2 = v[j, 0] * Lut(i, n)
    elem3 = n[i, 1] * u[j, 0] + n[i, 2] * u[j, 1]/sqrt(2.) 
    elem4 = u[j, 0] * Ldt(i, n)
    return [elem1, elem2, elem3, elem4]

def cNC( v, c ):
    if c == 0: return v[0]**2 + v[2]**2
    if c == 1: return v[1]**2 + v[3]**2
    if c == 2: return 2.*( v[0]*v[1] - v[2]*v[3] )
    if c == 3: return 2.*( v[1]*v[2] - v[0]*v[3] )
    if c == 4: return 2.* v[0]*v[2]
    if c == 5: return 2.* v[1]*v[3]

def gNN( i, j, n ):
    elem1 = n[i, 3] * n[j, 3] - n[i, 2] * n[j, 2]
    elem2 = Lut(i, n) * Lut(j, n)
    elem3 = Rut(i, n) * Rut(j, n)
    elem4 = Ldt(i, n) * Ldt(j, n)
    elem5 = Rdt(i, n) * Rdt(j, n)
    return [elem1, elem2, elem3, elem4, elem5]

def cNN( v, c ):
    if c == 0: return v[0]**2
    if c == 1: return v[1]**2 + v[2]**2
    if c == 2: return 2.*( Lu * v[0] * v[1] + Ru * v[0] * v[2] )
    if c == 3: return v[3]**2 + v[4]**2
    if c == 4: return 2.*( Ld * v[0] * v[3] + Rd * v[0] * v[4] )    


def gCC( i, j, u, v ):
    elem1 = u[i, 0] * u[j, 0]
    elem2 = v[i, 0] * v[j, 0]
    elem3 = u[i, 1] * u[j, 1]
    elem4 = v[i, 1] * v[j, 1]
    return [elem1, elem2, elem3, elem4]

def cCC( v, c ):
    if c == 0: return v[0]**2
    if c == 1: return 2. * v[0] * v[1]
    if c == 2: return 2. * v[0] * v[2]
    if c == 3: return 2. * v[0] * v[3]
    if c == 4: return v[1]**2
    if c == 5: return 2. * v[1] * v[2]
    if c == 6: return 2. * v[1] * v[3]
    if c == 7: return v[2]**2 + v[3]**2
    if c == 8: return 2. * v[2] * v[3]
    if c == 9: return 1.
    if c == 10: return v[0] - v[1] 
    if c == 11: return v[2] + v[3]



#####################################################
#####################################################
#####################################################

try:
    run_mode = sys.argv[1]
except: 
    run_mode = ''
possible_modes = ['CCsame', 'NNsame', 'C2+C1-', 'C2-C1+', 'NN', 'NC+', 'NC-']
if run_mode not in possible_modes:
    print 'the first argument has to be one of', possible_modes
    exit()

input_error = 0
sqrtS, nlo_string, error_switch, mass_string = '','','',''
if 'same' in run_mode:
    try:
        m1 = int(sys.argv[2])
        mQ = int(sys.argv[3])
        sqrtS = sys.argv[4]
        nlo_string = sys.argv[5]
        error_switch = sys.argv[6]
        mass_string = str(m1) +'_'+ str(mQ)
    except:
        input_error = 1
elif 'NC' in run_mode:
    try:
        mN = int(sys.argv[2])
        mC = int(sys.argv[3])    
        mQ = int(sys.argv[4])
        sqrtS = sys.argv[5]
        nlo_string = sys.argv[6]
        error_switch = sys.argv[7]
        mass_string = str(mN) +'_'+ str(mC) +'_'+ str(mQ)        
    except:
        input_error = 1
else:
    try:
        mL = int(sys.argv[2])
        mS = int(sys.argv[3])    
        mQ = int(sys.argv[4])
        sqrtS = sys.argv[5]
        nlo_string = sys.argv[6]
        error_switch = sys.argv[7]
        mass_string = str(mL) +'_'+ str(mS) +'_'+ str(mQ)        
    except:
        input_error = 1        

print 'sqrtS = ' + sqrtS
print 'nlo_string = ' + nlo_string
print 'error_switch = ' + error_switch

if input_error == 1:
    print "Input Error:"
    if 'same' in run_mode: 
        print "[run_mode] [m1] [mQ] [sqrtS (TeV)] [LO/NLO] [1/0 (with/without scale uncertaintly)]"
    elif 'NC' in run_mode: 
        print "[run_mode] [mN] [mC] [mQ] [sqrtS (TeV)] [LO/NLO] [1/0 (with/without scale uncertaintly)]"    
    else:
        print "[run_mode] [mL] [mS] [mQ] [sqrtS (TeV)] [LO/NLO] [1/0 (with/without scale uncertaintly)]"        
    exit()

if nlo_string not in ['LO','NLO']: 
    print 'ERROR: nlo_string should be LO or NLO'
    exit()
if error_switch not in ['0','1']: 
    print 'ERROR: error_switch should be 0 or 1 (LO or NLO)'
    exit()

outfile = "output_{run_mode}_{mass_string}_{sqrtS}_{nlo_string}_{error_switch}.dat".format(run_mode=run_mode, mass_string=mass_string, sqrtS=sqrtS, nlo_string=nlo_string, error_switch=error_switch)
f = open(outfile,'w')

#####################################################
#####################################################
#####################################################

mdec = 500000
mCha = [mdec, mdec]
mNeu = [mdec, mdec, mdec, mdec]    
nsamp = 12
if run_mode == 'NC+':
    i1, i2 = 0, 0
    i1_pros, i2_pros = i1 + 1, i2 + 5 ## C+ ##
    mCha = [mC, mdec]
    mNeu = [mN, mdec, mdec, mdec]        
if run_mode == 'NC-':
    i1, i2 = 0, 0
    i1_pros, i2_pros = i1 + 1, i2 + 7 ## C- ##
    mCha = [mC, mdec]
    mNeu = [mN, mdec, mdec, mdec]            
if run_mode == 'NN':
    i1, i2 = 0, 1
    i1_pros, i2_pros = i1+1, i2+1    
    mNeu = [mS, mL, mdec, mdec]

mass_lines  = "mC = " + get_line(mCha) + '\n'
mass_lines += "mN = " + get_line(mNeu) + '\n'

run_info = """
#######################################################
{run_mode}, {sqrtS}TeV, {nlo_string}  
""".format(run_mode=run_mode, nlo_string=nlo_string, sqrtS=sqrtS)
run_info += mass_lines
run_info += """mQ = {mQ} 
i1, i2 = {i1}, {i2}
i1, i2 (prospino) = {i1_pros}, {i2_pros}
#######################################################
""".format(mQ=str(mQ),i1=str(i1),i2=str(i2),i1_pros=str(i1_pros),i2_pros=str(i2_pros))
print run_info

# i1_pros, i2_pros = '5', '8'
# i1, i2 = 0, 1
# sqrtS = '8'
# nlo_string = 'LO'
# error_switch = '0'
# mQ = 500
# mNeu = [100, 200, 300, 400]
# mCha = [110, 110]

if error_switch == '0': mu_list = ['1.0']
if error_switch == '1': mu_list = ['0.5', '1.0', '2.0']
n_scale = len(mu_list)

if run_mode in ['NC+', 'NC-']:

    a0, a1, a2, a3, a4, a5 = [], [], [], [], [], []

    ##################################################################
    ###                a0
    ##################################################################

    U, V, N = get_zero_matrices()

    N[0,3], N[3,0] = r2, r2
    V[0,1], V[1,0] = 1., 1.

    xsec = get_xsec_from_prospino(mNeu, mCha, mQ, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)
    a0 = xsec
    #xsec = a0 = 1.58071
    #print gNC(0, 0, N, V, U)
    #print 'a0 =', a0

    ##################################################################
    ###                a1
    ##################################################################

    U, V, N = get_zero_matrices()

    n11 = 1./(cw*(2./3.) - (1./2. - sw**2 * (2./3.) )/cw)

    N[0,0] = n11
    V[0,0] = 1.

    xsec = get_xsec_from_prospino(mNeu, mCha, mQ, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)
    a1 = xsec
    #xsec = a1 = 0.0572713
    #print gNC(0, 0, N, V, U)
    print 'a1 =', a1

    ##################################################################
    ###                a2
    ##################################################################

    U, V, N = get_zero_matrices()
    x = random.uniform(-1, 1)

    N[0,0], N[0,1], N[1,0] = x, 1., 1.
    V[0,0] = 1.

    xsec = get_xsec_from_prospino(mNeu, mCha, mQ, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)
    gnc = gNC(0, 0, N, V, U)

    for i in xrange(n_scale):
        a2__ = (xsec[i] - a0[i] * cNC(gnc, 0) - a1[i] * cNC(gnc, 1))/cNC(gnc, 2)
        a2.append( a2__ )

    #print gNC(0, 0, N, V, U)
    print 'a2 =', a2

    ##################################################################
    ###                a3
    ##################################################################

    U, V, N = get_zero_matrices()
    x = random.uniform(-1, 1)

    n11 = 1./(-cw/3. - (-1./2. + sw**2/3.)/cw)

    N[0,0], N[0,3], N[3,0] = n11, x*r2, x*r2
    V[0,1], V[1,0] = 1., 1.
    U[0,0] = 1.

    xsec = get_xsec_from_prospino(mNeu, mCha, mQ, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)
    gnc = gNC(0, 0, N, V, U)

    for i in xrange(n_scale):
        a3__ = (xsec[i] - a0[i] * cNC(gnc, 0) - a1[i] * cNC(gnc, 1))/cNC(gnc, 3)
        a3.append(a3__)

    #print gNC(0, 0, N, V, U)
    print 'a3 =', a3

    ##################################################################
    ###                a4
    ##################################################################

    U, V, N = get_zero_matrices()
    x = random.uniform(-1, 1)

    N[0,2] = N[2,0] = r2
    N[0,3] = N[3,0] = -x*r2
    V[0,1] = V[1,0] = 1.
    U[0,1] = U[1,0] = 1.

    xsec = get_xsec_from_prospino(mNeu, mCha, mQ, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)
    gnc = gNC(0, 0, N, V, U)

    for i in xrange(n_scale):
        a4__ = (xsec[i] - a0[i] * cNC(gnc, 0))/cNC(gnc, 4)
        a4.append( a4__ )

    #print gNC(0, 0, N, V, U)
    print 'a4 =', a4

    ##################################################################
    ###                a5
    ##################################################################

    U, V, N = get_zero_matrices()
    x = random.uniform(-1, 1)

    n11 = 1./(-cw/3. - (-1./2. + sw**2/3.)/cw)

    N[0,0] = n11

    V[0,0] = 1.
    U[0,0] = x

    xsec = get_xsec_from_prospino(mNeu, mCha, mQ, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)
    gnc = gNC(0, 0, N, V, U)

    for i in xrange(n_scale):
        a5__ = (xsec[i] - a1[i] * cNC(gnc, 1))/cNC(gnc, 5)
        a5.append( a5__ )
    #print gNC(0, 0, N, V, U)
    print 'a5 =', a5

    ##################################################################
    ###                Results
    ##################################################################

    Coeff = []
    for i in xrange(n_scale):
        Coeff.append([a0[i], a1[i], a2[i], a3[i], a4[i]])

    result = get_line(mass_string.split('_')) + '\n'
    result += '---- result ----  \n'
    for i in range(len(mu_list)): 
        result += 'scale=' + mu_list[i] + '\n'
        result += get_line(Coeff[i]) + '\n'
    result += '---------------- \n'
    print result
    f.write(result)

    ##################################################################
    ###                test
    ##################################################################
    for itest in xrange(2):

        N = get_N_matrix()
        V, U = get_VU()
        xsec_pros = get_xsec_from_prospino(mNeu, mCha, mQ, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)

        gnc = gNC(0, 0, N, V, U)
        xsec_predic = []
        for i in xrange(n_scale):
            xsec =  a0[i] * cNC(gnc, 0) + \
                    a1[i] * cNC(gnc, 1) + \
                    a2[i] * cNC(gnc, 2) + \
                    a3[i] * cNC(gnc, 3) + \
                    a4[i] * cNC(gnc, 4) + \
                    a5[i] * cNC(gnc, 5) 
            xsec_predic.append(xsec)

        outline = '#---  Test  ---#'
        outline = 'iTest'.rjust(4) + str('prospino').rjust(15) + 'prediction'.rjust(25) + 'ratio'.rjust(25) + 'scale'.rjust(10) + '\n'
        for i in range(len(mu_list)):
            mu = mu_list[i]
            outline += str(itest).rjust(4) + str(xsec_pros[i]).rjust(15) + str(xsec_predic[i]).rjust(25) + str(xsec_predic[i]/xsec_pros[i]).rjust(25) + mu.rjust(10) + '\n'
        print outline
        f.write(outline)
    f.write(run_info)

exit()

##################################################################

#! /usr/bin/env python

import os, sys
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

def get_xsec_from_prospino(mNneu, mCha, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch):

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
    callOK = sbproc.check_call( ['sh', 'run_prospino.sh'] + input_argument )
    print "Subprocess finished with value " , callOK
    # sbproc.call( ['sh', 'run_prospino.sh'] + input_argument, stdout=FNULL ) # original Kazuki
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


def get_Mline(run_mode, i1, i2):

    if run_mode in ['CCsame', 'C2+C1-', 'C2-C1+']:
        data = {}
        V, U = get_VU()
        elem = []
        elem.append( V[0,i1] * V[0,i2] )
        elem.append( U[0,i1] * U[0,i2] )
        elem.append( V[1,i1] * V[1,i2] )
        elem.append( U[1,i1] * U[1,i2] )

        Mline = []
        for i in range(4):
            for j in range(i, 4):
                Mline.append( elem[i]*elem[j] )
        data['Mline'] = Mline
        data['V'] = V
        data['U'] = U
        data['N'] = np.identity(4)
        return data

    if run_mode in ['NNsame', 'NN']:

        data = {}
        N = get_N_matrix()

        elem = []
        elem.append( N[i1,0] * N[i2,0] )  # B-B
        elem.append( N[i1,0] * N[i2,1] + N[i1,1] * N[i2,0] )  # B-W        
        elem.append( N[i1,1] * N[i2,1] )  # W-W
        elem.append( N[i1,2] * N[i2,2] - N[i1,3] * N[i2,3] )  # Hd-Hd

        Mline = []
        nmax = len(elem)
        for i in range(nmax):
            for j in range(i, nmax):
                Mline.append( elem[i]*elem[j] )
        data['Mline'] = Mline
        data['V'] = np.identity(2)
        data['U'] = np.identity(2)
        data['N'] = N
        return data

    if run_mode in ['NC+', 'NC-']:

        data = {}

        V, U = get_VU()
        N = get_N_matrix()

        elem = []
        elem.append( N[i1,0] * V[0,i2] )  # B-W
        elem.append( N[i1,0] * U[0,i2] )  # B-W        
        elem.append( N[i1,1] * V[0,i2] )  # W-W
        elem.append( N[i1,1] * U[0,i2] )  # W-W
        elem.append( N[i1,2] * U[1,i2] )  # Hd-Hd        
        elem.append( N[i1,3] * V[1,i2] )  # Hu-Hu

        Mline = []
        nmax = len(elem)
        for i in range(nmax):
            for j in range(i, nmax):
                Mline.append( elem[i]*elem[j] )
        data['Mline'] = Mline
        data['V'] = V
        data['U'] = U
        data['N'] = N
        return data

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

mdec = 2000
msf = mQ
mCha = [mdec, mdec]
mNeu = [mdec, mdec, mdec, mdec]    
nsamp = 20
if run_mode == 'CCsame':
    i1, i2 = 0, 0   # <== C1, C1 production
    i1_pros, i2_pros = i1+5, i2+7
    mCha = [m1, mdec]
if run_mode == 'NNsame':
    i1, i2 = 0, 0
    i1_pros, i2_pros = i1+1, i2+1
    mNeu = [m1, mdec, mdec, mdec]    
if run_mode == 'NN':
    i1, i2 = 0, 1
    i1_pros, i2_pros = i1+1, i2+1    
    mNeu = [mS, mL, mdec, mdec]    
if run_mode == 'NC+':
    i1, i2 = 0, 0
    i1_pros, i2_pros = i1 + 1, i2 + 5 ## C+ ##
    mCha = [mC, mdec]
    mNeu = [mN, mdec, mdec, mdec]        
    nsamp = 50
if run_mode == 'NC-':
    i1, i2 = 0, 0
    i1_pros, i2_pros = i1 + 1, i2 + 7 ## C- ##
    mCha = [mC, mdec]
    mNeu = [mN, mdec, mdec, mdec]            
    nsamp = 50
if run_mode == 'C2+C1-':
    i1, i2 = 0, 1   
    i1_pros, i2_pros = i1+7, i2+5 # <== C1-, C2+ production
    mCha = [mS, mL]
if run_mode == 'C2-C1+':
    i1, i2 = 0, 1   
    i1_pros, i2_pros = i1+5, i2+7 # <== C1+, C2- production
    mCha = [mS, mL]

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
f.write(run_info)

Coeff = []
#(mC1, mQ) = (500, 3000)
#Coeff = np.array([0.00118747693324, 0.00177949336646, 0.00136938810036, 0.000942739583667, 0.00110269896945, 0.000873844048213, 0.00135227791823, 0.000498039005221, 0.000372399410175, 0.000538524411907, ])
#(mC1, mQ) = (300, -500, 1000)
#Coeff = np.array([0.000277325009582, 0.000123795840113, -0.000277325009582, -0.000123795840113, 0.000326157067479, -0.000123795840113, -0.000326157067479, 0.000277325009582, 0.000123795840113, 0.000326157067479, ])
if len(Coeff) == 0:

    isamp = 0
    samp = []
    Mat = []
    while len(samp) < nsamp:

        data = get_Mline(run_mode, i1, i2)
        Mat.append(data['Mline'])
        samp.append(data)

    Mat = np.array(Mat)
    #print np.shape(Mat)
    #print Mat

    ######################################
    #    Obtaining Cross Section
    ######################################
    xsec_05 = np.array([])
    xsec_10 = np.array([])
    xsec_20 = np.array([])    
    print 'Run Prospino'
    for i in range(nsamp):

        V = samp[i]['V']
        U = samp[i]['U']
        N = samp[i]['N']

        f.write('---- sample {isamp} ---- \n'.format(isamp=i))
        f.write('V = \n')
        np.savetxt(f, V)
        f.write('U = \n')        
        np.savetxt(f, U)
        f.write('N = \n')                
        np.savetxt(f, N)
        f.write('Mline = \n')                
        np.savetxt(f, samp[i]['Mline'])

        xval = get_xsec_from_prospino(mNeu, mCha, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)

        f.write('prospino xsec: \n')                        
        if error_switch == '1':
            xsec_05 = np.append(xsec_05, xval[0] )
            xsec_10 = np.append(xsec_10, xval[1] )
            xsec_20 = np.append(xsec_20, xval[2] )
            f.write(str(xval[0]) + '     --- 0.5 \n')                        
            f.write(str(xval[1]) + '     --- 1.0 \n')                        
            f.write(str(xval[2]) + '     --- 2.0 \n')                        
        if error_switch == '0':
            xsec_10 = np.append(xsec_10, xval[0] )
            f.write(str(xval[0]) + '     --- 1.0 \n')                                    
        print i, xval
    if error_switch == '1': xsec_list = [xsec_05, xsec_10, xsec_20]
    if error_switch == '0': xsec_list = [xsec_10]


    f.write('======================== \n')                        

    ######################################
    #    Obtaining Coeff
    ######################################
    #Minv = la.inv( Mat )
    #Coeff = np.dot(Minv, xsec)

    xsec_repro = []
    chi2_result = []
    Coeff = []
    iscale = 0
    for xsec in xsec_list:

        f.write('Mat {isc} = \n'.format(isc=iscale))
        np.savetxt(f, Mat)
        f.write('xsec {isc} = \n'.format(isc=iscale))
        np.savetxt(f, xsec)
        iscale += 1

        ######################################
        if 'NC' in run_mode:

            def chi2(F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21):
                __coeff = np.array([F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21])
                MR = (np.dot( Mat, __coeff ) - xsec) / xsec
                ML = np.transpose(MR)
                return np.dot(ML, MR)

            m = im.Minuit(chi2, 
                F1=0, F2=0, F3=0, F4=0, F5=0, F6=0, F7=0, F8=0, F9=0, F10=0,
                F11=0, F12=0, F13=0, F14=0, F15=0, F16=0, F17=0, F18=0, F19=0, F20=0,
                F21=0,
                error_F1=0.001, 
                error_F2=0.001, 
                error_F3=0.001, 
                error_F4=0.001, 
                error_F5=0.001, 
                error_F6=0.001, 
                error_F7=0.001, 
                error_F8=0.001, 
                error_F9=0.001, 
                error_F10=0.001, 
                error_F11=0.001, 
                error_F12=0.001, 
                error_F13=0.001, 
                error_F14=0.001, 
                error_F15=0.001, 
                error_F16=0.001, 
                error_F17=0.001, 
                error_F18=0.001, 
                error_F19=0.001, 
                error_F20=0.001, 
                error_F21=0.001, 
                print_level=1) 

            m.migrad()
            vals = m.values

            F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21 = vals['F1'], vals['F2'], vals['F3'], vals['F4'], vals['F5'], vals['F6'], vals['F7'], vals['F8'], vals['F9'], vals['F10'], vals['F11'], vals['F12'], vals['F13'], vals['F14'], vals['F15'], vals['F16'], vals['F17'], vals['F18'], vals['F19'], vals['F20'], vals['F21']
            chi2_result.append(chi2(F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21)) 
            __coeff = np.array([F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21])
        ######################################
        else:
            def chi2(F1, F2, F3, F4, F5, F6, F7, F8, F9, F10):
                __coeff = np.array([F1, F2, F3, F4, F5, F6, F7, F8, F9, F10])
                MR = (np.dot( Mat, __coeff ) - xsec) / xsec
                ML = np.transpose(MR)
                return np.dot(ML, MR)

            m = im.Minuit(chi2, 
                F1=0, F2=0, F3=0, F4=0, F5=0, F6=0, F7=0, F8=0, F9=0, F10=0,
                error_F1=0.001, 
                error_F2=0.001, 
                error_F3=0.001, 
                error_F4=0.001, 
                error_F5=0.001, 
                error_F6=0.001, 
                error_F7=0.001, 
                error_F8=0.001, 
                error_F9=0.001, 
                error_F10=0.001, 
                print_level=1) 

            m.migrad()
            vals = m.values

            F1, F2, F3, F4, F5, F6, F7, F8, F9, F10 = vals['F1'], vals['F2'], vals['F3'], vals['F4'], vals['F5'], vals['F6'], vals['F7'], vals['F8'], vals['F9'], vals['F10']
            chi2_result.append( chi2(F1, F2, F3, F4, F5, F6, F7, F8, F9, F10) )
            __coeff = np.array([F1, F2, F3, F4, F5, F6, F7, F8, F9, F10])
        ######################################

        Coeff.append( __coeff )
        xsec_repro.append( np.dot(Mat, __coeff) )

    ######################################
    ######################################

    if error_switch == '0': mu_list = ['1.0']
    if error_switch == '1': mu_list = ['0.5', '1.0', '2.0' ]

    result = get_line(mass_string.split('_')) + '\n'
    for i in range(len(mu_list)): 
        result += 'scale=' + mu_list[i] + '\n'
        result += get_line(Coeff[i]) + '\n'
    print result
    f.write(result)

    outline = ''
    for i in range(len(mu_list)):
        mu = mu_list[i]
        outline += '#######################################################' + '\n'
        if error_switch == '0': outline += 'LO' + '\n'
        if error_switch == '1': outline += 'NLO  ' + mu + '\n'          
        outline += "Chi2 = " + str(chi2_result[i]) + '\n'   
        outline += "Coeff = np.array([" + get_line(Coeff[i]) + "])" + '\n'
        outline += "xsec(pros) = " + get_line(xsec_list[i]) + '\n'   
        outline += "xsec(repro) = " + get_line(xsec_repro[i]) + '\n'
        outline += "ratio = " + get_line(xsec[i] / xsec_repro[i]) + '\n'
        outline += '#######################################################' + '\n'
    print outline
    f.write(outline)

#####################################################

ntest = 1
for itest in range(ntest):

    data = get_Mline(run_mode, i1, i2)
    U = data['U']    
    V = data['V']
    N = data['N']
    
    ######################
    # prospino
    xsec_pros = get_xsec_from_prospino(mNeu, mCha, U, V, N, i1_pros, i2_pros, sqrtS, nlo_string, error_switch)

    ######################
    # prediction
    Mline = np.array(data['Mline'])
    xsec_predic = []
    for i in range(len(mu_list)):
        xsec_predic.append( np.dot(Mline, Coeff[i]) )

    outline = '#---  Test  ---#'
    outline = 'iTest'.rjust(4) + str('prospino').rjust(15) + 'prediction'.rjust(30) + 'ratio'.rjust(30) + 'scale'.rjust(10) + '\n'
    for i in range(len(mu_list)):
        mu = mu_list[i]
        outline += str(itest).rjust(4) + str(xsec_pros[i]).rjust(15) + str(xsec_predic[i]).rjust(30) + str(xsec_predic[i]/xsec_pros[i]).rjust(30) + mu.rjust(10) + '\n'
    print outline
    f.write(outline)
    f.write(run_info)

exit()

##################################################################

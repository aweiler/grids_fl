
#! /usr/bin/env python

import os, sys


iniString1="""export PYTHONPATH=/space/heptools/iminuit/lib/python2.7/site-packages/ 
echo $JOB_ID $HOME $USER $JOB_ID  $JOB_NAME $HOSTNAME $TASK_ID   
scratchDIR="/scratch/${JOB_ID}.aweiler_scratch/"
echo $scratchDIR
mkdir $scratchDIR
cd $scratchDIR
cp -Rv /home/t30/wlr/gi23pah/Projects/EW_Fastlim/grid_generation .
cd grid_generation
"""

mainString = " 8 LO 1 \n"
#mainString = " 8 NLO 1 \n"

finalString="""pwd
ls
cp -v out*.dat /home/t30/wlr/gi23pah/Projects/EW_Fastlim/grid_results/
rm -rf $scratchDIR\n"""

#  1) main_CCsame.py  <==  grid_same.dat
#  2) main_NNsame.py  <==  grid_same.dat
#  3) main_NN.py          <==  grid_XX++.dat
#  4) main_NN.py          <==  grid_XX+-.dat
#  5) main_C2+C1-.py   <==  grid_XX++.dat
#  6) main_C2+C1-.py   <==  grid_XX+-.dat
#  7) main_C2-C1+.py   <==  grid_XX++.dat 
#  8) main_C2-C1+.py   <==  grid_XX+-.dat
#  9) main_NC+.py        <==  grid_CN++.dat
# 10) main_NC+.py        <==  grid_CN+-.dat
# 11) main_NC-.py         <==  grid_CN++.dat
# 12) main_NC-.py         <==  grid_CN+-.dat


fp = open("grid_same.dat")
for i, line in enumerate(fp):
    fname="/home/t30/wlr/gi23pah/Projects/EW_Fastlim/clust_submit/clustrun_CCsame_"+str(i)+"_LO.sh"
    runfile = open(fname,"w")
    runfile.write(iniString1)
    runfile.write("python main.py CCsame "+ line.rstrip()+ mainString)
    runfile.write(finalString)
    runfile.close()
    print "qsub ", fname
fp.close()

# fp = open("grid_same.dat")
# for i, line in enumerate(fp):
#     fname="/home/t30/wlr/gi23pah/Projects/EW_Fastlim/clust_submit/clustrun_NNsame_"+str(i)+".sh"
#     runfile = open(fname,"w")
#     runfile.write(iniString1)
#     runfile.write("python main.py NNsame "+ line.rstrip()+ mainString)
#     runfile.write(finalString)
#     runfile.close()
#     print "qsub ", fname
# fp.close()

# fp = open("grid_XX++.dat")
# for i, line in enumerate(fp):
#     fname="/home/t30/wlr/gi23pah/Projects/EW_Fastlim/clust_submit/clustrun_NN++_"+str(i)+".sh"
#     runfile = open(fname,"w")
#     runfile.write(iniString1)
#     runfile.write("python main.py NN "+ line.rstrip()+ mainString)
#     runfile.write(finalString)
#     runfile.close()
#     print "qsub ", fname
# fp.close()

# fp = open("grid_XX+-.dat")
# for i, line in enumerate(fp):    
#     fname="/home/t30/wlr/gi23pah/Projects/EW_Fastlim/clust_submit/clustrun_NN+-_"+str(i)+".sh"
#     runfile = open(fname,"w")
#     runfile.write(iniString1)
#     runfile.write("python main.py NN "+ line.rstrip()+ mainString)
#     runfile.write(finalString)
#     runfile.close()
#     print "qsub ", fname
# fp.close()

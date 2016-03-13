
#! /usr/bin/env python

import os, sys

#scriptDir="/home/t30/wlr/gi23pah/Projects/EW_Fastlim/clust_submit/"

scriptDir="/space/gi23pah/EW_Fastlim/clust_submit/"


iniString1="""export PYTHONPATH=/space/heptools/iminuit/lib/python2.7/site-packages/ 
SECONDS=0
echo $TMPDIR $JOB_ID $HOME $USER $JOB_ID  $JOB_NAME $HOSTNAME $TASK_ID   
cd $TMPDIR
cp -Rv /home/t30/wlr/gi23pah/Projects/EW_Fastlim/grid_generation .
cd grid_generation
"""

mainString = " 8 NLO 1 \n"
#mainString = " 8 NLO 1 \n"

queue = " -q longrun "

def finalString(fname):
	finalS="""pwd
ls
echo "[LOG] Time since start $SECONDS"
if [ $SECONDS -lt 300 ]
then
        echo "[LOG] Prospino didn't start! resubmitting again"
        ssh zuse "qsub """ + queue + fname + """ "
fi

echo "[LOG] total runtime $SECONDS"
cp -v out*.dat /home/t30/wlr/gi23pah/Projects/EW_Fastlim/grid_results/
\n"""
	return finalS


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
     fname=scriptDir+"clustrun_NNsame_"+str(i)+"_NLO.sh"
     runfile = open(fname,"w")
     runfile.write(iniString1)
     runfile.write("python main.py NNsame "+ line.rstrip()+ mainString)
     runfile.write(finalString(fname))
     runfile.close()
     print "qsub ",queue, fname
fp.close()

fp = open("grid_XX++.dat")
for i, line in enumerate(fp):
    fname=scriptDir+"clustrun_NN++_"+str(i)+"_NLO.sh"
    runfile = open(fname,"w")
    runfile.write(iniString1)
    runfile.write("python main.py NN "+ line.rstrip()+ mainString)
    runfile.write(finalString(fname))
    runfile.close()
    print "qsub ", queue, fname
fp.close()

fp = open("grid_XX+-.dat")
for i, line in enumerate(fp):    
    fname=scriptDir+"clustrun_NN+-_"+str(i)+"_NLO.sh"
    runfile = open(fname,"w")
    runfile.write(iniString1)
    runfile.write("python main.py NN "+ line.rstrip()+ mainString)
    runfile.write(finalString(fname))
    runfile.close()
    print "qsub ", queue , fname
fp.close()

fp = open("grid_same.dat")
for i, line in enumerate(fp):
    fname=scriptDir+"clustrun_CCsame_"+str(i)+"_NLO.sh"
    runfile = open(fname,"w")
    runfile.write(iniString1)
    runfile.write("python main.py CCsame "+ line.rstrip()+ mainString)
    runfile.write(finalString(fname))
    runfile.close()
    print "qsub ", queue, fname
fp.close()



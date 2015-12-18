#!/bin/bash
echo "RUNNING RUNPROSPINO.sh ..................................................... "
wkdir=`pwd`
#prospino_orig=$HOME'/Tools/prospino-2.1'
prospino_orig=/space/gi23pah/prospino_10_17_14
pwd
slha_path=$1
i1=$2
i2=$3
sqrtS=$4
nlo_string=$5
error_switch=$6

#if [[ $i2 == '' ]]; then
    echo '[SLHA file] [i1] [i2] [sqrtS (TeV)] [LO/NLO] [1/0 (with/without scale variation)]'
#    exit
#fi

#if [[ ! -e $slha_path ]]; then
#    echo $slha_path 'does not exist'
#    exit
#fi

##############################

prosDir=$wkdir/working
cp -r $prospino_orig $prosDir

prosDir=$prosDir

cp $slha_path $prosDir'/prospino.in.les_houches'
cd $prosDir

python $wkdir/write_prospino_main.py $i1 $i2 $sqrtS $nlo_string $error_switch > $wkdir/prospino_main
diff $wkdir/prospino_main $prosDir/prospino_main.f90 > /dev/null 2>&1
#if [[ $? -ne 0 ]]; then
    cp $wkdir/prospino_main $prosDir/prospino_main.f90
    make 
#fi

cp prospino.dat prospino.dat_old
# make > /dev/null
./prospino_2.run
cp prospino.dat $wkdir/prospino.out
cat prospino.dat

##############################
echo "DONE ......................................................................"
#rm -rf $prosDir
cd $wkdir
pwd

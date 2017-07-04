#! /usr/bin/env bash

# Defaults

fname="vel.dat"
tstart="11900000"
tend="12000000"
every="1"
bins="400"
unit="metal"
xlow="0"
xhigh="30"
ylow="0"
yhigh="30"
zlow="0"
zhigh="30"
name="None"


while getopts "f:t:T:r:N:u:x:X:y:Y:z:Z:n:" opt; do
    case $opt in
        f)
                fname="$OPTARG"
                ;;
        t)
                tstart="$OPTARG"
                ;;
        T)
                tend="$OPTARG"
                ;;
        r)
                every="$OPTARG"
                ;;
        N)
                bins="$OPTARG"
                ;;
        u)
                unit="$OPTARG"
                ;;
        x)
                xlow="$OPTARG"
                ;;
        X)
                xhigh="$OPTARG"
                ;;
        y)
                ylow="$OPTARG"
                ;;
        Y)
                yhigh="$OPTARG"
                ;;
        z)
                zlow="$OPTARG"
                ;;
        Z)
                zhigh="$OPTARG"
                ;;
        n)
                name="$OPTARG"
                ;;
	esac
done

if [ $unit='metal' ]
then
   dt=0.0005
elif [ $unit='real' ]
then
   dt=2
else
    echo 'Wrong unit'
fi

t_int=$((${tend}-${tstart}))

echo "The parameters used are:"
echo "Filename:" $fname
echo "Start time:" $tstart
echo "End time:" $tend
echo "Time interval:" $t_int
echo "Every:" $every
echo "Number of bins:" $bins
echo "Unit:" $unit
echo "Coordinates:" $xlow $xhigh $ylow $yhigh $zlow $zhigh
echo "The additional name is:" $name

vacf_cmd="dump_read_Cvv_spce.py -f ${fname} -tse ${tstart} ${tend} -re ${every} -dt ${dt} -N ${bins} -coord ${xlow} ${xhigh} ${ylow} ${yhigh} ${zlow} ${zhigh} -name ${name}"
$vacf_cmd
mv C_vv_${every}_${t_int}.dat C_vv_${every}_${t_int}_z${zlow}_${zhigh}.dat
mv C_vv_${every}_${t_int}_norm.dat C_vv_${every}_${t_int}_norm_z${zlow}_${zhigh}.dat
mv C_vv_y_${every}_${t_int}.dat C_vv_y_${every}_${t_int}_z${zlow}_${zhigh}.dat
diff_cmd="diffusion_quick.py -f C_vv_${name}_${every}_${t_int}_z${zlow}_${zhigh}.dat -u ${unit}"
echo $diff_cmd
$diff_cmd

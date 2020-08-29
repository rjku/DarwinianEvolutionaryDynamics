#!/bin/sh

clear

options="--delete -azv --rsh=ssh"
locale="/scratch/$USER/gridWorld"
remoto="riccardorao@mbp142.wireless.ias.edu:/Users/$USER/projects/fisheria/analysis/gridWorld/scratch"


echo "Pulling from <$remoto>."

echo "Pulling simulation outcome:"

rsync --dry-run $options $remoto"/" $locale

echo "Do you want to proceed?[Y/n]"

read -s -n 1 A

if [[ $A = "" ]]; then
	rsync $options $remoto"/" $locale
	echo "Pulling completed"
else
	echo "Pulling aborted"
fi


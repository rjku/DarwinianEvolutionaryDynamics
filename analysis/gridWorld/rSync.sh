#!/bin/sh

clear

HOSTNAME="mbp142.wireless.ias.edu"

options="--delete -azvP --rsh=ssh"
locale="/scratch/$USER/gridWorld"
remoto="$USER@$HOSTNAME:/Users/$USER/projects/fisheria/analysis/gridWorld/scratch"

echo "Syncronizing <$remoto> with <$locale>."
echo "Syncronization simulation outcome:"
rsync --dry-run $options $locale"/" $remoto

echo "Do you want to proceed?[Y/n]"
read -s -n 1 A
if [[ $A = "" ]]; then
	rsync $options $locale"/" $remoto
	echo "Syncronization completed"
else
	echo "Syncronization aborted"
fi


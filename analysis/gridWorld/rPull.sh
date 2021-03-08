#!/bin/sh

clear

HOSTNAME="selene.sns.ias.edu"

options="--delete -azvP --rsh=ssh"
locale="/Users/$USER/projects/fisheria/analysis/gridWorld/scratch"
remoto="$USER@$HOSTNAME:/scratch/$USER/gridWorld"


echo "Pulling from $remoto."
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


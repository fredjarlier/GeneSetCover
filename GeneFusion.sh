#!/bin/bash

SLEEP_MAX=300 #300 <=> 5min
SLEEP_START=1
CHECK_NEED=1
AUTO_DELETE=0
if [ -z "$1" ] #If no arg
then
	CHECK_NEED=1
else	#If at least one arg
	if [ -z "$2" ]	#If only one 
	then
		AUTO_DELETE=0
	else	#If two
		if [[ $2 == "y" ]]; then
			AUTO_DELETE=1
		fi
	fi

	#Ask for the configs ?
	if [[ $1 == "n" ]]; then
		CHECK_NEED=0
	fi
fi

#GeneFusion.sh
	echo "Starting bash GeneFusion"

#########################################################
#Parameters
	echo -n "Loading variables .."
#########################################################
#Custom variables
USR_REP=tmagalha
WR_DIR="/data/tmp/"$USR_REP
PWR_DIR="/mnt/fhgfs/"$USR_REP
DATE=$(date +%d_%m_%Y)
TIME=$(date +%H_%M_%S)
### - Custom variables

echo -n ".."
#Script
INIT_SCRIPT="init_ref_genes_files.sh"
### - Script

	echo -n ".."
#MPI
MPIRUN="/bioinfo/local/build/openmpi-1.8.3/bin/mpirun"
### - MPI

	echo -n ".."
#GeneFusion bin
BIN_NAME="/bioinfo/users/tmagalha/workspace/Release/src/psort"
### - GeneFusion bin

	echo -n ".."
#GeneFusion arguments
	#1 - SAM input
INPUT_FILE=$PWR_DIR"/file_57G.sam"
	#2 - Ouput directory
OUTPUT_DIR=$PWR_DIR"/OUTPUT_MPI/CLUSTERING"
	#3 - Read name
READ_NAME="HWI-ST"
#READ_NAME="HISEQ:"
	#4 - Host name
HOST_NAME="bi-calc"
	#5 - Gene reference file
CHR_REF=$PWR_DIR"/OUTPUT_MPI/CLUSTERING/RefGene/"
#CHR_REF="/bioinfo/users/tmagalha/workspace/Release/script/output/"
	#6 - Time log output
TIME_LOG=$OUTPUT_DIR"/time_"$DATE"_"$TIME".log"
	#7 - Family path
FAMILY_REF=$CHR_REF"family_gene.txt"
	#8 - Exon path
EXON_REF=$CHR_REF
	#9 - Optionnal
SUP_OPT="-p"
### - GeneFusion arguments

	echo -n ".."
#Torque
PPN=20
NODES=2
MEM=100
WALLTIME_H=00
WALLTIME_M=20
WALLTIME_S=00
QUEUE=fjarlier
STDOUT_DIR=$WR_DIR"/PBS_OUTPUT"
STDERR_DIR=$WR_DIR"/PBS_ERROR"
JOB_NAME="GeneFusion_$DATE"
### - Torque

	echo -n ".."
#Cmd
#gdb -ex run --batch --args 
GENE_FUSION_CMD="gdb -ex run --batch --args $BIN_NAME $INPUT_FILE $OUTPUT_DIR $READ_NAME $HOST_NAME $CHR_REF $TIME_LOG $FAMILY_REF $EXON_REF $SUP_OPT"
MPIRUN_CMD="$MPIRUN $GENE_FUSION_CMD"
QSUB_CMD="qsub -o $STDOUT_DIR -e $STDERR_DIR -N ${JOB_NAME} -q $QUEUE -l nodes=$NODES:ppn=$PPN:tengiga,mem=$MEMgb,walltime=$WALLTIME_H:$WALLTIME_M:$WALLTIME_S"
### - Cmd

	echo  " Done"
#########################################################
#Start commands
#########################################################

#########################################################
#Checks
#########################################################

#Ask for deleting existing files in output (Works only if OUTPUT_DIR is accessible)
if [ $AUTO_DELETE == 1 ]; then
	rm -f $OUTPUT_DIR/*
	echo "$OUPUT_DIR blanked"
else
	if [ "$(ls -A $OUTPUT_DIR)" ]; then
		echo "$OUTPUT_DIR is not empty"
		echo "Do you wish to delete existing output (recommended)"
		select result in "Yes|yes|YES|Y|y" "No|no|NO|N|n"
		do
			case $REPLY in
				1|"Y"|"y"|"Yes"|"yes"|"YES") 
					rm -f $OUTPUT_DIR/*
					echo "$OUPUT_DIR blanked"
					break;;
				2|"N"|"n"|"No"|"no"|"NO") 
					echo "$OUPUT_DIR not blanked"
					break;;
				*) echo "Unknown answer";;
			esac 
		done
	else
		echo "$OUTPUT_DIR is empty"
	fi
fi

if [ $CHECK_NEED == 1 ]; then
	#Ask if the chrNames.txt has been edited
	echo "Did you configure the chromosomes's names file ?"
	select result in "Yes|yes|YES|Y|y" "No|no|NO|N|n"
	do
		case $REPLY in
			1|"Y"|"y"|"Yes"|"yes"|"YES") 
				echo "Chromosomes's names configured."
				break;;
			2|"N"|"n"|"No"|"no"|"NO") 
				echo "Do it please.";;
			*) echo "Unknown answer";;
		esac 
	done

	#Ask if the $INIT_SCRIPT has been edited
	echo "Did you configure the $INIT_SCRIPT file ?"
	select result in "Yes|yes|YES|Y|y" "No|no|NO|N|n"
	do
		case $REPLY in
			1|"Y"|"y"|"Yes"|"yes"|"YES") 
				echo "$INIT_SCRIPT configured."
				break;;
			2|"N"|"n"|"No"|"no"|"NO") 
				echo "Do it please.";;
			*) echo "Unknown answer";;
		esac 
	done

	#Check the geneRef
	if [ "$(ls -A $CHR_REF)" ]; then
		echo "$CHR_REF is not empty. Script will consider it has been configured correctly."
	else
		echo "$CHR_REF is empty. Script will generate it."
		bash $INIT_SCRIPT
	fi
fi

#Launch qsub
echo -n "Submitting job to qsub ... "
echo $PBS_JOBID $MPIRUN_CMD $QSUB_CMD
JOB_NAME_ID=`echo "echo $PBS_JOBID; $MPIRUN_CMD" | $QSUB_CMD`
echo " Done"
JOB_ID=${JOB_NAME_ID%%.*}
echo "Created job $JOB_ID"
#End commands


#Wait for job end
QUEUED=1
Q_TIME=0
RUNNING=0
R_TIME=0
SUSPENDED=0
S_TIME=0
echo -n "Queued"
WATCH_CMD=$(qstat -f $JOB_ID | grep job_state | cut -f2 --delimiter='=')
while [[ $WATCH_CMD != *"C"* ]]
do
	echo -n -e "\r\e[K"
	if [[ $WATCH_CMD == *"Q"* ]]; then
		echo -e -n "Queued"
		QUEUED=$(($QUEUED+1))
		Q_TIME=$(($Q_TIME+$SLEEP_START))
	elif [[ $WATCH_CMD == *"R"* ]]; then
		echo -e -n "Running"
		RUNNING=$(($RUNNING+1))
		R_TIME=$(($R_TIME+$SLEEP_START))
	elif [[ $WATCH_CMD == *"S"* ]]; then
		echo -e -n "SUSPENDED !!"
		SUSPENDED=$(($SUSPENDED+1))
		S_TIME=$(($S_TIME+$SLEEP_START))
	fi
	echo -e -n " ($SLEEP_START)"
	sleep $SLEEP_START
	if [ $SLEEP_START -lt $SLEEP_MAX ]; then
		SLEEP_START=$(($SLEEP_START+1))
	fi
	WATCH_CMD=$(qstat -f $JOB_ID | grep job_state | cut -f2 --delimiter='=')
done
echo -n -e "\r\e[K"
echo -e "Complete:"
echo -e "\tQueued: $QUEUED times - $(($Q_TIME/3600)):$(($Q_TIME/60)):$(($Q_TIME%60))"
echo -e "\tSuspended: $SUSPENDED times - $(($S_TIME/3600)):$(($S_TIME/60)):$(($S_TIME%60))"
echo -e "\tRunning: $RUNNING times - $(($R_TIME/3600)):$(($R_TIME/60)):$(($R_TIME%60))"

echo "Script over!"
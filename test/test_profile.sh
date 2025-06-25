FILES=$1
N_VARS=$2
DIR=formulas/cnfs/$FILES
METHOD=$3
if [ -z "$3" ]
  then
    METHOD="highs-ipm"
fi
TYPE=$4
if [ -z "$4" ]
  then
    TYPE="optimization"
fi


for file in $(ls $DIR-clean);
	do 
		echo $file;
		bash scripts/profile.sh main.py $METHOD $DIR-clean/$file $TYPE $experiments/profiles/$file;
done

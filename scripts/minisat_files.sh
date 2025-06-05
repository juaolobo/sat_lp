N_VARS=$1
mkdir solutions/uf$N_VARS
for d in $(ls cnfs/uf$N_VARS-clean);
	do 
		echo $d;
		bash ./sat_lp/scripts/minisat.sh ./cnfs/uf$N_VARS-clean/$d ./solutions/uf$N_VARS/$d.out;
done

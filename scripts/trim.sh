N_VARS=$1
FILE=uuf$N_VARS

mkdir cnfs/uuf$N_VARS-clean;
# mkdir cnfs/uf$N_VARS-clean;
for d in $(ls cnfs/$FILE);
do
	head -n -3 cnfs/$FILE/$d > cnfs/$FILE-clean/$d
done

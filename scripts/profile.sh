INFILE=$1
OUTFILE=$5
if [ -z "$5" ]
  then
    OUTFILE="prof.out"
fi
kernprof -l $INFILE -m $2 -f $3 -t $4
python3 -m line_profiler $INFILE.lprof > /tmp/prof.out
awk '/Total time/ {print} 
	/Line #/ {print} 
	/self.solve_linear()/ {print} 
	/self.solve_boolean()/ {print}
	/resolved/' /tmp/prof.out  > $OUTFILE
rm /tmp/prof.out $INFILE.lprof

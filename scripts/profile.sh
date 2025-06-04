INFILE=$1
OUTFILE=$2
if [ -z "$2" ]
  then
    OUTFILE="prof.out"
fi
kernprof -l $INFILE 
python3 -m line_profiler $INFILE.lprof > /tmp/prof.out
awk '/Total time/ {print} 
	/Line #/ {print} 
	/self.solve_linear()/ {print} 
	/self.solve_boolean()/ {print}
	/resolved, formula/' /tmp/prof.out  > $OUTFILE
rm /tmp/prof.out $INFILE.lprof

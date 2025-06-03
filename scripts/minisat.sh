#!/bin/sh
SOLUTIONS=0
cp $1 /tmp/tmpsat
rm $2 || echo "creating file $2"
echo $1 $2
while :; do

  minisat -verb=0 /tmp/tmpsat /tmp/tmpout > /tmp/tmpmsg 2> /tmp/tmpmsg

  if [ `head -1 /tmp/tmpout` = UNSAT ]; then
    break
  fi
 SOLUTIONS=$((SOLUTIONS + 1))
 tail -1 /tmp/tmpout >> $2
  tail -1 /tmp/tmpout |
    awk '{
      for(i=1;i<NF;++i) { $i = -$i }
      print
    }' >> /tmp/tmpsat

done

echo There are $SOLUTIONS solutions.
rm -f /tmp/tmpsat
rm -f /tmp/tmpout
rm -f /tmp/tmpmsg

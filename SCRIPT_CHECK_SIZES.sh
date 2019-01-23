#!/bin/sh

threadID="$1"			# Thread ID
lFILE="./outputs/tmp/sizes_$threadID"
if test -s $lFILE
   then 
      echo "The file $lFILE exist,"
      echo "       and it contains the following: $(cat $lFILE)"
   else 
      R0_cnts=`grep "#(SMOOTH DIVISORS REACHED)/#(DYNAMIC FACTOR BASE)" outputs/logs/GET_SMOOTH_DIVISORS_USING_THREAD_$threadID.log | wc -l`
      FB_size=`grep -c "^$" outputs/tmp/FB_lc_$threadID`
      R0_size=`grep -c "^$" outputs/tmp/R0_lc_$threadID`

      echo "$R0_cnts $FB_size $R0_size" > $lFILE
      echo "The file $lFILE was successfully created,"
      echo "           and it contains the following: $(cat $lFILE)"
fi


#!/bin/bash

### r.sun mode 2 loop ###
 BEGIN=21
 END=365
 STEP=30
 NUM_CORES=4
 
 for DAY in `seq $BEGIN $STEP $END` ; do
    DAY_STR=`echo $DAY | awk '{printf("%.03d", $1)}'`
    echo "Processing day $DAY_STR at `date` ..."
 
    LINKE="`g.linke_by_day.py $DAY`"
 
    CMD="r.sun -s elevin=MORA_elev_NED_clip@ian day=$DAY lin=$LINKE step=0.5 \
         beam_rad=rad_beam.$DAY_STR diff_rad=rad_diffuse.$DAY_STR \
         refl_rad=rad_reflected.$DAY_STR glob_rad=rad_global.$DAY_STR \
         insol_time=rad_insol_time.$DAY_STR --quiet"
         
        
    # poor man's multi-threading for a multi-core CPU
    MODULUS=`echo "$DAY $STEP $NUM_CORES" | awk '{print $1 % ($2 * $3)}'`
    if [ "$MODULUS" = "$STEP" ] || [ "$DAY" = "$END" ] ; then
       # stall to let the background jobs finish
       $CMD
       sleep 2
       wait
       #while [ `pgrep -c r.sun` -ne 0 ] ; do
       #   sleep 5
       #done
    else
       $CMD &
    fi
 done
 wait   # wait for background jobs to finish to avoid race conditions
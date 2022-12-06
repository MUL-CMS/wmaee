#!/bin/bash
export LC_NUMERIC="en_US.UTF-8" # if not set, commas are written as "," and not as "."

##########################################################################################################

# This small script prepares the scalar-relativistic Green's function EMTO calculation step!
# And, It is adjusted for the mul-hpc!!
# In principle, it creates a bunch of calculations with varying Wigner-Seitz radius.
# If FCD is set to Y in the input.dat, the full charge density is saved.
# This can be then used to calulate the total energy in the last step, called kfcd!

##########################################################################################################

sws_list=`seq 2.47 0.01 2.69` # change here the desired range for the Wigner-Seitz radius

# create directories
for i in $sws_list; do
    mkdir SWS="$i"
done

file="dlm.dat" # change the name of the corresponding input file here

########## Options for the Exc functional ##########
# adjust IEX tag in the input.dat file
# many more options exist, but these two are the most common one!
# IEX...=  7 GGA-PBE
# IEX...=  4 LDA-PW
####################################################

for s in SWS*; do
    echo $s
    cd $s
	cp ../$file .
    mkdir chd pot # necessary for emto-cpa calculations
    sed -i "s/AAAA/${s: -4}/g" $file
    nohup /calc/asakic/software/emto/bin/kgrn < $file & # executes calculation. Path to kgrn executable needs to set appropriate!

    # The next 3 commands are used to interrupt (not kill!) the started calculation and execute the next one!
    PID=$! # Take the process ID of each calculation executed by nohup.
    sleep 2 # Wait for 2 seconds
    kill -s SIGINT $PID # This command doesn't kill the process. It just interrupts it. It is equivivalent to Ctrl+c
    cd $current_path
done

#!/bin/bash

# Variable initializations
benchmarks=bench/dgbsv.dat
target=bin/tp2poisson1D_direct

# Compute the mean and standard deviation for the run with the current number of queens
awk_prog='
    BEGIN {
        x = 0
        y = 0
    }
    {
        x += $1
        y += $1 ^ 2
    }
    END {
        mean = x / NR
        dev = sqrt(y / NR - mean ^ 2)
        var = (dev / mean) * 100
        printf "%lfs\t%lfs\t\t%lf%\n", mean, dev, var
    }
'

if [ ! -f $target ]; then
    printf "\t\033[1;32mCompiling\033[0m %s\n" "$target"
    make $target
    echo ""
fi

#--------------------#
# BENCHMAKRK SECTION #
#--------------------#
echo "Number of points: 1 000 000"
printf "Number of runs:   10\n\n"

printf "Mean time\tStandard deviation\tVariation percentage\n" > $benchmarks

start=$(date +%s.%N)
printf "Running dgbsv benchmark... "
for i in $(seq 0 10); do
	# Execute and parse output with awk to only retain the time measured
	taskset -c 7 bin/tp2poisson1D_direct | awk '/TIME/ {print $4}' >> tmp.dat
done

awk "$awk_prog" tmp.dat >> $benchmarks
rm -f tmp.dat

printf "\033[1;32m done\033[0m\n"
end=$(date +%s.%N)
diff=$(echo "scale=3; $end - $start" | bc -l)
printf "\nBenchmark finished in %.2f seconds\n\n" "$diff"

cat $benchmarks

exit 0

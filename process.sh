if [ -e $(date '+%Y_%m_%d') ]
then
        mkdir $(date '+%Y_%m_%d_2')
        cd $(date '+%Y_%m_%d_2')
else
        mkdir $(date '+%Y_%m_%d')
        cd $(date '+%Y_%m_%d')
fi

echo "Type the pressure at which the data was collected then press [ENTER]"
read pressure

mkdir $pressure
cd $pressure 

echo "Type the lowest anode voltage with which the data was collected then press [ENTER]"
read low_v

echo "Type the highest anode voltage with which the data was collected then press [ENTER]"
read high_v

echo "Type the intervals between the anode voltages then press [ENTER]"
read interval 

diff="$(($high_v - $low_v))"
#div="$(( ($diff/$interval) +1 ))"
#echo "$div"

for (( i = "$(($low_v))"; i < "$(($high_v+$interval))"; $((i+=$interval)) ))
do 
	mkdir $i
        cd $i
        mkdir raw
        mkdir processed
        cd raw
	echo "Type the last 5 digits of run number of interest for $i anode voltage"
	read digit
	scp -r lariat@lariat-daq06.fnal.gov:/daqdata/lariat_r0$(($digit))* .
	ls -d "$PWD"/* > signals.txt
	cd ../../
done  

for (( j = "$(($low_v))"; j < "$(($high_v+$interval))"; $((j+=$interval)) ))
do	
	cd $j
	cd processed
        lar -c multiple_input_output_slice_job.fcl -S ../raw/signals.txt
        ls -d "$PWD"/* > signals.txt
        cd ../../
done


cd ../

for (( c = "$(($low_v))"; c < "$(($high_v+$interval))"; $((c+=$interval)) ))
do
root -b -l <<EOF
.L ../getData.C
getData("$pressure/$c/processed/signals.txt","$pressure/$c/2atm.root");
.q
EOF
done



python outputter.py --lowv_anode "$(($low_v))" --highv_anode "$(($high_v))" --interval "$(($interval))" --pressure $pressure

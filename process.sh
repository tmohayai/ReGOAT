# create a directory with today's date
if [ -e $(date '+%Y_%m_%d') ]
then
        mkdir $(date '+%Y_%m_%d_2')
        cd $(date '+%Y_%m_%d_2')
else
        mkdir $(date '+%Y_%m_%d')
        cd $(date '+%Y_%m_%d')
fi
# make directories named after the anode voltage settings with which the data were collected
for c in {1250..1375..25}
do
        mkdir $c
        cd $c
        mkdir raw
        mkdir signal
        cd ../
done

# pause to scp the data files with various anode voltage settings from the lariat machine into the directories created above
read -p "Press enter to continue"

# continue with the data unpacking
for d in {1250..1375..25}
do
        
        cd $d
        cd raw
        ls -d "$PWD"/* > signals.txt
        cd ../
        cd signal
        lar -c multiple_input_output_slice_job.fcl -S ../raw/signals.txt
        ls -d "$PWD"/* > signals.txt
        cd ../../

done

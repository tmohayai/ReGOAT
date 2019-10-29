# run this in the same directory as the raw GOAT signal
ls -d "$PWD"/* > signals.txt
lar -c multiple_input_output_slice_job.fcl -S signals.txt
mkdir processed_signal
mv lariat_digit_*.root
cd processed_signal
ls -d "$PWD"/* > signals.txt

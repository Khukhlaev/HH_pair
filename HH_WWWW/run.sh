

input_file_lhe="../../../HH_WWWW_samples/WWWW_500k_events.lhe"
input_file_delphes="../../../HH_WWWW_samples/WWWW_500k_events.root"


root -l -b -q "../make_interface.C(\""$input_file_delphes"\")"
python ../turn_off_unused_branches.py
root -l -b -q "analyser_WWWW_huhlaev.C(\""$input_file_delphes"\", \""$input_file_lhe"\", false)"

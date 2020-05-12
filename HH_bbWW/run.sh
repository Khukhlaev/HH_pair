
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.12.04-13971/x86_64-slc6-gcc62-opt/ROOT-env.sh
export PATH=/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.12.04-13971/x86_64-slc6-gcc62-opt/bin:$PATH
source /cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0/x86_64-slc6/setup.sh

wdir=`pwd`
adir="results"
mkdir $adir

run_anal(){
  postfix=$1
  input_file_delphes=$2
  input_file_lhe=$3

  root -l -b -q "../make_interface.C(\""$input_file_delphes"\")"
  python ../turn_off_unused_branches.py
  root -l -b -q "analyser_WWbb.C(\""$input_file_delphes"\", \""$input_file_lhe"\")" > $adir/"log_"$postfix".log"
}

for index in 1 2 3 4 5 6 7 8 9 10 11 12; do
  run_anal "EFT_"$index "../../../test_bbWW/delphes_"$index".root" "../../../HH_lhe_samples/GF_HH_"$index"_cmsgrid_final.lhe"
done

run_anal "SM" "../../../test_bbWW/delphes_1.root" "../../../HH_lhe_samples/GF_HH_1_cmsgrid_final.lhe"


/*
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.12.04-13971/x86_64-slc6-gcc62-opt/ROOT-env.sh
export PATH=/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.12.04-13971/x86_64-slc6-gcc62-opt/bin:$PATH
source /cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0/x86_64-slc6/setup.sh
cd /beegfs/lfi.mipt.su/scratch/MadGraph

source /cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0/x86_64-slc6/setup.sh
source      /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.08.06-b32f2/x86_64-slc6-gcc62-opt/ROOT-env.sh
export PATH=/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.08.06-b32f2/x86_64-slc6-gcc62-opt/bin:$PATH

*/


#include "Delphes.C"
#include "../LHEReader.C"

TRandom rgen;

TLorentzVector make_jet(Delphes* reader, Long64_t index){
    TLorentzVector vec;
    vec.SetPtEtaPhiM(reader->Jet_PT[index], reader->Jet_Eta[index], reader->Jet_Phi[index], reader->Jet_Mass[index] );
    return vec;
}

TLorentzVector make_lepton(Delphes* reader, Long64_t index, bool electron) {
    TLorentzVector vec;
    if (electron)
        vec.SetPtEtaPhiM(
                reader->Electron_PT[index],
                reader->Electron_Eta[index],
                reader->Electron_Phi[index],
                0.000510998928);
    else
        vec.SetPtEtaPhiM(
                reader->Muon_PT[index],
                reader->Muon_Eta[index],
                reader->Muon_Phi[index],
                0.1056583715);
    return vec;
}

struct selected_lepton {
    selected_lepton () {}
    selected_lepton (Delphes *reader, ULong64_t index1, double pt1, Long64_t charge1, bool electron) {
        index = index1;
        pt = pt1;
        charge = charge1;
        is_electron = electron;
        vec = make_lepton(reader, index1, electron);
    }
    ULong64_t index;
    double pt;
    Long64_t charge = 0;
    bool is_electron; // true - electron, false - muon
    TLorentzVector vec;

    selected_lepton& operator= (const selected_lepton& lepton) {
        index = lepton.index;
        pt = lepton.pt;
        charge = lepton.charge;
        is_electron = lepton.is_electron;
        vec = lepton.vec;
        return *this;
    }
};

bool operator< (const selected_lepton& lhs, const selected_lepton& rhs) {
    return lhs.pt < rhs.pt;
}

bool operator!= (const selected_lepton& lhs, const selected_lepton& rhs) {
    return lhs.vec != rhs.vec;
}


void analyser_WWWW_huhlaev (
        string delphes_file = "../../../HH_WWWW_samples/WWWW_10k_events.root",
        string lhe_file = "../../../HH_WWWW_samples/WWWW_10k_events.lhe",
        bool lhe_format = false, string mode = "" ) {

    TRandom rgen;

    TFile *file = TFile::Open(delphes_file.c_str());
    TTree *tree = (TTree *) file->Get("Delphes");
    Delphes *reader = new Delphes(tree);

    LHEReader reader_lhe;
    if (lhe_format) {
        reader_lhe.weight_open_pattern = "wgt id=\"";
        reader_lhe.weight_exit_pattern = "</wgt>";
        reader_lhe.weight_middle_pattern = "\">";
    }
    reader_lhe.Init(lhe_file.c_str());

    TH1D *selections = new TH1D("selections", "selections", 100, 0, 100);
    selections->Fill("Total", 0);
    selections->Fill("Total_X_Weight", 0);
    selections->Fill("Selected", 0);
    selections->Fill("Selected_X_Weight", 0);

    vector<int> pdf_indexes;
    if (lhe_format) for (int i = 9; i < 40; i++) pdf_indexes.push_back(i);
    else for (int i = 45; i < 145; i++) pdf_indexes.push_back(i);


    selections->Fill("Pdf_Up", 0);
    selections->Fill("Pdf_Down", 0);
    selections->Fill("muR_05_muF_05", 0);
    selections->Fill("muR_05_muF_10", 0);
    selections->Fill("muR_10_muF_05", 0);
    selections->Fill("muR_10_muF_20", 0);
    selections->Fill("muR_20_muF_10", 0);
    selections->Fill("muR_20_muF_20", 0);
    selections->Fill("2l channel - number_leptons == 2", 0);
    selections->Fill("2l channel - TMath::Abs(sum_leptons_charge) != 2", 0);
    selections->Fill("2l channel - jets_indexes.size() < 3", 0);
    selections->Fill("2l channel - selected", 0);
    selections->Fill("Selected 3 leptons events", 0);
    selections->Fill("Summary mass of 4 leptons < 180 GeV", 0);
    selections->Fill("Summary mass of 4 leptons > 180 GeV", 0);
    selections->Fill("Selected events", 0);
	
    Long64_t entries = tree->GetEntries();
    double weight_sum = 0;

    for (Long64_t entry = 0; entry < entries; entry++) {
        reader->GetEntry(entry);
        selections->Fill("Total", 1);
//if (entry > 10000) break;
        if (entry % 5000 == 0) {
            cerr << entry << '/' << entries << endl;
        }

        reader_lhe.ReadEvent();
        LHEEvent *lhe_info = &reader_lhe.event;
        double weight = 1;
        if (lhe_format) weight = lhe_info->weights_v[0];
        else weight = lhe_info->weights_v[44];

        weight_sum += weight;

        // OBJECT SELECTIONS ==============================================

        selections->Fill("Total_X_Weight", weight);
        selections->Fill("Total (weighted)", weight);
        reader->GetEntry(entry);

        vector<ULong64_t> muons_indexes, electrons_indexes;
        for (ULong64_t i = 0; i < reader->Muon_; ++i) {
            if (reader->Muon_PT[i] <= 10) continue;
            if (TMath::Abs(reader->Muon_Eta[i]) > 2.5) continue;
            if(reader->Muon_IsolationVar[i] > 0.15) continue;
            muons_indexes.push_back(i);
            selections->Fill("Pre-selected muons", 1);
        }

        for (ULong64_t i = 0; i < reader->Electron_; ++i) {
            double eta = TMath::Abs(reader->Electron_Eta[i]);
            if (reader->Electron_PT[i] <= 10) continue;
            if (eta > 2.47 or (eta > 1.37 and eta < 1.52)) continue;
            if(reader->Electron_IsolationVar[i] > 0.15) continue;
            electrons_indexes.push_back(i);
            selections->Fill("Pre-selected electrons", 1);
        }

        selections->Fill("Minimum lepton p_t > 10 GeV", 1);

        vector <ULong64_t> jets_indexes;
        for (ULong64_t i = 0; i < reader->Jet_; ++i) {
            if (reader->Jet_PT[i] < 25) continue;
            if (TMath::Abs(reader->Jet_Eta[i]) > 2.5) continue;
            jets_indexes.push_back(i);
            selections->Fill("Pre-selected jets", 1);
        }

        // Remove jets which too close to leptons
        vector<ULong64_t> selected_jets;
        for (auto i : jets_indexes) {
            bool pass = true;
            TLorentzVector jet_vector = make_jet(reader, i);
            for (auto j : electrons_indexes) {
                TLorentzVector electron_vector = make_lepton(reader, j, true);
                if (jet_vector.DeltaR(electron_vector) > 0.2) continue;

                break;
            }
	    if(pass) selected_jets.push_back(i);
	}
        jets_indexes = selected_jets;
        selected_jets.clear();
        selections->Fill("Selected jets", jets_indexes.size());

        // Remove electrons which too close to jets
        vector<ULong64_t> selected_electrons;
        for(ULong64_t i : electrons_indexes) {
            TLorentzVector electron_vector = make_lepton(reader, i, true);
            bool pass = true;
            for(auto j : jets_indexes) {
                TLorentzVector jet_vector = make_jet(reader, j);
                if(jet_vector.DeltaR(electron_vector) > 0.4) continue;
                pass = false;
                break;
            }
            if(pass) selected_electrons.push_back(i);
        }
        electrons_indexes = selected_electrons;
        selected_electrons.clear();


        // Remove muons which too close to jets
        vector<ULong64_t> selected_muons;
        for (auto j : muons_indexes) {
            bool pass = true;
            TLorentzVector muon_vector = make_lepton(reader, j, false);
            for(auto i : jets_indexes) {
                TLorentzVector jet_vector = make_jet(reader, i);
                if (jet_vector.DeltaR(muon_vector) > TMath::Min(0.4, 0.04 + 10. / muon_vector.Pt())) continue;
                pass = false;
                break;
            }
            if(pass) selected_muons.push_back(j);
        }
        muons_indexes = selected_muons;
        selected_muons.clear();


        if (electrons_indexes.size() + muons_indexes.size() < 4) {
            for (int i = 0; i < electrons_indexes.size(); i++) {
                if (reader->Electron_IsolationVar[electrons_indexes[i]] > 0.06) {
                    electrons_indexes.erase(electrons_indexes.begin() + i);
                    selections->Fill("Removed electrons cause big isolationVar", 1);
                }
            }
            for (int i = 0; i < muons_indexes.size(); i++) {
                if (reader->Muon_IsolationVar[muons_indexes[i]] > 0.06) {
                    muons_indexes.erase(muons_indexes.begin() + i);
                    selections->Fill("Removed muons cause big isolationVar", 1);
                }
            }
         } 

        selections->Fill("Selected electrons", electrons_indexes.size());
        selections->Fill("Selected muons", muons_indexes.size());

        vector<selected_lepton> leptons;
        for (auto i : electrons_indexes)  {
            double pt = reader->Electron_PT[i];
            Long64_t charge = reader->Electron_Charge[i];
            leptons.push_back(selected_lepton(reader, i, pt, charge, true));
        }
        for (auto i : muons_indexes)  {
            double pt = reader->Muon_PT[i];
            Long64_t charge = reader->Muon_Charge[i];
            leptons.push_back(selected_lepton(reader, i, pt, charge, false));
        }
        sort(leptons.begin(), leptons.end()); // To know leading, subleading lepton, operator< has been overridden

        // EVENTS SELECTION =========================================================


        ULong64_t number_leptons = leptons.size();
        if (number_leptons < 2 || number_leptons > 4) continue;
        selections->Fill("At least 2 electrons or muons", 1);

        bool b_tag = false;
        for (auto i : jets_indexes) {
            if (reader->Jet_BTag[i]) {
                b_tag = true;
                break;
            }
        }
        if (b_tag) continue;
        selections->Fill("Events with no b tags", 1);

        vector<vector<selected_lepton>> sfos; // Same-flavor opposite sign (SFOS)
        Long64_t sum_leptons_charge = 0, min_pt = 1000, max_pt = 0;
        for (ULong64_t i = 0; i < leptons.size(); i++) {
            sum_leptons_charge += leptons[i].charge;
            // Search SFOS
            for (int j = i + 1; j < leptons.size(); j++) {
                if (leptons[i].charge != leptons[j].charge && leptons[i].is_electron == leptons[j].is_electron) {
                    vector<selected_lepton> pair = {leptons[i], leptons[j]};
                    sfos.push_back(pair);
                }
            }
        }
        int number_sfos = sfos.size();

        if (number_leptons == 2) {
            selections->Fill("2l channel - number_leptons == 2", 1);

            if (TMath::Abs(sum_leptons_charge) != 2) continue;
            selections->Fill("2l channel - TMath::Abs(sum_leptons_charge) != 2", 1);

            if (jets_indexes.size() < 3) continue;
            selections->Fill("2l channel - jets_indexes.size() < 3", 1);

            TLorentzVector sum_leptons = leptons[0].vec + leptons[1].vec;
	   
            double nearest_R_sub_leading = 100000000, nearest_R_leading = 100000000;
            ULong64_t index_nearest_jet_leading; // Compared to leading lepton (leading and sub_leading - about leptons)
            for (auto i : jets_indexes) {
                TLorentzVector vec_jet = make_jet(reader, i);
                double R_sub_leading = vec_jet.DeltaR(leptons[0].vec),
                        R_leading = vec_jet.DeltaR(leptons[1].vec);
                if (nearest_R_sub_leading > R_sub_leading)
                    nearest_R_sub_leading = R_sub_leading;
                if (nearest_R_leading > R_leading) {
                    nearest_R_leading = R_leading;
                    index_nearest_jet_leading = i;
                }
            }

            double sub_min_R = 1000000;
            ULong64_t sub_nearest_jet_index;
            for (auto i : jets_indexes) {
                TLorentzVector vec_jet = make_jet(reader, i);
                double R = vec_jet.DeltaR(leptons[1].vec);
                if (i != index_nearest_jet_leading && R < sub_min_R) {
                    sub_min_R = R;
                    sub_nearest_jet_index = i;
                }
            }
            TLorentzVector vec_leading_lepton_jets = leptons[1].vec + make_jet(reader, index_nearest_jet_leading) +
                                                    make_jet(reader, sub_nearest_jet_index);

	    if (leptons[0].pt <= 20 || leptons[1].pt <= 30) continue;
		
            if (leptons[0].is_electron && leptons[1].is_electron) {
                if (TMath::Abs(sum_leptons.M() - 91.1876) <= 10) continue;  // Comparison with M_Z
                if (sum_leptons.M() > 270 || sum_leptons.M() < 55) continue;
                if (nearest_R_leading > 1.15 || nearest_R_leading < 0.2) continue;
                if (nearest_R_sub_leading > 1.4 || nearest_R_sub_leading < 0.2) continue;
                if (vec_leading_lepton_jets.M() < 40 || vec_leading_lepton_jets.M() > 285) continue;
		selections->Fill("2l - 2 electrons", 1);
            }
            else if (!leptons[0].is_electron && !leptons[1].is_electron) {
                if (sum_leptons.M() > 250 || sum_leptons.M() < 60) continue;
                if (nearest_R_leading > 0.75 || nearest_R_leading < 0.2) continue;
                if (nearest_R_sub_leading > 1.05 || nearest_R_sub_leading < 0.2) continue;
                if (vec_leading_lepton_jets.M() < 30 || vec_leading_lepton_jets.M() > 310) continue;
		selections->Fill("2l - 2 muons", 1);
            }
            else { // If there are one electron and one muon
                if (sum_leptons.M() > 250 || sum_leptons.M() < 75) continue;
                if (nearest_R_leading > 0.8 || nearest_R_leading < 0.2) continue;
                if (nearest_R_sub_leading > 1.15 || nearest_R_sub_leading < 0.2) continue;
                if (vec_leading_lepton_jets.M() < 35 || vec_leading_lepton_jets.M() > 350) continue;
		selections->Fill("2l - muon and electoron", 1);
            } 

            selections->Fill("2l channel - selected", 1);
        }
        else if (number_leptons == 3) {
	    selections->Fill("3l events - ", 1);
            if (TMath::Abs(sum_leptons_charge) != 1) continue;
            selections->Fill("3l events, correct summary charge", 1);
            if (jets_indexes.size() < 2) continue;
	    selections->Fill("3l events, number_jets >= 2", 1);

            selected_lepton l1, l2, l3;
            double closest_to_Z_mass = 1000000;
            for (int i = 0; i < 3; i++) {
                for (int j = i + 1; j < 3; j++) {
                    if (leptons[i].is_electron == leptons[j].is_electron) {
                        TLorentzVector sum = leptons[i].vec + leptons[j].vec;
                        if (TMath::Abs(sum.M() - 91.1876) < TMath::Abs(closest_to_Z_mass - 91.1876))
                            closest_to_Z_mass = sum.M();
                    }
                    if (leptons[i].charge == leptons[j].charge) {
                        l2 = leptons[i];
                        l3 = leptons[j];
                    }
                }
            }

            if (TMath::Abs(closest_to_Z_mass - 91.1876) < 10) continue; // Conparison with M_Z

            for (selected_lepton i : leptons) {
                if (i != l2 && i != l3) {
                    l1 = i;
                    break;
                }
            }
	    
            if (l2.vec.DeltaR(l1.vec) > l3.vec.DeltaR(l1.vec)) {
                selected_lepton tmp = l2;
                l2 = l3;
                l3 = tmp;
            }
            if (l2.pt <= 20 || l3.pt <= 20) continue;

            bool too_low_mass = false;
            for (vector<selected_lepton> i : sfos) {
                TLorentzVector sum = i[0].vec + i[1].vec;
                if (sum.M() < 15) {
                    too_low_mass = true;
                    break;
                }
            }
            if (too_low_mass) continue;
	    
	    selections->Fill("3l - Normall SFOSs mass and pt", 1);
            double nearest_R = 100000000, sub_nearest_R = 100000000;
	    ULong64_t index_nearest_jet = 0, index_sub_nearest_jet = 0;
            for (auto i : jets_indexes) {
                TLorentzVector vec_jet = make_jet(reader, i);
                double R = vec_jet.DeltaR(l3.vec);
                if (R < nearest_R) {
                    nearest_R = R;
                    index_nearest_jet = i;
                }
            }
	    for (auto i : jets_indexes) {
		TLorentzVector vec_jet = make_jet(reader, i);
		double R = vec_jet.DeltaR(l3.vec);
		if (R < sub_nearest_R && i != index_nearest_jet) {
		    sub_nearest_R = R;
		    index_sub_nearest_jet = i;
		}
	    }       

            TLorentzVector sum23 = l2.vec + l3.vec;
            TLorentzVector sum3_jet = l3.vec + make_jet(reader, index_nearest_jet);
            TLorentzVector sum3_two_jets = sum3_jet;
	    if (jets_indexes.size() > 1) sum3_two_jets += make_jet(reader, index_sub_nearest_jet);
//cout << "3l channel, SFOS - " << number_sfos << ", sum23 mass = " << sum23.M() << ", DeltaR = " << l2.vec.DeltaR(l3.vec)<< ", sum_jet mass = " << sum3_jet.M() << ", sum_jetss mass = " << sum3_two_jets.M() << endl;
            if (number_sfos == 0) {
                if (l2.vec.DeltaR(l3.vec) < 2.47 || l2.vec.DeltaR(l3.vec) > 5.85) continue;
                if (sum23.M() < 10 || sum23.M() > 70) continue;
                if (sum3_two_jets.M() < 50 || sum3_two_jets.M() > 110) continue;
                if (sum3_jet.M() < 15 || sum3_jet.M() > 50) continue;

            }
            else { // If number_sfos = 1 or 2
                if (l2.vec.DeltaR(l3.vec) < 2.16 || l2.vec.DeltaR(l3.vec) > 3.5) continue;
                if (sum23.M() < 10 || sum23.M() > 70) continue;
                if (sum3_two_jets.M() < 50 || sum3_two_jets.M() > 115) continue;
                if (sum3_jet.M() < 15 || sum3_jet.M() > 45) continue;

            } 
            selections->Fill("Selected 3 leptons events", 1);

        }

        else {  //  FOUR LEPTONS CHANNEL
	selections->Fill("4l events", 1);
            if (TMath::Abs(sum_leptons_charge) != 0) continue;
            selections->Fill("4l, Correct summary charge", 1);
            if (jets_indexes.size() < 2) continue;
            if (leptons[3].pt <= 22) continue;
		selections->Fill("4l, jets >= 2 and highest pt > 22", 1);
            bool too_low_mass = false;
            for (auto i : sfos) {
                TLorentzVector sum = i[0].vec + i[1].vec;
                if (sum.M() <= 4) {
                    too_low_mass = true;
                    break;
                }
            }
            if (too_low_mass) continue;
		selections->Fill("4l - normal SFOSs mass", 1);
            selected_lepton l0, l1, l2, l3;
            if (number_sfos == 0) {

                double closest_to_Z_mass = 1000000;
                for (int i = 0; i < 4; i++) {
                    for (int j = i + 1; j < 4; j++) {
                        TLorentzVector sum = leptons[i].vec + leptons[j].vec;
                        if (leptons[i].is_electron != leptons[j].is_electron && leptons[i].charge != leptons[j].charge
                            && TMath::Abs(sum.M() - 91.1876) < TMath::Abs(closest_to_Z_mass - 91.1876)) {
                            closest_to_Z_mass = sum.M();
                            if (leptons[i].pt > leptons[j].pt) {
                                l2 = leptons[i];
                                l3 = leptons[j];
                            } else {
                                l2 = leptons[j];
                                l3 = leptons[i];
                            }
                        }
                    }
                }

                if (TMath::Abs(closest_to_Z_mass - 91.1876) <= 5) continue;

                for (int i = 0; i < 4; i++) {
                    if (leptons[i] != l2 && leptons[i] != l3 && l0.charge == 0) l0 = leptons[i];
                    else if (leptons[i] != l2 && leptons[i] != l3) l1 = leptons[i];
                }
                if (l1.pt > l0.pt) {
                    auto tmp = l0;
                    l0 = l1;
                    l1 = tmp;
                }
                TLorentzVector sum = l0.vec + l1.vec;
                if (sum.M() <= 10) continue;

                sum += l2.vec + l3.vec;
                if (sum.M() < 180) selections->Fill("Summary mass of 4 leptons < 180 GeV", 1);
                else selections->Fill("Summary mass of 4 leptons > 180 GeV", 1);
            }
            else { // If number_sfos = 1 or 2

                double closest_to_Z_mass = 1000000;
                for (auto i : sfos) {
                    TLorentzVector sum = i[0].vec + i[1].vec;
                    if (TMath::Abs(sum.M() - 91.1876) < TMath::Abs(closest_to_Z_mass - 91.1876)) {
                        closest_to_Z_mass = sum.M();
                        if (i[0].pt > i[1].pt) {
                            l2 = i[0];
                            l3 = i[1];
                        } else {
                            l2 = i[1];
                            l3 = i[0];
                        }
                    }
                }

                for (int i = 0; i < 4; i++) {
                    if (leptons[i] != l2 && leptons[i] != l3 && l0.charge == 0) l0 = leptons[i];
                    else if (leptons[i] != l2 && leptons[i] != l3) l1 = leptons[i];
                }
                if (l1.pt > l0.pt) {
                    auto tmp = l0;
                    l0 = l1;
                    l1 = tmp;
                }

                TLorentzVector sum_01 = l0.vec + l1.vec, sum_23 = l2.vec + l3.vec;
                TLorentzVector sum = sum_01 + sum_23;
                if (sum_01.M() <= 10) continue;

                if (number_sfos == 1) {
                    if (TMath::Abs(closest_to_Z_mass - 91.1876) <= 5) continue;

                    if (sum.M() < 180) selections->Fill("Summary mass of 4 leptons < 180 GeV", 1);
                    else selections->Fill("Summary mass of 4 leptons > 180 GeV", 1);
                }
                else { // If number_sfos = 2

                    if (sum_23.M() >= 70 && sum_23.M() <= 110) continue;

                    if (sum.M() < 180) {
                        if (l2.vec.DeltaPhi(l3.vec) >= 2.6) continue;
                        selections->Fill("Summary mass of 4 leptons < 180 GeV", 1);
                    }
                    else {
                        if (sum_01.M() >= 70 && sum_01.M() <= 110) continue;
                        selections->Fill("Summary mass of 4 leptons > 180 GeV", 1);
                    }
                }
            } 
        } 

        selections->Fill("Selected events", 1);
    }

    for(int i = 1; i < 100; ++i) {
        double passed = selections->GetBinContent(i);
        string label = selections->GetXaxis()->GetBinLabel(i);

        if(label.size() < 1) break;
        cout << "dic[\"" + label + "\"]" << "=[" << passed
            << ", " << passed  / entries << ", " << passed / weight_sum << "]" << endl;
    }
}

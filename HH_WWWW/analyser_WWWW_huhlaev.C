
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

bool operator != (const selected_lepton& lhs, const selected_lepton& rhs) {
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

    TH1D *selections = new TH1D("selections", "selections", 200, 0, 100);
    selections->Fill("Total", 0);
    selections->Fill("Total (weighted)", 0);
    selections->Fill("Selected events", 0);
    selections->Fill("Selected events (weighted)", 0);
    selections->Fill("Selected events Pdf_Up", 0);
    selections->Fill("Selected events Pdf_Down", 0);
    selections->Fill("Selected events muR_05_muF_05", 0);
    selections->Fill("Selected events muR_05_muF_10", 0);
    selections->Fill("Selected events muR_10_muF_05", 0);
    selections->Fill("Selected events muR_10_muF_20", 0);
    selections->Fill("Selected events muR_20_muF_10", 0);
    selections->Fill("Selected events muR_20_muF_20", 0);

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
    selections->Fill("Pre-selected muons", 0);
    selections->Fill("Pre-selected electrons", 0);
    selections->Fill("Minimum lepton p_t > 10 GeV", 0);
    selections->Fill("Pre-selected jets", 0);
    selections->Fill("Selected jets", 0);
    selections->Fill("Removed electrons cause big isolationVar", 0);
    selections->Fill("Removed muons cause big isolationVar", 0);
    selections->Fill("Selected electrons", 0);
    selections->Fill("Selected muons", 0);
    selections->Fill("At least 2 electrons or muons", 0);
    selections->Fill("Events with no b tags", 0);

    selections->Fill("2l channel - number_leptons == 2", 0);
    selections->Fill("2l channel - TMath::Abs(sum_leptons_charge) == 2", 0);
    selections->Fill("2l channel - jets_indexes.size() >= 3", 0);
    selections->Fill("2l channel - MET > 10", 0);
    selections->Fill("2l channel - selected", 0);
    selections->Fill("2l - 2 electrons", 0);
    selections->Fill("2l - 2 electrons (weighted)", 0);
    selections->Fill("2l - 2 electrons Pdf_Up", 0);
    selections->Fill("2l - 2 electrons Pdf_Down", 0);
    selections->Fill("2l - 2 electrons muR_05_muF_05", 0);
    selections->Fill("2l - 2 electrons muR_05_muF_10", 0);
    selections->Fill("2l - 2 electrons muR_10_muF_05", 0);
    selections->Fill("2l - 2 electrons muR_10_muF_20", 0);
    selections->Fill("2l - 2 electrons muR_20_muF_10", 0);
    selections->Fill("2l - 2 electrons muR_20_muF_20", 0);

    selections->Fill("2l - 2 muons", 0);
    selections->Fill("2l - 2 muons (weighted)", 0);
    selections->Fill("2l - 2 muons Pdf_Up", 0);
    selections->Fill("2l - 2 muons Pdf_Down", 0);
    selections->Fill("2l - 2 muons muR_05_muF_05", 0);
    selections->Fill("2l - 2 muons muR_05_muF_10", 0);
    selections->Fill("2l - 2 muons muR_10_muF_05", 0);
    selections->Fill("2l - 2 muons muR_10_muF_20", 0);
    selections->Fill("2l - 2 muons muR_20_muF_10", 0);
    selections->Fill("2l - 2 muons muR_20_muF_20", 0);

    selections->Fill("2l - muon and electoron", 0);
    selections->Fill("2l - muon and electoron (weighted)", 0);
    selections->Fill("2l - muon and electoron Pdf_Up", 0);
    selections->Fill("2l - muon and electoron Pdf_Down", 0);
    selections->Fill("2l - 2 muon and electoron muR_05_muF_05", 0);
    selections->Fill("2l - 2 muon and electoron muR_05_muF_10", 0);
    selections->Fill("2l - 2 muon and electoron muR_10_muF_05", 0);
    selections->Fill("2l - 2 muon and electoron muR_10_muF_20", 0);
    selections->Fill("2l - 2 muon and electoron muR_20_muF_10", 0);
    selections->Fill("2l - 2 muon and electoron muR_20_muF_20", 0);

    selections->Fill("3l events - ", 0);
    selections->Fill("3l events, correct summary charge", 0);
    selections->Fill("3l events, number_jets >= 2", 0);
    selections->Fill("3l channel - MET > 30", 0); // !
    selections->Fill("3l - Normall SFOSs mass and pt", 0);
    selections->Fill("3l - selected", 0);
    selections->Fill("3l - selected (weighted)", 0);
    selections->Fill("3l - selected Pdf_Up", 0);
    selections->Fill("3l - selected Pdf_Down", 0);

    selections->Fill("3l - selected SFOS=0", 0);
    selections->Fill("3l - selected SFOS=0 (weighted)", 0);
    selections->Fill("3l - selected SFOS=0 Pdf_Up", 0);
    selections->Fill("3l - selected SFOS=0 Pdf_Down", 0);
    selections->Fill("3l - selected SFOS=0 muR_05_muF_05", 0);
    selections->Fill("3l - selected SFOS=0 muR_05_muF_10", 0);
    selections->Fill("3l - selected SFOS=0 muR_10_muF_05", 0);
    selections->Fill("3l - selected SFOS=0 muR_10_muF_20", 0);
    selections->Fill("3l - selected SFOS=0 muR_20_muF_10", 0);
    selections->Fill("3l - selected SFOS=0 muR_20_muF_20", 0);

    selections->Fill("3l - selected SFOS=1,2", 0);
    selections->Fill("3l - selected SFOS=1,2 (weighted)", 0);
    selections->Fill("3l - selected SFOS=1,2 Pdf_Up", 0);
    selections->Fill("3l - selected SFOS=1,2 Pdf_Down", 0);
    selections->Fill("3l - selected SFOS=1,2 muR_05_muF_05", 0);
    selections->Fill("3l - selected SFOS=1,2 muR_05_muF_10", 0);
    selections->Fill("3l - selected SFOS=1,2 muR_10_muF_05", 0);
    selections->Fill("3l - selected SFOS=1,2 muR_10_muF_20", 0);
    selections->Fill("3l - selected SFOS=1,2 muR_20_muF_10", 0);
    selections->Fill("3l - selected SFOS=1,2 muR_20_muF_20", 0);

    selections->Fill("4l events", 0);
    selections->Fill("4l, Correct summary charge", 0);
    selections->Fill("4l, jets >= 2 and highest pt > 22", 0);
    selections->Fill("4l - normal SFOSs mass", 0);
    selections->Fill("4l - selected", 0);
    selections->Fill("4l - selected (weighted)", 0);
    selections->Fill("4l - selected Pdf_Up", 0);
    selections->Fill("4l - selected Pdf_Down", 0);

    selections->Fill("4l - sfos == 0, M_4l < 180 GeV", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 GeV (weighted)", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 GeV Pdf_Up", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 GeV Pdf_Down", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 muR_05_muF_05", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 muR_05_muF_10", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 muR_10_muF_05", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 muR_10_muF_20", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 muR_20_muF_10", 0);
    selections->Fill("4l - sfos == 0, M_4l < 180 muR_20_muF_20", 0);

    selections->Fill("4l - sfos == 0, M_4l > 180 GeV", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 GeV (weighted)", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 GeV Pdf_Up", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 GeV Pdf_Down", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 muR_05_muF_05", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 muR_05_muF_10", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 muR_10_muF_05", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 muR_10_muF_20", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 muR_20_muF_10", 0);
    selections->Fill("4l - sfos == 0, M_4l > 180 muR_20_muF_20", 0);

    selections->Fill("4l - sfos == 1, M_4l < 180 GeV", 0);
    selections->Fill("4l - sfos == 1, M_4l > 180 GeV", 0);

    selections->Fill("4l - sfos == 2, M_4l < 180 GeV", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 GeV (weighted)", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 GeV Pdf_Up", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 GeV Pdf_Down", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 muR_05_muF_05", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 muR_05_muF_10", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 muR_10_muF_05", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 muR_10_muF_20", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 muR_20_muF_10", 0);
    selections->Fill("4l - sfos == 2, M_4l < 180 muR_20_muF_20", 0);

    selections->Fill("4l - sfos == 2, M_4l > 180 GeV", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 GeV (weighted)", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 GeV Pdf_Up", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 GeV Pdf_Down", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 muR_05_muF_05", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 muR_05_muF_10", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 muR_10_muF_05", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 muR_10_muF_20", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 muR_20_muF_10", 0);
    selections->Fill("4l - sfos == 2, M_4l > 180 muR_20_muF_20", 0);

    Long64_t entries = tree->GetEntries();
    double weight_sum = 0;

    for (Long64_t entry = 0; entry < entries; entry++) {
        selections->Fill("Total", 1);

        if (entry % 5000 == 0) {
            cerr << entry << '/' << entries << endl;
        }

        reader_lhe.ReadEvent();
        LHEEvent *lhe_info = &reader_lhe.event;
        double weight = 1;
        if (lhe_format) weight = lhe_info->weights_v[0];
        else weight = lhe_info->weights_v[44];

        weight_sum += weight;

        double weight_pdf_up, weight_pdf_down;
        lhe_info->GetPDFErrors(weight_pdf_up, weight_pdf_down, weight, pdf_indexes);

        selections->Fill("Pdf_Up", weight_pdf_up);
        selections->Fill("Pdf_Down", weight_pdf_down);

        double weights_muR_05_muF_05, weights_muR_05_muF_10, weights_muR_10_muF_05, weights_muR_10_muF_20,
                weights_muR_20_muF_10, weights_muR_20_muF_20;
        if(lhe_format) {
            weights_muR_05_muF_05 = lhe_info->weights_v[8] ;
            weights_muR_05_muF_10 = lhe_info->weights_v[6] ;
            weights_muR_10_muF_05 = lhe_info->weights_v[2] ;
            weights_muR_10_muF_20 = lhe_info->weights_v[1] ;
            weights_muR_20_muF_10 = lhe_info->weights_v[3] ;
            weights_muR_20_muF_20 = lhe_info->weights_v[4] ;
        } else {
            weights_muR_05_muF_05 = lhe_info->weights_v[0] ;
            weights_muR_05_muF_10 = lhe_info->weights_v[5] ;
            weights_muR_10_muF_05 = lhe_info->weights_v[15] ;
            weights_muR_10_muF_20 = lhe_info->weights_v[24] ;
            weights_muR_20_muF_10 = lhe_info->weights_v[34] ;
            weights_muR_20_muF_20 = lhe_info->weights_v[39] ;
        }

        selections->Fill( "muR_05_muF_05", weights_muR_05_muF_05 );
        selections->Fill( "muR_05_muF_10", weights_muR_05_muF_10 );
        selections->Fill( "muR_10_muF_05", weights_muR_10_muF_05 );
        selections->Fill( "muR_10_muF_20", weights_muR_10_muF_20 );
        selections->Fill( "muR_20_muF_10", weights_muR_20_muF_10 );
        selections->Fill( "muR_20_muF_20", weights_muR_20_muF_20 );

        // OBJECT SELECTIONS ==============================================

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
            if(reader->Electron_IsolationVar[i] > 0.30) continue;
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
                pass = false;
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
                    selections->Fill("Removed electrons cause big isolationVar", 1);
                    electrons_indexes.erase(electrons_indexes.begin() + i);
                    i --;
                }
            }
            for (int i = 0; i < muons_indexes.size(); i++) {
                if (reader->Muon_IsolationVar[muons_indexes[i]] > 0.06) {
                    selections->Fill("Removed muons cause big isolationVar", 1);
                    muons_indexes.erase(muons_indexes.begin() + i);
                    i --;
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
        Long64_t sum_leptons_charge = 0;
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
            selections->Fill("2l channel - TMath::Abs(sum_leptons_charge) == 2", 1);

            if (jets_indexes.size() < 3) continue;
            selections->Fill("2l channel - jets_indexes.size() >= 3", 1);

            if ( reader->MissingET_MET[0] < 10 ) continue;
            selections->Fill("2l channel - MET > 10", 1);

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
                selections->Fill("2l - 2 electrons (weighted)", weight);
                selections->Fill("2l - 2 electrons Pdf_Up", weight_pdf_up);
                selections->Fill("2l - 2 electrons Pdf_Down", weight_pdf_down);
                selections->Fill("2l - 2 electrons muR_05_muF_05", weights_muR_05_muF_05);
                selections->Fill("2l - 2 electrons muR_05_muF_10", weights_muR_05_muF_10);
                selections->Fill("2l - 2 electrons muR_10_muF_05", weights_muR_10_muF_05);
                selections->Fill("2l - 2 electrons muR_10_muF_20", weights_muR_10_muF_20);
                selections->Fill("2l - 2 electrons muR_20_muF_10", weights_muR_20_muF_10);
                selections->Fill("2l - 2 electrons muR_20_muF_20", weights_muR_20_muF_20);
            }
            else if (!leptons[0].is_electron && !leptons[1].is_electron) {
                if (sum_leptons.M() > 250 || sum_leptons.M() < 60) continue;
                if (nearest_R_leading > 0.75 || nearest_R_leading < 0.2) continue;
                if (nearest_R_sub_leading > 1.05 || nearest_R_sub_leading < 0.2) continue;
                if (vec_leading_lepton_jets.M() < 30 || vec_leading_lepton_jets.M() > 310) continue;
                selections->Fill("2l - 2 muons", 1);
                selections->Fill("2l - 2 muons (weighted)", weight);
                selections->Fill("2l - 2 muons Pdf_Up", weight_pdf_up);
                selections->Fill("2l - 2 muons Pdf_Down", weight_pdf_down);
                selections->Fill("2l - 2 muons muR_05_muF_05", weights_muR_05_muF_05);
                selections->Fill("2l - 2 muons muR_05_muF_10", weights_muR_05_muF_10);
                selections->Fill("2l - 2 muons muR_10_muF_05", weights_muR_10_muF_05);
                selections->Fill("2l - 2 muons muR_10_muF_20", weights_muR_10_muF_20);
                selections->Fill("2l - 2 muons muR_20_muF_10", weights_muR_20_muF_10);
                selections->Fill("2l - 2 muons muR_20_muF_20", weights_muR_20_muF_20);
            }
            else { // If there are one electron and one muon
                if (sum_leptons.M() > 250 || sum_leptons.M() < 75) continue;
                if (nearest_R_leading > 0.8 || nearest_R_leading < 0.2) continue;
                if (nearest_R_sub_leading > 1.15 || nearest_R_sub_leading < 0.2) continue;
                if (vec_leading_lepton_jets.M() < 35 || vec_leading_lepton_jets.M() > 350) continue;
                selections->Fill("2l - muon and electoron", 1);
                selections->Fill("2l - muon and electoron (weighted)", weight);
                selections->Fill("2l - muon and electoron Pdf_Up", weight_pdf_up);
                selections->Fill("2l - muon and electoron Pdf_Down", weight_pdf_down);
                selections->Fill("2l - 2 muon and electoron muR_05_muF_05", weights_muR_05_muF_05);
                selections->Fill("2l - 2 muon and electoron muR_05_muF_10", weights_muR_05_muF_10);
                selections->Fill("2l - 2 muon and electoron muR_10_muF_05", weights_muR_10_muF_05);
                selections->Fill("2l - 2 muon and electoron muR_10_muF_20", weights_muR_10_muF_20);
                selections->Fill("2l - 2 muon and electoron muR_20_muF_10", weights_muR_20_muF_10);
                selections->Fill("2l - 2 muon and electoron muR_20_muF_20", weights_muR_20_muF_20);

            }

            selections->Fill("2l channel - selected", 1);
        }
        else if (number_leptons == 3) {
            selections->Fill("3l events - ", 1);
            if (TMath::Abs(sum_leptons_charge) != 1) continue;
            selections->Fill("3l events, correct summary charge", 1);
            if (jets_indexes.size() < 2) continue;
            selections->Fill("3l events, number_jets >= 2", 1);
            if ( reader->MissingET_MET[0] < 30 ) continue;
            selections->Fill("3l channel - MET > 30", 1);

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
            TLorentzVector sum3_jet      = l3.vec   + make_jet(reader, index_nearest_jet);
            TLorentzVector sum3_two_jets = sum3_jet + make_jet(reader, index_sub_nearest_jet);

            if (number_sfos == 0) {
                if (l2.vec.DeltaR(l3.vec) < 2.47 || l2.vec.DeltaR(l3.vec) > 5.85) continue;
                if (sum23.M() < 10 || sum23.M() > 70) continue;
                if (sum3_two_jets.M() < 50 || sum3_two_jets.M() > 110) continue;
                if (sum3_jet.M() < 15 || sum3_jet.M() > 50) continue;
                selections->Fill("3l - selected SFOS=0", 1);
                selections->Fill("3l - selected SFOS=0 (weighted)", weight);
                selections->Fill("3l - selected SFOS=0 Pdf_Up", weight_pdf_up);
                selections->Fill("3l - selected SFOS=0 Pdf_Down", weight_pdf_down);
                selections->Fill("3l - selected SFOS=0 muR_05_muF_05", weights_muR_05_muF_05);
                selections->Fill("3l - selected SFOS=0 muR_05_muF_10", weights_muR_05_muF_10);
                selections->Fill("3l - selected SFOS=0 muR_10_muF_05", weights_muR_10_muF_05);
                selections->Fill("3l - selected SFOS=0 muR_10_muF_20", weights_muR_10_muF_20);
                selections->Fill("3l - selected SFOS=0 muR_20_muF_10", weights_muR_20_muF_10);
                selections->Fill("3l - selected SFOS=0 muR_20_muF_20", weights_muR_20_muF_20);
            }
            else { // If number_sfos = 1 or 2
                if (l2.vec.DeltaR(l3.vec) < 2.16 || l2.vec.DeltaR(l3.vec) > 3.5) continue;
                if (sum23.M() < 10 || sum23.M() > 70) continue;
                if (sum3_two_jets.M() < 50 || sum3_two_jets.M() > 115) continue;
                if (sum3_jet.M() < 15 || sum3_jet.M() > 45) continue;
                selections->Fill("3l - selected SFOS=1,2", 1);
                selections->Fill("3l - selected SFOS=1,2 (weighted)", weight);
                selections->Fill("3l - selected SFOS=1,2 Pdf_Up", weight_pdf_up);
                selections->Fill("3l - selected SFOS=1,2 Pdf_Down", weight_pdf_down);
                selections->Fill("3l - selected SFOS=1,2 muR_05_muF_05", weights_muR_05_muF_05);
                selections->Fill("3l - selected SFOS=1,2 muR_05_muF_10", weights_muR_05_muF_10);
                selections->Fill("3l - selected SFOS=1,2 muR_10_muF_05", weights_muR_10_muF_05);
                selections->Fill("3l - selected SFOS=1,2 muR_10_muF_20", weights_muR_10_muF_20);
                selections->Fill("3l - selected SFOS=1,2 muR_20_muF_10", weights_muR_20_muF_10);
                selections->Fill("3l - selected SFOS=1,2 muR_20_muF_20", weights_muR_20_muF_20);
            }
            selections->Fill("3l - selected", 1);
            selections->Fill("3l - selected (weighted)", weight);
            selections->Fill("3l - selected Pdf_Up", weight_pdf_up);
            selections->Fill("3l - selected Pdf_Down", weight_pdf_down);
        }

        else {  //  FOUR LEPTONS CHANNEL
            selections->Fill("4l events", 1);
            if (TMath::Abs(sum_leptons_charge) != 0) continue;
            selections->Fill("4l, Correct summary charge", 1);

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

                TLorentzVector sum = l0.vec + l1.vec;
                if (sum.M() <= 10) continue;

                sum += l2.vec + l3.vec;
                if (sum.M() < 180) {
                    selections->Fill("4l - sfos == 0, M_4l < 180 GeV", 1);
                    selections->Fill("4l - sfos == 0, M_4l < 180 GeV (weighted)", weight);
                    selections->Fill("4l - sfos == 0, M_4l < 180 GeV Pdf_Up", weight_pdf_up);
                    selections->Fill("4l - sfos == 0, M_4l < 180 GeV Pdf_Down", weight_pdf_down);
                    selections->Fill("4l - sfos == 0, M_4l < 180 muR_05_muF_05", weights_muR_05_muF_05);
                    selections->Fill("4l - sfos == 0, M_4l < 180 muR_05_muF_10", weights_muR_05_muF_10);
                    selections->Fill("4l - sfos == 0, M_4l < 180 muR_10_muF_05", weights_muR_10_muF_05);
                    selections->Fill("4l - sfos == 0, M_4l < 180 muR_10_muF_20", weights_muR_10_muF_20);
                    selections->Fill("4l - sfos == 0, M_4l < 180 muR_20_muF_10", weights_muR_20_muF_10);
                    selections->Fill("4l - sfos == 0, M_4l < 180 muR_20_muF_20", weights_muR_20_muF_20);
                }
                else {
                    selections->Fill("4l - sfos == 0, M_4l > 180 GeV", 1);
                    selections->Fill("4l - sfos == 0, M_4l > 180 GeV (weighted)", weight);
                    selections->Fill("4l - sfos == 0, M_4l > 180 GeV Pdf_Up", weight_pdf_up);
                    selections->Fill("4l - sfos == 0, M_4l > 180 GeV Pdf_Down", weight_pdf_down);
                    selections->Fill("4l - sfos == 0, M_4l > 180 muR_05_muF_05", weights_muR_05_muF_05);
                    selections->Fill("4l - sfos == 0, M_4l > 180 muR_05_muF_10", weights_muR_05_muF_10);
                    selections->Fill("4l - sfos == 0, M_4l > 180 muR_10_muF_05", weights_muR_10_muF_05);
                    selections->Fill("4l - sfos == 0, M_4l > 180 muR_10_muF_20", weights_muR_10_muF_20);
                    selections->Fill("4l - sfos == 0, M_4l > 180 muR_20_muF_10", weights_muR_20_muF_10);
                    selections->Fill("4l - sfos == 0, M_4l > 180 muR_20_muF_20", weights_muR_20_muF_20);
                }
            }
            else { // If number_sfos > 0

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

                    if (sum.M() < 180) selections->Fill("4l - sfos == 1, M_4l < 180 GeV", 1);
                    else selections->Fill("4l - sfos == 1, M_4l > 180 GeV", 1);

                }
                else if (number_sfos == 2){ // If number_sfos = 2

                    if (sum_23.M() >= 70 && sum_23.M() <= 110) continue;

                    if (sum.M() < 180) {
                        if (TMath::Abs(l2.vec.DeltaPhi(l3.vec)) >= 2.6) continue;
                        selections->Fill("4l - sfos == 2, M_4l < 180 GeV", 1);
                        selections->Fill("4l - sfos == 2, M_4l < 180 GeV (weighted)", weight);
                        selections->Fill("4l - sfos == 2, M_4l < 180 GeV Pdf_Up", weight_pdf_up);
                        selections->Fill("4l - sfos == 2, M_4l < 180 GeV Pdf_Down", weight_pdf_down);
                        selections->Fill("4l - sfos == 2, M_4l < 180 muR_05_muF_05", weights_muR_05_muF_05);
                        selections->Fill("4l - sfos == 2, M_4l < 180 muR_05_muF_10", weights_muR_05_muF_10);
                        selections->Fill("4l - sfos == 2, M_4l < 180 muR_10_muF_05", weights_muR_10_muF_05);
                        selections->Fill("4l - sfos == 2, M_4l < 180 muR_10_muF_20", weights_muR_10_muF_20);
                        selections->Fill("4l - sfos == 2, M_4l < 180 muR_20_muF_10", weights_muR_20_muF_10);
                        selections->Fill("4l - sfos == 2, M_4l < 180 muR_20_muF_20", weights_muR_20_muF_20);
                    }
                    else {
                        if (sum_01.M() >= 70 && sum_01.M() <= 110) continue;
                        selections->Fill("4l - sfos == 2, M_4l > 180 GeV", 1);
                        selections->Fill("4l - sfos == 2, M_4l > 180 GeV (weighted)", weight);
                        selections->Fill("4l - sfos == 2, M_4l > 180 GeV Pdf_Up", weight_pdf_up);
                        selections->Fill("4l - sfos == 2, M_4l > 180 GeV Pdf_Down", weight_pdf_down);
                        selections->Fill("4l - sfos == 2, M_4l > 180 muR_05_muF_05", weights_muR_05_muF_05);
                        selections->Fill("4l - sfos == 2, M_4l > 180 muR_05_muF_10", weights_muR_05_muF_10);
                        selections->Fill("4l - sfos == 2, M_4l > 180 muR_10_muF_05", weights_muR_10_muF_05);
                        selections->Fill("4l - sfos == 2, M_4l > 180 muR_10_muF_20", weights_muR_10_muF_20);
                        selections->Fill("4l - sfos == 2, M_4l > 180 muR_20_muF_10", weights_muR_20_muF_10);
                        selections->Fill("4l - sfos == 2, M_4l > 180 muR_20_muF_20", weights_muR_20_muF_20);
                    }

                }
            }
            selections->Fill("4l - selected", 1);
            selections->Fill("4l - selected (weighted)", weight);
            selections->Fill("4l - selected Pdf_Up", weight_pdf_up);
            selections->Fill("4l - selected Pdf_Down", weight_pdf_down);
        }

        selections->Fill("Selected events", 1);
        selections->Fill("Selected events (weighted)", weight);
        selections->Fill("Selected events Pdf_Up", weight_pdf_up);
        selections->Fill("Selected events Pdf_Down", weight_pdf_down);
        selections->Fill("Selected events muR_05_muF_05", weights_muR_05_muF_05);
        selections->Fill("Selected events muR_05_muF_10", weights_muR_05_muF_10);
        selections->Fill("Selected events muR_10_muF_05", weights_muR_10_muF_05);
        selections->Fill("Selected events muR_10_muF_20", weights_muR_10_muF_20);
        selections->Fill("Selected events muR_20_muF_10", weights_muR_20_muF_10);
        selections->Fill("Selected events muR_20_muF_20", weights_muR_20_muF_20);

    }

    for(int i = 1; i < 200; ++i) {
        double passed = selections->GetBinContent(i);
        string label = selections->GetXaxis()->GetBinLabel(i);

        if(label.size() < 1) break;
        cout << "dic[\"" + label + "\"]" << "=[" << passed
             << ", " << passed  / entries << ", " << passed / weight_sum << "]" << endl;
    }
}


/*
dic["Total"]=[500000, 1, 76.0959]
dic["Total (weighted)"]=[6570.66, 0.0131413, 1]
dic["Selected events"]=[1700, 0.0034, 0.258726]
dic["Selected events (weighted)"]=[22.3431, 4.46862e-05, 0.00340044]
dic["Selected events Pdf_Up"]=[22.6407, 4.52813e-05, 0.00344572]
dic["Selected events Pdf_Down"]=[22.0455, 4.40911e-05, 0.00335515]
dic["Selected events muR_05_muF_05"]=[28.1847, 5.63694e-05, 0.00428948]
dic["Selected events muR_05_muF_10"]=[26.7414, 5.34828e-05, 0.00406982]
dic["Selected events muR_10_muF_05"]=[23.5514, 4.71027e-05, 0.00358432]
dic["Selected events muR_10_muF_20"]=[21.1876, 4.23753e-05, 0.00322458]
dic["Selected events muR_20_muF_10"]=[18.9599, 3.79199e-05, 0.00288555]
dic["Selected events muR_20_muF_20"]=[17.9781, 3.59563e-05, 0.00273612]
dic["Pdf_Up"]=[6672.35, 0.0133447, 1.01548]
dic["Pdf_Down"]=[6468.96, 0.0129379, 0.984523]
dic["muR_05_muF_05"]=[8271.83, 0.0165437, 1.2589]
dic["muR_05_muF_10"]=[7879.6, 0.0157592, 1.19921]
dic["muR_10_muF_05"]=[6898.45, 0.0137969, 1.04989]
dic["muR_10_muF_20"]=[6252.21, 0.0125044, 0.951535]
dic["muR_20_muF_10"]=[5566.64, 0.0111333, 0.847197]
dic["muR_20_muF_20"]=[5296.46, 0.0105929, 0.806078]
dic["Pre-selected muons"]=[166960, 0.33392, 25.4099]
dic["Pre-selected electrons"]=[148692, 0.297384, 22.6297]
dic["Minimum lepton p_t > 10 GeV"]=[500000, 1, 76.0959]
dic["Pre-selected jets"]=[2.0855e+06, 4.17099, 317.395]
dic["Selected jets"]=[2.0855e+06, 4.17099, 317.395]
dic["Removed electrons cause big isolationVar"]=[20173, 0.040346, 3.07016]
dic["Removed muons cause big isolationVar"]=[11111, 0.022222, 1.691]
dic["Selected electrons"]=[127985, 0.25597, 19.4783]
dic["Selected muons"]=[151860, 0.30372, 23.1118]
dic["At least 2 electrons or muons"]=[51788, 0.103576, 7.88171]
dic["Events with no b tags"]=[46682, 0.093364, 7.10461]
dic["2l channel - number_leptons == 2"]=[40732, 0.081464, 6.19907]
dic["2l channel - TMath::Abs(sum_leptons_charge) == 2"]=[12162, 0.024324, 1.85096]
dic["2l channel - jets_indexes.size() >= 3"]=[6653, 0.013306, 1.01253]
dic["2l channel - MET > 10"]=[6528, 0.013056, 0.993508]
dic["2l channel - selected"]=[1385, 0.00277, 0.210786]
dic["2l - 2 electrons"]=[304, 0.000608, 0.0462663]
dic["2l - 2 electrons (weighted)"]=[3.99547, 7.99094e-06, 0.000608078]
dic["2l - 2 electrons Pdf_Up"]=[4.04738, 8.09475e-06, 0.000615977]
dic["2l - 2 electrons Pdf_Down"]=[3.94357, 7.88713e-06, 0.000600178]
dic["2l - 2 electrons muR_05_muF_05"]=[5.0289, 1.00578e-05, 0.000765358]
dic["2l - 2 electrons muR_05_muF_10"]=[4.787, 9.574e-06, 0.000728542]
dic["2l - 2 electrons muR_10_muF_05"]=[4.19775, 8.3955e-06, 0.000638863]
dic["2l - 2 electrons muR_10_muF_20"]=[3.79954, 7.59909e-06, 0.000578259]
dic["2l - 2 electrons muR_20_muF_10"]=[3.38752, 6.77504e-06, 0.000515553]
dic["2l - 2 electrons muR_20_muF_20"]=[3.2212, 6.4424e-06, 0.00049024]
dic["2l - 2 muons"]=[442, 0.000884, 0.0672687]
dic["2l - 2 muons (weighted)"]=[5.80921, 1.16184e-05, 0.000884113]
dic["2l - 2 muons Pdf_Up"]=[5.88915, 1.17783e-05, 0.00089628]
dic["2l - 2 muons Pdf_Down"]=[5.72926, 1.14585e-05, 0.000871946]
dic["2l - 2 muons muR_05_muF_05"]=[7.34034, 1.46807e-05, 0.00111714]
dic["2l - 2 muons muR_05_muF_10"]=[6.94834, 1.38967e-05, 0.00105748]
dic["2l - 2 muons muR_10_muF_05"]=[6.1375, 1.2275e-05, 0.000934077]
dic["2l - 2 muons muR_10_muF_20"]=[5.49762, 1.09952e-05, 0.000836692]
dic["2l - 2 muons muR_20_muF_10"]=[4.93218, 9.86436e-06, 0.000750637]
dic["2l - 2 muons muR_20_muF_20"]=[4.66733, 9.33465e-06, 0.000710329]
dic["2l - muon and electoron"]=[639, 0.001278, 0.0972505]
dic["2l - muon and electoron (weighted)"]=[8.39838, 1.67968e-05, 0.00127816]
dic["2l - muon and electoron Pdf_Up"]=[8.51097, 1.70219e-05, 0.0012953]
dic["2l - muon and electoron Pdf_Down"]=[8.28578, 1.65716e-05, 0.00126103]
dic["2l - 2 muon and electoron muR_05_muF_05"]=[10.6116, 2.12232e-05, 0.00161499]
dic["2l - 2 muon and electoron muR_05_muF_10"]=[10.0428, 2.00857e-05, 0.00152844]
dic["2l - 2 muon and electoron muR_10_muF_05"]=[8.87486, 1.77497e-05, 0.00135068]
dic["2l - 2 muon and electoron muR_10_muF_20"]=[7.94669, 1.58934e-05, 0.00120942]
dic["2l - 2 muon and electoron muR_20_muF_10"]=[7.13189, 1.42638e-05, 0.00108542]
dic["2l - 2 muon and electoron muR_20_muF_20"]=[6.74785, 1.34957e-05, 0.00102697]
dic["3l events - "]=[5511, 0.011022, 0.838729]
dic["3l events, correct summary charge"]=[5507, 0.011014, 0.83812]
dic["3l events, number_jets >= 2"]=[3287, 0.006574, 0.500254]
dic["3l channel - MET > 30"]=[3051, 0.006102, 0.464337]
dic["3l - Normall SFOSs mass and pt"]=[1597, 0.003194, 0.24305]
dic["3l - selected"]=[32, 6.4e-05, 0.00487014]
dic["3l - selected (weighted)"]=[0.420576, 8.41152e-07, 6.40082e-05]
dic["3l - selected Pdf_Up"]=[0.426366, 8.52732e-07, 6.48894e-05]
dic["3l - selected Pdf_Down"]=[0.414786, 8.29572e-07, 6.3127e-05]
dic["3l - selected SFOS=0"]=[10, 2e-05, 0.00152192]
dic["3l - selected SFOS=0 (weighted)"]=[0.13143, 2.6286e-07, 2.00026e-05]
dic["3l - selected SFOS=0 Pdf_Up"]=[0.133192, 2.66384e-07, 2.02708e-05]
dic["3l - selected SFOS=0 Pdf_Down"]=[0.129668, 2.59336e-07, 1.97344e-05]
dic["3l - selected SFOS=0 muR_05_muF_05"]=[0.165466, 3.30932e-07, 2.51826e-05]
dic["3l - selected SFOS=0 muR_05_muF_10"]=[0.157558, 3.15115e-07, 2.3979e-05]
dic["3l - selected SFOS=0 muR_10_muF_05"]=[0.138033, 2.76067e-07, 2.10075e-05]
dic["3l - selected SFOS=0 muR_10_muF_20"]=[0.12501, 2.5002e-07, 1.90255e-05]
dic["3l - selected SFOS=0 muR_20_muF_10"]=[0.111376, 2.22752e-07, 1.69505e-05]
dic["3l - selected SFOS=0 muR_20_muF_20"]=[0.105932, 2.11865e-07, 1.6122e-05]
dic["3l - selected SFOS=1,2"]=[22, 4.4e-05, 0.00334822]
dic["3l - selected SFOS=1,2 (weighted)"]=[0.289146, 5.78292e-07, 4.40056e-05]
dic["3l - selected SFOS=1,2 Pdf_Up"]=[0.293174, 5.86347e-07, 4.46186e-05]
dic["3l - selected SFOS=1,2 Pdf_Down"]=[0.285118, 5.70237e-07, 4.33927e-05]
dic["3l - selected SFOS=1,2 muR_05_muF_05"]=[0.362335, 7.24669e-07, 5.51443e-05]
dic["3l - selected SFOS=1,2 muR_05_muF_10"]=[0.347236, 6.94472e-07, 5.28465e-05]
dic["3l - selected SFOS=1,2 muR_10_muF_05"]=[0.30175, 6.035e-07, 4.59239e-05]
dic["3l - selected SFOS=1,2 muR_10_muF_20"]=[0.276565, 5.53131e-07, 4.2091e-05]
dic["3l - selected SFOS=1,2 muR_20_muF_10"]=[0.244678, 4.89355e-07, 3.72379e-05]
dic["3l - selected SFOS=1,2 muR_20_muF_20"]=[0.234014, 4.68029e-07, 3.5615e-05]
dic["4l events"]=[439, 0.000878, 0.0668122]
dic["4l, Correct summary charge"]=[437, 0.000874, 0.0665078]
dic["4l, jets >= 2 and highest pt > 22"]=[436, 0.000872, 0.0663556]
dic["4l - normal SFOSs mass"]=[430, 0.00086, 0.0654424]
dic["4l - selected"]=[283, 0.000566, 0.0430703]
dic["4l - selected (weighted)"]=[3.71947, 7.43894e-06, 0.000566072]
dic["4l - selected Pdf_Up"]=[3.76679, 7.53359e-06, 0.000573275]
dic["4l - selected Pdf_Down"]=[3.67214, 7.34429e-06, 0.00055887]
dic["4l - sfos == 0, M_4l < 180 GeV"]=[14, 2.8e-05, 0.00213068]
dic["4l - sfos == 0, M_4l < 180 GeV (weighted)"]=[0.184002, 3.68004e-07, 2.80036e-05]
dic["4l - sfos == 0, M_4l < 180 GeV Pdf_Up"]=[0.186216, 3.72431e-07, 2.83405e-05]
dic["4l - sfos == 0, M_4l < 180 GeV Pdf_Down"]=[0.181788, 3.63577e-07, 2.76667e-05]
dic["4l - sfos == 0, M_4l < 180 muR_05_muF_05"]=[0.230502, 4.61005e-07, 3.50806e-05]
dic["4l - sfos == 0, M_4l < 180 muR_05_muF_10"]=[0.221127, 4.42254e-07, 3.36537e-05]
dic["4l - sfos == 0, M_4l < 180 muR_10_muF_05"]=[0.191815, 3.8363e-07, 2.91926e-05]
dic["4l - sfos == 0, M_4l < 180 muR_10_muF_20"]=[0.176151, 3.52301e-07, 2.68087e-05]
dic["4l - sfos == 0, M_4l < 180 muR_20_muF_10"]=[0.155608, 3.11217e-07, 2.36823e-05]
dic["4l - sfos == 0, M_4l < 180 muR_20_muF_20"]=[0.148962, 2.97924e-07, 2.26708e-05]
dic["4l - sfos == 0, M_4l > 180 GeV"]=[23, 4.6e-05, 0.00350041]
dic["4l - sfos == 0, M_4l > 180 GeV (weighted)"]=[0.302289, 6.04578e-07, 4.60059e-05]
dic["4l - sfos == 0, M_4l > 180 GeV Pdf_Up"]=[0.30626, 6.1252e-07, 4.66103e-05]
dic["4l - sfos == 0, M_4l > 180 GeV Pdf_Down"]=[0.298318, 5.96636e-07, 4.54015e-05]
dic["4l - sfos == 0, M_4l > 180 muR_05_muF_05"]=[0.381784, 7.63568e-07, 5.81044e-05]
dic["4l - sfos == 0, M_4l > 180 muR_05_muF_10"]=[0.361491, 7.22982e-07, 5.50159e-05]
dic["4l - sfos == 0, M_4l > 180 muR_10_muF_05"]=[0.31929, 6.38581e-07, 4.85934e-05]
dic["4l - sfos == 0, M_4l > 180 muR_10_muF_20"]=[0.286151, 5.72302e-07, 4.35498e-05]
dic["4l - sfos == 0, M_4l > 180 muR_20_muF_10"]=[0.256697, 5.13394e-07, 3.90672e-05]
dic["4l - sfos == 0, M_4l > 180 muR_20_muF_20"]=[0.242976, 4.85952e-07, 3.69789e-05]
dic["4l - sfos == 1, M_4l < 180 GeV"]=[0, 0, 0]
dic["4l - sfos == 1, M_4l > 180 GeV"]=[0, 0, 0]
dic["4l - sfos == 2, M_4l < 180 GeV"]=[68, 0.000136, 0.010349]
dic["4l - sfos == 2, M_4l < 180 GeV (weighted)"]=[0.893724, 1.78745e-06, 0.000136017]
dic["4l - sfos == 2, M_4l < 180 GeV Pdf_Up"]=[0.904998, 1.81e-06, 0.000137733]
dic["4l - sfos == 2, M_4l < 180 GeV Pdf_Down"]=[0.88245, 1.7649e-06, 0.000134302]
dic["4l - sfos == 2, M_4l < 180 muR_05_muF_05"]=[1.11619, 2.23238e-06, 0.000169875]
dic["4l - sfos == 2, M_4l < 180 muR_05_muF_10"]=[1.07517, 2.15034e-06, 0.000163632]
dic["4l - sfos == 2, M_4l < 180 muR_10_muF_05"]=[0.927873, 1.85575e-06, 0.000141215]
dic["4l - sfos == 2, M_4l < 180 muR_10_muF_20"]=[0.85859, 1.71718e-06, 0.00013067]
dic["4l - sfos == 2, M_4l < 180 muR_20_muF_10"]=[0.755157, 1.51031e-06, 0.000114929]
dic["4l - sfos == 2, M_4l < 180 muR_20_muF_20"]=[0.725442, 1.45088e-06, 0.000110406]
dic["4l - sfos == 2, M_4l > 180 GeV"]=[130, 0.00026, 0.0197849]
dic["4l - sfos == 2, M_4l > 180 GeV (weighted)"]=[1.70859, 3.41718e-06, 0.000260033]
dic["4l - sfos == 2, M_4l > 180 GeV Pdf_Up"]=[1.73065, 3.4613e-06, 0.00026339]
dic["4l - sfos == 2, M_4l > 180 GeV Pdf_Down"]=[1.68653, 3.37306e-06, 0.000256676]
dic["4l - sfos == 2, M_4l > 180 muR_05_muF_05"]=[2.15703, 4.31406e-06, 0.000328282]
dic["4l - sfos == 2, M_4l > 180 muR_05_muF_10"]=[2.04318, 4.08636e-06, 0.000310955]
dic["4l - sfos == 2, M_4l > 180 muR_10_muF_05"]=[1.80402, 3.60805e-06, 0.000274558]
dic["4l - sfos == 2, M_4l > 180 muR_10_muF_20"]=[1.61798, 3.23596e-06, 0.000246243]
dic["4l - sfos == 2, M_4l > 180 muR_20_muF_10"]=[1.45093, 2.90185e-06, 0.000220819]
dic["4l - sfos == 2, M_4l > 180 muR_20_muF_20"]=[1.37386, 2.74771e-06, 0.00020909]
*/
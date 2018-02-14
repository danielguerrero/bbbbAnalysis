/*
** class  : skim_ntuple.cpp
** author : L. Cadamuro (UF)
** date   : 31/12/2017
** brief  : transforms a NanoAOD into a bbbb ntuple for the subsequent plots/analysis
*/

#include <iostream>
#include <string>
#include <iomanip>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "CfgParser.h"
#include "NanoAODTree.h"
#include "EventInfo.h"

#include "SkimUtils.h"
namespace su = SkimUtils;

#include "OfflineProducerHelper.h"
namespace oph = OfflineProducerHelper;

#include "OutputTree.h"
#include "SkimEffCounter.h"

#include "TFile.h"

using namespace std;

// #include <boost/version.hpp>
// std::cout << "Boost version: " 
//       << BOOST_VERSION / 100000
//       << "."
//       << BOOST_VERSION / 100 % 1000
//       << "."
//       << BOOST_VERSION % 100 
//       << std::endl;
// return 0;

// skim_ntuple.exe --cfg config/skim.cfg --input inputFiles/Samples_80X/VBF_HH_4b_10gen2018.txt --output test_bbbb_tree.root --xs 10 --is-signal

int main(int argc, char** argv)
{
    cout << "[INFO] ... starting program" << endl;

    ////////////////////////////////////////////////////////////////////////
    // Decalre command line options
    ////////////////////////////////////////////////////////////////////////
    
    po::options_description desc("Skim options");
    desc.add_options()
        ("help", "produce help message")
        // required
        ("cfg"   , po::value<string>()->required(), "skim config")
        ("input" , po::value<string>()->required(), "input file list")
        ("output", po::value<string>()->required(), "output file LFN")
        // optional
        ("xs"       , po::value<float>(), "cross section [pb]")
        ("maxEvts"  , po::value<int>()->default_value(-1), "max number of events to process")
        // flags
        ("is-data",    po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a data sample (default is false)")
        ("is-signal",  po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a HH signal sample (default is false)")
        ("is-VBF-sig", po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a HH VBF signal sample (default is false)")
        ("save-p4",    po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "save the tlorentzvectors in the output")
    ;

    po::variables_map opts;
    try {
        po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), opts);
        if (opts.count("help")) {
            cout << desc << "\n";
            return 1;
        }
        po::notify(opts);
    }    
    catch (po::error& e) {
        // cerr << "ERROR: " << e.what() << endl << endl << desc << endl;
        cerr << "** [ERROR] " << e.what() << endl;
        return 1;
    }

    ////////////////////////////////////////////////////////////////////////
    // Read config and other cmd line options for skims
    ////////////////////////////////////////////////////////////////////////
    
    CfgParser config;
    if (!config.init(opts["cfg"].as<string>())) return 1;
    cout << "[INFO] ... using config file " << opts["cfg"].as<string>() << endl;

    const bool is_data = opts["is-data"].as<bool>();
    if (!is_data && !opts.count("xs")){
        cerr << "** [ERROR] please provide xs for this sample that is marked as MC" << endl;
        return 1;
    }
    cout << "[INFO] ... is a data sample? " << std::boolalpha << is_data << std::noboolalpha << endl;

    const bool is_signal = (is_data ? false : opts["is-signal"].as<bool>());
    cout << "[INFO] ... is a HH signal sample? " << std::boolalpha << is_signal << std::noboolalpha << endl;

    const bool is_VBF_sig = (is_data ? false : opts["is-VBF-sig"].as<bool>());
    cout << "[INFO] ... is a HH VBF sample? " << std::boolalpha << is_VBF_sig << std::noboolalpha << endl;

    const float xs = (is_data ? 1.0 : opts["xs"].as<float>());
    cout << "[INFO] ... cross section is : " << xs << " pb" << endl;

    const string bbbbChoice = config.readStringOpt("parameters::bbbbChoice");
    const oph::bbbbSelectionStrategy b_sel_strat = (
        bbbbChoice == "OneClosestToMh"  ? oph::bbbbSelectionStrategy::kOneClosestToMh  :
        bbbbChoice == "BothClosestToMh" ? oph::bbbbSelectionStrategy::kBothClosestToMh :
        bbbbChoice == "MostBackToBack"  ? oph::bbbbSelectionStrategy::kMostBackToBack  :
        throw std::runtime_error("cannot recognize bbbb pair choice strategy " + bbbbChoice)
    );
    cout << "[INFO] ... chosing bb bb jet pairs with strategy : " << bbbbChoice << endl;
    

    ////////////////////////////////////////////////////////////////////////
    // Prepare event loop
    ////////////////////////////////////////////////////////////////////////

    cout << "[INFO] ... opening file list : " << opts["input"].as<string>().c_str() << endl;
    if ( access( opts["input"].as<string>().c_str(), F_OK ) == -1 ){
        cerr << "** [ERROR] The input file list does not exist, aborting" << endl;
        return 1;        
    }

    TChain ch("Events");
    int nfiles = su::appendFromFileList(&ch, opts["input"].as<string>());
    
    if (nfiles == 0){
        cerr << "** [ERROR] The input file list contains no files, aborting" << endl;
        return 1;
    }
    cout << "[INFO] ... file list contains " << nfiles << " files" << endl;

    cout << "[INFO] ... creating tree reader" << endl;
    NanoAODTree nat (&ch, is_data);

    cout << "[INFO] ... loading the following triggers" << endl;
    for (auto trg : config.readStringListOpt("triggers::makeORof"))
        cout << "   - " << trg << endl;
    nat.triggerReader().setTriggers(config.readStringListOpt("triggers::makeORof"));

    ////////////////////////////////////////////////////////////////////////
    // Prepare the output
    ////////////////////////////////////////////////////////////////////////
    
    string outputFileName = opts["output"].as<string>();
    cout << "[INFO] ... saving output to file : " << outputFileName << endl;
    TFile outputFile(outputFileName.c_str(), "recreate");
    OutputTree ot(
        opts["save-p4"].as<bool>()
    );

    SkimEffCounter ec;

    // /// to declare user branches for studies
    // ot.declareUserIntBranch("integer1", 0);
    // ot.declareUserIntBranch("integer2", 1);
    // ot.declareUserFloatBranch("float1", 0.01);
    // ot.declareUserFloatBranch("float2", 1.001);

    ////////////////////////////////////////////////////////////////////////
    // Execute event loop
    ////////////////////////////////////////////////////////////////////////

    int maxEvts = opts["maxEvts"].as<int>();
    if (maxEvts >= 0)
        cout << "[INFO] ... running on : " << maxEvts << " events" << endl;

    for (int iEv = 0; true; ++iEv)
    {
        if (maxEvts >= 0 && iEv >= maxEvts)
            break;

        if (!nat.Next()) break;
        if (iEv % 10000 == 0) cout << "... processing event " << iEv << endl;

        ot.clear();
        EventInfo ei;

        // --- - --- - --- --- - --- - --- --- - --- - --- 
        // FIXME: compute the weights
        double w_fixme = 0.75;
        ec.updateProcessed(w_fixme);

        // --- - --- - --- --- - --- - --- --- - --- - --- 

        // cout << "---- my jets" << endl;
        // for (uint ij = 0; ij < *(nat.nJet); ++ij)
        // {
        //     Jet jet (ij, &nat);
        //     // cout << jet.P4().Px() << " " << jet.P4().Py() << " " << jet.P4().Pz() << " " << jet.P4().E() << " btag " << nat.Jet_btagCSVV2.At(ij) << endl;
        //     cout << " Pt  : "     << setw(10) << jet.P4().Pt()
        //          << " Eta : "     << setw(10) << jet.P4().Eta()
        //          << " Phi : "     << setw(10) << jet.P4().Phi()
        //          << " E   : "     << setw(10) << jet.P4().E()
        //          << " | CSV   : " << nat.Jet_btagCSVV2.At(ij)
        //          << endl;
        // }

        // if (!nat.triggerReader().getTrgOr()) continue;

        // if (!oph::select_bbbb_jets(nat, ei, oph::bbbbSelectionStrategy::kMostBackToBack)) continue;
        if (!oph::select_bbbb_jets(nat, ei, b_sel_strat)) continue;

        if (is_signal){
            oph::select_gen_HH(nat, ei);
            if (!oph::select_gen_bb_bb(nat, ei))
                return 1;
            
        }

        if (is_VBF_sig){
            bool got_gen_VBF = oph::select_gen_VBF_partons(nat, ei);
            if (!got_gen_VBF){
                cout << "Failed on iEv = " << iEv << " evt num = " << *(nat.event) << " run = " << *(nat.run) << endl;
                return 1;
            }
        }

        // --- - --- - --- --- - --- - --- --- - --- - --- 

        // --- use the user Floats
        // if (iEv == 5)
        //     ot.userInt("integer1") = 5;
        // if (iEv == 6)
        //     ot.userInt("integer2") = 66;
        // if (iEv == 10){
        //     ot.userFloat("float1") = 5.1;
        //     ot.userFloat("float2") = 5.2;
        //     ot.userInt("integer2") = 52;
        // }
        // wrong names cause an exception
        // ot.userInt("integer22") = 52;

        // cout << "EVENT " << *(nat.event) << endl;
        // cout << "Selected H1_b1 : " << ei.H1_b1->getIdx() << endl;
        // cout << "Selected H1_b2 : " << ei.H1_b2->getIdx() << endl;
        // cout << "Selected H2_b1 : " << ei.H2_b1->getIdx() << endl;
        // cout << "Selected H2_b2 : " << ei.H2_b2->getIdx() << endl;
        // cout << "And linked jets " << endl;
        // for (uint ij = 0; ij < *(nat.nJet); ++ij)
        // {
        //     cout << "** " <<  ij << " mothIdx " << endl; 
        // }

        // --- - --- - --- --- - --- - --- --- - --- - --- 
        // auto jets_all  = oph::make_jets(nat);
        // auto jets_pt50 = oph::make_jets(nat);
        // oph::filter_jets(jets_pt50,
        //     std::function<bool (Jet)> ([](Jet jet) -> bool {return jet.P4().Pt() > 50;})
        // );

        // auto jets_pt50 = oph::make_jets (
        //     nat,
        //     std::function<bool (Jet)> ([](Jet jet) -> bool {return jet.P4().Pt() > 50;})
        // );
        // auto jets_pt50_etapos = oph::make_jets(
        //     nat,
        //     std::function<bool (Jet)> ([](Jet jet) -> bool {return (jet.P4().Pt() > 50 && jet.P4().Eta() > 0);})
        // );

        // cout << "ALL - iEv " << iEv << endl;
        // for (auto& j : jets_all)
        //     cout << j.P4().Pt() << " " << j.P4().Eta() << endl;

        // cout << "PT 50 - iEv " << iEv << endl;
        // for (auto& j : jets_pt50)
        //     cout << j.P4().Pt() << " " << j.P4().Eta() << endl;
        
        // cout << "PT 50 - eta pos iEv " << iEv << endl;
        // for (auto& j : jets_pt50_etapos)
        //     cout << j.P4().Pt() << " " << j.P4().Eta() << endl;

        // --- - --- - --- --- - --- - --- --- - --- - --- 

        ec.updateSelected(w_fixme);
        su::fill_output_tree(ot, nat, ei);

    }

    outputFile.cd();
    ot.write();
    ec.write();
}
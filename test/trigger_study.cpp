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
#include "TMath.h"
#include "TFile.h"

using namespace std;
// trigger_study.exe --cfg config/skim.cfg --input inputFiles/Samples_80X/VBF_HH_4b_10gen2018_test.txt --maxEvts 1000

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
        // optional
        ("maxEvts"  , po::value<int>()->default_value(-1), "max number of events to process")
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
    // Read config and other cmd line options for triggers
    ////////////////////////////////////////////////////////////////////////
    
    CfgParser config;
    if (!config.init(opts["cfg"].as<string>())) return 1;
    cout << "[INFO] ... using config file " << opts["cfg"].as<string>() << endl;
    
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
    NanoAODTree nat (&ch);

    cout << "[INFO] ... reading all HLT paths" << endl;
    vector <string> trigger_paths;
    TObjArray* List = ch.GetListOfLeaves();
    for(int  j=0;  j < List->GetEntries();  j++){
              string name = List->At(j)->GetName();
              if(name.substr(0,4)=="HLT_"){
              trigger_paths.push_back(name);
              }
    }    

    cout << "[INFO] ... getting the HLT efficiencies" << endl;
    vector <int> ntrig(trigger_paths.size(),0);    
    int ntot=0;
    int maxEvts = opts["maxEvts"].as<int>();
    nat.triggerReader().setTriggers(trigger_paths);  
  
    ////////////////////////////////////////////////////////////////////////
    // Execute Event loop
    ////////////////////////////////////////////////////////////////////////        
    for (int iEv = 0; true; ++iEv)
        {
        if (maxEvts >= 0 && iEv >= maxEvts)break;
        if (!nat.Next()) break;
        if (iEv % 10000 == 0) cout << "... processing event " << iEv << endl;
        ////////////////////////////////////////////////////////////////////////
        // Execute HLT paths loop
        ////////////////////////////////////////////////////////////////////////  
        for (unsigned k=0; k<trigger_paths.size(); k++){        
          if (nat.triggerReader().getTrgResult(trigger_paths[k])) ntrig[k]=ntrig[k]+1;
        }
        ++ntot;
    }

    cout << "[INFO] ... saving the HLT efficiencies" << endl;
    ofstream outputFile;
    outputFile.open("efficiencies.txt");
    for (unsigned p=0; p<trigger_paths.size(); p++){        
           outputFile << setw(100)<< trigger_paths[p] << setw(20)<<(double)ntrig[p]/ntot<<endl;
    }
    outputFile.close();




/*
    cout << "[INFO] ... loading the best triggers" << endl;
    vector <string> trigger_path;
    for (auto trg : config.readStringListOpt("triggers::makeORof")){
        cout << "   - " << trg << endl;
        trigger_path.push_back(trg); 
    }
    nat.triggerReader().setTriggers(config.readStringListOpt("triggers::makeORof"));


    ////////////////////////////////////////////////////////////////////////
    // Execute event loop
    ////////////////////////////////////////////////////////////////////////
    int ntot = 0;
    int trig1 = 0;   
    int trig2 = 0;
    int trigOR = 0;   
    int maxEvts = opts["maxEvts"].as<int>();
    if (maxEvts >= 0)
        cout << "[INFO] ... running on : " << maxEvts << " events" << endl;

    for (int iEv = 0; true; ++iEv)
    {
        if (maxEvts >= 0 && iEv >= maxEvts)
            break;

        if (!nat.Next()) break;
        if (iEv % 10000 == 0) cout << "... processing event " << iEv << endl;

        ++ntot;

        ////////////////////////////////////////////////////////////////////////
        // Check individual and OR trigger status
        ////////////////////////////////////////////////////////////////////////
        if (nat.triggerReader().getTrgResult(trigger_path[0])) ++trig1;
        if (nat.triggerReader().getTrgResult(trigger_path[1])) ++trig2;  
        if (nat.triggerReader().getTrgOr()) ++trigOR;   

    }
 
    cout << "[INFO] Trigger report based on " << maxEvts <<" events: "<<endl;
    cout << "|Path       |"<< setw(30) << trigger_path[0] << setw(30) << trigger_path[1] << setw(30) <<"OR"<< endl;
    cout << "|Events     |"<< setw(30) << trig1 << setw(30) << trig2 << setw(30) << trigOR << endl;
    cout << "|Efficiency |"<< setw(30) << (double)trig1 / ntot << setw(30) << (double)trig2 / ntot << setw(30) << (double)trigOR / ntot << endl;


*/
}


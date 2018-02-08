#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <algorithm>
#include <cassert>
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
// USAGE:
// ./bin/trigger_study.exe --cfg config/skim.cfg --input inputFiles/Samples_80X/VBF_HH_4b_10gen2018_test.txt --maxEvts 1000


template <typename T>
void combination(const std::vector<T>& v, std::size_t count)
{
    assert(count <= v.size());
    std::vector<bool> bitset(v.size() - count, 0);
    bitset.resize(v.size(), 1);
    ofstream outputFile;
    string idx=to_string(count);
    outputFile.open(Form("efficiencies_%s.txt",idx.c_str()));  
    do {
        for (std::size_t i = 0; i != v.size(); ++i) {
            if (bitset[i]) {
                outputFile << v[i] << " ";
            }
        }   
        outputFile << std::endl;

    } while (std::next_permutation(bitset.begin(), bitset.end()));
    outputFile.close();

}

int main(int argc, char** argv)
{
    cout << "[INFO] ... starting program" << endl;

    ////////////////////////////////////////////////////////////////////////
    // Declare command line options
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

    cout << "[INFO] ... reading unprescaled HLT paths" << endl;
    vector <string> trigger_paths;
    vector <double> trigger_effs;
    double eff; string name;
    ifstream myfile("test/efficiencies_unprescaled.txt", ifstream::in);
    while(myfile >> name >> eff){ 
       trigger_paths.push_back(name);
       trigger_effs.push_back(eff);
    } 

    cout << "[INFO] ... creating double-path combinations for OR calculation" << endl;
    combination(trigger_paths, 2);
    vector <string> trigger_paths_two_1, trigger_paths_two_2;

    string name1,name2;
    ifstream myfile_2("efficiencies_2.txt", ifstream::in);
    while(myfile_2 >> name1 >> name2){ 
       trigger_paths_two_1.push_back(name1); 
       trigger_paths_two_2.push_back(name2);
    } 

    cout << "[INFO] ... creating triple-path combinations for OR calculation" << endl;
    combination(trigger_paths, 3);
    vector <string> trigger_paths_three_1, trigger_paths_three_2, trigger_paths_three_3;

    string name_three_1,name_three_2,name_three_3;
    ifstream myfile_3("efficiencies_3.txt", ifstream::in);
    while(myfile_3 >> name_three_1 >> name_three_2 >> name_three_3){ 
       trigger_paths_three_1.push_back(name_three_1); 
       trigger_paths_three_2.push_back(name_three_2);
       trigger_paths_three_3.push_back(name_three_3);       
    } 

    cout << "[INFO] ... creating quad-path combinations for OR calculation" << endl;
    combination(trigger_paths, 4);
    vector <string> trigger_paths_four_1, trigger_paths_four_2, trigger_paths_four_3, trigger_paths_four_4;

    string name_four_1,name_four_2,name_four_3,name_four_4;
    ifstream myfile_4("efficiencies_4.txt", ifstream::in);
    while(myfile_4 >> name_four_1 >> name_four_2 >> name_four_3 >> name_four_4){ 
       trigger_paths_four_1.push_back(name_four_1); 
       trigger_paths_four_2.push_back(name_four_2);
       trigger_paths_four_3.push_back(name_four_3);
       trigger_paths_four_4.push_back(name_four_4);              
    }  

    cout << "[INFO] ... getting the OR efficiencies" << endl;
    vector <int> ntrig_two(trigger_paths_two_1.size(),0); 
    vector <int> ntrig_three(trigger_paths_three_1.size(),0);   
    vector <int> ntrig_four(trigger_paths_four_1.size(),0);        
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

        // Double OR
        for (unsigned k=0; k<trigger_paths_two_1.size(); k++){        
          if ( nat.triggerReader().getTrgResult(trigger_paths_two_1[k]) 
            || nat.triggerReader().getTrgResult(trigger_paths_two_2[k]) )  {ntrig_two[k]=ntrig_two[k]+1;}
        }

        // Triple OR
        for (unsigned k=0; k<trigger_paths_three_1.size(); k++){        
          if ( nat.triggerReader().getTrgResult(trigger_paths_three_1[k]) 
            || nat.triggerReader().getTrgResult(trigger_paths_three_2[k]) 
            || nat.triggerReader().getTrgResult(trigger_paths_three_3[k]) )  {ntrig_three[k]=ntrig_three[k]+1;}
        }

        // Quadruple OR
        for (unsigned k=0; k<trigger_paths_four_1.size(); k++){        
          if ( nat.triggerReader().getTrgResult(trigger_paths_four_1[k]) 
            || nat.triggerReader().getTrgResult(trigger_paths_four_2[k]) 
            || nat.triggerReader().getTrgResult(trigger_paths_four_3[k])
            || nat.triggerReader().getTrgResult(trigger_paths_four_4[k]) )  { ntrig_four[k]=ntrig_four[k]+1;}
        }        
        ++ntot;
    }

    cout << "[INFO] ... making Trigger Optimization Report (efficiencies_report.txt)" << endl;
    ofstream outputReport;
    outputReport.open("efficiencies_report.txt");
    outputReport << "" << endl;outputReport << setw(100)<<"*** Optimal OR Trigger Efficiencies ***" << endl;outputReport << "" << endl;    
    outputReport << "" << endl;outputReport << "[RESULTS] Optimal OR Trigger Efficiencies using two HLT paths (eff>0.18)" << endl;outputReport << "" << endl;
    for (unsigned p=0; p<trigger_paths_two_1.size(); p++){ 
        if ( (double)ntrig_two[p]/ntot > 0.18){       
           outputReport << setw(50)<< trigger_paths_two_1[p] << setw(50)<<  trigger_paths_two_2[p] << setw(20)<< (double)ntrig_two[p]/ntot<<endl;
        }
    }
    outputReport << "" << endl;outputReport << "[RESULTS] Optimal OR Trigger Efficiencies using three HLT paths (eff>0.21)" << endl;outputReport << "" << endl;
    for (unsigned p=0; p<trigger_paths_three_1.size(); p++){ 
        if ( (double)ntrig_three[p]/ntot >= 0.21){       
           outputReport << setw(60)<< trigger_paths_three_1[p] << setw(60)<<  trigger_paths_three_2[p] << setw(70)<<  trigger_paths_three_3[p] << setw(20)<< (double)ntrig_three[p]/ntot<<endl;
        }
    }
    outputReport << "" << endl;outputReport << "[RESULTS] Optimal OR Trigger Efficiencies using four HLT paths (eff>0.24)" << endl;outputReport << "" << endl;
    for (unsigned p=0; p<trigger_paths_four_1.size(); p++){ 
        if ( (double)ntrig_four[p]/ntot >= 0.24){       
           outputReport << setw(60)<< trigger_paths_four_1[p] << setw(60)<<  trigger_paths_four_2[p] << setw(70)<<  trigger_paths_four_3[p] <<setw(70)<<  trigger_paths_four_4[p] << setw(20)<< (double)ntrig_four[p]/ntot<<endl;
        }
    }    
    outputReport.close();

}


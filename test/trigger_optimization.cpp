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
#include <tuple>
#include "SkimUtils.h"
namespace su = SkimUtils;
#include <iterator>
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
    vector <string> trigger_paths_allHT;
    vector <double> trigger_effs;  
    vector <double> trigger_effs_allHT;
    double eff; string name;
    ifstream myfile("efficiencies_unprescaled.txt", ifstream::in);
    while(myfile >> name >> eff){ 
         trigger_paths_allHT.push_back(name);
         trigger_effs_allHT.push_back(eff);
         if(name.substr(0,6) != "HLT_HT"){
         trigger_paths.push_back(name);
         trigger_effs.push_back(eff);
         cout << name << endl;
         }        
    } 
    trigger_paths.push_back("HLT_HT430_650");

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
    int nHT=0;        
    int ntot=0;
    int maxEvts = opts["maxEvts"].as<int>();
    nat.triggerReader().setTriggers(trigger_paths_allHT); 
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
        ++ntot;
        // HLT_HT OR
        if( nat.triggerReader().getTrgResult("HLT_HT430to450") || nat.triggerReader().getTrgResult("HLT_HT450to470") ||
            nat.triggerReader().getTrgResult("HLT_HT470to500") || nat.triggerReader().getTrgResult("HLT_HT500to550") ||
            nat.triggerReader().getTrgResult("HLT_HT550to650") || nat.triggerReader().getTrgResult("HLT_HT650") ){nHT++;}

        // Double OR
        for (unsigned k=0; k<trigger_paths_two_1.size(); k++)
        {     
          //Case 1 and case 2  
          if (trigger_paths_two_1[k]=="HLT_HT430_650" || trigger_paths_two_2[k]=="HLT_HT430_650")
          {
                //Case 1
                if( trigger_paths_two_1[k]=="HLT_HT430_650" ){
                     if( nat.triggerReader().getTrgResult(trigger_paths_two_2[k]) || 
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")) {ntrig_two[k]=ntrig_two[k]+1;}
                }
                if( trigger_paths_two_2[k]=="HLT_HT430_650" ){
                //Case2    
                     if( nat.triggerReader().getTrgResult(trigger_paths_two_1[k]) || 
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_two[k]=ntrig_two[k]+1;}
                }
          }
          //case 3
          else{  if( nat.triggerReader().getTrgResult(trigger_paths_two_1[k]) || 
                     nat.triggerReader().getTrgResult(trigger_paths_two_2[k]) )  {ntrig_two[k]=ntrig_two[k]+1;}
          }    
        }

        // Triple OR
        for (unsigned k=0; k<trigger_paths_three_1.size(); k++){ 
          //Case 1,2 and 3
          if(trigger_paths_three_1[k]=="HLT_HT430_650" || trigger_paths_three_2[k]=="HLT_HT430_650" || trigger_paths_three_3[k]=="HLT_HT430_650" )
          {
                //Case 1
                if(trigger_paths_three_1[k]=="HLT_HT430_650"){
                     if( nat.triggerReader().getTrgResult(trigger_paths_three_2[k]) || 
                         nat.triggerReader().getTrgResult(trigger_paths_three_3[k]) ||
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_three[k]=ntrig_three[k]+1;}
                }
                //Case 2
                if(trigger_paths_three_2[k]=="HLT_HT430_650"){
                     if( nat.triggerReader().getTrgResult(trigger_paths_three_1[k]) || 
                         nat.triggerReader().getTrgResult(trigger_paths_three_3[k]) ||
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_three[k]=ntrig_three[k]+1;}
                }
                //Case 3
                if(trigger_paths_three_3[k]=="HLT_HT430_650"){
                     if( nat.triggerReader().getTrgResult(trigger_paths_three_1[k]) || 
                         nat.triggerReader().getTrgResult(trigger_paths_three_2[k]) ||
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_three[k]=ntrig_three[k]+1;}
                }  
          }
          //Case 4
          else
          {
               if( nat.triggerReader().getTrgResult(trigger_paths_three_1[k]) || 
                   nat.triggerReader().getTrgResult(trigger_paths_three_2[k]) || 
                   nat.triggerReader().getTrgResult(trigger_paths_three_3[k]) ){ntrig_three[k]=ntrig_three[k]+1;}
          }

        }

        // Quadruple OR
        for (unsigned k=0; k<trigger_paths_four_1.size(); k++){
          //Case 1,2,3  
          if (trigger_paths_four_1[k]=="HLT_HT430_650" || trigger_paths_four_2[k]=="HLT_HT430_650" ||
              trigger_paths_four_3[k]=="HLT_HT430_650" || trigger_paths_four_4[k]=="HLT_HT430_650")
          {
                //Case 1
                if(trigger_paths_four_1[k]=="HLT_HT430_650"){
                     if( nat.triggerReader().getTrgResult(trigger_paths_four_2[k]) || 
                         nat.triggerReader().getTrgResult(trigger_paths_four_3[k]) ||
                         nat.triggerReader().getTrgResult(trigger_paths_four_4[k]) ||
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_four[k]=ntrig_four[k]+1;}
                }
                //Case 2
                if(trigger_paths_four_2[k]=="HLT_HT430_650"){
                     if( nat.triggerReader().getTrgResult(trigger_paths_four_1[k]) || 
                         nat.triggerReader().getTrgResult(trigger_paths_four_3[k]) ||
                         nat.triggerReader().getTrgResult(trigger_paths_four_4[k]) ||
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_four[k]=ntrig_four[k]+1;}
                }
                //Case 3
                if(trigger_paths_four_3[k]=="HLT_HT430_650"){
                     if( nat.triggerReader().getTrgResult(trigger_paths_four_1[k]) || 
                         nat.triggerReader().getTrgResult(trigger_paths_four_2[k]) ||
                         nat.triggerReader().getTrgResult(trigger_paths_four_4[k]) ||
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_four[k]=ntrig_four[k]+1;}
                }
                //Case 4
                if(trigger_paths_four_4[k]=="HLT_HT430_650"){
                     if( nat.triggerReader().getTrgResult(trigger_paths_four_1[k]) || 
                         nat.triggerReader().getTrgResult(trigger_paths_four_2[k]) ||
                         nat.triggerReader().getTrgResult(trigger_paths_four_3[k]) ||
                         nat.triggerReader().getTrgResult("HLT_HT430to450") ||
                         nat.triggerReader().getTrgResult("HLT_HT450to470") ||
                         nat.triggerReader().getTrgResult("HLT_HT470to500") ||
                         nat.triggerReader().getTrgResult("HLT_HT500to550") ||
                         nat.triggerReader().getTrgResult("HLT_HT550to650") ||
                         nat.triggerReader().getTrgResult("HLT_HT650")){ntrig_four[k]=ntrig_four[k]+1;}
                }
          }             
          else
          {
           if ( nat.triggerReader().getTrgResult(trigger_paths_four_1[k]) 
              || nat.triggerReader().getTrgResult(trigger_paths_four_2[k]) 
              || nat.triggerReader().getTrgResult(trigger_paths_four_3[k])
              || nat.triggerReader().getTrgResult(trigger_paths_four_4[k])){ ntrig_four[k]=ntrig_four[k]+1;}
          }

       }   
   }
    double ht_eff = double(nHT)/ntot;
    trigger_effs.push_back(ht_eff);

    cout << "[INFO] ... making Trigger Optimization Report (efficiencies_report.txt)" << endl;
    vector<tuple<double,string>>  result_1;
    for (unsigned p=0; p<trigger_paths.size(); p++){
    //if (trigger_paths[p]=="HLT_HT430_650") continue;
    result_1.push_back(tuple<double, string>(trigger_effs[p],trigger_paths[p]));    
    }
    std::sort( std::rbegin( result_1 ), std::rend( result_1 ) );
    vector<tuple<double,string,string>>  result_2;
    for (unsigned p=0; p<trigger_paths_two_1.size(); p++){
    //if (trigger_paths_two_1[p]=="HLT_HT430_650" || trigger_paths_two_2[p]=="HLT_HT430_650") continue;   
    result_2.push_back(tuple<double, string, string>((double)ntrig_two[p]/ntot,trigger_paths_two_1[p],trigger_paths_two_2[p])); 
    }
    std::sort( std::rbegin( result_2 ), std::rend( result_2 ) );
    vector<tuple<double,string,string,string>>  result_3;
    for (unsigned p=0; p<trigger_paths_three_1.size(); p++){
    //if (trigger_paths_three_1[p]=="HLT_HT430_650" || trigger_paths_three_2[p]=="HLT_HT430_650" || trigger_paths_three_3[p]=="HLT_HT430_650") continue;         
    result_3.push_back(tuple<double,string, string, string>((double)ntrig_three[p]/ntot,trigger_paths_three_1[p],trigger_paths_three_2[p],trigger_paths_three_3[p])); 
    }
    std::sort( std::rbegin( result_3 ), std::rend( result_3 ) );
    vector<tuple<double,string,string,string,string>>  result_4;
    for (unsigned p=0; p<trigger_paths_four_1.size(); p++){
    //if (trigger_paths_four_1[p]=="HLT_HT430_650" || trigger_paths_four_2[p]=="HLT_HT430_650" 
     //|| trigger_paths_four_3[p]=="HLT_HT430_650" || trigger_paths_four_4[p]=="HLT_HT430_650") continue;  
    result_4.push_back(tuple<double,string, string, string, string>((double)ntrig_four[p]/ntot,trigger_paths_four_1[p],trigger_paths_four_2[p],trigger_paths_four_3[p],trigger_paths_four_4[p]));
    }
    std::sort( std::rbegin( result_4 ), std::rend( result_4 ) );

    ofstream outputReport;
    outputReport.open("efficiencies_report_HT.txt");
    outputReport << "" << endl;outputReport << setw(100)<<"*** Optimal OR Trigger Efficiencies ***" << endl;outputReport << "" << endl;  
    outputReport << "" << endl;outputReport << "[RESULTS] Optimal OR Trigger Efficiencies using one HLT path" << endl;outputReport << "" << endl;
    for ( const auto &i : result_1 )
    {
    //double eff1=get<0>(i);   
//    if (eff1 < 0.8) continue;
    cout << get<0>(i) << get<1>(i) << endl;   
    outputReport <<"Eff = "<<setw(6)<< get<0>(i)*100<<'%'<<setw(15)<<"PATH:"<< setw(70)<< get<1>(i)<<endl;
    }

    outputReport << "" << endl;outputReport << "[RESULTS] Optimal OR Trigger Efficiencies using two HLT paths (eff>=0.17)" << endl;outputReport << "" << endl;
    for ( const auto &i : result_2 )
    {
    double eff2=get<0>(i);   
    if (eff2 < 0.17) continue;
    cout << get<0>(i) << get<1>(i) << get<2>(i) << endl;
    outputReport <<"Eff = "<<setw(6)<< get<0>(i)*100<<'%'<<setw(15)<<"PATHS:"<< setw(70)<< get<1>(i) << setw(70)<< get<2>(i) <<endl;    
    }
    outputReport << "" << endl;outputReport << "[RESULTS] Optimal OR Trigger Efficiencies using three HLT paths (eff>=0.24)" << endl;outputReport << "" << endl;//0.22
    for ( const auto &i : result_3 )
    {
    double eff3=get<0>(i);   
    if (eff3 < 0.24) continue;
    cout << get<0>(i) << get<1>(i) << get<2>(i) << get<3>(i) << endl;
    outputReport <<"Eff = "<<setw(6)<< get<0>(i)*100<<'%'<<setw(15)<<"PATHS:"<< setw(70)<< get<1>(i) << setw(70)<< get<2>(i) << setw(70)<< get<3>(i) <<endl;    
    }
    outputReport << "" << endl;outputReport << "[RESULTS] Optimal OR Trigger Efficiencies using four HLT paths (eff>=0.30)" << endl;outputReport << "" << endl;//0.24
    for ( const auto &i : result_4 )
    {
    double eff4=get<0>(i);   
    if (eff4 < 0.30) continue;
    cout << get<0>(i) << get<1>(i) << get<2>(i) << get<3>(i) << get<4>(i) << endl;
    outputReport <<"Eff = "<<setw(6)<< get<0>(i)*100<<'%'<<setw(15)<<"PATHS:"<< setw(70)<< get<1>(i) << setw(70)<< get<2>(i) << setw(70)<< get<3>(i) << setw(70)<< get<4>(i) <<endl;    
    }
    outputReport.close();
    
}

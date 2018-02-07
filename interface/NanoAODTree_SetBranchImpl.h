#ifndef NANOAODTREE_SETBRANCHIMPL
#define NANOAODTREE_SETBRANCHIMPL

/**
** class  : NanoAODTree_SetBranchImpl
** author : L. Cadamuro (UF)
** date   : 31/01/2018
** brief  : interface to the Nano AOD Event ntuple.
**        : implements the read via setbranchaddress instructions
** NOTE1  : you can do SetBranchStatus(name, 0) on some branches to speed uo
**        : but if you will try to access them in future, there is no warning that they are inactive and undefined values are returned!
**        :
** NOTE2  : IMPORTANT: NanoAOD use plain arrays instead of vectors to store variables (che minchiata colossale)
**        : so the vector where one stores the variable must have a size >= than the array to be stored itself
**        : the NanoAODTree makes a sort of check after calling GetEntry(), but at that point the memory has already been written and anuthing bad can happen
**        : so CHECK CAREFULLY the size of the ReaderArrayPatch below and increase it if necessary
**        : (the only way to avoid the problem could be the ReaderArray implementation, but fails due to a bug in ROOT as of Feb 01 2018)
**/

#include <vector>
#include <iostream>
#include "TChain.h"
#include "ReaderValuePatch.h"
#include "ReaderArrayPatch.h"
#include "TriggerReader_SetBranchImpl.h"

class NanoAODTree_SetBranchImpl {
    public:
        NanoAODTree_SetBranchImpl  (TChain* chain, bool isData = false):
            chain_(chain),
            trg_reader_(chain),
            old_tree_nr_(-1),
            nEv_(-1),
            is_data_ (isData)
        {}
        ~NanoAODTree_SetBranchImpl (){};

        TriggerReader_SetBranchImpl& triggerReader(){return trg_reader_;}

        // bool Next() {return fReader.Next();}
        bool Next();

        TChain* chain_;

        // the trigger reader - a special care is needed for branches that can disappear
        TriggerReader_SetBranchImpl trg_reader_;

        int old_tree_nr_; // the current tree
        Long64_t nEv_; // the current event
        bool is_data_;

        // declare readers
        // you can get the following list with the script generate_NanoAODTree_SetBranchImpl_vars.py
        // watch out for the max size of the array
        ReaderValuePatch<UInt_t>      run                                 {chain_, "run"};
        ReaderValuePatch<UInt_t>      luminosityBlock                     {chain_, "luminosityBlock"};
        ReaderValuePatch<ULong64_t>   event                               {chain_, "event"};
        ReaderValuePatch<Float_t>     CaloMET_phi                         {chain_, "CaloMET_phi"};
        ReaderValuePatch<Float_t>     CaloMET_pt                          {chain_, "CaloMET_pt"};
        ReaderValuePatch<Float_t>     CaloMET_sumEt                       {chain_, "CaloMET_sumEt"};
        ReaderValuePatch<UInt_t>      nElectron                           {chain_, "nElectron"};
        ReaderArrayPatch<Float_t>     Electron_deltaEtaSC                 {chain_, "Electron_deltaEtaSC", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_dr03EcalRecHitSumEt        {chain_, "Electron_dr03EcalRecHitSumEt", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_dr03HcalDepth1TowerSumEt   {chain_, "Electron_dr03HcalDepth1TowerSumEt", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_dr03TkSumPt                {chain_, "Electron_dr03TkSumPt", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_dxy                        {chain_, "Electron_dxy", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_dxyErr                     {chain_, "Electron_dxyErr", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_dz                         {chain_, "Electron_dz", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_dzErr                      {chain_, "Electron_dzErr", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_eCorr                      {chain_, "Electron_eCorr", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_eInvMinusPInv              {chain_, "Electron_eInvMinusPInv", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_energyErr                  {chain_, "Electron_energyErr", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_eta                        {chain_, "Electron_eta", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_hoe                        {chain_, "Electron_hoe", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_ip3d                       {chain_, "Electron_ip3d", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_mass                       {chain_, "Electron_mass", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_miniPFRelIso_all           {chain_, "Electron_miniPFRelIso_all", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_miniPFRelIso_chg           {chain_, "Electron_miniPFRelIso_chg", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_mvaSpring16GP              {chain_, "Electron_mvaSpring16GP", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_mvaSpring16HZZ             {chain_, "Electron_mvaSpring16HZZ", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_pfRelIso03_all             {chain_, "Electron_pfRelIso03_all", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_pfRelIso03_chg             {chain_, "Electron_pfRelIso03_chg", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_phi                        {chain_, "Electron_phi", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_pt                         {chain_, "Electron_pt", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_r9                         {chain_, "Electron_r9", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_sieie                      {chain_, "Electron_sieie", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_sip3d                      {chain_, "Electron_sip3d", nElectron, 100};
        ReaderArrayPatch<Float_t>     Electron_mvaTTH                     {chain_, "Electron_mvaTTH", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_charge                     {chain_, "Electron_charge", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_cutBased                   {chain_, "Electron_cutBased", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_cutBased_HLTPreSel         {chain_, "Electron_cutBased_HLTPreSel", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_jetIdx                     {chain_, "Electron_jetIdx", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_pdgId                      {chain_, "Electron_pdgId", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_photonIdx                  {chain_, "Electron_photonIdx", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_tightCharge                {chain_, "Electron_tightCharge", nElectron, 100};
        ReaderArrayPatch<Int_t>       Electron_vidNestedWPBitmap          {chain_, "Electron_vidNestedWPBitmap", nElectron, 100};
        ReaderArrayPatch<Bool_t>      Electron_convVeto                   {chain_, "Electron_convVeto", nElectron, 100};
        ReaderArrayPatch<Bool_t>      Electron_cutBased_HEEP              {chain_, "Electron_cutBased_HEEP", nElectron, 100};
        ReaderArrayPatch<Bool_t>      Electron_isPFcand                   {chain_, "Electron_isPFcand", nElectron, 100};
        ReaderArrayPatch<UChar_t>     Electron_lostHits                   {chain_, "Electron_lostHits", nElectron, 100};
        ReaderArrayPatch<Bool_t>      Electron_mvaSpring16GP_WP80         {chain_, "Electron_mvaSpring16GP_WP80", nElectron, 100};
        ReaderArrayPatch<Bool_t>      Electron_mvaSpring16GP_WP90         {chain_, "Electron_mvaSpring16GP_WP90", nElectron, 100};
        ReaderArrayPatch<Bool_t>      Electron_mvaSpring16HZZ_WPL         {chain_, "Electron_mvaSpring16HZZ_WPL", nElectron, 100};
        ReaderValuePatch<UInt_t>      nFatJet                             {chain_, "nFatJet"};
        ReaderArrayPatch<Float_t>     FatJet_area                         {chain_, "FatJet_area", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_btagCMVA                     {chain_, "FatJet_btagCMVA", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_btagCSVV2                    {chain_, "FatJet_btagCSVV2", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_btagDeepB                    {chain_, "FatJet_btagDeepB", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_btagHbb                      {chain_, "FatJet_btagHbb", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_eta                          {chain_, "FatJet_eta", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_mass                         {chain_, "FatJet_mass", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_msoftdrop                    {chain_, "FatJet_msoftdrop", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_n2b1                         {chain_, "FatJet_n2b1", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_n3b1                         {chain_, "FatJet_n3b1", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_phi                          {chain_, "FatJet_phi", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_pt                           {chain_, "FatJet_pt", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_tau1                         {chain_, "FatJet_tau1", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_tau2                         {chain_, "FatJet_tau2", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_tau3                         {chain_, "FatJet_tau3", nFatJet, 100};
        ReaderArrayPatch<Float_t>     FatJet_tau4                         {chain_, "FatJet_tau4", nFatJet, 100};
        ReaderArrayPatch<Int_t>       FatJet_subJetIdx1                   {chain_, "FatJet_subJetIdx1", nFatJet, 100};
        ReaderArrayPatch<Int_t>       FatJet_subJetIdx2                   {chain_, "FatJet_subJetIdx2", nFatJet, 100};
        ReaderValuePatch<UInt_t>      nGenJetAK8                          {chain_, "nGenJetAK8"};
        ReaderArrayPatch<Float_t>     GenJetAK8_eta                       {chain_, "GenJetAK8_eta", nGenJetAK8, 100};
        ReaderArrayPatch<Float_t>     GenJetAK8_mass                      {chain_, "GenJetAK8_mass", nGenJetAK8, 100};
        ReaderArrayPatch<Float_t>     GenJetAK8_phi                       {chain_, "GenJetAK8_phi", nGenJetAK8, 100};
        ReaderArrayPatch<Float_t>     GenJetAK8_pt                        {chain_, "GenJetAK8_pt", nGenJetAK8, 100};
        ReaderValuePatch<UInt_t>      nGenJet                             {chain_, "nGenJet"};
        ReaderArrayPatch<Float_t>     GenJet_eta                          {chain_, "GenJet_eta", nGenJet, 100};
        ReaderArrayPatch<Float_t>     GenJet_mass                         {chain_, "GenJet_mass", nGenJet, 100};
        ReaderArrayPatch<Float_t>     GenJet_phi                          {chain_, "GenJet_phi", nGenJet, 100};
        ReaderArrayPatch<Float_t>     GenJet_pt                           {chain_, "GenJet_pt", nGenJet, 100};
        ReaderValuePatch<UInt_t>      nGenPart                            {chain_, "nGenPart"};
        ReaderArrayPatch<Float_t>     GenPart_eta                         {chain_, "GenPart_eta", nGenPart, 400};
        ReaderArrayPatch<Float_t>     GenPart_mass                        {chain_, "GenPart_mass", nGenPart, 400};
        ReaderArrayPatch<Float_t>     GenPart_phi                         {chain_, "GenPart_phi", nGenPart, 400};
        ReaderArrayPatch<Float_t>     GenPart_pt                          {chain_, "GenPart_pt", nGenPart, 400};
        ReaderArrayPatch<Int_t>       GenPart_genPartIdxMother            {chain_, "GenPart_genPartIdxMother", nGenPart, 400};
        ReaderArrayPatch<Int_t>       GenPart_pdgId                       {chain_, "GenPart_pdgId", nGenPart, 400};
        ReaderArrayPatch<Int_t>       GenPart_status                      {chain_, "GenPart_status", nGenPart, 400};
        ReaderArrayPatch<Int_t>       GenPart_statusFlags                 {chain_, "GenPart_statusFlags", nGenPart, 400};
        ReaderValuePatch<Float_t>     Generator_scalePDF                  {chain_, "Generator_scalePDF"};
        ReaderValuePatch<Float_t>     Generator_x1                        {chain_, "Generator_x1"};
        ReaderValuePatch<Float_t>     Generator_x2                        {chain_, "Generator_x2"};
        ReaderValuePatch<Float_t>     Generator_xpdf1                     {chain_, "Generator_xpdf1"};
        ReaderValuePatch<Float_t>     Generator_xpdf2                     {chain_, "Generator_xpdf2"};
        ReaderValuePatch<Int_t>       Generator_id1                       {chain_, "Generator_id1"};
        ReaderValuePatch<Int_t>       Generator_id2                       {chain_, "Generator_id2"};
        ReaderValuePatch<UInt_t>      nGenVisTau                          {chain_, "nGenVisTau"};
        ReaderArrayPatch<Float_t>     GenVisTau_eta                       {chain_, "GenVisTau_eta", nGenVisTau, 100};
        ReaderArrayPatch<Float_t>     GenVisTau_mass                      {chain_, "GenVisTau_mass", nGenVisTau, 100};
        ReaderArrayPatch<Float_t>     GenVisTau_phi                       {chain_, "GenVisTau_phi", nGenVisTau, 100};
        ReaderArrayPatch<Float_t>     GenVisTau_pt                        {chain_, "GenVisTau_pt", nGenVisTau, 100};
        ReaderArrayPatch<Int_t>       GenVisTau_charge                    {chain_, "GenVisTau_charge", nGenVisTau, 100};
        ReaderArrayPatch<Int_t>       GenVisTau_genPartIdxMother          {chain_, "GenVisTau_genPartIdxMother", nGenVisTau, 100};
        ReaderArrayPatch<Int_t>       GenVisTau_status                    {chain_, "GenVisTau_status", nGenVisTau, 100};
        ReaderValuePatch<Float_t>     genWeight                           {chain_, "genWeight"};
        ReaderValuePatch<Float_t>     LHEWeight_originalXWGTUP            {chain_, "LHEWeight_originalXWGTUP"};
        ReaderValuePatch<UInt_t>      nLHEPdfWeight                       {chain_, "nLHEPdfWeight"};
        ReaderArrayPatch<Float_t>     LHEPdfWeight                        {chain_, "LHEPdfWeight", nLHEPdfWeight, 100};
        ReaderValuePatch<UInt_t>      nLHEScaleWeight                     {chain_, "nLHEScaleWeight"};
        ReaderArrayPatch<Float_t>     LHEScaleWeight                      {chain_, "LHEScaleWeight", nLHEScaleWeight, 100};
        ReaderValuePatch<UInt_t>      nJet                                {chain_, "nJet"};
        ReaderArrayPatch<Float_t>     Jet_area                            {chain_, "Jet_area", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_btagCMVA                        {chain_, "Jet_btagCMVA", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_btagCSVV2                       {chain_, "Jet_btagCSVV2", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_btagDeepB                       {chain_, "Jet_btagDeepB", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_btagDeepC                       {chain_, "Jet_btagDeepC", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_chEmEF                          {chain_, "Jet_chEmEF", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_chHEF                           {chain_, "Jet_chHEF", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_eta                             {chain_, "Jet_eta", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_mass                            {chain_, "Jet_mass", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_neEmEF                          {chain_, "Jet_neEmEF", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_neHEF                           {chain_, "Jet_neHEF", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_phi                             {chain_, "Jet_phi", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_pt                              {chain_, "Jet_pt", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_qgl                             {chain_, "Jet_qgl", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_rawFactor                       {chain_, "Jet_rawFactor", nJet, 100};
        ReaderArrayPatch<Float_t>     Jet_bReg                            {chain_, "Jet_bReg", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_electronIdx1                    {chain_, "Jet_electronIdx1", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_electronIdx2                    {chain_, "Jet_electronIdx2", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_jetId                           {chain_, "Jet_jetId", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_muonIdx1                        {chain_, "Jet_muonIdx1", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_muonIdx2                        {chain_, "Jet_muonIdx2", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_nConstituents                   {chain_, "Jet_nConstituents", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_nElectrons                      {chain_, "Jet_nElectrons", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_nMuons                          {chain_, "Jet_nMuons", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_puId                            {chain_, "Jet_puId", nJet, 100};
        ReaderValuePatch<Float_t>     LHE_HT                              {chain_, "LHE_HT"};
        ReaderValuePatch<Float_t>     LHE_HTIncoming                      {chain_, "LHE_HTIncoming"};
        ReaderValuePatch<Float_t>     LHE_Vpt                             {chain_, "LHE_Vpt"};
        ReaderValuePatch<UChar_t>     LHE_Njets                           {chain_, "LHE_Njets"};
        ReaderValuePatch<UChar_t>     LHE_Nb                              {chain_, "LHE_Nb"};
        ReaderValuePatch<UChar_t>     LHE_Nc                              {chain_, "LHE_Nc"};
        ReaderValuePatch<UChar_t>     LHE_Nuds                            {chain_, "LHE_Nuds"};
        ReaderValuePatch<UChar_t>     LHE_Nglu                            {chain_, "LHE_Nglu"};
        ReaderValuePatch<UChar_t>     LHE_NpNLO                           {chain_, "LHE_NpNLO"};
        ReaderValuePatch<UChar_t>     LHE_NpLO                            {chain_, "LHE_NpLO"};
        ReaderValuePatch<Float_t>     GenMET_phi                          {chain_, "GenMET_phi"};
        ReaderValuePatch<Float_t>     GenMET_pt                           {chain_, "GenMET_pt"};
        ReaderValuePatch<Float_t>     MET_MetUnclustEnUpDeltaX            {chain_, "MET_MetUnclustEnUpDeltaX"};
        ReaderValuePatch<Float_t>     MET_MetUnclustEnUpDeltaY            {chain_, "MET_MetUnclustEnUpDeltaY"};
        ReaderValuePatch<Float_t>     MET_covXX                           {chain_, "MET_covXX"};
        ReaderValuePatch<Float_t>     MET_covXY                           {chain_, "MET_covXY"};
        ReaderValuePatch<Float_t>     MET_covYY                           {chain_, "MET_covYY"};
        ReaderValuePatch<Float_t>     MET_phi                             {chain_, "MET_phi"};
        ReaderValuePatch<Float_t>     MET_pt                              {chain_, "MET_pt"};
        ReaderValuePatch<Float_t>     MET_significance                    {chain_, "MET_significance"};
        ReaderValuePatch<Float_t>     MET_sumEt                           {chain_, "MET_sumEt"};
        ReaderValuePatch<UInt_t>      nMuon                               {chain_, "nMuon"};
        ReaderArrayPatch<Float_t>     Muon_dxy                            {chain_, "Muon_dxy", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_dxyErr                         {chain_, "Muon_dxyErr", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_dz                             {chain_, "Muon_dz", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_dzErr                          {chain_, "Muon_dzErr", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_eta                            {chain_, "Muon_eta", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_ip3d                           {chain_, "Muon_ip3d", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_mass                           {chain_, "Muon_mass", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_miniPFRelIso_all               {chain_, "Muon_miniPFRelIso_all", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_miniPFRelIso_chg               {chain_, "Muon_miniPFRelIso_chg", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_pfRelIso03_all                 {chain_, "Muon_pfRelIso03_all", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_pfRelIso03_chg                 {chain_, "Muon_pfRelIso03_chg", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_pfRelIso04_all                 {chain_, "Muon_pfRelIso04_all", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_phi                            {chain_, "Muon_phi", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_pt                             {chain_, "Muon_pt", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_ptErr                          {chain_, "Muon_ptErr", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_segmentComp                    {chain_, "Muon_segmentComp", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_sip3d                          {chain_, "Muon_sip3d", nMuon, 100};
        ReaderArrayPatch<Float_t>     Muon_mvaTTH                         {chain_, "Muon_mvaTTH", nMuon, 100};
        ReaderArrayPatch<Int_t>       Muon_charge                         {chain_, "Muon_charge", nMuon, 100};
        ReaderArrayPatch<Int_t>       Muon_jetIdx                         {chain_, "Muon_jetIdx", nMuon, 100};
        ReaderArrayPatch<Int_t>       Muon_nStations                      {chain_, "Muon_nStations", nMuon, 100};
        ReaderArrayPatch<Int_t>       Muon_nTrackerLayers                 {chain_, "Muon_nTrackerLayers", nMuon, 100};
        ReaderArrayPatch<Int_t>       Muon_pdgId                          {chain_, "Muon_pdgId", nMuon, 100};
        ReaderArrayPatch<Int_t>       Muon_tightCharge                    {chain_, "Muon_tightCharge", nMuon, 100};
        ReaderArrayPatch<UChar_t>     Muon_highPtId                       {chain_, "Muon_highPtId", nMuon, 100};
        ReaderArrayPatch<Bool_t>      Muon_isPFcand                       {chain_, "Muon_isPFcand", nMuon, 100};
        ReaderArrayPatch<Bool_t>      Muon_mediumId                       {chain_, "Muon_mediumId", nMuon, 100};
        ReaderArrayPatch<Bool_t>      Muon_softId                         {chain_, "Muon_softId", nMuon, 100};
        ReaderArrayPatch<Bool_t>      Muon_tightId                        {chain_, "Muon_tightId", nMuon, 100};
        ReaderValuePatch<UInt_t>      nPhoton                             {chain_, "nPhoton"};
        ReaderArrayPatch<Float_t>     Photon_eCorr                        {chain_, "Photon_eCorr", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_energyErr                    {chain_, "Photon_energyErr", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_eta                          {chain_, "Photon_eta", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_hoe                          {chain_, "Photon_hoe", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_mass                         {chain_, "Photon_mass", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_mvaID                        {chain_, "Photon_mvaID", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_pfRelIso03_all               {chain_, "Photon_pfRelIso03_all", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_pfRelIso03_chg               {chain_, "Photon_pfRelIso03_chg", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_phi                          {chain_, "Photon_phi", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_pt                           {chain_, "Photon_pt", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_r9                           {chain_, "Photon_r9", nPhoton, 100};
        ReaderArrayPatch<Float_t>     Photon_sieie                        {chain_, "Photon_sieie", nPhoton, 100};
        ReaderArrayPatch<Int_t>       Photon_charge                       {chain_, "Photon_charge", nPhoton, 100};
        ReaderArrayPatch<Int_t>       Photon_cutBased                     {chain_, "Photon_cutBased", nPhoton, 100};
        ReaderArrayPatch<Int_t>       Photon_electronIdx                  {chain_, "Photon_electronIdx", nPhoton, 100};
        ReaderArrayPatch<Int_t>       Photon_jetIdx                       {chain_, "Photon_jetIdx", nPhoton, 100};
        ReaderArrayPatch<Int_t>       Photon_pdgId                        {chain_, "Photon_pdgId", nPhoton, 100};
        ReaderArrayPatch<Int_t>       Photon_vidNestedWPBitmap            {chain_, "Photon_vidNestedWPBitmap", nPhoton, 100};
        ReaderArrayPatch<Bool_t>      Photon_electronVeto                 {chain_, "Photon_electronVeto", nPhoton, 100};
        ReaderArrayPatch<Bool_t>      Photon_mvaID_WP80                   {chain_, "Photon_mvaID_WP80", nPhoton, 100};
        ReaderArrayPatch<Bool_t>      Photon_mvaID_WP90                   {chain_, "Photon_mvaID_WP90", nPhoton, 100};
        ReaderArrayPatch<Bool_t>      Photon_pixelSeed                    {chain_, "Photon_pixelSeed", nPhoton, 100};
        ReaderValuePatch<Int_t>       Pileup_nPU                          {chain_, "Pileup_nPU"};
        ReaderValuePatch<Int_t>       Pileup_nTrueInt                     {chain_, "Pileup_nTrueInt"};
        ReaderValuePatch<Float_t>     PuppiMET_phi                        {chain_, "PuppiMET_phi"};
        ReaderValuePatch<Float_t>     PuppiMET_pt                         {chain_, "PuppiMET_pt"};
        ReaderValuePatch<Float_t>     PuppiMET_sumEt                      {chain_, "PuppiMET_sumEt"};
        ReaderValuePatch<Float_t>     RawMET_phi                          {chain_, "RawMET_phi"};
        ReaderValuePatch<Float_t>     RawMET_pt                           {chain_, "RawMET_pt"};
        ReaderValuePatch<Float_t>     RawMET_sumEt                        {chain_, "RawMET_sumEt"};
        ReaderValuePatch<Float_t>     fixedGridRhoFastjetAll              {chain_, "fixedGridRhoFastjetAll"};
        ReaderValuePatch<Float_t>     fixedGridRhoFastjetCentralCalo      {chain_, "fixedGridRhoFastjetCentralCalo"};
        ReaderValuePatch<Float_t>     fixedGridRhoFastjetCentralNeutral   {chain_, "fixedGridRhoFastjetCentralNeutral"};
        ReaderValuePatch<UInt_t>      nGenDressedLepton                   {chain_, "nGenDressedLepton"};
        ReaderArrayPatch<Float_t>     GenDressedLepton_eta                {chain_, "GenDressedLepton_eta", nGenDressedLepton, 100};
        ReaderArrayPatch<Float_t>     GenDressedLepton_mass               {chain_, "GenDressedLepton_mass", nGenDressedLepton, 100};
        ReaderArrayPatch<Float_t>     GenDressedLepton_phi                {chain_, "GenDressedLepton_phi", nGenDressedLepton, 100};
        ReaderArrayPatch<Float_t>     GenDressedLepton_pt                 {chain_, "GenDressedLepton_pt", nGenDressedLepton, 100};
        ReaderArrayPatch<Int_t>       GenDressedLepton_pdgId              {chain_, "GenDressedLepton_pdgId", nGenDressedLepton, 100};
        ReaderValuePatch<UInt_t>      nSoftActivityJet                    {chain_, "nSoftActivityJet"};
        ReaderArrayPatch<Float_t>     SoftActivityJet_eta                 {chain_, "SoftActivityJet_eta", nSoftActivityJet, 100};
        ReaderArrayPatch<Float_t>     SoftActivityJet_phi                 {chain_, "SoftActivityJet_phi", nSoftActivityJet, 100};
        ReaderArrayPatch<Float_t>     SoftActivityJet_pt                  {chain_, "SoftActivityJet_pt", nSoftActivityJet, 100};
        ReaderValuePatch<Float_t>     SoftActivityJetHT                   {chain_, "SoftActivityJetHT"};
        ReaderValuePatch<Float_t>     SoftActivityJetHT10                 {chain_, "SoftActivityJetHT10"};
        ReaderValuePatch<Float_t>     SoftActivityJetHT2                  {chain_, "SoftActivityJetHT2"};
        ReaderValuePatch<Float_t>     SoftActivityJetHT5                  {chain_, "SoftActivityJetHT5"};
        ReaderValuePatch<Int_t>       SoftActivityJetNjets10              {chain_, "SoftActivityJetNjets10"};
        ReaderValuePatch<Int_t>       SoftActivityJetNjets2               {chain_, "SoftActivityJetNjets2"};
        ReaderValuePatch<Int_t>       SoftActivityJetNjets5               {chain_, "SoftActivityJetNjets5"};
        ReaderValuePatch<UInt_t>      nSubJet                             {chain_, "nSubJet"};
        ReaderArrayPatch<Float_t>     SubJet_btagCMVA                     {chain_, "SubJet_btagCMVA", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_btagCSVV2                    {chain_, "SubJet_btagCSVV2", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_btagDeepB                    {chain_, "SubJet_btagDeepB", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_eta                          {chain_, "SubJet_eta", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_mass                         {chain_, "SubJet_mass", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_n2b1                         {chain_, "SubJet_n2b1", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_n3b1                         {chain_, "SubJet_n3b1", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_phi                          {chain_, "SubJet_phi", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_pt                           {chain_, "SubJet_pt", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_tau1                         {chain_, "SubJet_tau1", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_tau2                         {chain_, "SubJet_tau2", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_tau3                         {chain_, "SubJet_tau3", nSubJet, 100};
        ReaderArrayPatch<Float_t>     SubJet_tau4                         {chain_, "SubJet_tau4", nSubJet, 100};
        ReaderValuePatch<UInt_t>      nTau                                {chain_, "nTau"};
        ReaderArrayPatch<Float_t>     Tau_chargedIso                      {chain_, "Tau_chargedIso", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_dxy                             {chain_, "Tau_dxy", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_dz                              {chain_, "Tau_dz", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_eta                             {chain_, "Tau_eta", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_footprintCorr                   {chain_, "Tau_footprintCorr", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_leadTkDeltaEta                  {chain_, "Tau_leadTkDeltaEta", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_leadTkDeltaPhi                  {chain_, "Tau_leadTkDeltaPhi", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_leadTkPtOverTauPt               {chain_, "Tau_leadTkPtOverTauPt", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_mass                            {chain_, "Tau_mass", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_neutralIso                      {chain_, "Tau_neutralIso", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_phi                             {chain_, "Tau_phi", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_photonsOutsideSignalCone        {chain_, "Tau_photonsOutsideSignalCone", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_pt                              {chain_, "Tau_pt", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_puCorr                          {chain_, "Tau_puCorr", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_rawAntiEle                      {chain_, "Tau_rawAntiEle", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_rawIso                          {chain_, "Tau_rawIso", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_rawMVAnewDM                     {chain_, "Tau_rawMVAnewDM", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_rawMVAoldDM                     {chain_, "Tau_rawMVAoldDM", nTau, 100};
        ReaderArrayPatch<Float_t>     Tau_rawMVAoldDMdR03                 {chain_, "Tau_rawMVAoldDMdR03", nTau, 100};
        ReaderArrayPatch<Int_t>       Tau_charge                          {chain_, "Tau_charge", nTau, 100};
        ReaderArrayPatch<Int_t>       Tau_decayMode                       {chain_, "Tau_decayMode", nTau, 100};
        ReaderArrayPatch<Int_t>       Tau_jetIdx                          {chain_, "Tau_jetIdx", nTau, 100};
        ReaderArrayPatch<Int_t>       Tau_rawAntiEleCat                   {chain_, "Tau_rawAntiEleCat", nTau, 100};
        ReaderArrayPatch<UChar_t>     Tau_idAntiEle                       {chain_, "Tau_idAntiEle", nTau, 100};
        ReaderArrayPatch<UChar_t>     Tau_idAntiMu                        {chain_, "Tau_idAntiMu", nTau, 100};
        ReaderArrayPatch<Bool_t>      Tau_idDecayMode                     {chain_, "Tau_idDecayMode", nTau, 100};
        ReaderArrayPatch<Bool_t>      Tau_idDecayModeNewDMs               {chain_, "Tau_idDecayModeNewDMs", nTau, 100};
        ReaderArrayPatch<UChar_t>     Tau_idMVAnewDM                      {chain_, "Tau_idMVAnewDM", nTau, 100};
        ReaderArrayPatch<UChar_t>     Tau_idMVAoldDM                      {chain_, "Tau_idMVAoldDM", nTau, 100};
        ReaderArrayPatch<UChar_t>     Tau_idMVAoldDMdR03                  {chain_, "Tau_idMVAoldDMdR03", nTau, 100};
        ReaderValuePatch<Float_t>     TkMET_phi                           {chain_, "TkMET_phi"};
        ReaderValuePatch<Float_t>     TkMET_pt                            {chain_, "TkMET_pt"};
        ReaderValuePatch<Float_t>     TkMET_sumEt                         {chain_, "TkMET_sumEt"};
        ReaderValuePatch<UInt_t>      nTrigObj                            {chain_, "nTrigObj"};
        ReaderArrayPatch<Float_t>     TrigObj_pt                          {chain_, "TrigObj_pt", nTrigObj, 100};
        ReaderArrayPatch<Float_t>     TrigObj_eta                         {chain_, "TrigObj_eta", nTrigObj, 100};
        ReaderArrayPatch<Float_t>     TrigObj_phi                         {chain_, "TrigObj_phi", nTrigObj, 100};
        ReaderArrayPatch<Float_t>     TrigObj_l1pt                        {chain_, "TrigObj_l1pt", nTrigObj, 100};
        ReaderArrayPatch<Float_t>     TrigObj_l1pt_2                      {chain_, "TrigObj_l1pt_2", nTrigObj, 100};
        ReaderArrayPatch<Float_t>     TrigObj_l2pt                        {chain_, "TrigObj_l2pt", nTrigObj, 100};
        ReaderArrayPatch<Int_t>       TrigObj_id                          {chain_, "TrigObj_id", nTrigObj, 100};
        ReaderArrayPatch<Int_t>       TrigObj_l1iso                       {chain_, "TrigObj_l1iso", nTrigObj, 100};
        ReaderArrayPatch<Int_t>       TrigObj_l1charge                    {chain_, "TrigObj_l1charge", nTrigObj, 100};
        ReaderArrayPatch<Int_t>       TrigObj_filterBits                  {chain_, "TrigObj_filterBits", nTrigObj, 100};
        ReaderValuePatch<Int_t>       genTtbarId                          {chain_, "genTtbarId"};
        ReaderValuePatch<UInt_t>      nOtherPV                            {chain_, "nOtherPV"};
        ReaderArrayPatch<Float_t>     OtherPV_z                           {chain_, "OtherPV_z", nOtherPV, 100};
        ReaderValuePatch<Float_t>     PV_ndof                             {chain_, "PV_ndof"};
        ReaderValuePatch<Float_t>     PV_x                                {chain_, "PV_x"};
        ReaderValuePatch<Float_t>     PV_y                                {chain_, "PV_y"};
        ReaderValuePatch<Float_t>     PV_z                                {chain_, "PV_z"};
        ReaderValuePatch<Float_t>     PV_chi2                             {chain_, "PV_chi2"};
        ReaderValuePatch<Float_t>     PV_score                            {chain_, "PV_score"};
        ReaderValuePatch<Int_t>       PV_npvs                             {chain_, "PV_npvs"};
        ReaderValuePatch<UInt_t>      nSV                                 {chain_, "nSV"};
        ReaderArrayPatch<Float_t>     SV_dlen                             {chain_, "SV_dlen", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_dlenSig                          {chain_, "SV_dlenSig", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_pAngle                           {chain_, "SV_pAngle", nSV, 100};
        ReaderArrayPatch<Int_t>       Electron_genPartIdx                 {chain_, "Electron_genPartIdx", nElectron, 100};
        ReaderArrayPatch<UChar_t>     Electron_genPartFlav                {chain_, "Electron_genPartFlav", nElectron, 100};
        ReaderArrayPatch<Int_t>       GenJetAK8_partonFlavour             {chain_, "GenJetAK8_partonFlavour", nGenJetAK8, 100};
        ReaderArrayPatch<UChar_t>     GenJetAK8_hadronFlavour             {chain_, "GenJetAK8_hadronFlavour", nGenJetAK8, 100};
        ReaderArrayPatch<Int_t>       GenJet_partonFlavour                {chain_, "GenJet_partonFlavour", nGenJet, 100};
        ReaderArrayPatch<UChar_t>     GenJet_hadronFlavour                {chain_, "GenJet_hadronFlavour", nGenJet, 100};
        ReaderArrayPatch<Int_t>       Jet_genJetIdx                       {chain_, "Jet_genJetIdx", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_hadronFlavour                   {chain_, "Jet_hadronFlavour", nJet, 100};
        ReaderArrayPatch<Int_t>       Jet_partonFlavour                   {chain_, "Jet_partonFlavour", nJet, 100};
        ReaderArrayPatch<Int_t>       Muon_genPartIdx                     {chain_, "Muon_genPartIdx", nMuon, 100};
        ReaderArrayPatch<UChar_t>     Muon_genPartFlav                    {chain_, "Muon_genPartFlav", nMuon, 100};
        ReaderArrayPatch<Int_t>       Photon_genPartIdx                   {chain_, "Photon_genPartIdx", nPhoton, 100};
        ReaderArrayPatch<UChar_t>     Photon_genPartFlav                  {chain_, "Photon_genPartFlav", nPhoton, 100};
        ReaderValuePatch<Float_t>     MET_fiducialGenPhi                  {chain_, "MET_fiducialGenPhi"};
        ReaderValuePatch<Float_t>     MET_fiducialGenPt                   {chain_, "MET_fiducialGenPt"};
        ReaderArrayPatch<UChar_t>     Electron_cleanmask                  {chain_, "Electron_cleanmask", nElectron, 100};
        ReaderArrayPatch<UChar_t>     Jet_cleanmask                       {chain_, "Jet_cleanmask", nJet, 100};
        ReaderArrayPatch<UChar_t>     Muon_cleanmask                      {chain_, "Muon_cleanmask", nMuon, 100};
        ReaderArrayPatch<UChar_t>     Photon_cleanmask                    {chain_, "Photon_cleanmask", nPhoton, 100};
        ReaderArrayPatch<UChar_t>     Tau_cleanmask                       {chain_, "Tau_cleanmask", nTau, 100};
        ReaderArrayPatch<Float_t>     SV_chi2                             {chain_, "SV_chi2", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_eta                              {chain_, "SV_eta", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_mass                             {chain_, "SV_mass", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_ndof                             {chain_, "SV_ndof", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_phi                              {chain_, "SV_phi", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_pt                               {chain_, "SV_pt", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_x                                {chain_, "SV_x", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_y                                {chain_, "SV_y", nSV, 100};
        ReaderArrayPatch<Float_t>     SV_z                                {chain_, "SV_z", nSV, 100};
        ReaderArrayPatch<Int_t>       Tau_genPartIdx                      {chain_, "Tau_genPartIdx", nTau, 100};
        ReaderArrayPatch<UChar_t>     Tau_genPartFlav                     {chain_, "Tau_genPartFlav", nTau, 100};
};

#endif
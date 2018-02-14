#include "OfflineProducerHelper.h"
#include <iostream>
#include <iomanip>
#include <cmath>

#include "CompositeCandidate.h"
#include "Jet.h"
#include "GenPart.h"

using namespace std;

std::vector<Jet> OfflineProducerHelper::make_jets(NanoAODTree& nat, const std::function<bool (Jet)>& presel_function)
{
    std::vector<Jet> jets;
    for (uint ij = 0; ij < *(nat.nJet); ++ij)
    {
        Jet jet (ij, &nat);
        bool pass = presel_function ? presel_function(jet) : true;
        if (pass)
            jets.push_back(jet);
    }
    return jets;
}

void OfflineProducerHelper::filter_jets(std::vector<Jet>& jets, const std::function<bool (Jet)>& filter_function)
{
    jets.erase(remove_if(jets.begin(), jets.end(), std::not1(filter_function)), jets.end());
    return;
}

bool OfflineProducerHelper::select_bbbb_jets(NanoAODTree& nat, EventInfo& ei, bbbbSelectionStrategy strat)
{
    // cout << "FIXME : select_bbbb_jets : this is just a dummy function " << endl;

    if (*(nat.nJet) < 4)
        return false;

    std::vector<Jet> jets;
    jets.reserve(*(nat.nJet));
    for (uint ij = 0; ij < *(nat.nJet); ++ij)
        jets.push_back(Jet(ij, &nat));

    stable_sort(jets.begin(), jets.end(), [](const Jet & a, const Jet & b) -> bool
    {
        return ( get_property(a, Jet_btagCSVV2) < get_property(b, Jet_btagCSVV2) );
    }); // sort by CSV (lowest to highest)

    // now need to pair the jets
    std::array<Jet, 4> presel_jets = {{
        *(jets.rbegin()+0),
        *(jets.rbegin()+1),
        *(jets.rbegin()+2),
        *(jets.rbegin()+3)
    }};

    std::array<Jet, 4> ordered_jets = presel_jets;

    // cout << "IN (pt) " << endl;
    // for (auto& j : presel_jets)
    //     cout << j.P4().Pt() << endl;

    // pair the jets
    if (strat == bbbbSelectionStrategy::kOneClosestToMh)
        ordered_jets = bbbb_jets_idxs_OneClosestToMh(presel_jets, 115.);
    
    else if (strat == bbbbSelectionStrategy::kBothClosestToMh)
        ordered_jets = bbbb_jets_idxs_BothClosestToMh(presel_jets, 115.);
    
    else if (strat == bbbbSelectionStrategy::kMostBackToBack)
        ordered_jets = bbbb_jets_idxs_MostBackToBack(presel_jets);
    else
        cout << "** [WARNING] : select_bbbb_jets : did not recognize bbbb selection strategy, pairing by pt" << endl;

    // cout << "OUT (pt) " << endl;
    // for (auto& j : ordered_jets)
    //     cout << j.P4().Pt() << endl;

    order_by_pT(ordered_jets.at(0), ordered_jets.at(1));
    order_by_pT(ordered_jets.at(2), ordered_jets.at(3));

    // ei.H1 = CompositeCandidate(ei.H1_b1.get(), ei.H1_b2.get());
    // ei.H2 = CompositeCandidate(ei.H2_b1.get(), ei.H2_b2.get());

    ei.H1 = CompositeCandidate(ordered_jets.at(0), ordered_jets.at(1));
    ei.H2 = CompositeCandidate(ordered_jets.at(2), ordered_jets.at(3));

    // order H1, H2 by pT: pT(H1) > pT (H2)
    bool swapped = order_by_pT(ei.H1.get(), ei.H2.get());

    if (!swapped)
    {
        ei.H1_b1 = ordered_jets.at(0);
        ei.H1_b2 = ordered_jets.at(1);
        ei.H2_b1 = ordered_jets.at(2);
        ei.H2_b2 = ordered_jets.at(3);
    }
    else
    {
        ei.H1_b1 = ordered_jets.at(2);
        ei.H1_b2 = ordered_jets.at(3);
        ei.H2_b1 = ordered_jets.at(0);
        ei.H2_b2 = ordered_jets.at(1);        
    }

    return true;
}

// one pair is closest to the Higgs mass, the other follows
std::array<Jet, 4> OfflineProducerHelper::bbbb_jets_idxs_OneClosestToMh(std::array<Jet, 4> presel_jets, float targetmH)
{
    // 12, 13, 14, 23, 24, 34    
    std::vector<pair <float, pair<int,int>>> mHs_and_jetIdx; // each entry is the mH and the two idxs of the pair
    
    for (uint i = 0; i < presel_jets.size(); ++i)
        for (uint j = i+1; j < presel_jets.size(); ++j)
        {
            TLorentzVector p4sum = (presel_jets.at(i).P4() + presel_jets.at(j).P4());
            float dmh = fabs(p4sum.Mag() - targetmH);
            mHs_and_jetIdx.push_back(make_pair(dmh, make_pair(i,j)));
        }

    // sort to get the pair closest to mH
    stable_sort (mHs_and_jetIdx.begin(), mHs_and_jetIdx.end());

    std::array<Jet, 4> outputJets = presel_jets;
    int ij0 = mHs_and_jetIdx.begin()->second.first;
    int ij1 = mHs_and_jetIdx.begin()->second.second;

    // get the other two jets. The following creates a vector with idxs 0/1/2/3, and then removes ij1 and ij2
    std::vector<int> vres;
    for (uint i = 0; i < presel_jets.size(); ++i) vres.push_back(i);
    vres.erase(std::remove(vres.begin(), vres.end(), ij0), vres.end());
    vres.erase(std::remove(vres.begin(), vres.end(), ij1), vres.end());
    
    int ij2 = vres.at(0);
    int ij3 = vres.at(1);

    outputJets.at(0) = presel_jets.at(ij0);
    outputJets.at(1) = presel_jets.at(ij1);
    outputJets.at(2) = presel_jets.at(ij2);
    outputJets.at(3) = presel_jets.at(ij3);

    return outputJets;
}

// minimize the distance from (targetmH, targetmH) in the 2D plane
std::array<Jet, 4> OfflineProducerHelper::bbbb_jets_idxs_BothClosestToMh(std::array<Jet, 4> presel_jets, float targetmH)
{
    // there are three possible pairings:
    // 12 34
    // 13 24
    // 14 23
    
    // 0   1   2   3   4   5
    // 12, 13, 14, 23, 24, 34    
    std::vector<float> mHs;
    
    for (uint i = 0; i < presel_jets.size(); ++i)
        for (uint j = i+1; j < presel_jets.size(); ++j)
        {
            TLorentzVector p4sum = (presel_jets.at(i).P4() + presel_jets.at(j).P4());
            mHs.push_back(p4sum.Mag());
        }

    std::pair<float, float> m_12_34 = make_pair (mHs.at(0), mHs.at(5));
    std::pair<float, float> m_13_24 = make_pair (mHs.at(1), mHs.at(4));
    std::pair<float, float> m_14_23 = make_pair (mHs.at(2), mHs.at(3));

    float r12_34 = sqrt (pow(m_12_34.first - targetmH, 2) + pow(m_12_34.second - targetmH, 2) );
    float r13_24 = sqrt (pow(m_13_24.first - targetmH, 2) + pow(m_13_24.second - targetmH, 2) );
    float r14_23 = sqrt (pow(m_14_23.first - targetmH, 2) + pow(m_14_23.second - targetmH, 2) );

    float the_min = std::min({r12_34, r13_24, r14_23});

    std::array<Jet, 4> outputJets = presel_jets;

    if (the_min == r12_34){
        outputJets.at(0) = presel_jets.at(1 - 1);
        outputJets.at(1) = presel_jets.at(2 - 1);
        outputJets.at(2) = presel_jets.at(3 - 1);
        outputJets.at(3) = presel_jets.at(4 - 1);
    }

    else if (the_min == r13_24){
        outputJets.at(0) = presel_jets.at(1 - 1);
        outputJets.at(1) = presel_jets.at(3 - 1);
        outputJets.at(2) = presel_jets.at(2 - 1);
        outputJets.at(3) = presel_jets.at(4 - 1);
    }

    else if (the_min == r14_23){
        outputJets.at(0) = presel_jets.at(1 - 1);
        outputJets.at(1) = presel_jets.at(4 - 1);
        outputJets.at(2) = presel_jets.at(2 - 1);
        outputJets.at(3) = presel_jets.at(3 - 1);   
    }

    else
        cout << "** [WARNING] : bbbb_jets_idxs_BothClosestToMh : something went wrong with finding the smallest radius" << endl;

    return outputJets;
}

// make the pairs that maximize their dR separation
std::array<Jet, 4> OfflineProducerHelper::bbbb_jets_idxs_MostBackToBack(std::array<Jet, 4> presel_jets)
{
    std::pair<CompositeCandidate,CompositeCandidate> H1H2 = get_two_best_jet_pairs (
        presel_jets,
        std::function<float (std::pair<CompositeCandidate,CompositeCandidate>)> ([](pair<CompositeCandidate,CompositeCandidate> x)-> float {return x.first.P4().DeltaR(x.second.P4());}),
        false);

    std::array<Jet, 4> output_jets;
    output_jets.at(0) = static_cast<Jet&> (H1H2.first.getComponent1());
    output_jets.at(1) = static_cast<Jet&> (H1H2.first.getComponent2());
    output_jets.at(2) = static_cast<Jet&> (H1H2.second.getComponent1());
    output_jets.at(3) = static_cast<Jet&> (H1H2.second.getComponent2());

    return output_jets;
}


// std::vector<std::pair<CompositeCandidate, CompositeCandidate>> OfflineProducerHelper::make_jet_combinations (std::vector<Jet> jets)
// {
//     // make all possible jet pairs (CompositeCandidates)
//     std::vector<CompositeCandidate> pairs;
//     for (uint i = 0; i < jets.size(); ++i)
//         for (uint j = i+1; j < jets.size(); ++j)
//         {
//             CompositeCandidate cc (jets.at(i), jets.at(j));
//             pairs.push_back(cc);
//         }

//     // make all possible pairs of CompositeCandiates, ensuring that each jet is taken only once
//     std::vector<pair<CompositeCandidate, CompositeCandidate>> jet_combinations;
//     for (uint i = 0; i < pairs.size(); ++i)
//         for (uint j = i+1; j < pairs.size(); ++j)
//         {
//             auto& cc1 = pairs.at(i);
//             auto& cc2 = pairs.at(j);
//             if (cc1.sharesComponentWith(cc2)) continue;
//             jet_combinations.push_back(make_pair(cc1, cc2));
//         }

//     return jet_combinations;
// }

bool OfflineProducerHelper::select_gen_HH (NanoAODTree& nat, EventInfo& ei)
{
    bool all_ok = true;
    for (uint igp = 0; igp < *(nat.nGenPart); ++igp)
    {
        GenPart gp (igp, &nat);
        if (abs(get_property(gp, GenPart_pdgId)) != 25) continue;
        // bool isFirst = checkBit(get_property(gp, GenPart_statusFlags), 12); // 12: isFirstCopy
        // if (!isFirst) continue;
        if (gp.isFirstCopy())
        {
            if (!assign_to_uninit(gp, {&ei.gen_H1, &ei.gen_H2} )) {
                cout << "** [WARNING] : select_gen_HH : more than two Higgs found" << endl;
                all_ok = false;
            }
            
            // // assign
            // if      (!ei.gen_H1) ei.gen_H1 = gp;
            // else if (!ei.gen_H2) ei.gen_H2 = gp;
            // else{
            //     cout << "** [WARNING] : select_gen_HH : more than two Higgs found" << endl;
            //     return false;
            // }
        }
        if (gp.isLastCopy())
        {
            if (!assign_to_uninit(gp, {&ei.gen_H1_last, &ei.gen_H2_last} )) {
                cout << "** [WARNING] : select_gen_HH : more than two Higgs found (last copy)" << endl;
                all_ok = false;
            }

            // assign
            // if      (!ei.gen_H1_last) ei.gen_H1_last = gp;
            // else if (!ei.gen_H2_last) ei.gen_H2_last = gp;
            // else{
            //     cout << "** [WARNING] : select_gen_HH : more than two Higgs found (last copy)" << endl;
            //     return false;
            // }
        }
    }

    if (!ei.gen_H1 || !ei.gen_H2){
        cout << "** [WARNING] : select_gen_HH : didn't find two Higgs : "
             << std::boolalpha
             << "H1 :" << ei.gen_H1.is_initialized()
             << "H2 :" << ei.gen_H2.is_initialized()
             << std::noboolalpha
             << endl;
        all_ok = false;
    }
    return all_ok;
}

bool OfflineProducerHelper::select_gen_bb_bb (NanoAODTree& nat, EventInfo& ei)
{
    if (!ei.gen_H1 || !ei.gen_H2)
    {
        cout << "** [WARNING] : select_gen_bb_bb : you need to search for H1 and H2 before" << endl;
        return false;
    }

    bool all_ok = true;
    for (uint igp = 0; igp < *(nat.nGenPart); ++igp)
    {
        GenPart gp (igp, &nat);
        if (abs(get_property(gp, GenPart_pdgId)) != 5) continue;
        // bool isFirst = checkBit(get_property(gp, GenPart_statusFlags), 12); // 12: isFirstCopy
        // if (!isFirst) continue;
        if (!gp.isFirstCopy()) continue;
        
        int idxMoth = get_property(gp, GenPart_genPartIdxMother);
        if (idxMoth == ei.gen_H1_last->getIdx())
        {
            if (!assign_to_uninit(gp, {&ei.gen_H1_b1, &ei.gen_H1_b2} )) {
                cout << "** [WARNING] : select_gen_bb_bb : more than two b from H1 found" << endl;
                all_ok = false;                
            }
            // if      (!ei.gen_H1_b1) ei.gen_H1_b1 = gp;
            // else if (!ei.gen_H1_b2) ei.gen_H1_b2 = gp;
            // else{
            //     cout << "** [WARNING] : select_gen_bb_bb : more than two b from H1 found" << endl;
            //     all_ok = false;
            // }
        }
        else if (idxMoth == ei.gen_H2_last->getIdx())
        {
            if (!assign_to_uninit(gp, {&ei.gen_H2_b1, &ei.gen_H2_b2} )) {
                cout << "** [WARNING] : select_gen_bb_bb : more than two b from H2 found" << endl;
                all_ok = false;                
            }

            // if      (!ei.gen_H2_b1) ei.gen_H2_b1 = gp;
            // else if (!ei.gen_H2_b2) ei.gen_H2_b2 = gp;
            // else{
            //     cout << "** [WARNING] : select_gen_bb_bb : more than two b from H2 found" << endl;
            //     all_ok = false;
            // }            
        }
        else
        {
            // cout << "** [WARNING] : select_gen_bb_bb : found a b quark of idx " << gp.getIdx() << " and mother " << idxMoth
            //      << " that does not match last H1 H2 idx " << ei.gen_H1_last->getIdx() << " " << ei.gen_H2_last->getIdx()
            //      << endl;
            //      all_ok = false;
            // possibly something that has b --> b + g --> b + (other stuff)
        }
    }
    if (!all_ok)
    {
        cout << "** [DEBUG] : select_gen_bb_bb : something went wrong, dumping gen parts" << endl;
        dump_gen_part(nat, true);
    }
    return all_ok;
}

bool OfflineProducerHelper::select_gen_VBF_partons (NanoAODTree& nat, EventInfo& ei)
{

    bool all_ok = true;

    for (uint igp = 0; igp < *(nat.nGenPart); ++igp)
    {
        GenPart gp (igp, &nat);
        if (abs(get_property(gp, GenPart_pdgId)) >= 6) continue; // only light quarks + b
        // auto flags = get_property(gp, GenPart_statusFlags);
        // search for ingoing partons
        // if (!) continue; // 12: isFirstCopy
        // search for outgoing partons

        if (get_property(gp, GenPart_status) == 21) // incoming VBF parton
        {
            // if      (!ei.gen_q1_in) ei.gen_q1_in = gp;
            // else if (!ei.gen_q2_in) ei.gen_q2_in = gp;
            // else{
            //     cout << "** [WARNING] : select_gen_VBF_partons : more than two incoming partons found" << endl;
            //     all_ok = false;
            // }    

            if (!assign_to_uninit(gp, {&ei.gen_q1_in, &ei.gen_q2_in} )) {
                cout << "** [WARNING] : select_gen_VBF_partons : more than two incoming partons found" << endl;
                all_ok = false;
            }
        }

        else
        {
            // it appears that in the VBF sample, the first objects in the list are the VBF ingoing partons
            // and that the next ones are the outgoing partons, 1st copy
            // so checking that idxMother == 0 should be enough for the VBF sample in use (checked on 80X on 19/01/2018)
            
            // FIXME: how to ensure the correspondence in_1 --> out_1, in_2 --> out_2 ?
            // unfortunately there is no valid gen info associated
            if (get_property(gp, GenPart_genPartIdxMother) == 0)
            {
                // if      (!ei.gen_q1_out) ei.gen_q1_out = gp;
                // else if (!ei.gen_q2_out) ei.gen_q2_out = gp;
                // else {
                //     cout << "** [WARNING] : select_gen_VBF_partons : more than two outgoing partons found" << endl;
                //     all_ok = false;
                // }
                if (!assign_to_uninit(gp, {&ei.gen_q1_out, &ei.gen_q2_out} )) {
                    cout << "** [WARNING] : select_gen_VBF_partons : more than two outgoing partons found" << endl;
                    all_ok = false;
                }
            }

            // the last copy of the outgoing parton can be checked by verifying two conditions:
            // 1. is last copy
            // 2. its mother is one of the outgoing VBF partons at ME. As they are before in the gen part list, at this point they should already have been identified
            // 23/01/2018: commenting out the outgoing partons, as it is not clear if this information makes sense
            else
            {
                // bool good_last = gp.isLastCopy();
                // if (ei.gen_q1_out && ei.gen_q2_out)
                // {
                //     // auto moth_idx = get_property(gp, GenPart_genPartIdxMother);
                //     auto moth_idx = recursive_gen_mother_idx(gp); // make a recursive search, if intermediate copies where radiated
                //     good_last = good_last && ((moth_idx == ei.gen_q1_out->getIdx()) || (moth_idx == ei.gen_q2_out->getIdx()));
                // }
                // else
                //     good_last = false;
                
                // // now check and save
                // if (good_last)
                // {
                //     if      (!ei.gen_q1_out_lastcopy) ei.gen_q1_out_lastcopy = gp;
                //     else if (!ei.gen_q2_out_lastcopy) ei.gen_q2_out_lastcopy = gp;
                //     else {
                //         cout << "** [WARNING] : select_gen_VBF_partons : more than two outgoing partons (last copy) found" << endl;
                //         all_ok = false;
                //     }
                // }
            }
        }
    }

    // quality checks

    if (!ei.gen_q1_in || !ei.gen_q2_in){
    cout << "** [WARNING] : select_gen_VBF_partons : didn't find two incoming partons : "
         << std::boolalpha
         << " q1_in :" << ei.gen_q1_in.is_initialized()
         << " q2_in :" << ei.gen_q2_in.is_initialized()
         << std::noboolalpha
         << endl;
        all_ok = false;
    }
    
    if (!ei.gen_q1_out || !ei.gen_q2_out){
    cout << "** [WARNING] : select_gen_VBF_partons : didn't find two outgoing partons : "
         << std::boolalpha
         << " q1_out :" << ei.gen_q1_out.is_initialized()
         << " q2_out :" << ei.gen_q2_out.is_initialized()
         << std::noboolalpha
         << endl;
        all_ok = false;
    }    
    
    // if (!ei.gen_q1_out_lastcopy || !ei.gen_q2_out_lastcopy){
    // cout << "** [WARNING] : select_gen_VBF_partons : didn't find two outgoing partons (last copy) : "
    //      << std::boolalpha
    //      << " q1_out_lastcopy :" << ei.gen_q1_out_lastcopy.is_initialized()
    //      << " q2_out_lastcopy :" << ei.gen_q2_out_lastcopy.is_initialized()
    //      << std::noboolalpha
    //      << endl;
    //     all_ok = false;
    // }    

    // // print a debug
    // if (!all_ok)
    // {
    //     cout << "** [DEBUG] select_gen_VBF_partons : something went wrong, dumping gen parts" << endl;
    //     for (uint igp = 0; igp < *(nat.nGenPart); ++igp)
    //     {
    //         GenPart gp (igp, &nat);
    //         // if (abs(get_property(gp, GenPart_pdgId)) >= 6) continue; // only light quarks + b
    //         cout << boolalpha;
    //         cout << igp << " -- "
    //              << " pdgId: "      << setw(4)  << get_property(gp, GenPart_pdgId)
    //              << " pt: "         << setw(10) << get_property(gp, GenPart_pt)
    //              << " eta: "        << setw(10) << get_property(gp, GenPart_eta)
    //              << " phi: "        << setw(10) << get_property(gp, GenPart_phi)
    //              << " status: "     << setw(4)  << get_property(gp, GenPart_status)
    //              << " moth_idx: "   << setw(4)  << get_property(gp, GenPart_genPartIdxMother)
    //              << endl;
    //         cout << "    . Flags :" << endl
    //             << "       isPrompt :                           " << gp.isPrompt() << endl
    //             << "       isDecayedLeptonHadron :              " << gp.isDecayedLeptonHadron() << endl
    //             << "       isTauDecayProduct :                  " << gp.isTauDecayProduct() << endl
    //             << "       isPromptTauDecayProduct :            " << gp.isPromptTauDecayProduct() << endl
    //             << "       isDirectTauDecayProduct :            " << gp.isDirectTauDecayProduct() << endl
    //             << "       isDirectPromptTauDecayProduct :      " << gp.isDirectPromptTauDecayProduct() << endl
    //             << "       isDirectHadronDecayProduct :         " << gp.isDirectHadronDecayProduct() << endl
    //             << "       isHardProcess :                      " << gp.isHardProcess() << endl
    //             << "       fromHardProcess :                    " << gp.fromHardProcess() << endl
    //             << "       isHardProcessTauDecayProduct :       " << gp.isHardProcessTauDecayProduct() << endl
    //             << "       isDirectHardProcessTauDecayProduct : " << gp.isDirectHardProcessTauDecayProduct() << endl
    //             << "       fromHardProcessBeforeFSR :           " << gp.fromHardProcessBeforeFSR() << endl
    //             << "       isFirstCopy :                        " << gp.isFirstCopy() << endl
    //             << "       isLastCopy :                         " << gp.isLastCopy() << endl
    //             << "       isLastCopyBeforeFSR :                " << gp.isLastCopyBeforeFSR() << endl
    //             << endl;
    //         cout << noboolalpha;
    //     }
    //     cout << endl << endl;
    // } // end of debug block

    return all_ok;
}

void OfflineProducerHelper::dump_gen_part (NanoAODTree& nat, bool printFlags)
{
    // print a debug
    cout << "** [DEBUG] dumping gen parts of event " << *(nat.event) << endl;
    cout << boolalpha;
    for (uint igp = 0; igp < *(nat.nGenPart); ++igp)
    {
        GenPart gp (igp, &nat);
        // if (abs(get_property(gp, GenPart_pdgId)) >= 6) continue; // only light quarks + b
        cout << igp << " -- "
             << " pdgId: "      << setw(4)  << get_property(gp, GenPart_pdgId)
             << " pt: "         << setw(10) << get_property(gp, GenPart_pt)
             << " eta: "        << setw(10) << get_property(gp, GenPart_eta)
             << " phi: "        << setw(10) << get_property(gp, GenPart_phi)
             << " status: "     << setw(4)  << get_property(gp, GenPart_status)
             << " moth_idx: "   << setw(4)  << get_property(gp, GenPart_genPartIdxMother)
             << endl;
        if (printFlags)
        {
            cout << "    . Flags :" << endl
                << "       isPrompt :                           " << gp.isPrompt() << endl
                << "       isDecayedLeptonHadron :              " << gp.isDecayedLeptonHadron() << endl
                << "       isTauDecayProduct :                  " << gp.isTauDecayProduct() << endl
                << "       isPromptTauDecayProduct :            " << gp.isPromptTauDecayProduct() << endl
                << "       isDirectTauDecayProduct :            " << gp.isDirectTauDecayProduct() << endl
                << "       isDirectPromptTauDecayProduct :      " << gp.isDirectPromptTauDecayProduct() << endl
                << "       isDirectHadronDecayProduct :         " << gp.isDirectHadronDecayProduct() << endl
                << "       isHardProcess :                      " << gp.isHardProcess() << endl
                << "       fromHardProcess :                    " << gp.fromHardProcess() << endl
                << "       isHardProcessTauDecayProduct :       " << gp.isHardProcessTauDecayProduct() << endl
                << "       isDirectHardProcessTauDecayProduct : " << gp.isDirectHardProcessTauDecayProduct() << endl
                << "       fromHardProcessBeforeFSR :           " << gp.fromHardProcessBeforeFSR() << endl
                << "       isFirstCopy :                        " << gp.isFirstCopy() << endl
                << "       isLastCopy :                         " << gp.isLastCopy() << endl
                << "       isLastCopyBeforeFSR :                " << gp.isLastCopyBeforeFSR() << endl
                << endl;
        }
    }
    cout << noboolalpha;
    cout << endl << endl;
}

bool OfflineProducerHelper::checkBit (int number, int bitpos)
{
    bool res = number & (1 << bitpos);
    return res;
}

int OfflineProducerHelper::recursive_gen_mother_idx(const GenPart& gp, bool stop_at_moth_zero)
{
    int imoth = get_property(gp, GenPart_genPartIdxMother);
    if (imoth < 0)
        return gp.getIdx();
    if (imoth == 0 && stop_at_moth_zero)
        return gp.getIdx();
    GenPart gp_moth (imoth, gp.getNanoAODTree());
    return recursive_gen_mother_idx(gp_moth, stop_at_moth_zero);
}
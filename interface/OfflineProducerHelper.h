#ifndef OFFLINEPRODUCERHELPER_H
#define OFFLINEPRODUCERHELPER_H

/**
** class  : OfflineProducerHelper
** author : L. Cadamuro (UF)
** date   : 10/01/2018
** brief  : A namespace containing a set of helper functions to analyze the events and process the skims
** note   : the prototype of the general function accepts as first input the NanoAODTree 
**          and as the second a reference to EventInto.
**          It's a flexible way to let the function know what it should know on the event, and possibly set extra information
**/

#include "NanoAODTree.h"
#include "EventInfo.h"
#include "CompositeCandidate.h"
#include "Jet.h"
#include <array>
#include <utility>
#include <functional>
#include <initializer_list>
#include <boost/optional.hpp>

namespace OfflineProducerHelper {

    // build a collection of jets in the event. If a prese;_function is given, preselect the jets
    std::vector<Jet> make_jets(NanoAODTree& nat, const std::function<bool (Jet)>& presel_function = nullptr);
    // filter the jet collection by keeping only the elements that pass filter_function
    void filter_jets(std::vector<Jet>& jets, const std::function<bool (Jet)>& filter_function);

    // --- - --- - --- - --- - --- - --- - --- - --- - --- - --- - --- - --- - --- - --- - 

    enum class bbbbSelectionStrategy{
        kOneClosestToMh,
        kBothClosestToMh,
        kMostBackToBack,
    };

    // functions that act on the EventInfo
    bool select_bbbb_jets (NanoAODTree& nat, EventInfo& ei, bbbbSelectionStrategy strat);

    // functions to pair a preselected set of four jets. They all shuffle the input set of jets and return them as (H1_b1, H1_b2, H2_b1, H2_b2)
    //
    // pick the pair that is closest to the Higgs mass, the other two remaining jets form the other pair
    std::array<Jet, 4> bbbb_jets_idxs_OneClosestToMh(std::array<Jet, 4> presel_jets, float targetmH = 115.);
    // minimize the distance of both pairs from the (targetmH, targetmH) point in a 2D plane
    std::array<Jet, 4> bbbb_jets_idxs_BothClosestToMh(std::array<Jet, 4> presel_jets, float targetmH = 115.);
    // by pair that is most back-to-back
    std::array<Jet, 4> bbbb_jets_idxs_MostBackToBack(std::array<Jet, 4> presel_jets);

    // combines a collection of type C of jets (either std::vector or std::array) into a collection of H H possible combinations
    // (i.e. all possible H1 = (jA, jB) and H2 = (jC, jD) choices)
    template <typename C>
    std::vector<std::pair<CompositeCandidate, CompositeCandidate>> make_jet_combinations (C jets);

    // given a collection of jets, returns the best two pairs.
    // the choice is made passing a sort_par_func that is evaluated on pairs of CompositeCandidates cc1 and cc2
    // e.g. if sort_par_func(cc1, cc2) is cc1.DeltaR(cc2), then the two CompositeCandidates that are less back-to-back is selected
    // using get_smallest = false returns the largest element instead
    // template <typename T>
    // std::pair<CompositeCandidate,CompositeCandidate> get_two_best_jet_pairs (
    //     std::vector<Jet> jets,
    //     std::function<T (std::pair<CompositeCandidate,CompositeCandidate>)> sort_par_func,
    //     bool get_smallest = true);

    template <typename T_sortpar, typename T_coll>
    std::pair<CompositeCandidate,CompositeCandidate> get_two_best_jet_pairs (
        T_coll jets,
        std::function<T_sortpar (std::pair<CompositeCandidate,CompositeCandidate>)> sort_par_func,
        bool get_smallest = true);


    bool select_gen_HH    (NanoAODTree& nat, EventInfo& ei);
    bool select_gen_bb_bb (NanoAODTree& nat, EventInfo& ei);
    bool select_gen_VBF_partons (NanoAODTree& nat, EventInfo& ei);

    // other utilities
    void dump_gen_part (NanoAODTree& nat, bool printFlags = true);
    bool checkBit(int number, int bitpos);
    int  recursive_gen_mother_idx(const GenPart& gp, bool stop_at_moth_zero = true); // returns the uppermost ancestor in the gen particle chain
                                                                                    // if stop_at_moth_zero = true, it returns when a GenPart is found and his mother has idx 0
                                                                                    // instead of returning the GenPart # 0 itself
    
    // loops on targets, and assigns value to the first element of target that is found to be uninitialized
    // returns false if none could be assigned, else return true
    template <typename T>
    bool assign_to_uninit(T value, std::initializer_list<boost::optional<T>*> targets);

    // given to Candidates val1 and val2, order them by pt
    // the function returns true if the objects were swapped, or false if they were left unchanged
    // if max_first = true, pt (val1) > pt (val2), if max_first = false the opposite
    template <typename T>
    bool order_by_pT(T& val1, T& val2, bool max_first = true);
};

template <typename C>
std::vector<std::pair<CompositeCandidate, CompositeCandidate>> OfflineProducerHelper::make_jet_combinations (C jets)
{
    const static bool debug = false;
    if (debug)
    {
        std::cout << "DEBUG make_jet_combinations : my jets in input are (pT)" << std::endl;
        for (auto& j : jets)
            std::cout << j.P4().Pt() << std::endl;        
    }

    // make all possible jet pairs (CompositeCandidates)
    std::vector<CompositeCandidate> pairs;
    for (uint i = 0; i < jets.size(); ++i)
        for (uint j = i+1; j < jets.size(); ++j)
        {
            CompositeCandidate cc (jets.at(i), jets.at(j));
            pairs.push_back(cc);
        }

    if (debug)
    {
        std::cout << "DEBUG make_jet_combinations : my pairs are (pT)" << std::endl;
        for (auto& j : pairs)
            std::cout << "(" << j.getComponent1().P4().Pt() << ", " << j.getComponent2().P4().Pt() << ")" << std::endl;        
    }

    // make all possible pairs of CompositeCandiates, ensuring that each jet is taken only once
    std::vector<std::pair<CompositeCandidate, CompositeCandidate>> jet_combinations;
    for (uint i = 0; i < pairs.size(); ++i)
        for (uint j = i+1; j < pairs.size(); ++j)
        {
            auto& cc1 = pairs.at(i);
            auto& cc2 = pairs.at(j);
            if (cc1.sharesComponentWith(cc2)) continue;
            jet_combinations.push_back(std::make_pair(cc1, cc2));
        }

    if (debug)
    {
        std::cout << "DEBUG : my combinations are (pT, pT)" << std::endl;
        for (auto& p : jet_combinations)
            std::cout << "(" << p.first.getComponent1().P4().Pt()  << ", " << p.first.getComponent2().P4().Pt()  << ") . "
                      << "(" << p.second.getComponent1().P4().Pt() << ", " << p.second.getComponent2().P4().Pt() << ")"
                      << std::endl;        
    }


    return jet_combinations;
}


// the function is templated on the return type of the sort_par_func, which must be sortable.
// all jet pairs (== Higgs cand) are made, and they are grouped into exclusive pairs (H1 and H2, not made from the same jet)
// e.g. 5 jets : A B C D E
// --> build: AB AC AD AE BC BD BE CD CE DE
// --> build all possible H1 H2 cands:
// AB CD
// AB DE
// AC BD
// etc...
// --> compute the sort parameter for each one, and sort it keeping track of the original (H1, H2) index
// --> return the best (H1, H2) pair
template <typename T_sortpar, typename T_coll>
std::pair<CompositeCandidate,CompositeCandidate> OfflineProducerHelper::get_two_best_jet_pairs (T_coll jets, std::function<T_sortpar (std::pair<CompositeCandidate,CompositeCandidate>)> sort_par_func, bool get_smallest)
{
    const static bool debug = false;

    if (debug)
    {
        std::cout << "DEBUG get_two_best_jet_pairs : my jets in input are (pT)" << std::endl;
        for (auto& j : jets)
            std::cout << j.P4().Pt() << std::endl;        
    }

    // make all (H1, H2) candidates
    auto combs = make_jet_combinations (jets);

    if (debug)
    {
        std::cout << "DEBUG : my combinations are (pT, pT)" << std::endl;
        for (auto& p : combs)
            std::cout << "(" << p.first.getComponent1().P4().Pt()  << ", " << p.first.getComponent2().P4().Pt()  << ") . "
                      << "(" << p.second.getComponent1().P4().Pt() << ", " << p.second.getComponent2().P4().Pt() << ")"
                      << std::endl;        
    }

    // compute sort parameter
    std::vector<std::pair<T_sortpar,int>> sort_pars;
    for (unsigned int idx = 0; idx < combs.size(); ++idx)
    {
        auto& jet = combs.at(idx);
        T_sortpar sort_par = sort_par_func(jet);
        sort_pars.push_back(std::make_pair(sort_par, idx));
    }

    // make the sort and get best two results
    std::stable_sort(sort_pars.begin(), sort_pars.end());
    if (get_smallest)
        return combs.at(sort_pars.begin()->second);
    
    else
        return combs.at(sort_pars.rbegin()->second);
}

template <typename T>
bool OfflineProducerHelper::assign_to_uninit(T value, std::initializer_list<boost::optional<T>*> targets)
{
    for (boost::optional<T>* tar : targets)
    {
        if (!(*tar))
        {
            *tar = value;
            return true;
        }
    }
    return false;
}

template <typename T>
bool OfflineProducerHelper::order_by_pT(T& val1, T& val2, bool max_first)
{
    float pt1 = val1.P4().Pt();
    float pt2 = val2.P4().Pt();

    if (max_first && pt1 > pt2)
        return false;
    if (!max_first && pt1 <= pt2)
        return false;
    std::swap(val1, val2);
    return true;
}

#endif
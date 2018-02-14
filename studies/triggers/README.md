# trigger

## General

Some utilities to :
1. filter the triggers and get the most efficient ones
2. verify if the triggers were prescaled

*NOTE*: some useful links to get trigger information are:
* the summary of the menu: http://fwyzard.web.cern.ch/fwyzard/hlt/2016/summary (can replace 2016 with any year)
* code and tools to access the prescales: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HLTStandAlonePrescaleInformation
* google doc with 2016 menu v4.2 V6 https://docs.google.com/spreadsheets/d/1bII_92pCrgk20A9FMIIHsOsYYj3f-lLjsoRkP_ZNQW4/edit#gid=1040120160
* google doc with the full 2017 menu (gives the prescaled/uprescaled status of each path) - only for 2017 menu :
https://docs.google.com/spreadsheets/d/1SqSKATiHcZdSB-t7rMb0SRSUL4Ot1jr2TExjnMrDUpA/edit#gid=0

Some quicker info is also found from WBM: https://cmswbm.cern.ch

In particular:
* find a fill --> find a run that is representative (e.g. 2016H, long fill, high PU --> run 283408)
* check the HLT paths and prescales (on RunSummary --> trigger mode)
* check the L1 seeds and prescales (on RunSummary --> other services : prescale set)

## The code in this folder
* filterTriggers.py : select trigger paths by min efficiency and check the unprescaled ones
* parseEfficiencies.py : get a list of trg path and eff from efficiencies.txt
* parsePrescalePage.py : parse the HTML file (to copy from internet with "view source" of the WBM page). It's a quick and dirty code to parse:
  * the HLT prescale columns page
  * the L1 prescale columns page
  * the L1 initial/final prescale in the run page (from RunSummary --> L1 key)
  * the google doc of the menu to get the prescaled/unprescaled status of a path
* the HTML files for run 283408 and the 2017 gdoc

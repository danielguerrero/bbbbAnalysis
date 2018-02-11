#!/bin/bash

#####Calculate all HLT trigger efficiencies (OUTPUT: efficiencies.txt)

#./bin/trigger_study.exe --cfg config/skim.cfg --input inputFiles/Samples_80X/VBF_HH_4b_10gen2018.txt --maxEvts 10000
#mv efficiencies.txt /studies/triggers

#####Filter all unprescaled trigger (OUTPUT: efficiencies_unprescaled.txt)

cd studies/triggers
python filterTriggers.py 
cd ..
cd ..
mv studies/triggers/efficiencies_unprescaled.txt .


#####Optimize trigger efficiencies (OUTPUT: efficiencies_report.txt)

./bin/trigger_optimization.exe --cfg config/skim.cfg --input inputFiles/Samples_80X/VBF_HH_4b_10gen2018.txt --maxEvts 10000
rm efficiencies_2.txt efficiencies_3.txt efficiencies_4.txt
#mv efficiencies_report.txt studies/triggers

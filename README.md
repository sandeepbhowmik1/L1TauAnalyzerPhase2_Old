# L1TauAnalyzerPhase2_Old


cmsrel CMSSW_10_6_1_patch2

cd CMSSW_10_6_1_patch2/src

cmsenv

git cms-init


# To get L1 HPS Tau


git remote add sandeepbhowmik1 https://github.com/sandeepbhowmik1/cmssw

git fetch sandeepbhowmik1 from-CMSSW_10_6_1_patch2

git cms-merge-topic -u sandeepbhowmik1:Phase2_L1HPSPFTau_v8

scram b -j 8



# To analyze L1 Tau


git clone https://github.com/sandeepbhowmik1/L1TauAnalyzerPhase2_Old $CMSSW_BASE/src/L1TauAnalyzerPhase2


scram b -j 8


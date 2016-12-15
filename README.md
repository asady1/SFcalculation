Instructions for running SFs for Double B with Double Muon selection

1. To download:

cmsrel CMSSW_7_6_4

cd CMSSW_7_6_4/src

cmsenv

git clone git://github.com/kskovpen/CFIT.git

cd CFIT

make

git clone git://github.com/asady/SFcalculation.git

cd SFcalculation

make

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../

2. To calculate the +/- templates systematic, edit runSFtemplates.sh to run over the +/- templates - it’s recommended to edit test.cxx so that it only prints out the nominal SF value (https://github.com/asady1/SFcalculation/blob/master/test.cxx#L897) and make sure you have inclSYS=0 (https://github.com/asady1/SFcalculation/blob/master/test.cxx#L7). Note that every time you make a change to test.cxx you have to run make again.

3. To calculate the 5 template systematic, edit runSFcalc.sh to run over your nominal templates, and edit test.cxx so that it doesn’t combine templates (comment out https://github.com/asady1/SFcalculation/blob/master/test.cxx#L234-L245, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L275-L278, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L281-L290, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L522-L532, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L547-L550, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L553-L561). Make sure inclSYS is still turned off, and don’t forget to uncomment what you commented out after you calculate.

4. Once you have these values for all pt bins and all WP, arrange them as an array for each systematic so that each ith component of each array corresponds to a systematic for the same WP and pt bin, i.e. the 0th and input these arrays into test.cxx (https://github.com/asady1/SFcalculation/blob/master/test.cxx#L786-L850).

5. Now you run SFcalc.sh, making sure that it reads the appropriate ith component for each WP and pt bin, making sure inclSYS=1, and making sure it outputs the information you need to make a .tex file (https://github.com/asady1/SFcalculation/blob/master/test.cxx#L397, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L406, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L710, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L716, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L874-L882, https://github.com/asady1/SFcalculation/blob/master/test.cxx#L894-L895). This will give you the percentage error for each systematic, and all the information you need for the SF table.

6. Once you have the nominal SFs, the systematic error up and down, and the systematic + statistical error up and down, along with the pT binning, you can input all of this into SFComp.C to produce plots for each WP, and you can use your output from step 5 inside a .tex table to make tables of percentage systematic error for each WP, and a table of SF information for all the WP.

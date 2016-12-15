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

2. To calculate the +/- templates systematic, edit runSFtemplates.sh to run over the +/- templates - it’s recommended to edit test.cxx so that it only prints out the nominal SF value ()

3. To calculate the 5 template systematic, edit runSFcalc.sh to run over your nominal templates, and edit test.cxx so that it doesn’t combine templates (comment out ) don’t forget to uncomment what you commented out after you calculate

4. Once you have these values for all pt bins and all WP, arrange them as an array for each systematic so that each ith component of each array corresponds to a systematic for the same WP and pt bin, i.e. the 0th and input these arrays into test.cxx ()

5. Now you run SFcalc.sh, making sure that it reads the appropriate ith component for each WP and pt bin, and making sure it outputs the information you need to make a .tex file ()

6. Once you have the nominal SFs, the systematic error up and down, and the systematic + statistical error up and down, along with the pT binning, you can input all of this into SFComp.C to produce plots for each WP.

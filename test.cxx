#include "cfit.h"
#include <iostream>
#include <vector>
#include <sstream>

//whether to include systematics 
bool inclSYS = 1;
//number of points to consider for stat error
int nStat = 100;
bool doSFc = 0;
std::string var = "JP";

void getResults(CFIT::cfit *cf,float *par,float *err);
void process(CFIT::cfit *cf,float *par_notag,float *err_notag,float *par_tag,float *err_tag,std::string sysVar,std::string t1Name,std::string t2Name,std::string t3Name,std::string t4Name,std::string t5Name,int statVariation);
double errfMC(double v1,double ve1,double v2,double ve2);
double errfDATA_UNCOR(double v1,double ve1,double v2,double ve2);
std::string getPtHeader(std::string pt);

std::vector<int> run;

std::vector<float> par_notag_bfromg;
std::vector<float> par_notag_b;
std::vector<float> par_notag_cfromg;
std::vector<float> par_notag_c;
std::vector<float> par_notag_l;

std::vector<float> err_notag_bfromg;
std::vector<float> err_notag_b;
std::vector<float> err_notag_cfromg;
std::vector<float> err_notag_c;
std::vector<float> err_notag_l;

std::vector<float> par_tag_bfromg;
std::vector<float> par_tag_b;
std::vector<float> par_tag_cfromg;
std::vector<float> par_tag_c;
std::vector<float> par_tag_l;

std::vector<float> err_tag_bfromg;
std::vector<float> err_tag_b;
std::vector<float> err_tag_cfromg;
std::vector<float> err_tag_c;
std::vector<float> err_tag_l;

std::vector<float> chi2_notag;
std::vector<float> chi2_tag;

std::vector<int> ndof_notag;
std::vector<int> ndof_tag;

std::vector<float> ndata_notag;
std::vector<float> ndata_tag;

std::vector<float> nmc_bfromg_notag;
std::vector<float> nmc_bfromg_tag;
std::vector<float> nmc_b_notag;
std::vector<float> nmc_b_tag;
std::vector<float> nmc_cfromg_notag;
std::vector<float> nmc_cfromg_tag;
std::vector<float> nmc_c_notag;
std::vector<float> nmc_c_tag;
std::vector<float> nmc_l_notag;
std::vector<float> nmc_l_tag;

std::vector<float> nmc_bfromgErr_notag;
std::vector<float> nmc_bfromgErr_tag;
std::vector<float> nmc_bErr_notag;
std::vector<float> nmc_bErr_tag;
std::vector<float> nmc_cfromgErr_notag;
std::vector<float> nmc_cfromgErr_tag;
std::vector<float> nmc_cErr_notag;
std::vector<float> nmc_cErr_tag;
std::vector<float> nmc_lErr_notag;
std::vector<float> nmc_lErr_tag;

std::vector<float> nmc_notag;
std::vector<float> nmc_tag;

std::vector<float> fr_bfromg_notag;
std::vector<float> fr_bfromg_tag;
std::vector<float> fr_b_notag;
std::vector<float> fr_b_tag;
std::vector<float> fr_cfromg_notag;
std::vector<float> fr_cfromg_tag;
std::vector<float> fr_c_notag;
std::vector<float> fr_c_tag;
std::vector<float> fr_l_notag;
std::vector<float> fr_l_tag;

std::vector<float> effMC;
std::vector<float> effMCErr;
std::vector<float> effDATA;
std::vector<float> effDATAErr;
std::vector<float> SF;
std::vector<float> SFErr;

int main(int argc,char *argv[])
{
   if( argc < 6 )
     {
	std::cout << "Please specify input parameters" << std::endl;
	exit(1);
     }
   //input file
   std::string fname = std::string(argv[1]);
   //pt bin
   std::string pt = std::string(argv[2]);
   //WP
   std::string tag = std::string(argv[3]);
   //selection (double muon, double muon tight, single muon)
   std::string runName = std::string(argv[4]);
   //directory with output directories
   std::string outdir = std::string(argv[5]);
   //output directory
   std::string dirname = std::string(argv[6]);
   //number corresponding to systematic for +/- templates and 5 template syst
   std::string sfnoST = std::string(argv[7]);
   //int sfno = 0; 
   int sfno = std::stoi(sfnoST);

   std::string tagMod = "";
   if( tag == "DoubleBL" ) tagMod = "DoubleBL";
   if( tag == "DoubleBM1" ) tagMod = "DoubleBM1";
   if( tag == "DoubleBM2" ) tagMod = "DoubleBM2";
   if( tag == "DoubleBH" ) tagMod = "DoubleBT";
   
   std::string foutName = outdir+"/"+dirname+"/"+pt+"_"+tagMod+".csv";
   std::ofstream fout(foutName.c_str());

   fout << "ptMin,ptMax,tag,run,\
nDATA,nDATATag,nMCb,nMCbTag,nMCc,\
nMCcTag,nMCl,nMClTag,pb,pbTag,frb,\
frbTag,frc,frcTag,frl,frlTag,chi2,\
chi2Tag,effMC,effData,SF,\
effMCErr,effDataErr,SFErr,pbErr,pbTagErr,\
pc,pcErr,pcTag,pcTagErr,pl,plErr,plTag,plTagErr,ndof,ndofTag,p4,p4Err,p4Tag,p4ErrTag,p5,p5Err,p5Tag,p5ErrTag,nMC4,nMC4Tag,nMC5,nMC5Tag\
\n";
   
   float par_notag[100];
   float err_notag[100];

   float par_tag[100];
   float err_tag[100];

   int nSYS = 0;
//   if( inclSYS ) nSYS = 6;
   if( inclSYS ) {
     nSYS = 7;
     //ALICE: this prevents us from running into trouble with the JES for doubleBH pt 300-400 for double muon
     //     if( runName == "DoubleT" )nSYS = 6;
     // else nSYS = 7;
   }
   std::string sys_down[100];
   std::string sys_up[100];
   std::string sys[100];

   //systematics read from input file
   if( inclSYS )
     {	
       	 sys[0] = "JES";
	 sys[1] = "NTRACKS";
	 sys[2] = "BFRAG";
	 sys[3] = "CFRAG";
	 sys[4] = "CD";
	 sys[5] = "K0L";
	 sys[6] = "PU";
     }

   CFIT::cfit *cf = new CFIT::cfit("JP discriminator");
   cf->SetVerbose(0);
   cf->ProducePlots(1);
   cf->SetLegendHeader(getPtHeader(pt));
   cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
   cf->SetCovarianceMode(COV_MAX);
//   cf->SetOptimization(OPT_NOCORR);
   cf->SetMorphing(OPTMORPH_CUTOFF,0.5);
//   cf->SetMorphing(OPTMORPH_GEOMETRIC);
   
   cf->SetInputFile(fname.c_str());

   for(int is=0;is<nSYS;is++)
     {
	sys_down[is] = "_"+sys[is]+"down";
	sys_up[is] = "_"+sys[is]+"up";
	cf->AddSys(sys[is].c_str(),sys_down[is].c_str(),sys_up[is].c_str());
     }
   
   std::string mname = "matrix_"+pt;
   cf->SetMatrixName(mname.c_str());
   cf->SetMatrixOption("WRITE");

   std::string t1Name = "g #rightarrow b#bar{b}";
   std::string t2Name = "b";
   std::string t3Name = "g #rightarrow c#bar{c}";
   std::string t4Name = "c";
   std::string t5Name = "light";

   std::string hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt";
   cf->SetData(hname_data.c_str());   
   std::string hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+tag+"pass_"+pt+"_data_opt";
   cf->SetDataTag(hname_data_tag.c_str());
   std::string hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+tag+"fail_"+pt+"_data_opt";
   cf->SetDataUntag(hname_data_untag.c_str());
   
   std::string hname_bfromg = "UNWEIGHTED__QCD__FatJet_"+var+"_all_"+pt+"_bfromg_opt";
   std::string hname_b = "UNWEIGHTED__QCD__FatJet_"+var+"_all_"+pt+"_b_opt";
   std::string hname_cfromg = "UNWEIGHTED__QCD__FatJet_"+var+"_all_"+pt+"_cfromg_opt";
   std::string hname_c = "UNWEIGHTED__QCD__FatJet_"+var+"_all_"+pt+"_c_opt";
   std::string hname_l = "UNWEIGHTED__QCD__FatJet_"+var+"_all_"+pt+"_l_opt";

   std::string hname_b_cfromg = "UNWEIGHTED__QCD__FatJet_"+var+"_all_"+pt+"_b_cfromg_opt";
   std::string hname_c_l = "UNWEIGHTED__QCD__FatJet_"+var+"_all_"+pt+"_c_l_opt";
   //defining pretag templates
   if( runName == "Double" || runName == "" || runName == "Single" )
     {	
       //defining all 5 templates: gbb, gcc, b, c, l
	cf->AddTemplate(t1Name,hname_bfromg.c_str(),65);
	cf->AddTemplate(t2Name,hname_b.c_str(),213);
	cf->AddTemplate(t3Name,hname_cfromg.c_str(),208);
	cf->AddTemplate(t4Name,hname_c.c_str(),206);
	cf->AddTemplate(t5Name,hname_l.c_str(),212);
     }
   if( runName == "DoubleT" )
     {	
       //defining templates already added together gbb, gcc + b, c + l
       //used if one of the 5 separate template histos is empty
        t2Name = "b + g #rightarrow c#bar{c}";
	t3Name = "c + light";
	cf->AddTemplate(t1Name,hname_bfromg.c_str(),65);
	cf->AddTemplate(t2Name,hname_b_cfromg.c_str(),221);
	cf->AddTemplate(t3Name,hname_c_l.c_str(),93);
     }
   
   if( runName == "Double" )
     {	
       //adding some of the 5 templates so that we now only have three: gbb, gcc + b, c + l 
	std::vector<std::string> nameV;
	nameV.push_back(t2Name);
	nameV.push_back(t3Name);
      	cf->GlueTemplates(nameV,"b + g #rightarrow c#bar{c}",221);
	nameV.clear();
	nameV.push_back(t4Name);
	nameV.push_back(t5Name);
	cf->GlueTemplates(nameV,"c + light",93);
	}
    
   std::string hname_bfromg_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_bfromg_opt";
   std::string hname_b_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_b_opt";
   std::string hname_cfromg_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_cfromg_opt";
   std::string hname_c_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_c_opt";
   std::string hname_l_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_l_opt";

   std::string hname_b_cfromg_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_b_cfromg_opt";
   std::string hname_c_l_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_c_l_opt";
   std::string hname_b_cfromg_c_l_tag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"pass_"+pt+"_b_cfromg_c_l_opt";
   //defining post tag templates
     if( runName == "Double" || runName == "" || runName == "Single" )
     {
       //defining all 5 templates
	cf->AddTemplateTag(t1Name,hname_bfromg_tag.c_str(),65);
	cf->AddTemplateTag(t2Name,hname_b_tag.c_str(),213);
	cf->AddTemplateTag(t3Name,hname_cfromg_tag.c_str(),208);
	cf->AddTemplateTag(t4Name,hname_c_tag.c_str(),206);
	cf->AddTemplateTag(t5Name,hname_l_tag.c_str(),212);
     }
     if( runName == "DoubleT" )
     {	
       //defining the 3 condensed templates, adding them so that now we have gbb, gcc + b + l + c
        t2Name = "b + g #rightarrow c#bar{c}";
	t3Name = "c + light";
	cf->AddTemplateTag(t1Name,hname_bfromg_tag.c_str(),65);
	cf->AddTemplateTag(t2Name,hname_b_cfromg_tag.c_str(),221);
	cf->AddTemplateTag(t3Name,hname_c_l_tag.c_str(),93);
	
	std::vector<std::string> nameV;
	nameV.push_back(t2Name);
	nameV.push_back(t3Name);
	cf->GlueTemplatesTag(nameV,"other flavours",42);
	 }
     
     if( runName == "Double" )
     {	
       //adding templates so that now we have gbb, gcc + b + l + c 
	std::vector<std::string> nameV;
	nameV.push_back(t2Name);
	nameV.push_back(t3Name);
	nameV.push_back(t4Name);
	nameV.push_back(t5Name);
	cf->GlueTemplatesTag(nameV,"other flavours",42);
	}
     
      if( runName == "Single" )
     {
       if( tagMod == "DoubleBL" || tagMod == "DoubleBM1" )
	 {
	   std::vector<std::string> nameV;
	   nameV.push_back(t2Name);
	   nameV.push_back(t3Name);
	   cf->GlueTemplatesTag(nameV,"b + g #rightarrow c#bar{c}",221);
	   nameV.clear();
	   nameV.push_back(t4Name);
	   nameV.push_back(t5Name);
	   cf->GlueTemplatesTag(nameV,"c + light",93);
	 }
        if( tagMod == "DoubleBM2" || tagMod == "DoubleBT" )
	 {
	   std::vector<std::string> nameV;
	   nameV.push_back(t2Name);
	   nameV.push_back(t3Name);
	   nameV.push_back(t4Name);
	   nameV.push_back(t5Name);
	   cf->GlueTemplatesTag(nameV,"other flavours",42);
	   }
     }
   std::string hname_bfromg_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_bfromg_opt";
   std::string hname_b_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_b_opt";
   std::string hname_cfromg_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_cfromg_opt";
   std::string hname_c_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_c_opt";
   std::string hname_l_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_l_opt";

   std::string hname_b_cfromg_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_b_cfromg_opt";
   std::string hname_c_l_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_c_l_opt";
   std::string hname_b_cfromg_c_l_untag = "UNWEIGHTED__QCD__FatJet_"+var+"_"+tag+"fail_"+pt+"_b_cfromg_c_l_opt";
   //defining untag templates
   if( runName == "Double" || runName == "" || runName == "Single" )
     {	
	cf->AddTemplateUntag(t1Name,hname_bfromg_untag,7);
	cf->AddTemplateUntag(t2Name,hname_b_untag,213);
	cf->AddTemplateUntag(t3Name,hname_cfromg_untag,42);
	cf->AddTemplateUntag(t4Name,hname_c_untag,8);
	cf->AddTemplateUntag(t5Name,hname_l_untag,4);
     }   
   if( runName == "DoubleT" )
     {	
        t2Name = "b + g #rightarrow c#bar{c}";
      	t3Name = "c + light";
	cf->AddTemplateUntag(t1Name,hname_bfromg_untag.c_str(),65);
	cf->AddTemplateUntag(t2Name,hname_b_cfromg_untag.c_str(),221);
	cf->AddTemplateUntag(t3Name,hname_c_l_untag.c_str(),93);
     }   
   std::cout << "Run systematic variations" << std::endl;
   
   if( !inclSYS )
     {	
	process(cf,par_notag,err_notag,par_tag,err_tag,"",t1Name,t2Name,t3Name,t4Name,t5Name,-1);
	
	run.push_back(0);
	
	par_notag_bfromg.push_back(par_notag[0]);
	par_notag_b.push_back(par_notag[1]);
	par_notag_cfromg.push_back(par_notag[2]);
	par_notag_c.push_back(par_notag[3]);
	par_notag_l.push_back(par_notag[4]);
	
	err_notag_bfromg.push_back(err_notag[0]);
	err_notag_b.push_back(err_notag[1]);
	err_notag_cfromg.push_back(err_notag[2]);
	err_notag_c.push_back(err_notag[3]);
	err_notag_l.push_back(err_notag[4]);
	
	par_tag_bfromg.push_back(par_tag[0]);
	par_tag_b.push_back(par_tag[1]);
	par_tag_cfromg.push_back(par_tag[2]);
	par_tag_c.push_back(par_tag[3]);
	par_tag_l.push_back(par_tag[4]);
	
	err_tag_bfromg.push_back(err_tag[0]);
	err_tag_b.push_back(err_tag[1]);
	err_tag_cfromg.push_back(err_tag[2]);
	err_tag_c.push_back(err_tag[3]);
	err_tag_l.push_back(err_tag[4]);
     }
   else
     {
	// nosys and sys
	for(int is=0;is<=2*nSYS;is++)
	  {
	     // if( is != 0 ) continue;
	     //	if( is != 0 && is != 1 ) continue;
	     int iss = is;			
	     
	     std::string sysVar = "";

	     if( iss > nSYS )
	       {
		  iss = -(iss-nSYS);
	       }
	     if( iss != 0 )
	       {
		  cf->SetMatrixOption("READ");
		  
		  if( iss > 0 )
		    sysVar = sys_up[iss-1];
		  else
		    sysVar = sys_down[abs(iss)-1];
		  
		 std::cout << "Systematics variation " << sysVar << std::endl;
	       } 
	     //	    if( sysVar != "_NTRACKSup"){
	     if( iss != 0 ) cf->ProducePlots(0);
	     
	     process(cf,par_notag,err_notag,par_tag,err_tag,sysVar,t1Name,t2Name,t3Name,t4Name,t5Name,-1);

	     run.push_back(iss);

	     std::cout << iss << std::endl;
	     par_notag_bfromg.push_back(par_notag[0]);
	     par_notag_b.push_back(par_notag[1]);
	     par_notag_cfromg.push_back(par_notag[2]);
	     par_notag_c.push_back(par_notag[3]);
	     par_notag_l.push_back(par_notag[4]);
	     
	     err_notag_bfromg.push_back(err_notag[0]);
	     err_notag_b.push_back(err_notag[1]);
	     err_notag_cfromg.push_back(err_notag[2]);
	     err_notag_c.push_back(err_notag[3]);
	     err_notag_l.push_back(err_notag[4]);
	     
	     par_tag_bfromg.push_back(par_tag[0]);
	     par_tag_b.push_back(par_tag[1]);
	     par_tag_cfromg.push_back(par_tag[2]);
	     par_tag_c.push_back(par_tag[3]);
	     par_tag_l.push_back(par_tag[4]);
	
	     err_tag_bfromg.push_back(err_tag[0]);
	     err_tag_b.push_back(err_tag[1]);
	     err_tag_cfromg.push_back(err_tag[2]);
	     err_tag_c.push_back(err_tag[3]);
	     err_tag_l.push_back(err_tag[4]);
	     }
	  }
   //        }
   
   //std::cout << "Run statistical variations" << std::endl;
   
   for(int is=0;is<nStat;is++)
     {	
	int isys = 666+is;
	cf->SetMatrixOption("READ");
	cf->ProducePlots(0);

	//	std::cout << "Statistics variation " << isys << std::endl;
	
	process(cf,par_notag,err_notag,par_tag,err_tag,"",t1Name,t2Name,t3Name,t4Name,t5Name,isys);

	run.push_back(isys);
	
	par_notag_bfromg.push_back(par_notag[0]);
	par_notag_b.push_back(par_notag[1]);
	par_notag_cfromg.push_back(par_notag[2]);
	par_notag_c.push_back(par_notag[3]);
	par_notag_l.push_back(par_notag[4]);
	
	err_notag_bfromg.push_back(err_notag[0]);
	err_notag_b.push_back(err_notag[1]);
	err_notag_cfromg.push_back(err_notag[2]);
	err_notag_c.push_back(err_notag[3]);
	err_notag_l.push_back(err_notag[4]);
	
	par_tag_bfromg.push_back(par_tag[0]);
	par_tag_b.push_back(par_tag[1]);
	par_tag_cfromg.push_back(par_tag[2]);
	par_tag_c.push_back(par_tag[3]);
	par_tag_l.push_back(par_tag[4]);
	
	err_tag_bfromg.push_back(err_tag[0]);
	err_tag_b.push_back(err_tag[1]);
	err_tag_cfromg.push_back(err_tag[2]);
	err_tag_c.push_back(err_tag[3]);
	err_tag_l.push_back(err_tag[4]);
     }
   
   // additional systematics variations
   for(int i=100;i<=104;i++)
     {
	// recalculate correlation matrix according to optimisation settings
	cf->SetMatrixOption("WRITE");
	cf->ProducePlots(0);
	
	if( i == 100 ) cf->SetOptimization(OPT_NOCORR);
	if( i == 101 ) cf->SetMorphing(OPTMORPH_CUTOFF,0.25);
	if( i == 102 ) cf->SetMorphing(OPTMORPH_CUTOFF,0.75);
//	if( i == 103 ) cf->SetMorphing(OPTMORPH_CUTOFF,0.5);
	if( i == 103 ) cf->SetMorphing(OPTMORPH_GEOMETRIC);
	if( i == 104 )
	  {
	     delete cf;
	     cf = new CFIT::cfit();
	     
	     cf->SetVerbose(0);
	     cf->ProducePlots(0);
	     cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
	     cf->SetCovarianceMode(COV_MAX);
	     cf->SetMorphing(OPTMORPH_CUTOFF,0.5);
//	     cf->SetMorphing(OPTMORPH_GEOMETRIC);
	     
	     cf->SetInputFile(fname.c_str());
   
	     std::string mname = "matrix_"+pt;
	     cf->SetMatrixName(mname.c_str());
	     cf->SetMatrixOption("WRITE");

	     cf->SetData(hname_data.c_str());   
	     cf->SetDataTag(hname_data_tag.c_str());
	     cf->SetDataUntag(hname_data_untag.c_str());
	     //defining pretag, tag, untag: copy what's already done above
	     if( runName == "Double" || runName == "" || runName == "Single" )
	       {		  
		  cf->AddTemplate(t1Name,hname_bfromg.c_str(),65);
		  cf->AddTemplate(t2Name,hname_b.c_str(),213);
		  cf->AddTemplate(t3Name,hname_cfromg.c_str(),208);
		  cf->AddTemplate(t4Name,hname_c.c_str(),206);
		  cf->AddTemplate(t5Name,hname_l.c_str(),212);
	       }	     
	     if( runName == "DoubleT" )
	       {		  
		  cf->AddTemplate(t1Name,hname_bfromg.c_str(),65);
		  cf->AddTemplate(t2Name,hname_b_cfromg.c_str(),221);
		  cf->AddTemplate(t3Name,hname_c_l.c_str(),93);
	       }	     
	     
	     if( runName == "Double" )
	       {	
		  std::vector<std::string> nameV;
		  nameV.push_back(t2Name);
		  nameV.push_back(t3Name);
		  cf->GlueTemplates(nameV,"b + g #rightarrow c#bar{c}",221);
		  nameV.clear();
		  nameV.push_back(t4Name);
		  nameV.push_back(t5Name);
		  cf->GlueTemplates(nameV,"c + light",93);
		  }
      	     if( runName == "Double" || runName == "" || runName == "Single" )
	       {		  
		  cf->AddTemplateTag(t1Name,hname_bfromg_tag.c_str(),65);
		  cf->AddTemplateTag(t2Name,hname_b_tag.c_str(),213);
		  cf->AddTemplateTag(t3Name,hname_cfromg_tag.c_str(),208);
		  cf->AddTemplateTag(t4Name,hname_c_tag.c_str(),206);
		  cf->AddTemplateTag(t5Name,hname_l_tag.c_str(),212);
	       }	     
	     if( runName == "DoubleT" )
	       {		  
		  cf->AddTemplateTag(t1Name,hname_bfromg_tag.c_str(),65);
		  cf->AddTemplateTag(t2Name,hname_b_cfromg_tag.c_str(),221);
		  cf->AddTemplateTag(t3Name,hname_c_l_tag.c_str(),93);
		  
		  std::vector<std::string> nameV;
		  nameV.push_back(t2Name);
		  nameV.push_back(t3Name);
		  cf->GlueTemplatesTag(nameV,"other flavours",42);
		  }
	     
	     if( runName == "Double" )
	       {	
		  std::vector<std::string> nameV;
		  nameV.push_back(t2Name);
		  nameV.push_back(t3Name);
		  nameV.push_back(t4Name);
		  nameV.push_back(t5Name);
		  cf->GlueTemplatesTag(nameV,"other flavours",42);
		  }
	     
	     if( runName == "Single" )
	       {
		 if( tagMod == "DoubleBL" || tagMod == "DoubleBM1" )
		   {
		     std::vector<std::string> nameV;
		     nameV.push_back(t2Name);
		     nameV.push_back(t3Name);
		     cf->GlueTemplatesTag(nameV,"b + g #rightarrow c#bar{c}",221);
		     nameV.clear();
		     nameV.push_back(t4Name);
		     nameV.push_back(t5Name);
		     cf->GlueTemplatesTag(nameV,"c + light",93);
		   }
		 if( tagMod == "DoubleBM2" || tagMod == "DoubleBT" )
		   {
		     std::vector<std::string> nameV;
		     nameV.push_back(t2Name);
		     nameV.push_back(t3Name);
		     nameV.push_back(t4Name);
		     nameV.push_back(t5Name);
		     cf->GlueTemplatesTag(nameV,"other flavours",42);
		     }
	       }
	     
	     if( runName == "Double" || runName == "" || runName == "Single" )
	       {		  
		  cf->AddTemplateUntag(t1Name,hname_bfromg_untag,7);
		  cf->AddTemplateUntag(t2Name,hname_b_untag,60);
		  cf->AddTemplateUntag(t3Name,hname_cfromg_untag,42);
		  cf->AddTemplateUntag(t4Name,hname_c_untag,8);
		  cf->AddTemplateUntag(t5Name,hname_l_untag,4);
	       }	     
	     if( runName == "DoubleT" )
	       {		  
		  cf->AddTemplateUntag(t1Name,hname_bfromg_untag.c_str(),65);
		  cf->AddTemplateUntag(t2Name,hname_b_cfromg_untag.c_str(),221);
		  cf->AddTemplateUntag(t3Name,hname_c_l_untag.c_str(),93);		  
	       }
	  }

	process(cf,par_notag,err_notag,par_tag,err_tag,"",t1Name,t2Name,t3Name,t4Name,t5Name,-1);
	
	run.push_back(i);
	
	par_notag_bfromg.push_back(par_notag[0]);
	par_notag_b.push_back(par_notag[1]);
	par_notag_cfromg.push_back(par_notag[2]);
	par_notag_c.push_back(par_notag[3]);
	par_notag_l.push_back(par_notag[4]);
	
	err_notag_bfromg.push_back(err_notag[0]);
	err_notag_b.push_back(err_notag[1]);
	err_notag_cfromg.push_back(err_notag[2]);
	err_notag_c.push_back(err_notag[3]);
	err_notag_l.push_back(err_notag[4]);
	
	par_tag_bfromg.push_back(par_tag[0]);
	par_tag_b.push_back(par_tag[1]);
	par_tag_cfromg.push_back(par_tag[2]);
	par_tag_c.push_back(par_tag[3]);
	par_tag_l.push_back(par_tag[4]);
	
	err_tag_bfromg.push_back(err_tag[0]);
	err_tag_b.push_back(err_tag[1]);
	err_tag_cfromg.push_back(err_tag[2]);
	err_tag_c.push_back(err_tag[3]);
	err_tag_l.push_back(err_tag[4]);
	
	// restore initial settings
	cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
	cf->SetMorphing(OPTMORPH_CUTOFF,0.5);
//	cf->SetMorphing(OPTMORPH_GEOMETRIC);
     }
      
   std::stringstream ss(pt);
   std::string item;
   int idx = 0;
   std::string ptMin = "";
   std::string ptMax = "";
//   while( std::getline(ss,item,'t') )
//     {
//	if( idx == 0 ) ptMin = item;
//	if( idx == 2 ) ptMax = item;
//	idx++;
//     }
   //defining pT regions
   if( pt == "pt5" ) {ptMin = "300"; ptMax = "350";}
   else if( pt == "pt4" ) {ptMin = "400"; ptMax = "500";}
   else if( pt == "pt2" ) {ptMin = "500"; ptMax = "600";}
   else if( pt == "pt3" ) {ptMin = "600"; ptMax = "700";}
   else if( pt == "pt0" ) {ptMin = "400"; ptMax = "450";}
   else if( pt == "pt1" ) {ptMin = "450"; ptMax = "500";}
   else if( pt == "pt6" ) {ptMin = "500"; ptMax = "700";}
   else if( pt == "pt7" ) {ptMin = "250"; ptMax = "300";}
   else if( pt == "ptall" ) {ptMin = "300"; ptMax = "700";}
   
   std::string dl = ",";
   int nRun = SF.size();
   //ALICE: initializing variables needed to calculate stat and syst errors
   float sf_stat_sum = 0;
   float sigma_syst_up = 0;
   float sigma_syst_down = 0;

   for(int i=0;i<nRun;i++)
     {
	float nmc = 0;
	float nmc_tag = 0;
	
	nmc = (nmc_bfromg_notag[i]+nmc_b_notag[i]+nmc_cfromg_notag[i]+nmc_c_notag[i]+nmc_l_notag[i]);	
	nmc_tag = (nmc_bfromg_tag[i]+nmc_b_tag[i]+nmc_cfromg_tag[i]+nmc_c_tag[i]+nmc_l_tag[i]);

	float fr_bfromg_notag_new = 0;
	float fr_b_notag_new = 0;
	float fr_cfromg_notag_new = 0;
	float fr_c_notag_new = 0;
	float fr_l_notag_new = 0;
	
	fr_bfromg_notag_new = nmc/ndata_notag[i]*par_notag_bfromg[i]*fr_bfromg_notag[i];
	fr_b_notag_new = nmc/ndata_notag[i]*par_notag_b[i]*fr_b_notag[i];
	fr_cfromg_notag_new = nmc/ndata_notag[i]*par_notag_cfromg[i]*fr_cfromg_notag[i];
	fr_c_notag_new = nmc/ndata_notag[i]*par_notag_c[i]*fr_c_notag[i];
	fr_l_notag_new = nmc/ndata_notag[i]*par_notag_l[i]*fr_l_notag[i];

	float fr_bfromg_tag_new = 0;
	float fr_b_tag_new = 0;
	float fr_cfromg_tag_new = 0;
	float fr_c_tag_new = 0;
	float fr_l_tag_new = 0;

	fr_bfromg_tag_new = nmc_tag/ndata_tag[i]*par_tag_bfromg[i]*fr_bfromg_tag[i];
	fr_b_tag_new = nmc_tag/ndata_tag[i]*par_tag_b[i]*fr_b_tag[i];
	fr_cfromg_tag_new = nmc_tag/ndata_tag[i]*par_tag_cfromg[i]*fr_cfromg_tag[i];
	fr_c_tag_new = nmc_tag/ndata_tag[i]*par_tag_c[i]*fr_c_tag[i];
	fr_l_tag_new = nmc_tag/ndata_tag[i]*par_tag_l[i]*fr_l_tag[i];

	std::string strOut = "";
	//calculating syst errors up and down, here we add all of them together, taken from combine.C
	if ( run[i] < 666 ){
	  float delta = SF[0] - SF[i];
	  // std::cout << SF[i] << std::endl;
	  bool isUp =(delta < 0);
	  float err = delta*delta;
	  // std::cout << "run " << run[i] << " : " << err << std::endl;
	  if ( run[i] < 100 ){
	    if( isUp ) sigma_syst_up += err;
	    else sigma_syst_down += err;
	    //print out percentage for sys
	    std::cout << run[i] << " & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor((delta/SF[0])*1000+0.5)/1000 << std::endl ;
	   }
	  if ( run[i] >= 100 ){
	    sigma_syst_up += pow(delta, 2);
	    sigma_syst_down += pow(delta, 2);
	    //print out percentage for sys 
	    std::cout << run[i] << " & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor((delta/SF[0])*1000+0.5)/1000 << std::endl ;
	   }
	}
	//calculating average SF for stat errors
	else{
	  sf_stat_sum += SF[i];
	}
	strOut = ptMin+dl+ptMax+dl+tag+dl+
	  std::to_string(run[i])+dl+
	  std::to_string(ndata_notag[i])+dl+
	  std::to_string(ndata_tag[i])+dl+
	  std::to_string(nmc_bfromg_notag[i])+dl+
	  std::to_string(nmc_bfromg_tag[i])+dl+
	  std::to_string(nmc_b_notag[i])+dl+
	  std::to_string(nmc_b_tag[i])+dl+
	  std::to_string(nmc_cfromg_notag[i])+dl+
	  std::to_string(nmc_cfromg_tag[i])+dl+
	  std::to_string(par_notag_bfromg[i])+dl+
	  std::to_string(par_tag_bfromg[i])+dl+
	  std::to_string(fr_bfromg_notag_new)+dl+
	  std::to_string(fr_bfromg_tag_new)+dl+
	  std::to_string(fr_b_notag_new)+dl+
	  std::to_string(fr_b_tag_new)+dl+
	  std::to_string(fr_cfromg_notag_new)+dl+
	  std::to_string(fr_cfromg_tag_new)+dl+
	  std::to_string(chi2_notag[i])+dl+
	  std::to_string(chi2_tag[i])+dl+
	  std::to_string(effMC[i])+dl+
	  std::to_string(effDATA[i])+dl+
	  std::to_string(SF[i])+dl+
	  std::to_string(effMCErr[i])+dl+
	  std::to_string(effDATAErr[i])+dl+
	  std::to_string(SFErr[i])+dl+
	  std::to_string(err_notag_bfromg[i])+dl+
	  std::to_string(err_tag_bfromg[i])+dl+
	  std::to_string(par_notag_b[i])+dl+
	  std::to_string(err_notag_b[i])+dl+
	  std::to_string(par_tag_b[i])+dl+
	  std::to_string(err_tag_b[i])+dl+
	  std::to_string(par_notag_cfromg[i])+dl+
	  std::to_string(err_notag_cfromg[i])+dl+
	  std::to_string(par_tag_cfromg[i])+dl+
	  std::to_string(err_tag_cfromg[i])+dl+
	  std::to_string(ndof_notag[i])+dl+
	  std::to_string(ndof_tag[i])+dl+
	  std::to_string(par_notag_c[i])+dl+
	  std::to_string(err_notag_c[i])+dl+
	  std::to_string(par_tag_c[i])+dl+
	  std::to_string(err_tag_c[i])+dl+
	  std::to_string(par_notag_l[i])+dl+
	  std::to_string(err_notag_l[i])+dl+
	  std::to_string(par_tag_l[i])+dl+
	  std::to_string(err_tag_l[i])+dl+
	  std::to_string(nmc_c_notag[i])+dl+
	  std::to_string(nmc_c_tag[i])+dl+
	  std::to_string(nmc_l_notag[i])+dl+
	  std::to_string(nmc_l_tag[i])+
	  "\n";
	fout << strOut.c_str();
	}   
   //calculating average SF for stat errors
   sf_stat_sum /= 100;
   float sigma_stat = 0;
   //calculating stat errors
   for (int i=0;i<nRun;i++){
     if( run[i] > 600 ) sigma_stat += pow(sf_stat_sum - SF[i],2)/100;
   }
   sigma_stat = sqrt(sigma_stat); 
  

       //hard coded up and down for all four background templates to add systematic for these - once calculated, the sfno fed to the command line should correspond to the i in the array so that the correct SF_sys is paired with the correct SF_nom
   /*   float SF_b_up[16] = {0.91, 0.76, 0.9, 0.91, 0.81, 0.84, 1.03, 0.8, 0.71, 0.95, 1.05, 0.69, 0.94, 0.86, 0.87, 1.09};
   float SF_b_down[16] = {0.9, 0.78, 0.89, 0.9, 0.8, 0.83, 1.02, 0.84, 0.72, 0.97, 1.07, 0.69, 0.94, 0.85, 0.87, 1.14};
   float SF_c_up[16] = {0.91, 0.77, 0.9, 0.91, 0.8, 0.84, 1.03, 0.81, 0.72, 0.96, 1.05, 0.69, 0.94, 0.86, 0.87, 1.13};
   float SF_c_down[16] = {0.9, 0.77, 0.9, 0.9, 0.81, 0.83, 1.02, 0.83, 0.71, 0.95, 1.06, 0.69, 0.94, 0.86, 0.87, 1.1};
   float SF_l_up[16] = {0.91, 0.79, 0.9, 0.91, 0.81, 0.84, 1.02, SF[0], 0.72, SF[0], SF[0], SF[0], 0.93, 0.85, SF[0], SF[0]};
   float SF_l_down[16] = {0.9, 0.76, 0.89, 0.91, 0.81, 0.83, 1.03, SF[0], 0.71, SF[0], SF[0], SF[0], 0.94, 0.86, SF[0],SF[0]};
   float SF_cfg_up[16] = {0.91, 0.76, 0.89, 0.9, 0.8, 0.84, 1.02, SF[0], 0.71, SF[0], SF[0], SF[0],0.93, 0.85, SF[0], SF[0]};
   float SF_cfg_down[16] = {0.9, 0.81, 0.92, 0.91, 0.85, 0.86, 1.03, SF[0], 0.73, SF[0], SF[0], SF[0], 0.94, 0.86, SF[0], SF[0]};
   float SF_5_temp[16] = {0.92, 0.73, 0.9, 0.91, 0.81, 0.83, 1.05, 0.81, 0.72, SF[0], 1.05, 0.75, 0.95, 0.85, 0.87, 1.08};*/
   
   /*   float SF_b_up[8] = {0.89, 0.95, 0.83, 0.89, 0.83, 0.92, 1.08, 0.68};
   float SF_b_down[8] = {0.91, 0.97, 0.88, 0.96, 0.83, 0.93, 1.07, 0.71};
   float SF_c_up[8] = {0.89, 0.96, 0.84, 0.94, 0.83, 0.92, 1.07, 0.68};
   float SF_c_down[8] = {0.9, 0.96, 0.86, 0.95, 0.83, 0.93, 1.07, 0.71};
   float SF_l_up[8] = {0.89, 0.96, 0.85, 0.94, 0.83, 0.93, 1.08, 0.7};
   float SF_l_down[8] = {0.9, 0.96, 0.86, 0.95, 0.83, 0.93, 1.07, 0.69};
   float SF_cfg_up[8] = {0.93, 0.99, 0.92, 0.95, 0.84, 0.95, 1.07, 0.74};
   float SF_cfg_down[8] = {0.79, 0.95, 0.75, 0.93, 0.75, 0.92, 0.96, 0.68};
   float SF_5_temp[8] = {0.83, 1.02, 0.76, 0.94, 0.92, 0.91, 1.04, 0.67};
   //std::cout << SF_b_up[sfno] << " b up SF" << std::endl;*/

   //FINAL SFs UP DOWN TEMPLATES FOR pT1 + pt6, run G
   /*   float SF_b_up[8] = {0.88, 0.99, 0.88, 0.94, 0.89, 0.89, 0.84, 0.87};
   float SF_b_down[8] = {0.88, 1.00, 0.89, 0.95, 0.89, 0.89, 0.84, 0.87};
   float SF_c_up[8] = {0.88, 0.99, 0.89, 0.95, 0.89, 0.89, 0.84, 0.87};
   float SF_c_down[8] = {0.88, 0.99, 0.89, 0.95, 0.89, 0.89, 0.84, 0.87};
   float SF_l_up[8] = {0.88, 1.00, 0.89, 0.95, 0.90, 0.89, 0.84, 0.87};
   float SF_l_down[8] = {0.88, 1.00, 0.89, 0.95, 0.89, 0.89, 0.84, 0.87};
   float SF_cfg_up[8] = {0.88, 1.00, 0.90, 0.95, 0.90, 0.90, 0.84, 0.87};
   float SF_cfg_down[8] = {0.88, 0.99, 0.88, 0.95, 0.90, 0.89, 0.84, 0.87};
   float SF_5_temp[8] = {0.88, 0.98, 0.88, 0.94, 0.90, 0.88, 0.84, 0.86};*/

   //FINAL SFs UP DOWN TEMPLATES FOR pt5 (350-400) + pt0
   /*   float SF_b_up[8] = {1.01, 0.97, 0.96, 0.95, 0.97, 0.93, 0.85, 0.86};
   float SF_b_down[8] = {1.01, 0.98, 0.96, 0.95, 0.97, 0.94, 0.86, 0.85};
   float SF_c_up[8] = {1.01, 0.98, 0.96, 0.95, 0.97, 0.94, 0.85, 0.86};
   float SF_c_down[8] = {1.01, 0.97, 0.96, 0.95, 0.97, 0.94, 0.85, 0.86};
   float SF_l_up[8] = {1.01, 0.97, 0.96, 0.95, 0.97, 0.93, 0.85, 0.85};
   float SF_l_down[8] = {1.01, 0.97, 0.96, 0.95, 0.96, 0.94, 0.85, 0.86};
   float SF_cfg_up[8] = {1.01, 0.97, 0.96, 0.94, 0.96, 0.93, 0.86, 0.85};
   float SF_cfg_down[8] = {1.01, 0.98, 0.96, 0.96, 0.97, 0.94, 0.85, 0.86};
   float SF_5_temp[8] = {1.01, 0.89, 0.96, 0.93, 0.94, 0.93, 0.83, 0.81};
   */
   //FINAL SFs UP DOWN TEMPLATES FOR pt5 (300-350) + pt7
   /*   float SF_b_up[8] = {0.97, 0.95, 0.98, 0.96, 1.04, 0.91, 0.88, 0.91};
   float SF_b_down[8] = {0.97, 0.95, 0.98, 0.95, 1.04, 0.9, 0.88, 0.87};
   float SF_c_up[8] = {0.97, 0.95, 0.98, 0.95, 1.04, 0.9, 0.88, 0.89};
   float SF_c_down[8] = {0.97, 0.95, 0.97, 0.95, 1.04, 0.9, 0.87, 0.89};
   float SF_l_up[8] = {0.97, 0.95, 0.98, 0.95, 1.04, 0.91, SF[0], 0.89};
   float SF_l_down[8] = {0.97, 0.95, 0.98, 0.95, 1.04, 0.9, SF[0], 0.89};
   float SF_cfg_up[8] = {0.97, 0.95, 0.98, 0.94, 1.04, 0.89, SF[0], 0.88};
   float SF_cfg_down[8] = {0.98, 0.97, 0.98, 0.98, 1.05, 0.93, SF[0], 0.93};
   float SF_5_temp[8] = {0.99, 0.95, 1.05, 0.93, 1.04, 0.88, 0.87,0.86};*/

   //FINAL SFs UP DOWN TEMPLATES FOR pt4 
   float SF_b_up[8] = {0.97, 0.96, 0.97, 0.9};
   float SF_b_down[8] = {0.96, 0.95, 0.97, 0.9};
   float SF_c_up[8] = {0.96, 0.95, 0.97, 0.9};
   float SF_c_down[8] = {0.96, 0.95, 0.97, 0.9};
   float SF_l_up[8] = {0.97, 0.95, 0.97, 0.9};
   float SF_l_down[8] = {0.96, 0.95, 0.97, 0.9};
   float SF_cfg_up[8] = {0.96, 0.95, 0.97, 0.9};
   float SF_cfg_down[8] = {0.97, 0.96, 0.98, 0.9};
   float SF_5_temp[8] = {0.99, 0.95, 0.97, 0.9};
   
   //calculate sys from 5 template and up/down sys
   if( inclSYS ){
   if ( SF[0] < SF_b_up[sfno]) sigma_syst_up += (SF[0] - SF_b_up[sfno])*(SF[0] - SF_b_up[sfno]);
   else if (SF[0] > SF_b_up[sfno] ) sigma_syst_down += (SF[0] - SF_b_up[sfno])*(SF[0] - SF_b_up[sfno]);
   if ( SF[0] < SF_b_down[sfno]) sigma_syst_up += (SF[0] - SF_b_down[sfno])*(SF[0] - SF_b_down[sfno]);
   else if (SF[0] > SF_b_down[sfno]) sigma_syst_down += (SF[0] - SF_b_down[sfno])*(SF[0] - SF_b_down[sfno]);
   if ( SF[0] < SF_c_up[sfno]) sigma_syst_up += (SF[0] - SF_c_up[sfno])*(SF[0] - SF_c_up[sfno]);
   else if (SF[0] > SF_c_up[sfno]) sigma_syst_down += (SF[0] - SF_c_up[sfno])*(SF[0] - SF_c_up[sfno]);
   if ( SF[0] < SF_c_down[sfno]) sigma_syst_up += (SF[0] - SF_c_down[sfno])*(SF[0] - SF_c_down[sfno]);
   else if (SF[0] > SF_c_down[sfno]) sigma_syst_down += (SF[0] - SF_c_down[sfno])*(SF[0] - SF_c_down[sfno]);
   if ( SF[0] < SF_cfg_up[sfno]) sigma_syst_up += (SF[0] - SF_cfg_up[sfno])*(SF[0] - SF_cfg_up[sfno]);
   else if (SF[0] > SF_cfg_up[sfno]) sigma_syst_down += (SF[0] - SF_cfg_up[sfno])*(SF[0] - SF_cfg_up[sfno]);
   if ( SF[0] < SF_cfg_down[sfno]) sigma_syst_up += (SF[0] - SF_cfg_down[sfno])*(SF[0] - SF_cfg_down[sfno]);
   else if (SF[0] > SF_cfg_down[sfno]) sigma_syst_down += (SF[0] - SF_cfg_down[sfno])*(SF[0] - SF_cfg_down[sfno]);
   if ( SF[0] < SF_l_up[sfno]) sigma_syst_up += (SF[0] - SF_l_up[sfno])*(SF[0] - SF_l_up[sfno]);
   else if (SF[0] > SF_l_up[sfno]) sigma_syst_down += (SF[0] - SF_l_up[sfno])*(SF[0] - SF_l_up[sfno]);
   if ( SF[0] < SF_l_down[sfno]) sigma_syst_up += (SF[0] - SF_l_down[sfno])*(SF[0] - SF_l_down[sfno]);
   else if (SF[0] > SF_l_down[sfno] ) sigma_syst_down += (SF[0] - SF_l_down[sfno])*(SF[0] - SF_l_down[sfno]);
   if ( SF[0] < SF_5_temp[sfno]) sigma_syst_up += (SF[0] - SF_5_temp[sfno])*(SF[0] - SF_5_temp[sfno]);
   else if (SF[0] > SF_5_temp[sfno]) sigma_syst_down += (SF[0] - SF_5_temp[sfno])*(SF[0] - SF_5_temp[sfno]);
   }
   //print out SF percentage sys
   std::cout << "b + 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_b_up[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "b - 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_b_down[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "c + 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_c_up[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "c - 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_c_down[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "l + 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_l_up[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "l - 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_l_down[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "cfromg + 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_cfg_up[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "cfromg - 50% & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_cfg_down[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   std::cout << "unmerged templates & " <<  tagMod << " & " << getPtHeader(pt) << " & " << floor(((SF[0] - SF_5_temp[sfno])/SF[0])*1000+0.5)/1000 << std::endl;
   

   //taking square of syst errors and then calculating total up and down errors including syst and stat
   sigma_syst_up = sqrt(sigma_syst_up);
   sigma_syst_down = sqrt(sigma_syst_down);
   //   std::cout << "syst down " << sigma_syst_down << std::endl;
   // std::cout << "syst up " << sigma_syst_up << std::endl;
   float sigma_total_up = sqrt(sigma_stat*sigma_stat+sigma_syst_up*sigma_syst_up);
   float sigma_total_down = sqrt(sigma_stat*sigma_stat+sigma_syst_down*sigma_syst_down);
   //std::cout << "up " << sigma_total_up << " down " << sigma_total_down << std::endl;
   //printing out information for .tex file - uncomment for final run
   std::cout << runName << " & pre-tag & " << getPtHeader(pt) << " & " << floor(chi2_notag[0]*100+0.5)/100 << " & -- & -- & -- & -- & -- & -- \\" << std::endl;  
   std::cout << runName << " & " << tagMod << " & " << getPtHeader(pt) << " & " << floor(chi2_tag[0]*100+0.5)/100 << " & " << floor(effMC[0]*100+0.5)/100 << " +/- " << floor(effMCErr[0]*1000+0.5)/1000 << " & " << floor(effDATA[0]*100+0.5)/100 << " +/- " << floor(effDATAErr[0]*1000+0.5)/1000 << " & " << floor(SF[0]*100+0.5)/100 << " +/- " << floor(sigma_total_up*100+0.5)/100 << "/" << floor(sigma_total_down*100+0.5)/100 << " & " << floor(sigma_stat*1000+0.5)/1000 << " & " << floor(sigma_syst_up*1000+0.5)/1000 << " & " << floor(sigma_syst_down*1000+0.5)/1000 << std::endl;
   //printing out just SF - uncomment for +/- template run and 5 template run
   // std::cout << floor(SF[0]*100+0.5)/100 << std::endl;
   fout.close();
   
   delete cf;
}

void process(CFIT::cfit *cf,float *par_notag,float *err_notag,float *par_tag,float *err_tag,std::string sysVar,
	     std::string t1Name,std::string t2Name,std::string t3Name,std::string t4Name,std::string t5Name,
	     int statVariation)
{
   // pre-tag fit
   cf->SetSysVariation(sysVar);
   if( statVariation >= 0 ) cf->SetStatVariation(statVariation);
   cf->Run();   
   getResults(cf,par_notag,err_notag);
   chi2_notag.push_back(cf->GetChisq());
   ndof_notag.push_back(cf->GetNDOF());
   ndata_notag.push_back(cf->GetNData());
   nmc_bfromg_notag.push_back(cf->GetNTemplate(t1Name));
   nmc_bfromgErr_notag.push_back(cf->GetErrTemplate(t1Name));
   nmc_b_notag.push_back(cf->GetNTemplate(t2Name));
   nmc_bErr_notag.push_back(cf->GetErrTemplate(t2Name));
   nmc_cfromg_notag.push_back(cf->GetNTemplate(t3Name));
   nmc_c_notag.push_back(cf->GetNTemplate(t4Name));
   nmc_l_notag.push_back(cf->GetNTemplate(t5Name));
   nmc_notag.push_back(cf->GetNTemplate(t1Name)+cf->GetNTemplate(t2Name)+cf->GetNTemplate(t3Name)+cf->GetNTemplate(t4Name)+cf->GetNTemplate(t5Name));
   
   // post-tag fit
   cf->SetSysVariation(sysVar);
   if( statVariation >= 0 ) cf->SetStatVariation(statVariation);
   cf->Run("tag");   
   getResults(cf,par_tag,err_tag);
   chi2_tag.push_back(cf->GetChisq());
   ndof_tag.push_back(cf->GetNDOF());
   ndata_tag.push_back(cf->GetNData());
   nmc_bfromg_tag.push_back(cf->GetNTemplate(t1Name));
   nmc_bfromgErr_tag.push_back(cf->GetErrTemplate(t1Name));
   nmc_b_tag.push_back(cf->GetNTemplate(t2Name));
   nmc_bErr_tag.push_back(cf->GetErrTemplate(t2Name));
   nmc_cfromg_tag.push_back(cf->GetNTemplate(t3Name));
   nmc_c_tag.push_back(cf->GetNTemplate(t4Name));
   nmc_l_tag.push_back(cf->GetNTemplate(t5Name));
   nmc_tag.push_back(cf->GetNTemplate(t1Name)+cf->GetNTemplate(t2Name)+cf->GetNTemplate(t3Name)+cf->GetNTemplate(t4Name)+cf->GetNTemplate(t5Name));

   // std::cout << "notag chi2=" << chi2_notag.back() << std::endl;

   for(int i=0;i<5;i++)
     {
       //std::cout << par_notag[i] << " +- " << err_notag[i] << std::endl;
     }   
   
   // std::cout << "tag chi2=" << chi2_tag.back() << std::endl;

   for(int i=0;i<5;i++)
     {
       //	std::cout << par_tag[i] << " +- " << err_tag[i] << std::endl;
     }   

   fr_bfromg_notag.push_back(nmc_bfromg_notag.back()/nmc_notag.back());
   fr_bfromg_tag.push_back(nmc_bfromg_tag.back()/nmc_tag.back());

   fr_b_notag.push_back(nmc_b_notag.back()/nmc_notag.back());
   fr_b_tag.push_back(nmc_b_tag.back()/nmc_tag.back());

   fr_cfromg_notag.push_back(nmc_cfromg_notag.back()/nmc_notag.back());
   fr_cfromg_tag.push_back(nmc_cfromg_tag.back()/nmc_tag.back());

   fr_c_notag.push_back(nmc_c_notag.back()/nmc_notag.back());
   fr_c_tag.push_back(nmc_c_tag.back()/nmc_tag.back());
	
   fr_l_notag.push_back(nmc_l_notag.back()/nmc_notag.back());
   fr_l_tag.push_back(nmc_l_tag.back()/nmc_tag.back());
   
   float errMCErr_val = errfMC(nmc_bfromg_tag.back(),nmc_bfromgErr_tag.back(),
			       nmc_bfromg_notag.back(),nmc_bfromgErr_notag.back());
   if( doSFc ) errMCErr_val = errfMC(nmc_b_tag.back(),nmc_bErr_tag.back(),
				     nmc_b_notag.back(),nmc_bErr_notag.back());
   
   if( !doSFc ) effMC.push_back(nmc_bfromg_tag.back()/nmc_bfromg_notag.back());
   else effMC.push_back(nmc_b_tag.back()/nmc_b_notag.back());

//   effDATA.push_back(ndata_tag.back()/ndata_notag.back()*par_tag[0]/par_notag[0]*frb_tag.back()/frb_notag.back());

   if( !doSFc ) effDATA.push_back(nmc_tag.back()/nmc_notag.back()*
				  par_tag[0]/par_notag[0]*
				  fr_bfromg_tag.back()/fr_bfromg_notag.back());
   else effDATA.push_back(nmc_tag.back()/nmc_notag.back()*
				  par_tag[1]/par_notag[1]*
				  fr_b_tag.back()/fr_b_notag.back());
   
   float effMCErr_val = errfMC(nmc_bfromg_tag.back(),nmc_bfromgErr_tag.back(),
			       nmc_bfromg_notag.back(),nmc_bfromgErr_notag.back());
   if( doSFc ) effMCErr_val = errfMC(nmc_b_tag.back(),nmc_bErr_tag.back(),
				     nmc_b_notag.back(),nmc_bErr_notag.back());

   float effDATAErr_val = errfDATA_UNCOR(par_tag[0],err_tag[0],
					 par_notag[0],err_notag[0]);
   if( doSFc ) effDATAErr_val = errfDATA_UNCOR(par_tag[1],err_tag[1],
					       par_notag[1],err_notag[1]);

   effDATAErr_val *= nmc_tag.back()/nmc_notag.back();
   if( !doSFc ) effDATAErr_val *= fr_bfromg_tag.back()/fr_bfromg_notag.back();
   else effDATAErr_val *= fr_b_tag.back()/fr_b_notag.back();
   
   effMCErr.push_back(effMCErr_val);
   effDATAErr.push_back(effDATAErr_val);
   
   // std::cout << "effMC = " << effMC.back() << " +- " << effMCErr.back() << std::endl;
   //std::cout << "effDATA = " << effDATA.back() << " +- " << effDATAErr.back() << std::endl;

   float SFErr_val = errfDATA_UNCOR(effDATA.back(),effDATAErr.back(),
				    effMC.back(),effMCErr.back());

   SF.push_back(effDATA.back()/effMC.back());
   SFErr.push_back(SFErr_val);
   
   //   std::cout << "sf = " << SF.back() << " +- " << SFErr.back() << std::endl;
}

void getResults(CFIT::cfit *cf,float *par,float *err)
{   
   float chi2 = cf->GetChisq();
   int nPar = cf->GetNPar();
   for(int i=0;i<nPar;i++)
     {
	par[i] = cf->GetPar(i);
	err[i] = cf->GetParErr(i);	
     }
}

double errfMC(double v1,double ve1,double v2,double ve2)
{
   if( v2 == 0 ) return -666;
   
   double err = pow(v2-v1,2)*ve1*ve1/pow(v2,4) +
     v1*v1*(ve2*ve2-ve1*ve1)/pow(v2,4);
   
   err = sqrt(err);
   
   return err;
}

double errfDATA(double v1,double ve1,double v2,double ve2)
{
   if( v2 == 0 ) return -666;
   
   double err = ve1*ve1/v2/v2 + v1*v1*ve2*ve2/v2/v2/v2/v2;
   
   err = sqrt(err);
   
   return err;
}

double errfDATA_UNCOR(double v1,double ve1,double v2,double ve2)
{
   if( v2 == 0 ) return -666;
   
   double err = ve1*ve1/v1/v1+ve2*ve2/v2/v2;
   
   err = sqrt(err)*v1/v2;
   
   return err;
}

std::string getPtHeader(std::string pt)
{
   std::string ptHeader = "";
   
   if( pt == "pt0" ) ptHeader = "p_{T} 400-450 GeV";
   else if( pt == "pt1" ) ptHeader = "p_{T} 450-500 GeV";
   else if( pt == "pt2" ) ptHeader = "p_{T} 500-600 GeV";
   else if( pt == "pt3" ) ptHeader = "p_{T} 600-700 GeV";
   else if( pt == "pt4" ) ptHeader = "p_{T} 400-500 GeV";
   else if( pt == "pt5" ) ptHeader = "p_{T} 300-350 GeV";
   else if( pt == "pt6" ) ptHeader = "p_{T} 500-700 GeV"; 
   else if( pt == "pt7" ) ptHeader = "p_{T} 250-300 GeV"; 
   else if( pt == "ptall" ) ptHeader = "p_{T} 300-700 GeV";
   else
     {
	std::cout << "pT label is unknown" << std::endl;
	exit(1);
     }
   return ptHeader;
}

#include "TStyle.h"

void SetStyle();
TStyle* Style();

void addbin(TH1D *h);

void SFComp()
{
   gROOT->SetBatch();
   gROOT->SetStyle("Plain");
   
   gStyle->SetOptStat(0);

   SetStyle();
   
   TCanvas *c1 = new TCanvas("c1","c1",0,0,600,500);
   c1->Draw();
   c1->cd();
   
   gStyle->SetHistTopMargin(0);
   
   c1->SetLogx(0);
   c1->SetGrid(1);
   // cMVAv2L ///////////////////////////////////////////////////////////////
   
   // LT
   const int LT_xbinsN_cMVAv2L = 4;
   //middle of bin
   double LT_xbins_cMVAv2L[LT_xbinsN_cMVAv2L] = {275,350,450,600};
   //width of bin to either side of middle
   double LT_xWidth_cMVAv2L[LT_xbinsN_cMVAv2L] = {25, 50, 50, 100};
   //SF values
   double LT_SF_cMVAv2L[LT_xbinsN_cMVAv2L] = {0.91, 0.86, 0.78, 0.91};
   //sys up and down along with sys+stat box drawn around sys error
   double LT_SF_sysUp_cMVAv2L[LT_xbinsN_cMVAv2L] = {1.025,0.914,0.847,0.977};
   double LT_SF_sysDown_cMVAv2L[LT_xbinsN_cMVAv2L] = {0.795,0.806,0.713,0.843};   
   double LT_SF_sysstatx_cMVAv2L[19] = {250,275,300,350,400,450,500,600,700,700,600,500,450,400,350,300,275,250,250};
   double LT_SF_sysstaty_cMVAv2L[19] = {1.18, 1.18, 1.1, 1.01, 0.96, 0.85, 0.97, 1.08, 1.08, 0.83, 0.83, 0.75, 0.69, 0.73, 0.78, 0.73,0.68,0.68, 1.18};
   TPolyLine *sysTot_SFL = new TPolyLine(19,LT_SF_sysstatx_cMVAv2L,LT_SF_sysstaty_cMVAv2L);
   double LT_SF_sysUpRel_cMVAv2L[LT_xbinsN_cMVAv2L];
   double LT_SF_sysDownRel_cMVAv2L[LT_xbinsN_cMVAv2L];   
   for(int i=0;i<LT_xbinsN_cMVAv2L;i++)
     {
	LT_SF_sysUpRel_cMVAv2L[i] = LT_SF_sysUp_cMVAv2L[i]-LT_SF_cMVAv2L[i];
	LT_SF_sysDownRel_cMVAv2L[i] = LT_SF_cMVAv2L[i]-LT_SF_sysDown_cMVAv2L[i];
     }   

   TGraphAsymmErrors *gLT_SF_cMVAv2L = new TGraphAsymmErrors(LT_xbinsN_cMVAv2L,
							     LT_xbins_cMVAv2L,
							     LT_SF_cMVAv2L,
							     LT_xWidth_cMVAv2L,
							     LT_xWidth_cMVAv2L,
							     LT_SF_sysDownRel_cMVAv2L,
							     LT_SF_sysUpRel_cMVAv2L);

   /*   // 2TagCounting
   const int TC_xbinsN_cMVAv2L = 3;
   double TC_xbins_cMVAv2L[TC_xbinsN_cMVAv2L] = {45,90,220};
   double TC_xWidth_cMVAv2L[TC_xbinsN_cMVAv2L] = {15,30,100};
   double TC_SF_cMVAv2L[TC_xbinsN_cMVAv2L] = {0.995,0.995,0.964};
   double TC_SF_sysUpRel_cMVAv2L[TC_xbinsN_cMVAv2L] = {0.054,0.033,0.044};
   double TC_SF_sysDownRel_cMVAv2L[TC_xbinsN_cMVAv2L] = {0.054,0.033,0.044};

   TGraphAsymmErrors *gTC_SF_cMVAv2L = new TGraphAsymmErrors(TC_xbinsN_cMVAv2L,
							     TC_xbins_cMVAv2L,
							     TC_SF_cMVAv2L,
							     TC_xWidth_cMVAv2L,
							     TC_xWidth_cMVAv2L,
							     TC_SF_sysDownRel_cMVAv2L,
							     TC_SF_sysUpRel_cMVAv2L);
   */
   // cMVAv2M ///////////////////////////////////////////////////////////////
   
   const int LT_xbinsN_cMVAv2M1 = 4;
   double LT_xbins_cMVAv2M1[LT_xbinsN_cMVAv2M1] = {275, 350, 450, 600};
   double LT_xWidth_cMVAv2M1[LT_xbinsN_cMVAv2M1] = {25,50,50,100};
   double LT_SF_cMVAv2M1[LT_xbinsN_cMVAv2M1] = {0.64,0.84, 0.83, 0.85};
   double LT_SF_sysUp_cMVAv2M1[LT_xbinsN_cMVAv2M1] = {0.754,0.9, 0.922, 0.938};
   double LT_SF_sysDown_cMVAv2M1[LT_xbinsN_cMVAv2M1] = {0.526,0.78,0.738,0.762};
   double LT_SF_sysstatx_cMVAv2M1[19] = {250,275,300,350,400,450,500,600,700,700,600,500,450,400,350,300,275,250,250};
   double LT_SF_sysstaty_cMVAv2M1[19] = {1.48, 1.48, 1.2, 1.06, 1., 0.93, 0.94, 0.95, 0.95, 0.75, 0.75, 0.73, 0.7, 0.73, 0.75, 0.4,0.11,0.11, 1.48};
   TPolyLine *sysTot_SFM1 = new TPolyLine(19,LT_SF_sysstatx_cMVAv2M1,LT_SF_sysstaty_cMVAv2M1);
   double LT_SF_sysUpRel_cMVAv2M1[LT_xbinsN_cMVAv2M1];
   double LT_SF_sysDownRel_cMVAv2M1[LT_xbinsN_cMVAv2M1];   
   for(int i=0;i<LT_xbinsN_cMVAv2M1;i++)
     {
	LT_SF_sysUpRel_cMVAv2M1[i] = LT_SF_sysUp_cMVAv2M1[i]-LT_SF_cMVAv2M1[i];
	LT_SF_sysDownRel_cMVAv2M1[i] = LT_SF_cMVAv2M1[i]-LT_SF_sysDown_cMVAv2M1[i];
     }   

   TGraphAsymmErrors *gLT_SF_cMVAv2M1 = new TGraphAsymmErrors(LT_xbinsN_cMVAv2M1,
							     LT_xbins_cMVAv2M1,
							     LT_SF_cMVAv2M1, 
							     LT_xWidth_cMVAv2M1,
							     LT_xWidth_cMVAv2M1,
							     LT_SF_sysDownRel_cMVAv2M1,
							     LT_SF_sysUpRel_cMVAv2M1);

   const int LT_xbinsN_cMVAv2M2 = 4;
   double LT_xbins_cMVAv2M2[LT_xbinsN_cMVAv2M2] = {275,350,450,600};
   double LT_xWidth_cMVAv2M2[LT_xbinsN_cMVAv2M2] = {25,50,50,100};
   double LT_SF_cMVAv2M2[LT_xbinsN_cMVAv2M2] = {1.02,0.84,0.82,0.70};
   double LT_SF_sysUp_cMVAv2M2[LT_xbinsN_cMVAv2M2] = {1.307,0.931,0.922,0.877};
   double LT_SF_sysDown_cMVAv2M2[LT_xbinsN_cMVAv2M2] = {0.733,0.749,0.718,0.523};
   double LT_SF_sysstatx_cMVAv2M2[19] = {250,275,300,350,400,450,500,600,700,700,600,500,450,400,350,300,275,250,250};
   double LT_SF_sysstaty_cMVAv2M2[19] = {1.4, 1.4, 1.44, 1.47, 1.2, 0.93, 0.91, 0.89, 0.89, 0.52, 0.52, 0.6, 0.7, 0.6, 0.57, 0.54,0.51,0.51, 1.4};
   TPolyLine *sysTot_SFM2 = new TPolyLine(19,LT_SF_sysstatx_cMVAv2M2,LT_SF_sysstaty_cMVAv2M2);   
   double LT_SF_sysUpRel_cMVAv2M2[LT_xbinsN_cMVAv2M2];
   double LT_SF_sysDownRel_cMVAv2M2[LT_xbinsN_cMVAv2M2];   
   for(int i=0;i<LT_xbinsN_cMVAv2M2;i++)
     {
	LT_SF_sysUpRel_cMVAv2M2[i] = LT_SF_sysUp_cMVAv2M2[i]-LT_SF_cMVAv2M2[i];
	LT_SF_sysDownRel_cMVAv2M2[i] = LT_SF_cMVAv2M2[i]-LT_SF_sysDown_cMVAv2M2[i];
     }   

   TGraphAsymmErrors *gLT_SF_cMVAv2M2 = new TGraphAsymmErrors(LT_xbinsN_cMVAv2M2,
							     LT_xbins_cMVAv2M2,
							     LT_SF_cMVAv2M2, 
							     LT_xWidth_cMVAv2M2,
							     LT_xWidth_cMVAv2M2,
							     LT_SF_sysDownRel_cMVAv2M2,
							     LT_SF_sysUpRel_cMVAv2M2);

   /*   // 2TagCounting
   const int TC_xbinsN_cMVAv2M = 3;
   double TC_xbins_cMVAv2M[TC_xbinsN_cMVAv2M] = {45,90,220};
   double TC_xWidth_cMVAv2M[TC_xbinsN_cMVAv2M] = {15,30,100};
   double TC_SF_cMVAv2M[TC_xbinsN_cMVAv2M] = {0.966,0.971,0.927};
   double TC_SF_sysUpRel_cMVAv2M[TC_xbinsN_cMVAv2M] = {0.032,0.017,0.028};
   double TC_SF_sysDownRel_cMVAv2M[TC_xbinsN_cMVAv2M] = {0.032,0.017,0.028};

   TGraphAsymmErrors *gTC_SF_cMVAv2M = new TGraphAsymmErrors(TC_xbinsN_cMVAv2M,
							     TC_xbins_cMVAv2M,
							     TC_SF_cMVAv2M,
							     TC_xWidth_cMVAv2M,
							     TC_xWidth_cMVAv2M,
							     TC_SF_sysDownRel_cMVAv2M,
							     TC_SF_sysUpRel_cMVAv2M);*/

   // cMVAv2T ///////////////////////////////////////////////////////////////
   
   // LT
   const int LT_xbinsN_cMVAv2T = 4;
   double LT_xbins_cMVAv2T[LT_xbinsN_cMVAv2T] = {275,350,450,600};
   double LT_xWidth_cMVAv2T[LT_xbinsN_cMVAv2T] = {25,50,50,100};
   double LT_SF_cMVAv2T[LT_xbinsN_cMVAv2T] = {0.94,0.88,1.0,0.72};
   double LT_SF_sysUp_cMVAv2T[LT_xbinsN_cMVAv2T] = {1.365,0.998,1.125,0.98};
   double LT_SF_sysDown_cMVAv2T[LT_xbinsN_cMVAv2T] = {0.515,0.762,0.875,0.46};
   double LT_SF_sysstatx_cMVAv2T[19] = {250,275,300,350,400,450,500,600,700,700,600,500,450,400,350,300,275,250,250};
   double LT_SF_sysstaty_cMVAv2T[19] = {1.8, 1.8, 1.4, 1.12, 1.15, 1.18, 1.09, 1., 1., 0.45,0.45, 0.65, 0.86, 0.77, 0.7, 0.4,0.17,0.17, 1.8};
   TPolyLine *sysTot_SFT = new TPolyLine(19,LT_SF_sysstatx_cMVAv2T,LT_SF_sysstaty_cMVAv2T);   
   double LT_SF_sysUpRel_cMVAv2T[LT_xbinsN_cMVAv2T];
   double LT_SF_sysDownRel_cMVAv2T[LT_xbinsN_cMVAv2T];   
   for(int i=0;i<LT_xbinsN_cMVAv2T;i++)
     {
	LT_SF_sysUpRel_cMVAv2T[i] = LT_SF_sysUp_cMVAv2T[i]-LT_SF_cMVAv2T[i];
	LT_SF_sysDownRel_cMVAv2T[i] = LT_SF_cMVAv2T[i]-LT_SF_sysDown_cMVAv2T[i];
     }   

   TGraphAsymmErrors *gLT_SF_cMVAv2T = new TGraphAsymmErrors(LT_xbinsN_cMVAv2T,
							     LT_xbins_cMVAv2T,
							     LT_SF_cMVAv2T,
							     LT_xWidth_cMVAv2T,
							     LT_xWidth_cMVAv2T,
							     LT_SF_sysDownRel_cMVAv2T,
							     LT_SF_sysUpRel_cMVAv2T);

   // 2TagCounting
   /*   const int TC_xbinsN_cMVAv2T = 3;
   double TC_xbins_cMVAv2T[TC_xbinsN_cMVAv2T] = {45,90,220};
   double TC_xWidth_cMVAv2T[TC_xbinsN_cMVAv2T] = {15,30,100};
   double TC_SF_cMVAv2T[TC_xbinsN_cMVAv2T] = {0.940,0.959,0.917};
   double TC_SF_sysUpRel_cMVAv2T[TC_xbinsN_cMVAv2T] = {0.037,0.020,0.032};
   double TC_SF_sysDownRel_cMVAv2T[TC_xbinsN_cMVAv2T] = {0.037,0.020,0.032};

   TGraphAsymmErrors *gTC_SF_cMVAv2T = new TGraphAsymmErrors(TC_xbinsN_cMVAv2T,
							     TC_xbins_cMVAv2T,
							     TC_SF_cMVAv2T,
							     TC_xWidth_cMVAv2T,
							     TC_xWidth_cMVAv2T,
							     TC_SF_sysDownRel_cMVAv2T,
							     TC_SF_sysUpRel_cMVAv2T);*/
   
   // cMVAv2L
     {   	
       //plot for each WP
	gLT_SF_cMVAv2L->SetMarkerSize(1.2);
	gLT_SF_cMVAv2L->SetMarkerStyle(24);
	gLT_SF_cMVAv2L->SetLineWidth(3);
	gLT_SF_cMVAv2L->SetMarkerColor(46);
	gLT_SF_cMVAv2L->SetLineColor(46);
	gLT_SF_cMVAv2L->Draw("AP"); 
	 
	gLT_SF_cMVAv2L->GetXaxis()->SetMoreLogLabels();
	gLT_SF_cMVAv2L->GetXaxis()->SetNoExponent();
	gLT_SF_cMVAv2L->SetMinimum(0.);
	gLT_SF_cMVAv2L->SetMaximum(2.0);
	gLT_SF_cMVAv2L->GetXaxis()->SetTitle("p_{T} [GeV]");
	gLT_SF_cMVAv2L->GetYaxis()->SetTitle("SF_{double b}");
	sysTot_SFL->SetFillColor(7);
	sysTot_SFL->Draw("PS");
	//sysTot_SFL->Draw();
	  
	
	/*	gTC_SF_cMVAv2L->SetMarkerColor(38);
	gTC_SF_cMVAv2L->SetLineColor(38);
	gTC_SF_cMVAv2L->SetMarkerStyle(22);
	gTC_SF_cMVAv2L->SetMarkerSize(1.2);
	gTC_SF_cMVAv2L->SetLineWidth(3);
	gTC_SF_cMVAv2L->Draw("PS");
	*/
	TLegend *leg_cMVAv2L;
	leg_cMVAv2L = new TLegend(0.35,0.35,0.70,0.20);
	leg_cMVAv2L->SetFillColor(253);
	leg_cMVAv2L->SetBorderSize(0);
	
	TLatex *legl_cMVAv2L = new TLatex();
	legl_cMVAv2L->SetNDC();
	legl_cMVAv2L->SetTextAlign(22);
	legl_cMVAv2L->SetTextFont(63);
	legl_cMVAv2L->SetTextSizePixels(30);
	legl_cMVAv2L->DrawLatex(0.55,0.82,"DoubleBL");
	
	leg_cMVAv2L->AddEntry(gLT_SF_cMVAv2L,"LT","lp");
	//	leg_cMVAv2L->AddEntry(gTC_SF_cMVAv2L,"2TagCounting","lp");
	
	//	leg_cMVAv2L->Draw();
	
	c1->RedrawAxis("g");
	
	c1->Print("pics/SFComp_cMVAv2L.eps");
	c1->Clear();
     }

   // cMVAv2M
     {   	
	gLT_SF_cMVAv2M1->SetMarkerSize(1.2);
	gLT_SF_cMVAv2M1->SetMarkerStyle(24);
	gLT_SF_cMVAv2M1->SetLineWidth(3);
	gLT_SF_cMVAv2M1->SetMarkerColor(46);
	gLT_SF_cMVAv2M1->SetLineColor(46);
	gLT_SF_cMVAv2M1->Draw("AP");   
	
	gLT_SF_cMVAv2M1->GetXaxis()->SetMoreLogLabels();
	gLT_SF_cMVAv2M1->GetXaxis()->SetNoExponent();
	
	gLT_SF_cMVAv2M1->SetMinimum(0.);
	gLT_SF_cMVAv2M1->SetMaximum(2.0);
	gLT_SF_cMVAv2M1->GetXaxis()->SetTitle("p_{T} [GeV]");
	gLT_SF_cMVAv2M1->GetYaxis()->SetTitle("SF_{double b}");
	sysTot_SFM1->SetFillColor(7);
	sysTot_SFM1->Draw("PS");
	
	/*	gTC_SF_cMVAv2M1->SetMarkerColor(38);
	gTC_SF_cMVAv2M1->SetLineColor(38);
	gTC_SF_cMVAv2M1->SetMarkerStyle(22);
	gTC_SF_cMVAv2M1->SetMarkerSize(1.2);
	gTC_SF_cMVAv2M1->SetLineWidth(3);
	gTC_SF_cMVAv2M1->Draw("PS");
	*/
	TLegend *leg_cMVAv2M1;
	leg_cMVAv2M1 = new TLegend(0.35,0.35,0.70,0.20);
	leg_cMVAv2M1->SetFillColor(253);
	leg_cMVAv2M1->SetBorderSize(0);
	
	TLatex *legl_cMVAv2M1 = new TLatex();
	legl_cMVAv2M1->SetNDC();
	legl_cMVAv2M1->SetTextAlign(22);
	legl_cMVAv2M1->SetTextFont(63);
	legl_cMVAv2M1->SetTextSizePixels(30);
	legl_cMVAv2M1->DrawLatex(0.55,0.82,"DoubleBM1");
	
	leg_cMVAv2M1->AddEntry(gLT_SF_cMVAv2M1,"LT","lp");
	//	leg_cMVAv2M->AddEntry(gTC_SF_cMVAv2M,"2TagCounting","lp");
	
	//	leg_cMVAv2M1->Draw();
	
	c1->RedrawAxis("g");
	
	c1->Print("pics/SFComp_cMVAv2M1.eps");
	c1->Clear();
     }

   // cMVAv2M
     {   	
	gLT_SF_cMVAv2M2->SetMarkerSize(1.2);
	gLT_SF_cMVAv2M2->SetMarkerStyle(24);
	gLT_SF_cMVAv2M2->SetLineWidth(3);
	gLT_SF_cMVAv2M2->SetMarkerColor(46);
	gLT_SF_cMVAv2M2->SetLineColor(46);
	gLT_SF_cMVAv2M2->Draw("AP");   
	
	gLT_SF_cMVAv2M2->GetXaxis()->SetMoreLogLabels();
	gLT_SF_cMVAv2M2->GetXaxis()->SetNoExponent();
	
	gLT_SF_cMVAv2M2->SetMinimum(0.);
	gLT_SF_cMVAv2M2->SetMaximum(2.0);
	gLT_SF_cMVAv2M2->GetXaxis()->SetTitle("p_{T} [GeV]");
	gLT_SF_cMVAv2M2->GetYaxis()->SetTitle("SF_{double b}");
        sysTot_SFM2->SetFillColor(7);
	sysTot_SFM2->Draw("PS");
	/*	gTC_SF_cMVAv2M2->SetMarkerColor(38);
	gTC_SF_cMVAv2M2->SetLineColor(38);
	gTC_SF_cMVAv2M2->SetMarkerStyle(22);
	gTC_SF_cMVAv2M2->SetMarkerSize(1.2);
	gTC_SF_cMVAv2M2->SetLineWidth(3);
	gTC_SF_cMVAv2M2->Draw("PS");
	*/
	TLegend *leg_cMVAv2M2;
	leg_cMVAv2M2 = new TLegend(0.35,0.35,0.70,0.20);
	leg_cMVAv2M2->SetFillColor(253);
	leg_cMVAv2M2->SetBorderSize(0);
	
	TLatex *legl_cMVAv2M2 = new TLatex();
	legl_cMVAv2M2->SetNDC();
	legl_cMVAv2M2->SetTextAlign(22);
	legl_cMVAv2M2->SetTextFont(63);
	legl_cMVAv2M2->SetTextSizePixels(30);
	legl_cMVAv2M2->DrawLatex(0.55,0.82,"DoubleBM2");
	
	leg_cMVAv2M2->AddEntry(gLT_SF_cMVAv2M2,"LT","lp");
	//	leg_cMVAv2M->AddEntry(gTC_SF_cMVAv2M,"2TagCounting","lp");
	
	//	leg_cMVAv2M2->Draw();
	
	c1->RedrawAxis("g");
	
	c1->Print("pics/SFComp_cMVAv2M2.eps");
	c1->Clear();
     }

   // cMVAv2T
     {   	
	gLT_SF_cMVAv2T->SetMarkerSize(1.2);
	gLT_SF_cMVAv2T->SetMarkerStyle(24);
	gLT_SF_cMVAv2T->SetLineWidth(3);
	gLT_SF_cMVAv2T->SetMarkerColor(46);
	gLT_SF_cMVAv2T->SetLineColor(46);
	gLT_SF_cMVAv2T->Draw("AP");   
	
	gLT_SF_cMVAv2T->GetXaxis()->SetMoreLogLabels();
	gLT_SF_cMVAv2T->GetXaxis()->SetNoExponent();
	
	gLT_SF_cMVAv2T->SetMinimum(0.);
	gLT_SF_cMVAv2T->SetMaximum(2.0);
	gLT_SF_cMVAv2T->GetXaxis()->SetTitle("p_{T} [GeV]");
	gLT_SF_cMVAv2T->GetYaxis()->SetTitle("SF_{double b}");
        sysTot_SFT->SetFillColor(7);
        sysTot_SFT->Draw("PS");
	
	/*	gTC_SF_cMVAv2T->SetMarkerColor(38);
	gTC_SF_cMVAv2T->SetLineColor(38);
	gTC_SF_cMVAv2T->SetMarkerStyle(22);
	gTC_SF_cMVAv2T->SetMarkerSize(1.2);
	gTC_SF_cMVAv2T->SetLineWidth(3);
	gTC_SF_cMVAv2T->Draw("PS");
	*/
	TLegend *leg_cMVAv2T;
	leg_cMVAv2T = new TLegend(0.35,0.35,0.70,0.20);
	leg_cMVAv2T->SetFillColor(253);
	leg_cMVAv2T->SetBorderSize(0);
	
	TLatex *legl_cMVAv2T = new TLatex();
	legl_cMVAv2T->SetNDC();
	legl_cMVAv2T->SetTextAlign(22);
	legl_cMVAv2T->SetTextFont(63);
	legl_cMVAv2T->SetTextSizePixels(30);
	legl_cMVAv2T->DrawLatex(0.55,0.82,"DoubleBT");
	
	leg_cMVAv2T->AddEntry(gLT_SF_cMVAv2T,"LT","lp");
	//leg_cMVAv2T->AddEntry(gTC_SF_cMVAv2T,"2TagCounting","lp");
	
	//	leg_cMVAv2T->Draw();
	
	c1->RedrawAxis("g");
	
	c1->Print("pics/SFComp_cMVAv2T.eps");
	c1->Clear();
     }
   
   gApplication->Terminate();
}

void SetStyle()
{
  static TStyle* style = 0;
//  std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
  if( style==0 ) style = Style();
  gROOT->SetStyle("STYLE");
  gROOT->ForceStyle();
}

TStyle* Style() 
{
  TStyle *style = new TStyle("STYLE","User style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  style->SetFrameBorderMode(icol);
  style->SetFrameFillColor(icol);
  style->SetCanvasBorderMode(icol);
  style->SetCanvasColor(icol);
  style->SetPadBorderMode(icol);
  style->SetPadColor(icol);
  style->SetStatColor(icol);
  //style->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  style->SetPaperSize(20,26);

  // set margin sizes
  style->SetPadTopMargin(0.05);
  style->SetPadRightMargin(0.05);
  style->SetPadBottomMargin(0.16);
  style->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  style->SetTitleXOffset(1.4);
  style->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  style->SetTextFont(font);

  style->SetTextSize(tsize);
  style->SetLabelFont(font,"x");
  style->SetTitleFont(font,"x");
  style->SetLabelFont(font,"y");
  style->SetTitleFont(font,"y");
  style->SetLabelFont(font,"z");
  style->SetTitleFont(font,"z");
  
  style->SetLabelSize(tsize,"x");
  style->SetTitleSize(tsize,"x");
  style->SetLabelSize(tsize,"y");
  style->SetTitleSize(tsize,"y");
  style->SetLabelSize(tsize,"z");
  style->SetTitleSize(tsize,"z");

  // use bold lines and markers
  style->SetMarkerStyle(20);
  style->SetMarkerSize(1.2);
  style->SetHistLineWidth(2.);
  style->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //style->SetErrorX(0.001);
  // get rid of error bar caps
  style->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  style->SetOptTitle(0);
  //style->SetOptStat(1111);
  style->SetOptStat(0);
  //style->SetOptFit(1111);
  style->SetOptFit(0);

  // put tick marks on top and RHS of plots
  style->SetPadTickX(1);
  style->SetPadTickY(1);

  return style;
}

void addbin(TH1D *h)
{   
   // Add overflow and underflow bins
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
   h->SetBinError(1,TMath::Sqrt(pow(h->GetBinError(0),2)+pow(h->GetBinError(1),2)));
   h->SetBinContent(x_nbins,h->GetBinContent(x_nbins)+h->GetBinContent(x_nbins+1));
   h->SetBinError(x_nbins,TMath::Sqrt(pow(h->GetBinError(x_nbins),2)+
				      pow(h->GetBinError(x_nbins+1),2)));
   // Set overflow and underflow bins to 0
   h->SetBinContent(0,0.);
   h->SetBinError(0,0.);
   h->SetBinContent(x_nbins+1,0.);
   h->SetBinError(x_nbins+1,0.);
}

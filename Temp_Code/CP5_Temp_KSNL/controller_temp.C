void controller_temp()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(2);	
  gStyle->SetLabelOffset(0.002,"Y");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetLabelSize(0.04, "X");
  gStyle->SetNdivisions(510,"X");
  gStyle->SetTitleSize(0.05, "X" );
  gStyle->SetTitleSize(0.06, "Y" );
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetLabelFont(22,"X");
  gStyle->SetLabelFont(22,"Y");
  gStyle->SetTitleOffset(0.9,"X");
  gStyle->SetTitleOffset(0.84,"Y");
  
  TCanvas *plot = new TCanvas("plot","",1200,1200);
  TPad *pad1 = new TPad("pad1", "",0.0,0.0,1.0,1.0);
  pad1->Draw();
  pad1->cd();
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.03);
  pad1->SetBottomMargin(0.1);
  pad1->SetLeftMargin(0.10);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameBorderMode(0);
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetGridx(1);
  pad1->SetGridy(1);

  TH2F *frame = new TH2F("frame","",10, 0, 500, 10, 28, 40); 
  //TH2F *frame = new TH2F("frame","",10, 423, 470, 10, 28, 40); 
  frame->GetXaxis()->SetTitle("Days");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("Controller temp (^{o}C)");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleOffset(0.84);
  frame->Draw("");

  //ECGe PPC 1500
  
  const int data = 66604;
  char err_1[data][20];
  
  int days_1[data], month_1[data],year_1[data], hours_1[data], min_1[data], sec_1[data];
  float coolar_watt_1[data], dip_temp_1[data], dip_set_1[data], coldhead_temp_1[data],  compressor_temp_1[data],  controller_temp_1[data]; 

  float Xaxis_1[data];
  
  FILE *file1 =fopen("TEST_70X70.txt","r");
  
  for(int ll =0; ll<data; ll++)
    {
      fscanf(file1, "%d %d %d %d %d %d", &days_1[ll],&month_1[ll],&year_1[ll],&hours_1[ll], &min_1[ll], &sec_1[ll]);
      fscanf(file1, "%f %f %f  %f %f %f \n", &coolar_watt_1[ll],  &dip_temp_1[ll],  &dip_set_1[ll],  &coldhead_temp_1[ll], &compressor_temp_1[ll], &controller_temp_1[ll]);
      
      //printf( "%d %d %d %d %d %d", days_1[ll],month_1[ll],year_1[ll],hours_1[ll],min_1[ll],sec_1[ll]);
      //printf("\t%f %f %f  %f %f %f \n", coolar_watt_1[ll],  dip_temp_1[ll],  dip_set_1[ll],  coldhead_temp_1[ll], compressor_temp_1[ll], controller_temp_1[ll]);

      Xaxis_1[ll] = ((double)ll*10.0/(60.0*24.0)); 
    }
  printf("file 1st Over\n");

  //ECGe PPC 500

  const int range = 24210;
  char err_2[range][20];
  
  int days_2[range], month_2[range],year_2[range], hours_2[range], min_2[range], sec_2[range];
  float coolar_watt_2[range], dip_temp_2[range], dip_set_2[range], coldhead_temp_2[range],  compressor_temp_2[range],  controller_temp_2[range]; 

  float Xaxis_2[range];
  
  FILE *file2 =fopen("TEST_50X50.txt","r");
  
  for(int ll =0; ll<range; ll++)
    {
      fscanf(file2, "%d %d %d %d %d %d", &days_2[ll],&month_2[ll],&year_2[ll],&hours_2[ll], &min_2[ll], &sec_2[ll]);
      fscanf(file2, "%f %f %f  %f %f %f \n", &coolar_watt_2[ll],  &dip_temp_2[ll],  &dip_set_2[ll],  &coldhead_temp_2[ll], &compressor_temp_2[ll], &controller_temp_2[ll]);
      
      //printf( "%d %d %d %d %d %d", days_2[ll],month_2[ll],year_2[ll],hours_2[ll],min_2[ll],sec_2[ll]);
      //printf("\t%f %f %f  %f %f %f \n", coolar_watt_2[ll],  dip_temp_2[ll],  dip_set_2[ll],  coldhead_temp_2[ll], compressor_temp_2[ll], controller_temp_2[ll]);

       if(ll>=1614 && ll<=2359)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0; 
	}
      else if(ll>=2360 && ll<=2454)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0; 
	}
      else if(ll>=2455 && ll<=2799)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0; 
	}
      else if(ll>=2800 && ll<=2946)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0; 
	}
      else if(ll>=2947 && ll<=4443)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0; 
	}
      else if(ll>=4444 && ll<=5111)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0; 
	}
      else if(ll>=5112 && ll<=8243)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+5.0; 
	}
      else if(ll>=8244 && ll<=8728)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+4.6; 
	}
      else if(ll>=8729 && ll<=13139)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+7.0+4.0+2.6+5.0;  //<--Peak--Fluctuation
	}
      else if(ll>=13140 && ll<=18276)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0;
	}
      else if(ll>=18277 && ll<=18644)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0;
	}
      else if(ll>=18645 && ll<=18988)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0;
	}
      else if(ll>=18989 && ll<=19433)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0;
	}
      else if(ll>=19434 && ll<=19606)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0;
	}
      else if(ll>=19607 && ll<=19848)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0;
	}
      else if(ll>=19849 && ll<=19952)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0;
	}
      else if(ll>=19953 && ll<=20126)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0;
	}
      else if(ll>=20127 && ll<=20204)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0;
	}
      else if(ll>=20205 && ll<=20967)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0;
	}
      else if(ll>=20968 && ll<=21218)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0;
	}
      else if(ll>=21219 && ll<=21597)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0;
	}
      else if(ll>=21598 && ll<=23183)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0;
	}
      else if(ll>=23184 && ll<=23623)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0;
	}
      else if(ll>=23624 && ll<=23633)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0+8.0;
	}
      else if(ll>=23634 && ll<=23967)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0+8.0+5.0;
	}
      else if(ll>=23968 && ll<=24210)
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0+8.0+5.0+2.0;
	}
      else
	{
	  Xaxis_2[ll] = 71.0 + ((double)ll*10.0/(60.0*24.0)); 
	}
       
      printf("%f  %d\n", Xaxis_2[ll], ll);
      
      // Xaxis_2[ll] = ((double)ll*10.0/(60.0*24.0)); 
    }
  printf("file 2nd Over\n");

  TGraph *gr_1 = new TGraph(data, Xaxis_1, controller_temp_1);
  gr_1->SetLineColor(1);
  gr_1->SetLineWidth(3);
  gr_1->Draw("l");  
  
  TGraph *gr_2 = new TGraph(range, Xaxis_2, controller_temp_2);
  gr_2->SetLineColor(2);
  gr_2->SetLineWidth(3);
  gr_2->Draw("l");  
  
  TLatex *tex = new TLatex(0.48,0.91,"ECGe 70#times70");
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(22);
  tex->SetTextSize(0.055);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex *tex0 = new TLatex(0.48,0.85,"ECGe 50#times50");
  tex0->SetNDC();
  tex0->SetTextColor(2);
  tex0->SetTextFont(22);
  tex0->SetTextSize(0.055);
  tex0->SetLineWidth(2);
  tex0->Draw();

  /*
  TLatex *tex1 = new TLatex(0.11,0.71,"HV ON");
  tex1->SetNDC();
  tex1->SetTextColor(1);
  tex1->SetTextFont(22);
  tex1->SetTextSize(0.055);
  tex1->SetLineWidth(2);
  tex1->Draw();
  
  TLatex *tex2 = new TLatex(0.71,0.61,"-195^{o}C");
  tex2->SetNDC();
  tex2->SetTextColor(2);
  tex2->SetTextFont(22);
  tex2->SetTextSize(0.036);
  tex2->SetLineWidth(2);
  tex2->Draw();
  
  TLatex *tex_3 = new TLatex(0.11,0.61,"Temp. changed (-190^{o}C -> -200^{o}C)");
  tex_3->SetNDC();
  tex_3->SetTextColor(1);
  tex_3->SetTextFont(22);
  tex_3->SetTextSize(0.036);
  tex_3->SetLineWidth(2);
  tex_3->Draw();
  
  TLatex *tex_4 = new TLatex(0.61,0.41,"Temp. changed (-200^{o}C -> -195^{o}C)");
  tex_4->SetNDC();
  tex_4->SetTextColor(1);
  tex_4->SetTextFont(22);
  tex_4->SetTextSize(0.036);
  tex_4->SetLineWidth(2);
  tex_4->Draw();
  
  TLatex *tex11 = new TLatex(0.61,0.31,"Power OFF");
  tex11->SetNDC();
  tex11->SetTextColor(1);
  tex11->SetTextFont(22);
  tex11->SetTextSize(0.055);
  tex11->SetLineWidth(2);
  tex11->Draw();
  
  TLatex *tex112 = new TLatex(0.31,0.21,"-200.0^{o}C");
  tex112->SetNDC();
  tex112->SetTextColor(1);
  tex112->SetTextFont(22);
  tex112->SetTextSize(0.055);
  tex112->SetLineWidth(2);
  tex112->Draw();

  TLatex *tex1122 = new TLatex(0.31,0.41,"-190.0^{o}C");
  tex1122->SetNDC();
  tex1122->SetTextColor(2);
  tex1122->SetTextFont(22);
  tex1122->SetTextSize(0.055);
  tex1122->SetLineWidth(2);
  tex1122->Draw();
  */

}

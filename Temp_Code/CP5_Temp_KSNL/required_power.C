void required_power()
{
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetFrameLineWidth(2);	
   gStyle->SetLabelOffset(0.007,"Y");
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
   gStyle->SetTitleOffset(0.88,"Y");
  
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
   pad1->SetLeftMargin(0.11);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderMode(0);
   pad1->SetTickx(1);
   pad1->SetTicky(1);
   pad1->SetGridx(1);
   pad1->SetGridy(1);

   TH2F *frame = new TH2F("frame","",10, 0.0, 500, 10, 80, 180);
   //TH2F *frame = new TH2F("frame","",10, 423.0, 470, 10, 80, 180);
   frame->GetXaxis()->SetTitle("Days");
   frame->GetXaxis()->CenterTitle();
   frame->GetYaxis()->SetTitle("Cooler Power (watts)");
   frame->GetYaxis()->CenterTitle();
   frame->GetYaxis()->SetTitleSize(0.06);
   frame->GetYaxis()->SetTitleOffset(0.88);
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
      fscanf(file1, "%d  %d  %d  %d  %d  %d  %f  %f  %f  %f  %f  %f \n",   &days_1[ll],  &month_1[ll],  &year_1[ll],  &hours_1[ll], &min_1[ll], &sec_1[ll], &coolar_watt_1[ll],  &dip_temp_1[ll],  &dip_set_1[ll],  &coldhead_temp_1[ll], &compressor_temp_1[ll], &controller_temp_1[ll]);

      
      //  printf("%d %d %d %d %d %d ", days_1[ll],month_1[ll],year_1[ll],hours_1[ll],min_1[ll],sec_1[ll]);
      //  printf("\t %f %f %f  %f %f %f \n", coolar_watt_1[ll],  dip_temp_1[ll],  dip_set_1[ll],  coldhead_temp_1[ll], compressor_temp_1[ll], controller_temp_1[ll]);   

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
  
  for(int jj =0; jj<range; jj++)
    {
      fscanf(file2, "%d  %d  %d  %d  %d  %d  %f  %f  %f  %f  %f  %f \n", &days_2[jj],&month_2[jj],&year_2[jj],&hours_2[jj], &min_2[jj], &sec_2[jj], &coolar_watt_2[jj],  &dip_temp_2[jj],  &dip_set_2[jj],  &coldhead_temp_2[jj], &compressor_temp_2[jj], &controller_temp_2[jj]);
      
      //  printf( "%d %d %d %d %d %d", days_2[jj],month_2[jj],year_2[jj],hours_2[jj],min_2[jj],sec_2[jj]);
      // printf("\t %f %f %f  %f %f %f \n", coolar_watt_2[jj],  dip_temp_2[jj],  dip_set_2[jj],  coldhead_temp_2[jj], compressor_temp_2[jj], controller_temp_2[jj]);
      

      if(jj>=1614 && jj<=2359)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0; 
	}
      else if(jj>=2360 && jj<=2454)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0; 
	}
      else if(jj>=2455 && jj<=2799)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0; 
	}
      else if(jj>=2800 && jj<=2946)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0; 
	}
      else if(jj>=2947 && jj<=4443)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0; 
	}
      else if(jj>=4444 && jj<=5111)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0; 
	}
      else if(jj>=5112 && jj<=8243)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+5.0; 
	}
      else if(jj>=8244 && jj<=8728)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+4.6; 
	}
      else if(jj>=8729 && jj<=13139)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+7.0+4.0+2.6+5.0;  //<--Peak--Fluctuation
	}
      else if(jj>=13140 && jj<=18276)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0;
	}
      else if(jj>=18277 && jj<=18644)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0;
	}
      else if(jj>=18645 && jj<=18988)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0;
	}
      else if(jj>=18989 && jj<=19433)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0;
	}
      else if(jj>=19434 && jj<=19606)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0;
	}
      else if(jj>=19607 && jj<=19848)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0;
	}
      else if(jj>=19849 && jj<=19952)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0;
	}
      else if(jj>=19953 && jj<=20126)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0;
	}
      else if(jj>=20127 && jj<=20204)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0;
	}
      else if(jj>=20205 && jj<=20967)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0;
	}
      else if(jj>=20968 && jj<=21218)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0;
	}
      else if(jj>=21219 && jj<=21597)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0;
	}
      else if(jj>=21598 && jj<=23183)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0;
	}
      else if(jj>=23184 && jj<=23623)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0;
	}
      else if(jj>=23624 && jj<=23633)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0+8.0;
	}
      else if(jj>=23634 && jj<=23967)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0+8.0+5.0;
	}
      else if(jj>=23968 && jj<=24210)
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)) +1.0+1.0+4.0+4.0+9.0+5.0+3.0+7.0+11.0+24.0+5.0+38.0+11.0+2.0+4.0+20.0+18.0+20.0+9.0+3.0+3.0+6.0+8.0+5.0+2.0;
	}

      else
	{
	  Xaxis_2[jj] = 71.0 + ((double)jj*10.0/(60.0*24.0)); 
	}
      printf("%f  %d\n", Xaxis_2[jj], jj);
    }
  printf("file 2nd Over\n");
  
  
  TGraph *gr_1 = new TGraph(data, Xaxis_1, coolar_watt_1);
  gr_1->SetLineColor(1);
  gr_1->SetLineWidth(3);
  gr_1->Draw("l");  
  
  TGraph *gr_2 = new TGraph(range, Xaxis_2, coolar_watt_2);
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
  
  TLatex *tex2 = new TLatex(0.11,0.61,"Temp. changed (-190^{o}C -> -200^{o}C)");
  tex2->SetNDC();
  tex2->SetTextColor(1);
  tex2->SetTextFont(22);
  tex2->SetTextSize(0.036);
  tex2->SetLineWidth(2);
  tex2->Draw();
  
  TLatex *tex_4 = new TLatex(0.61,0.41,"Temp. changed (-200^{o}C -> -195^{o}C)");
  tex_4->SetNDC();
  tex_4->SetTextColor(1);
  tex_4->SetTextFont(22);
  tex_4->SetTextSize(0.036);
  tex_4->SetLineWidth(2);
  tex_4->Draw();
  
  TLatex *tex_2 = new TLatex(0.31,0.51,"T_{Coldtip} = -190^{o}C");
  tex_2->SetNDC();
  tex_2->SetTextColor(2);
  tex_2->SetTextFont(22);
  tex_2->SetTextSize(0.036);
  tex_2->SetLineWidth(2);
  tex_2->Draw();

  TLatex *tex_3 = new TLatex(0.41,0.61,"T_{Coldtip} = -195^{o}C");
  tex_3->SetNDC();
  tex_3->SetTextColor(2);
  tex_3->SetTextFont(22);
  tex_3->SetTextSize(0.036);
  tex_3->SetLineWidth(2);
  tex_3->Draw();
  
  TLatex *tex11 = new TLatex(0.61,0.31,"T_{Coldtip} = -200^{o}C");
  tex11->SetNDC();
  tex11->SetTextColor(1);
  tex11->SetTextFont(22);
  tex11->SetTextSize(0.055);
  tex11->SetLineWidth(2);
  tex11->Draw();
  
  TLatex *tex112 = new TLatex(0.31,0.21,"Stable power (105.6 - 110.0) watts");
  tex112->SetNDC();
  tex112->SetTextColor(1);
  tex112->SetTextFont(22);
  tex112->SetTextSize(0.055);
  tex112->SetLineWidth(2);
  //tex112->Draw();

  TLatex *tex1122 = new TLatex(0.31,0.51,"Stable power (127.0 - 132.32) watts");
  tex1122->SetNDC();
  tex1122->SetTextColor(2);
  tex1122->SetTextFont(22);
  tex1122->SetTextSize(0.055);
  tex1122->SetLineWidth(2);
  //tex1122->Draw();
  */
}

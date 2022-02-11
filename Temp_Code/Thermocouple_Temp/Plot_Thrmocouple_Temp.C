void Plot_Thrmocouple_Temp()
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(2);	
  gStyle->SetLabelOffset(0.007,"Y");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetLabelSize(0.04, "X");
  gStyle->SetNdivisions(512,"X");
  gStyle->SetTitleSize(0.05, "X" );
  gStyle->SetTitleSize(0.06, "Y" );
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetLabelFont(22,"X");
  gStyle->SetLabelFont(22,"Y");
  gStyle->SetTitleOffset(0.9,"X");
  gStyle->SetTitleOffset(0.40,"Y");
  
  TCanvas *plot = new TCanvas("plot","",1500,900);
  TPad *pad2 = new TPad("pad2", "",0.0,0.0,1.0,1.0);
  pad2->Draw();
  pad2->cd();
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(2);
  pad2->SetRightMargin(0.03);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.1);
  pad2->SetLeftMargin(0.09);
  pad2->SetFrameBorderMode(0);
  pad2->SetFrameBorderMode(0);
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  
  TH2F *frame2 = new TH2F("frame2","",10, -0.05, 500, 10, 10, 40);
  //TH2F *frame2 = new TH2F("frame2","",10, -0.05, 80, 10, 10, 40);
  frame2->GetXaxis()->SetTitle("Days");
  frame2->GetXaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitle("AC Temperature (^{o}C)");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitleSize(0.06);
  frame2->GetYaxis()->SetTitleOffset(0.6);
  frame2->Draw("");

    
  //Ge1500
  
  const int data = 2000;
  char err_1[data][20];
  
  int FILE_NO1[data],FILE_NO2[data],FILE_NO3[data],FILE_NO4[data],FILE_NO5[data],FILE_NO6[data],FILE_NO7[data],FILE_NO8[data],FILE_NO9[data];
  int FILE_NO10[data],FILE_NO11[data],FILE_NO12[data],FILE_NO13[data],FILE_NO14[data],FILE_NO15[data],FILE_NO16[data],FILE_NO17[data],FILE_NO18[data];
  int FILE_NO19[data],FILE_NO20[data],FILE_NO21[data],FILE_NO22[data],FILE_NO23[data],FILE_NO24[data];
  
  double T_Inner_Shield1[data], Broken_Channel1[data], T_Outer_Shield1[data], T_Controller1[data], Run_End_TIME1[data]; 
  double T_Inner_Shield2[data], Broken_Channel2[data], T_Outer_Shield2[data], T_Controller2[data], Run_End_TIME2[data];
  double T_Inner_Shield3[data], Broken_Channel3[data], T_Outer_Shield3[data], T_Controller3[data], Run_End_TIME3[data];
  double T_Inner_Shield4[data], Broken_Channel4[data], T_Outer_Shield4[data], T_Controller4[data], Run_End_TIME4[data];
  double T_Inner_Shield5[data], Broken_Channel5[data], T_Outer_Shield5[data], T_Controller5[data], Run_End_TIME5[data];
  double T_Inner_Shield6[data], Broken_Channel6[data], T_Outer_Shield6[data], T_Controller6[data], Run_End_TIME6[data];
  double T_Inner_Shield7[data], Broken_Channel7[data], T_Outer_Shield7[data], T_Controller7[data], Run_End_TIME7[data];
  double T_Inner_Shield8[data], Broken_Channel8[data], T_Outer_Shield8[data], T_Controller8[data], Run_End_TIME8[data];
  double T_Inner_Shield9[data], Broken_Channel9[data], T_Outer_Shield9[data], T_Controller9[data], Run_End_TIME9[data];
  double T_Inner_Shield10[data], Broken_Channel10[data], T_Outer_Shield10[data], T_Controller10[data], Run_End_TIME10[data];
  double T_Inner_Shield11[data], Broken_Channel11[data], T_Outer_Shield11[data], T_Controller11[data], Run_End_TIME11[data];
  double T_Inner_Shield12[data], Broken_Channel12[data], T_Outer_Shield12[data], T_Controller12[data], Run_End_TIME12[data];
  double T_Inner_Shield13[data], Broken_Channel13[data], T_Outer_Shield13[data], T_Controller13[data], Run_End_TIME13[data];
  double T_Inner_Shield14[data], Broken_Channel14[data], T_Outer_Shield14[data], T_Controller14[data], Run_End_TIME14[data];
  double T_Inner_Shield15[data], Broken_Channel15[data], T_Outer_Shield15[data], T_Controller15[data], Run_End_TIME15[data];
  double T_Inner_Shield16[data], Broken_Channel16[data], T_Outer_Shield16[data], T_Controller16[data], Run_End_TIME16[data];
  double T_Inner_Shield17[data], Broken_Channel17[data], T_Outer_Shield17[data], T_Controller17[data], Run_End_TIME17[data];
  double T_Inner_Shield18[data], Broken_Channel18[data], T_Outer_Shield18[data], T_Controller18[data], Run_End_TIME18[data];
  double T_Inner_Shield19[data], Broken_Channel19[data], T_Outer_Shield19[data], T_Controller19[data], Run_End_TIME19[data];
  double T_Inner_Shield20[data], Broken_Channel20[data], T_Outer_Shield20[data], T_Controller20[data], Run_End_TIME20[data];
  double T_Inner_Shield21[data], Broken_Channel21[data], T_Outer_Shield21[data], T_Controller21[data], Run_End_TIME21[data];
  double T_Inner_Shield22[data], Broken_Channel22[data], T_Outer_Shield22[data], T_Controller22[data], Run_End_TIME22[data];
  double T_Inner_Shield23[data], Broken_Channel23[data], T_Outer_Shield23[data], T_Controller23[data], Run_End_TIME23[data];
  double T_Inner_Shield24[data], Broken_Channel24[data], T_Outer_Shield24[data], T_Controller24[data], Run_End_TIME24[data];


  
  double Xaxis_1[data],Xaxis_2[data],Xaxis_3[data],Xaxis_4[data],Xaxis_5[data],Xaxis_6[data],Xaxis_7[data],Xaxis_8[data],Xaxis_9[data],Xaxis_10[data],Xaxis_11[data];
  double Xaxis_12[data],Xaxis_13[data],Xaxis_14[data],Xaxis_15[data],Xaxis_16[data],Xaxis_17[data],Xaxis_18[data],Xaxis_19[data],Xaxis_20[data];
  double Xaxis_21[data],Xaxis_22[data],Xaxis_23[data],Xaxis_24[data];
  
  FILE *file1 =fopen("200911_1753_PhysOn_Det70x70NaI3ThrdP685_Det50x50NaI5Thrd0_FADC_RAW_Data.txt","r");
  FILE *file2 =fopen("200925_2042_PhysOn_Det70x70_3020HV_200C_114SA_Thrd685_Det50x50NaI5Thrd0_FADC_RAW_Data.txt","r");
  FILE *file3 =fopen("201102_1158_PhysOn_Det70x70_3020HV_Tri_Out1_2Copy_No_Out2_ThrdP690_Det50x50NaI5Thrd0_FADC_RAW_Data.txt","r");
  FILE *file4 =fopen("201109_1412_PhysOn_Det70x70_3020HV_Tri_newDataSettings_ThrdM490_Det50x50NaI5Thrd0_FADC_RAW_Data.txt","r");
  FILE *file5 =fopen("201116_1914_PhysOn_Det70x70_3020HV_Tri_newDataSettings_ThrdM490_Det50x50NaI5Thrd0_FADC_RAW_Data.txt","r");
  FILE *file6 =fopen("201117_1320_PhysOn_Det70x70_3020HV_Tri_onlyOneouterSA_ThrdM450_Det50x50NaI5Thrd0_FADC_RAW_Data.txt","r");
  FILE *file7 =fopen("201119_1316_PhysOn_Det70x70_3020HV_Tri_out1SALE_out2SA_n_SABacktoTA_ThrdM455_DetFADC_RAW_Data.txt","r");
  FILE *file8 =fopen("201201_1346_PhysOn_Det70x70_SA_3Ch_TA_1Ch_ThrdP2295_Det50x50_OldSett_ThrdM560FADC_RAW_Data.txt","r");
  FILE *file9 =fopen("201231_1200_PhysOn_Det70x70_Out1_SA_Ch01_Out2_Copies_SAFastCh2_TACh0_HE60M3_ThrdM522_Det50x50_OldSett_SApt5usCh5onwrds_ThrdM568_FADC_RAW_Data.txt","r");
  FILE *file10 =fopen("210125_1411_PhysOn_PrePostTrg_Det70x70_ThrdM532_Det50x50_ThrdM570_FADC_RAW_Data.txt","r");
  FILE *file11 =fopen("210226_1257_PhysOn_Det70x70_ThrdM529_Det50x50_ThrdM566_FADC_RAW_Data.txt","r");
  FILE *file12 =fopen("210331_1354_PhysOn_Det70x70_ThrdM526_Det50x50_ThrdM800_FADC_RAW_Data.txt","r");
  FILE *file13 =fopen("210424_1445_PhysOn_Det70x70_ThrdM520_Det50x50_ThrdM570_FADC_RAW_Data.txt","r");
  FILE *file14 =fopen("210527_1200_PhysOn_Det70x70_ThrdM503_Det50x50_ThrdM556_FADC_RAW_Data.txt","r");
  FILE *file15 =fopen("210630_1010_PhysOn_Det70x70_ThrdM504_Det50x50_ThrdM558_FADC_RAW_Data.txt","r");
  FILE *file16 =fopen("210727_1228_HiThrd_PhysOFF_Det70x70_ThrdM485_Det50x50_ThrdM535_FADC_RAW_Data.txt","r");
  FILE *file17 =fopen("210827_1555_PhysOFF_Det70x70_ThrdM520_Det50x50_ThrdM570_FADC_RAW_Data.txt","r");
  FILE *file18 =fopen("210925_1313_PhysOFF_Det70x70_ThrdM520_Det50x50_ThrdM575_FADC_RAW_Data.txt","r");
  FILE *file19 =fopen("211027_1051_PhysOFF_Det70x70_ThrdM518_Det50x50_ThrdM573_FADC_RAW_Data.txt","r");
  FILE *file20 =fopen("211130_1129_PhysOFF_Det70x70_ThrdM515_Det50x50_ThrdM570_FADC_RAW_Data.txt","r");
  FILE *file21 =fopen("211207_1238_PhysOFF_Det70x70_ThrdM515_Det50x50_ThrdM570_FADC_RAW_Data.txt","r");
  FILE *file22 =fopen("211228_1159_PhysOFF_Det70x70_ThrdM510_Det50x50_ThrdM557_FADC_RAW_Data.txt","r");
  FILE *file23 =fopen("211230_1453_PhysOFF_Det70x70_ThrdM522_Det50x50_Temp188_ThrdM565_FADC_RAW_Data.txt","r");
  FILE *file24 =fopen("220105_1221_PhysOFF_Det70x70_ThrdM522_Det50x50_Temp188_ThrdM565_FADC_RAW_Data.txt","r");
  
  for(int ll =0; ll<data; ll++)
    {
      fscanf(file1, "%d %lf %lf %lf %lf %lf \n", &FILE_NO1[ll], &T_Inner_Shield1[ll], &Broken_Channel1[ll], &T_Outer_Shield1[ll], &T_Controller1[ll], &Run_End_TIME1[ll]);
      fscanf(file2, "%d %lf %lf %lf %lf %lf \n", &FILE_NO2[ll], &T_Inner_Shield2[ll], &Broken_Channel2[ll], &T_Outer_Shield2[ll], &T_Controller2[ll], &Run_End_TIME2[ll]);
      fscanf(file3, "%d %lf %lf %lf %lf %lf \n", &FILE_NO3[ll], &T_Inner_Shield3[ll], &Broken_Channel3[ll], &T_Outer_Shield3[ll], &T_Controller3[ll], &Run_End_TIME3[ll]);
      fscanf(file4, "%d %lf %lf %lf %lf %lf \n", &FILE_NO4[ll], &T_Inner_Shield4[ll], &Broken_Channel4[ll], &T_Outer_Shield4[ll], &T_Controller4[ll], &Run_End_TIME4[ll]);
      fscanf(file5, "%d %lf %lf %lf %lf %lf \n", &FILE_NO5[ll], &T_Inner_Shield5[ll], &Broken_Channel5[ll], &T_Outer_Shield5[ll], &T_Controller5[ll], &Run_End_TIME5[ll]);
      fscanf(file6, "%d %lf %lf %lf %lf %lf \n", &FILE_NO6[ll], &T_Inner_Shield6[ll], &Broken_Channel6[ll], &T_Outer_Shield6[ll], &T_Controller6[ll], &Run_End_TIME6[ll]);
      fscanf(file7, "%d %lf %lf %lf %lf %lf \n", &FILE_NO7[ll], &T_Inner_Shield7[ll], &Broken_Channel7[ll], &T_Outer_Shield7[ll], &T_Controller7[ll], &Run_End_TIME7[ll]);
      fscanf(file8, "%d %lf %lf %lf %lf %lf \n", &FILE_NO8[ll], &T_Inner_Shield8[ll], &Broken_Channel8[ll], &T_Outer_Shield8[ll], &T_Controller8[ll], &Run_End_TIME8[ll]);
      fscanf(file9, "%d %lf %lf %lf %lf %lf \n", &FILE_NO9[ll], &T_Inner_Shield9[ll], &Broken_Channel9[ll], &T_Outer_Shield9[ll], &T_Controller9[ll], &Run_End_TIME9[ll]);
      fscanf(file10, "%d %lf %lf %lf %lf %lf \n", &FILE_NO10[ll], &T_Inner_Shield10[ll], &Broken_Channel10[ll], &T_Outer_Shield10[ll], &T_Controller10[ll], &Run_End_TIME10[ll]);
      fscanf(file11, "%d %lf %lf %lf %lf %lf \n", &FILE_NO11[ll], &T_Inner_Shield11[ll], &Broken_Channel11[ll], &T_Outer_Shield11[ll], &T_Controller11[ll], &Run_End_TIME11[ll]);
      fscanf(file12, "%d %lf %lf %lf %lf %lf \n", &FILE_NO12[ll], &T_Inner_Shield12[ll], &Broken_Channel12[ll], &T_Outer_Shield12[ll], &T_Controller12[ll], &Run_End_TIME12[ll]);
      fscanf(file13, "%d %lf %lf %lf %lf %lf \n", &FILE_NO13[ll], &T_Inner_Shield13[ll], &Broken_Channel13[ll], &T_Outer_Shield13[ll], &T_Controller13[ll], &Run_End_TIME13[ll]);
      fscanf(file14, "%d %lf %lf %lf %lf %lf \n", &FILE_NO14[ll], &T_Inner_Shield14[ll], &Broken_Channel14[ll], &T_Outer_Shield14[ll], &T_Controller14[ll], &Run_End_TIME14[ll]);
      fscanf(file15, "%d %lf %lf %lf %lf %lf \n", &FILE_NO15[ll], &T_Inner_Shield15[ll], &Broken_Channel15[ll], &T_Outer_Shield15[ll], &T_Controller15[ll], &Run_End_TIME15[ll]);
      fscanf(file16, "%d %lf %lf %lf %lf %lf \n", &FILE_NO16[ll], &T_Inner_Shield16[ll], &Broken_Channel16[ll], &T_Outer_Shield16[ll], &T_Controller16[ll], &Run_End_TIME16[ll]);
      fscanf(file17, "%d %lf %lf %lf %lf %lf \n", &FILE_NO17[ll], &T_Inner_Shield17[ll], &Broken_Channel17[ll], &T_Outer_Shield17[ll], &T_Controller17[ll], &Run_End_TIME17[ll]);
      fscanf(file18, "%d %lf %lf %lf %lf %lf \n", &FILE_NO18[ll], &T_Inner_Shield18[ll], &Broken_Channel18[ll], &T_Outer_Shield18[ll], &T_Controller18[ll], &Run_End_TIME18[ll]);
      fscanf(file19, "%d %lf %lf %lf %lf %lf \n", &FILE_NO19[ll], &T_Inner_Shield19[ll], &Broken_Channel19[ll], &T_Outer_Shield19[ll], &T_Controller19[ll], &Run_End_TIME19[ll]);
      fscanf(file20, "%d %lf %lf %lf %lf %lf \n", &FILE_NO20[ll], &T_Inner_Shield20[ll], &Broken_Channel20[ll], &T_Outer_Shield20[ll], &T_Controller20[ll], &Run_End_TIME20[ll]);
      fscanf(file21, "%d %lf %lf %lf %lf %lf \n", &FILE_NO21[ll], &T_Inner_Shield21[ll], &Broken_Channel21[ll], &T_Outer_Shield21[ll], &T_Controller21[ll], &Run_End_TIME21[ll]);
      fscanf(file22, "%d %lf %lf %lf %lf %lf \n", &FILE_NO22[ll], &T_Inner_Shield22[ll], &Broken_Channel22[ll], &T_Outer_Shield22[ll], &T_Controller22[ll], &Run_End_TIME22[ll]);
      fscanf(file23, "%d %lf %lf %lf %lf %lf \n", &FILE_NO23[ll], &T_Inner_Shield23[ll], &Broken_Channel23[ll], &T_Outer_Shield23[ll], &T_Controller23[ll], &Run_End_TIME23[ll]);
      fscanf(file24, "%d %lf %lf %lf %lf %lf \n", &FILE_NO24[ll], &T_Inner_Shield24[ll], &Broken_Channel24[ll], &T_Outer_Shield24[ll], &T_Controller24[ll], &Run_End_TIME24[ll]);

      
      if(FILE_NO1[ll]==1)
       	{
       	  Xaxis_1[ll] = 0.0;
	}
      else
       	{
	  // Xaxis_1[ll] = ((Run_End_TIME[ll]-Run_End_TIME[ll-1]));   // 
	  Xaxis_1[ll] = ((Run_End_TIME1[ll]-Run_End_TIME1[0])/(60.0*60.0*24.0));                     
	}

      if(FILE_NO2[ll]==1)
       	{
       	  Xaxis_2[ll] = 15+0.0;
	}
      else
       	{
	  Xaxis_2[ll] = 15.0 + ((Run_End_TIME2[ll]-Run_End_TIME2[0])/(60.0*60.0*24.0));                     
	}

      if(FILE_NO3[ll]==1)
       	{
       	  Xaxis_3[ll] = 37.0+15+0.0;
	}
      else
       	{
	  Xaxis_3[ll] = 37.0+15.0 + ((Run_End_TIME3[ll]-Run_End_TIME3[0])/(60.0*60.0*24.0));                     
	}
      
      if(FILE_NO4[ll]==1)
       	{
       	  Xaxis_4[ll] = 7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_4[ll] = 7.0+37.0+15.0 + ((Run_End_TIME4[ll]-Run_End_TIME4[0])/(60.0*60.0*24.0));                     
	}

      if(FILE_NO5[ll]==1)
       	{
       	  Xaxis_5[ll] = 7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_5[ll] = 7.0+7.0+37.0+15.0 + ((Run_End_TIME5[ll]-Run_End_TIME5[0])/(60.0*60.0*24.0));                     
	}

      if(FILE_NO6[ll]==1)
       	{
       	  Xaxis_6[ll] = 1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_6[ll] = 1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME6[ll]-Run_End_TIME6[0])/(60.0*60.0*24.0));                     
	}
      
      if(FILE_NO7[ll]==1)
       	{
       	  Xaxis_7[ll] = 2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_7[ll] = 2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME7[ll]-Run_End_TIME7[0])/(60.0*60.0*24.0));                     
	}
      
      if(FILE_NO8[ll]==1)
       	{
       	  Xaxis_8[ll] = 12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_8[ll] = 12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME8[ll]-Run_End_TIME8[0])/(60.0*60.0*24.0));                     
	}
      
      if(FILE_NO9[ll]==1)
       	{
       	  Xaxis_9[ll] = 31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_9[ll] = 31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME9[ll]-Run_End_TIME9[0])/(60.0*60.0*24.0));                     
	}
      
      if(FILE_NO10[ll]==1)
       	{
       	  Xaxis_10[ll] = 30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_10[ll] = 30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME10[ll]-Run_End_TIME10[0])/(60.0*60.0*24.0));                     
	}

      if(FILE_NO11[ll]==1)
       	{
       	  Xaxis_11[ll] = 30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_11[ll] = 30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME11[ll]-Run_End_TIME11[0])/(60.0*60.0*24.0));                     
	}
      
       if(FILE_NO12[ll]==1)
       	{
       	  Xaxis_12[ll] = 35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_12[ll] = 35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME12[ll]-Run_End_TIME12[0])/(60.0*60.0*24.0));                     
	}
       
       if(FILE_NO13[ll]==1)
       	{
       	  Xaxis_13[ll] = 30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_13[ll] = 30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME13[ll]-Run_End_TIME13[0])/(60.0*60.0*24.0));                     
	}
       
	  
       if(FILE_NO14[ll]==1)
       	{
       	  Xaxis_14[ll] = 30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
      else
       	{
	  Xaxis_14[ll] = 30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME14[ll]-Run_End_TIME14[0])/(60.0*60.0*24.0));                     
	}

       if(FILE_NO15[ll]==1)
       	{
       	  Xaxis_15[ll] = 30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	}
       else
	 {
	   Xaxis_15[ll] = 30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME15[ll]-Run_End_TIME15[0])/(60.0*60.0*24.0));                     
	 }
      
       if(FILE_NO16[ll]==1)
	 {
	   Xaxis_16[ll] = 30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_16[ll] = 30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME16[ll]-Run_End_TIME16[0])/(60.0*60.0*24.0));                     
	 }

       if(FILE_NO17[ll]==1)
	 {
	   Xaxis_17[ll] = 30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_17[ll] = 30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME17[ll]-Run_End_TIME17[0])/(60.0*60.0*24.0));                     
	 }

       if(FILE_NO18[ll]==1)
	 {
	   Xaxis_18[ll] = 30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_18[ll] = 30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME18[ll]-Run_End_TIME18[0])/(60.0*60.0*24.0));                     
	 }

       if(FILE_NO19[ll]==1)
	 {
	   Xaxis_19[ll] = 30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_19[ll] = 30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME19[ll]-Run_End_TIME19[0])/(60.0*60.0*24.0));                     
	 }
       
       if(FILE_NO20[ll]==1)
	 {
	   Xaxis_20[ll] = 30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_20[ll] = 30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME20[ll]-Run_End_TIME20[0])/(60.0*60.0*24.0));                     
	 }

       if(FILE_NO21[ll]==1)
	 {
	   Xaxis_21[ll] = 7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_21[ll] = 7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME21[ll]-Run_End_TIME21[0])/(60.0*60.0*24.0));                     
	 }

       if(FILE_NO22[ll]==1)
	 {
	   Xaxis_22[ll] = 23.0+7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_22[ll] = 23.0+7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME22[ll]-Run_End_TIME22[0])/(60.0*60.0*24.0));                     
	 }
       
       if(FILE_NO23[ll]==1)
	 {
	   Xaxis_23[ll] = 6.0+7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_23[ll] = 6.0+7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME23[ll]-Run_End_TIME23[0])/(60.0*60.0*24.0));                     
	 }

       if(FILE_NO23[ll]==1)
	 {
	   Xaxis_23[ll] = 3.0+7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15+0.0;
	 }
       else
	 {
	   Xaxis_23[ll] = 3.0+7.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+30.0+35.0+30.0+30.0+31.0+12.0+2.0+1.0+7.0+7.0+37.0+15.0 + ((Run_End_TIME23[ll]-Run_End_TIME23[0])/(60.0*60.0*24.0));                     
	 }


       
       // printf( "%d  %f  %f  %f  %f  %f  %f\n", FILE_NO1[ll], T_Inner_Shield1[ll], Broken_Channel1[ll], T_Outer_Shield1[ll], T_Controller1[ll], Run_End_TIME1[ll], Xaxis_1[ll]);
       printf( "%d  %f  %f  %f  %f  %f  %f\n", FILE_NO2[ll], T_Inner_Shield2[ll], Broken_Channel2[ll], T_Outer_Shield2[ll], T_Controller2[ll], Run_End_TIME2[ll], Xaxis_2[ll]);
       
      
    }
  printf("file 1st Over\n");


   TGraph *gr_11 = new TGraph(243, Xaxis_1, T_Inner_Shield1);
   gr_11->SetLineColor(1);
   gr_11->SetLineWidth(3);
   gr_11->SetMarkerColor(1);
   gr_11->SetMarkerSize(3);
   gr_11->SetMarkerStyle(21);
   gr_11->Draw("p");

   TGraph *gr_22 = new TGraph(450, Xaxis_2, T_Inner_Shield2);
   gr_22->SetLineColor(1);
   gr_22->SetLineWidth(3);
   gr_22->SetMarkerColor(1);
   gr_22->SetMarkerSize(3);
   gr_22->SetMarkerStyle(21);
   gr_22->Draw("p");
   
   TGraph *gr_33 = new TGraph(289, Xaxis_3, T_Inner_Shield3);
   gr_33->SetLineColor(1);
   gr_33->SetLineWidth(3);
   gr_33->SetMarkerColor(1);
   gr_33->SetMarkerSize(3);
   gr_33->SetMarkerStyle(21);
   gr_33->Draw("p");
   
   TGraph *gr_44 = new TGraph(450, Xaxis_4, T_Inner_Shield4);
   gr_44->SetLineColor(1);
   gr_44->SetLineWidth(3);
   gr_44->SetMarkerColor(1);
   gr_44->SetMarkerSize(3);
   gr_44->SetMarkerStyle(21);
   gr_44->Draw("p");

   TGraph *gr_55 = new TGraph(450, Xaxis_5, T_Inner_Shield5);
   gr_55->SetLineColor(1);
   gr_55->SetLineWidth(3);
   gr_55->SetMarkerColor(1);
   gr_55->SetMarkerSize(3);
   gr_55->SetMarkerStyle(21);
   gr_55->Draw("p");
   
   TGraph *gr_66 = new TGraph(450, Xaxis_6, T_Inner_Shield6);
   gr_66->SetLineColor(1);
   gr_66->SetLineWidth(3);
   gr_66->SetMarkerColor(1);
   gr_66->SetMarkerSize(3);
   gr_66->SetMarkerStyle(21);
   gr_66->Draw("p");
   
   TGraph *gr_77 = new TGraph(450, Xaxis_7, T_Inner_Shield7);
   gr_77->SetLineColor(1);
   gr_77->SetLineWidth(3);
   gr_77->SetMarkerColor(1);
   gr_77->SetMarkerSize(3);
   gr_77->SetMarkerStyle(21);
   gr_77->Draw("p");
   
   TGraph *gr_88 = new TGraph(450, Xaxis_8, T_Inner_Shield8);
   gr_88->SetLineColor(1);
   gr_88->SetLineWidth(3);
   gr_88->SetMarkerColor(1);
   gr_88->SetMarkerSize(3);
   gr_88->SetMarkerStyle(21);
   gr_88->Draw("p");

   TGraph *gr_99 = new TGraph(450, Xaxis_9, T_Inner_Shield9);
   gr_99->SetLineColor(1);
   gr_99->SetLineWidth(3);
   gr_99->SetMarkerColor(1);
   gr_99->SetMarkerSize(3);
   gr_99->SetMarkerStyle(21);
   gr_99->Draw("p");
   
   TGraph *gr_1010 = new TGraph(450, Xaxis_10, T_Inner_Shield10);
   gr_1010->SetLineColor(1);
   gr_1010->SetLineWidth(3);
   gr_1010->SetMarkerColor(1);
   gr_1010->SetMarkerSize(3);
   gr_1010->SetMarkerStyle(21);
   gr_1010->Draw("p");

   TGraph *gr_1111 = new TGraph(450, Xaxis_11, T_Inner_Shield11);
   gr_1111->SetLineColor(1);
   gr_1111->SetLineWidth(3);
   gr_1111->SetMarkerColor(1);
   gr_1111->SetMarkerSize(3);
   gr_1111->SetMarkerStyle(21);
   gr_1111->Draw("p");
   
   TGraph *gr_1212 = new TGraph(450, Xaxis_12, T_Inner_Shield12);
   gr_1212->SetLineColor(1);
   gr_1212->SetLineWidth(3);
   gr_1212->SetMarkerColor(1);
   gr_1212->SetMarkerSize(3);
   gr_1212->SetMarkerStyle(21);
   gr_1212->Draw("p");
   
   TGraph *gr_1313 = new TGraph(450, Xaxis_13, T_Inner_Shield13);
   gr_1313->SetLineColor(1);
   gr_1313->SetLineWidth(3);
   gr_1313->SetMarkerColor(1);
   gr_1313->SetMarkerSize(3);
   gr_1313->SetMarkerStyle(21);
   gr_1313->Draw("p");
   
   TGraph *gr_1414 = new TGraph(450, Xaxis_14, T_Inner_Shield14);
   gr_1414->SetLineColor(1);
   gr_1414->SetLineWidth(3);
   gr_1414->SetMarkerColor(1);
   gr_1414->SetMarkerSize(3);
   gr_1414->SetMarkerStyle(21);
   gr_1414->Draw("p");

   TGraph *gr_1515 = new TGraph(450, Xaxis_15, T_Inner_Shield15);
   gr_1515->SetLineColor(1);
   gr_1515->SetLineWidth(3);
   gr_1515->SetMarkerColor(1);
   gr_1515->SetMarkerSize(3);
   gr_1515->SetMarkerStyle(21);
   gr_1515->Draw("p");

   TGraph *gr_1616 = new TGraph(450, Xaxis_16, T_Inner_Shield16);
   gr_1616->SetLineColor(1);
   gr_1616->SetLineWidth(3);
   gr_1616->SetMarkerColor(1);
   gr_1616->SetMarkerSize(3);
   gr_1616->SetMarkerStyle(21);
   gr_1616->Draw("p");

   TGraph *gr_1717 = new TGraph(450, Xaxis_17, T_Inner_Shield17);
   gr_1717->SetLineColor(1);
   gr_1717->SetLineWidth(3);
   gr_1717->SetMarkerColor(1);
   gr_1717->SetMarkerSize(3);
   gr_1717->SetMarkerStyle(21);
   gr_1717->Draw("p");

   TGraph *gr_1818 = new TGraph(450, Xaxis_18, T_Inner_Shield18);
   gr_1818->SetLineColor(1);
   gr_1818->SetLineWidth(3);
   gr_1818->SetMarkerColor(1);
   gr_1818->SetMarkerSize(3);
   gr_1818->SetMarkerStyle(21);
   gr_1818->Draw("p");

   TGraph *gr_1919 = new TGraph(450, Xaxis_19, T_Inner_Shield19);
   gr_1919->SetLineColor(1);
   gr_1919->SetLineWidth(3);
   gr_1919->SetMarkerColor(1);
   gr_1919->SetMarkerSize(3);
   gr_1919->SetMarkerStyle(21);
   gr_1919->Draw("p");

   TGraph *gr_2020 = new TGraph(450, Xaxis_20, T_Inner_Shield20);
   gr_2020->SetLineColor(1);
   gr_2020->SetLineWidth(3);
   gr_2020->SetMarkerColor(1);
   gr_2020->SetMarkerSize(3);
   gr_2020->SetMarkerStyle(21);
   gr_2020->Draw("p");

   TGraph *gr_2021 = new TGraph(450, Xaxis_21, T_Inner_Shield21);
   gr_2021->SetLineColor(1);
   gr_2021->SetLineWidth(3);
   gr_2021->SetMarkerColor(1);
   gr_2021->SetMarkerSize(3);
   gr_2021->SetMarkerStyle(21);
   gr_2021->Draw("p");
   
   TGraph *gr_2022 = new TGraph(450, Xaxis_22, T_Inner_Shield22);
   gr_2022->SetLineColor(1);
   gr_2022->SetLineWidth(3);
   gr_2022->SetMarkerColor(1);
   gr_2022->SetMarkerSize(3);
   gr_2022->SetMarkerStyle(21);
   gr_2022->Draw("p");
   
   TGraph *gr_2023 = new TGraph(450, Xaxis_23, T_Inner_Shield23);
   gr_2023->SetLineColor(1);
   gr_2023->SetLineWidth(3);
   gr_2023->SetMarkerColor(1);
   gr_2023->SetMarkerSize(3);
   gr_2023->SetMarkerStyle(21);
   gr_2023->Draw("p");
   
   TGraph *gr_2024 = new TGraph(450, Xaxis_24, T_Inner_Shield24);
   gr_2024->SetLineColor(1);
   gr_2024->SetLineWidth(3);
   gr_2024->SetMarkerColor(1);
   gr_2024->SetMarkerSize(3);
   gr_2024->SetMarkerStyle(21);
   gr_2024->Draw("p");
   
   //======================== Controller =========================================


   TGraph *CR_11 = new TGraph(243, Xaxis_1, T_Controller1);
   CR_11->SetLineColor(kGreen+2);
   CR_11->SetLineWidth(3);
   CR_11->SetMarkerColor(kGreen+2);
   CR_11->SetMarkerSize(3);
   CR_11->SetMarkerStyle(21);
   CR_11->Draw("p");

   TGraph *CR_22 = new TGraph(450, Xaxis_2, T_Controller2);
   CR_22->SetLineColor(kGreen+2);
   CR_22->SetLineWidth(3);
   CR_22->SetMarkerColor(kGreen+2);
   CR_22->SetMarkerSize(3);
   CR_22->SetMarkerStyle(21);
   CR_22->Draw("p");
   
   TGraph *CR_33 = new TGraph(289, Xaxis_3, T_Controller3);
   CR_33->SetLineColor(kGreen+2);
   CR_33->SetLineWidth(3);
   CR_33->SetMarkerColor(kGreen+2);
   CR_33->SetMarkerSize(3);
   CR_33->SetMarkerStyle(21);
   CR_33->Draw("p");
   
   TGraph *CR_44 = new TGraph(450, Xaxis_4, T_Controller4);
   CR_44->SetLineColor(kGreen+2);
   CR_44->SetLineWidth(3);
   CR_44->SetMarkerColor(kGreen+2);
   CR_44->SetMarkerSize(3);
   CR_44->SetMarkerStyle(21);
   CR_44->Draw("p");

   TGraph *CR_55 = new TGraph(450, Xaxis_5, T_Controller5);
   CR_55->SetLineColor(kGreen+2);
   CR_55->SetLineWidth(3);
   CR_55->SetMarkerColor(kGreen+2);
   CR_55->SetMarkerSize(3);
   CR_55->SetMarkerStyle(21);
   CR_55->Draw("p");
   
   TGraph *CR_66 = new TGraph(450, Xaxis_6, T_Controller6);
   CR_66->SetLineColor(kGreen+2);
   CR_66->SetLineWidth(3);
   CR_66->SetMarkerColor(kGreen+2);
   CR_66->SetMarkerSize(3);
   CR_66->SetMarkerStyle(21);
   CR_66->Draw("p");
   
   TGraph *CR_77 = new TGraph(450, Xaxis_7, T_Controller7);
   CR_77->SetLineColor(kGreen+2);
   CR_77->SetLineWidth(3);
   CR_77->SetMarkerColor(kGreen+2);
   CR_77->SetMarkerSize(3);
   CR_77->SetMarkerStyle(21);
   CR_77->Draw("p");
   
   TGraph *CR_88 = new TGraph(450, Xaxis_8, T_Controller8);
   CR_88->SetLineColor(kGreen+2);
   CR_88->SetLineWidth(3);
   CR_88->SetMarkerColor(kGreen+2);
   CR_88->SetMarkerSize(3);
   CR_88->SetMarkerStyle(21);
   CR_88->Draw("p");

   TGraph *CR_99 = new TGraph(450, Xaxis_9, T_Controller9);
   CR_99->SetLineColor(kGreen+2);
   CR_99->SetLineWidth(3);
   CR_99->SetMarkerColor(kGreen+2);
   CR_99->SetMarkerSize(3);
   CR_99->SetMarkerStyle(21);
   CR_99->Draw("p");
   
   TGraph *CR_1010 = new TGraph(450, Xaxis_10, T_Controller10);
   CR_1010->SetLineColor(kGreen+2);
   CR_1010->SetLineWidth(3);
   CR_1010->SetMarkerColor(kGreen+2);
   CR_1010->SetMarkerSize(3);
   CR_1010->SetMarkerStyle(21);
   CR_1010->Draw("p");

   TGraph *CR_1111 = new TGraph(450, Xaxis_11, T_Controller11);
   CR_1111->SetLineColor(kGreen+2);
   CR_1111->SetLineWidth(3);
   CR_1111->SetMarkerColor(kGreen+2);
   CR_1111->SetMarkerSize(3);
   CR_1111->SetMarkerStyle(21);
   CR_1111->Draw("p");
   
   TGraph *CR_1212 = new TGraph(450, Xaxis_12, T_Controller12);
   CR_1212->SetLineColor(kGreen+2);
   CR_1212->SetLineWidth(3);
   CR_1212->SetMarkerColor(kGreen+2);
   CR_1212->SetMarkerSize(3);
   CR_1212->SetMarkerStyle(21);
   CR_1212->Draw("p");
   
   TGraph *CR_1313 = new TGraph(450, Xaxis_13, T_Controller13);
   CR_1313->SetLineColor(kGreen+2);
   CR_1313->SetLineWidth(3);
   CR_1313->SetMarkerColor(kGreen+2);
   CR_1313->SetMarkerSize(3);
   CR_1313->SetMarkerStyle(21);
   CR_1313->Draw("p");
   
   TGraph *CR_1414 = new TGraph(450, Xaxis_14, T_Controller14);
   CR_1414->SetLineColor(kGreen+2);
   CR_1414->SetLineWidth(3);
   CR_1414->SetMarkerColor(kGreen+2);
   CR_1414->SetMarkerSize(3);
   CR_1414->SetMarkerStyle(21);
   CR_1414->Draw("p");

   TGraph *CR_1515 = new TGraph(450, Xaxis_15, T_Controller15);
   CR_1515->SetLineColor(kGreen+2);
   CR_1515->SetLineWidth(3);
   CR_1515->SetMarkerColor(kGreen+2);
   CR_1515->SetMarkerSize(3);
   CR_1515->SetMarkerStyle(21);
   CR_1515->Draw("p");

   TGraph *CR_1616 = new TGraph(450, Xaxis_16, T_Controller16);
   CR_1616->SetLineColor(kGreen+2);
   CR_1616->SetLineWidth(3);
   CR_1616->SetMarkerColor(kGreen+2);
   CR_1616->SetMarkerSize(3);
   CR_1616->SetMarkerStyle(21);
   CR_1616->Draw("p");

   TGraph *CR_1717 = new TGraph(450, Xaxis_17, T_Controller17);
   CR_1717->SetLineColor(kGreen+2);
   CR_1717->SetLineWidth(3);
   CR_1717->SetMarkerColor(kGreen+2);
   CR_1717->SetMarkerSize(3);
   CR_1717->SetMarkerStyle(21);
   CR_1717->Draw("p");

   TGraph *CR_1818 = new TGraph(450, Xaxis_18, T_Controller18);
   CR_1818->SetLineColor(kGreen+2);
   CR_1818->SetLineWidth(3);
   CR_1818->SetMarkerColor(kGreen+2);
   CR_1818->SetMarkerSize(3);
   CR_1818->SetMarkerStyle(21);
   CR_1818->Draw("p");

   TGraph *CR_1919 = new TGraph(450, Xaxis_19, T_Controller19);
   CR_1919->SetLineColor(kGreen+2);
   CR_1919->SetLineWidth(3);
   CR_1919->SetMarkerColor(kGreen+2);
   CR_1919->SetMarkerSize(3);
   CR_1919->SetMarkerStyle(21);
   CR_1919->Draw("p");

   TGraph *CR_2020 = new TGraph(450, Xaxis_20, T_Controller20);
   CR_2020->SetLineColor(kGreen+2);
   CR_2020->SetLineWidth(3);
   CR_2020->SetMarkerColor(kGreen+2);
   CR_2020->SetMarkerSize(3);
   CR_2020->SetMarkerStyle(21);
   CR_2020->Draw("p");
   
   TGraph *CR_2021 = new TGraph(450, Xaxis_21, T_Controller21);
   CR_2021->SetLineColor(kGreen+2);
   CR_2021->SetLineWidth(3);
   CR_2021->SetMarkerColor(kGreen+2);
   CR_2021->SetMarkerSize(3);
   CR_2021->SetMarkerStyle(21);
   CR_2021->Draw("p");
   
   TGraph *CR_2022 = new TGraph(450, Xaxis_22, T_Controller22);
   CR_2022->SetLineColor(kGreen+2);
   CR_2022->SetLineWidth(3);
   CR_2022->SetMarkerColor(kGreen+2);
   CR_2022->SetMarkerSize(3);
   CR_2022->SetMarkerStyle(21);
   CR_2022->Draw("p");
   
   TGraph *CR_2023 = new TGraph(450, Xaxis_23, T_Controller23);
   CR_2023->SetLineColor(kGreen+2);
   CR_2023->SetLineWidth(3);
   CR_2023->SetMarkerColor(kGreen+2);
   CR_2023->SetMarkerSize(3);
   CR_2023->SetMarkerStyle(21);
   CR_2023->Draw("p");
   
   TGraph *CR_2024 = new TGraph(450, Xaxis_24, T_Controller24);
   CR_2024->SetLineColor(kGreen+2);
   CR_2024->SetLineWidth(3);
   CR_2024->SetMarkerColor(kGreen+2);
   CR_2024->SetMarkerSize(3);
   CR_2024->SetMarkerStyle(21);
   CR_2024->Draw("p");
   
   //=======================================Outer Shielding ======================

   
   TGraph *OUS_11 = new TGraph(243, Xaxis_1, T_Outer_Shield1);
   OUS_11->SetLineColor(4);
   OUS_11->SetLineWidth(3);
   OUS_11->SetMarkerColor(4);
   OUS_11->SetMarkerSize(3);
   OUS_11->SetMarkerStyle(21);
   OUS_11->Draw("p");

   TGraph *OUS_22 = new TGraph(450, Xaxis_2, T_Outer_Shield2);
   OUS_22->SetLineColor(4);
   OUS_22->SetLineWidth(3);
   OUS_22->SetMarkerColor(4);
   OUS_22->SetMarkerSize(3);
   OUS_22->SetMarkerStyle(21);
   OUS_22->Draw("p");
   
   TGraph *OUS_33 = new TGraph(289, Xaxis_3, T_Outer_Shield3);
   OUS_33->SetLineColor(4);
   OUS_33->SetLineWidth(3);
   OUS_33->SetMarkerColor(4);
   OUS_33->SetMarkerSize(3);
   OUS_33->SetMarkerStyle(21);
   OUS_33->Draw("p");
   
   TGraph *OUS_44 = new TGraph(450, Xaxis_4, T_Outer_Shield4);
   OUS_44->SetLineColor(4);
   OUS_44->SetLineWidth(3);
   OUS_44->SetMarkerColor(4);
   OUS_44->SetMarkerSize(3);
   OUS_44->SetMarkerStyle(21);
   OUS_44->Draw("p");

   TGraph *OUS_55 = new TGraph(450, Xaxis_5, T_Outer_Shield5);
   OUS_55->SetLineColor(4);
   OUS_55->SetLineWidth(3);
   OUS_55->SetMarkerColor(4);
   OUS_55->SetMarkerSize(3);
   OUS_55->SetMarkerStyle(21);
   OUS_55->Draw("p");
   
   TGraph *OUS_66 = new TGraph(450, Xaxis_6, T_Outer_Shield6);
   OUS_66->SetLineColor(4);
   OUS_66->SetLineWidth(3);
   OUS_66->SetMarkerColor(4);
   OUS_66->SetMarkerSize(3);
   OUS_66->SetMarkerStyle(21);
   OUS_66->Draw("p");
   
   TGraph *OUS_77 = new TGraph(450, Xaxis_7, T_Outer_Shield7);
   OUS_77->SetLineColor(4);
   OUS_77->SetLineWidth(3);
   OUS_77->SetMarkerColor(4);
   OUS_77->SetMarkerSize(3);
   OUS_77->SetMarkerStyle(21);
   OUS_77->Draw("p");
   
   TGraph *OUS_88 = new TGraph(450, Xaxis_8, T_Outer_Shield8);
   OUS_88->SetLineColor(4);
   OUS_88->SetLineWidth(3);
   OUS_88->SetMarkerColor(4);
   OUS_88->SetMarkerSize(3);
   OUS_88->SetMarkerStyle(21);
   OUS_88->Draw("p");

   TGraph *OUS_99 = new TGraph(450, Xaxis_9, T_Outer_Shield9);
   OUS_99->SetLineColor(4);
   OUS_99->SetLineWidth(3);
   OUS_99->SetMarkerColor(4);
   OUS_99->SetMarkerSize(3);
   OUS_99->SetMarkerStyle(21);
   OUS_99->Draw("p");
   
   TGraph *OUS_1010 = new TGraph(450, Xaxis_10, T_Outer_Shield10);
   OUS_1010->SetLineColor(4);
   OUS_1010->SetLineWidth(3);
   OUS_1010->SetMarkerColor(4);
   OUS_1010->SetMarkerSize(3);
   OUS_1010->SetMarkerStyle(21);
   OUS_1010->Draw("p");

   TGraph *OUS_1111 = new TGraph(450, Xaxis_11, T_Outer_Shield11);
   OUS_1111->SetLineColor(4);
   OUS_1111->SetLineWidth(3);
   OUS_1111->SetMarkerColor(4);
   OUS_1111->SetMarkerSize(3);
   OUS_1111->SetMarkerStyle(21);
   OUS_1111->Draw("p");
   
   TGraph *OUS_1212 = new TGraph(450, Xaxis_12, T_Outer_Shield12);
   OUS_1212->SetLineColor(4);
   OUS_1212->SetLineWidth(3);
   OUS_1212->SetMarkerColor(4);
   OUS_1212->SetMarkerSize(3);
   OUS_1212->SetMarkerStyle(21);
   OUS_1212->Draw("p");
   
   TGraph *OUS_1313 = new TGraph(450, Xaxis_13, T_Outer_Shield13);
   OUS_1313->SetLineColor(4);
   OUS_1313->SetLineWidth(3);
   OUS_1313->SetMarkerColor(4);
   OUS_1313->SetMarkerSize(3);
   OUS_1313->SetMarkerStyle(21);
   OUS_1313->Draw("p");
   
   TGraph *OUS_1414 = new TGraph(450, Xaxis_14, T_Outer_Shield14);
   OUS_1414->SetLineColor(4);
   OUS_1414->SetLineWidth(3);
   OUS_1414->SetMarkerColor(4);
   OUS_1414->SetMarkerSize(3);
   OUS_1414->SetMarkerStyle(21);
   OUS_1414->Draw("p");

   TGraph *OUS_1515 = new TGraph(450, Xaxis_15, T_Outer_Shield15);
   OUS_1515->SetLineColor(4);
   OUS_1515->SetLineWidth(3);
   OUS_1515->SetMarkerColor(4);
   OUS_1515->SetMarkerSize(3);
   OUS_1515->SetMarkerStyle(21);
   OUS_1515->Draw("p");

   TGraph *OUS_1616 = new TGraph(450, Xaxis_16, T_Outer_Shield16);
   OUS_1616->SetLineColor(4);
   OUS_1616->SetLineWidth(3);
   OUS_1616->SetMarkerColor(4);
   OUS_1616->SetMarkerSize(3);
   OUS_1616->SetMarkerStyle(21);
   OUS_1616->Draw("p");

   TGraph *OUS_1717 = new TGraph(450, Xaxis_17, T_Outer_Shield17);
   OUS_1717->SetLineColor(4);
   OUS_1717->SetLineWidth(3);
   OUS_1717->SetMarkerColor(4);
   OUS_1717->SetMarkerSize(3);
   OUS_1717->SetMarkerStyle(21);
   OUS_1717->Draw("p");

   TGraph *OUS_1818 = new TGraph(450, Xaxis_18, T_Outer_Shield18);
   OUS_1818->SetLineColor(4);
   OUS_1818->SetLineWidth(3);
   OUS_1818->SetMarkerColor(4);
   OUS_1818->SetMarkerSize(3);
   OUS_1818->SetMarkerStyle(21);
   OUS_1818->Draw("p");

   TGraph *OUS_1919 = new TGraph(450, Xaxis_19, T_Outer_Shield19);
   OUS_1919->SetLineColor(4);
   OUS_1919->SetLineWidth(3);
   OUS_1919->SetMarkerColor(4);
   OUS_1919->SetMarkerSize(3);
   OUS_1919->SetMarkerStyle(21);
   OUS_1919->Draw("p");

   TGraph *OUS_2020 = new TGraph(450, Xaxis_20, T_Outer_Shield20);
   OUS_2020->SetLineColor(4);
   OUS_2020->SetLineWidth(3);
   OUS_2020->SetMarkerColor(4);
   OUS_2020->SetMarkerSize(3);
   OUS_2020->SetMarkerStyle(21);
   OUS_2020->Draw("p");
   
   TGraph *OUS_2021 = new TGraph(450, Xaxis_21, T_Outer_Shield21);
   OUS_2021->SetLineColor(4);
   OUS_2021->SetLineWidth(3);
   OUS_2021->SetMarkerColor(4);
   OUS_2021->SetMarkerSize(3);
   OUS_2021->SetMarkerStyle(21);
   OUS_2021->Draw("p");

   TGraph *OUS_2022 = new TGraph(450, Xaxis_22, T_Outer_Shield22);
   OUS_2022->SetLineColor(4);
   OUS_2022->SetLineWidth(3);
   OUS_2022->SetMarkerColor(4);
   OUS_2022->SetMarkerSize(3);
   OUS_2022->SetMarkerStyle(21);
   OUS_2022->Draw("p");
   
   TGraph *OUS_2023 = new TGraph(450, Xaxis_23, T_Outer_Shield23);
   OUS_2023->SetLineColor(4);
   OUS_2023->SetLineWidth(3);
   OUS_2023->SetMarkerColor(4);
   OUS_2023->SetMarkerSize(3);
   OUS_2023->SetMarkerStyle(21);
   OUS_2023->Draw("p");
   
   TGraph *OUS_2024 = new TGraph(450, Xaxis_24, T_Outer_Shield24);
   OUS_2024->SetLineColor(4);
   OUS_2024->SetLineWidth(3);
   OUS_2024->SetMarkerColor(4);
   OUS_2024->SetMarkerSize(3);
   OUS_2024->SetMarkerStyle(21);
   OUS_2024->Draw("p");
   
   //===============================Names ===================

   TLatex *tex = new TLatex(0.48,0.34,"Inside the Shielding (Near AC)");
   tex->SetNDC();
   tex->SetTextColor(1);
   tex->SetTextFont(22);
   tex->SetTextSize(0.050);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(0.48,0.28,"Control room");
   tex->SetNDC();
   tex->SetTextColor(4);
   tex->SetTextFont(22);
   tex->SetTextSize(0.050);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(0.48,0.22,"Near Outlet of Compressor");
   tex->SetNDC();
   tex->SetTextColor(kGreen+2);
   tex->SetTextFont(22);
   tex->SetTextSize(0.050);
   tex->SetLineWidth(2);
   tex->Draw();

   
   
}

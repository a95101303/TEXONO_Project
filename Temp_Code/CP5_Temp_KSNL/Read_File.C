#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TDatime.h"

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <string>
#include <filesystem>
#include <ctype.h>
namespace fs = std::filesystem;

double Number_Trans(string tmp_string)
{
    vector<char> vec_str(tmp_string.begin(), tmp_string.end());
    int size=vec_str.size();
    string Number_char="";
    for(int Loop=0; Loop<size-2; Loop++){Number_char = Number_char + vec_str[Loop];}
    Number_char = Number_char + '.';
    for(int Loop=size-2; Loop<size; Loop++){Number_char = Number_char + vec_str[Loop];}

    return stod(Number_char);
}

/*
void Read_File_1()
{
    cout << Number_Trans("3854");
vector<TDatime*> Date_check;
TDatime *t1 = new TDatime(2002,1,12,12,24,59);
Date_check.push_back(t1);
cout << Date_check[0]->GetYear();
}
*/
void Read_File()
{
    
    string filename1("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Temp_Code/CP5_Temp_KSNL/Original_Files/KS_201117_50x50_CP5.txt");
    string filename2("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Temp_Code/CP5_Temp_KSNL/Original_Files/KS_200908_70x70_CP5.txt");
    //cout << "filename: " << filename << endl;
    vector<TGraph*> Cooler_Power_Watts_50_70;
    vector<TGraph*> D_Coldtip_Temp_50_70;
    vector<TGraph*> D_Coldtip_Setp_50_70;
    vector<TGraph*> D_Coldhead_Temp_50_70;
    vector<TGraph*> D_Compressor_Temp_50_70;
    vector<TGraph*> D_Controller_Temp_50_70;
    //vector<TGraph*> D_Coldtip_Temp_50_70;
    TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.04,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);

    
    for(int File=0; File<2; File++)
    {
            string filename;
            if(File==0)filename = filename1;
            if(File==1)filename = filename2;
            ifstream input_file(filename);

            char byte = 0;int  Mark_Set = 0;
            int Start_Counting_Date=0;int Start_Counting_Underscore=0;
            vector<char> Bytes;vector<char> Bytes_Digit_char;Bytes.push_back(0);Bytes_Digit_char.push_back(0);
            string Date="";int Set_Started=0;
            vector<char> Line_Set;
            vector<TDatime*> Date_check;
            vector<double>   D_Cooler_Power_Watts;
            vector<double>   D_Coldtip_Temp;vector<double>   D_Coldtip_Setp;vector<double>   D_Coldhead_Temp;vector<double>   D_Compressor_Temp;vector<double>  D_Controller_Temp;
            while (input_file.get(byte))
            {
                
                  if(byte!=' ')
                  {
                      cout << byte ;
                      Mark_Set = Mark_Set + 1;
                      Bytes.push_back(byte);
                      //cout << "======" << endl;
                      //cout << "byte: " << byte << endl;
                      //cout << "======" << endl;
                  }
                  if(byte!=' ' and Bytes.size()>3)
                  {
                      Date  = Date + (Bytes[Mark_Set-3]);
                      Date  = Date + (Bytes[Mark_Set-2]);
                      Date  = Date + (Bytes[Mark_Set-1]);
                      Date  = Date + (Bytes[Mark_Set]);
                      //cout << "Date" << Date << endl;
                      //Mark the Date and start a new set
                      if(Date=="Date")
                      {Start_Counting_Date = Start_Counting_Date + 1;}
                      //Use - as a sign to start the number
                      if(Bytes[Mark_Set]=='-' and Bytes[Mark_Set-1]!='-' and isdigit(Bytes[Mark_Set-1])==0 and Bytes[Mark_Set-1]!='o')
                      {Start_Counting_Underscore = Start_Counting_Underscore + 1;}
                      //cout << "Start_Counting_Date: " << Start_Counting_Date << "Start_Counting_Underscore: " << Start_Counting_Underscore << endl;
                      if(Start_Counting_Date==1 and Start_Counting_Underscore==1 and (isdigit(Bytes[Mark_Set])==1 and isdigit(Bytes[Mark_Set-1])==0) )
                      {
                          //cout <<  byte ;
                          Set_Started=1;
                          //cout << "===S===" << endl;

                          /*
                          cout << "======" << endl;
                          cout << "byte: " << byte << endl;;
                          cout << "isdigit(Bytes[Mark_Set]): " << isdigit(Bytes[Mark_Set]) << endl;
                          cout << "Bytes[Mark_Set-1]: " << Bytes[Mark_Set-1] << endl;
                          cout << "======" << endl;
                           */
                           
                      }
                      //cout << "Start_Counting_Date: " << Start_Counting_Date << endl;
                      //cout << "Start_Counting_Underscore: " << Start_Counting_Underscore << endl;
        if(Start_Counting_Date==1 and Start_Counting_Underscore==1 and Set_Started==1 and (( isdigit(Bytes[Mark_Set])==1 or Bytes[Mark_Set]=='-') or (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N')) )
                      {
                          cout  << "OK: " ;
                          Line_Set.push_back(byte);
                          string Date_individual;string Month_individual ;string Year_individual;
                          string Hour_individual;string Minute_individual;string Second_individual;
                          string Cooler_Power_Watts;
                          string Coldtip_Temp;string Coldtip_Setp;string Coldhead_Temp;string Compressor_Temp;string Controller_Temp;
                          //cout << byte;
                          if( (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N') )
                          {
                              //Define the Date/Month/Year Hour/Minute/Second
                              int Minus_Number=0;
                              for(int Index=0; Index<Line_Set.size(); Index++)
                              {
                                  //cout  << Line_Set[Index] ;
                                  if(Line_Set[Index]=='-')Minus_Number=Minus_Number+1;
                              }
                              Date_individual    = Date_individual   + Line_Set[0]+Line_Set[1];
                              Month_individual   = Month_individual  + Line_Set[2]+Line_Set[3];
                              Year_individual    = Year_individual   + Line_Set[4]+Line_Set[5]+Line_Set[6]+Line_Set[7];
                              Hour_individual    = Hour_individual   + Line_Set[8]+Line_Set[9];
                              Minute_individual  = Minute_individual + Line_Set[10]+Line_Set[11];
                              if(Minus_Number==2)Second_individual  = Second_individual + Line_Set[12]+Line_Set[13];
                              else{Second_individual  = Second_individual + '8' + '0';}//Weird point
                              cout << "Minus_Number: " << Minus_Number << endl;

                              cout << "Date_individual: "  << Date_individual.c_str() << endl;
                              cout << "Month_individual:" << Month_individual.c_str() << endl;
                              cout << "Year_individual: "  << Year_individual.c_str() << endl;
                              cout << "Hour_individual: "  << Hour_individual.c_str() << endl;
                              cout << "Minute_individual: "  << Minute_individual.c_str() << endl;
                              cout << "Second_individual: "  << Second_individual.c_str() << endl;
          TDatime *t1 = new TDatime(stoi(Year_individual),stoi(Month_individual),stoi(Date_individual),stoi(Hour_individual),stoi(Minute_individual),stoi(Second_individual));
                              cout << "GetDate: " << t1->GetDate() << endl;
                              cout << "GetTime: " << t1->GetTime() << endl;
                              cout << "stoi(Second_individual): " << stoi(Second_individual) << endl;
                              
                              //(Cooler_Power, Temperature
                              if(stoi(Second_individual)!=80 and Minus_Number==2)
                              {
                                  Date_check.push_back(t1);
                                  int First_Minus_Index=0;int Second_Minus_Index=0;
                                  int First_Zero_Index=0;int Second_Zero_Index=0;int Third_Zero_Index=0;int Fourth_Zero_Index=0;
                                  for(int Index=14; Index<Line_Set.size(); Index++)
                                  {
                                      
                                      cout  << Line_Set[Index] ;
                                      if(First_Minus_Index ==0 and Line_Set[Index]=='-'){First_Minus_Index =Index;continue;}
                                      if(First_Minus_Index >0 and Second_Minus_Index==0 and Line_Set[Index]=='-'){Second_Minus_Index=Index;continue;}
            if(Second_Minus_Index>0 and First_Zero_Index==0 and Line_Set[Index-1]=='0' and Line_Set[Index]=='0' and Line_Set[Index+1]!='0'){First_Zero_Index =Index;continue;}
                                      Second_Zero_Index=First_Zero_Index+4;
                                      Third_Zero_Index=Second_Zero_Index+4;
                                      Fourth_Zero_Index=Third_Zero_Index+4;
            //if(First_Zero_Index>0  and Second_Zero_Index==0 and Line_Set[Index]=='0' and Line_Set[Index+1]!='0' and Line_Set[Index+1]!='.'){Second_Zero_Index=Index;continue;}
            //if(Second_Zero_Index>0 and Third_Zero_Index==0 and Line_Set[Index]=='0' and Line_Set[Index+1]!='0' and Line_Set[Index+1]!='.'){Third_Zero_Index =Index;continue;}
            //if(Third_Zero_Index>0  and Fourth_Zero_Index==0 and Line_Set[Index]=='0' and Line_Set[Index+1]=='N' and Line_Set[Index+1]!='.'){Fourth_Zero_Index=Index;continue;}
                                  }
                                  cout << "First_Minus_Index: " << Line_Set[First_Minus_Index+1] << endl;
                                  cout << "Second_Minus_Index: " << Line_Set[Second_Minus_Index+1] << endl;
                                  cout << "First_Zero_Index: " << Line_Set[First_Zero_Index+1] << endl;
                                  cout << "Second_Zero_Index: " << Line_Set[Second_Zero_Index+1] << endl;
                                  cout << "Third_Zero_Index: " << Line_Set[Third_Zero_Index+1] << endl;

                                  for(int i=14;i<First_Minus_Index;i++){Cooler_Power_Watts = Cooler_Power_Watts + Line_Set[i];}
                                  for(int i=First_Minus_Index+1;i<Second_Minus_Index;i++){cout << Line_Set[i];Coldtip_Temp = Coldtip_Temp + Line_Set[i];}
                                  for(int i=Second_Minus_Index+1;i<First_Zero_Index+1;i++){cout << Line_Set[i];Coldtip_Setp = Coldtip_Setp + Line_Set[i];}
                                  for(int i=First_Zero_Index+1;i<Second_Zero_Index+1;i++){cout << Line_Set[i];Coldhead_Temp = Coldhead_Temp + Line_Set[i];}
                                  for(int i=Second_Zero_Index+1;i<Third_Zero_Index+1;i++){cout << Line_Set[i];Compressor_Temp = Compressor_Temp + Line_Set[i];}
                                  for(int i=Third_Zero_Index+1;i<Fourth_Zero_Index+1;i++){cout << Line_Set[i];Controller_Temp = Controller_Temp + Line_Set[i];}

                                  cout << "Cooler_Power_Watts: " << Number_Trans(Cooler_Power_Watts) << endl;
                                  cout << "Coldtip_Temp: " << Number_Trans(Coldtip_Temp) << endl;
                                  cout << "Coldtip_Setp: " << Number_Trans(Coldtip_Setp) << endl;
                                  cout << "Coldhead_Temp: " << Number_Trans(Coldhead_Temp) << endl;
                                  cout << "Compressor_Temp: " << Number_Trans(Compressor_Temp) << endl;
                                  cout << "Controller_Temp: " << Number_Trans(Controller_Temp) << endl;

                                  D_Cooler_Power_Watts.push_back(Number_Trans(Cooler_Power_Watts));
                                  D_Coldtip_Temp.push_back(Number_Trans(Coldtip_Temp));
                                  D_Coldtip_Setp.push_back(Number_Trans(Coldtip_Setp));
                                  D_Coldhead_Temp.push_back(Number_Trans(Coldhead_Temp));
                                  D_Compressor_Temp.push_back(Number_Trans(Compressor_Temp));
                                  D_Controller_Temp.push_back(Number_Trans(Controller_Temp));
                                  
                                  //cout << "D_Controller_Temp[0]: " << D_Controller_Temp[0] << endl;

                              }//if(stoi(Second_individual)!=80 and Minus_Number==2)
                              Line_Set.clear();
                              Date_individual="";Month_individual="" ;Year_individual="";
                              Hour_individual="";Minute_individual="";Second_individual="";
                              Cooler_Power_Watts="";Coldtip_Temp="";Coldtip_Setp="";
                              Coldhead_Temp="";Compressor_Temp="";Controller_Temp="";
                          }//if( (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N') )
                      }//if(Start_Counting_Date==1 and Start_Counting_Underscore==1 and Set_Started==1 and ( isdigit(Bytes[Mark_Set])==1 or Bytes[Mark_Set]=='-') or (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N') )
                      Date="";
                  }//if(byte!=' ' and Bytes.size()>3)
                if( (Start_Counting_Date>0 and (Start_Counting_Date==Start_Counting_Underscore) and byte=='*') or byte=='*')
                {
                    cout << "End!" << endl;
                    Start_Counting_Date = 0;Start_Counting_Underscore=0;
                    Set_Started=0;
                }
                
            }//while (input_file.get(byte))
        
            vector<double> Time_Bef_Aft;
            cout << "===================================" << endl;
            for(int Index=0; Index<Date_check.size()-1; Index++)
            {
                /*
                cout << "Index: " << Index << endl;
                cout << "GetDate: " << Date_check[Index]->GetDate() << endl;
                cout << "GetTime: " << Date_check[Index]->GetTime() << endl;
                cout << "D_Cooler_Power_Watts: " << D_Cooler_Power_Watts[Index] << endl;
                cout << "D_Coldtip_Temp: " << D_Coldtip_Temp[Index] << endl;
                cout << "D_Coldtip_Setp: " << D_Coldtip_Setp[Index] << endl;
                cout << "D_Coldhead_Temp: " << D_Coldhead_Temp[Index] << endl;
                cout << "D_Compressor_Temp: " << D_Compressor_Temp[Index] << endl;
                cout << "Controller_Temp: " << D_Controller_Temp[Index] << endl;
                 */
                //cout << "GetDate: " << Date_check[Index]->GetDate() << endl;
                if(D_Cooler_Power_Watts[Index]>300.)
                {   cout << "GetDate: " << Date_check[Index]->GetDate() << endl;
                    cout << "GetTime: " << Date_check[Index]->GetTime() << endl;
                    cout << "D_Cooler_Power_Watts: " << D_Cooler_Power_Watts[Index] << endl;}
                //cout << "D_Coldtip_Setp: " << D_Coldtip_Setp[Index] << endl;
                /*
                if(D_Cooler_Power_Watts[Index]>200. or D_Cooler_Power_Watts[Index]<100.) cout << "D_Cooler_Power_Watts[Index]: " << D_Cooler_Power_Watts[Index] << endl;
                if(D_Coldtip_Temp[Index]<180. or D_Coldtip_Temp[Index]>210.) cout << "D_Coldtip_Temp[Index]: " << D_Coldtip_Temp[Index] << endl;
                if(D_Coldtip_Setp[Index]<180. or D_Coldtip_Setp[Index]>210.) cout << "D_Coldtip_Setp[Index]: " << D_Coldtip_Setp[Index] << endl;
                if(D_Coldhead_Temp[Index]>50. or D_Coldhead_Temp[Index]<0.) cout << "D_Coldhead_Temp[Index]: " << D_Coldhead_Temp[Index] << endl;
                if(D_Compressor_Temp[Index]>50. or D_Compressor_Temp[Index]<0.) cout << "D_Compressor_Temp[Index]: " << D_Compressor_Temp[Index] << endl;
                if(D_Controller_Temp[Index]>100.)cout << "Controller_Temp: " << D_Controller_Temp[Index] << endl;
                 */
                D_Coldtip_Temp[Index] = -D_Coldtip_Temp[Index];
                D_Coldtip_Setp[Index] = -D_Coldtip_Setp[Index];
                Time_Bef_Aft.push_back(Index);
            }
            cout << "===================================" << endl;


        TGraph *TGaph_D_Cooler_Power_Watts = new TGraph(Date_check.size()-1,&Time_Bef_Aft[0],&D_Cooler_Power_Watts[0]);//Total Cross Section(TCS), Electron Number(EN)
        Cooler_Power_Watts_50_70.push_back(TGaph_D_Cooler_Power_Watts);
        //TGaph_D_Cooler_Power_Watts->Draw("AP");
        TGraph *TGaph_D_Coldtip_Temp = new TGraph(Date_check.size()-1,&Time_Bef_Aft[0],&D_Coldtip_Temp[0]);//Total Cross Section(TCS), Electron Number(EN)
        D_Coldtip_Temp_50_70.push_back(TGaph_D_Coldtip_Temp);
        //TGaph_D_Coldtip_Temp->Draw("AP");
        TGraph *TGaph_D_Coldtip_Setp = new TGraph(Date_check.size()-1,&Time_Bef_Aft[0],&D_Coldtip_Setp[0]);//Total Cross Section(TCS), Electron Number(EN)
        D_Coldtip_Setp_50_70.push_back(TGaph_D_Coldtip_Setp);
        //TGaph_D_Coldtip_Setp->Draw("AP");
        TGraph *TGaph_D_Coldhead_Temp = new TGraph(Date_check.size()-1,&Time_Bef_Aft[0],&D_Coldhead_Temp[0]);//Total Cross Section(TCS), Electron Number(EN)
        D_Coldhead_Temp_50_70.push_back(TGaph_D_Coldhead_Temp);
        //TGaph_D_Coldhead_Temp->Draw("AP");
        TGraph *TGaph_D_Compressor_Temp = new TGraph(Date_check.size()-1,&Time_Bef_Aft[0],&D_Compressor_Temp[0]);//Total Cross Section(TCS), Electron Number(EN)
        D_Compressor_Temp_50_70.push_back(TGaph_D_Compressor_Temp);
        //TGaph_D_Compressor_Temp->Draw("AP");
        TGraph *TGaph_D_Controller_Temp = new TGraph(Date_check.size()-1,&Time_Bef_Aft[0],&D_Controller_Temp[0]);//Total Cross Section(TCS), Electron Number(EN)
        D_Controller_Temp_50_70.push_back(TGaph_D_Controller_Temp);
        //TGaph_D_Controller_Temp->Draw("AP");
        Bytes.clear();Bytes_Digit_char.clear();Date_check.clear();
        D_Cooler_Power_Watts.clear();D_Coldtip_Temp.clear();D_Coldtip_Setp.clear();D_Coldhead_Temp.clear();D_Compressor_Temp.clear();D_Controller_Temp.clear();

    }
        char fout_name[100];
        sprintf(fout_name,"/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Temp_Code/Temp.root");
        TFile *fout=new TFile(fout_name,"recreate");
        Cooler_Power_Watts_50_70[0]->Write("Cooler_Power_Watts_50");
        Cooler_Power_Watts_50_70[1]->Write("Cooler_Power_Watts_70");
        D_Coldtip_Temp_50_70[0]->Write("Coldtip_Temp_50");
        D_Coldtip_Temp_50_70[1]->Write("Coldtip_Temp_70");
        D_Coldtip_Setp_50_70[0]->Write("Coldtip_Setp_50");
        D_Coldtip_Setp_50_70[1]->Write("Coldtip_Setp_70");
        D_Coldhead_Temp_50_70[0]->Write("Coldhead_Temp_50");
        D_Coldhead_Temp_50_70[1]->Write("Coldhead_Temp_70");
        D_Compressor_Temp_50_70[0]->Write("Compressor_Temp_50");
        D_Compressor_Temp_50_70[1]->Write("Compressor_Temp_70");
        D_Controller_Temp_50_70[0]->Write("Controller_Temp_50");
        D_Controller_Temp_50_70[1]->Write("Controller_Temp_70");

        /*
        TGraph *TGaph_D_Cooler_Power_Watts = new TGraph(Date_check.size()-1,&Time_Bef_Aft[0],&D_Cooler_Power_Watts[0]);//Total Cross Section(TCS), Electron Number(EN)
        TGaph_D_Cooler_Power_Watts->SetMarkerStyle(20);
        TGaph_D_Cooler_Power_Watts->SetMarkerColor(3);
        TGaph_D_Cooler_Power_Watts->SetLineColor(3);
        TGaph_D_Cooler_Power_Watts->SetLineWidth(5);
        TGaph_D_Cooler_Power_Watts->GetXaxis()->SetTitle("Date");
        TGaph_D_Cooler_Power_Watts->GetYaxis()->SetTitle("Cooler Power[W]");
        TGaph_D_Cooler_Power_Watts->GetXaxis()->SetRangeUser(0.,Time_Bef_Aft.size()+1);
        TGaph_D_Cooler_Power_Watts->GetYaxis()->SetRangeUser(0.,200.);
        TGaph_D_Cooler_Power_Watts->Draw("apl");
     */
}
/*
 
 TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
 gStyle->SetOptStat(0);
 gStyle->SetTitleSize(0.04,"XY");
 gStyle->SetTitleFont(62,"XY");
 gStyle->SetLegendFont(62);
 
 c1->SetLogy();
 c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Temp_Code/Cooler.pdf");

 string AAA="2300";
 cout << "AAA: " << AAA << endl;
 int Number = stoi(AAA);
 cout << "Number: " << Number << endl;
*/

/*
     string filename("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Temp_Code/CP5_Temp_KSNL/Original_Files/KS_200908_70x70_CP5.txt");
    //cout << "filename: " << filename << endl;
    ifstream input_file(filename);//Read the texts in the file

    char byte = 0;int  Mark_Set = 0;
    int Start_Counting_Date=0;int Start_Counting_Underscore=0;
    vector<char> Bytes;vector<char> Bytes_Digit_char;Bytes.push_back(0);Bytes_Digit_char.push_back(0);
    string Date="";int Set_Started=0;
    vector<char> Line_Set;
    vector<TDatime*> Date_check;
    vector<double>   D_Cooler_Power_Watts;
    vector<double>   D_Coldtip_Temp;vector<double>   D_Coldtip_Setp;vector<double>   D_Coldhead_Temp;vector<double>   D_Compressor_Temp;vector<double>  D_Controller_Temp;
    while (input_file.get(byte))
    {
        
          if(byte!=' ')
          {
              Mark_Set = Mark_Set + 1;
              Bytes.push_back(byte);
              //cout << "======" << endl;
              //cout << "byte: " << byte << endl;
              //cout << "======" << endl;
          }
          if(byte!=' ' and Bytes.size()>3)
          {
              Date  = Date + (Bytes[Mark_Set-3]);
              Date  = Date + (Bytes[Mark_Set-2]);
              Date  = Date + (Bytes[Mark_Set-1]);
              Date  = Date + (Bytes[Mark_Set]);
              //cout << "Date" << Date << endl;
              //Mark the Date and start a new set
              if(Date=="Date")
              {Start_Counting_Date = Start_Counting_Date + 1;}
              //Use - as a sign to start the number
              if(Bytes[Mark_Set]=='-' and Bytes[Mark_Set-1]!='-' and isdigit(Bytes[Mark_Set-1])==0 and Bytes[Mark_Set-1]!='o')
              {Start_Counting_Underscore = Start_Counting_Underscore + 1;}
              //cout << "Start_Counting_Date: " << Start_Counting_Date << "Start_Counting_Underscore: " << Start_Counting_Underscore << endl;
              if(Start_Counting_Date==1 and Start_Counting_Underscore==1 and (isdigit(Bytes[Mark_Set])==1 and isdigit(Bytes[Mark_Set-1])==0) )
              {
                  //cout <<  byte ;
                  Set_Started=1;
                  //cout << "===S===" << endl;

                  
                  //cout << "======" << endl;
                  //cout << "byte: " << byte << endl;;
                  //cout << "isdigit(Bytes[Mark_Set]): " << isdigit(Bytes[Mark_Set]) << endl;
                  //cout << "Bytes[Mark_Set-1]: " << Bytes[Mark_Set-1] << endl;
                  //cout << "======" << endl;
                   
                   
              }
              //cout << "Start_Counting_Date: " << Start_Counting_Date << endl;
              //cout << "Start_Counting_Underscore: " << Start_Counting_Underscore << endl;
if(Start_Counting_Date==1 and Start_Counting_Underscore==1 and Set_Started==1 and (( isdigit(Bytes[Mark_Set])==1 or Bytes[Mark_Set]=='-') or (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N')) )
              {
                  //cout  << "OK: " << byte ;
                  Line_Set.push_back(byte);
                  string Date_individual;string Month_individual ;string Year_individual;
                  string Hour_individual;string Minute_individual;string Second_individual;
                  string Cooler_Power_Watts;
                  string Coldtip_Temp;string Coldtip_Setp;string Coldhead_Temp;string Compressor_Temp;string Controller_Temp;
                  //cout << byte;
                  if( (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N') )
                  {
                      //Define the Date/Month/Year Hour/Minute/Second
                      int Minus_Number=0;
                      for(int Index=0; Index<Line_Set.size(); Index++)
                      {
                          //cout  << Line_Set[Index] ;
                          if(Line_Set[Index]=='-')Minus_Number=Minus_Number+1;
                      }
                      //cout << "Minus_Number: " << Minus_Number << endl;
                      
                      Date_individual    = Date_individual   + Line_Set[0]+Line_Set[1];
                      Month_individual   = Month_individual  + Line_Set[2]+Line_Set[3];
                      Year_individual    = Year_individual   + Line_Set[4]+Line_Set[5]+Line_Set[6]+Line_Set[7];
                      Hour_individual    = Hour_individual   + Line_Set[8]+Line_Set[9];
                      Minute_individual  = Minute_individual + Line_Set[10]+Line_Set[11];
                      if(Minus_Number==2)Second_individual  = Second_individual + Line_Set[12]+Line_Set[13];
                      else{Second_individual  = Second_individual + '8' + '0';}//Weird point
                      cout << "Date_individual: "  << Date_individual.c_str() << endl;
                      cout << "Month_individual:" << Month_individual.c_str() << endl;
                      cout << "Year_individual: "  << Year_individual.c_str() << endl;
                      cout << "Hour_individual: "  << Hour_individual.c_str() << endl;
                      cout << "Minute_individual: "  << Minute_individual.c_str() << endl;
                      cout << "Second_individual: "  << Second_individual.c_str() << endl;
  TDatime *t1 = new TDatime(stoi(Year_individual),stoi(Month_individual),stoi(Date_individual),stoi(Hour_individual),stoi(Minute_individual),stoi(Second_individual));
                      cout << "GetDate: " << t1->GetDate() << endl;
                      cout << "GetTime: " << t1->GetTime() << endl;
                      cout << "stoi(Second_individual): " << stoi(Second_individual) << endl;
                      Date_check.push_back(t1);
                      
                      //(Cooler_Power, Temperature
                      if(stoi(Second_individual)!=80 and Minus_Number==2)
                      {
                          int First_Minus_Index=0;int Second_Minus_Index=0;
                          int First_Zero_Index=0;int Second_Zero_Index=0;int Third_Zero_Index=0;int Fourth_Zero_Index=0;
                          for(int Index=14; Index<Line_Set.size(); Index++)
                          {
                              cout  << Line_Set[Index] ;
                              if(First_Minus_Index ==0 and Line_Set[Index]=='-'){First_Minus_Index =Index;continue;}
                              if(First_Minus_Index >0 and Second_Minus_Index==0 and Line_Set[Index]=='-'){Second_Minus_Index=Index;continue;}
                              if(Second_Minus_Index>0 and First_Zero_Index==0 and Line_Set[Index]=='0' and Line_Set[Index+1]!='0'){First_Zero_Index =Index;continue;}
                              if(First_Zero_Index>0   and Second_Zero_Index==0 and Line_Set[Index]=='0' and Line_Set[Index+1]!='0'){Second_Zero_Index=Index;continue;}
                              if(Second_Zero_Index>0  and Third_Zero_Index==0 and Line_Set[Index]=='0' and Line_Set[Index+1]!='0'){Third_Zero_Index =Index;continue;}
                              if(Third_Zero_Index>0   and Fourth_Zero_Index==0 and Line_Set[Index]=='0' and Line_Set[Index+1]=='N'){Fourth_Zero_Index=Index;continue;}
                          }
                          //cout << "First_Minus_Index: " << First_Minus_Index << endl;
                          //cout << "Second_Minus_Index: " << Second_Minus_Index << endl;
                          //cout << "First_Zero_Index: " << First_Zero_Index << endl;
                          //cout << "Second_Zero_Index: " << Second_Zero_Index << endl;
                          //cout << "Third_Zero_Index: " << Third_Zero_Index << endl;

                          for(int i=First_Minus_Index-5;i<First_Minus_Index;i++){Cooler_Power_Watts = Cooler_Power_Watts + Line_Set[i];}
                          for(int i=First_Minus_Index+1;i<Second_Minus_Index;i++){cout << Line_Set[i];Coldtip_Temp = Coldtip_Temp + Line_Set[i];}
                          for(int i=Second_Minus_Index+1;i<First_Zero_Index+1;i++){cout << Line_Set[i];Coldtip_Setp = Coldtip_Setp + Line_Set[i];}
                          for(int i=First_Zero_Index+1;i<Second_Zero_Index+1;i++){cout << Line_Set[i];Coldhead_Temp = Coldhead_Temp + Line_Set[i];}
                          for(int i=Second_Zero_Index+1;i<Third_Zero_Index+1;i++){cout << Line_Set[i];Compressor_Temp = Compressor_Temp + Line_Set[i];}
                          for(int i=Third_Zero_Index+1;i<Fourth_Zero_Index+1;i++){cout << Line_Set[i];Controller_Temp = Controller_Temp + Line_Set[i];}
                          cout << "Cooler_Power_Watts: " << Number_Trans(Cooler_Power_Watts) << endl;
                          cout << "Coldtip_Temp: " << Number_Trans(Coldtip_Temp) << endl;
                          cout << "Coldtip_Setp: " << Number_Trans(Coldtip_Setp) << endl;
                          cout << "Coldhead_Temp: " << Number_Trans(Coldhead_Temp) << endl;
                          cout << "Compressor_Temp: " << Number_Trans(Compressor_Temp) << endl;
                          cout << "Controller_Temp: " << Number_Trans(Controller_Temp) << endl;

                          D_Cooler_Power_Watts.push_back(Number_Trans(Cooler_Power_Watts));
                          D_Coldtip_Temp.push_back(Number_Trans(Coldtip_Temp));
                          D_Coldtip_Setp.push_back(Number_Trans(Coldtip_Setp));
                          D_Coldhead_Temp.push_back(Number_Trans(Coldhead_Temp));
                          D_Compressor_Temp.push_back(Number_Trans(Compressor_Temp));
                          D_Controller_Temp.push_back(Number_Trans(Controller_Temp));
                          
                          //cout << "D_Controller_Temp[0]: " << D_Controller_Temp[0] << endl;

                      }//if(stoi(Second_individual)!=80 and Minus_Number==2)
                      Line_Set.clear();
                      Date_individual="";Month_individual="" ;Year_individual="";
                      Hour_individual="";Minute_individual="";Second_individual="";
                      Cooler_Power_Watts="";Coldtip_Temp="";Coldtip_Setp="";
                      Coldhead_Temp="";Compressor_Temp="";Controller_Temp="";
                  }//if( (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N') )
              }//if(Start_Counting_Date==1 and Start_Counting_Underscore==1 and Set_Started==1 and ( isdigit(Bytes[Mark_Set])==1 or Bytes[Mark_Set]=='-') or (Bytes[Mark_Set-1]=='0' and  Bytes[Mark_Set]=='N') )
              Date="";
          }//if(byte!=' ' and Bytes.size()>3)
        if( Start_Counting_Date>0 and (Start_Counting_Date==Start_Counting_Underscore) and byte=='*')
        {
            cout << "End!" << endl;
            Start_Counting_Date = 0;Start_Counting_Underscore=0;
            Set_Started=0;
        }
        
    }//while (input_file.get(byte))
*/

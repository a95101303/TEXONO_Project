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

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <string>
#include <filesystem>
namespace fs = std::filesystem;


void Read_File()
{
        string filename("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Temp_Code/CP5_Temp_KSNL/Original_Files/KS_201117_50x50_CP5.txt");
        cout << "filename: " << filename << endl;
        ifstream input_file(filename);//Read the texts in the file
    
        char byte = 0;int  Mark = 0;
        vector<char> bytes;bytes.push_back(0);
        string Year="";
        while (input_file.get(byte))
        {
              if(byte!=' ')cout << "byte: " << byte << endl;
              Mark = Mark + 1;
              bytes.push_back(byte);
              if(Mark>3) Year  = Year + (bytes[Mark-3] + bytes[Mark-2] + bytes[Mark-1] + bytes[Mark]);
              cout << "Year: " << Year.c_str() << endl;
              //if(stoi(Year)==2020 or stoi(Year)==2021)cout << "FUCK!!!" << endl;
        }
         
         

}

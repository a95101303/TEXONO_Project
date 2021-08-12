#include "DataJune_LongRange/DataJune_0_06GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_07GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_08GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_09GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_10GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_12GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_20GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_25GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_30GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_40GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_50GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_60GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_70GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_80GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_0_90GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_1_0GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_1_3GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_1_5GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_2_0GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_2_5GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_3_0GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_3_5GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_4_0GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_5_0GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_10_0GeV_d1_dcs.h"
#include "DataJune_LongRange/DataJune_20_0GeV_d1_dcs.h"


double LongRangeDcs_June(double T, double mass)
{
  double result;
  double a, b, c, d, y;
  double TKeV[DataBin], DCS[DataBin];

  //if(mass != 0.06 || mass != 0.07  ||  mass != 0.08  ||
  //   mass != 0.09 || mass != 0.10  ||  mass != 0.12  ||
  //   mass != 0.20 || mass != 0.25  ||  mass != 0.30  ||
  //   mass != 0.40 || mass != 0.50  ||  mass != 0.60  ||
  //   mass != 0.70 || mass != 0.80  ||  mass != 0.90  ||
  //   mass != 1.00 || mass != 1.30  ||  mass != 1.50  ||
  //   mass != 2.00 || mass != 2.50  ||  mass != 3.00  ||
  //   mass != 3.50 || mass != 4.00  ||  mass != 5.00  ||
  //   mass != 10.00 ||mass != 20.00 )
  //  {
  //    printf("This is Mass is not Calculated \n");
  //    return ;
  //  }
  //
  if(mass == 0.06)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_06GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_06GeV_d1_dcs[i][1];
	}
    }

  if(mass == 0.07)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_07GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_07GeV_d1_dcs[i][1];
	}
    }

  if(mass == 0.08)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_08GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_08GeV_d1_dcs[i][1];
	}
    }



  if(mass == 0.09)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_09GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_09GeV_d1_dcs[i][1];
	}
    }


if(mass == 0.10)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_10GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_10GeV_d1_dcs[i][1];
	}
    }



if(mass == 0.12)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_12GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_12GeV_d1_dcs[i][1];
	}
    }




if(mass == 0.20)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_20GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_20GeV_d1_dcs[i][1];
	}
    }



if(mass == 0.25)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_25GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_25GeV_d1_dcs[i][1];
	}
    }



if(mass == 0.30)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_30GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_30GeV_d1_dcs[i][1];
	}
    }




if(mass == 0.40)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_40GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_40GeV_d1_dcs[i][1];
	}
    }




 if(mass == 0.50)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_50GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_50GeV_d1_dcs[i][1];
	}
    }


if(mass == 0.60)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_60GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_60GeV_d1_dcs[i][1];
	}
    }



 if(mass == 0.70)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_70GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_70GeV_d1_dcs[i][1];
	}
    }


if(mass == 0.80)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_80GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_80GeV_d1_dcs[i][1];
	}
    }


if(mass == 0.90)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_0_90GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_0_90GeV_d1_dcs[i][1];
	}
    }



  if(mass == 1.00)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_1_0GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_1_0GeV_d1_dcs[i][1];
	}
    }


if(mass == 1.30)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_1_3GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_1_3GeV_d1_dcs[i][1];
	}
    }


if(mass == 1.50)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_1_5GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_1_5GeV_d1_dcs[i][1];
	}
    }


if(mass == 2.00)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_2_0GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_2_0GeV_d1_dcs[i][1];
	}
    }


if(mass == 2.50)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_2_5GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_2_5GeV_d1_dcs[i][1];
	}
    }



if(mass == 3.00)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_3_0GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_3_0GeV_d1_dcs[i][1];
	}
    }


if(mass == 3.50)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_3_5GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_3_5GeV_d1_dcs[i][1];
	}
    }


if(mass == 4.00)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_4_0GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_4_0GeV_d1_dcs[i][1];
	}
    }


if(mass == 5.00)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_5_0GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_5_0GeV_d1_dcs[i][1];
	}
    }


if(mass == 10.00)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_10_0GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_10_0GeV_d1_dcs[i][1];
	}
    }


if(mass == 20.00)
    {

      for(int i=0;i<DataBin;i++)
	{
	  TKeV[i] = DataJune_20_0GeV_d1_dcs[i][0];
	  DCS[i] =  DataJune_20_0GeV_d1_dcs[i][1];
	}
    }



  
  if(T<=TKeV[0])
    {
      result = DCS[0];
    }
  
  for(int i=0;i<(DataBin-1);i++)
    {
      if((T>TKeV[i])&&(T<=TKeV[i+1]))
	{
	  a = T-TKeV[i];
	  b = TKeV[i+1]-TKeV[i];
	  c = DCS[i+1] - DCS[i];
	  d = DCS[i];
	  result = ((a/b)*c) + d;
	   
	}
    }
  
   
  return result;


}

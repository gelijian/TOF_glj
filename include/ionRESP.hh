#include <fstream>
#include <cmath>
#include <string>
using namespace std;


double CECILresp(double E,string pname)
{
	double Eion=E/1000;     //in unit MeV
	double PHLight;
	
	if(pname == "proton")
	{
		PHLight = 0.95*Eion-8*(1-exp(-0.1*pow(Eion,0.9)));          // NIM 161(1979)439-447
		if(PHLight<0)
			PHLight = 0;
		return PHLight*1000;     // in keV
	}
	if(pname == "e-" || pname == "e+")
	{
		return Eion*1000;     // in keV
	}
	if(pname == "alpha")
	{
		PHLight = 0.41*Eion-5.9*(1.0-exp(-0.065*pow(Eion,1.01)));     // NIM 161(1979)439-447
		if(PHLight<0)
			PHLight = 0;
		return PHLight*1000;     // in keV
	}
	
	G4cout<<"No Light Output Data for Ions: "<<pname<<" exists!!!"<<G4endl;
	return 0;
}

double Songresp(double E,string pname)
{
	double Eion=E/1000;     //in unit MeV
	double PHLight;
	
	if(pname== "e-" || pname== "e+")
	{
		return Eion*1000;     // in MeV
	}

	if(pname == "proton")
	{
		if(Eion>=0&&Eion<=0.45)
			PHLight = 0.302*pow(Eion,5)+0.019*pow(Eion,4)-0.198*pow(Eion,3)+0.183*pow(Eion,2)+0.037*Eion+0.001;
		if(Eion>0.45&&Eion<0.55)
			PHLight = (0.302*pow(Eion,5)+0.019*pow(Eion,4)-0.198*pow(Eion,3)+0.183*pow(Eion,2)+0.037*Eion+0.001)*(0.55-Eion)/0.1+(-0.021*pow(Eion,5)+0.14*pow(Eion,4)-0.367*pow(Eion,3)+0.54*pow(Eion,2)-0.181*Eion+0.044)*(Eion-0.45)/0.1;
		if(Eion>=0.55&&Eion<=1.9)
			PHLight = -0.021*pow(Eion,5)+0.14*pow(Eion,4)-0.367*pow(Eion,3)+0.54*pow(Eion,2)-0.181*Eion+0.044;
		if(Eion>1.9&&Eion<2.1)
			PHLight = (-0.021*pow(Eion,5)+0.14*pow(Eion,4)-0.367*pow(Eion,3)+0.54*pow(Eion,2)-0.181*Eion+0.044)*(2.1-Eion)/0.2+(-1.4E-4*pow(Eion,5)+0.00452*pow(Eion,4)-0.05405*pow(Eion,3)+0.32042*pow(Eion,2)-0.46104*Eion+0.49415)*(Eion-1.9)/0.2;
		if(Eion>=2.1&&Eion<=4.9)
			PHLight = -1.4E-4*pow(Eion,5)+0.00452*pow(Eion,4)-0.05405*pow(Eion,3)+0.32042*pow(Eion,2)-0.46104*Eion+0.49415;
		if(Eion>4.9&&Eion<5.1)
			PHLight = (-1.4E-4*pow(Eion,5)+0.00452*pow(Eion,4)-0.05405*pow(Eion,3)+0.32042*pow(Eion,2)-0.46104*Eion+0.49415)*(5.1-Eion)/0.2+(2.9E-8*pow(Eion,6)-3.68E-6*pow(Eion,5)+0.0001808183*pow(Eion,4)-0.004413601*pow(Eion,3)+0.0601558495*pow(Eion,2)+0.1523790656*Eion-0.0306757377)*(Eion-4.9)/0.2;
		if(Eion>=5.1&&Eion<=40.0)
			PHLight = 2.9E-8*pow(Eion,6)-3.68E-6*pow(Eion,5)+0.0001808183*pow(Eion,4)-0.004413601*pow(Eion,3)+0.0601558495*pow(Eion,2)+0.1523790656*Eion-0.0306757377;
		if(Eion>40.0)
			PHLight = 2.31E-10*pow(Eion,6)-1.23896E-7*pow(Eion,5)+1.7645469E-5*pow(Eion,4)-0.001015910879*pow(Eion,3)+0.029063489921*Eion*Eion+0.250762837876*Eion-0.076982286193;
		if(PHLight<0)
			PHLight = 0;
		return PHLight*1000;     // in keV
	}

	if(pname == "alpha")
	{
		if(Eion>=0&&Eion<=0.5)
			PHLight = 0.1677922184*pow(Eion,6)-0.2200822987*pow(Eion,5)+0.0787625366*pow(Eion,4)+0.0042364096*pow(Eion,3)-0.000065703*Eion*Eion+0.0147137075*Eion+0.0001578641;
		if(Eion>0.5&&Eion<=4.0)
			PHLight = 0.0002546146*pow(Eion,6)-0.0034573049*pow(Eion,5)+0.0180969426*pow(Eion,4)-0.0455116438*pow(Eion,3)+0.0645135764*Eion*Eion-0.0201880666*Eion+0.0073383595;
		if(Eion>4.0)
			PHLight = -8.6E-9*pow(Eion,6)+1.1229E-6*pow(Eion,5)-5.01551E-5*pow(Eion,4)+6.158314E-4*pow(Eion,3)+0.0176338525*Eion*Eion-0.0889066524*Eion+0.2239099841;
		if(PHLight<0)
			PHLight = 0;
		return PHLight*1000;     // in keV
	}

	if(pname == "C12"||pname == "C13[0.0]" || pname == "C13")
	{
		PHLight = 3.3E-9*pow(Eion,5)-3.464E-7*pow(Eion,4)+1.32441E-5*pow(Eion,3)-1.343176E-4*pow(Eion,2)+5.6415155E-3*Eion+8.14417E-4;
		if(PHLight<0)
			PHLight = 0;
		return PHLight*1000;     // in keV
	}

	G4cout<<"No Light Output Data for Ions: "<<pname<<" exists!!!"<<G4endl;
	return 0;
}

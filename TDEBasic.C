/*
TDEBasic takes an array/waveform created by python and 
*/
#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<Python.h>
using namespace std;

int iLength(int[]);
int fLength(double[]);
void waveformMaker(double[],double[]);
double averager(double[],int);
void dtCoorelator(double[],double[],double[],int[],int,int);

int main()
{
  double waveformA[24]={};
  double waveformB[24]={};
  double Rsquared[6]={};
  int starts[6]={};
  int choice=0;
  int count,divisions,smallDT,dt,coorInitial,beginner;
  double aveA,aveB;
  cout<<fLength(waveformA)<<endl<<iLength(starts)<<endl;
  count=divisions=smallDT=coorInitial=beginner=0;
  dt=4;
  cout<<"Do you have python.h yet? 1 for yes, 0 for no:"<<endl;
  cin>>choice;
  if(choice==0)
    waveformMaker(waveformA,waveformB);
  if (fLength(waveformA)==fLength(waveformB))
    {
      aveA=averager(waveformA,fLength(waveformA));
      aveB=averager(waveformB,fLength(waveformB));
      for(int n=0;n<fLength(waveformA);n++)
	{
	  waveformA[n]=waveformA[n]-aveA;
	  waveformB[n]=waveformB[n]-aveB;
	}
      count=fLength(waveformA);
      divisions=count/dt;
      if(count%dt!=0)
        {
	  smallDT=count%dt;
          divisions=divisions+1;
          cout<<"\nLast cell has a dt of "<<smallDT<<endl;
        }
      for(int j=coorInitial;j<divisions;j++)
	{
	Rsquared[j]=0;
	starts[j]=beginner;
	beginner=beginner+dt;
	}
      starts[divisions-1]=count;
      dtCoorelator(waveformA,waveformB,Rsquared,starts,dt,smallDT);
      //cout<<Rsquared
    }
  
  return 0;
}

int iLength(int R[])
{
  return (sizeof(R)/sizeof(R[0]));
}

int fLength(double R[])
{
  return (sizeof(R)/sizeof(R[0]));
}

void waveformMaker(double A[24],double B[24])
{
  srand(time(0));
  double pi=atan(1)*4;
  int Period=rand()%24+1;
  double f=1.0/Period;
  double maxNoise,SNR;
  int Amp,shift;
  shift=rand()%24;
  Amp=(rand()%81+20);SNR=(rand()%9)*((rand()%100+1)/100.0)+2;
  maxNoise=Amp/SNR;
  
  for (int m=0;m<24;m++)
    {
      A[m]=Amp*sin(2*pi*f*m)+(maxNoise*(rand()%100)/100.0)-(maxNoise*(rand()%100/100.0));
      B[m]=Amp*sin(2*pi*f*m+shift)+(maxNoise*(rand()%100)/100.0)-(maxNoise*(rand()%100/100.0));
    }
  for (int p=0;p<24;p++)
    cout<<A[p]<<" ";
  cout<<endl;
  for (int q=0;q<24;q++)
    cout<<B[q]<<" ";
  cout<<endl<<"Amp="<<Amp<<" period="<<Period<<" shift="<<shift<<" SNR="<<SNR<<endl;
  return;
}

double averager(double X[],int size)
{
  double sum=0;
  for(int k=0;k<size;k++)
    sum=sum+X[k];
  return (sum/size);
}

void dtCoorelator(double A[],double B[],double R[],int starter[],int dt,int smallDT)
{
  int cell=0;
  int sum=0;
  int l=0;
  cout<<fLength(R);
  for(l=0;l<fLength(R)-1;l++)
    {
   for(int i=starter[l];i<=starter[l]+dt;i++)
    {
      sum=sum+(A[i+cell]*B[i+cell]);
      cell++;
    }
      R[l]=sum;
      cout<<R[l]<<endl;
      cell=sum=0;
    }
  for(int o=starter[l];o<=starter[l]+smallDT;o++)
    {
      sum=sum+(A[starter[l]+cell]*B[starter[l]+cell]);
      cell++;
    }
  return;
}

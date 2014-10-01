
#include <iostream>
#include <iomanip>
#include <bitset>
#include "quadlods.h"

using namespace std;

quadlods quads;

void testcoverage()
/* Cycle sizes: 408, 571, 377.
 * Total cycle size: 87828936 (LCM).
 * Remaining factor: 1.
 */
{
  int i,j,x,y,z,inx;
  quadlods cov;
  bitset<87828936> *histo;
  vector<mpq_class> point;
  cov.init(3,370);
  histo=new bitset<87828936>;
  for (i=0;i<87828936;i++)
  {
    point=cov.gen();
    x=mpz_class(point[0]*408).get_si();
    y=mpz_class(point[1]*571).get_si();
    z=mpz_class(point[2]*377).get_si();
    inx=x+408*(y+571*z);
    histo->set(inx);
    if (i%878290==0)
    {
      cout<<setw(2)<<i/878290<<"%\b\b\b";
      cout.flush();
    }
  }
  cout<<histo->count()<<endl;
  if (histo->count()==87828936)
    cout<<"Complete coverage"<<endl;
  delete histo;
}

int main(int argc,char **argv)
{
  int i,j;
  vector<double> point;
  quads.init(5,1e10);
  quads.advance(-1);
  for (i=0;i<5;i++)
    cout<<quads.num[i]<<'/'<<quads.denom[i]<<' '<<quads.acc[i]<<endl;
  for (i=0;i<30;i++)
  {
    point=quads.dgen();
    for (j=0;j<5;j++)
      cout<<point[j]<<' ';
    cout<<endl;
  }
  testcoverage();
  return 0;
}

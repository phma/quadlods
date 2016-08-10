/* Copyright 2014,2016 Pierre Abbat.
 * This file is part of the Quadlods program.
 * 
 * The Quadlods program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Quadlods is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License and Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and Lesser General Public License along with Quadlods. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <bitset>
#include "quadlods.h"
#include "ps.h"

using namespace std;

quadlods quads;

void plotxy(quadlods& quad,int xdim,int ydim)
{
  int i;
  double x,y;
  vector<double> point;
  startpage();
  setscale(0,0,1,1);
  // (2,0) and (3,0) look splotchy at 30000, but fill in well at 100000.
  for (i=0;i<30000;i++)
  {
    point=quad.dgen();
    x=point[xdim];
    y=point[ydim];
    dot(x,y);
  }
  endpage();
}

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
  psopen("quadlods.ps");
  psprolog();
  quads.init(5,1e10);
  quads.setjumble(QL_JUMBLE_GRAY);
  quads.advance(-1);
  for (i=0;i<5;i++)
    cout<<quads.getnum(i)<<'/'<<quads.getdenom(i)<<' '<<quads.getacc(i)<<endl;
  for (i=0;i<30;i++)
  {
    point=quads.dgen();
    for (j=0;j<5;j++)
      cout<<point[j]<<' ';
    cout<<endl;
  }
  //testcoverage();
  for (i=0;i<5;i++)
    for (j=0;j<i;j++)
      plotxy(quads,i,j);
  pstrailer();
  psclose();
  return 0;
}

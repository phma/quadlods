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
#include <set>
#include <map>
#include <ctime>
#include <cmath>
#include "quadlods.h"
#include "ps.h"
#include "discrepancy.h"

using namespace std;

quadlods quads;
char rpint[][2]=
{ // relatively prime integers whose sum of squares is less than 100
  {1,1},{1,2},{1,3},{1,4},{2,3},{1,5},{1,6},{2,5},
  {3,4},{1,7},{3,5},{1,8},{2,7},{4,5},{2,9},{3,8},
  {4,7},{5,6},{5,7},{4,9},{5,8},{6,7}
};

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

void testdiscrepancy(int dims,double res,int npoints)
{
  int i;
  quadlods q;
  double star;
  set<int> plotpoints;
  vector<vector<double> > seq;
  if (npoints<1)
    npoints=1;
  for (i=0;i<=300;i++)
    plotpoints.insert(rint(pow(npoints,i/300.))+1);
  q.init(dims,res);
  for (i=0;i<npoints;i++)
  {
    seq.push_back(q.dgen());
  }
  for (i=0;i<10;i++)
  {
    star=stardiscrepancy(seq);
    cout<<"Trial "<<i<<" discrepancy "<<star<<endl;
  }
}

void findclosequad()
{
  map<double,int> singleq;
  map<double,int>::iterator j;
  vector<double> qlist;
  int i,p,lastp,closep0,closep1;
  double q,lastq,closeness,closeq0,closeq1;
  for (i=0;i<QL_MAX_DIMS;i++)
  {
    p=nthprime(i);
    q=nthquad(i);
    singleq[q]=p;
  }
  for (j=singleq.begin();j!=singleq.end();j++)
    qlist.push_back(j->first);
  for (lastq==-3,closeness=1,j=singleq.begin();j!=singleq.end();j++)
  {
    p=j->second;
    q=j->first;
    if (q-lastq<closeness)
    {
      closeness=q-lastq;
      closep1=p;
      closep0=lastp;
      closeq1=q;
      closeq0=lastq;
    }
    lastq=q;
    lastp=p;
  }
  cout<<"The closest quads are "<<closeq0<<" ("<<closep0<<") and "<<closeq1<<" ("<<closep1<<")."<<endl;
}

/* Commands for interactive mode, which can be used as a server:
 * init n s res: Initialize generator #n with s dimensions and resolution res.
 * form n dec/hex/flo/rat: Set format to decimal/hexadecimal/floating point/rational.
 * gene n i: Generate i points from generator n.
 * seed n: Seed generator n with random numbers.
 */

int main(int argc,char **argv)
{
  int i,j,seedlen;
  char *seedbuf;
  time_t now;
  vector<double> point;
  vector<int> badprimes; // none of the primes is bad per se, they just have very close q values
  badprimes.push_back(65029);
  badprimes.push_back(65027);
  badprimes.push_back(64013);
  badprimes.push_back(64007);
  badprimes.push_back(62003);
  psopen("quadlods.ps");
  psprolog();
  quads.init(badprimes,1e10);
  quads.setjumble(QL_JUMBLE_GRAY);
  quads.advance(-1);
  for (i=0;i<5;i++)
    cout<<quads.getnum(i)<<'/'<<quads.getdenom(i)<<' '<<quads.getacc(i)<<endl;
  seedlen=quads.seedsize();
  cout<<seedlen<<" bytes needed to seed"<<endl;
  seedbuf=new char[seedlen];
  now=time(NULL);
  for (i=0;i<sizeof(time_t) && i<seedlen;i++)
    seedbuf[i]=now>>(8*i);
  for (;i<seedlen;i++)
    seedbuf[i]=seedbuf[i-1]+seedbuf[i-sizeof(time_t)];
  quads.seed(seedbuf,seedlen);
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
  //testdiscrepancy(5,1e10,1000);
  findclosequad();
  psclose();
  delete[] seedbuf;
  return 0;
}

/* Copyright 2014,2016,2018 Pierre Abbat.
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
#include <cassert>
#include <cmath>
#include "quadlods.h"
#include "ps.h"
#include "discrepancy.h"
#include "circletest.h"
#include "filltest.h"
#include "contfrac.h"
#include "ldecimal.h"
#include "histogram.h"

using namespace std;

struct command
{
  std::string word;
  void (*fun)();
  std::string desc;
  command(std::string w,void (*f)(),std::string d)
  {
    word=w;
    fun=f;
    desc=d;
  }
};

vector<command> commands;
quadlods quads,cirquads;
char rpint[][2]=
{ // relatively prime integers whose sum of squares is less than 100
  {1,1},{1,2},{1,3},{1,4},{2,3},{1,5},{1,6},{2,5},
  {3,4},{1,7},{3,5},{1,8},{2,7},{4,5},{2,9},{3,8},
  {4,7},{5,6},{5,7},{4,9},{5,8},{6,7}
};
map<double,int> singleq;
int minlargerp;
PostScript ps;

void plotxy(quadlods& quad,int xdim,int ydim)
{
  int i;
  double x,y;
  vector<double> point;
  ps.startpage();
  ps.setscale(0,0,1,1);
  // (2,0) and (3,0) look splotchy at 30000, but fill in well at 100000.
  for (i=0;i<30000;i++)
  {
    point=quad.dgen();
    x=point[xdim];
    y=point[ydim];
    ps.dot(x,y);
  }
  ps.endpage();
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
  vector<int> p235; // since sorting the primes, the list begins 5,3,7
  p235.push_back(2);
  p235.push_back(3);
  p235.push_back(5);
  cov.init(p235,370);
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

void test2quads(double q0,double q1,double toler)
{
  int i,a=0,b=0;
  int p0,p1;
  double diff;
  for (i=0;i<sizeof(rpint)/sizeof(rpint[0]);i++)
  {
    a=rpint[i][0];
    b=rpint[i][1];
    diff=a*q0-b*q1;
    diff-=rint(diff);
    if (fabs(diff)<toler)
      break;
    a=-a;
    diff=a*q0-b*q1;
    diff-=rint(diff);
    if (fabs(diff)<toler)
      break;
    if (i==0)
      continue;
    swap(a,b);
    diff=a*q0-b*q1;
    diff-=rint(diff);
    if (fabs(diff)<toler)
      break;
    b=-b;
    diff=a*q0-b*q1;
    diff-=rint(diff);
    if (fabs(diff)<toler)
      break;
    a=b=0;
  }
  if (a)
  {
    p0=singleq[q0];
    p1=singleq[q1];
    cout<<"Close quads are "<<a<<'*'<<q0<<" ("<<p0<<") and "<<b<<'*'<<q1<<" ("<<p1<<")."<<endl;
    if (p0<p1)
      swap(p0,p1);
    if (p0<minlargerp)
      minlargerp=p0;
  }
}

void findclosequad()
{
  map<double,int>::iterator j;
  vector<double> qlist;
  int i,k,p,lastp,closep0,closep1;
  double q,lastq,closeness,closeq0,closeq1;
  singleq.clear();
  minlargerp=65536;
  for (i=0;i<QL_MAX_DIMS;i++)
  {
    p=nthprime(i);
    q=nthquad(i);
    singleq[q]=p;
  }
  for (j=singleq.begin();j!=singleq.end();j++)
    qlist.push_back(j->first);
  for (i=0;i<QL_MAX_DIMS;i++)
    for (k=0;k<i;k++)
      test2quads(qlist[i],qlist[k],2e-9);
  cout<<"Use primes below "<<minlargerp<<endl;
}

void testBadPrimes()
{
  int i,j;
  vector<int> badprimes; // none of the primes is bad per se, they just have very close q values
  badprimes.push_back(65027);
  badprimes.push_back(62003);
  badprimes.push_back(59051);
  badprimes.push_back(56171);
  badprimes.push_back(50627);
  badprimes.push_back(47963);
  badprimes.push_back(65029);
  badprimes.push_back(64013);
  badprimes.push_back(60029);
  badprimes.push_back(59053);
  ps.open("quadlods.ps");
  ps.prolog();
  quads.init(badprimes,1e10);
  quads.setjumble(QL_JUMBLE_GRAY);
  quads.advance(-1);
  for (i=0;i<quads.size();i++)
    cout<<quads.getnum(i)<<'/'<<quads.getdenom(i)<<' '<<quads.getacc(i)<<endl;
  for (i=0;i<quads.size();i++)
    for (j=0;j<i;j++)
      plotxy(quads,i,j);
  ps.trailer();
  //testdiscrepancy(5,1e10,1000);
  cirquads.init(badprimes,1e17,QL_JUMBLE_NONE);
  circletest(cirquads);
  filltest(cirquads);
  ps.close();
}

void testGoodPrimes()
{
  int i,j;
  ps.open("quadlods.ps");
  ps.prolog();
  quads.init(10,1e17);
  quads.setjumble(QL_JUMBLE_GRAY);
  quads.advance(-1);
  for (i=0;i<quads.size();i++)
    cout<<quads.getnum(i)<<'/'<<quads.getdenom(i)<<' '<<quads.getacc(i)<<endl;
  for (i=0;i<quads.size();i++)
    for (j=0;j<i;j++)
      plotxy(quads,i,j);
  ps.trailer();
  //testdiscrepancy(5,1e10,1000);
  cirquads.init(10,1e17,QL_JUMBLE_NONE);
  circletest(cirquads);
  filltest(cirquads);
  ps.close();
}

void testSeed()
{
  int i,j,seedlen;
  vector<double> point;
  char *seedbuf;
  time_t now;
  quads.init(5,1e10);
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
    for (j=0;j<point.size();j++)
      cout<<point[j]<<' ';
    cout<<endl;
  }
  delete[] seedbuf;
}

void writeshort(ostream &file,unsigned short i)
{
  char buf[2];
  *(unsigned short *)buf=i;
  file.write(buf,2);
}

void sortPrimes()
{
  int i,j;
  const double logKhinchin=log(2.685452001065306);
  quads.init(0,0);
  ofstream primeFile("primes.dat",ios::binary);
  ofstream primeText("primes.txt");
  histogram hist(0,6.235); // log scale, 1 to 510
  PostScript ps;
  ContinuedFraction cf;
  set<PrimeContinuedFraction> pcf;
  set<PrimeContinuedFraction>::iterator k;
  PrimeContinuedFraction pcf0;
  for (i=0;i<QL_MAX_DIMS;i++)
  {
    //cout<<nthprime(i)<<':';
    cf=contFrac(nthquadQi(i));
    for (j=0;false && j<cf.terms.size();j++)
    {
      if (j+cf.period==cf.terms.size())
	cout<<'(';
      else
	cout<<' ';
      cout<<cf.terms[j];
    }
    //cout<<')'<<cf.averageTerm()<<endl;
    pcf0.prime=nthprime(i);
    pcf0.cf=cf;
    pcf.insert(pcf0);
    hist<<log(pcf0.cf.averageTerm());
  }
  for (k=pcf.begin();k!=pcf.end();k++)
  {
    writeshort(primeFile,k->prime);
    writeshort(primeFile,k->cf.terms.size());
    writeshort(primeFile,k->cf.period);
    for (i=0;i<k->cf.terms.size();i++)
      writeshort(primeFile,k->cf.terms[i]);
    primeText<<k->prime<<' '<<ldecimal(k->cf.averageTerm());
    for (i=j=0;i<40;i++,j++)
    {
      if (j>=k->cf.terms.size())
	j-=k->cf.period;
      primeText<<' '<<k->cf.terms[j];
    }
    primeText<<" ...\n";
  }
  ps.open("primes.ps");
  ps.prolog();
  ps.setpaper(a4land,0);
  ps.startpage();
  hist.plot(ps,HISTO_LOG);
  ps.setcolor(1,0,0);
  ps.line2p(xy(logKhinchin*3/6.235,2),xy(logKhinchin*3/6.235,2.1));
}

/* Commands for interactive mode, which can be used as a server:
 * init n s res: Initialize generator #n with s dimensions and resolution res.
 * form n dec/hex/flo/rat: Set format to decimal/hexadecimal/floating point/rational.
 * gene n i: Generate i points from generator n.
 * seed n: Seed generator n with random numbers.
 */

int main(int argc,char **argv)
{
  string arg1;
  int cmd,i;
  if (argc>1)
    arg1=argv[1];
  commands.push_back(command("sortprimes",sortPrimes,"Sort primes by average continued fraction term"));
  commands.push_back(command("coverage",testcoverage,"Test coverage of small 3D generator"));
  commands.push_back(command("goodprimes",testGoodPrimes,"Test primes with low CF terms"));
  commands.push_back(command("badprimes",testBadPrimes,"Test primes with low CF terms"));
  for (cmd=-1,i=0;i<commands.size();i++)
    if (commands[i].word==arg1)
      cmd=i;
  if (cmd>=0)
    commands[cmd].fun();
  //testBadPrimes();
  //testSeed();
  //findclosequad();
  return 0;
}

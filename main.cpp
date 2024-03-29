/******************************************************/
/*                                                    */
/* main.cpp - main program                            */
/*                                                    */
/******************************************************/
/* Copyright 2014,2016,2018-2021,2023 Pierre Abbat.
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
#include <ctime>
#include <cassert>
#include <cmath>
#include <boost/program_options.hpp>
#include "main.h"
#include "config.h"
#include "hstep.h"
#include "ps.h"
#include "threads.h"
#include "random.h"
#include "circletest.h"
#include "filltest.h"
#include "contfrac.h"
#include "ldecimal.h"
#include "histogram.h"
#include "manysum.h"
#include "flowertest.h"
#include "fourier.h"
#include "plot.h"
#include "pairpoint.h"
#include "matrix.h"
#include "interact.h"
#include "discrepancy.h"

#define tassert(x) testfail|=(!(x))

using namespace std;
using namespace quadlods;
namespace po=boost::program_options;

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
Quadlods cirquads;
map<int,Quadlods> quads;
bool testfail=false;
char rpint[][2]=
{ // relatively prime integers whose sum of squares is less than 100
  {1,1},{1,2},{1,3},{1,4},{2,3},{1,5},{1,6},{2,5},
  {3,4},{1,7},{3,5},{1,8},{2,7},{4,5},{2,9},{3,8},
  {4,7},{5,6},{5,7},{4,9},{5,8},{6,7}
};
map<double,int> singleq;
int minlargerp;
PostScript ps;
string resstr;
double resolution;
string primestr;
vector<int> primelist;
string scramblestr;
int ndims,niter,scramble;
string filename;
bool useMinMax=false,disc2d=false;
int maxStairStep;

void listCommands()
{
  int i;
  cout<<"Quadlods version "<<VERSION<<" © "<<COPY_YEAR<<" Pierre Abbat\n"<<
    "Quadratic-irrational (Richtmyer) and Halton low-discrepancy sequences\n"<<
    "Distributed under GPL v3 or later; library distributed under LGPL.\n"<<
    "This is free software with no warranty.\n";
  cout<<"Commands:\n";
  for (i=0;i<commands.size();i++)
    cout<<setw(15)<<left<<commands[i].word<<commands[i].desc<<'\n';
}

int parseInt(string intstr)
{
  int ret;
  size_t idx;
  ret=stoi(intstr,&idx);
  if (idx<intstr.length())
    throw(invalid_argument("parseInt"));
  return ret;
}

double parseDouble(string intstr)
{
  double ret;
  size_t idx;
  ret=stod(intstr,&idx);
  if (idx<intstr.length())
    throw(invalid_argument("parseDouble"));
  return ret;
}

bool isValidPrime(int n)
{
  int i;
  bool ret=n>1 && n<65522 && (n==2 || (n&1));
  for (i=3;ret && i*i<=n;i+=2)
    if (n>i && (n%i==0))
      ret=false;
  return ret;
}

void testContinuedFraction()
{
  ContinuedFraction cf2,cf3,cf5,cf7;
  cout<<"Continued fraction unit test\n";
  cf2.terms.push_back(1);
  cf2.terms.push_back(2);
  cf2.period=1;
  cf3.terms.push_back(1);
  cf3.terms.push_back(1);
  cf3.terms.push_back(2);
  cf3.period=2;
  cf5.terms.push_back(1);
  cf5.period=1;
  cf7.terms.push_back(2);
  cf7.terms.push_back(1);
  cf7.terms.push_back(1);
  cf7.terms.push_back(1);
  cf7.terms.push_back(4);
  cf7.period=4;
  tassert(fabs(cf3.averageTerm()-cf7.averageTerm())<1e-15);
  tassert(cf3.maximumTerm()<cf7.maximumTerm());
  tassert(cf2.averageTerm()==2);
  tassert(cf2.maximumTerm()==2);
  tassert(cf5.averageTerm()==1);
  tassert(cf5.maximumTerm()==1);
}

void testReverseScramble()
{
  int i,j;
  cout<<"Reverse scramble unit test\n";
  fillReverseScrambleTable(5,QL_SCRAMBLE_NONE);
  fillReverseScrambleTable(5,QL_SCRAMBLE_POWER);
  tassert(reverseScrambleTable[(QL_SCRAMBLE_NONE<<16)+5][10646]==5638);
  tassert(reverseScrambleTable[(QL_SCRAMBLE_POWER<<16)+5][10646]==5642);
  for (i=2;i<12;i++)
  {
    for (j=0;j<i;j++)
      cout<<faureperm(j,i)<<' ';
    cout<<endl;
  }
}

vector<int> factor(int n)
{
  vector<int> ret;
  int i,p;
  for (i=0;i<QL_MAX_DIMS && n>1;i++)
  {
    p=nthprime(i);
    while (n%p==0)
    {
      ret.push_back(p);
      n/=p;
    }
  }
  return ret;
}

bool parsePrimeList()
{
  string numstr;
  int num;
  size_t pos;
  bool ret=true;
  while (primestr.length())
  {
    pos=primestr.find_first_of(", ");
    numstr=primestr.substr(0,pos);
    if (pos<=primestr.length())
      pos++;
    primestr.erase(0,pos);
    if (numstr.length())
    {
      if (numstr.find_first_not_of("0123456789")>numstr.length())
      {
	try
	{
	  num=parseInt(numstr);
	  if (isValidPrime(num))
	    primelist.push_back(num);
	  else
	  {
	    cerr<<num<<" is not a prime less than 65535\n";
	    ret=false;
	  }
	}
	catch (...)
	{
	  cerr<<numstr<<" is out of range for an int\n";
	  ret=false;
	}
      }
      else
      {
	cerr<<numstr<<" is not a positive integer\n";
	ret=false;
      }
    }
  }
  return ret;
}

array<short,676> digraphs(string word)
{
  int l0,l1,i;
  array<short,676> ret;
  for (i=0;i<676;i++)
    ret[i]=0;
  for (i=0;i+1<word.length();i++)
  {
    l0=tolower(word[i]);
    l1=tolower(word[i+1]);
    if (l0>='a' && l0<='z' && l1>='a' && l1<='z')
      ret[(l0-'a')*26+(l1-'a')]++;
  }
  return ret;
}

int parseScramble(string scramblestr)
{
  vector<array<short,676> > digs;
  array<short,676> digj=digraphs(scramblestr);
  vector<int> match;
  int i,j,ret,maxmatch=-1,nmatch=0;
  digs.push_back(digraphs("None"));
  digs.push_back(digraphs("Third"));
  digs.push_back(digraphs("Thue-Morse"));
  digs.push_back(digraphs("Gray"));
  digs.push_back(digraphs("Power"));
  digs.push_back(digraphs("Faure"));
  digs.push_back(digraphs("tipwitch"));
  digs.push_back(digraphs("default"));
  match.resize(digs.size());
  for (i=0;i<676;i++)
    for (j=0;j<digs.size();j++)
      match[j]+=digs[j][i]*digj[i];
  for (i=0;i<match.size();i++)
    if (match[i]>maxmatch)
      maxmatch=match[i];
  if (match[0]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_NONE;
  }
  if (match[1]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_THIRD;
  }
  if (match[2]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_THUEMORSE;
  }
  if (match[3]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_GRAY;
  }
  if (match[4]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_POWER;
  }
  if (match[5]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_FAURE;
  }
  if (match[6]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_TIPWITCH;
  }
  if (match[7]==maxmatch)
  {
    nmatch++;
    ret=QL_SCRAMBLE_DEFAULT;
  }
  if (nmatch>1)
  {
    ret=-1;
  }
  return ret;
}

double parseResolution(string resstr)
{
  size_t rest;
  double ret;
  try
  {
    ret=stod(resstr,&rest);
    if (rest<resstr.length())
      ret=-1;
  }
  catch (...)
  {
    int i;
    string lower;
    for (i=0;i<resstr.length();i++)
      lower+=tolower(resstr[i]);
    if (lower==string("halton").substr(0,lower.length()) && lower.length()>0)
      ret=0;
    else
      ret=-1;
  }
  return ret;
}

void plotxy(Quadlods& quad,int xdim,int ydim,int inx,int allinx)
{
  int i;
  double x,y;
  time_t now,then;
  char buf[24];
  vector<double> point;
  Quadlods sel2;
  PairCompressor pc;
  vector<int> pinx2;
  pinx2.push_back(xdim);
  pinx2.push_back(ydim);
  sel2=select(quad,pinx2);
  ps.startpage();
  ps.setscale(0,0,1,1);
  // (2,5) and (7,2) ((2,0) and (5,2)) look splotchy at 30000, but fill in well at 100000.
  sprintf(buf,"%d %d",sel2.getprime(0),sel2.getprime(1));
  ps.write(0,1,buf);
  for (i=0;i<niter;i++)
  {
    point=sel2.dgen();
    x=point[0];
    y=point[1];
    pc.insert(xy(x,y));
    now=time(nullptr);
    if (now!=then)
    {
      cout<<rint((inx+(double)i/niter)/allinx*100)<<"% \r";
      cout.flush();
      then=now;
    }
  }
  ps.draw(pc);
  ps.endpage();
}

void testcoverage()
/* Cycle sizes: 408, 571, 377.
 * Total cycle size: 87828936 (LCM).
 * Remaining factor: 1.
 */
{
  int i,x,y,z,inx;
  Quadlods cov;
  bitset<87828936> *histo;
  vector<mpq_class> point;
  vector<int> p235; // since sorting the primes, the list begins 5,3,2
  cout<<"Testing coverage\n";
  p235.push_back(2);
  p235.push_back(3);
  p235.push_back(5);
  cov.init(p235,370);
  cov.setscramble(scramble);
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
    cout<<"Complete coverage\n";
  else
  {
    testfail=true;
    cout<<"Coverage test failed\n";
  }
  delete histo;
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
  int i,k,p;
  double q;
  singleq.clear();
  minlargerp=65536;
  for (i=0;i<QL_MAX_DIMS;i++)
  {
    p=nthprime(i);
    q=nthquad(i);
    singleq[q]=p;
  }
  for (j=singleq.begin();j!=singleq.end();++j)
    qlist.push_back(j->first);
  for (i=0;i<QL_MAX_DIMS;i++)
    for (k=0;k<i;k++)
      test2quads(qlist[i],qlist[k],2e-9);
  cout<<"Use primes below "<<minlargerp<<endl;
}

/* Some bad pairs of primes for Richtmyer:
 * 65027,65029
 * 59053,59051
 * 13691,13693
 * A set that is somewhat bad:
 * 5,4933,42533,57713 (CF expansions all begin n;1,1,1,1,1,1,1,1)
 * 42533,57713 is bad, but not as bad as 65027,65029.
 * Twin primes, 2 and 4 more than a square:
 * 3,5 227,229 1091,1093 2027,2029 3251,3253
 * 13691,13693 21611,21613 59051,59053 65027,65029.
 * All these have CF expansions of the form a;(a,b) and c;(a)
 * where b=2a and a=2c-1.
 * Any large enough pair of twin primes is bad for Halton, and even
 * when Faure-scrambled, such pairs as 21839,21841 (close to 65536/3)
 * are bad. So are pairs whose ratio is close to 2, such as 63277,31643.
 *
 * 44741,63277 is a good pair for Faure; their ratio is very nearly √2.
 */

void testCircle()
{
  if (niter>0 && (ndims>0 || primelist.size()))
  {
    ps.open(filename.length()?filename:"circle.ps");
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    quads[0].advance(-1);
    circletest(quads[0],niter,ps);
    ps.trailer();
    ps.close();
  }
  if (niter<=0)
    cerr<<"Please specify number of iterations with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void testScatter()
{
  int i,j,inx,allinx;
  time_t now,then;
  if (niter>0 && (ndims>0 || primelist.size()))
  {
    ps.open(filename.length()?filename:"scatter.ps");
    ps.prolog();
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    quads[0].advance(-1);
    allinx=quads[0].size()*(quads[0].size()-1)/2;
    //for (i=0;resolution && i<quads[0].size();i++)
      //cout<<quads[0].getnum(i)<<'/'<<quads[0].getdenom(i)<<' '<<quads[0].getacc(i)<<endl;
    for (i=0;i<quads[0].size();i++)
      for (j=0;j<i;j++)
      {
	inx=i*(i-1)/2+j;
	plotxy(quads[0],i,j,inx,allinx);
	now=time(nullptr);
	if (now!=then)
	{
	  cout<<rint((double)inx/allinx*100)<<"% \r";
	  cout.flush();
	  then=now;
	}
      }
    ps.trailer();
    ps.close();
  }
  if (niter<=0)
    cerr<<"Please specify number of iterations with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void testFill()
{
  if (niter>0 && (ndims>0 || primelist.size()))
  {
    ps.open(filename.length()?filename:"fill.ps");
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    quads[0].advance(-1);
    filltest(quads[0],niter,ps);
    ps.trailer();
    ps.close();
  }
  if (niter<=0)
    cerr<<"Please specify number of iterations with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void testFlower()
{
  if (niter>0 && (ndims>0 || primelist.size()))
  {
    ps.open(filename.length()?filename:"flower.ps");
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    quads[0].advance(-1);
    flowertest(quads[0],niter,ps,disc2d);
    ps.trailer();
    ps.close();
  }
  if (niter<=0)
    cerr<<"Please specify number of iterations with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void testFourier()
{
  if (niter>0 && (ndims>0 || primelist.size()))
  {
    ps.open(filename.length()?filename:"fourier.ps");
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    quads[0].advance(-1);
    fouriertest(quads[0],niter,ps);
    ps.trailer();
    ps.close();
  }
  if (niter<=0)
    cerr<<"Please specify number of iterations with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void quadPlot()
{
  ps.open(filename.length()?filename:"quadplot.ps");
  quads[0].init(ndims,resolution);
  quads[0].init(primelist,resolution);
  quads[0].setscramble(scramble);
  quads[0].advance(-1);
  quadplot(ps);
  ps.trailer();
  ps.close();
}

void testRandom()
{
  int hist[256],i;
  int done=0,max,min,maxstep=0;
  memset(hist,0,sizeof(hist));
  while (!done)
  {
    hist[rng.ucrandom()]++;
    for (max=i=0,min=16777777;i<256;i++)
    {
      if (hist[i]>max)
	max=hist[i];
      if (hist[i]<min)
	min=hist[i];
    }
    if (max>16777215)
      done=-1;
    if (max-min>1.1*pow(max,1/3.) && max-min<0.9*pow(max,2/3.))
      done=1;
    if (max-maxstep>=16384)
    {
      maxstep=max;
      //cout<<max<<' '<<min<<endl;
    }
  }
  tassert(done==1);
  cout<<"Random test: max "<<max<<" min "<<min<<endl;
  done=maxstep=0;
  memset(hist,0,sizeof(hist));
  while (!done)
  {
    hist[rng.rangerandom(251)]++;
    for (max=i=0,min=16777777;i<251;i++)
    {
      if (hist[i]>max)
	max=hist[i];
      if (hist[i]<min)
	min=hist[i];
    }
    if (max>16777215)
      done=-1;
    if (max-min>1.1*pow(max,1/3.) && max-min<0.9*pow(max,2/3.))
      done=1;
    if (max-maxstep>=16384)
    {
      maxstep=max;
      //cout<<max<<' '<<min<<endl;
    }
  }
  tassert(done==1);
  cout<<"Range random test: max "<<max<<" min "<<min<<endl;
  done=maxstep=0;
  memset(hist,0,sizeof(hist));
  while (!done)
  {
    hist[rng.frandom(mpq_class(1,3))]++;
    for (max=i=0,min=16777777;i<2;i++)
    {
      if (hist[i]*(1+i)>max)
	max=hist[i]*(1+i);
      if (hist[i]*(1+i)<min)
	min=hist[i]*(1+i);
    }
    if (max>16777215)
      done=-1;
    if (min>256 && max-min>1.1*pow(max,1/3.) && max-min<0.9*pow(max,2/3.))
      done=1;
    if (max-maxstep>=16384)
    {
      maxstep=max;
      //cout<<max<<' '<<min<<endl;
    }
  }
  tassert(done==1);
  cout<<"Fractional random test: "<<hist[0]<<" zeros "<<hist[1]<<" ones\n";
}

void testSeed()
{
  int n,i,j,seedlen;
  vector<double> point;
  char *seedbuf;
  time_t now;
  quads[0].init(5,1e10);
  quads[1].init(5,0);
  for (n=0;n<2;n++)
  {
    seedlen=quads[n].seedsize();
    cout<<seedlen<<" bytes needed to seed"<<endl;
    seedbuf=new char[seedlen];
    now=time(NULL);
    for (i=0;i<sizeof(time_t) && i<seedlen;i++)
      seedbuf[i]=now>>(8*i);
    for (;i<seedlen;i++)
      seedbuf[i]=seedbuf[i-1]+seedbuf[i-sizeof(time_t)];
    quads[n].seed(seedbuf,seedlen);
    for (i=0;i<30;i++)
    {
      point=quads[n].dgen();
      for (j=0;j<point.size();j++)
	cout<<point[j]<<' ';
      cout<<endl;
    }
    delete[] seedbuf;
  }
}

void testNegativeHalton()
{
  int i,j;
  vector<double> point,sum;
  quads[0].init(5,0,QL_SCRAMBLE_NONE);
  quads[0].advance(-16);
  sum.resize(5,0);
  cout<<"Negative Halton test\n";
  for (i=0;i<30;i++)
  {
    point=quads[0].dgen();
    for (j=0;j<point.size();j++)
    {
      sum[j]+=point[j];
      cout<<point[j]<<' ';
    }
    cout<<endl;
  }
  for (j=0;j<sum.size();j++)
  {
    tassert(fabs(sum[j]-15)<1e-12);
  }
}

array<double,3> spherePoint(double x0,double x1)
{
  array<double,3> ret;
  double r;
  ret[2]=2*x1-1;
  r=sqrt(1-sqr(ret[2]));
  ret[1]=r*sin(x0*2*M_PI);
  ret[0]=r*cos(x0*2*M_PI);
  return ret;
}

double upperDiff(int n)
/* For a six-dimensional integration, the error in the sum should be bounded
 * by the sixth power of the logarithm. For small n, though, it isn't.
 */
{
  return pow(log(n),6)/5e4+pow(log(n),4)/180+pow(log(n),3)/128+log(n)*2.22;
}

void testuvmatrix()
// Calculate the distribution of the determinant of three unit 3-vectors
{
  int i,j,k;
  double det,reldiff,maxreldiff=0;
  time_t now,then;
  bool analyze=false;
  vector<double> allsqdet;
  double absdiff,minabsdiff,maxabsdiff;
  double hi=-INFINITY,lo=INFINITY,scale;
  manysum sumsqdet;
  histogram hist(-1,1);
  matrix mat(3,3);
  array<double,3> row;
  array<double,4> vline;
  vector<array<double,4> > vlines;
  vector<double> point;
  if (niter==0)
  {
    analyze=true;
    niter=7776000;
  }
  quads[0].init(6,resolution);
  quads[0].setscramble(scramble);
  cout<<"Calculating determinants of 3×3 unit vector matrices\n";
  for (i=0;i<niter;i++)
  {
    point=quads[0].dgen();
    row=spherePoint(point[0],point[1]);
    for (j=0;j<3;j++)
      mat[0][j]=row[j];
    row=spherePoint(point[2],point[3]);
    for (j=0;j<3;j++)
      mat[1][j]=row[j];
    row=spherePoint(point[4],point[5]);
    for (j=0;j<3;j++)
      mat[2][j]=row[j];
    det=mat.determinant();
    now=time(nullptr);
    if (now!=then)
    {
      cout<<rint((double)i/niter*100)<<"% \r";
      cout.flush();
      reldiff=fabs(sumsqdet.total()-i*6./27.)/pow(log(i),6);
      if (reldiff>maxreldiff)
	maxreldiff=reldiff;
      then=now;
    }
    hist<<det;
    sumsqdet+=det*det;
    if (analyze)
      allsqdet.push_back(det*det);
  }
  ps.open(filename.length()?filename:"uvmatrix.ps");
  ps.setpaper(a4land,0);
  ps.prolog();
  ps.startpage();
  ps.setscale(0,-1,3,1);
  hist.plot(ps,HISTO_LINEAR);
  ps.endpage();
  cout<<niter<<" iterations, average square determinant is "<<sumsqdet.total()/niter<<endl;
  //cout<<"compare "<<sumsqdet.total()-niter*6./27.<<' '<<upperDiff(niter)<<endl;
  if (niter<2)
    cout<<"Please increase the number of iterations with -n\n";
  else if (fabs(sumsqdet.total()-niter*6./27.)>upperDiff(niter))
  {
    cout<<"Test failed, average square determinant should be 6/27\n";
    // 6/27 is 3!/3^3
    testfail=true;
  }
  //cout<<"Maximum relative difference is "<<maxreldiff<<endl;
  if (analyze)
  {
    for (i=2;i<2;i++)
    {
      minabsdiff=INFINITY;
      maxabsdiff=0;
      if (niter%i==0)
      {
	for (j=0;j<niter;j+=i)
	{
	  sumsqdet.clear();
	  for (k=0;k<i;k++)
	    sumsqdet+=allsqdet[j+k];
	  absdiff=fabs(sumsqdet.total()-i*6./27.);
	  if (absdiff<minabsdiff)
	    minabsdiff=absdiff;
	  if (absdiff>maxabsdiff)
	    maxabsdiff=absdiff;
	}
	cout<<i<<' '<<minabsdiff<<'-'<<maxabsdiff<<" < "<<upperDiff(i)<<endl;
	vline[0]=log(i);
	vline[1]=log(minabsdiff);
	vline[2]=maxabsdiff;
	vline[3]=upperDiff(i);
	vlines.push_back(vline);
	for (k=1;k<4;k++)
	{
	  if (vline[k]>hi)
	    hi=vline[k];
	  if (vline[k]<lo)
	    lo=vline[k];
	}
      }
    }
    if (vlines.size())
    {
      ps.startpage();
      ps.setscale(0,-1,3,1);
      scale=(hi-lo)/2;
      ps.startline();
      ps.lineto(0,-1);
      ps.lineto(3,-1);
      ps.lineto(3,1);
      ps.lineto(0,1);
      ps.endline(true);
      ps.setcolor(0,0,1);
      for (i=0;i<vlines.size();i++)
	ps.line2p(xy(vlines[i][0]/log(niter)*3,(vlines[i][1]-lo)/scale-1),
		  xy(vlines[i][0]/log(niter)*3,(vlines[i][2]-lo)/scale-1));
      ps.setcolor(1,0,0);
      ps.startline();
      for (i=0;i<vlines.size();i++)
	ps.lineto(vlines[i][0]/log(niter)*3,(vlines[i][3]-lo)/scale-1);
      ps.endline();
      ps.endpage();
    }
    niter=0;
  }
  ps.close();
}

void outHaltonAccumulators(vector<unsigned short> terHacc,vector<unsigned short> quinHacc,vector<unsigned short> septHacc,bool sign)
{
  int i;
  for (i=terHacc.size()-1;i>=0;i--)
    cout<<terHacc[i]<<' ';
  cout<<endl<<haccValue(terHacc,59049,sign)<<endl;
  for (i=quinHacc.size()-1;i>=0;i--)
    cout<<quinHacc[i]<<' ';
  cout<<endl<<haccValue(quinHacc,15625,sign)<<endl;
  for (i=septHacc.size()-1;i>=0;i--)
    cout<<septHacc[i]<<' ';
  cout<<endl<<haccValue(septHacc,16807,sign)<<endl;
}

void testHaltonAccumulator()
{
  vector<unsigned short> terHacc,quinHacc,septHacc;
  int i;
  mpz_class inc;
  bool sign=false,newsign;
  cout<<"Halton accumulator test\n";
  newsign=incHacc(terHacc,59049,496125,0,sign);
  tassert(!newsign);
  newsign=incHacc(quinHacc,15625,496125,0,sign);
  tassert(!newsign);
  newsign=incHacc(septHacc,16807,496125,0,sign);
  tassert(!newsign);
  outHaltonAccumulators(terHacc,quinHacc,septHacc,newsign);
  sign=newsign;
  newsign=incHacc(terHacc,59049,-524288,0,sign);
  tassert(newsign);
  newsign=incHacc(quinHacc,15625,-524288,0,sign);
  tassert(newsign);
  newsign=incHacc(septHacc,16807,-524288,0,sign);
  tassert(newsign);
  outHaltonAccumulators(terHacc,quinHacc,septHacc,newsign);
  sign=newsign;
  for (inc=1,i=0;i<41;i++)
    inc*=41;
  newsign=incHacc(terHacc,59049,inc,sign);
  tassert(!newsign);
  newsign=incHacc(quinHacc,15625,inc,sign);
  tassert(!newsign);
  newsign=incHacc(septHacc,16807,inc,sign);
  tassert(!newsign);
  outHaltonAccumulators(terHacc,quinHacc,septHacc,newsign);
  tassert(terHacc[0]==41670);
  tassert(terHacc.size()==14);
  tassert(terHacc.back()==12542);
  tassert(quinHacc[0]==3603);
  tassert(quinHacc.size()==16);
  tassert(quinHacc.back()==1647);
  tassert(septHacc[0]==4036);
  tassert(septHacc.size()==16);
  tassert(septHacc.back()==551);
}

void test1AreaInCircle(double minx,double miny,double maxx,double maxy,double rightArea)
{
  double area=areaInCircle(minx,miny,maxx,maxy);
  cout<<"Area="<<ldecimal(area);
  if (rightArea>=0)
  {
    cout<<", should be "<<ldecimal(rightArea);
    tassert(fabs(area-rightArea)<1e-9);
  }
  cout<<endl;
}

void test4AreaInCircle(double x0,double x1,double x2,double y0,double y1,double y2)
{
  double area[5];
  double diff;
  int i;
  if (x0>x1)
    swap(x0,x1);
  if (x1>x2)
    swap(x1,x2);
  if (x0>x1)
    swap(x0,x1);
  if (y0>y1)
    swap(y0,y1);
  if (y1>y2)
    swap(y1,y2);
  if (y0>y1)
    swap(y0,y1);
  area[4]=areaInCircle(x0,y0,x2,y2);
  area[0]=areaInCircle(x0,y0,x1,y1);
  area[1]=areaInCircle(x0,y1,x1,y2);
  area[2]=areaInCircle(x1,y0,x2,y1);
  area[3]=areaInCircle(x1,y1,x2,y2);
  diff=area[4]-area[0]-area[1]-area[2]-area[3];
  if (fabs(diff)>1e-9)
  {
    cout<<"x "<<x0<<' '<<x1<<' '<<x2<<" y "<<y0<<' '<<y1<<' '<<y2<<endl;
    cout<<"Areas differ by "<<diff<<endl;
    testfail=true;
  }
  for (i=0;i<5;i++)
  {
    tassert(area[i]>=0);
    tassert(area[i]<355/113.);
  }
}

double pythFrac[]=
{ // In double, (5,12,13) is just greater than 1, while (9,40,41) is just less than 1.
  -1,
  -40/41.,-24/25.,-35/37.,-12/13.,-15/17.,-4/5.,-21/29.,
  -20/29.,-3/5.,-8/17.,-5/13.,-12/37.,-7/25.,-9/41.,
  0,
  9/41.,7/25.,12/37.,5/13.,8/17.,3/5.,20/29.,
  21/29.,4/5.,15/17.,12/13.,35/37.,24/25.,40/41.,
  1
};

void testAreaInCircle()
{
  int i,j,k,l;
  vector<double> point;
  test1AreaInCircle(-0.8,-0.8,0.8,0.8,2.487588218416656);
  test1AreaInCircle(-1,-1,-0.8,-0.8,0);
  test1AreaInCircle(0.8,0.6,1,1,0);
  quads[0].init(2,resolution);
  quads[0].setscramble(scramble);
  for (i=1;i<31;i++)
    for (j=0;j<i;j++)
      for (k=1;k<31;k++)
	for (l=0;l<k;l++)
	{
	  point=quads[0].dgen();
	  point[0]=2*point[0]-1;
	  point[1]=2*point[1]-1;
	  test4AreaInCircle(pythFrac[i],point[0],pythFrac[j],
			    pythFrac[k],point[1],pythFrac[l]);
	}
}

void runTests()
{
  testContinuedFraction();
  testReverseScramble();
  testNegativeHalton();
  testRandom();
  testSeed();
  testHaltonAccumulator();
  testAreaInCircle();
}

void runLongTests()
{
  testcoverage();
  testuvmatrix();
}

void textOutput()
{
  int i,j;
  ostream *out;
  vector<double> point;
  if (niter>0 && (ndims>0 || primelist.size()))
  {
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    if (filename.length())
      out=new ofstream(filename);
    else
      out=&cout;
    for (i=0;i<niter;i++)
    {
      point=quads[0].dgen();
      for (j=0;j<point.size();j++)
      {
	if (j)
	  *out<<' ';
	*out<<ldecimal(point[j]);
      }
      *out<<endl;
    }
    if (filename.length())
      delete out;
  }
  if (niter<=0)
    cerr<<"Please specify number of iterations with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void computeDiscrepancy()
{
  int i;
  vector<vector<double> > points;
  if (niter>1 && (ndims>0 || primelist.size()))
  {
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    for (i=0;i<niter;i++)
      points.push_back(quads[0].dgen());
    cout<<ldecimal(discrepancy(points))<<endl;
  }
  if (niter<2)
    cerr<<"Please specify number of points (at least 2) with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void plotDiscrepancy()
{
  int i;
  set<int> halfsteps=hsteps(1,niter);
  set<int>::iterator it;
  vector<vector<double> > points;
  vector<int> npts;
  double lastdisc;
  vector<double> disc;
  if (niter>1 && (ndims>0 || primelist.size()))
  {
    ps.open(filename.length()?filename:"discrepancy.ps");
    quads[0].init(ndims,resolution);
    quads[0].init(primelist,resolution);
    quads[0].setscramble(scramble);
    for (i=0;i<niter;i++)
    {
      points.push_back(quads[0].dgen());
      if (halfsteps.count(i+1))
      {
	npts.push_back(i+1);
	cout<<i+1<<' ';
	lastdisc=discrepancy(points,true);
	disc.push_back(log(lastdisc));
	cout<<'\r';
      }
    }
    cout<<endl;
    ps.setpaper(a4land,0);
    ps.prolog();
    logLogPlot(ps,halfsteps,disc);
  }
  if (niter<2)
    cerr<<"Please specify number of points (at least 2) with -n\n";
  if (ndims<=0 && primelist.size()==0)
    cerr<<"Please specify number of dimensions with -d or primes with -p\n";
}

void writeshort(ostream &file,unsigned short i)
{
  char buf[2];
  *(unsigned short *)buf=i;
  file.write(buf,2);
}

int facRev(int n,const vector<int> &factors, map<int,vector<int> > scram)
{
  int i;
  vector<int> parts;
  for (i=0;i<factors.size();i++)
  {
    parts.push_back(n%factors[i]);
    n/=factors[i];
  }
  for (i=0;i<factors.size();i++)
  {
    n=n*factors[i]+scram[factors[i]][parts[i]];
  }
  return n;
}

vector<int> stepStairs(int p,int s)
{
  int i,q,diff;
  vector<int> ret;
  set<int> sret;
  set<int>::iterator j;
  q=lrint((double)p/s);
  for (i=0;i<=s;i++)
    ret.push_back(i*q);
  diff=p-1-ret[s];
  for (i=0;i<diff;i++)
    ret[s-i]+=diff-i;
  for (i=0;i<-diff;i++)
    ret[s-i]+=i+diff;
  for (i=0;i<=s;i++)
    sret.insert(scrambledig(ret[i],p));
  ret.clear();
  for (j=sret.begin();j!=sret.end();++j)
    ret.push_back(*j);
  return ret;
}

void writeStairs(ostream &file,vector<int> stairs)
{
  int i;
  for (i=1;i<stairs.size();i++)
  {
    file.put(stairs[i]-stairs[i-1]);
    if (stairs[i]-stairs[i-1]>maxStairStep)
      maxStairStep=stairs[i]-stairs[i-1];
    if (stairs[i]-stairs[i-1]>255)
      cout<<stairs[i]-stairs[i-1]<<" too big\n";
  }
}

void writePerm(ostream &file,vector<int> perm)
{
  int i,last;
  for (i=last=0;i<perm.size();i++)
  {
    writeshort(file,perm[i]-last);
    last=perm[i];
  }
}

void newScramble()
/* This is a sequence of scrambling permutations defined recursively
 * as follows:
 * P(2)=(0,1)
 * P(3)=(0,1,2)
 * P(p) where p>4 is (0,a-S(a),p-1) interleaved with b-1-S(b)
 * where S(n)[i] is computed by factoring n, expressing i as if the factors
 * were a base, reversing the expression, applying the previously computed
 * P(q) where q are the factors to the "digits", and turning it back
 * into one number. a is 2 less than the index of p (p=5,a=0; p=7,a=1; p=11,a=2;
 * etc.), a+b+2=p, and the interleaving uses the power scrambling function.
 */
{
  map<int,vector<int> > scram,submap;
  vector<int> row,perm0,perm1;
  int p,i,j,last,sz,stsz,inx;
  vector<int> factors,stairs,skipStairs;
  ofstream permuteFile("permute.dat",ios::binary);
  ofstream permuteText("permute.txt");
  PostScript ps;
  scram[2].push_back(0);
  scram[2].push_back(1);
  scram[3].push_back(0);
  scram[3].push_back(1);
  scram[3].push_back(2);
  maxStairStep=0;
  for (p=5,inx=0;p<65536;p+=2)
  {
    factors=factor(p);
    if (factors.size()>1)
      continue;
    cout<<p<<"     \r";
    cout.flush();
    inx++;
    stairs=stepStairs(p,inx);
    writeStairs(permuteFile,stairs);
    skipStairs.clear();
    for (i=j=0;i<p;i++)
      if (i==stairs[j])
	j++;
      else
	skipStairs.push_back(i);
    stsz=skipStairs.size();
    factors=factor(stsz);
    sz=factors.size();
    for (i=0;i*2<sz;i++)
      swap(factors[i],factors[sz-1-i]);
    if (p>66000)
      cout<<p<<' '<<factors[0]<<endl;
    row.resize(p);
    perm0.clear();
    perm1.clear();
    submap.clear();
    for (i=0;i<sz;i++)
      submap[factors[i]]=scram[factors[i]];
    for (i=0;i<stsz;i++)
    {
      perm0.push_back(stsz-1-facRev(i,factors,submap));
      row[skipStairs[i]]=skipStairs[perm0.back()];
    }
    writePerm(permuteFile,perm0);
    stsz=stairs.size()-2;
    factors=factor(stsz);
    sz=factors.size();
    for (i=0;i*2<sz;i++)
      swap(factors[i],factors[sz-1-i]);
    submap.clear();
    for (i=0;i<sz;i++)
      submap[factors[i]]=scram[factors[i]];
    for (i=0;i<stsz;i++)
    {
      perm1.push_back(stsz-1-facRev(i,factors,submap));
      row[stairs[i+1]]=stairs[perm1.back()+1];
    }
    row[0]=0;
    row[p-1]=p-1;
    writePerm(permuteFile,perm1);
    if (p<128)
    {
      for (i=0;i<p;i++)
	permuteText<<(i?' ':'(')<<row[i];
      permuteText<<")\n";
    }
    last=0;
    scram[p]=row;
  }
  ps.open("scramble.ps");
  ps.setpaper(a4land,0);
  ps.prolog();
  ps.startpage();
  ps.setscale(0,0,65026,65026);
  for (i=0;i<65027;i++)
    ps.dot(i,scram[65027][i]);
  ps.endpage();
  ps.startpage();
  ps.setscale(0,0,65028,65028);
  for (i=0;i<65029;i++)
    ps.dot(i,scram[65029][i]);
  ps.endpage();
  ps.close();
  cout<<"Maximum stair step: "<<maxStairStep<<'\n';
}

quadirr findMinMaxQuad(int p)
{
  int i,j,k,l;
  quadirr q,ret;
  QuadMax eqc;
  int minmaxterm=100000;
  ContinuedFraction cf;
  int triesSince=0;
  bool cont=true;
  for (i=1;cont;i=i+k+1)
    for (k=0;cont && i>0;i--,k++)
      for (j=1;j<=i;j++)
	if (gcd(i,j)==1)
	  for (l=0;l<k;l++)
	    if (gcd(k,l)==1)
	    {
	      q=quadirr(l,k,j,i,p);
	      try
	      {
		eqc=equivClass(q); // can throw OVERFLOW
		triesSince++;
		if (eqc.max<minmaxterm)
		{
		  minmaxterm=eqc.max;
		  ret=q;
		  cf=contFrac(eqc.qi);
		  triesSince=0;
		}
	      }
	      catch (...)
	      {
	      }
	      if (triesSince>256) // 4096
		cont=false;
	    }
  /*cout<<ret.stringval()<<' ';
  for (i=0;i<cf.terms.size();i++)
  {
    if (i==cf.terms.size()-cf.period)
      cout<<'(';
    else
      cout<<' ';
    cout<<cf.terms[i];
  }
  if (cf.period)
    cout<<')';
  cout<<endl;*/
  return ret;
}

void checkEquivClasses(int p)
{
  int i,j,k,l;
  map<double,quadirr> classes;
  map<double,vector<quadirr> > originals;
  map<double,quadirr>::iterator it;
  quadirr q;
  QuadMax eqc;
  double rv;
  int minmaxterm=100000;
  ContinuedFraction cf;
  for (i=1;i<=p;i++)
  {
    cout<<i<<'\r';
    cout.flush();
    for (j=1;j<=i;j++)
      if (gcd(i,j)==1)
	for (k=1;k<=p;k++)
	  for (l=0;l<k;l++)
	    if (gcd(k,l)==1)
	    {
	      q=quadirr(l,k,j,i,p);
	      if (q.realval())
	      {
		try
		{
		  eqc=equivClass(q); // can throw OVERFLOW
		  classes[eqc.qi.realval()]=eqc.qi;
		  originals[eqc.qi.realval()].push_back(q);
		}
		catch (...)
		{
		}
	      }
	    }
  }
  //cout<<classes.size()<<" equivalence classes:\n";
  for (it=classes.begin();it!=classes.end();++it)
  {
    cf=contFrac(q=it->second);
    if (cf.maximumTerm()<minmaxterm)
      minmaxterm=cf.maximumTerm();
  }
  for (it=classes.begin();it!=classes.end();++it)
  {
    cf=contFrac(q=it->second);
    if (cf.maximumTerm()==minmaxterm)
    {
      cout<<ldecimal(rv=it->first)<<" = "<<q.stringval()<<endl;
      for (i=0;i<cf.terms.size();i++)
      {
	if (i==cf.terms.size()-cf.period)
	  cout<<'(';
	else
	  cout<<' ';
	cout<<cf.terms[i];
      }
      if (cf.period)
	cout<<')';
      cout<<endl;
      for (i=0;i<originals[rv].size();i++)
      {
	if (i)
	  cout<<',';
	cout<<originals[rv][i].stringval();
      }
      cout<<endl;
    }
  }
}

void sortPrimes()
{
  int i,j,n;
  int min[512],max[512];
  const double logKhinchin=log(2.685452001065306);
  quads[0].init(0,0);
  ofstream primeFile("primes.dat",ios::binary);
  ofstream primeText("primes.txt");
  ofstream bFile("b322289.txt");
  //ofstream minMaxFile("minmax.txt");
  histogram hist(0,6.235); // log scale, 1 to 510
  PostScript ps;
  ContinuedFraction cf;
  map<int,quadirr> qi;
  set<PrimeContinuedFraction> pcf;
  set<PrimeContinuedFraction>::iterator k;
  PrimeContinuedFraction pcf0;
  for (i=0;i<512;i++)
  {
    min[i]=65536;
    max[i]=0;
  }
  for (i=0;i<QL_MAX_DIMS;i++)
  {
    cout<<nthprime(i)<<"     \r";
    cout.flush();
    if (useMinMax)
      qi[nthprime(i)]=findMinMaxQuad(nthprime(i));
    else
      qi[nthprime(i)]=nthquadQi(i);
    cf=contFrac(qi[nthprime(i)]);
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
  for (n=1,k=pcf.begin();k!=pcf.end();++n,++k)
  {
    writeshort(primeFile,k->prime);
    writeshort(primeFile,k->cf.terms.size());
    writeshort(primeFile,k->cf.period);
    for (i=0;i<k->cf.terms.size();i++)
      writeshort(primeFile,k->cf.terms[i]);
    if (n<=4228)
      bFile<<n<<' '<<k->prime<<endl;
    primeText<<k->prime<<' '<<qi[k->prime].stringval()<<'='<<ldecimal(qi[k->prime].realval())<<' '<<ldecimal(k->cf.averageTerm())<<' '<<k->cf.maximumTerm();
    for (i=j=0;i<40;i++,j++)
    {
      if (j>=k->cf.terms.size())
	j-=k->cf.period;
      primeText<<(i?' ':'\n')<<k->cf.terms[j];
    }
    primeText<<" ...\n";
    if (k->prime>max[k->cf.maximumTerm()])
      max[k->cf.maximumTerm()]=k->prime;
    if (k->prime<min[k->cf.maximumTerm()])
      min[k->cf.maximumTerm()]=k->prime;
  }
  //for (i=1;i<511;i++)
    //minMaxFile<<i<<' '<<min[i]<<' '<<max[i]<<'\n';
  ps.open("primes.ps");
  ps.setpaper(a4land,0);
  ps.prolog();
  ps.startpage();
  hist.plot(ps,HISTO_LOG);
  ps.setcolor(1,0,0);
  ps.line2p(xy(logKhinchin*3/6.235,2),xy(logKhinchin*3/6.235,2.1));
}

void primesqrt2()
/* 2 2
 * 2 3
 * 5 7
 * 29 41
 * 691 977
 * 773 1093
 * 971 1373
 * 1009 1427
 * 2617 3701
 * 6079 8597
 * 15731 22247
 * 27073 38287
 * 37879 53569
 * 43789 61927
 */
{
  int i=0,j=0;
  double close=1,diff;
  while (j<=QL_MAX_DIMS)
  {
    diff=double(primes[j])/primes[i]-M_SQRT2;
    if (fabs(diff)<close)
    {
      cout<<primes[i]<<' '<<primes[j]<<endl;
      close=fabs(diff);
    }
    if (diff>0)
      ++i;
    else
      ++j;
  }
}

int main(int argc,char **argv)
/* Commands:
 * sortprimes	Write a list of primes sorted by average continued fraction term
 * scramblehalton Compute permutations for scrambling Halton
 * test		Test the program
 * longtest	Run long tests
 * quadplot	Plot quadratic irrationals mod 1
 * scatter	Draw scatterplots of pairs of dimensions of a sequence
 * circle	Test how well pairs of dimensions estimate the area of a circle
 * fill		Graph how well a sequence fills space
 * flower	Draw flower plots
 * fourier	Plot Fourier transform of each dimension
 * discrepancy	Generate points and compute discrepancy
 * discplot	Plot discrepancy versus number of points
 * textout	Output a text file for the star_discrepancy program
 * interact	Run interactively
 * Options:
 * -d n		Use the first n primes in sortprimes order (d means dimensions)
 * -p p1,p2,p3	Use the specified primes
 * -r n		Set the resolution to n (default 1e17)
 * -s x		Set the scrambling option (none, third, morse, gray)
 * -n n		Output n lines (textout) or run n iterations (testprimes)
 * -o fname	Write to the specified file
 */
{
  int cmd,i;
  string cmdstr;
  bool validArgs,validCmd=true;
  int nthreads;
  po::options_description generic("Options");
  po::options_description hidden("Hidden options");
  po::options_description cmdline_options;
  po::positional_options_description p;
  po::variables_map vm;
  generic.add_options()
    ("dimensions,d",po::value<int>(&ndims),"Number of dimensions")
    ("primes,p",po::value<string>(&primestr),"List of primes")
    ("resolution,r",po::value<string>(&resstr)->default_value("1e17"),"Resolution (H for Halton)")
    ("scramble,s",po::value<string>(&scramblestr)->default_value("default"),"Scrambling: none, third, Thue-Morse, Gray, power, Faure, tipwitch, default")
    ("niter,n",po::value<int>(&niter),"Number of iterations or lines of output")
    ("threads,t",po::value<int>(&nthreads)->default_value(thread::hardware_concurrency()),"Number of threads")
    ("disc","Compute discrepancy of plot")
    ("output,o",po::value<string>(&filename),"Output file");
  hidden.add_options()
    ("command",po::value<string>(&cmdstr),"Command");
  p.add("command",1);
  cmdline_options.add(generic).add(hidden);
  commands.push_back(command("sortprimes",sortPrimes,"Sort primes by average continued fraction term"));
  commands.push_back(command("scramblehalton",newScramble,"Compute permutations for scrambling Halton"));
  commands.push_back(command("test",runTests,"Run unit and other short tests"));
  commands.push_back(command("longtest",runLongTests,"Run long tests"));
  commands.push_back(command("quadplot",quadPlot,"Plot quadratic irrationals mod 1"));
  commands.push_back(command("scatter",testScatter,"Scatter plot pairs of primes"));
  commands.push_back(command("circle",testCircle,"Test pairs of primes by estimating area of circle"));
  commands.push_back(command("fill",testFill,"Graph how well sequence fills space"));
  commands.push_back(command("flower",testFlower,"Graph each dimension separately"));
  commands.push_back(command("fourier",testFourier,"Plot log spectrum of each dimension"));
  commands.push_back(command("discrepancy",computeDiscrepancy,"Generate points and compute discrepancy"));
  commands.push_back(command("discplot",plotDiscrepancy,"Plot discrepancy versus number of points"));
  commands.push_back(command("textout",textOutput,"Output a text stream of points"));
  commands.push_back(command("interact",interact,"Enter interactive mode"));
  try
  {
    po::store(po::command_line_parser(argc,argv).options(cmdline_options).positional(p).run(),vm);
    po::notify(vm);
    scramble=parseScramble(scramblestr);
    if (scramble<0)
      cerr<<"Unrecognized or ambiguous scrambling: "<<scramblestr<<endl;
    resolution=parseResolution(resstr);
    if (!(resolution>=0))
      cerr<<"Invalid resolution (should be positive number or H): "<<resstr<<endl;
    disc2d=vm.count("disc")>0;
    validArgs=parsePrimeList() && scramble>=0 && resolution>=0;
  }
  catch (exception &e)
  {
    cerr<<e.what()<<endl;
    validCmd=false;
  }
  for (cmd=-1,i=0;i<commands.size();i++)
    if (commands[i].word==cmdstr)
      cmd=i;
  if (nthreads<0)
    nthreads=1;
  /* Fuzzing should be done with 0 threads, but calculating discrepancy
   * hangs if done with 0 threads.
   */
  if (cmd==9 && nthreads==0) // discrepancy
    nthreads=thread::hardware_concurrency();
  startThreads(nthreads);
  waitForThreads(TH_RUN);
  switch (nthprime(0)) // This initializes the list of primes.
  {
    case 2:
      cout<<"primes.dat has not been installed\n";
      break;
    case 5:
      break;
    default:
      cout<<"primes.dat is non-standard\n";
  }
  if (cmd>=0)
    if (validArgs)
      commands[cmd].fun();
    else;
  else
  {
    if (cmdstr.length())
    {
      cerr<<"Unrecognized command: "<<cmdstr<<endl;
      validCmd=false;
    }
    listCommands();
    cout<<generic<<endl;
  }
  for (i=0;i<0;i++)
    findMinMaxQuad(primes[i]);
  waitForThreads(TH_STOP);
  joinThreads();
  destroyPlans();
  if (testfail)
    cerr<<"Test failed\n";
  return !validArgs || !validCmd || testfail;
}

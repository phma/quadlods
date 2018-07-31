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
#include <boost/program_options.hpp>
#include "quadlods.h"
#include "ps.h"
#include "discrepancy.h"
#include "circletest.h"
#include "filltest.h"
#include "contfrac.h"
#include "ldecimal.h"
#include "histogram.h"

using namespace std;
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
double resolution;
string primestr;
vector<int> primelist;
string jumblestr;
int ndims,niter,jumble;
string filename;

void listCommands()
{
  int i;
  cout<<"Commands:\n";
  for (i=0;i<commands.size();i++)
    cout<<setw(15)<<left<<commands[i].word<<commands[i].desc<<'\n';
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
	num=stoi(numstr);
	if (isValidPrime(num))
	  primelist.push_back(num);
	else
	{
	  cerr<<num<<" is not a prime less than 65535\n";
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
  for (i=0;i<word.length()-1;i++)
  {
    l0=tolower(word[i]);
    l1=tolower(word[i+1]);
    if (l0>='a' && l0<='z' && l1>='a' && l1<='z')
      ret[(l0-'a')*26+(l1-'a')]++;
  }
  return ret;
}

bool parseJumble()
{
  array<short,676> dig0=digraphs("None");
  array<short,676> dig1=digraphs("Third");
  array<short,676> dig2=digraphs("Thue-Morse");
  array<short,676> dig3=digraphs("Gray");
  array<short,676> digj=digraphs(jumblestr);
  int i,match0=0,match1=0,match2=0,match3=0,maxmatch,nmatch=0;
  for (i=0;i<676;i++)
  {
    match0+=dig0[i]*digj[i];
    match1+=dig1[i]*digj[i];
    match2+=dig2[i]*digj[i];
    match3+=dig3[i]*digj[i];
  }
  maxmatch=match0;
  if (match1>maxmatch)
    maxmatch=match1;
  if (match2>maxmatch)
    maxmatch=match2;
  if (match3>maxmatch)
    maxmatch=match3;
  if (match0==maxmatch)
  {
    nmatch++;
    jumble=QL_JUMBLE_NONE;
  }
  if (match1==maxmatch)
  {
    nmatch++;
    jumble=QL_JUMBLE_THIRD;
  }
  if (match2==maxmatch)
  {
    nmatch++;
    jumble=QL_JUMBLE_THUEMORSE;
  }
  if (match3==maxmatch)
  {
    nmatch++;
    jumble=QL_JUMBLE_GRAY;
  }
  if (nmatch>1)
  {
    cerr<<"Unrecognized or ambiguous jumbling: "<<jumblestr<<endl;
    jumble=-1;
  }
  return jumble>=0;
}

void plotxy(quadlods& quad,int xdim,int ydim)
{
  int i;
  double x,y;
  vector<double> point;
  ps.startpage();
  ps.setscale(0,0,1,1);
  // (2,0) and (3,0) look splotchy at 30000, but fill in well at 100000.
  for (i=0;i<niter;i++)
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

/* Some bad pairs of primes:
 * 65027,65029
 * 59053,59051
 * 13691,13693
 */

void testCircle()
{
  int i,j;
  ps.open(filename.length()?filename:"circletest.ps");
  quads.init(ndims,resolution);
  quads.init(primelist,resolution);
  quads.setjumble(jumble);
  quads.advance(-1);
  circletest(quads,niter,ps);
  ps.trailer();
  ps.close();
}

void testScatter()
{
  int i,j;
  ps.open(filename.length()?filename:"scattertest.ps");
  ps.prolog();
  quads.init(ndims,resolution);
  quads.init(primelist,resolution);
  quads.setjumble(jumble);
  quads.advance(-1);
  for (i=0;i<quads.size();i++)
    cout<<quads.getnum(i)<<'/'<<quads.getdenom(i)<<' '<<quads.getacc(i)<<endl;
  for (i=0;i<quads.size();i++)
    for (j=0;j<i;j++)
      plotxy(quads,i,j);
  ps.trailer();
  ps.close();
}

void testFill()
{
  int i,j;
  ps.open(filename.length()?filename:"filltest.ps");
  quads.init(ndims,resolution);
  quads.init(primelist,resolution);
  quads.setjumble(jumble);
  quads.advance(-1);
  filltest(quads,niter,ps);
  ps.trailer();
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
  ps.setpaper(a4land,0);
  ps.prolog();
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
/* Commands:
 * sortprimes	Write a list of primes sorted by average continued fraction term
 * coverage	Test generator with small numbers 2, 3, 5
 * scatter	Draw scatterplots of pairs of dimensions of a sequence
 * circle	Test how well pairs of dimensions estimate the area of a circle
 * fill		Graph how well a sequence fills space
 * textout	Output a text file for the star_discrepancy program
 * interact	Run interactively
 * Options:
 * -d n		Use the first n primes in sortprimes order (d means dimensions)
 * -p p1,p2,p3	Use the specified primes
 * -r n		Set the resolution to n (default 1e17)
 * -j x		Set the jumbling option (none, third, morse, gray)
 * -n n		Output n lines (textout) or run n iterations (testprimes) (default 1048576)
 * -o fname	Write to the specified file
 */
{
  int cmd,i;
  string cmdstr;
  bool validArgs,validCmd=true;
  po::options_description generic("Options");
  po::options_description hidden("Hidden options");
  po::options_description cmdline_options;
  po::positional_options_description p;
  po::variables_map vm;
  generic.add_options()
    ("dimensions,d",po::value<int>(&ndims),"Number of dimensions")
    ("primes,p",po::value<string>(&primestr),"List of primes")
    ("resolution,r",po::value<double>(&resolution)->default_value(1e17),"Resolution")
    ("jumble,j",po::value<string>(&jumblestr)->default_value("Gray"),"Jumbling: none, third, Thue-Morse, Gray")
    ("niter,n",po::value<int>(&niter),"Number of iterations or lines of output")
    ("output,o",po::value<string>(&filename),"Output file");
  hidden.add_options()
    ("command",po::value<string>(&cmdstr),"Command");
  p.add("command",1);
  cmdline_options.add(generic).add(hidden);
  commands.push_back(command("sortprimes",sortPrimes,"Sort primes by average continued fraction term"));
  commands.push_back(command("coverage",testcoverage,"Test coverage of small 3D generator"));
  commands.push_back(command("scatter",testScatter,"Scatter plot pairs of primes"));
  commands.push_back(command("circle",testCircle,"Test pairs of primes by estimating area of circle"));
  commands.push_back(command("fill",testFill,"Graph how well sequence fills space"));
  po::store(po::command_line_parser(argc,argv).options(cmdline_options).positional(p).run(),vm);
  po::notify(vm);
  validArgs=parsePrimeList()&parseJumble();
  for (cmd=-1,i=0;i<commands.size();i++)
    if (commands[i].word==cmdstr)
      cmd=i;
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
  return !validArgs || !validCmd;
}

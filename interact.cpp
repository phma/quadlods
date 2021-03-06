/******************************************************/
/*                                                    */
/* interact.cpp - interactive mode                    */
/*                                                    */
/******************************************************/
/* Copyright 2019,2020 Pierre Abbat.
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
/* On 2020-10-01 I fuzzed interact for 16 hours; afl++ found 23 hangs
 * and no crashes. All the hangs fell in one of three classes:
 * • Initialize a generator, then generate a huge number of tuples.
 * • Initialize a generator with a huge number of dimensions, then generate.
 * • Compute a continued fraction with a huge period.
 * Examples:
 * init 0 5 1e17 gray
 * gene 0 2010709
 * init 1 52222222 1e17 gray
 * gene 1 20
 * cfra 0 1 1 1 622222222 20
 * These validly take a long time; they are not bugs.
 */
#include <iostream>
#include <cassert>
#include <boost/locale.hpp>
#include "config.h"
#include "interact.h"
#include "ldecimal.h"
#include "main.h"
#include "random.h"
#include "contfrac.h"
using namespace std;
using namespace boost::locale;

map<int,int> formats;

int commandInt(string &command)
/* Removes the first word of command and returns it as an int.
 * The commands are four letters each.
 */
{
  int i,ret=0,ch;
  for (i=0;i<command.length();i++)
  {
    ch=command[i];
    if (isspace(ch))
      break;
    ch=toupper(ch);
    if (i<4)
      ret+=(ch&0xff)<<(8*(3-i));
  }
  command.erase(0,i+1);
  return ret;
}

string firstWord(string &command)
{
  int i,ch;
  string ret;
  for (i=0;i<command.length();i++)
  {
    ch=command[i];
    if (isspace(ch))
      break;
    ret+=(char)ch;
  }
  command.erase(0,i+1);
  return ret;
}

void reply(int code,bool done,string text)
{
  assert(code>=100 && code<1000);
  cout<<code<<(done?' ':'-')<<text<<endl;
}

void cmdInit(string command)
{
  int n,s,scram;
  double res;
  string scramStr;
  int replyCode=200;
  string replyText=boost::locale::gettext("OK");
  try
  {
    n=parseInt(firstWord(command));
    s=parseInt(firstWord(command));
    res=parseResolution(firstWord(command));
  }
  catch (...)
  {
    replyCode=420;
    replyText=boost::locale::gettext("Parse error");
  }
  scramStr=firstWord(command);
  if (s<0)
  {
    replyCode=400;
    replyText=boost::locale::gettext("Number of dimensions must be nonnegative");
  }
  if (res<0) 
  {
    replyCode=401;
    replyText=boost::locale::gettext("Resolution must be positive or Halton");
  }
  if (replyCode<300)
  {
    scram=parseScramble(scramStr);
    if (scram<0)
    {
      replyCode=402;
      replyText=boost::locale::gettext("Unrecognized scrambling method");
    }
  }
  if (replyCode<300)
  {
    try
    {
      quads[n].init(s,res);
      if (formats[n]==0)
	formats[n]=10;
      quads[n].setscramble(scram);
    }
    catch (...)
    {
      replyCode=500;
      replyText=boost::locale::gettext("Internal service error");
    }
  }
  reply(replyCode,true,replyText);
}

string toString(vector<mpq_class> tuple,int format)
{
  string ret;
  double dval;
  int i;
  char buf[32];
  for (i=0;i<tuple.size();i++)
  {
    if (i)
      ret+=' ';
    if (format&0xff00)
      ret+=tuple[i].get_str(format&0xff);
    else
    {
      dval=tuple[i].get_d();
      if (format==10)
	ret+=ldecimal(dval);
      else
      {
	snprintf(buf,31,"%a",dval);
	ret+=buf;
      }
    }
  }
  return ret;
}

void cmdGene(string command)
{
  int n,i,j;
  int replyCode=200;
  string replyText=boost::locale::gettext("OK");
  vector<mpq_class> tuple;
  try
  {
    n=parseInt(firstWord(command));
    i=parseInt(firstWord(command));
  }
  catch (...)
  {
    replyCode=420;
    replyText=boost::locale::gettext("Parse error");
  }
  if (replyCode<300 && quads.count(n)==0)
  {
    replyCode=410;
    replyText=boost::locale::gettext("Generator is uninitialized");
  }
  if (replyCode<300)
  {
    if (i>0)
      replyCode=220;
    if (i<0)
    {
      replyCode=402;
      replyText=boost::locale::gettext("Number of tuples must be positive");
    }
  }
  if (replyCode==220)
    for (j=0;j<i;j++)
    {
      tuple=quads[n].gen();
      replyText=toString(tuple,formats[n]);
      reply(replyCode,!(i-j-1),replyText);
    }
  else
    reply(replyCode,true,replyText);
}

size_t findsubseq(string haystack,string needle)
{
  size_t i,j,ret=ULLONG_MAX;
  for (i=j=0;i<haystack.length() && j<needle.length();i++)
    if (tolower(haystack[i])==tolower(needle[j]))
      if (++j==needle.length())
	ret=i;
  return ret;
}

void cmdSeed(string command)
{
  int n,i,b;
  int replyCode=200;
  string replyText=boost::locale::gettext("OK");
  vector<char> bytes;
  try
  {
    n=parseInt(firstWord(command));
  }
  catch (...)
  {
    replyCode=420;
    replyText=boost::locale::gettext("Parse error");
  }
  if (replyCode<300 && quads.count(n)==0)
  {
    replyCode=410;
    replyText=boost::locale::gettext("Generator is uninitialized");
  }
  if (replyCode<300)
  {
    b=quads[n].seedsize();
    for (i=0;i<b;i++)
      bytes.push_back(rng.ucrandom());
    quads[n].seed(&bytes[0],b);
  }
  reply(replyCode,true,replyText);
}

int parseFormat(string fmt)
{
  size_t pos,lastpos=ULLONG_MAX;
  int ret=-1;
  pos=findsubseq(fmt,"dec");
  if (pos<lastpos)
  {
    lastpos=pos;
    ret=10;
  }
  pos=findsubseq(fmt,"hex");
  if (pos<lastpos)
  {
    lastpos=pos;
    ret=16;
  }
  pos=findsubseq(fmt,"flo");
  if (pos<lastpos)
  {
    lastpos=pos;
    ret=0;
  }
  pos=findsubseq(fmt,"rat");
  if (pos<lastpos)
  {
    lastpos=pos;
    ret=256;
  }
  return ret;
}

void cmdForm(string command)
{
  int n,fmt,fmt1,replyCode=200;
  string replyText=boost::locale::gettext("OK");
  try
  {
    n=parseInt(firstWord(command));
  }
  catch (...)
  {
    replyCode=420;
    replyText=boost::locale::gettext("Parse error");
  }
  if (replyCode<300 && quads.count(n)==0)
  {
    replyCode=410;
    replyText=boost::locale::gettext("Generator is uninitialized");
  }
  if (replyCode<300)
    fmt=formats[n];
  while (replyCode<300 && command.length())
  {
    fmt1=parseFormat(firstWord(command));
    if (fmt1<0)
    {
      replyCode=421;
      replyText=boost::locale::gettext("Unknown format");
    }
    if (fmt1&0xff)
      fmt=(fmt&0xff00)+fmt1;
    else
      fmt=(fmt&0xff)+fmt1;
  }
  if (replyCode<300)
    formats[n]=fmt;
  reply(replyCode,true,replyText);
}

void cmdCfra(string command)
{
  int replyCode=220;
  int a,b,c,d,p,i;
  quadirr q;
  string replyText;
  ContinuedFraction cf;
  try
  {
    a=parseInt(firstWord(command));
    b=parseInt(firstWord(command));
    c=parseInt(firstWord(command));
    d=parseInt(firstWord(command));
    p=parseInt(firstWord(command));
  }
  catch (...)
  {
    replyCode=420;
    replyText=boost::locale::gettext("Parse error");
  }
  if (replyCode<300)
  {
    try
    {
      q=quadirr(a,b,c,d,p);
      cf=contFrac(q);
    }
    catch (int e)
    {
      replyCode=410;
      replyText="Error "+to_string(e);
      replyText=boost::locale::gettext(replyText.c_str());
#if 0
      gettext("Error 1"); // overflow
      gettext("Error 2"); // zero divide
      gettext("Error 3"); // imaginary
#endif
    }
  }
  if (replyCode<300)
  {
    replyText=q.stringval()+"=[";
    for (i=0;i<cf.terms.size();i++)
    {
      if (i)
	replyText+=(i>1)?',':';';
      if (i==cf.terms.size()-cf.period)
	replyText+='(';
      replyText+=to_string(cf.terms[i]);
    }
    if (cf.period)
      replyText+=')';
    replyText+=']';
  }
  reply(replyCode,true,replyText);
}

void cmdHelp(string command)
{
  int replyCode=220;
  string replyText=boost::locale::gettext("help");
  string line;
  size_t pos;
  while ((pos=replyText.find('\n'))<=replyText.length())
  {
    line=replyText.substr(0,pos);
    reply(replyCode,false,line);
    replyText.erase(0,pos+1);
  }
  reply(replyCode,true,replyText);
}

/* Commands for interactive mode, which can be used as a server:
 * init n s res scram: Initialize generator #n with s dimensions and resolution res.
 * form n dec/hex/flo/rat: Set format to decimal/hexadecimal/floating point/rational.
 * gene n i: Generate i points from generator n.
 * seed n: Seed generator n with random numbers.
 */
void interact()
{
  bool cont=true;
  string command;
  int opcode;
  generator gen;
  gen.add_messages_path(SHARE_DIR);
  gen.add_messages_domain("interact");
  locale::global(gen("en_US.UTF-8")); // UTF-8 because of √ sign
  reply(220,true,string("Quadlods version ")+VERSION+" ready");
  while (cont)
  {
    getline(cin,command);
    if (!cin)
    {
      cont=false;
      continue;
    }
    opcode=commandInt(command);
    switch (opcode)
    {
      case 0x45584954:
      case 0x51554954:
	cont=false;
	reply(221,true,"Quadlods exiting");
	break;
      case 0x494e4954:
	cmdInit(command);
	break;
      case 0x47454e45:
	cmdGene(command);
	break;
      case 0x464f524d:
	cmdForm(command);
	break;
      case 0x53454544:
	cmdSeed(command);
	break;
      case 0x43465241:
	cmdCfra(command);
	break;
      case 0x48454c50:
	cmdHelp(command);
	break;
      default:
	reply(400,true,"Invalid command");
    }
  }
}

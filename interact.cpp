/******************************************************/
/*                                                    */
/* interact.cpp - interactive mode                    */
/*                                                    */
/******************************************************/
/* Copyright 2019 Pierre Abbat.
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
#include <cassert>
#include "config.h"
#include "interact.h"
#include "main.h"
using namespace std;

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
  int n,s;
  double res;
  string scramStr;
  int replyCode=200;
  string replyText="OK";
  n=stoi(firstWord(command));
  s=stoi(firstWord(command));
  res=stod(firstWord(command));
  scramStr=firstWord(command);
  if (s<0)
  {
    replyCode=400;
    replyText="Number of dimensions must be nonnegative";
  }
  if (res<=0) // res==0 means Halton, but that's not finished yet
  {
    replyCode=401;
    replyText="Resolution must be positive";
  }
  if (replyCode<300)
  {
    try
    {
      quads[n].init(s,res);
    }
    catch (...)
    {
      replyCode=500;
      replyText="Internal service error";
    }
  }
  reply(replyCode,true,replyText);
}

string toString(vector<mpq_class> tuple)
{
  string ret;
  int i;
  for (i=0;i<tuple.size();i++)
  {
    if (i)
      ret+=' ';
    ret+=tuple[i].get_str();
  }
  return ret;
}

void cmdGene(string command)
{
  int n,i,j;
  int replyCode=200;
  string replyText="OK";
  vector<mpq_class> tuple;
  try
  {
    n=stoi(firstWord(command));
    i=stoi(firstWord(command));
  }
  catch (...)
  {
    replyCode=420;
    replyText="Parse error";
  }
  if (replyCode<300 && quads.count(n)==0)
  {
    replyCode=410;
    replyText="Generator is uninitialized";
  }
  if (replyCode<300)
  {
    if (i>0)
      replyCode=220;
    if (i<0)
    {
      replyCode=402;
      replyText="Number of tuples must be positive";
    }
  }
  if (replyCode==220)
    for (j=0;j<i;j++)
    {
      tuple=quads[n].gen();
      replyText=toString(tuple);
      reply(replyCode,!(i-j-1),replyText);
    }
  else
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
  reply(220,true,string("Quadlods version ")+VERSION+" ready");
  while (cont)
  {
    getline(cin,command);
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
      default:
	reply(400,true,"Invalid command");
    }
  }
}

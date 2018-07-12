/******************************************************/
/*                                                    */
/* ps.cpp - PostScript output                         */
/*                                                    */
/******************************************************/
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <iomanip>
#include "ps.h"
#include "ldecimal.h"
using namespace std;

#define PAPERRES 0.004

char rscales[]={10,12,15,20,25,30,40,50,60,80};
const double PSPoint=25.4/72;
papersize a4land={297,210};

PostScript::PostScript()
{
  oldr=oldg=oldb=NAN;
  paperx=210;
  papery=297;
  scale=1;
  pageorientation=orientation=pages=0;
  indocument=inpage=inlin=false;
  psfile=nullptr;
}

PostScript::~PostScript()
{
  if (psfile)
    close();
}

void PostScript::setpaper(papersize pap,int ori)
/* ori is 0 for no rotation, 1 for 90° rotation, making portrait
 * into landscape and vice versa. Do this before each page,
 * or before calling prolog if all pages are the same size.
 */
{
  paperx=pap.width;
  papery=pap.height;
  pageorientation=ori;
}

double PostScript::aspectRatio()
// Returns >1 for landscape, <1 for portrait.
{
  if (pageorientation&1)
    return papery/paperx;
  else
    return paperx/papery;
}

void PostScript::open(string psfname)
{
  if (psfile)
    close();
  psfile=new ofstream(psfname);
}

void PostScript::prolog()
{
  if (psfile && !indocument)
  {
    *psfile<<"%!PS-Adobe-3.0\n%%BeginProlog\n%%%%Pages: (atend)"<<endl;
    *psfile<<"%%BoundingBox: 0 0 "<<rint(paperx/PSPoint)<<' '<<rint(papery/PSPoint)<<endl;
    *psfile<<"\n/. % ( x y )\n{ newpath 0.1 0 360 arc fill } bind def\n\n";
    *psfile<<"/- % ( x1 y1 x2 y2 )\n\n{ newpath moveto lineto stroke } bind def\n\n";
    *psfile<<"/c. % ( str )\n{ dup stringwidth -2 div exch -2 div exch\n"<<
            "3 2 roll 2 index 2 index rmoveto show rmoveto } bind def\n\n";
    *psfile<<"/mmscale { 720 254 div dup scale } bind def\n";
    *psfile<<"%%EndProlog"<<endl;
    indocument=true;
    pages=0;
  }
}

void PostScript::startpage()
{
  if (psfile && indocument && !inpage)
  {
    ++pages;
    *psfile<<"%%Page: "<<pages<<' '<<pages<<"\ngsave mmscale 0.1 setlinewidth\n";
    *psfile<<paperx/2<<' '<<papery/2<<' '<<"translate ";
    *psfile<<(pageorientation&3)*90<<" rotate ";
    *psfile<<paperx/-2<<' '<<papery/-2<<' '<<"translate"<<endl;
    *psfile<<"/Helvetica findfont 3 scalefont setfont"<<endl;
    oldr=oldg=oldb=NAN;
    inpage=true;
  }
}

void PostScript::endpage()
{
  if (psfile && indocument && inpage)
  {
    *psfile<<"grestore showpage"<<endl;
    inpage=false;
  }
}

void PostScript::trailer()
{
  if (inpage)
    endpage();
  if (psfile && indocument)
  {
    *psfile<<"%%BeginTrailer\n%%Pages: "<<pages<<"\n%%EndTrailer"<<endl;
    indocument=false;
  }
}

void PostScript::close()
{
  if (indocument)
    trailer();
  delete(psfile);
  psfile=nullptr;
}

double PostScript::xscale(double x)
{
  return scale*(x-centerx)+paperx/2;
}

double PostScript::yscale(double y)
{
  return scale*(y-centery)+papery/2;
}

string PostScript::escape(string text)
{
  string ret;
  int ch;
  while (text.length())
  {
    ch=text[0];
    if (ch=='(' || ch==')')
      ret+='\\';
    ret+=ch;
    text.erase(0,1);
  }
  return ret;
}

void PostScript::setcolor(double r,double g,double b)
{
  if (r!=oldr || g!=oldg || b!=oldb)
  {
    *psfile<<fixed<<setprecision(3)<<r<<' '<<g<<' '<<b<<" setrgbcolor"<<endl;
    oldr=r;
    oldg=g;
    oldb=b;
  }
}

void PostScript::setscale(double minx,double miny,double maxx,double maxy,int ori)
{
  double xsize,ysize,papx,papy;
  int i;
  orientation=ori;
  centerx=(minx+maxx)/2;
  centery=(miny+maxy)/2;
  xsize=fabs(minx-maxx);
  ysize=fabs(miny-maxy);
  papx=paperx;
  papy=papery;
  if (pageorientation&1)
    swap(papx,papy);
  for (scale=1;scale*xsize/10<papx && scale*ysize/10<papy;scale*=10);
  for (;scale*xsize/80>papx*0.9 || scale*ysize/80>papy*0.9;scale/=10);
  for (i=0;i<9 && (scale*xsize/rscales[i]>papx*0.9 || scale*ysize/rscales[i]>papy*0.9);i++);
  scale/=rscales[i];
  *psfile<<"% minx="<<minx<<" miny="<<miny<<" maxx="<<maxx<<" maxy="<<maxy<<" scale="<<scale<<endl;
}

void PostScript::dot(double x,double y)
{
  assert(psfile);
  *psfile<<fixed<<setprecision(2)<<xscale(x)<<' '<<yscale(y)<<" .";
  *psfile<<endl;
}

void PostScript::line2p(xy pnt1,xy pnt2)
{
  //pnt1=turn(pnt1,orientation);
  //pnt2=turn(pnt2,orientation);
  if (isfinite(pnt1.getx()) && isfinite(pnt1.gety()) && isfinite(pnt2.getx()) && isfinite(pnt2.gety()))
    *psfile<<ldecimal(xscale(pnt1.getx()),PAPERRES)<<' '<<ldecimal(yscale(pnt1.gety()),PAPERRES)
    <<' '<<ldecimal(xscale(pnt2.getx()),PAPERRES)<<' '<<ldecimal(yscale(pnt2.gety()),PAPERRES)<<" -"<<endl;
}

void PostScript::startline()
{
  assert(psfile);
  *psfile<<"newpath"<<endl;
}

void PostScript::lineto(double x,double y)
{
  assert(psfile);
  *psfile<<fixed<<setprecision(2)<<xscale(x)<<' '<<yscale(y)<<(inlin?" lineto":" moveto");
  *psfile<<endl;
  inlin=true;
}

void PostScript::endline(bool closed)
{
  assert(psfile);
  if (closed)
    *psfile<<"closepath ";
  *psfile<<"stroke"<<endl;
  inlin=false;
}

void PostScript::write(double x,double y,string text)
{
  *psfile<<fixed<<setprecision(2)<<xscale(x)<<' '<<yscale(y)
  <<" moveto ("<<text<<") show"<<endl;
}

void PostScript::centerWrite(xy pnt,string text)
{
  //pnt=turn(pnt,orientation);
  *psfile<<ldecimal(xscale(pnt.getx()),PAPERRES)<<' '<<ldecimal(yscale(pnt.gety()),PAPERRES)
  <<" moveto ("<<escape(text)<<") c."<<endl;
}

void PostScript::comment(string text)
{
  *psfile<<'%'<<text<<endl;
}

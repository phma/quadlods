/******************************************************/
/*                                                    */
/* ps.cpp - PostScript output                         */
/*                                                    */
/******************************************************/
/* Copyright 2014 Pierre Abbat.
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
#include <unistd.h>

FILE *psfile;
int pages;
double scale=1; // paper size is in millimeters, but model space is in meters
int orientation;
double paperx(210),papery(297),centerx,centery;
char rscales[]={10,12,15,20,25,30,40,50,60,80};

void setscale(double minx,double miny,double maxx,double maxy)
{double xsize,ysize;
 int i;
 centerx=(minx+maxx)/2;
 centery=(miny+maxy)/2;
 xsize=fabs(minx-maxx);
 ysize=fabs(miny-maxy);
 for (scale=1;scale*xsize/10<paperx && scale*ysize/10<papery;scale*=10);
 for (;scale*xsize/80>paperx*0.9 || scale*ysize/80>papery*0.9;scale/=10);
 for (i=0;i<9 && (scale*xsize/rscales[i]>paperx*0.9 || scale*ysize/rscales[i]>papery*0.9);i++);
 scale/=rscales[i];
 printf("scale=%f\n",scale);
 //sleep(3);
 }

void widen(double factor)
{fprintf(psfile,"currentlinewidth %f mul setlinewidth\n",factor);
 }

void setcolor(double r,double g,double b)
{fprintf(psfile,"%f %f %f setrgbcolor\n",r,g,b);
 }

void psopen(const char * psfname)
{psfile=fopen(psfname,"w");
 }

void psclose()
{fclose(psfile);
 printf("scale=%f\n",scale);
 sleep(3);
 }

void psprolog()
{fprintf(psfile,"%%!PS-Adobe-3.0\n\
%%%%BeginProlog\n\
%%%%Pages: (atend)\n\
%%%%BoundingBox: 0 0 596 843\n\
%% A4 paper.\n\
\n\
/. %% ( x y )\n\
{ newpath 0.3 0 360 arc fill } bind def\n\
\n\
/- %% ( x1 y1 x2 y2 )\n\
{ newpath moveto lineto stroke } bind def\n\
\n\
/mmscale { 720 254 div dup scale } bind def\n\
%%%%EndProlog\n\
");
 pages=0;
 //scale=10;
 fflush(psfile);
 }

void pstrailer()
{fprintf(psfile,"%%%%BeginTrailer\n\
%%%%Pages: %d\n\
%%%%EndTrailer\n\
",pages);
 }

double xscale(double x)
{return scale*(x-centerx)+105;
 }

double yscale(double y)
{return scale*(y-centery)+148.5;
 }

void startpage()
{++pages;
 fprintf(psfile,"%%%%Page: %d %d\ngsave mmscale 0.1 setlinewidth\n\
/Helvetica findfont 3 scalefont setfont\n",pages,pages);
 }

void endpage()
{fputs("grestore showpage\n",psfile);
 fflush(psfile);
 }

void dot(double x,double y)
{
  double r,g,b;
  fprintf(psfile,"%7.3f %7.3f .\n",
          xscale(x),yscale(y));
}

int fibmod3(int n)
{
  int i,a,b;
  for (i=a=0,b=1;a<n;i++)
  {
    b+=a;
    a=b-a;
  }
  return (a==n)?(i%3):-1;
}

void line2p(double x1,double y1,double x2,double y2)
{fprintf(psfile,"%7.3f %7.3f %7.3f %7.3f -\n",
        xscale(x1),yscale(y1),xscale(x2),yscale(y2));
 }

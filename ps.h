/******************************************************/
/*                                                    */
/* ps.h - PostScript output                           */
/*                                                    */
/******************************************************/
/* Copyright 2014,2016-2018,2023 Pierre Abbat.
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
#ifndef PS_H
#define PS_H
#include <string>
#include <iostream>
#include <fstream>
#include "xy.h"
#include "polyline.h"
#include "pairpoint.h"

struct papersize
{
  double width,height;
};

extern papersize a4land,a4port;

class PostScript
{
protected:
  std::ostream *psfile;
  int pages;
  bool indocument,inpage,inlin;
  double scale; // paper size is in millimeters
  double paperx,papery,centerx,centery;
  int orientation,pageorientation;
  double oldr,oldg,oldb;
public:
  PostScript();
  ~PostScript();
  void setpaper(papersize pap,int ori);
  double aspectRatio();
  void open(std::string psfname);
  void prolog();
  void startpage();
  void endpage();
  void trailer();
  void close();
  double xscale(double x);
  double yscale(double y);
  std::string escape(std::string text);
  void setcolor(double r,double g,double b);
  void setscale(double minx,double miny,double maxx,double maxy,int ori=0);
  void dot(double x,double y);
  void line2p(xy pnt1,xy pnt2);
  void startline();
  void lineto(double x,double y);
  void endline(bool closed=false);
  void plot(polyline pl,bool fill=false);
  void draw(PairCompressor pnts);
  void write(double x,double y,std::string text);
  void centerWrite(xy pnt,std::string text);
  void comment(std::string text);
};
#endif

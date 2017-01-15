/******************************************************/
/*                                                    */
/* ps.h - PostScript output                           */
/*                                                    */
/******************************************************/
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
#include <string>
#include <iostream>
#include <fstream>

struct papersize
{
  double width,height;
};

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
  void setcolor(double r,double g,double b);
  void setscale(double minx,double miny,double maxx,double maxy,int ori=0);
  void dot(double x,double y);
  void startline();
  void lineto(double x,double y);
  void endline();
  void comment(std::string text);
};

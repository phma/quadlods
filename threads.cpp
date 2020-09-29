/******************************************************/
/*                                                    */
/* threads.cpp - multithreading                       */
/*                                                    */
/******************************************************/
/* Copyright 2020 Pierre Abbat.
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
#include <queue>
#include <iostream>
#include "threads.h"
#include "random.h"
#include "manysum.h"
#include "discrepancy.h"
using namespace std;
namespace cr=std::chrono;

mutex startMutex;
int threadCommand;
vector<thread> threads;
vector<int> threadStatus; // Bit 8 indicates whether the thread is sleeping.
vector<double> sleepTime,sleepFraction;
int currentAction;

cr::steady_clock clk;

double busyFraction()
{
  int i,numBusy=0;
  for (i=0;i<threadStatus.size();i++)
    if ((threadStatus[i]&256)==0)
      numBusy++;
  return (double)numBusy/i;
}

void startThreads(int n)
{
  int i;
  threadCommand=TH_WAIT;
  sleepTime.resize(n);
  sleepFraction.resize(n);
  for (i=0;i<n;i++)
  {
    sleepFraction[i]=0.5;
    sleepTime[i]=1000;
    threads.push_back(thread(QuadThread(),i));
    this_thread::sleep_for(chrono::milliseconds(10));
  }
}

void joinThreads()
{
  int i;
  for (i=0;i<threads.size();i++)
    threads[i].join();
}

void sleep(int thread)
{
  sleepTime[thread]*=1.125;
  if (sleepTime[thread]>32768)
    sleepTime[thread]*=0.875;
  this_thread::sleep_for(chrono::microseconds(lrint(sleepTime[thread])));
}

void unsleep(int thread)
{
  sleepTime[thread]*=0.875;
  if (sleepTime[thread]<1)
    sleepTime[thread]*=1.125;
  if (sleepTime[thread]<0 || std::isnan(sleepTime[thread]))
    sleepTime[thread]=1;
}

double maxSleepTime()
{
  int i;
  double max=0;
  for (i=0;i<sleepTime.size();i++)
    if (sleepTime[i]>max)
      max=sleepTime[i];
  return max;
}

void setThreadCommand(int newStatus)
{
  threadCommand=newStatus;
  //cout<<statusNames[newStatus]<<endl;
}

int getThreadStatus()
/* Returns aaaaaaaaaabbbbbbbbbbcccccccccc where
 * aaaaaaaaaa is the status all threads should be in,
 * bbbbbbbbbb is 0 if all threads are in the same state, and
 * cccccccccc is the state the threads are in.
 * If all threads are in the commanded state, but some may be asleep and others awake,
 * getThreadStatus()&0x3ffbfeff is a multiple of 1048577.
 */
{
  int i,oneStatus,minStatus=-1,maxStatus=0;
  for (i=0;i<threadStatus.size();i++)
  {
    oneStatus=threadStatus[i];
    maxStatus|=oneStatus;
    minStatus&=oneStatus;
  }
  return (threadCommand<<20)|((minStatus^maxStatus)<<10)|(minStatus&0x3ff);
}

void waitForThreads(int newStatus)
// Waits until all threads are in the commanded status.
{
  int i,n;
  threadCommand=newStatus;
  do
  {
    for (i=n=0;i<threadStatus.size();i++)
      if ((threadStatus[i]&255)!=threadCommand)
	n++;
    this_thread::sleep_for(chrono::milliseconds(n));
  } while (n);
}

void QuadThread::operator()(int thread)
{
  startMutex.lock();
  if (threadStatus.size()!=thread)
  {
    cout<<"Starting thread "<<threadStatus.size()<<", was passed "<<thread<<endl;
    thread=threadStatus.size();
  }
  threadStatus.push_back(0);
  startMutex.unlock();
  while (threadCommand!=TH_STOP)
  {
    if (threadCommand==TH_RUN)
    {
      threadStatus[thread]=TH_RUN;
      if (countAnyBlock())
	unsleep(thread);
      else
	sleep(thread);
    }
    if (threadCommand==TH_PAUSE)
    { // The job is ongoing, but has to pause to write out the files.
      threadStatus[thread]=TH_PAUSE;
      sleep(thread);
    }
    if (threadCommand==TH_WAIT)
    { // There is no job. The threads are waiting for a job.
      threadStatus[thread]=TH_WAIT;
      sleep(thread);
    }
  }
  threadStatus[thread]=TH_STOP;
}

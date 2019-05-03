# -*- coding: utf-8 -*-
"""
Released on April 29, 2019

@author: Yeoh Jing Wui <bchyjw@nus.edu.sg>; Ivan Ng Kai Boon; Poh Chueh Loo <poh.chuehloo@nus.edu.sg>

The code is part of BMSS software.

Copyright (c) 2019, National University of Singapore.

"""

### This module is to get the name for the BMSS-generated files ###

import datetime

def gettxtfilename():
    timenow = datetime.datetime.now()

    # print (timenow.year, timenow.month, timenow.day, timenow.hour, timenow.minute, timenow.second)

    year = str(timenow.year % 100)
    month = str(timenow.month).zfill(2)
    day = str(timenow.day).zfill(2)
    hour = str(timenow.hour).zfill(2)
    minute = str(timenow.minute).zfill(2)
    second = str(timenow.second).zfill(2)

    Txtfilename = year+month+day+'_'+hour+minute+'.txt'
    DateTime = day + '-' + month + '-' + str(timenow.year) + ', ' + hour + ':' + minute + ':' + second

    return(Txtfilename, DateTime)

def getcsvfilename():
    timenow = datetime.datetime.now()

    # print (timenow.year, timenow.month, timenow.day, timenow.hour, timenow.minute, timenow.second)

    year = str(timenow.year % 100)
    month = str(timenow.month).zfill(2)
    day = str(timenow.day).zfill(2)
    hour = str(timenow.hour).zfill(2)
    minute = str(timenow.minute).zfill(2)
    #second = str(timenow.second).zfill(2)

    CSVfilename = year+month+day+'_'+hour+minute+'.csv'

    return (CSVfilename)

def getxmlfilename():
    timenow = datetime.datetime.now()

    # print (timenow.year, timenow.month, timenow.day, timenow.hour, timenow.minute, timenow.second)

    year = str(timenow.year % 100)
    month = str(timenow.month).zfill(2)
    day = str(timenow.day).zfill(2)
    hour = str(timenow.hour).zfill(2)
    minute = str(timenow.minute).zfill(2)
    #second = str(timenow.second).zfill(2)

    XMLfilename = year+month+day+'_'+hour+minute+'.xml'

    return (XMLfilename)



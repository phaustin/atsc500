"""
Read atmospheric soundings from a soundings file screenscraped from
http://weather.uwyo.edu/upperair/sounding.html using their "Text List" option
"""
from __future__ import print_function

from  datetime import timezone as tz
import re
from a405thermo.thermlib import esat
from a405thermo.thermlib  import constants as con
import sys
from collections import OrderedDict
import pandas as pd
import datetime

month_nums=OrderedDict([('Jan', 1), ('Feb', 2), ('Mar', 3), ('Apr', 4), ('May', 5), ('Jun', 6), ('Jul', 7),
                        ('Aug', 8), ('Sep', 9), ('Oct', 10), ('Nov', 11), ('Dec', 12)])

columns="PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV"
columns=columns.lower().split()


def parsekey(key):
    oldTup=key.split(' ')
    hour,day,month,year=oldTup
    hour=hour[:-1]  #strip the Z
    theDay=int(day)
    theHour=int(hour)
    theYear=int(year)
    theMonth=month_nums[month]
    theDate=datetime.datetime(theYear,theMonth,theDay,theHour,0,0,0,tz.utc)
    return theDate


def get_segment(filename):
    """
      72365 ABQ Albuquerque Observations at 00Z 01 Jul 2010
    is parsed to return
      ('72365 ABQ Albuquerque Observations at', '00Z', '01', 'Jul', '2010')
    """
    findchunk=re.compile('(.*?Observations at)(.*?\d{4,4})')
    with open(filename) as f:
        firstline=next(f)
        thematch=findchunk.match(firstline)
        if thematch:
            out=thematch.groups()
            time,day,month,year=out[1].split()
        else:
            raise Exception('broken')
        return out[0],time,day,month,year


def readsound(inFileName = 'soundings_july.txt'):
    """ Read San Diego DYCOMS soundings from file soundings.txt.

    Returns a dictionary allsounds[key] which contains all soundings
    in the file. The key is of the form '12Z 18 Jul 2001'.
    allsounds[key] is an array [level, variable], where level is the
    index of the level and variable is
         [,0] PRES(hPa); [,1] HGHT(m); [,2] TEMP[C];
         [,3] DWPT(C); [,4] RELH(%); [,5] MIXR(g/kg);
         [,6] DRCT(deg); [,7] SKNT(knot); [,8] THTA(K);
         [,9] THTE(K); [,10] THTV(K)
    """
    markerText,time,day,month,year=get_segment(inFileName)
    #
    # Open input file and read all lines into thelines.
    with open(inFileName) as f:
        thelines = f.read()
        thelines=thelines.strip()
    #
    # Each sounding starts with a line that contains the string
    # 'NKX San Diego Observations at'. Split the lines into individual
    # soundings.
    marker = re.compile(markerText)
    soundings = marker.split(thelines)
    #
    # Expressions that are used in the following loop.
    lineend = re.compile('\n')
    whitespace = re.compile('\s+')
    #
    # Initialise a dictionary that will hold the soundings. The soundings
    # will be accessible using date/time as key, e.g. '12Z 01 Jul 2001'.
    allsounds = {}

    #
    # Loop over all soundings: Extract the lines that contain the data and
    # convert the data from strings to floats.
    for itemcount,item in enumerate(soundings):
        newsplit = lineend.split(item)
        count = 0
        theLine = newsplit[0]
        #
        # The data lines start with a header. Search for this header.
        #
        try:
            while theLine.find('PRES   HGHT   TEMP   DWPT') < 0:
                count += 1
                theLine = newsplit[count]
        except IndexError:
            if itemcount > 0:  #item[0] is blank
                print("no column header line found in sounding" \
                      " %d\ncontents: %s" % (itemcount, item))
            continue
        count += 1
        theLine = newsplit[count]
        if theLine.find('hPa     m      C      C      %') < 0:
            raise IOError("no units line found in file %s" % theLine)
        count += 1
        theLine = newsplit[count]
        if theLine.find('-------------') < 0:
            raise IOError("expected ---- separator")
        count += 1
        theLine = newsplit[count]
        #
        # lineSave will hold all levels of data (one line in the file
        # corresponds to one level); header will become the key for the
        # sounding in allsounds, e.g. '12Z 01 Jul 2001'.
        lineSave = []
        header = newsplit[0].strip()

        while True:
            try:
                #
                # For each line in the sounding convert the string with the
                # data to an array of floats.
                # Data fields: [0] PRES(hPa); [1] HGHT(m); [2] TEMP[C];
                #              [3] DWPT(C); [4] RELH(%); [5] MIXR(g/kg);
                #              [6] DRCT(deg); [7] SKNT(knot); [8] THTA(K);
                #              [9] THTE(K); [10] THTV(K)
                #
                # at elevation, the bottom lines may be blank:
                #
                #
                #  1000.0     76                                                               
                #   925.0    777                                                               
                #   850.0   1525                                                               
                #   841.0   1620   28.8    0.8     16   4.85     85      9  317.3  333.1  318.2
                #
                dataFields = re.split(whitespace, theLine.strip())
                if len(dataFields) == 11:
                    try:
                        dataFields = [float(number) for number in dataFields]
                    except VauleError:
                        print('trouble converting dataFields to float')
                        print(dataFields)
                        sys.exit(1)
                    #
                    # Recompute the mixing ratio. This is done because due to only
                    # two digits accuracy in sounding.txt, the mixing ratio becomes
                    # 0.00 from time to time.
                    es = esat(dataFields[3] + 273.15)*0.01  #get vapor pressure from dewpoint in hPa
                    dataFields[5] = (con.eps*es/(dataFields[0] - es))*1.e3   #g/kg
                    #
                    # Only store complete lines with 11 elements.
                    lineSave.append(dataFields)
                    count += 1
                    theLine = newsplit[count]
                #
                # get the next line
                #
                count += 1
                theLine = newsplit[count]
            except IndexError:
                # out of lines, go to next sounding
                break
            
        if len(lineSave) < 5:
            #
            # The sounding should have more than 4 levels.
            print("trouble with item %d %d" % (itemcount,len(lineSave)))
        else:
            #
            # Store the sounding to allsounds as a pandas dataframe
            the_date = parsekey(header)
            allsounds[the_date] = pd.DataFrame.from_records(lineSave,columns=columns)
        the_keys=list(allsounds.keys())
        the_keys.sort()
        out_dict=OrderedDict()
        for key in the_keys:
            out_dict[key]=allsounds[key]
    return out_dict


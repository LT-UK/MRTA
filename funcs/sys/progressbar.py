#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:50:41 2020

Print progress bar in console panel. 

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""
import time

# =============================================================================
# 
# =============================================================================
class colors:
    RED = '\x1b[0;31;40m'
    GREEN = '\x1b[0;32;40m'
    ORANGE = '\x1b[0;33;40m'
    BLUE = '\x1b[0;34;40m'
    VIOLET = '\x1b[0;35;40m'
    CYAN = '\x1b[0;36;40m'
    GREY = '\x1b[0;37;40m'
    END = '\x1b[0m'
    
    
# =============================================================================
# 
# =============================================================================
# Print iterations progress
def printProgressBar(iteration, total, time_start, 
                     color = colors.CYAN,
                     prefix = 'TA Calculating:', suffix = 'Complete', 
                     decimals = 2, length = 50, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        time_start  - Required  : progress start time
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    dt = time.time() - time_start
    hours = int(dt//3600)
    mins = int(dt%3600//60)
    min_1 = int(mins//10)
    min_2 = int(mins%10)
    secs = int(dt%3600%60//1)
    sec_1 = int(secs//10)
    sec_2 = int(secs%10)
    print('\r%s %s %s%% %s' %
          (prefix, 
           color + bar + colors.END, 
           color + percent + colors.END, 
           suffix+
           ", Time Used: " + 
           color + 
           str(hours)+':'+str(min_1)+str(min_2)+':'+str(sec_1)+str(sec_2) + 
           colors.END + '  '),
          end ='\r')


# =============================================================================
# 
# =============================================================================
print("----- progressbar.py is loaded -----")

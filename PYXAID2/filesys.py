#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os.path
import os

def check(prefix,suffix,mini,maxi):
# The function creates the continuity regions of existing files
# prefix - prefix of the file to be found
# suffix - suffix of the file to be found
# mini - minimal index of the file to be found
# maxi - maximal index of the file to be found


    res = []
    lst = []
    stop = 0
    i = mini
    while i<=maxi:
        file_path = prefix + str(i) + suffix
        if os.path.exists(file_path):
            lst.append(i)
        else:
            if len(lst)>0:
                res.append(lst)
            lst = []

        i = i + 1

    for lst in res:
        print "range("+str(lst[0])+","+str(lst[-1]+1)+")       nelts = "+str(lst[-1]+1-lst[0])



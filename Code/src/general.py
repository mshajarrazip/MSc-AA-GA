""" Store general functions that may be needed anywhere in the program.

Methods(s):
    disp(*arg) - Prints arg to the console. Use(s): Mostly for progress messages.
    disp_dict(table) - Display values in dictionaries in a nice layout.

"""
from random import random
from time import time, localtime


def disp(*arg):
    msg = ''
    for item in arg:
        if type(item) == list:
            str1 = "|"
            for w in item:
                str1 += (str(w) + "|")
            item = str1
        msg += str(item) + " "
    print(msg)

def disp_head(*arg):
    msg = ''
    for item in arg:
        msg += str(item) + " "
    trailing = ">>>>>"
    msg = trailing + msg
    print(msg)

def disp_dict(table, detail = False):
    """Prints the keys and values in a dictionary

    Keyword argument(s):
        table
            A dictionary.
        detail
            If the value of dictionary is a list, we check if it is a list of lists.
            Value(s):
            > False (Default) - Display only the size of the value list.
            > True - Display all the lists in value in separate line.

    """
    for k in table:
        val = table[k]
        if (type(val) == list):
            if (detail == True):
                if (all(isinstance(elem, list) for elem in val)):
                    str1 = ''.join("\t"+str(row)+"\n" for row in val)
                    val = str1
            else:
                val = "{0} of size {1}".format(type(val), len(val))
        disp(k,':', val)

def disp_status(*arg):
    msg = ''
    for item in arg:
        msg += str(item) + " "
    trailing = "............."
    msg = trailing + msg
    print(msg)

def random_decision(Pr, N=1):
    '''Return True if success based on probability Pr for N events.

    Parameter(s):
        Pr:float
            - Probability of success.
        N:int (Default = 100)
            - Number of events

    '''
    decision = False
    event = N
    heads = 0
    if (Pr > 0):
        heads = 0
        for t in range(event):
            if (random() < Pr):
                heads += 1
    if (heads >= 1):
        decision = True
    return decision

def timestamp():
    sec = time()
    T = localtime(sec)
    tstamp = ""
    for i in range(6):
        t = str(T[i])
        if (len(t)==1):
            t = "0" + t
        tstamp +=  t
    return (tstamp)


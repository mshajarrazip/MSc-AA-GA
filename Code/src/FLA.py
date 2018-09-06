"""Fitness Landscape Analyser"""

from json import load as jload
from math import sqrt, pow

from GAcomp import GApool, chromosome
from problem import problem


def FDC(pool, MIPdata, prob):
    """Fitness distance correlator.

    Calculates the FDC between the population and s', and returns the coefficient r.
    Based on Finger(2002)

    Parameter(s):
        pool:GApool
            - The population pool
        MIPdata:dict
            - The MIPdata. We assume that the MIP JSON file has been parsed, so we only pass ths resulting dictionary into this.
        prob:dict
            - dict of problem info data.
    """

    #converts MIP to optimal string
    OPT = chromosome(prob)
    OPT.genotype = MIPdata["x"]
    OPT.ID = "opt"
    F = dict()
    D = dict()

    for chr in pool:
        chr.setHD(OPT)
        D[chr.ID] = chr.HD["opt"]
        chr.setFD(OPT)
        F[chr.ID] = chr.fitness
        #print ("{0} - {1}".format(chr.HD["opt"], chr.FD["opt"]))
    '''
    count = 0
    for kw in F:
        count += 1
        print("F{1}:{0}".format(F[kw], count))

    count = 0
    for kw in D:
        count += 1
        print("D{1}:{0}".format(D[kw], count))
    '''

    uF = sum([float(F[kw]) for kw in F]) / float(len(F))
    uD = sum([float(D[kw]) for kw in D]) / float(len(D))
    #print(uF, uD)

    C_FD = sum([(float(F[kw])-uF)*(float(D[kw])-uD) for kw in F])/ (float(len(F))-1)
    #print(C_FD)

    sF = sqrt(sum([pow(float(F[kw]) - float(uF),2) for kw in F]) / (float(len(F))-1))
    sD = sqrt(sum([pow(float(D[kw]) - float(uD),2) for kw in D]) / (float(len(D))-1))
    #print(sF, sD)

    r = C_FD/(sF*sD)
    print("r is {0}".format(r))

    del(OPT)
    del(F)
    del(D)
    return r
'''
prob = problem()
prob.load("C:\Users\Hajar\Documents\UTP II\MASTER\Python_Codes\MUTP_Exp_AAGA-remastered\ORLib\p4\scp41.txt")

MIPpath = "C:\Users\Hajar\Documents\UTP II\MASTER\Python_Codes\MUTP_Exp_AAGA-remastered\MIP_JSON\scp41_MIP_f0_v429.0_20170225133001_JSON.txt"

with open(MIPpath) as json_data:
    jdict = jload(json_data)


FDC(1, jdict, prob)
'''
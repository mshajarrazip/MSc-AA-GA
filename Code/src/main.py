"""Main src file.


example:
source = '../ORLib/scp51.txt' # Set the source file
GA = GA() # Create an instance of GA
GA.setProblem(source) # Bind source data file to GA instance
GA.setParameter(N=20, gen=10000, mut_param = (10, 200, 2.0), init_method="AA")
    # Set GA parameters
#if init_method = AAR, then AAparam = ("AAR", 0.5, 1)
GA.solve()

"""

# currently running without N-1
from os import listdir, makedirs
from os.path import join, isfile,basename, exists
from random import seed as randoseed
from json import load as jload
from csv import DictReader as csvread
from sys import exit, argv
from time import time

from general import timestamp, disp_head
from GAsolver import GA
from data_manager import solution


if (len(argv) < 5):
    print("Usage: python main.py <<name of main problem data folder>> <<subfolder>> <<population data file>> << poptag >>")
    print("Exiting ... ")
    exit()

datasubfile = argv[1]
popdatafile = argv[3]
poptag = argv[4]
#for probsetname in ["p4"]:#,"p5", "p6", "pa", "pb", "pc", "pd", "pe"]:
probsetfile = argv[2] # pass the name of subfile at command line
for probsetname in [probsetfile]:

    # parameter settings

    # >>>>>>>> mutation rate: in /ORLib, edit mut_rate.csv to modify settings
    with open(join("..",datasubfile,"mut_rate.csv")) as mut_file:
        mut_data = csvread(mut_file)
        m_ = [row for row in mut_data if (row["problem"] == probsetname)][0]
	print(m_)
    MUT_PARAM = (float(m_["mf"]),float(m_["mc"]),float(m_["mg"]))

    # >>>>>>>>> Initialization method
    INIT_METHOD = "LoadPopulation"
    #poptag = "f" # input file would be named scp#_<tag>... , please insert what is in teh place of <tag>
    MaxMode = True

    # >>>>>>>>> Problem data source file settings
    probsrc = join("..", datasubfile, probsetname)
    probsrclist = [join(probsrc, f) for f in listdir(probsrc) if isfile(join(probsrc, f))]
    probnamelist = [basename(f)[:basename(f).index(".txt")] for f in probsrclist]
    #probnamelist = ["scp43"];

    # >>>>>>>>> for load population
    inputdir = join("..", "..", "Input")
    if (exists(inputdir) == False):
        makedirs(inputdir)
        print("Input folder missing \n I've created one for you. \n Please put population data folder there. \n")
        exit()

    popdatadir = ""
    if (INIT_METHOD == "LoadPopulation"):
        popdatadir = popdatafile
        popdir = join(inputdir, popdatadir)
        if (exists(popdir) == False):
            print("Population data folder is missing. \n Please add folder in Input.")
            exit()
        popdatalist = [join(popdir, f) for f in listdir(popdir) if isfile(join(popdir, f))]

    # parameter settings for experiment
    trial = 10
    popsize = 100
    genmax = 100000 # termination condition (100000)
    gen_lim = 10000 # Determines if a GA has already converged (e.g. avg fitness stuck at certain value for max how long?)(10000)
    nseed = 0.0 # set seed size [0.0~1.0]

    # other settings
    MIPstat = False # if False, MIP file will not be read for FDC, optval set to 0
    BFStoJSON = False# if true, the best GA solution in t trials is printed to a JSON file
    optval = 0


    # >>>>>>>>> Output directory setup
    outputdir = join("..", "..", "Output")
    if (exists(outputdir) == False):
        makedirs(outputdir)


    TS = timestamp()
    logfolder = "LOG_GA_RUN_{0}_{1}_{2}_s{3}_{4}".format(INIT_METHOD, popdatadir, popsize, nseed, TS)
    logdir = join(outputdir, logfolder)
    print("Log files will be saved in: {0}".format(logdir))
    makedirs(logdir)

    # Output for best BFS OBJ
    if (INIT_METHOD == "Random"):
        popdatadir = "Random"
    outOBJ = join(logdir, "mainlog_{3}_OBJ_{0}_s{1}_{2}.csv".format(popdatadir, genmax,
                                                                             str(TS), probsetname))
    # Output for best BFS first gen found
    outGEN = join(logdir, "mainlog_{3}_GEN_{0}_s{1}_{2}.csv".format(popdatadir, genmax,
                                                                         str(TS), probsetname))
    # Output for seeds for each run
    outSEED = join(logdir, "mainlog_{3}_SEED_{0}_s{1}_{2}.csv".format(popdatadir, genmax,
                                                                             str(TS), probsetname))

    # make header for output file
    header_str = "Problem"
    for j in range(trial):
        header_str += ",{0}".format("t{0}".format(str(j+1)))
    header_str += ",BEST,AVG\n"

    with open(outOBJ, "a+") as f:
        f.write(header_str)
    with open(outGEN, "a+") as f:
        f.write(header_str)
    with open(outSEED, "a+") as f:
        f.write(header_str)

    for pr in probnamelist:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        disp_head("Problem: {0}".format(pr))
        prdatafile = [f for f in probsrclist if (pr + ".txt") in f][0]
        print("Problem data file: {0}".format(prdatafile))
        if (INIT_METHOD == "LoadPopulation"):
            popdatafile = []
            try:
                popdatafile = [f for f in popdatalist if (pr + "_") in f][0]
            except IndexError:
                print("POP DATA FILE NOT FOUND (%s): Proceeding to next problem ..." % pr)
                continue
            print("Population data file: {0}".format(popdatafile))
        else:
            popdatafile = ""

        mipdata = {}
        if (MIPstat == True):
            MIPsrc = join("..", "GA_BFS_H")
            MIPdatafile = [join(MIPsrc, f) for f in listdir(MIPsrc) if ((isfile(join(MIPsrc, f)) and ((pr + "_") in f)))][0]
            print("MIP JSON file: {0}".format(MIPdatafile))
            # get exact solution value
            with open(MIPdatafile) as json_data:
                mipdata = jload(json_data)
            optval = mipdata['OBJ_VAL']
            print("\t OPT val: {0}".format(optval))
            print("Output logs:")
            print("\t OBJ: {0}".format(outOBJ))
            print("\t GEN: {0}".format(outGEN))
            print("\t SEED: {0}".format(outSEED))

        listOBJ = list()
        listGEN = list()
        # chromosome list to keep track of BFS in each trial
        listBFS = list()

        status = 0

        for i in range(trial):
            print(">>>Trial {0}\n", i)

            seedval = time()
            randoseed(seedval)

            Ga = GA()
            Ga.setProblem(prdatafile)
            Ga.setParameter(N=popsize, N_seed= nseed, gen=genmax, gen_max = gen_lim, opt_val = optval,
                            mut_param=MUT_PARAM, init_method=INIT_METHOD,
                            popsource=popdatafile,
                            MIPdata = mipdata,
                            logDir = logdir,
                            poptags = poptag)

            # if init_method = AAR, then AAparam = ("AAR", 0.5, 1)
            status = Ga.solve()

            if (status > 0):
                print("NOTICE: Population size not sufficient for seeding. Skipping problem {0}".format(pr))
                break

            if (i == 0):
                with open(outOBJ, "a+") as f:
                    f.write("{0},".format(pr))
                with open(outGEN, "a+") as f:
                    f.write("{0},".format(pr))
                with open(outSEED, "a+") as f:
                    f.write("{0},".format(pr))

            listOBJ.append(float(Ga.population.BFS.fitness))
            listGEN.append(float(Ga.population.BFS.ID[Ga.population.BFS.ID.index("g")+1:Ga.population.BFS.ID.index("ge")]))
            listBFS.append(Ga.population.BFS)

            with open(outOBJ, "a+") as f:
                f.write("{0},".format(Ga.population.BFS.fitness))
            with open(outGEN, "a+") as f:
                f.write("{0},".format(Ga.population.BFS.ID[Ga.population.BFS.ID.index("g")+1:Ga.population.BFS.ID.index("ge")]))
            with open(outSEED, "a+") as f:
                f.write("{0},".format(seedval))

            del(Ga)

        if (BFStoJSON == True):
            # get fittest chromosome in all trials
            GAbest = listBFS[0]
            for f in listBFS:
                if (f.fitness < GAbest.fitness):
                    GAbest = f
            print("... BFS in {0} trials is {1}".format(trial, GAbest.ID))

            # print GAbest to JSON file
            gbdata = dict()
            mode = "GA_BFS"
            gbdata["OBJ_VAL"] = GAbest.fitness
            gbdata["X"] = GAbest.genotype
            gbdata["feasible"] = GAbest.status
            GAsol = solution(pr, mode, 0)
            GAsol.update(gbdata)
            GAsol.display()
            GAsol.writeasjson(outpdir=logdir, filename="{0}_{1}_{2}_{3}.txt".format(pr, mode, trial, TS))

        if (status == 0):
            with open(outOBJ, "a+") as f:
                bestOBJ = min(listOBJ)
                avgOBJ = sum(listOBJ)/float(len(listOBJ))
                f.write("{0},{1}\n".format(bestOBJ, avgOBJ))
            with open(outGEN, "a+") as f:
                bestGEN = min(listGEN)
                avgGEN = sum(listGEN)/float(len(listGEN))
                f.write("{0},{1}\n".format(bestGEN, avgGEN))
            with open(outSEED, "a+") as f:
                f.write("\n")


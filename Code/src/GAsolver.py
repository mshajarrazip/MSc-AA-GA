"""Genetic Algorithm Solver

Class(es):
    GA - GA class.

"""


from random import choice, randint
from math import ceil, exp
from os import makedirs, path


from problem import problem, info
from general import *
from GAcomp import *
from FLA import FDC as fla_FDC


class GA:
    """GA class.

    Attribute(s):
        prob:problem.info
            - Problem instance information.
        N:int
            - Size of population.
        init_method:str
            - Method of initial population generation.
        gen:int
            - Number of generations (before termination).

    Method(s):
        __init__ - Constructor.
        setProblem - Bind problem from source to GA.
    """
    def __init__(self):
        """Constructor.

        Initialize self.prob to an instance of problem.problem (Class that holds problem data).

        """

        self.prob = problem()

    def setProblem(self, source):
        """Load the data from source.

        Bind problem info to GA (Store problem info in self.prob).

        Parameter(s):
            source:str
                - Problem data file.

        """

        self.prob.load(source)
        pass

    def setParameter(self, N=0, init_method = "Random",
                     gen=0, gen_max = False, opt_val = None,
                     sol_log = '',
                     mut_param = (0,0,0),
                     AA_param = (),
                     popsource = "",
                     MIPdata = {},
                     logDir = "..",
                     N_seed = 0,
                     poptags = ""):

        """Parameter setting.

        Parameter(s):
            N:int
                - Size of population.
            init_method:str
                - Method of initial population generation.
                    "Random" - Random generation
                    "AA" - Approximation algorithm. Pass parameters to AA_param for further settings.
            gen:int
                - Number of generations (before termination).
            gen_max
                - Max generations where BFS remains unchanged
            con_log:str
                - File path for convergence data log.
            mut_param:(mf:int, mc:int, mg:int) (Default = (0,0,0))
                - See method: mutate
            AA_param:(Mode:str, factor:float, Randomizer:int, (Obj:list(int), A:list([{0-1}])
                - Parameters for initialization by AA.
                    Mode: Mode of solution generation (AAR/ AAD) (Default = "AAD")
                        *note that this setting only creates 1 child if second argument is 0
                        because AAD with no factor is deterministic.
                    factor: Range of randomization. (Default = 0)
                        * suggested 0.5
                    Randomizer: Mode of randomization. (1:Pseudorandom, 2:Sobol etc) (Default = 1)
            logDir: String
                - Where the output files go
            N_seed:int (Default = 0, range [0.0 - 1.0])
                - If N_seed > 0, then int(N_seed*seed) will be loaded from LoadPopulation, and the rest will be randomized
        """
        self.poptag = poptags
        self.N = N
        self.N_seed = N_seed
        self.init_method = init_method
        self.opt_val = opt_val
        self.gen = gen
        self.gen_max = gen_max
        self.mut_param = mut_param
        if init_method == "AA":
            if AA_param == ():
                AA_param = ("AAD", 0, 1)
            self.AA_param = AA_param
        elif init_method == "LoadPopulation":
            self.popsource = popsource
        self.MIPdata = MIPdata


        # Set name for output log file, and write to log file the headers
        def_rep = path.join(logDir, "IterationStatus")
        if not path.exists(def_rep):
            makedirs(def_rep)


        if sol_log == '':
            filename = ""
            if init_method == "Random":
                filename = "sol_log_{1}.csv".format(self.prob.name,
                                                                    timestamp())
            elif init_method == "AA":
                filename = "sol_log_{1}_{2}_{3}.csv".format(self.prob.name,
                                                                    self.AA_param[1],
                                                                    self.AA_param[2],
                                                                    timestamp())
            elif init_method == "LoadPopulation":
                print(popsource)
                filename = "sol_log_{0}_{1}_{2}_{3}_{4}{5}.csv".format(self.prob.name,
                                                            self.popsource[popsource.index("scp"):
                                                                           popsource.index("_{0}".format(self.poptag))],
                                                            mut_param[0],
                                                            mut_param[1],
                                                            mut_param[2],
                                                            str(time()))
            filename = path.join(def_rep, filename)
            disp("Log file path: ", filename)
            self.sol_log = filename
        else:
            self.sol_log = sol_log

        out_log = self.sol_log
        with open(out_log, "a+") as f:
            f.write("Gen, Avg_Fitness, Mut_Rate, BFS_ID, BFS_val\n")

        # output message
        disp_head("Parameters set:")
        disp("\tPopulation size:", self.N)
        disp("\tSeed rate:", self.N_seed)
        disp("\tInitialization method:", self.init_method)
        if init_method == "AA":
            disp("\t\tAA parameters:{0},{1},{2}".format(AA_param[0], AA_param[1], AA_param[2]))
        disp("\tGen limit:", self.gen)
        disp("\tMutation parameters:", self.mut_param)

    def solve(self):
        """Solve the GA.

        Initialize the population according to init_method (Random/AA).


        """

        prob = self.prob
        init_method = self.init_method
        pop = population()
        if init_method == "Random":
            pop.initialize(prob, N=self.N)
        elif init_method == "AA":
            AA_param = self.AA_param
            pop.initialize(prob, N=self.N, init_method = "AA", AA_param = AA_param)
        elif init_method == "LoadPopulation":
            popsource = self.popsource
            status = pop.initialize(prob, init_method = "LoadPopulation", PopSource = popsource, N_seed = self.N_seed)
            if (status > 0):
                return status
            disp("Source: {0}".format(popsource))
            disp("Population size: {0}".format(len(pop.pool)))

        if (self.MIPdata != {}):
            self.FDC = fla_FDC(pop.pool, self.MIPdata, self.prob)

        avg_fitness = []

        gen_pool = GApool()
        gn = len(gen_pool)
        self.gen_n = gn
        gen_lim = self.gen_max
        optval = self.opt_val
        bfscount = 0
        bfstrack = None
        disp_status("Begin evolution")
        while ((gn < self.gen) and (bfscount < gen_lim)):

            # disp("current generation: ", gn)
            parent = self.tournament(pop)
            self.gen_n = gn
            child = self.crossover(parent)
            if ((gen_pool.ismember(child) != True) and (pop.pool.ismember(child) != True)):
                child.setID(gn)
                child.printinfo()
                pop.insert(child)
                gen_pool.append(child)
                avg_fitness.append(pop.avg_fitness)
                if (pop.BFS.ID == bfstrack):
                    bfscount += 1
                else:
                    bfstrack = pop.BFS.ID
                    bfscount = 0
                if (pop.BFS.fitness <= optval):
                    break
            gn = len(gen_pool)

            #format of console output
            dvder = 1000
            dvder2 = 100
            if ((gn%dvder == 0) | (bfscount%dvder2 == 0)):
                disp_status("Gen pool: {0} Same BFS: {1}".format(gn, bfscount))

            with open(self.sol_log, "a+") as f:
                f.write("{0},{1},{2},{3},{4}\n".format(gn, pop.avg_fitness, self.mut_rate*10, pop.BFS.ID, pop.BFS.fitness))

        self.population = pop
        disp_status("done")
        disp("Best solution found: ", pop.BFS.ID, pop.BFS.fitness)
        return 0

    def tournament(self, pop, T=2):
        """Return a tuple of two best parents from a torunament selection of size T.

        Parameter(s):
            pop:population
                - The population.
            T:int (Default = 2)
                - The size of the tournament pools.
        """

        #disp_head("TOURNAMENT SELECTION")
        #disp("T:", T)

        tpool = list(pop.pool)

        T1 = []
        T2 = []

        while len(T1) < 2:
            parent = choice(tpool)
            T1.append(parent)
            tpool.remove(parent)
        while len(T2) < 2:
            parent = choice(tpool)
            T2.append(parent)
            tpool.remove(parent)

        parent1 = best_chr(T1)
        parent2 = best_chr(T2)

        #disp("T1 pool:")
        for i in T1:
            fitness = i.fitness
            if i == parent1:
                fitness = str(fitness) + "*"
            #disp(fitness)
        #disp("T2 pool:")
        for i in T2:
            fitness = i.fitness
            if i == parent2:
                fitness = str(fitness) + "*"
            #disp(fitness)

        return (parent1, parent2)

    def crossover(self, parent):
        """Return a child from a crossover of both parents using the 'fusion operator'.

        Parameter(s):
            parent:tuple(chromosome(),chromosome())
                - A tuple of two chromosomes.

        """
        #disp_head("PERFORMING CROSSOVER")

        p1 = parent[0].genotype
        p2 = parent[1].genotype
        f1 = parent[0].fitness
        f2 = parent[1].fitness

        fp1 = 0
        fp2 = 0

        C = [0] * len(p1)
        for j in range(len(C)):
            if p1[j] == p2[j]:
                C[j] = p1[j]
                continue
            P = float(f2) / (f1 + f2)
            Q = 1 - P
            flag = random_decision(P)
            if flag == True:
                C[j] = p1[j]
                fp1 += 1
            else:
                C[j] = p2[j]
                fp2 += 1

        C = self.mutate_BeChu(C, self.mut_param)

        child = chromosome(parent[0].info) #fix here
        child.genotype = C
        #disp("Child produced: S {0} | genotype {1} | fitness {2} | feasibility {3}".format(len(child.S), child.genotype.count(1), child.fitness, child.status))
        #disp("{0} from P1 {1} from P2".format(fp1, fp2))
        return child

    def mutate_BeChu(self, C, mut_param = (0,0,0)):
        """Mutate k randomly selected columns in C where k is determined by the variable mutation schedule

        Parameter(s):
            C:list({0,1})
                - A genotype.
            mut_param:(mf:int, mc:int, mg:int) (Default = (0,0,0))
                Mf:int (Default = 0, no mutation applied by default)
                    - The final stable mutation rate (user-defined).
                Mc:int (Default = 0)
                    - The number of child solutions generated at which mutation rate of mf/2 is reached
                Mg:int (Default = 0)
                    - The gradient at t=mc
            *These values are provided in Beasley's 1996 paper

        """
        Mf = mut_param[0]

        if (Mf != 0):
            Mc = mut_param[1]
            Mg = mut_param[2]
            t  = self.gen_n # the number of child solutions that have been generated
            M = ((-4*Mg)*(t-Mc))/Mf
            X = 1+exp(M)
            rate = ceil( Mf/X)
            self.mut_rate = rate

            Ji = self.prob.Ji
            I = self.prob.I
            k_elite = 5

            S_elite = {Ji[i][j] for i in I for j in range(k_elite)}
            S_elite = sorted(list(S_elite)) # Bits to mutate only among this set

            mut_bits = []

            while len(mut_bits) < rate:
                mut_bits.append(S_elite.pop(randint(0, len(S_elite)-1)))

            for j in mut_bits:
                C[j] -= 1
                if C[j] < 0:
                    C[j] = 0
        else:
            self.mut_rate = 0

        return C


"""
Performs the parallel execution of the QM calls.
"""
import time
import os
from multiprocessing import Process
from typing import Any, Dict, List, Union
from concurrent.futures import ProcessPoolExecutor

from debug.censo_test.qm_processor import ProcessorFactory, QmProc
from debug.censo_test.tm_processor import TmProc
from debug.censo_test.orca_processor import OrcaProc
from censo_test.utilities import print
from censo_test.cfg import ENVIRON, WARNLEN
from censo_test.datastructure import MoleculeData
from censo_test.settings import CensoSettings


class ProcessHandler:
    
    # https://stackoverflow.com/questions/1540822/dumping-a-multiprocessing-queue-into-a-list/1541117#1541117
    
    def __init__(self, settings: CensoSettings, conformers: List[MoleculeData]):
        """
        creates instance with 
        also sets the number of processes (nprocs) and cores per processes (omp) to optimized vaules if possible
        """
        # stores the conformers
        self.conformers: List[MoleculeData] = conformers
        
        # stores the processor
        self._processor = None
        
        # get number of cores
        try:
            self.ncores = os.cpu_count()
        except AttributeError:
            raise AttributeError("Could not determine number of cores.")
        
        self.nprocs = None
        self.omp = None
        if (
            getattr(settings.settings_current.byname("balance"), "value", False) 
            or not self._balance_load()
        ):
            print("\nCould not readjust maxprocs and omp.\n")
            # get number of processes
            self.nprocs = getattr(
                settings.settings_current.byname("maxprocs"), 
                "value",
                None
            )
            
            # get number of cores per process
            self.omp = getattr(
                settings.settings_current.byname("omp"),
                "value",
                None
            )
    
    
    def _balance_load(self) -> bool:
        """
        distribute computational load between cores
        keeping the number of processes below number of conformers
        trying to utilize as many cores as possible
        """
        # check if load can be rebalanced
        if (
            not self.ncores is None 
            and not (
                self.nprocs is None or self.omp is None
            )
            and (
                self.ncores > len(self.conformers) 
                or self.ncores > self.nprocs * self.omp
            )
        ):
            # you want to create a number of parallel processes that in the best case is a divisor of nconf
            # the number of processes should allow for utilization of all cores
            try:
                self.nprocs = max([i for i in range(self.ncores) if self.ncores % i and i <= self.jobq.qsize()])
            except ValueError:
                print("There was an error while determining the number of processes in load balancing.") # TODO
                return False
            
            # should always be divisible, int casting only as safeguard
            self.omp = int(self.ncores / self.nprocs)
            
            return True
        else:
            return False
                

    def execute(self, jobtype: List[str], prog: str, instructions: Dict[str, Any]) -> List:
        """
        creates and executes the processes
        returns the results
        """
        res = []
        
        # initialize the processor for the respective program (depends on part)
        # and set the jobtype as well as instructions
        self._processor = ProcessorFactory.create_processor(prog)
        self._processor.jobtype = jobtype
        self._processor.instructions = instructions
        
        # execute jobs for conformers in queue
        # check if nprocs is a divisor of nconf
        # TODO - if this is not the case split jobq into two parts, rebalance load for the second part
        with ProcessPoolExecutor(max_workers=self.nprocs) as executor:
            res = [item for item in executor.map(self._processor.run, self.conformers)]

        return res


    def run_in_parallel(
        config,
        q,
        resultq,
        job,
        maxthreads,
        omp,
        loopover,
        instructdict,
        balance=False,
        foldername="",
    ):
        """Run jobs in parallel
        q = queue to put assemble tasks
        resultq = queue to retrieve results
        job = information which kind of job is to be performed tm_job , orca_job
        loopover is list of qm_class objects
        instructdict example : {'jobtype': 'prep', 'chrg': args.chrg}
        foldername is for existing objects to change the workdir
        results = list of qm_class objects with results from calculations
        
        
        how does it work?
        - takes jobs as input, balances load first
        - creates a number of processes, with run=execute_data, which gets tasks from the q and runs them,
        putting the results in resultq
        - waits for q to be finished, then puts the results from resultq into a list and returns the list
        - results are the jobs themselves with the actual results being stored in attributes of them (?!???!)
        """
        omp_initial = omp
        if instructdict.get("jobtype", None) is None:
            raise KeyError("jobtype is missing in instructdict!")
        if len(loopover) != 0:
            nconf = len(loopover)
        else:
            balance = False
        maxthreads, omp, changed = balance_load(maxthreads, omp, nconf, balance)
        if all(isinstance(x, QmProc) for x in loopover):
            for item in loopover:
                if isinstance(item, TmProc) and job == OrcaProc:
                    item.__class__ = job
                elif isinstance(item, OrcaProc) and job == TmProc:
                    item.__class__ = job
                elif isinstance(item, QmProc) and job != QmProc:
                    item.__class__ = job
                item.job["workdir"] = os.path.normpath(
                    os.path.join(config.cwd, "CONF" + str(item.id), foldername)
                )
                # update instructions
                item.job.update(instructdict)
                item.job["omp"] = omp
                # put item on queue
                q.put(item)
                time.sleep(0.02)
        time.sleep(0.02)
        njobs = q.qsize()
        if changed:
            ENVIRON["PARNODES"] = str(omp)
        if njobs != nconf:
            print(f"{'WARNING':{WARNLEN}}Number of conformers submitted to queue ({nconf}) and in queue ({njobs}) do not match,"
                f"{'':{WARNLEN}}but this can also be due to an approximate number of conformers from q.qsize() in python!"
                )
        if instructdict.get("onlyread", False):
            print(f"\nReading data from {njobs} conformers calculated in previous run.")
        else:
            response = {
                "prep": f"\nPreparing {njobs} calculations.",
                "sp": f"\nStarting {njobs} single-point calculations.",
                "xtb_sp": f"\nStarting {njobs} xTB - single-point calculations.",
                "lax_sp": f"\nStarting {njobs} lax-single-point calculations.",
                "cosmors": f"\nStarting {njobs} COSMO-RS-Gsolv calculations.",
                "gbsa_gsolv": f"\nStarting {njobs} GBSA-Gsolv calculations",
                "alpb_gsolv": f"\nStarting {njobs} ALPB-Gsolv calculations",
                "smd_gsolv": f"\nStarting {njobs} SMD-Gsolv calculations",
                "rrhoxtb": f"\nStarting {njobs} G_RRHO calculations.",
                "rrhoorca": f"\nStarting {njobs} G_RRHO calculations.",
                "rrhotm": f"\nStarting {njobs} G_RRHO calculations.",
                "opt": f"\nStarting {njobs} optimizations.",
                "xtbopt": f"\nStarting {njobs} optimizations.",
                "couplings": f"\nStarting {njobs} coupling constants calculations",
                "couplings_sp": f"\nStarting {njobs} coupling constants calculations",
                "shieldings": f"\nStarting {njobs} shielding constants calculations",
                "shieldings_sp": f"\nStarting {njobs} shielding constants calculations",
                "genericoutput": f"\nWriting {njobs} generic outputs.",
                "opt-rot": f"\nStarting {njobs} optical-rotation calculations.",
                "opt-rot_sp": f"\nStarting {njobs} optical-rotation calculations.",
            }
            if instructdict["jobtype"] in response:
                print(response[instructdict["jobtype"]])
        time.sleep(0.04)
        # start working in parallel
        for _ in range(int(maxthreads)):
            worker = Process(target=execute_data, args=(q, resultq))
            worker.daemon = True
            worker.start()
            # NOBODY IS ALLOWED TO TOUCH THIS SLEEP STATEMENT!!!!
            time.sleep(0.05)  # sleep is important don't remove it!!!
            # seriously don't remove it!
        q.join()

        if not instructdict.get("onlyread", False):
            print("Tasks completed!\n")
        else:
            print("Reading data from previous run completed!\n")

        # Get results
        results = []
        while not resultq.empty():
            results.append(resultq.get())
            if getattr(results[-1], "hugeERROR", False):
                print(getattr(results[-1], "tb"))
                raise getattr(results[-1], "hugeERROR")
            time.sleep(0.02)  # sleep is important don't remove it!!!
            # seriously don't remove it!

        time.sleep(0.02)
        results.sort(key=lambda x: int(x.id))
        if nconf != len(results):
            print(f"{'ERROR':{WARNLEN}}number of confs before ({nconf}) and after ({len(results)}) run_in_parallel do not match !")
        elif njobs != len(results):
            print(f"{'WARNING':{WARNLEN}}number some conformers do not match ({njobs} / {len(results)}) it appears as if some conformers were lost,"
                f"{'':{WARNLEN}}but this can also be due to an approximate number of conformers from q.qsize() in python!"
                )
        if changed:
            ENVIRON["PARNODES"] = str(omp_initial)
        return results

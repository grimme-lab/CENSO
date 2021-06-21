"""
Performs the parallel execution of the QM calls.
"""
import time
import os
import traceback
from multiprocessing import Process
from .qm_job import QmJob
from .tm_job import TmJob
from .orca_job import OrcaJob
from .utilities import print
from .cfg import ENVIRON


def balance_load(P, O, nconf, do_change):
    """Balance calculation load between threads (P) and number of cores per 
    thread (O)
    """
    changed = False
    max_cores = P * O
    if do_change:
        if nconf < P:
            try:
                P = nconf
                O = 1
                while True:
                    if P * O <= max_cores:
                        if P * (O + 1) <= max_cores:
                            O += 1
                        else:
                            break
                    else:
                        break
                changed = True
            except Exception:
                pass
        if changed:
            print(
                f"Adjusting the number of threads (P) = {P} and number of cores per thread (O) = {O}"
            )
    return P, O, changed


def execute_data(q, resultq):
    """
    code that the worker has to execute
    """
    while True:
        if q.empty():
            break
        task = q.get()
        try:
            task.execute()
        except Exception as e:
            print(e)
            task.hugeERROR = e
            task.tb = traceback.format_exc()
            resultq.put(task)
            q.task_done()
            if q.empty():
                break
            else:
                task = q.get()
                resultq.put(task)
                q.task_done()
        resultq.put(task)
        time.sleep(0.02)
        q.task_done()
        time.sleep(0.02)


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
    """
    omp_initial = omp
    if instructdict.get("jobtype", None) is None:
        raise KeyError("jobtype is missing in instructdict!")
    if len(loopover) != 0:
        nconf = len(loopover)
    else:
        balance = False
    maxthreads, omp, changed = balance_load(maxthreads, omp, nconf, balance)
    if all(isinstance(x, QmJob) for x in loopover):
        for item in loopover:
            if isinstance(item, TmJob) and job == OrcaJob:
                item.__class__ = job
            elif isinstance(item, OrcaJob) and job == TmJob:
                item.__class__ = job
            elif isinstance(item, QmJob) and job != QmJob:
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
    njobs = q.qsize()
    time.sleep(0.02)
    if changed:
        ENVIRON["PARNODES"] = str(omp)

    if instructdict.get("onlyread", False):
        print(f"\nReading data from {njobs} conformers calculated in previous run.")
    else:
        response = {
            "prep": f"\nPreparing {q.qsize()} calculations.",
            "sp": f"\nStarting {q.qsize()} single-point calculations.",
            "xtb_sp": f"\nStarting {q.qsize()} xTB - single-point calculations.",
            "lax_sp": f"\nStarting {q.qsize()} lax-single-point calculations.",
            "cosmors": f"\nStarting {q.qsize()} COSMO-RS-Gsolv calculations.",
            "gbsa_gsolv": f"\nStarting {q.qsize()} GBSA-Gsolv calculations",
            "alpb_gsolv": f"\nStarting {q.qsize()} ALPB-Gsolv calculations",
            "smd_gsolv": f"\nStarting {q.qsize()} SMD-Gsolv calculations",
            "rrhoxtb": f"\nStarting {q.qsize()} G_RRHO calculations.",
            "rrhoorca": f"\nStarting {q.qsize()} G_RRHO calculations.",
            "rrhotm": f"\nStarting {q.qsize()} G_RRHO calculations.",
            "opt": f"\nStarting {q.qsize()} optimizations.",
            "xtbopt": f"\nStarting {q.qsize()} optimizations.",
            "couplings": f"\nStarting {q.qsize()} coupling constants calculations",
            "couplings_sp": f"\nStarting {q.qsize()} coupling constants calculations",
            "shieldings": f"\nStarting {q.qsize()} shielding constants calculations",
            "shieldings_sp": f"\nStarting {q.qsize()} shielding constants calculations",
            "genericoutput": f"\nWriting {q.qsize()} generic outputs.",
            "opt-rot": f"\nStarting {q.qsize()} optical-rotation calculations.",
            "opt-rot_sp": f"\nStarting {q.qsize()} optical-rotation calculations.",
        }
        if instructdict["jobtype"] in response:
            print(response[instructdict["jobtype"]])

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
        time.sleep(0.01)  # sleep is important don't remove it!!!
        # seriously don't remove it!

    time.sleep(0.02)
    results.sort(key=lambda x: int(x.id))
    if njobs != len(results):
        print(f"ERROR some conformers were lost!")

    if changed:
        ENVIRON["PARNODES"] = str(omp_initial)
    return results

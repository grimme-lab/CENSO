"""
Utility functions which are used in the CENSO modules. From creating folders to
printout routines.
"""
import os
import sys
import shutil
import math
import hashlib
import time
import subprocess
from copy import deepcopy
from builtins import print as print_orig
from .cfg import ENVIRON, CODING, AU2J, AU2KCAL, BOHR2ANG, KB, WARNLEN


def print(*args, **kwargs):
    """
    patch print to always flush
    """
    sep = " "
    end = "\n"
    file = None
    flush = True
    for key, value in kwargs.items():
        if key == "sep":
            sep = value
        elif key == "end":
            end = value
        elif key == "file":
            file = value
        elif key == "flush":
            flush = value
    print_orig(*args, sep=sep, end=end, file=file, flush=flush)


def frange(start, end, step=1):
    """
    range with floats
    """
    try:
        start = float(start)
        end = float(end)
        step = float(step)
    except (ValueError, TypeError):
        raise
    if start > end:
        tmp = start
        start = end
        end = tmp
    count = 0
    while True:
        temp = float(start + count * step)
        if temp >= end:
            break
        yield temp
        count += 1


def mkdir_p(path):
    """
    create mkdir -p like behaviour
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if not os.path.isdir(path):
            raise e


def print_block(strlist, width=80):
    """Print all elements of strlist in block mode
    e.g. within 80 characters then newline
    - width [int] width of block
    """
    length = 0
    try:
        maxlen = max([len(str(x)) for x in strlist])
    except (ValueError, TypeError):
        maxlen = 12
    for item in strlist:
        length += maxlen + 2
        if length <= width:
            if not item == strlist[-1]:  # works only if item only once in list!
                print("{:>{digits}}, ".format(str(item), digits=maxlen), end="")
            else:
                print("{:>{digits}}".format(str(item), digits=maxlen), end="")
        else:
            print("{:>{digits}}".format(str(item), digits=maxlen))
            length = 0
    if length != 0:
        print("\n")


def t2x(path, writexyz=False, outfile="original.xyz"):
    """convert TURBOMOLE coord file to xyz data and/or write *.xyz ouput

     - path [abs. path] does not need to include the filename coord
     - writexyz [bool] default=False, directly write to outfile
     - outfile [filename] default = 'original.xyz' filename of xyz file which
                        is written into the same directory as
     returns:
     - coordxyz --> list of strings including atom x y z information
     - number of atoms
    """
    if not os.path.basename(path) == "coord":
        if os.path.isfile(path):
            with open(path, "r", encoding=CODING, newline=None) as f:
                if not '$coord' in f.readline():
                    path = os.path.join(path, "coord")
        else:
            path = os.path.join(path, "coord")
    with open(path, "r", encoding=CODING, newline=None) as f:
        coord = f.readlines()
    x = []
    y = []
    z = []
    atom = []
    for line in coord[1:]:
        if "$" in line:  # stop at $end ...
            break
        x.append(float(line.split()[0]) * BOHR2ANG)
        y.append(float(line.split()[1]) * BOHR2ANG)
        z.append(float(line.split()[2]) * BOHR2ANG)
        atom.append(str(line.split()[3].lower()))
    # natoms = int(len(coord[1:-1])) # unused
    coordxyz = []
    for i in range(len(x)):
        coordxyz.append(
            "{:3} {: .10f}  {: .10f}  {: .10f}".format(
                atom[i][0].upper() + atom[i][1:], x[i], y[i], z[i]
            )
        )
    if writexyz:
        with open(
            os.path.join(os.path.split(path)[0], outfile),
            "w",
            encoding=CODING,
            newline=None,
        ) as out:
            out.write(str(len(coordxyz)) + "\n\n")
            for line in coordxyz:
                out.write(line + "\n")
    return coordxyz, int(len(coordxyz))


def x2t(path, infile="inp.xyz"):
    """convert file inp.xyz to TURBOMOLE coord file"""
    if ".xyz" not in os.path.basename(path):
        path = os.path.join(path, infile)
    with open(path, "r", encoding=CODING, newline=None) as f:
        xyz = f.readlines()
        atom = []
        x = []
        y = []
        z = []
        for line in xyz[2:]:
            atom.append(str(line.split()[0].lower()))
            x.append(float(line.split()[1]) / BOHR2ANG)
            y.append(float(line.split()[2]) / BOHR2ANG)
            z.append(float(line.split()[3]) / BOHR2ANG)
        coordxyz = []
        for i in range(len(x)):
            coordxyz.append(f"{x[i]: .14f} {y[i]: .14f}  {z[i]: .14f}  {atom[i]}")
        with open(
            os.path.join(os.path.split(path)[0], "coord"), "w", newline=None
        ) as coord:
            coord.write("$coord\n")
            for line in coordxyz:
                coord.write(line + "\n")
            coord.write("$end\n")


def write_trj(
    results, cwd, outpath, optfolder, nat, attribute, overwrite=False, *args, **kwargs
):
    """
    Write trajectory (multiple xyz geometries) to file.
    """
    if overwrite and os.path.isfile(outpath):
        os.remove(outpath)
    for key, value in kwargs.items():
        if key == "rrho":
            rrho = value
        elif key == "energy":
            energy = value
    try:
        rrho
    except NameError:
        rrho = None
    try:
        energy
    except NameError:
        energy = None
    try:
        with open(outpath, "a", encoding=CODING, newline=None) as out:
            for conf in results:
                conf_xyz, nat = t2x(os.path.join(cwd, "CONF" + str(conf.id), optfolder))
                ### coordinates in xyz
                out.write("  {}\n".format(nat))
                xtbfree = conf.calc_free_energy(
                    e=energy, solv=None, rrho=rrho, out=True
                )
                if xtbfree is not None:
                    xtbfree = f"{xtbfree:20.8f}"
                out.write(
                    f"G(CENSO)= {getattr(conf, attribute):20.8f}"
                    f"  G(xTB)= {xtbfree}"
                    f"        !CONF{str(conf.id)}\n"
                )
                for line in conf_xyz:
                    out.write(line + "\n")
    except (FileExistsError, ValueError):
        print(f"{'WARNING:':{WARNLEN}}Could not write trajectory: "
              f"{last_folders(outpath, 1)}."
              )


def check_for_float(line):
    """ Go through line and check for float, return first float"""
    elements = line.strip().split()
    value = None
    for element in elements:
        try:
            value = float(element)
            found = True
        except ValueError:
            found = False
            value = None
        if found:
            break
    return value


def last_folders(path, number=1):
    """
    Return string of last folder or last two folders of path, depending on number
    """
    if number not in (1, 2, 3):
        number = 1
    if number == 1:
        folder = os.path.basename(path)
    if number == 2:
        folder = os.path.join(
            os.path.basename(os.path.dirname(path)), os.path.basename(path)
        )
    if number == 3:
        basename = os.path.basename(path)
        dirname = os.path.basename(os.path.dirname(path))
        predirname = os.path.basename(os.path.split(os.path.split(path)[0])[0])
        folder = os.path.join(predirname, dirname, basename)
    return folder


def get_energy_from_ensemble(path, config, conformers):
    """
    Get energies from the ensemble inputfile and assign xtb_energy and
    rel_xtb_energy
    """
    with open(path, "r", encoding=CODING, newline=None) as inp:
        data = inp.readlines()
    if config.maxconf * (config.nat + 2) > len(data):
        print(
            f"{'ERROR:':{WARNLEN}}Either the number of conformers ({config.nconf}) "
            f"or the number of atoms ({config.nat}) is wrong!"
        )
    # calc energy and rel energy:
    e = {}
    conformers.sort(key=lambda x: int(x.id))
    for conf in conformers:
        e[conf.id] = check_for_float(data[(conf.id - 1) * (config.nat + 2) + 1])
    try:
        lowest = float(min([i for i in e.values() if i is not None]))
    except (ValueError, TypeError):
        print(f"{'WARNING:':{WARNLEN}}Can't calculate rel_xtb_energy!")
        return
    for conf in conformers:
        try:
            conf.xtb_energy = e[conf.id]
            conf.rel_xtb_energy = (e[conf.id] - lowest) * AU2KCAL
            # print(f"CONF{conf.id} {conf.xtb_energy} {conf.rel_xtb_energy}")
        except (ValueError, TypeError) as e:
            print(e)
    return conformers


def ensemble2coord(config, foldername, conflist, store_confs, save_errors):
    """
    read ensemble file: e.g. 'crest_conformers.xyz' and write coord files into
    designated folders

    - path [abs path] to ensemble file
    - nat  [int] number of atoms in molecule
    - nconf [int] number of considered conformers
    - cwd [path] path to current working directory
    - foldername [str] name of folder into which the coord file is to be written
    - conflist [list with conf object] all conf objects

    returns list with conformer objects
    """
    if not os.path.isfile(config.ensemblepath):
        print(f"ERROR: File {os.path.basename(config.ensemblepath)} does not exist!")
    with open(config.ensemblepath, "r", encoding=CODING, newline=None) as inp:
        data = inp.readlines()
    if config.maxconf * (config.nat + 2) > len(data):
        print(
            f"ERROR: Either the number of conformers ({config.nconf}) "
            f"or the number of atoms ({config.nat}) is wrong!"
        )
    for conf in conflist:
        i = conf.id
        atom = []
        x = []
        y = []
        z = []
        start = (i - 1) * (config.nat + 2) + 2
        end = i * (config.nat + 2)
        for line in data[start:end]:
            atom.append(str(line.split()[0].lower()))
            x.append(float(line.split()[1]) / BOHR2ANG)
            y.append(float(line.split()[2]) / BOHR2ANG)
            z.append(float(line.split()[3]) / BOHR2ANG)
        coordxyz = []
        for j in range(len(x)):
            coordxyz.append(f"{x[j]: .14f} {y[j]: .14f}  {z[j]: .14f}  {atom[j]}")
        outpath = os.path.join(config.cwd, "CONF" + str(conf.id), foldername, "coord")
        if not os.path.isfile(outpath):
            # print(f"Write new coord file in {last_folders(outpath)}")
            with open(outpath, "w", newline=None) as coord:
                coord.write("$coord\n")
                for line in coordxyz:
                    coord.write(line + "\n")
                coord.write("$end")
    return conflist, store_confs, save_errors


def splitting(item):
    """
    Used in move_recursively.

    """
    try:
        return int(item.rsplit(".", 1)[1])
    except ValueError:
        return 0


def move_recursively(path, filename):
    """
    Check if file or file.x exists and move them to file.x+1 ignores e.g.
    file.save
    """
    files = [
        f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
    ]  # list of all files in directory
    newfiles = []  # list of all files in directory that contain filename and '.'
    for item in files:
        if filename + "." in item:
            newfiles.append(item)
    newfiles.sort(key=splitting, reverse=True)
    for item in newfiles:
        try:
            data = item.rsplit(".", 1)  # splits only at last '.'
            int(data[1])
        except ValueError:
            continue
        tmp_from = os.path.join(path, item)
        newfilename = str(data[0]) + "." + str(int(data[1]) + 1)
        tmp_to = os.path.join(path, newfilename)
        # print("Backing up {} to {}.".format(item, newfilename))
        shutil.move(tmp_from, tmp_to)

    if filename in files:
        print("Backing up {} to {}.".format(filename, filename + ".1"))
        shutil.move(os.path.join(path, filename), os.path.join(path, filename + ".1"))


def calc_boltzmannweights(confs, property, T):
    """
    Calculate Boltzmannweights:
    - confs [list] list with conformer objects
    - property [str] e.g. free_energy of conformer
    - T [float] temperature at which the Boltzmann weight has to be evaluated

    returns confs
    """
    if len(confs) == 1:
        confs[0].bm_weight = 1.0
        return confs
    try:
        T = float(T)
    except ValueError:
        T = 298.15  # K
        print(f"{'WARNING:':{WARNLEN}}Temperature can not be converted and is therfore set to T = {T} K.")
    if T == 0:
        T += 0.00001  # avoid division by zero
    try:
        minfree = min(
            [
                getattr(conf, property, None)
                for conf in confs
                if getattr(conf, property, None) is not None
            ]
        )
    except ValueError:
        print(f"{'ERROR:':{WARNLEN}}Boltzmann weight can not be calculated!")
    bsum = 0.0
    for item in confs:
        bsum += getattr(item, "gi", 1.0) * math.exp(
            -((item.free_energy - minfree) * AU2J) / (KB * T)
        )
    for item in confs:
        item.bm_weight = (
            getattr(item, "gi", 1.0)
            * math.exp(-((item.free_energy - minfree) * AU2J) / (KB * T))
            / bsum
        )
    return confs


def new_folders(cwd, conflist, foldername, save_errors, store_confs, silent=False):
    """ 
    create folders for all conformers in conflist
    """

    for conf in conflist:
        tmp_dir = os.path.join(cwd, "CONF" + str(conf.id), foldername)
        try:
            mkdir_p(tmp_dir)
        except Exception as e:
            print(e)
            if not os.path.isdir(tmp_dir):
                print(f"{'ERROR:':{WARNLEN}}Could not create folder for CONF{conf.id}!")
                print(f"{'ERROR:':{WARNLEN}}CONF{conf.id} is removed, because IO failed!")
                save_errors.append(f"{'ERROR:':{WARNLEN}}CONF{conf.id} was removed, because IO failed!")
                store_confs.append(conflist.pop(conflist.index(conf)))
    if not silent:
        print("Constructed folders!")
    return save_errors, store_confs, conflist


def check_for_folder(path, conflist, foldername, debug=False):
    """
    Check if folders exist (of conformers calculated in previous run)
    """
    error_logical = False
    for i in conflist:
        tmp_dir = os.path.join(path, "CONF" + str(i), foldername)
        if not os.path.exists(tmp_dir):
            print(
                f"{'ERROR:':{WARNLEN}}directory of {last_folders(tmp_dir, 2)} does not exist, although "
                "it was calculated before!"
            )
            error_logical = True
    if error_logical and not debug:
        print(f"{'ERROR:':{WARNLEN}}One or multiple directories are missing.\n")
    return error_logical


def do_md5(path):
    """
    Calculate md5 of file to identifly if restart happend on the same file!
    Input is buffered into smaller sizes to ease on memory consumption.
    """
    BUF_SIZE = 65536
    md5 = hashlib.md5()
    if os.path.isfile(path):
        with open(path, "rb") as f:
            while True:
                data = f.read(BUF_SIZE)
                if not data:
                    break
                md5.update(data)
        return md5.hexdigest()
    else:
        raise FileNotFoundError


def rank_simple(vector):
    """
    needed to rank vectors
    """
    return sorted(range(len(vector)), key=vector.__getitem__)


def rankdata(a):
    """
    rank vectors like in numpy
    """
    n = len(a)
    ivec = rank_simple(a)
    svec = [a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0] * n
    for i in range(n):
        sumranks += i
        dupcount += 1
        if i == n - 1 or svec[i] != svec[i + 1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i - dupcount + 1, i + 1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newarray


def pearson(A, B):
    """
    Calculate pearson correlation coefficient
    """
    if len(A) != len(B):
        print("ERROR in PEARSON lists are not of equal length!")
    n = float(len(A))
    muA = sum(A) / n
    muB = sum(B) / n
    diffA = map(lambda x: x - muA, A)
    diffB = map(lambda x: x - muB, B)
    stdA = math.sqrt((1 / (n - 1)) * sum([d * d for d in diffA]))
    stdB = math.sqrt((1 / (n - 1)) * sum([d * d for d in diffB]))
    try:
        return (sum([A[i] * B[i] for i in range(int(n))]) - n * muA * muB) / (
            (n - 1) * stdA * stdB
        )
    except ZeroDivisionError as e:
        print(f"{'WARNING:':{WARNLEN}}{e}")
        return 0.0


def spearman(A, B):
    """
    Calculate spearman correlation coefficient
    """
    return pearson(rankdata(A), rankdata(B))


def printout(
    outputpath,
    columncall,
    columnheader,
    columndescription,
    columnformat,
    calculate,
    minfree,
    columndescription2=None,
):
    """
    Create printout which is printed to stdout and file.
    """
    if columndescription2 is None:
        columndescription2 = []
    calculate.sort(key=lambda x: int(x.id))
    if not any(
        [
            len(i) == len(columncall)
            for i in (columnheader, columndescription, columnformat)
        ]
    ):
        print("Lists of uneqal length!")
    collength = []
    columnheaderprint = []
    columndescriptionprint = []
    columndescriptionprint2 = []
    if not columndescription2:
        columndescription2 = ["" for _ in range(len(columncall))]
    # split on "["   eg. COSMORS[B97-3c/def2-TZVP]

    for i in range(len(columndescription)):
        if "[" in columndescription[i] and columndescription[i] not in (
            "[Eh]",
            "[kcal/mol]",
            "[a.u.]",
        ):
            columndescription2[i] = "[" + str(columndescription[i]).split("[")[1]
            columndescription[i] = str(columndescription[i]).split("[")[0]
    try:
        for j in range(len(columncall)):
            if columnformat[j]:
                collength.append(
                    max(
                        [
                            len(str(f"{i:{columnformat[j][0]}.{columnformat[j][1]}f}"))
                            for i in map(columncall[j], calculate)
                        ]
                    )
                )
            else:
                collength.append(max([len(i) for i in map(columncall[j], calculate)]))
            if (
                max(
                    len(i)
                    for i in [
                        columndescription[j],
                        columnheader[j],
                        columndescription2[j],
                    ]
                )
                > collength[j]
            ):
                collength[j] = max(
                    len(i)
                    for i in [
                        columndescription[j],
                        columnheader[j],
                        columndescription2[j],
                    ]
                )
    except (ValueError, TypeError) as e:
        print(f"\n\nERRROR {e}")
        for j in range(len(columncall)):
            collength.append(12)

    for i in range(len(columncall)):
        columnheaderprint.append(f"{columnheader[i]:>{collength[i]}}")
        columndescriptionprint.append(f"{columndescription[i]:>{collength[i]}}")
        if columndescription2:
            columndescriptionprint2.append(f"{columndescription2[i]:>{collength[i]}}")
    with open(outputpath, "w", newline=None) as out:
        line = " ".join(columnheaderprint)
        print(line)
        out.write(line + "\n")
        line = " ".join(columndescriptionprint)
        print(line)
        out.write(line + "\n")
        if columndescription2:
            line = " ".join(columndescriptionprint2)
            print(line)
            out.write(line + "\n")
        for conf in calculate:
            columncallprint = []
            for i in range(len(columncall)):
                if columnformat[i]:
                    columncallprint.append(
                        f"{columncall[i](conf):{collength[i]}.{columnformat[i][1]}f}"
                    )
                else:
                    columncallprint.append(f"{columncall[i](conf):{collength[i]}}")
            if conf.free_energy != minfree:
                line = " ".join(columncallprint)
                print(line)
                out.write(line + "\n")
            else:
                line = " ".join(columncallprint + [f"    <------"])
                print(line)
                out.write(line + "\n")


def crest_routine(config, conformers, func, store_confs, prev_calculated=None):
    """
    check if two conformers are rotamers of each other,
    this check is always performed, but removing conformers depends on crestcheck
    returns conformers
    returns store_confs
    returns prev_calculated
    """
    if prev_calculated is None:
        prev_calculated = []
    dirn = "conformer_rotamer_check"  ### directory name
    fn = "conformers.xyz"  ### file name

    print("\nChecking for identical structures in ensemble with CREGEN!\n")
    # create folder for comparison
    if not os.path.isdir(os.path.join(config.cwd, dirn)):
        mkdir_p(os.path.join(config.cwd, dirn))
    # delete conformers.xyz file if it already exists
    if os.path.isfile(os.path.join(config.cwd, dirn, fn)):
        os.remove(os.path.join(config.cwd, dirn, fn))
    # delete coord file if exists
    if os.path.isfile(os.path.join(config.cwd, dirn, "coord")):
        os.remove(os.path.join(config.cwd, dirn, "coord"))

    allconfs = deepcopy(conformers)
    allconfs.extend(deepcopy(prev_calculated))

    ### sort conformers according to energy of optimization
    allconfs.sort(key=lambda conf: float(getattr(conf, "optimization_info")["energy"]))
    # write coord:
    try:
        shutil.copy(
            os.path.join(config.cwd, "CONF" + str(allconfs[0].id), func, "coord"),
            os.path.join(config.cwd, dirn, "coord"),
        )
    except Exception as e:
        print(f"{'ERROR:':{WARNLEN}}{e}")

    # write conformers.xyz file
    with open(
        os.path.join(config.cwd, dirn, fn), "w", encoding=CODING, newline=None
    ) as out:
        for conf in allconfs:
            conf_xyz, nat = t2x(os.path.join(config.cwd, "CONF" + str(conf.id), func))
            out.write("  {}\n".format(nat))  ### number of atoms
            out.write(
                "{:20.8f}        !{}\n".format(
                    getattr(conf, "optimization_info")["energy"], "CONF" + str(conf.id)
                )
            )
            for line in conf_xyz:
                out.write(line + "\n")
        for conf in allconfs:
            conf_xyz, nat = t2x(os.path.join(config.cwd, "CONF" + str(conf.id), func))
            out.write("  {}\n".format(nat))  ### number of atoms
            out.write(
                "{:20.8f}        !{}\n".format(
                    getattr(conf, "optimization_info")["energy"], "CONF" + str(conf.id)
                )
            )
            for line in conf_xyz:
                out.write(line + "\n")
    time.sleep(0.01)

    crestcall = [
        config.external_paths["crestpath"],
        "coord",
        "-cregen",
        fn,
        "-ethr",
        "0.15",
        "-rthr",
        "0.175",
        "-bthr",
        "0.03",
        "-enso",
    ]

    with open(
        os.path.join(config.cwd, dirn, "crest.out"), "w", newline=None, encoding=CODING
    ) as outputfile:
        subprocess.call(
            crestcall,
            shell=False,
            stdin=None,
            stderr=subprocess.STDOUT,
            universal_newlines=False,
            cwd=os.path.join(config.cwd, dirn),
            stdout=outputfile,
            env=ENVIRON,
        )
        time.sleep(0.05)
        try:
            with open(
                os.path.join(config.cwd, dirn, "enso.tags"),
                "r",
                encoding=CODING,
                newline=None,
            ) as inp:
                store = inp.readlines()
        except (Exception) as e:
            print(f"{'ERROR:':{WARNLEN}}{e}")
            print(f"{'ERROR:':{WARNLEN}}output file (enso.tags) of CREST routine does not exist!")
        keep = []
    if config.crestcheck:
        try:
            for line in store:
                keep.append(line.split()[1][1:])
            for conf in list(conformers):
                if "CONF" + str(conf.id) not in keep:
                    conf.optimization_info["info"] = "calculated"
                    conf.optimization_info["cregen_sort"] = "removed"
                    print(
                        f"!!!! Removing CONF{conf.id} because it is sorted "
                        "out by CREGEN."
                    )
                    store_confs.append(conformers.pop(conformers.index(conf)))
            for conf in list(prev_calculated):
                if "CONF" + str(conf.id) not in keep:
                    conf.optimization_info["info"] = "calculated"
                    conf.optimization_info["cregen_sort"] = "removed"
                    print(
                        f"!!!! Removing CONF{conf.id} because it is sorted "
                        "out by CREGEN."
                    )
                    store_confs.append(prev_calculated.pop(prev_calculated.index(conf)))
        except (NameError, Exception) as e:
            print(f"{'ERROR:':{WARNLEN}}{e}")
    return conformers, prev_calculated, store_confs


def format_line(key, value, options, optionlength=70, dist_to_options=30):
    """
    used in print_parameters
    """
    # limit printout of possibilities
    if len(str(options)) > optionlength:
        length = 0
        reduced = []
        for item in options:
            length += len(item) + 2
            if length < optionlength:
                reduced.append(item)
        reduced.append("...")
        options = reduced
    line = "{}: {:{digits}} # {} \n".format(
        key, str(value), options, digits=dist_to_options - len(key)
    )
    return line


def check_tasks(results, check=False, thresh=0.25):
    """
    Check if too many tasks failed and exit if so!
    """
    # Check if preparation failed too often:
    counter = 0
    for item in results:
        if not item.job["success"]:
            counter += 1
    try:
        fail_rate = float(counter) / float(len(results))
    except ZeroDivisionError:
        print(f"{'ERROR:':{WARNLEN}}Too many calculations failed!" "\nGoing to exit!")
        sys.exit(1)
    if fail_rate >= thresh and check:
        print(
            f"{'ERROR:':{WARNLEN}}{fail_rate*100} % of the calculations failed!" "\nGoing to exit!"
        )
        sys.exit(1)
    elif fail_rate >= thresh:
        print(f"{'WARNING:':{WARNLEN}}{fail_rate*100} % of the calculations failed!")


def isclose(value_a, value_b, rel_tol=1e-9, abs_tol=0.0):
    """
    Replace function if not available from math module (exists since python 3.5)
    """
    return abs(value_a - value_b) <= max(
        rel_tol * max(abs(value_a), abs(value_b)), abs_tol
    )


def calc_std_dev(data):
    """
    Calculate standard deviation
    """
    n = len(data)
    if len(data) != 0:
        mean = sum(data) / n
        variance = sum([(x - mean) ** 2 for x in data]) / (n - 1)
        std_dev = math.sqrt(variance)
    else:
        std_dev = 0.0
    return std_dev


def calc_weighted_std_dev(data, weights=None):
    """
    Calculate standard deviation
    """
    if weights is None:
        weights = []
    n = len(data)
    if n == 0:
        return 0.0
    if not weights or len(weights) < n:
        weights = [1.0 for _ in range(n)]
    w_mean = sum([data[i] * weights[i] for i in range(n)]) / sum(weights)
    m = 0
    for i in weights:
        if i != 0.0:
            m += 1
    variance = sum([weights[i] * (data[i] - w_mean) ** 2 for i in range(n)]) / (
        (m - 1) * sum(weights) / m
    )
    std_dev = math.sqrt(variance)
    return std_dev

def print_errors(line, save_errors):
    """print line and append to list save_errors"""
    print(line)
    try:
        save_errors.append(line)
    except Exception:
        pass
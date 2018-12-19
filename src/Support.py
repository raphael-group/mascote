import sys
import math
import random
import datetime

import Support as sp

def bins(size, step):
    return zip(xrange(0, size-step, step), xrange(step, size, step)) + [(int(math.floor((size - 1) / step)) * step, size)]

def cumsum(lis):
    total = 0
    for x in lis:
        total += x
        yield total

def cumsumbins(lis):
    total = 0
    for x in lis:
        total += x[1] - x[0]
        yield total

def startendbins(bins, start, size):
    cumsum = list(cumsumbins(bins))
    startbin = -1

    for index, item in enumerate(cumsum):
        if item > start:
            startbin = index
            break

    if startbin == -1:
        raise ValueError('The start bin must be inside the given sequence of bins!')

    endbin = -1
    for index, item in enumerate(cumsum[startbin:]):
        if item > start+size:
            endbin = index
            break

    return startbin, startbin+endbin+1 if endbin != -1 else len(bins)

def wchoice(elements, weights):
    assert(len(elements) == len(weights))
    total = int(sum(weights))
    pick = random.randint(1, total)
    upto = 1
    for el, w in zip(elements, weights):
        if upto + w >= pick:
            return el
        upto += w
    assert(False)


def parseSNPList(filename):
    if filename != None:
        snplist = {}
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip().split()
                if len(line) > 2:
                    pos = int(line[1])
                    if not line[0] in snplist:
                        snplist[line[0]] = {}
                    if not pos in snplist[line[0]]:
                        snplist[line[0]][pos] = set(char.upper() for rest in line[2:] for char in rest if char.upper() in ('A', 'C', 'G', 'T'))
        return snplist
    else:
        return None


def log(msg, level=None, lock=None):
    timestamp = '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
    if level == "STEP":
        if lock is None:
            sys.stderr.write("{}{}[{}]{}{}".format(bcolors.BOLD, bcolors.HEADER, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}{}[{}]{}{}".format(bcolors.BOLD, bcolors.HEADER, timestamp, msg, bcolors.ENDC))
    elif level == "INFO":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.BBLUE, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.OKGREEN, timestamp, msg, bcolors.ENDC))
    elif level == "WARN":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.WARNING, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.WARNING, timestamp, msg, bcolors.ENDC))
    elif level == "PROGRESS":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.OKGREEN, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.BBLUE, timestamp, msg, bcolors.ENDC))
    elif level == "ERROR":
        if lock is None:
            sys.stderr.write("{}[{}]{}{}".format(bcolors.FAIL, timestamp, msg, bcolors.ENDC))
        else:
            with lock: sys.stderr.write("{}[{}]{}{}".format(bcolors.FAIL, timestamp, msg, bcolors.ENDC))
    else:
        if lock is None:
            sys.stderr.write("{}".format(msg))
        else:
            with lock: sys.stderr.write("{}".format(msg))


def error(msg):
    return "{}{}{}".format(bcolors.FAIL, msg, bcolors.ENDC)


def close(msg):
    log(msg=msg, level="WARN")
    sys.exit(0)


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    BBLUE = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def basesize(isize):
    size = 0
    try:
        if isize[-2:] == "kb":
            size = int(isize[:-2]) * 1000
        elif isize[-2:] == "Mb":
            size = int(isize[:-2]) * 1000000
        else:
            size = int(isize)
        return size
    except:
        raise ValueError(sp.error("Size must be a number, optionally ending with either \"kb\" or \"Mb\"!"))

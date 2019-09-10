# -*- coding:utf-8 -*-
import ctypes
import os

import util


# port directly from C code just to be sure
def eed_c(hyp: list,
          ref: list,
          alpha: float,
          deletion: float,
          insertion: float,
          substitution: float,
          rho: float,
          norm: int):
    # float EED(const std::vector<unsigned long long> hyp,
    #     const std::vector<unsigned long long> ref,
    #     const float alpha,
    #     const float deletion,
    #     const float insertion,
    #     const float substitution,
    #     const float rho,
    #     const int norm){

    # std::vector<int> lj(hyp.size() + 1, 0); //Coverage Tracker
    lj = [0] * (len(hyp) + 1)

    # std::vector<float> row(hyp.size() + 1, 1);
    row = [1.0] * (len(hyp) + 1)

    # row[0] = float(0); //CDER initialisation 0,0 = 0 rest 1
    row[0] = 0.0

    # std::vector<float> nextRow(hyp.size() + 1, std::numeric_limits<float>::max());
    next_row = [float('inf')] * (len(hyp) + 1)

    # vector<vector<float>> m(ref.size()+1);
    m = [[] for _ in range(len(ref) + 1)]

    # int minInd; //Min index
    min_ind = float('nan')

    # int minimum; //Min value
    minimum = float('nan')

    # for(int w = 1; w < ref.size() + 1; ++w){
    for w in range(1, len(ref) + 1):

        # for(int i = 0; i < hyp.size() + 1; ++i){
        for i in range(len(hyp) + 1):

            # if(i > 0){
            if i > 0:
                # nextRow[i] = std::min({nextRow[i-1] + deletion,
                next_row[i] = min([next_row[i - 1] + deletion,
                                   # row[i-1] + substitution*(ref[w-1] != hyp[i-1]),
                                   row[i - 1] + substitution * (ref[w - 1] != hyp[i - 1]),
                                   # row[i]+ insertion});
                                   row[i] + insertion])
            # }
            # else{
            else:
                # nextRow[i] = row[i]+ 1.0;
                next_row[i] = row[i] + 1.0
            # }

        # }

        # minimum = nextRow[0];
        minimum = next_row[0]

        # minInd = 0;
        min_ind = 0

        # for(int i = 1; i < nextRow.size(); ++i){
        for i in range(1, len(next_row)):

            # if(nextRow[i] < minimum){
            if next_row[i] < minimum:
                # minInd = i;
                min_ind = i

                # minimum = nextRow[minInd];
                minimum = next_row[min_ind]
            # }
        # }
        assert min_ind == next_row.index(min(next_row))

        # lj[minInd] = lj[minInd] + 1;
        lj[min_ind] = lj[min_ind] + 1

        # //Long Jumps 32 is the magic number for white spaces
        # if(ref[w-1] == 32){
        if ref[w - 1].isspace():

            # float longJump = alpha + nextRow[minInd];
            long_jump = alpha + next_row[min_ind]

            # for(int i = 0; i < nextRow.size(); ++i){
            for i in range(len(next_row)):

                # if(nextRow[i] > longJump){
                if next_row[i] > long_jump:
                    # nextRow[i] = longJump;
                    next_row[i] = long_jump
                # }
            # }
        # }

        # m[w] = row;
        m[w] = row

        # row = nextRow;
        row = next_row

        # nextRow.assign(nextRow.size() ,std::numeric_limits<float>::max());
        next_row = [float('inf')] * len(next_row)
    # }

    # float coverage = 0;
    coverage = 0

    # for(int i = 1; i < lj.size(); ++i){
    for i in range(1, len(lj)):

        # if(lj[i] > 1){
        if lj[i] > 1:
            # coverage += lj[i];
            coverage += lj[i]
            # }
        # }

    assert coverage == sum([x if x > 1 else 0 for x in lj])

    # float errors = row[row.size()-1];
    errors = row[-1]

    # return (errors + rho*coverage)/(norm + rho*coverage);
    return (errors + rho * coverage) / (norm + rho * coverage)
    # }


def eed(hyp: list, ref: list):
    hyp.insert(0, " ")
    hyp.append(" ")
    ref.insert(0, " ")
    ref.append(" ")

    # hyp_c = [bytes_to_int(x.encode('utf-8')) for x in hyp]
    # ref_c = [bytes_to_int(x.encode('utf-8')) for x in ref]

    alpha = 2.0
    deletion = 0.2
    insertion = 1.0
    substitution = 1.0
    rho = 0.3
    norm = len(ref)

    result = eed_c(hyp, ref, alpha, deletion, insertion, substitution, rho, norm)

    return min(1.0, result)


# Python wrpaper for the C++ EED implementation
def eed_orig(hyp, ref):
    _eed = ctypes.CDLL(os.path.dirname(os.path.abspath(__file__)) + '/libEED.so')
    # print(os.path.dirname(os.path.abspath(__file__)))
    _eed.wrapper.restype = ctypes.c_float
    hyp.insert(0, " ")
    hyp.append(" ")
    ref.insert(0, " ")
    ref.append(" ")
    hyp_c = (ctypes.c_ulonglong * len(hyp))()
    hyp_c[:] = [bytes_to_int(x.encode('utf-8')) for x in hyp]
    ref_c = (ctypes.c_ulonglong * len(ref))()
    ref_c[:] = [bytes_to_int(x.encode('utf-8')) for x in ref]
    alpha = 2.0
    deletion = 0.2
    insertion = 1.0
    substitution = 1.0
    rho = 0.3
    norm = len(ref_c)
    result = _eed.wrapper(hyp_c, ref_c, len(hyp_c), len(ref_c),
                          ctypes.c_float(alpha), ctypes.c_float(deletion),
                          ctypes.c_float(insertion), ctypes.c_float(substitution), ctypes.c_float(rho),
                          norm)
    return min(1.0, result)


def bytes_to_int(bytes):
    result = 0
    for b in bytes:
        result = result * 256 + int(b)
    return result


# Distance measure used for substitutions/identity operation
def distance(refWord, hypWord):
    if refWord == hypWord:
        return 0
    else:
        return 1


# Python Implementation of EED, ~30x slower than the C++ one
def eed_python(hyp, ref):
    hyp.insert(0, " ")
    hyp.append(" ")
    ref.insert(0, " ")
    ref.append(" ")

    alpha = 2.0
    deletion = 0.2
    # substitutions are implemented via the distance function
    insertion = 1.0
    rho = 0.3

    lj = [0] * (len(hyp) + 1)
    row = [1] * (len(hyp) + 1)  # row[i] stores cost of cheapest path from (0,0) to (i,l) in CDER aligment grid.
    row[0] = 0  # CDER initialisation 0,0 = 0 rest 1
    nextRow = [float('inf')] * (len(hyp) + 1)
    for w in range(1, len(ref) + 1):
        for i in range(0, len(hyp) + 1):
            if i > 0:
                nextRow[i] = min(nextRow[i - 1] + deletion, row[i - 1] + distance(ref[w - 1], hyp[i - 1]),
                                 row[i] + insertion)
            else:
                nextRow[i] = row[i] + 1.0
        minInd = nextRow.index(min(nextRow))
        lj[minInd] += 1
        # Long Jumps
        if ref[w - 1] == " ":
            longJump = alpha + nextRow[minInd]
            nextRow = [x if x < longJump else longJump for x in nextRow]
        row = nextRow
        nextRow = [float('inf')] * (len(hyp) + 1)

    coverage = rho * sum([x if x > 1 else 0 for x in lj])

    return min(1, (row[-1] + coverage) / (float(len(ref)) + coverage))


# Provides System scoring with full preprocessing in the case where EED is used as an import
def score(hypIn, refIn):
    import codecs
    hyp = [util.preprocess(x) for x in open(hypIn, 'r', 'utf-8').readlines()]
    ref = [util.preprocess(x) for x in open(refIn, 'r', 'utf-8').readlines()]
    scores = []
    for (h, r) in zip(hyp, ref):
        h, r = list(h), list(r)
        score = eed(h, r)
        scores.append(score)
    return sum(scores) / len(scores)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Extended Edit Distance Metric for Machine Translation')
    parser.add_argument('-ref', '--reference', help='Reference file', required=True)
    parser.add_argument('-hyp', '--hypothesis', help='Input(test) file', required=True)
    parser.add_argument('-v', '--verbose', help='Show scores of each sentence.',
                        action='store_true', default=False)
    return parser.parse_args()


def main():
    import sys
    import codecs
    args = parse_args()
    hlines = [util.preprocess(x, "en") for x in codecs.open(args.hypothesis, 'r', 'utf-8').readlines()]
    rlines = [util.preprocess(x, "en") for x in codecs.open(args.reference, 'r', 'utf-8').readlines()]
    if len(hlines) != len(rlines):
        print("Error: input file has {0} lines, but reference has {1} lines.".format(len(hlines), len(rlines)))
        sys.exit(1)
    scores = []
    for lineno, (hline, rline) in enumerate(zip(hlines, rlines), start=1):
        rline, hline = list(rline), list(hline)
        score = eed(hline, rline)
        scores.append(score)
        if args.verbose:
            print("Sentence {0}: {1:.4f}".format(lineno, score))

    average = sum(scores) / len(scores)

    print("System Score={0:.4f}".format(average))
    sys.exit(0)


if __name__ == '__main__':
    # main()
    hline = list('The relationship between Obama and Netanyahu has been strained for years.')
    rline = list('Relations between Obama and Netanyahu have been strained for years.')
    print(eed(hline, rline))

    hline = list('The relationship between Obama and Netanyahu has been strained for years.')
    rline = list('Relations between Obama and Netanyahu have been strained for years.')
    print(eed(rline, hline))

    hline = list('The relationship between Obama and Netanyahu has been strained for years.')
    rline = list('Relations between Obama and Netanyahu have been strained for years.')
    print(eed_python(hline, rline))

    hline = list('The relationship between Obama and Netanyahu has been strained for years.')
    rline = list('Relations between Obama and Netanyahu have been strained for years.')
    print(eed_python(rline, hline))

import util


def eed(hyp, ref, deletion=0.2, insertion=1.0, substitution=1.0, jump=2.0, rho=0.3):
    """
    Extended Edit Distance (https://www.aclweb.org/anthology/W19-5359)
    Calculates the (non-symmetric) EED from hypothesis to reference string
    Note: EED builds on  CDER (https://www.aclweb.org/anthology/E06-1031)

    with default settings, the algorithm seems to:
    - skip the first char in each hyp word, even if it matches the ref
    - repeat the last hyp character until it sees a space in the ref
    - but it doesn't skip the first char if there's no jump
    - the c code stores the entire DP table for no reason

    other notes:
    - skipping the first char seems to be solved by setting jump cost to 0.9
    - setting jump cost higher than insertion cost tends to insert space as first word char
    - setting jump cost higher than substitution cost tends to insert random chars instead
    - rho changes nothing since it's computed at the end
    - attempts to integrate rho into the cost function failed miserably

    :param hyp: hypothesis sentence
    :type hyp: str
    :param ref: reference sentence
    :type ref: str
    :param deletion: deletion cost
    :type deletion: float
    :param insertion: insertion cost
    :type insertion: float
    :param substitution: substitution cost
    :type substitution: float
    :param jump: jump cost
    :type jump: float
    :param rho: coverage cost weight
    :type rho: float
    :return: EED score of hyp given ref (not symmetric)
    :rtype: float
    """

    # start and end with whitespace to facilitate jumps to front/end
    # works better when you use the most common whitespace (usually spaces)
    # this step is defined as part of the algorithm
    hyp = ' ' + hyp + ' '
    ref = ' ' + ref + ' '

    # only for debugging
    debug_table = []  # store full DP table
    debug_str = []  # matched string
    debug_cost = []  # costs
    debug_idx = []  # indices of matched string

    # coverage: count how many times each char is visited
    visit_coverage = [-1.0] * (len(hyp) + 1)

    # the i-th row stores cost of cheapest path from (0,0) to (i,l) in CDER alignment grid
    row = [0.0] + [1.0] * len(hyp)  # CDER initial row

    for ref_idx, ref_char in enumerate(ref):
        next_row = [float('inf')] * (len(hyp) + 1)

        # add 1 to the cost per row (same as edit distance)
        next_row[0] = row[0] + 1.0

        # do the normal edit distance calculation for the hyp sentence
        for hyp_idx, hyp_char in enumerate(hyp):
            next_row[hyp_idx + 1] = min([next_row[hyp_idx] + deletion,
                                         row[hyp_idx] + (0.0 if ref_char == hyp_char else substitution),
                                         row[hyp_idx + 1] + insertion])

        # this is the next char to be visited according to the EED algo
        min_cost, min_cost_idx = min((cost, idx) for idx, cost in enumerate(next_row))

        # increment the visit count
        visit_coverage[min_cost_idx] += 1.0

        # long jump allowed only if ref char is whitespace
        # the original algo only checks if ord(ref_char) == 32 (ascii space)
        if ref_char.isspace():
            long_jump_cost = jump + min_cost
            next_row = [min(x, long_jump_cost) for x in next_row]

        # for debug
        debug_table.append(row)
        debug_str.append(hyp[min_cost_idx - 1])
        debug_cost.append(min_cost)
        debug_idx.append(min_cost_idx)

        row = next_row

    # overall error == final cell of final row
    errors = row[-1]
    weighted_coverage = rho * sum(x if x >= 0.0 else 1.0 for x in visit_coverage)
    # weighted_coverage = rho * sum(1.0 for x in visit_coverage if x != 1.0)  # shouldn't this be the correct impl
    result = (errors + weighted_coverage) / (len(ref) + weighted_coverage)

    # debug
    debug_table.append(row)
    # print(hyp)
    # print(''.join(debug_str))
    # print(ref)
    # print(list(zip(debug_str, debug_idx)))
    for row in debug_table:
        # print(row)
        pass

    return min(1.0, result)


# Distance measure used for substitutions/identity operation
def distance(refWord, hypWord):
    return 0 if refWord == hypWord else 1


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
    hyp = [util.preprocess(x) for x in open(hypIn, mode='rt', encoding='utf-8').readlines()]
    ref = [util.preprocess(x) for x in open(refIn, mode='rt', encoding='utf-8').readlines()]
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
    hyp_sent = 'The relationship between Obama and Netanyahu has been strained for years.'
    ref_sent = 'Relations between Obama and Netanyahu have been strained for years.'

    print(eed_python(list(hyp_sent), list(ref_sent)))
    print(eed_python(list(ref_sent), list(hyp_sent)))

    print(eed(hyp_sent, ref_sent))
    print(eed(ref_sent, hyp_sent))

    hyp_sent = 'super calloused fragile mystic hexed by haliotosis'
    ref_sent = 'supercalifragilisticexpialidocious'

    print(eed_python(list(hyp_sent), list(ref_sent)))
    print(eed_python(list(ref_sent), list(hyp_sent)))

    print(eed(hyp_sent, ref_sent))
    print(eed(ref_sent, hyp_sent))

    hyp_sent = 'qwerty say asdfg hello world'
    ref_sent = 'hello world asdfg say qwerty'

    print(eed(hyp_sent, ref_sent))
    print(eed_python(list(hyp_sent), list(ref_sent)))

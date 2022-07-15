import sys
from Levenshtein import distance as lev

GAP = '-'

def fileToString(originals_file, aligned_file):
    origs = open(originals_file, "r")
    aligned = open(aligned_file, "r")
    origs_lines = origs.readlines()
    aligned_lines = aligned.readlines()

    if len(origs_lines) != len(aligned_lines):
        return None, None, None, None

    refs = []
    copies = []
    refs_aligned = []
    copies_aligned = []

    for i in range(0, (len(origs_lines)-(len(origs_lines)%3)+1), 3):
        refs.append(origs_lines[i].strip())
        copies.append(origs_lines[i + 1].strip())
        refs_aligned.append(aligned_lines[i].strip())
        copies_aligned.append(aligned_lines[i + 1].strip())

    '''print(refs)
    print(copies)
    print("______________________________________")
    print(refs_aligned)
    print(copies_aligned)'''

    return refs, copies, refs_aligned, copies_aligned


def levenshtein(ref, copy, ref_aligned, copy_aligned):
    lev_score = lev(ref, copy)
    misses_num = 0
    for i in range(len(ref_aligned)):
        if ref_aligned[i] != copy_aligned[i] or (ref_aligned[i] == "-" and copy_aligned[i] == "-"):
            misses_num = misses_num + 1
    return misses_num - lev_score


def SCSMetric(ref, copy, ref_aligned, copy_aligned):
    m = len(ref)
    n = len(copy)
    scs = SCS(ref, copy, m, n)
    subs = numOfSubstitutions(ref_aligned, copy_aligned)
    metric_ref = len(ref_aligned) - scs + subs

    return metric_ref


def numOfSubstitutions(X, Y):
    m = len(X)
    n = len(Y)
    if n != m:
        print("not aligned")
        print(X)
        print(Y)
    count = 0

    for i in range(m):
        if X[i] != GAP and Y[i] != GAP and X[i] != Y[i]:
            count += 1

    return count


def SCS(X, Y, m, n):
    l = lcs(X, Y, m, n)

    # Result is sum of input string
    # lengths - length of lcs
    return m + n - l


def lcs(X, Y, m, n):
    L = [[0] * (n + 2) for i in range(m + 2)]

    # Following steps build L[m + 1][n + 1]
    # in bottom up fashion. Note that L[i][j]
    # contains length of LCS of X[0..i - 1]
    # and Y[0..j - 1]
    for i in range(m + 1):

        for j in range(n + 1):

            if (i == 0 or j == 0):
                L[i][j] = 0

            elif (X[i - 1] == Y[j - 1]):
                L[i][j] = L[i - 1][j - 1] + 1

            else:
                L[i][j] = max(L[i - 1][j],
                              L[i][j - 1])

    # L[m][n] contains length of
    # LCS for X[0..n - 1] and Y[0..m - 1]
    return L[m][n]


def lcs_algorithm(strand_1: str, strand_2: str):
    m = len(strand_1)
    n = len(strand_2)
    matrix = [[0 for d in range(n + 1)] for e in range(m + 1)]

    # Building the matrix in bottom-up way
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                matrix[i][j] = 0
            elif strand_1[i - 1] == strand_2[j - 1]:
                matrix[i][j] = matrix[i - 1][j - 1] + 1
            else:
                matrix[i][j] = max(matrix[i - 1][j], matrix[i][j - 1])
    index = matrix[m][n]
    size = index
    lcs_algo = [""] * (index + 1)
    lcs_algo[index] = ""

    i = m
    j = n
    while i > 0 and j > 0:
        t = strand_1[i - 1]
        r = strand_2[j - 1]
        if strand_1[i - 1] == strand_2[j - 1]:
            lcs_algo[index - 1] = strand_1[i - 1]
            i -= 1
            j -= 1
            index -= 1

        elif matrix[i - 1][j] > matrix[i][j - 1]:
            i -= 1
        else:
            j -= 1
    return lcs_algo, size


def num_matches(strand_1: str, strand_2: str):
    count = 0
    length = min(len(strand_1), len(strand_2))
    for i in range(length):
        if strand_1[i] == strand_2[i]:
            count += 1

    return count


def lcs_metric(strand_1: str, strand_2: str, align_st1: str, align_st2: str):
    lcs, size = lcs_algorithm(strand_1, strand_2)
    num_of_matches = num_matches(align_st1, align_st2)
    diff = size - num_of_matches
    return diff


def prob_metric(ref_a, copy_a, deletion_prob, insertion_prob, sub_prob):
    prob = 1
    count =0
    for i in range(len(ref_a)):
        if ref_a[i] == '-':
            prob *= float(insertion_prob)
            count += 1
        elif copy_a[i] == '-':
            prob *= float(deletion_prob)
            count += 1
        elif copy_a[i] != ref_a[i]:
            prob *= float(sub_prob)
            count += 1

    return -prob*count


# python metrics.py simple_format.txt output.txt num
# del ins s0ubs
#num = number of strands couples
#del = probability of deletion error
#ins = probability of insertion error
#subs = probability of substitution error

def metrics_runner(args):
    sum_scs=0
    sum_lcs=0
    sum_lev=0
    refs, copies, refs_aligned, copies_aligned = fileToString(args[0], args[1])
    if refs is None:
        print("input error")
    else:
        metrics_file = open(args[3], "w")
        for i in range(len(refs)):
            SCS_metric = SCSMetric(refs[i], copies[i], refs_aligned[i], copies_aligned[i])
            LCS_metric = lcs_metric(refs[i], copies[i], refs_aligned[i], copies_aligned[i])
            levenshtein_metric = levenshtein(refs[i], copies[i], refs_aligned[i], copies_aligned[i])

            sum_scs += SCS_metric
            sum_lcs += LCS_metric
            sum_lev += levenshtein_metric

            # probability_metric = prob_metric(refs_aligned[i], copies_aligned[i], args[4], args[5], args[6])
            metrics_file.write("SCS metric is: " + str(SCS_metric) + "\n")
            metrics_file.write("LCS metric is: " + str(LCS_metric)+ "\n")
            metrics_file.write("levenshtein metric is: " + str(levenshtein_metric) + "\n")
            # print("probabilities function is: " + str(probability_metric))
            metrics_file.write("\n")

        metrics_file.write("To conclude : We have " + str(args[2]) + "strands" + "\n")
        avg_scs = sum_scs / int(args[2])
        avg_lcs = sum_lcs / int(args[2])
        avg_lev = sum_lev / int(args[2])
        metrics_file.write("The SCS average is:" + str(avg_scs)+"\n")
        metrics_file.write("The LCS average is:" + str(avg_lcs)+"\n")
        metrics_file.write("The Levenstein average is:" + str(avg_lev)+"\n")

        return avg_scs, avg_lcs, avg_lev

    return 0,0,0


    '''ref = 'ABCDEFG'
    copy = 'BBCDEFG'
    ref_aligned = 'ABCDEFG' 
    copy_aligned = 'BBCDEFG'
    metric = SCSMetric(ref, copy, ref_aligned, copy_aligned)
    print("SCS metric is: " + str(metric))'''




# See PyCharm help at https://www.jetbrains.com/help/pycharm/

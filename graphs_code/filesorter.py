import sys


def simulatorOutToLists(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()
    i = 0
    refs = []
    copies = []

    while i != len(lines):
        if i + 2 >= len(lines):
            break
        if lines[i + 2] == '\n':
            i += 4
            continue
        else:
            refs.append(lines[i].rstrip())
            copies.append(lines[i + 2].rstrip())
            i += 5
    return refs, copies


def ListsToSimpleFormat(refs, copies):
    simple_in = open("orig_strands.txt", 'w')
    for i in range(len(refs)):
        simple_in.write(refs[i])
        simple_in.write("\n")
        simple_in.write(copies[i])
        if (i == len(refs)-1):
            continue
        simple_in.write("\n")
        simple_in.write("\n")


if __name__ == '__main__':
    args = sys.argv[1:]
    refs, copies = simulatorOutToLists(args[0])
    print(len(refs))
    print(len(copies))

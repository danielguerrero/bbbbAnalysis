def getEfficiencies(filename = 'efficiencies.txt'):

    f = open(filename)
    vals = []
    for l in f:
        l = l.strip()
        if not l: continue
        xx = float(l.split()[1]), l.split()[0]
        vals.append(xx)
    # print vals
    vals.sort()
    vals.reverse()

    # for l in vals:
    #     print '{:<100} {:}'.format(l[1],l[0])

    print '... parsed', len(vals), 'eff paths'
    f.close()
    return vals
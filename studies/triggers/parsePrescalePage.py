def getPrescaleDict(ppgName = 'prescalePage.html'):

    ppg = open(ppgName)
    precales = {}


    # <TR><TD><B>454</B></TD><TD><A HREF=HLTPath?PATHID=2067392><B>HLT_PFMET110_PFMHT110_IDTight_v8</B></A> (2067392)</TD>
    # <TD>0</TD>
    # <TD>0</TD>
    # <TD>0</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>1</TD>
    # <TD>L1_ETM50 OR L1_ETM60 OR L1_ETM70 OR L1_ETM75 OR L1_ETM80 OR L1_ETM85 OR L1_ETM90 OR L1_ETM95 OR L1_ETM100  OR L1_ETM105  OR L1_ETM110  OR L1_ETM115 OR L1_ETM120 OR L1_ETM75_Jet60_dPhi_Min0p4 OR L1_DoubleJetC60_ETM60</TD></TR>

    ncols = 17
    removeVersion = True

    for l in ppg:
        if 'HLT_' in l: ## starts a line of a new trg path
            # cols = (-1,) * ncols
            name = 'HLT_' + l.split('HLT_')[1]
            name = name.split('</B>')[0]
            if removeVersion:
                name = name.rsplit('_', 1)[0]
            # print name

            cols = []
            for idx in range(0, ncols):
                l = next(ppg)
                x = l
                x = x.strip('<TD>')
                x = x.split('<')[0]
                # ps = int(l.strip('<TD>').strip('</TD>'))
                ps = int(x)
                cols.append(ps)
            l = next(ppg)
            l = l.split('<TD>')[1]
            l = l.split('</TD></TR>')[0]
            l = l.split(' OR ')
            l1seeds = [s for s in l if s != '&nbsp;'] ## there is a symbol meaning no L1 at all (for calib triggers)
            cols.append(l1seeds) ## the L1 seeds
            prescs = tuple(cols)
            precales[name] = prescs

    ppg.close()
    # print precales
    # print precales['HLT_DoubleJet90_Double30_DoubleBTagCSV_p087']

    print '... got', len(precales), 'paths'
    return precales

def getInitRunPrescaleDictL1(ppgName = 'L1prescalePage.html'):
    ppg = open(ppgName)
    prescalesL1 = {}
    for l in ppg:
        if not 'ALIGN=CENTER' in l:
            continue
        if not '</B></TD><TD' in l:
            continue ## some strange lines with no L1 seed associated
        tokens = l.split('</B></TD><TD')
        # print "I hve" , len(tokens)
        # print tokens
        init_prescale = tokens[-3]
        init_prescale = init_prescale.split('>')[-1]
        init_prescale = int(init_prescale)
        seedname = tokens[1]
        seedname = seedname.strip('><B>')
        prescalesL1[seedname] = init_prescale
        # print seedname, '===>', init_prescale
    ppg.close()
    return prescalesL1

import re
def getPrescaleDictL1(ppgName = 'L1prescalePageColumns.html'):
    ppg = open(ppgName)
    prescalesL1 = {}
    # ncols = 17
    for l in ppg:
        if 'L1_' in l:
            l = l.strip()
            # print l
            algoname = l.split('<B>')[-1]
            algoname = algoname.split('<')[0]
            l = next(ppg)
            # tokens = l.split('</B></TD><TD><B>')
            # tokens[0] = tokens[0].split('>')[-1]
            # tokens[-1] = tokens[-1].split('<')[0]
            # print tokens
            tokens = l.split('>')
            tokens = [re.sub("[^0-9]", "", x) for x in tokens]
            tokens = [int(x) for x in tokens if x]
            # print tokens
            cols = []
            for val in tokens:
                cols.append(int(val))
            # print algoname, cols
            prescalesL1[algoname] = cols
    ppg.close()
    return prescalesL1

import csv
def getHLTMenuStatusFromGDoc(csvGDocName = 'Final_Menu_2017.tsv', removeVersion=True):
    ### to get the inoput file, go to https://docs.google.com/spreadsheets/d/1SqSKATiHcZdSB-t7rMb0SRSUL4Ot1jr2TExjnMrDUpA/edit#gid=0
    ## end export it as a tsv spreadsheet

    # gdoc = open(csvGDocName)
    res = {}
    with open(csvGDocName) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in reader:
            # print row
            if not 'HLT_' in row[2]: continue
            pathName = row[2]
            if removeVersion:
                pathName = pathName.rsplit('_', 1)[0]
            pathStatus = row[6]   
            res[pathName] = pathStatus                  
    return res
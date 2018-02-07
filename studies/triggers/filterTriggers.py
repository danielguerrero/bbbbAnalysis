import parsePrescalePage
import parseEfficiencies

prescales          = parsePrescalePage.getPrescaleDict('prescalePage.html') ## a dictionary
prescalesL1InitRun = parsePrescalePage.getInitRunPrescaleDictL1('L1prescalePage.html') ## a dictionary
prescalesL1        = parsePrescalePage.getPrescaleDictL1('L1prescalePageColumns.html') ## a dictionary
effs               = parseEfficiencies.getEfficiencies('efficiencies.txt')  ## a list

## this is very good, but I only have the spreadsheet for 2017
# status = parsePrescalePage.getHLTMenuStatusFromGDoc(csvGDocName = 'Final_Menu_2017.tsv')
# print status

effs_names  = [x[1] for x in effs]
presc_names = prescales.keys()

presc_not_in_eff  = [x for x in presc_names if x not in effs_names]
effs_not_in_presc = [x for x in effs_names if x not in presc_names]

# print prescales
# print '\n\n\n'
# print prescalesL1

# the two collections nto nece

print "=== Trgs found in WBM page but not in MC ==="
for a in presc_not_in_eff:
    print a
print ""

print "=== Trg found in MC but not in WBM ==="
for a in effs_not_in_presc:
    print a

print "\n\n\n"

####### filter the triggers

min_eff = 0.02 ## at least 2%
found = []
for eff,path in effs:
    if eff > min_eff:
        
        ### ucnomment for the HLT + L1 from WBM page
        if not path in prescales:
            if eff > min_eff:
                print "** path: {:<100}".format(path) , " has an eff of {:.1f}%".format(100.*eff) , "but not found in WBM"
            continue
        if prescales[path][1] == 1: ## col 1 must be unprescaled

            ## also check if at least one seed is unprescaled
            one_unpresc = False
            # print path
            for l1seed in prescales[path][-1]:
                if prescalesL1[l1seed][1] == 1: ## check only col 1 (not 0 that would be emergency)
                    one_unpresc = True
                    break
            
            # print "UNPRESCALED: ", path, " . eff =", 100.*eff , "%"
            if one_unpresc:
                found.append((eff,path))

        # ## check from the gdoc spreacsheet
        # if not path in status:
        #     print "** path: " , path , " has an eff of", 100.*eff , "but not found in the GDOC"
        #     continue
        # the_status = status[path]
        # if the_status == 'unprescaled':
        #     found.append((eff,path))

print "\n\nI found", len(found), "triggers UNPRESCALED and eff > ", 100.*min_eff , "%"
for eff, path in found:
    print "{:<100} with eff = {:.1f}%".format(path, 100.*eff)
    # print "UNPRESCALED: ", path, " . eff =", 100.*eff , "%"

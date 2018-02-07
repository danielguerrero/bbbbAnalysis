#### a script that generates the code that defines the member variables of NanoAODTree_SetBranchImpl
# practical to use as both MakeSelector nor MakeClass are capable of making this kind of code
# usage: python generate... <MC_filename> <output_name (optional)>
# python generate_NanoAODTree_SetBranchImpl_vars.py root://cmsxrootd.fnal.gov//store/user/lcadamur/NanoTest/VBFHHTo4B_CV_1_C2V_1_C3_1_13TeV-madgraph/VBFHbbNanoAOD_10Gen2018/180110_005708/0000/test80X_NANO_1.root thecode.txt

import ROOT
import os
import sys
import fnmatch

if len (sys.argv) < 2:
    print "please provide a file"
    sys.exit()

fname = sys.argv[1]
oname = sys.argv[2] if len(sys.argv) > 2 else None

## the parrents to skip - do not create anything for these
patt_to_skip = [
    "HLT_*",
    "Flag_*",
    "L1simulation_step",
    "HLTriggerFirstPath",
    "HLTriggerFinalPath"
]

# ### after running on VBF HH > 4b, I get
# Fond these maximum sizes
#            nElectron  --> 8
#              nFatJet  --> 6
#           nGenJetAK8  --> 6
#              nGenJet  --> 25
#             nGenPart  --> 153
#           nGenVisTau  --> 3
#        nLHEPdfWeight  --> 101
#      nLHEScaleWeight  --> 9
#                 nJet  --> 29
#                nMuon  --> 6
#              nPhoton  --> 10
#    nGenDressedLepton  --> 2
#     nSoftActivityJet  --> 6
#              nSubJet  --> 10
#                 nTau  --> 6
#             nTrigObj  --> 31
#             nOtherPV  --> 3
#                  nSV  --> 16

max_sizes = {
    "nElectron"         : 100,
    "nFatJet"           : 100,
    "nGenJetAK8"        : 100,
    "nGenJet"           : 100,
    "nGenPart"          : 400,
    "nGenVisTau"        : 100,
    "nLHEPdfWeight"     : 100,
    "nLHEScaleWeight"   : 100,
    "nJet"              : 100,
    "nMuon"             : 100,
    "nPhoton"           : 100,
    "nGenDressedLepton" : 100,
    "nSoftActivityJet"  : 100,
    "nSubJet"           : 100,
    "nTau"              : 100,
    "nTrigObj"          : 100,
    "nOtherPV"          : 100,
    "nSV"               : 100,
    "DEFAULT"           : 100,
}

## ask root to generate a MakeClass, so I can easily get the types
fscript = open('the_generate_script.txt', 'w')
fscript.write('Events->MakeClass("the_generate_class")\n')
fscript.write('.q\n')
fscript.close()
os.system("root -l %s < the_generate_script.txt" % fname )
os.system('rm the_generate_script.txt')

## parse the MakeClass .h file and redefine as I need the branches
lines = []

print "... parsing the .h file"
fMakeClass = open('the_generate_class.h')
started = False
for l in fMakeClass:
    if r'// Declaration of leaf types' in l: # this comment marks the beginning of the declaration block
        started = True
        continue
    if r'// List of branches' in l: # this comment marks the end of the declaration
        break
    
    if started: # I am parsing the vars
        l = l.strip()
        if not l: continue
        tokens = l.split()
        tipo = tokens[0]
        var  = tokens[1].replace(';', '') # no semicolon at the end
        isArray = True if '[' in var and ']' in var else False
        if isArray:
            var = var.split('[', 1)[0] ## trim everything in the parenthesis

        to_skip = False
        for pts in patt_to_skip:
            to_skip = (to_skip or fnmatch.fnmatch(var, pts))
        if to_skip:
            continue

        lines.append((tipo, var, isArray))
fMakeClass.close()
os.system('rm the_generate_class.h')
os.system('rm the_generate_class.C')


## build the output

## make a pass to get the column size
cols0 = [len(l[0]) for l in lines] 
cols1 = [len(l[1]) for l in lines]
maxcol0 = max(cols0)
maxcol1 = max(cols1)

if oname:
    ofile = open(oname, 'w')

val_proto = 'ReaderValuePatch<{0}> {2} {1} {3} {{tree_, "{1}"}};'
arr_proto = 'ReaderArrayPatch<{0}> {4} {1} {5} {{tree_, "{1}", {2}, {3}}};'
sbnames = []
for l in lines:
    
    filler0 = ' ' * (maxcol0 - len(l[0]) + 1)
    filler1 = ' ' * (maxcol1 - len(l[1]) + 1)
    
    if not l[2]:
        oline = val_proto.format(l[0], l[1], filler0, filler1)
    else:
        ## I need to look for the size branch. all are in the format n<NOME>_...
        sbname = 'n' + l[1].split('_',1)[0]
        maxsize = max_sizes[sbname] if sbname in max_sizes else max_sizes["DEFAULT"]
        oline  = arr_proto.format(l[0], l[1], sbname, maxsize, filler0, filler1)
        if not sbname in sbnames: sbnames.append(sbname)

    if not oname:
        print oline
    else:
        ofile.write(oline + '\n')

if oname:
    ofile.close()
    print "Please verify that the max_size is enough for the following branches"
    print "You can check the size with getMaxBranchSize.exe <filelist_MC>"
    for sbn in sbnames:
        tipo = [l for l in lines if l[1] == sbn][0][0]
        print sbn
        if tipo != 'UInt_t':
            print "** WARNING: this type is", tipo, "instead of UInt_t"
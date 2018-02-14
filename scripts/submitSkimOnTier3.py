import os
import sys
import argparse
import getpass

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def splitInBlocks (l, n):
    """split the list l in n blocks of equal size"""
    k = len(l) / n
    r = len(l) % n

    i = 0
    blocks = []
    while i < len(l):
        if len(blocks)<r:
            blocks.append(l[i:i+k+1])
            i += k+1
        else:
            blocks.append(l[i:i+k])
            i += k

    return blocks

def formOutName (outputDir, tag, outputName):
    oput = outputDir
    if oput[-1] != '/':
        oput += '/'
    oput += tag + '/' + outputName 
    return oput

def parseInputFileList (fileName) :
    filelist = []
    with open (fileName) as fIn:
        for line in fIn:
            line = (line.split("#")[0]).strip()
            if line:
                filelist.append(line)
    return filelist

def writeln(f, line):
    f.write(line + '\n')

###########

parser = argparse.ArgumentParser(description='Command line parser of skim options')
parser.add_argument('--input',      dest='input',  help='input filelist',   required = True)
parser.add_argument('--tag',        dest='tag',    help='production tag',   required = True)
parser.add_argument('--njobs',      dest='njobs',  help='njobs', type=int,  default=500)
##
############################################################
## all the following are for "expert" use, not meant to be used by default
##
parser.add_argument('--outputName', dest='oname',  help='the name of the output (if not given, auto from filelist)', default = None)
parser.add_argument('--outputDir',  dest='odir',   help='the base EOS output directory. Use a {0} for username placeholder, or give it explicitely', default = "root://cmseos.fnal.gov//store/user/{0}/bbbb_ntuples/")
## if leaving a {0} argument, odir will be formatted with the current username
## otherwise a manual fixed path can also be given
##
parser.add_argument('--no-tar',         dest='tar',        help='do not tar the executable folder',       action='store_false', default=True)
parser.add_argument('--no-xrdcp-tar',   dest='xrdcptar',   help='do not xrdcp the tar to EOS',       action='store_false', default=True)
parser.add_argument('--no-xrdcp-flist', dest='xrdcpflist', help='do not xrdcp the filelist to EOS',  action='store_false', default=True)
parser.add_argument('--dry-run',        dest='dryrun',     help='dry run without launching',         action='store_true',  default=False)
parser.add_argument('--verbose',        dest='verbose',    help='set verbose mode',                  action='store_true',  default=False)
parser.add_argument('--force',          dest='force',  help='force submission (delete work folder if existing)', action='store_true', default=False)
##
# parser.add_argument('--cfg',       dest='skim_cfg',   help='the name of the skim cfg',          default='skim.cfg')

# args = parser.parse_args()
args, unknown = parser.parse_known_args()

#### the output is saved in
#### odir / tag / oname
##### e.g.
#### root://cmseos.fnal.gov//store/user/lcadamur/bbbb_ntuples/' + TAG + '/' + ONAME

##################
executable = 'bin/skim_ntuple.exe'
##################

username = getpass.getuser()
print "... Welcome", username

oname = args.oname
if not oname:
    # print "A:" , args.input
    # print "B:" , args.input.rsplit(r'/', 1)[-1]
    # print "C:" , args.input.rsplit(r'/', 1)[0].rsplit('.', 1)
    oname = args.input.rsplit(r'/', 1)[-1].rsplit('.', 1)[0]
oname = 'SKIM_' + oname

odir = args.odir.format(username)
oLFN_base = formOutName(odir, args.tag, oname)
if oLFN_base[-1] == '/':
    oLFN_base = oLFN_base[:-1]

jobsDir                = 'jobs_' + args.tag + '/' + oname
outListNameBareProto   = 'filelist_{0}.txt'
outScriptNameBareProto = 'job_{0}.sh'
outListNameProto       = (jobsDir + '/' + outListNameBareProto)
outScriptNameProto     = (jobsDir + '/' + outScriptNameBareProto)
EOSfilelistBase        = oLFN_base + '/filelist'
EOSfilelistProto       = EOSfilelistBase + '/' + outListNameBareProto
outFileNameProto       = 'bbbbNtuple_{0}.root'

cmssw_base    = os.environ['CMSSW_BASE']
cmssw_version = os.environ['CMSSW_VERSION']
scram_arch    = os.environ['SCRAM_ARCH']



tarName      = 'bbbbAnalysis.tar.gz' #%s_tar.tgz' % cmssw_version
bbbbWorkDir  = os.getcwd()
tarLFN       = bbbbWorkDir + '/tars/' + tarName

tarEOSdestLFN         = ('root://cmseos.fnal.gov//store/user/{0}/bbbb_analysis_tar/' + tarName).format(username)
filelistEOSdestLFNdir = ('root://cmseos.fnal.gov//store/user/{0}/bbbb_analysis_tar/' + tarName).format(username)

inputfiles = parseInputFileList (args.input)    ## parse input list
njobs      = args.njobs if args.njobs <= len (inputfiles) else len (inputfiles)
fileblocks = splitInBlocks (inputfiles, njobs)
if njobs != len(fileblocks):
    print "** ERROR: length of file lists and njobs do not match, something went wrong"
    sys.exit()

## build the full command to launch
## all to format on the job number

## come commands MUST NOT be forwarded - so check and throw in case
if '--input'  in unknown:
    print "** you MUST NOT pass a --input as it is generated by the script, exiting..."
    sys.exit()
if '--output' in unknown:
    print "** you MUST NOT pass a --output as it is generated by the script, exiting..."
    sys.exit()

commands = [
    executable,
    '--input %s'  % outListNameBareProto,
    '--output %s' % outFileNameProto,
    ]
full_command = ' '.join(commands)
## now forward all the other commands to full_command
full_command += ' ' + ' '.join(unknown)

#######

print "** Skim of filelist      :", args.input
print "** Production tag        :", args.tag
print "** Njobs                 :", njobs
print "** Working in            :", jobsDir
print "** Saving output to      :", oLFN_base
print "** Saving executables in :", tarEOSdestLFN
print "** The script will execute :", full_command

if os.path.isdir(jobsDir):
    if not args.force:
        print "... working folder", jobsDir, " already exists, exit"
        sys.exit()
    else:
        print "... working folder", jobsDir, " already exists, FORCING ITS DELETION "
        cmd = 'rm -r ' + jobsDir
        print ">>> ", cmd
        os.system(cmd)        
os.system('mkdir -p ' + jobsDir)

## filelist
for n in range(0, njobs):
    outListName = outListNameProto.format(n)
    jobfilelist = open(outListName, 'w')
    for f in fileblocks[n]: jobfilelist.write(f+"\n")
    jobfilelist.close()

################################################################
## prepare the tarball of the executables
if args.tar:
    # tar -zcvf <Tar output name> <folder to tar>
    # note -v is verbose option
    ## tp exlcude something:
    ### command = 'tar --exclude="{0}" --exclude="{1}" -zcf {2} -C {3} {4}'.format(excludePath, excludePath2, tarLFN, cmsswWorkDir, cmssw_version)

    ### NOTE: I must be in bbbb
    to_include = [
        'bin/',
        'lib/',
        'config/',
    ]

    command = 'tar -zcf {0} '.format(tarLFN)
    for ti in to_include:
        command += ti + ' '
    # command += ' -zcf {0} -C {1} {2}'.format(tarLFN, cmsswWorkDir, cmssw_version)

    print '** INFO: Going to tar executable folder into', tarName
    if args.verbose: print "** INFO: executing:", command
    os.system(command)
    print '** INFO: tar finished and saved in:', tarLFN
else:
    print '** INFO: Not going to tar executable folder, using', tarLFN

# script
for n in range(0, njobs):
    this_full_command = full_command.format(n)
    outListName     = outListNameProto.format(n)
    outListNameBare = outListNameBareProto.format(n)
    # outputFileName = 'ntuple_%i.root' % n
    outputFileName = outFileNameProto.format(str(n))
    outputEOSName  = '%s/output/%s' % (oLFN_base, outputFileName)
    outScriptName  = outScriptNameProto.format(n)
    outScript      = open(outScriptName, 'w')
    writeln(outScript, '#!/bin/bash')
    writeln(outScript, '{') ## start of redirection..., keep stderr and stdout in a single file, it's easier
    writeln(outScript, 'echo "... starting job on " `date` #Date/time of start of job')
    writeln(outScript, 'echo "... running on: `uname -a`" #Condor job is running on this node')
    writeln(outScript, 'echo "... system software: `cat /etc/redhat-release`" #Operating System on that node')
    writeln(outScript, 'source /cvmfs/cms.cern.ch/cmsset_default.sh')
    writeln(outScript, 'export SCRAM_ARCH=%s' % scram_arch)
    # writeln(outScript, 'cd %s/src/' % cmssw_version)
    writeln(outScript, 'eval `scramv1 project CMSSW %s`' % cmssw_version)
    writeln(outScript, 'cd %s/src' % cmssw_version)
    writeln(outScript, 'eval `scramv1 runtime -sh`')
    writeln(outScript, 'echo "... retrieving bbbb executables tarball"')
    writeln(outScript, 'xrdcp -f -s %s .' % tarEOSdestLFN) ## force overwrite CMSSW tar
    writeln(outScript, 'echo "... uncompressing bbbb executables tarball"')
    writeln(outScript, 'tar -xzf %s' % tarName)
    writeln(outScript, 'rm %s' % tarName)
    writeln(outScript, 'echo "... retrieving filelist"')
    writeln(outScript, 'xrdcp -f -s %s .' % EOSfilelistProto.format(n)) ## force overwrite file list
    # writeln(outScript, 'echo "... listing files"')
    # writeln(outScript, 'ls -altrh')
    writeln(outScript, 'export CPP_BOOST_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt')
    writeln(outScript, 'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:./lib:${CPP_BOOST_PATH}/lib')
    # writeln(outScript, 'echo "Arguments passed to this script are: for 1: $1, and for 2: $2"')
    # writeln(outScript, 'echo $LD_LIBRARY_PATH')
    # writeln(outScript, 'ls /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/lib')
    writeln(outScript, 'echo "... executing the command:"')
    writeln(outScript, 'echo "%s"' % this_full_command)
    writeln(outScript, 'echo "... starting execution"')
    # writeln(outScript, 'cmsRun %s inputFiles_load=%s outputFile=%s' % (executable, outListNameBare, outputFileName))
    writeln(outScript, '%s' % this_full_command)
    writeln(outScript, 'echo "... execution finished with status $?"')
    writeln(outScript, 'echo "... copying output file %s to EOS in %s"' % (outputFileName, outputEOSName))
    writeln(outScript, 'xrdcp -s %s %s' % (outputFileName, outputEOSName)) ## no not force overwrite output in destination
    writeln(outScript, 'echo "... copy done with status $?"')
    # writeln(outScript, 'remove the input and output files if you dont want it automatically transferred when the job ends')
    # writeln(outScript, 'rm nameOfOutputFile.root')
    # writeln(outScript, 'rm Filename1.root')
    # writeln(outScript, 'rm Filename2.root')
    writeln(outScript, 'cd ${_CONDOR_SCRATCH_DIR}')
    writeln(outScript, 'rm -rf %s' % cmssw_version)
    writeln(outScript, 'echo "... job finished with status $?"')
    writeln(outScript, 'echo "... finished job on " `date`')
    writeln(outScript, 'echo "... exiting script"')
    writeln(outScript, '} 2>&1') ## end of redirection
    outScript.close()


##############################
#### Ship the tarball and submit the jobs

if args.xrdcptar:
    print "** INFO: copying executables tarball to:", tarEOSdestLFN
    command = 'xrdcp -f -s %s %s' % (tarLFN, tarEOSdestLFN)
    os.system(command)
    if args.verbose: print "** INFO: executing:", command
else:
    print "** INFO: not going to xrdcp the CMSSW tarball to EOS, assuming it exists at", tarEOSdestLFN

if args.xrdcpflist:    
    
    print "** INFO: copying input filelists to:", EOSfilelistProto.format('*')
    command = 'eos root://cmseos.fnal.gov mkdir -p %s' % EOSfilelistBase.replace('root://cmseos.fnal.gov/', '/eos/uscms')
    if args.verbose: print "** INFO: executing:", command
    os.system(command)
    
    command = 'xrdcp -f -s %s %s' % (outListNameProto.format('*'), EOSfilelistBase)
    if args.verbose: print "** INFO: executing:", command
    os.system(command)

    # for n in range(0, njobs):
    #     command = 'xrdcp -f -s %s %s' % (outListNameProto.format(n), EOSfilelistProto.format(n))
    #     if args.verbose: print "** INFO: executing:", command
    #     os.system(command)
else:
    print "** INFO: not going to xrdcp the filelistsi to EOS, they exists at", EOSfilelistProto.format('*')


## set directory to job directory, so that logs will be saved there
os.chdir(jobsDir)
for n in range(0, njobs):
    command = "%s/scripts/t3submit %s" % (bbbbWorkDir, outScriptNameBareProto.format(n))
    if not args.dryrun:
        if args.verbose: print "** INFO: submit job with command", command
        os.system(command)
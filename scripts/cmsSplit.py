#https://github.com/ICBPHCMS/BPHParkingAnalysis/blob/master/NtupleProducer/scripts/cmsSplit.py

###!/usr/bin/python                                                                                                                                     
import sys
import os
import commands
import optparse
import fnmatch
import time
import math
import re

### parsing input options 
def parseOptions():
    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-a', '--anType', dest='ANALYSISTYPE', type='string', default='',help='cmsRun, script, scriptAndJOBID')
    parser.add_option('-c', '--cfg', dest='CONFIGFILE', type='string', default='',help='CMSSW config template name, if empty string the deafult one will be used')
    parser.add_option('-t', '--tag', dest='TAGNAME', type='string', default='',help='kind of job')
    parser.add_option('', '--dasquery', action='store_true', dest='DASquery',  default=False, help='retrieve list on DAS')
    parser.add_option('', '--folderquery', action='store_true', dest='FOLDERquery',  default=False, help='retrieve list on DAS')
    parser.add_option('', '--listquery', action='store_true', dest='LISTquery',  default=False, help='retrieve list from .txt')
    parser.add_option('', '--das', dest='DASCOMMAND', type='string', default='', help='set input dataset on das')
    parser.add_option('', '--folder', dest='FOLDERNAME', type='string', default='', help='set input dataset in a folder')
    parser.add_option('-j', '--jobs',  dest='NJOBS',  type=int,       default=-1,  help='total number of jobs')
    parser.add_option('-f', '--filesperjob', dest='FILESPERJOB', type=int, default=-1,   help='number of events per job')
    parser.add_option('-s', '--storeArea', dest='STOREAREA', type='string', default='', help='set output storage folder')
    parser.add_option('-i', '--inputList', dest='inputLIST', type='string', default='', help='give input .txt')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()
        
    # make sure CMSSW is set up
    if 'CMSSW_BASE' not in os.environ:
        print 'ERROR: CMSSW does not seem to be set up. Exiting...'
        sys.exit()
        
    # set the configfile or raise problem
    if (opt.ANALYSISTYPE != 'cmsRun' and opt.ANALYSISTYPE != 'script' and opt.ANALYSISTYPE != 'scriptAndJOBID'):
        print 'Need to specify the type of analysis. Exiting...'
        sys.exit()

    # set the configfile or raise problem
    if (opt.CONFIGFILE == ''):
        print 'Need to specify an input cfg.  Exiting...'
        sys.exit()

    if (opt.DASquery==False and opt.FOLDERquery==False and opt.LISTquery==False):
        print 'Need to specify an query place.  Exiting...'
        sys.exit()

    if (opt.NJOBS==-1 and opt.FILESPERJOB==-1):
        print 'Need to specify a splitting criteria... NJOBS or FILESPERJOB.  Exiting...'
        sys.exit()

    if (opt.NJOBS!=-1 and opt.FILESPERJOB!=-1):
        print 'Pleas make sure that splitting criteria are consistent... NJOBS or FILESPERJOB.  Exiting...'
        sys.exit()
       
       
def printOptions():
    print ' NJOBS = ', opt.NJOBS
    print ' FILESPERJOB = ', opt.FILESPERJOB
    print ' CONFIGFILE = ', opt.CONFIGFILE
    print ' TAGNAME = ', opt.TAGNAME
    print ' STOREAREA = ', opt.STOREAREA
    if opt.DASquery:
        print ' DAS query is option '
        print ' command => ', opt.DASCOMMAND
    if opt.FOLDERNAME:
        print ' FOLDER query is option '
        print ' command => ', opt.FOLDERNAME
    if opt.LISTquery:
        print ' LIST query is option '
        print ' command => ', opt.inputLIST



def getInputFileList(DASquery):
    inputList = []
    if not DASquery:
        pattern = '*.root'
        inputList = [f for f in os.listdir(opt.FOLDERNAME) if (os.path.isfile(os.path.join(opt.FOLDERNAME, f)) and (fnmatch.fnmatch(f, pattern)))]
    else:
        cmd = 'dasgoclient -query="file dataset='+opt.DASCOMMAND+' instance=prod/phys03"'
        #cmd = 'dasgoclient -query="file dataset='+opt.DASCOMMAND+'"'
        print 'cmd = ', cmd
        status, thisoutput = commands.getstatusoutput(cmd)
        if status !=0:
            print "Error in processing command: "+cmd
            print "Did you forget running voms-proxy-init?"
            sys.exit(1)
        inputList=thisoutput.split()
    return inputList




### submission of GSD/RECO production
def submitProduction():
    # parse the arguments and options
    global opt, args, particles
    parseOptions()


    printOptions()


    nJOBS = opt.NJOBS
    nFILESPERJOB = opt.FILESPERJOB


    print '[Submitting jobs]'
    jobCount=0


    DASquery=False
    if opt.DASquery:
        DASquery=True

    inputFilesList = ""
    if opt.LISTquery:
        file = open(opt.inputLIST, 'r')
        inputFilesList = file.read().split()
    else:
        inputFilesList = getInputFileList(DASquery)
        if len(inputFilesList) == 0:
            print 'list empty...'
            sys.exit()

    if (nFILESPERJOB == -1):
        nFILESPERJOB = int(math.ceil(float(len(inputFilesList))/float(nJOBS)))
    else:
        nJOBS = int(math.ceil(float(len(inputFilesList))/float(nFILESPERJOB)))


    # read the template file in a single string
    f_template= open(opt.CONFIGFILE, 'r')
    template = f_template.read()
    f_template.close()

    currentDir = os.getcwd()
    SCRAM_ARCH = os.getenv('SCRAM_ARCH')

    storeArea = opt.STOREAREA+'/'+opt.TAGNAME
    if not os.path.exists(storeArea):
        os.makedirs(storeArea)
    #os.system('mkdir %s' %(storeArea));

    outDir = currentDir+"/"+opt.TAGNAME;
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    #os.system('mkdir %s' %(outDir));

    outLancia = open(currentDir+'/launch_'+opt.TAGNAME+'.sh', 'w')

    for job in range(1,int(nJOBS)+1):
        subDir = outDir+"/JOB_"+str(job)
        if not os.path.exists(subDir):
            os.makedirs(subDir)
        #os.system('mkdir %s' %(subDir));

        s_template=template
        basename=opt.TAGNAME+'_'+str(job)
    
        inputFilesListPerJob = inputFilesList[(job-1)*nFILESPERJOB:(job)*nFILESPERJOB]
        if len(inputFilesListPerJob)==0: continue

        inputFiles = '"' + '", "'.join([str(f) for f in inputFilesListPerJob]) + '"'
            
        if not DASquery:
            recoInputPrefix = 'file:'+opt.FOLDERNAME+'/'
            inputFiles = '"' + '", "'.join([recoInputPrefix+str(f) for f in inputFilesListPerJob]) + '"'            
            
        if opt.LISTquery:
            inputFiles = '"'.join([str(f) for f in inputFilesListPerJob])

        if(opt.LISTquery and nFILESPERJOB != 1 and opt.ANALYSISTYPE != "cmsRun"):
            inputFiles = open(subDir+'/tempFileList_'+str(job)+'.txt', 'w')
            for f in inputFilesListPerJob:
                #print 'single file', str(f)
                inputFiles.write(str(f) + "\n")

        if(opt.LISTquery and nFILESPERJOB != 1 and opt.ANALYSISTYPE != "cmsRun"):
            s_template=s_template.replace('DUMMYINPUTFILELIST', inputFiles.name)
        else:
            s_template=s_template.replace('DUMMYINPUTFILELIST',inputFiles)

        if not opt.LISTquery:
            s_template=s_template.replace('DUMMYEVTSPERJOB',str(-1))

        outfile = storeArea+"/"+basename +'.root'
        s_template=s_template.replace('DUMMYOUTFILENAME',outfile)

        if(opt.ANALYSISTYPE == "scriptAndJOBID"):
            s_template=s_template.replace('DUMMYJOBID',str(job))
        
        if(opt.ANALYSISTYPE == "cmsRun"):
            cfgfile = basename +'.py'
        if(opt.ANALYSISTYPE == "script" or opt.ANALYSISTYPE == "scriptAndJOBID"):
            cfgfile = basename +'.sh'
        write_template = open(subDir+'/'+cfgfile, 'w')
        write_template.write(s_template)
        write_template.close()

        if(opt.ANALYSISTYPE == "cmsRun"):
            shFile = 'bjob_'+str(job)+'.sh'
            outScript = open(subDir+'/'+shFile,'w');
            outScript.write('#!/bin/sh \n');
            outScript.write('cd %s \n' %(subDir));
            outScript.write('export %s \n' %(SCRAM_ARCH));
            outScript.write('source /cvmfs/cms.cern.ch/cmsset_default.sh \n');    
            outScript.write('eval `scramv1 ru -sh` \n');
            outScript.write('export X509_USER_PROXY=~/myVoms/x509up_u1282033 \n');        
            outScript.write('cmsRun  %s \n' %(cfgfile) )
            outScript.write( '\n ' );
            os.system('chmod 777 %s/%s' %(subDir, shFile));

        if(opt.ANALYSISTYPE == "script" or opt.ANALYSISTYPE == "scriptAndJOBID"):
            shFile = 'bjob_'+str(job)+'.sh'
            outScript = open(subDir+'/'+shFile,'w');
            outScript.write('#!/bin/sh \n');
            outScript.write('cd %s \n' %(subDir));
            outScript.write('export %s \n' %(SCRAM_ARCH));
            outScript.write('source /cvmfs/cms.cern.ch/cmsset_default.sh \n');    
            outScript.write('eval `scramv1 ru -sh` \n');
            outScript.write('export X509_USER_PROXY=~/myVoms/x509up_u1282033 \n');
            outScript.write('source  %s \n' %(cfgfile) )
            outScript.write( '\n ' );
            os.system('chmod 777 %s/%s' %(subDir, shFile));


        # qsub -q hep.q -l h_rt=10800  (in sec) pathTosh
        cmd = ('qsub -q hep.q -l h_rt=10800 %s/%s ' %(subDir, shFile))
        outLancia.write(' qsub -q hep.q -l h_rt=10800 %s/%s \n' %(subDir, shFile));





### run the submitProduction() as main
if __name__ == "__main__":
submitProduction()

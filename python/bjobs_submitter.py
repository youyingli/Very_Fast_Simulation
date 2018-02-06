#!/usr/bin/env python

from optparse import OptionParser
import os, sys
import random

def Option_Parser(argv):

    usage='usage: %prog [options] arg\n'
    usage+='The script handles the MC generation by parallel calculation.\n'
    usage+='It can automatically produce some scripts for given MCTool and submit to batch system.\n'
    usage+='For more information, please see README !!!\n'
    parser = OptionParser(usage=usage)

    parser.add_option('-i', '--input',
            type='str', dest='input', default='myMCgenerated',
            help='Input is .txt file contained mc sample names you want to generate. These names must be the same as cards'
            )
    parser.add_option('-o', '--outdir',
            type='str', dest='outdir', default='/afs/cern.ch/user/<Y>/<YOURNAME>',
            help='Output directory'
            )
    parser.add_option('-t', '--tag',
            type='str', dest='tag', default='Standard',
            help='tag output MC samples'
            )
    parser.add_option('-q', '--queue',
            type='str', dest='queue', default='8nh',
            help='Queue for batch system'
            )
    parser.add_option('--madgraph5_pythia8',
            action='store_true', dest='madgraph5_pythia8',
            help='Run madgraph5_pythia8'
            )
    parser.add_option('--isNLO',
            action='store_true', dest='isNLO',
            help='Run madgraph5_pythia8 at NLO'
            )
    parser.add_option('--FastSim',
            action='store_true', dest='FastSim',
            help='Run the detector simulation by Delphes3'
            )
    parser.add_option('-d', '--delphes_card',
            type='str', dest='delphes_card', default='delphes_card_CMS.tcl',
            help='Datacard needed when running the detector simulation by Delphes3'
            )

    (options, args) = parser.parse_args(argv)
    return options

def decode_txt (filename) :

    file = open(filename, 'r')
    return file.readlines()

def batch_job_submitter (script_list, queue) :
    for script in script_list :
        os.system('bsub -q %s -o %s %s' % (queue, script.replace('.sh', '.log'), script))

def bjobs_setup (argv):

    package_path = os.environ.get('VFSIM_PACKAGE_PATH')
    options = Option_Parser(argv)

    #Madgraph5 + Pythia8
    if options.madgraph5_pythia8 :
        random.seed()
        bjob_list = list()
        import madgraph5_pythia8_plugin
        bjob_scripts_path = '%s/bjobs/madgraph5_pythia8/%s' % (package_path, options.tag)

        for mcline in decode_txt(options.input):
            mclist = mcline.split()
            mc_job = bjob_scripts_path + '/' + mclist[0]
            os.system('mkdir -p %s' % mc_job)
            for i in range(int(mclist[2])) :
                output = options.outdir + '/' + options.tag + '/' + mclist[0]
                outputdir = '%s/output%d' % (output, int(i))
                os.system('mkdir -p ' + outputdir)
                mc_base = package_path + '/Cards/MadGraph5_aMCatNLO/' + mclist[0]
                if os.path.isfile(mc_base + '/cript.dat') :
                    print 'Cannot find %s \n' % (mc_base + '/script.dat')
                    exit()
                USER = os.environ.get('USER')
                tmp_dir = '/tmp/%s/%s_%s_tmp' % (USER, mclist[0], options.tag)
                mc_script = '%s/madgraph5_script%d.dat' % (mc_job, int(i))
                os.system("sed -e 's/MYPATH/%s/" % tmp_dir.replace('/', '\/') + "g' %s > %s" % (mc_base + '/script.dat', mc_script))

                #Cars (If these cards are not present, it will use default cards)
                run_card = (os.popen("ls %s | grep 'run'" % mc_base).readlines())
                pythia8_card = (os.popen("ls %s | grep 'pythia8'" % mc_base).readlines())
                shower_card = (os.popen("ls %s | grep 'shower'" % mc_base).readlines())
                madspin_card = (os.popen("ls %s | grep 'madspin'" % mc_base).readlines())

                if pythia8_card != [] :
                    os.system("echo '%s' >> %s" %  (mc_base + '/' + pythia8_card[0], mc_script))
                if shower_card != [] :
                    os.system("echo '%s' >> %s" %  (mc_base + '/' + shower_card[0], mc_script))
                if madspin_card != [] :
                    os.system("echo '%s' >> %s" %  (mc_base + '/' + madspin_card[0], mc_script))
                os.system("echo 'set nevents %d' >> %s" % (int(mclist[1]), mc_script))
                os.system("echo 'set iseed %d' >> %s" % (random.randint(1, 1000000), mc_script))

                with open('%s/run_job%d.sh' % (mc_job, int(i)), 'w') as bjob_script :

                    bjob_script.write(madgraph5_pythia8_plugin.produce_batch_script(mc_script))
                    bjob_script.write('\n')
                    if (options.FastSim) :
                        if (options.isNLO) :
                            bjob_script.write('gunzip %s/Events/run_01_decayed_1/events_PYTHIA8_0.hepmc.gz\n' % tmp_dir)
                            bjob_script.write("$HEPTOOL/delphes_install/bin/DelphesHepMC $VFS_PACKAGE_PATH/%s %s/delphes_output_%d.root %s/Events/run_01_decayed_1/events_PYTHIA8_0.hepmc\n\n" % (options.delphes_card, output, int(i), tmp_dir))
                        else :
                            bjob_script.write('gunzip %s/Events/run_01/tag_1_pythia8_events.hepmc.gz\n' % tmp_dir)
                            bjob_script.write("$HEPTOOL/delphes_install/bin/DelphesHepMC $VFS_PACKAGE_PATH/%s %s/delphes_output_%d.root %s/Events/run_01/tag_1_pythia8_events.hepmc\n\n" % (options.delphes_card, output, int(i), tmp_dir))
                    for filetype in ['txt', 'dat', 'log', 'gnuplot', 'HwU'] :
                        bjob_script.write('cp %s/Events/run_01/*.%s %s/\n' % (tmp_dir, filetype, outputdir))
                        bjob_script.write('cp %s/Events/run_01_decayed_1/*.%s %s/\n' % (tmp_dir, filetype, outputdir))
                    bjob_script.write('rm -rf %s' % tmp_dir)
                os.system('chmod 755 %s/run_job%d.sh' % (mc_job, int(i)))
                bjob_list.append('%s/run_job%d.sh' % (mc_job, int(i)))
        batch_job_submitter (bjob_list, options.queue)

if __name__ == '__main__' :
    sys.exit(bjobs_setup(sys.argv[1:]))

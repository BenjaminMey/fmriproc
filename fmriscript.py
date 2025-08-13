import os
import sys
import glob
import argparse
import fmriproc

timesubdir = r'TIME_FILES/'

def make_argument_parser():
    p = argparse.ArgumentParser(
        prog='fmriproc',
        usage='Given a directory with functional and anatomical data, this script will choose a reference volume per' \
              'functional task and affine align the other volumes stripping those that exceed the cost function. ' \
              'Then performing regression with the time and motion parameters, and then aligning the functional data ' \
              'to the anatomical data.',
        epilog='Script default taskname functional files: SYNONYMS_5_MIN_MB, VISUAL_5_MIN_MB, MOTOR_6_MIN_MB, ' \
               'VERBSSILENTGRADIENTMB, STORIESSILENTGRADIENTMB.'
    )
    p.add_argument('--version',
                   action='version',
                   version='%(prog) v 1.0')
    p.add_argument('-d', '--datadir',
                   nargs=1,
                   help='Provide data directory.',
                   default='/data')
    p.add_argument('-a', '--anat',
                   nargs=1,
                   help="Provide anatomical image name, default is 'SAGT13DMPRAGE.nii'",
                   default='SAGT13DMPRAGE.nii')
    p.add_argument('-t', '--taskname',
                   nargs='*',
                   help='Add space separated taskname(s) to process, default is all.') #accomplished in main
    p.add_argument('-c', '--timecoding',
                   nargs='*',
                   help='Add space separated time course regression codings for each taskname to process, default for built-in tasks')
    p.add_argument('--test',
                   action='store_true',
                   help='Test workflow on provided dataset')
    return p

def main():
    try:
        parser = make_argument_parser()
        args = parser.parse_args()

        timedir = os.path.join(os.getcwd(), timesubdir)
        if args.test:
            datadir = os.path.join(os.getcwd(), r'test/')
            anatname = 'SAG_T1_3D_MPRAGE.nii'
        else:
            datadir = args.datadir[0]
            anatname = args.anat[0]

        #Change to io directory
        os.chdir(os.path.abspath(os.path.expanduser(datadir)))

        if args.taskname is None: #run all processing for standard tasknames
            # processall function does not seem to refer to all tasknames, but to all steps of processing
            fmriproc.processall('SYNONYMS_5_MIN_MB', timedir + '5MinutesMB.time', anatname)
            fmriproc.processall('VISUAL_5_MIN_MB', timedir + '5MinutesMB.time', anatname)

            fmriproc.processmotor('MOTOR_6_MIN_MB', timedir + 'LeftHandMB.time', timedir + 'RightHandMB.time',
                                  args.anat)

            #TODO acquire Verbs log example for testing
            #logfile = glob.glob("*Verbs*.log")
            #logfile = logfile[0]
            #fmriproc.verbslogtotime(logfile, 'Verbs.time')
            #fmriproc.hushprocessall('VERBSSILENTGRADIENTMB', 'Verbs.time', 4, anatname)

            fmriproc.hushprocessall('STORIESSILENTGRADIENTMB', timedir + 'StoriesMB.time', 3, anatname)
        else:
            for i, task in enumerate(args.taskname):
                match task:
                    case 'SYNONYMS5MINMB':
                        fmriproc.processall('SYNONYMS_5_MIN_MB', timedir + '5MinutesMB.time', anatname)
                    case 'VISUAL5MINMB':
                        fmriproc.processall('VISUAL_5_MIN_MB', timedir + '5MinutesMB.time', anatname)
                    case 'MOTOR6MINMB':
                        fmriproc.processmotor('MOTOR_6_MIN_MB', timedir + 'LeftHandMB.time', timedir + 'RightHandMB.time', anatname)
                    case 'VERBSSILENTGRADIENTMB':
                        logfile = glob.glob("*Verbs*.log")
                        logfile = logfile[0]
                        fmriproc.verbslogtotime(logfile, 'Verbs.time')
                        fmriproc.hushprocessall('VERBSSILENTGRADIENTMB', 'Verbs.time', 4, anatname)
                    case 'STORIESSILENTGRADIENTMB':
                        fmriproc.hushprocessall('STORIESSILENTGRADIENTMB', timedir + 'StoriesMB.time', 3, anatname)
                    case _:
                        #Other tasks with matching time course codings processing. For now not needed.
                        #fmriproc.processall(task, args.timecoding[i], anatname)
                        raise ValueError('Taskname %s unknown, no associated time course files' % (task))
    except KeyboardInterrupt:
        print('\nKeyboardInterrupt detected. Ending script.')
        # As of 202508, all file opens are handled with(open, 'mode') as file:, so should exit gracefully.
    finally:
        print('Script finished.')

if __name__ == '__main__':
    sys.exit(main())

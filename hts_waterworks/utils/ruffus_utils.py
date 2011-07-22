
from optparse import OptionParser
import logging
import logging.handlers
import os
import sys
from collections import defaultdict
import subprocess
from time import ctime

from ruffus import pipeline_printout, pipeline_printout_graph, pipeline_run
from ruffus.proxy_logger import make_shared_logger_and_proxy


def ruffus_opt_parser():
    'creates a ruffus optparse opts object'
    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %progs [options]")
    #   general options: verbosity / logging
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", default=3,
                      help="Print more verbose messages for each additional verbose level.")
    #   pipeline
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        default=4,
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously. (%default by default)")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", default=False,
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")
    #   Less common pipeline options
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true", default=False,
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", '-f', dest="forced_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")
    parser.add_option("--config_file", '-c', dest="config_file",
                        default='pipeline.cfg',
                        type="string",
                        help="Configuration file for the pipeline")
    return parser

class DefaultLog:
    log_file = 'ruffus.log'
    verbose = 4

def ruffus_logger(options=None, module_name='pipeline'):
    'creates a shared logger and mutex'
    if options is None:
        options = DefaultLog()
    logger = logging.getLogger(module_name)
    _setup_std_logging(logger, options.log_file, options.verbose)
    def get_logger (logger_name, args):
        return logger
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger, module_name, {})
    logger_proxy.log_file = options.log_file
    return logger_proxy, logging_mutex

def _setup_std_logging (logger, log_file, verbose):
    """
    set up logging using programme options
    """
    class NullHandler(logging.Handler):
        """
        for when there is no logging
        """
        def emit(self, record):
            pass
    # We are interested in all messages
    logger.setLevel(logging.DEBUG)
    has_handler = False

    # log to file if that is specified
    if log_file:
        handler = logging.FileHandler(log_file, delay=False)
        handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
        handler.setLevel(logging.DEBUG)
        logger.addHandler(handler)
        has_handler = True

    # log to stderr if verbose
    if verbose:
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
        stderrhandler.setLevel(logging.DEBUG)
        logger.addHandler(stderrhandler)
        has_handler = True

    # no logging
    if not has_handler:
        logger.addHandler(NullHandler())

main_logger, main_mutex = ruffus_logger()

def ruffus_main(options, args):
    'Main entry point for ruffus pipelines'
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            verbose=options.verbose)
    elif options.flowchart:
        pipeline_printout_graph (   open(options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    options.target_tasks,
                                    options.forced_tasks,
                                    no_key_legend   = not options.key_legend_in_graph)
    else:
        pipeline_run(options.target_tasks, options.forced_tasks,
                            multiprocess    = options.jobs,
                            logger          = main_logger,
                            verbose         = options.verbose)


def sys_call(cmd, logger=main_logger, log_mutex=main_mutex, file_log=True):
    """Fork a job and write the results to the log file for this process
    
    stderr and stdout from the forked job are available in LOGFILE.PID
    """
    logfile = '%s.%s.log' % (logger.log_file, os.getpid())
    if file_log:
        cmd += ' 2>&1 | tee -a %s ' % logfile
    with open(logfile, 'a') as logout:
        logout.write('\n****** %s ******\n%s\n\n' % (ctime(), cmd))
    main_logger.debug(cmd)
    process = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    output, _ = process.communicate()
    
    retcode = process.poll()
    if retcode:
        with log_mutex:
            logger.error(output)
        raise subprocess.CalledProcessError(retcode, cmd)
    else:
        with log_mutex:
            logger.debug(output)

def touch(fname, times = None):
    'touch the given file (update timestamp or create the file)'
    with file(fname, 'a'):
        os.utime(fname, times)

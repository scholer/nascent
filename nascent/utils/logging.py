# -*- coding: utf-8 -*-
##    Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
##
##    This file is part of Nascent.
##
##    Nascent is free software: you can redistribute it and/or modify
##    it under the terms of the GNU Affero General Public License as
##    published by the Free Software Foundation, either version 3 of the
##    License, or (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU Affero General Public License for more details.
##
##    You should have received a copy of the GNU Affero General Public License
##    along with this program. If not, see <http://www.gnu.org/licenses/>.

# pylint: disable=W0142,C0103,C0301,W0141

"""

Module for setting up logging.

"""

from __future__ import absolute_import, print_function
import os
import logging
import logging.handlers
logger = logging.getLogger(__name__)




def init_logging(args=None, logfilepath=None, logdir=None):
    """
    Set up standard logging system based on parameters in args, e.g. loglevel and testing.
    :args:          dict with keys to configure logging setup: 'loglevel', 'basic_logging', 'rotating', etc.
    :logfilepath:   filename to output logging output to.
    :logdir:        If logfilepath is not provided, generate a filename in this directory.
                    If neither logfilepath or logdir is provided, will try to file the default
                    Nascent application home/data root directory and place the logfile here.
    After initializing/setting up the logging system it can be further configured through direct access
    to logging.root properties, e.g. logging.root.handlers.

    """
    if args is None:
        args = {}
    if logfilepath is None:
        appname = "Nascent"
        if logdir is None:
            try:
                import appdirs
                logdir = appdirs.user_log_dir(appname)
            except ImportError:
                if os.environ.get('APPDATA'):
                    logdir = os.path.join(os.environ['APPDATA'], appname, "Logs")
                elif sys.platform == 'darwin':
                    logdir = os.path.join(os.path.expanduser("~"), "Library", "Logs", appname)
                else:
                    logdir = os.path.join(os.path.expanduser("~"), "."+appname, "logs")
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        logfilepath = os.path.join(logdir, appname+".log")

    ## We want different output formatting for file vs console logging output.
    ## File logs should be simple and easy to regex; console logs should be short and nice on the eyes
    logfilefmt = '%(asctime)s %(levelname)-6s - %(name)s:%(lineno)s - %(funcName)s() - %(message)s'
    logdatefmt = "%Y%m%d-%H:%M:%S.%f"
    loguserfmt = "%(asctime)s %(levelname)-5s %(name)30s:%(lineno)-4s%(funcName)16s() %(message)s"
    loguserfmt = "%(asctime)s %(levelname)-5s %(module)30s:%(lineno)-4s%(funcName)16s() %(message)s"
    logtimefmt = "%H:%M:%S" # Nice for output to user in console and testing.
    # See https://docs.python.org/3/library/logging.html#logrecord-attributes for full list of attributes

    # Loglevel (for console messages)
    if args.get('loglevel'):
        try:
            loglevel = int(args['loglevel'])
        except (TypeError, ValueError):
            loglevel = getattr(logging, args['loglevel'].upper())
    else:
        loglevel = logging.DEBUG if args.get('testing') else logging.WARNING

    if args.get('basic_logging', False):
        logging.basicConfig(level=loglevel,
                            format=loguserfmt,
                            datefmt=logtimefmt,
                            filename=logfilename)
        logger.debug("Logging system initialized with loglevel %s", loglevel)
    else:

        # Set up custom logger:
        logging.root.setLevel(logging.DEBUG)  # Ensure that root logger accepts all DEBUG messages.

        # Add a rotating file handler:
        if args.get('rotating', False):
            logfilehandler = logging.handlers.RotatingFileHandler(logfilepath, maxBytes=2*2**20, backupCount=2)
        else:
            logfilehandler = logging.FileHandler(logfilepath)
        logfileformatter = logging.Formatter(fmt=logfilefmt, datefmt=logdatefmt)
        logfilehandler.setFormatter(logfileformatter)
        logging.root.addHandler(logfilehandler)
        print("Logging to file:", logfilepath)

        # Add a custom StreamHandler for outputting to the console (default level is 0 = ANY)
        logstreamhandler = logging.StreamHandler() # default stream is sys.stderr
        logging.root.addHandler(logstreamhandler)
        logstreamformatter = logging.Formatter(loguserfmt, logtimefmt)
        logstreamhandler.setFormatter(logstreamformatter)

        # Set filter for debugging:
        if args.get('debug_modules'):
            debug_modules = args['debug_modules']
            def module_debug_filter(record):
                """
                All Filters attached to a logger or handler are asked.
                The record is discarted if any of the attached Filters return False.
                """
                return any(record.name.startswith(modstr) for modstr in args['debug_modules']) \
                    or record.levelno >= loglevel
            logstreamhandler.addFilter(module_debug_filter)
            # Default level is 0, which is appropriate when using module_debug_filter
        else:
            # only set a min level if we are not using module_debug_filter. (Level is an additional filter.)
            logstreamhandler.setLevel(loglevel)
    logger.info("Logging system initialized...")

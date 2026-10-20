#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import json
import sys
import os
import atexit
import resource
import platform

from checkgms_utils import *
from checkgms_stable import *

def print_rusage():
    rusage = resource.getrusage(resource.RUSAGE_SELF)
    print()
    print("checkgms resource usage:")
    print()
    print("real     %12.3f sec" % (rusage.ru_utime + rusage.ru_stime,))
    print("user     %12.3f sec" % (rusage.ru_utime,))
    print("sys      %12.3f sec" % (rusage.ru_stime,))
    print("maxrss   %12d KiByte" % (rusage.ru_maxrss,))

"""
If --json_create flag is passed to trigger the creation of new *.json validation files then it will look or log
files containg the extension specified by file_extension to parse in order to generate the new *.json validation file.

TODO - make file_extension a flag in checkgms.py

"""
file_extension = ".VALIDATION"

"""
In some cases we may have a *.json validation file that we wish to protect from being over-written even when the --json_create
flag is passed in. A way to allow this is by renaming the file so that it ends with what is specified by json_protect.

Meaning if --json_protect is passed in then the script will skip all input files that match the filename pattern *project.json
"""
json_protect = "protect.json"

run_arguments = parse_arguments()

if run_arguments["mode"] == "none":
    print("Exiting on failure!")
    print(
        "Need to specify testing mode, -m {validation,regression,both}, --mode {validation,regression,both}")
    sys.exit(1)

if int(run_arguments["ncpus"]) > 1:
    if len(run_arguments["filter_folder"]) < 1:
        run_arguments["filter_folder"] = "parallel"
    else:
        run_arguments["filter_folder"] += ",parallel"

# Populate the log_file_paths array
log_file_paths = []
if run_arguments["json_create"]:
    log_file_paths = get_log_file_paths(
        filepath_string_match=run_arguments["filter_filepath"],
        folder_string_match=run_arguments["filter_folder"],
        file_string_match=run_arguments["filter_file"],
        folder_string_skip=run_arguments["skip_folder"],
        file_string_skip=run_arguments["skip_file"],
        file_extension=file_extension,
        test_path=run_arguments["test_path"])
else:
    log_file_paths = get_log_file_paths(
        filepath_string_match=run_arguments["filter_filepath"],
        folder_string_match=run_arguments["filter_folder"],
        file_string_match=run_arguments["filter_file"],
        folder_string_skip=run_arguments["skip_folder"],
        file_string_skip=run_arguments["skip_file"],
        test_path=run_arguments["test_path"])

if run_arguments["print_resource_usage"] and "Windows" not in platform.system():
    atexit.register(print_rusage)

# Loop through the log_file_paths array and validate
for filenum, log_file_path in enumerate(log_file_paths, start=1):

    if run_arguments["json_create"]:
        validation_result = checkgms(
            filenum=filenum,
            log_file_path=log_file_path,
            log_file_count=len(log_file_paths),
            run_arguments=run_arguments,
            file_extension=file_extension,
            json_protect=json_protect)
    else:
        validation_result = checkgms(
            filenum=filenum,
            log_file_path=log_file_path,
            log_file_count=len(log_file_paths),
            run_arguments=run_arguments,
        )

    if "skip" in str(validation_result):
        continue

    if "pass" not in str(validation_result):
        if run_arguments["fail_rename"]:
            os.rename(
                log_file_path,
                log_file_path.replace(
                    ".log",
                    run_arguments["fail_extension"]))
    else:
        if run_arguments["pass_delete"]:
            os.remove(log_file_path)

    if "empty" in str(validation_result):
        print("Empty Log File",log_file_path)

    if "unknown" in str(validation_result):
        print("Parser Exception",log_file_path)

    if run_arguments["exit_on_fail"]:
        if "pass" not in str(validation_result):
            if not run_arguments["no_dump"]:
                os.system("cat " + log_file_path)
            print("Exiting on failure!")
            print("File", log_file_path)
            sys.exit(1)

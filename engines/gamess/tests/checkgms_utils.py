import argparse
import json
import os
import sys
import datetime


class ColorX:
    reset = '\033[0m'
    bold = '\033[01m'
    disable = '\033[02m'
    underline = '\033[04m'
    reverse = '\033[07m'
    strikethrough = '\033[09m'
    invisible = '\033[08m'
    black = '\033[30m'
    red = '\033[31m'
    green = '\033[32m'
    orange = '\033[33m'
    blue = '\033[34m'
    purple = '\033[35m'
    cyan = '\033[36m'
    lightgrey = '\033[37m'
    darkgrey = '\033[90m'
    lightred = '\033[91m'
    lightgreen = '\033[92m'
    yellow = '\033[93m'
    lightblue = '\033[94m'
    pink = '\033[95m'
    lightcyan = '\033[96m'

# Rename to Color if you want clean text without the color markup.


class Color:
    reset = ''
    bold = ''
    disable = ''
    underline = ''
    reverse = ''
    strikethrough = ''
    invisible = ''
    black = ''
    red = ''
    green = ''
    orange = ''
    blue = ''
    purple = ''
    cyan = ''
    lightgrey = ''
    darkgrey = ''
    lightred = ''
    lightgreen = ''
    yellow = ''
    lightblue = ''
    pink = ''
    lightcyan = ''

# Formatted printing for numerical values with adjustable precision.
# Print width = 20
# Must be able to typecast input_string into a float or int.


def n_(input_string, precision=10):
    if precision > 0:
        return '{0:20.{1}f}'.format(
            float(input_string), precision
        )
    else:
        return '{0:20}'.format(
            int(input_string)
        )

# Formatted printing for numerical values with adjustable precision.
# Print width = 10
# Must be able to typecast input_string into a float or int.


def n_small(input_string, precision=10):
    if precision > 0:
        return '{0:10.{1}f}'.format(
            float(input_string), precision
        )
    else:
        return '{0:20}'.format(
            int(input_string)
        )

# Formatted printing with left-justification
# Print width = 20


def ls_(input_string, precision=10):
    return '{0:<20}'.format(
        str(input_string)
    )

# Formatted printing with center-justification
# Print width = 20


def cs_(input_string, precision=10):
    return '{0:^20}'.format(
        str(input_string)
    )

# Formatted printing with right-justification
# Print width = 20


def rs_(input_string, precision=10):
    return '{0:>20}'.format(
        str(input_string)
    )

# Formated printing with left-justification
# Print width = 40
# Color = light grey


def l_(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:<40}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing with center-justification
# Print width = 40
# Color = light grey


def c_(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:^40}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing with right-justification
# Print width = 40
# Color = light grey


def r_(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:>40}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with left-justification
# Print width = 20
# Color = light grey


def l_small(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:<20}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with center-justification
# Print width = 20
# Color = light grey


def c_small(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:^20}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with right-justification
# Print width = 20
# Color = light grey


def r_small(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:>20}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing with left-justification
# Print width = 60
# Color = light grey


def l_large(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:<60}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing with center-justification
# Print width = 60
# Color = light grey


def c_large(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:^60}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing with right-justification
# Print width = 60
# Color = light grey


def r_large(*args):
    input_string = ' '.join(map(str, args))
    return '{0}{1:>60}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing with left-justification
# Print width = custom width w/ default of 84
# Color = light grey


def l_custom(*args, width=84):
    input_string = ' '.join(map(str, args))
    return '{0}{1:<{3}}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
        width
    )

# Formated printing with center-justification
# Print width = custom width w/ default of 84
# Color = light grey


def c_custom(*args, width=84):
    input_string = ' '.join(map(str, args))
    return '{0}{1:^{3}}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
        width
    )

# Formated printing with right-justification
# Print width = custom width w/ default of 84
# Color = light grey


def r_custom(*args, width=84):
    input_string = ' '.join(map(str, args))
    return '{0}{1:>{3}}{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
        width
    )

# Formated printing within [ ] with left-justification
# Print width = 42
# Color = light grey


def l_box(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:<40}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with center-justification
# Print width = 40 (42 if you include '[' ']')
# Color = light grey


def c_box(*args):
    """
    Color for general values
    """
    input_string = ' '.join(map(str, args))
    return '{0}[{1:^40}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with right-justification
# Print width = 42
# Color = light grey


def r_box(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:>40}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with left-justification
# Print width = 62
# Color = light grey


def l_box_large(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:<60}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with center-justification
# Print width = 62
# Color = light grey


def c_box_large(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:^60}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with right-justification
# Print width = 62
# Color = light grey


def r_box_large(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:>60}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with left-justification
# Print width = custom width w/ default of 84
# Color = light grey


def l_box_custom(*args, width=84):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:<{3}}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
        width
    )

# Formated printing within [ ] with center-justification
# Print width = custom width w/ default of 84
# Color = light grey


def c_box_custom(*args, width=84):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:^{3}}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
        width
    )

# Formated printing within [ ] with right-justification
# Print width = custom width w/ default of 84
# Color = light grey


def r_box_custom(*args, width=84):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:>{3}}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
        width
    )

# Formated printing within [ ] with left-justification
# Print width = 22
# Color = light grey


def l_box_small(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:<20}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with center-justification
# Print width = 22
# Color = light grey


def c_box_small(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:^20}]{2}'.format(
        Color.lightgrey,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with right-justification
# Print width = 22
#Color = red


def r_box_small(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:>20}]{2}'.format(
        Color.red,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with left-justification
# Print width = 42
# Color = light blue


def lb_box(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:<40}]{2}'.format(
        Color.lightblue,
        input_string,
        Color.reset,
    )

# Formated printing within [ ] with left-justification
# Print width = 42
# Color = light red


def lr_box(*args):
    input_string = ' '.join(map(str, args))
    return '{0}[{1:<40}]{2}'.format(
        Color.lightred,
        input_string,
        Color.reset,
    )

# Formated printing of file progress [#### / ####]
# Print width = 13
#Color = default


def file_progress(*args):
    width = max(len(str(args[0])), len(str(args[1])))
    return '{0}[{1:>{4}} / {2:>{4}}]{3}'.format(
        Color.lightblue,
        int(args[0]),
        int(args[1]),
        Color.reset,
        width,
    )

# Formated printing of [Passed]
# Print width = 8
# Color = light green


def pass_box():
    return '{0}[Passed]{1}'.format(
        Color.lightgreen,
        Color.reset,
    )

# Formated printing of [Skipped]
# Print width = 10
# Color = light red


def skip_box():
    return '{0}[Skipped]{1}'.format(
        Color.lightred,
        Color.reset,
    )

# Formated printing of [Failed]
# Print width = 8
# Color = light red


def fail_box():
    return '{0}[Failed]{1}'.format(
        Color.lightred,
        Color.reset,
    )

# Formated printing of validation header
# Print width = 93
#Color = default


def validation_header():
    return '{0:^42} {1:^20} {2:^20} {3:^8}'.format(
        "Parameter",
        "Value",
        "Error",
        "Result"
    )

# Formated printing of "="*84
# Print width = 84
# Color = light grey


def seperator(columns=84):
    return '{0}'.format(Color.lightgrey) + '=' * \
        columns + '{0}'.format(Color.reset)

# Formated printing of "-"*84
# Print width = 84
# Color = light grey


def seperator_(columns=84):
    return '{0}'.format(Color.lightgrey) + '-' * \
        columns + '{0}'.format(Color.reset)

# Parse input argument


def parse_arguments(checkgms=True):
    run_arguments = {}

    if checkgms:
        # Command-line arguments used by $GMS_DIR/tests/checkgms.py
        parser = argparse.ArgumentParser(description='GAMESS Test Validation')
        parser.add_argument(
            '-r',
            '--regression',
            help='performance regression testing',
            action="store_true")
        parser.add_argument(
            '-m',
            '--mode',
            help='testing mode, default: validation',
            default="validation",
            choices=[
                'validation',
                'regression',
                'both'])
        parser.add_argument(
            '-a',
            '--arch',
            help='architecture string used for performance regression',
            default=None)
        parser.add_argument(
            '-s',
            '--system',
            help='HPC system string used for performance regression',
            default=None)
        parser.add_argument(
            '-n',
            '--ncpus',
            help='number of GAMESS compute processes for performance regression',
            default="1")
        parser.add_argument(
            '-g',
            '--naccelerators',
            help='number of accelerators for performance regression',
            default="0")
        parser.add_argument(
            '--dryrun',
            help='cycles through filelist without parsing',
            action="store_true")
        parser.add_argument(
            '--file',
            help='process file(s) containing substring',
            default="")
        parser.add_argument(
            '--folder',
            help='process folder(s) containing substring',
            default="")
        parser.add_argument(
            '--filepath',
            help='process filepath(s) containing substring',
            default="")
        parser.add_argument(
            '--test_path',
            help='process file(s) in the specified path',
            default=".")
        parser.add_argument(
            '--json_create',
            help='force the creation of JSON validation files',
            action="store_true")
        parser.add_argument(
            '--array',
            help='print out array values',
            action="store_true")
        parser.add_argument(
            '-d',
            '--debug',
            help='debug print control',
            action="store_true")
        parser.add_argument(
            '-i',
            '--print_to_stdout',
            help='print_to_stdout',
            action="store_true")
        parser.add_argument(
            '-e',
            '--exit_on_fail',
            help='exit on first failed validation',
            action="store_true")
        parser.add_argument(
            '--group',
            help='print group header for values',
            action="store_true")
        parser.add_argument(
            '-p',
            '--verbose_parsing',
            help='verbose printing during parsing',
            action="store_true")
        parser.add_argument(
            '-v',
            '--verbose_validation',
            help='verbose printing during validation',
            action="store_true")
        parser.add_argument(
            '--skip_file',
            help='skip file(s) containing substring',
            default="")
        parser.add_argument(
            '--skip_folder',
            help='skip folder(s) containing substring',
            default="")
        parser.add_argument(
            '--skip_json_create',
            help='skip creation of new JSON validation files',
            action="store_true")
        parser.add_argument(
            '--test_type',
            help='test input type: small, medium, large, msucc',
            default="")
        parser.add_argument(
            '--fail_rename',
            help='if validation fail, replaces ".log" with arg value',
            default="")
        parser.add_argument(
            '--pass_delete',
            help='delete logfile if validation passes',
            action="store_true")
        parser.add_argument(
            '--no_dump',
            help='suppress logfile output when exiting on failure',
            action="store_true")
        parser.add_argument(
            '--supress',
            help='limit output',
            action="store_true")
        parser.add_argument(
            '--print_resource_usage',
            help='on Linux systems print resource usage after validation',
            action="store_true")
    else:
        # Command-line arguments used by $GMS_DIR/tests/runtests.py
        parser = argparse.ArgumentParser(description='GAMESS Test Launch')
        parser.add_argument(
            '--pre',
            help='command to run before each test',
            default="")
        parser.add_argument(
            '--post',
            help='command to run after each test',
            default="")
        parser.add_argument(
            '--dryrun',
            help='cycles through filelist without parsing',
            action="store_true")
        parser.add_argument(
            '--no_counter',
            help='suppress the file counter',
            action="store_true")
        parser.add_argument(
            '-i',
            '--print_to_stdout',
            help='print_to_stdout',
            action="store_true")
        parser.add_argument(
            '--file',
            help='process file(s) containing substring',
            default="")
        parser.add_argument(
            '--folder',
            help='process folder(s) containing substring',
            default="")
        parser.add_argument(
            '--filepath',
            help='process filepath(s) containing substring',
            default="")
        parser.add_argument(
            '--path',
            help='Use a different GAMESS path',
            default="")
        parser.add_argument(
            '--version',
            help='Use a GAMESS version',
            default="00")
        parser.add_argument(
            '-d',
            '--debug',
            help='debug print control',
            action="store_true")
        parser.add_argument(
            '-n',
            '--ncpus',
            help='number of GAMESS compute processes',
            default="1")
        parser.add_argument(
            '-t',
            '--threads',
            help='number of parallel GAMESS computing',
            default="1")
        parser.add_argument('--mpi', help='use mpi', action="store_true")
        parser.add_argument('--hpe-cray-ex', help='use hpe-cray-ex', action="store_true")
        parser.add_argument('--hpe-cray-cs', help='use hpe-cray-cs', action="store_true")
        parser.add_argument('--cray-xt', help='use cray-xt', action="store_true")
        parser.add_argument('--cray-xc', help='use cray-xc', action="store_true")
        parser.add_argument(
            '--output_extension',
            help='extension to use for output files default(".log")',
            default=".log")
        parser.add_argument(
            '-q',
            '--queue',
            help='specify SLURM queue',
            default="compute")
        parser.add_argument(
            '--skip_file',
            help='skip file(s) containing substring',
            default="")
        parser.add_argument(
            '--skip_folder',
            help='skip folder(s) containing substring',
            default="")
        parser.add_argument(
            '--skip_log',
            help='skip file if log file exists',
            action="store_true")
        parser.add_argument(
            '--no_skip',
            help='ignore TRAVIS-CI SKIP tags during input filtering',
            action="store_true")
        parser.add_argument(
            '--test_type',
            help='test input type: small, medium, large, msucc',
            default="")
        parser.add_argument(
            '--stderr',
            help='print to stderr',
            action="store_true")
        parser.add_argument(
            '-c',
            '--coverage',
            help='run code coverage with GNU gcov',
            action="store_true")
        parser.add_argument(
            '--breakdown',
            help='breakdown coverage information per input',
            action="store_true")
        parser.add_argument(
            '--no_accumulate',
            help='do not accumulate coverage statistics across inputs',
            action="store_true")
        parser.add_argument(
            '--sourcefile',
            help='limit coverage data to source file substring, only valid with --breakdown',
            default="")
        parser.add_argument(
            '--bwrap',
            help='execute tests in a bwrap sandbox to avoid leftover semaphores',
            action="store_true")
        parser.add_argument(
            '--supress',
            help='limit output',
            action="store_true")

    args = parser.parse_args()

    # Used by both $GMS_DIR/tests/checkgms.py and $GMS_DIR/tests/runtests.py
    run_arguments["dryrun"] = args.dryrun
    run_arguments["filter_filepath"] = args.filepath
    run_arguments["filter_file"] = args.file
    run_arguments["filter_folder"] = args.folder
    run_arguments["debug"] = args.debug
    run_arguments["skip_file"] = args.skip_file
    run_arguments["skip_folder"] = args.skip_folder
    run_arguments["test_type"] = args.test_type
    run_arguments["ncpus"] = args.ncpus
    run_arguments["supress"] = args.supress

    # Used by $GMS_DIR/tests/runtests.py
    if not checkgms:
        run_arguments["pre"] = args.pre
        run_arguments["post"] = args.post
        run_arguments["no_counter"] = args.no_counter
        run_arguments["path"] = args.path
        run_arguments["version"] = args.version
        run_arguments["coverage"] = args.coverage
        run_arguments["print_to_stdout"] = args.print_to_stdout
        run_arguments["breakdown"] = args.breakdown
        run_arguments["no_accumulate"] = args.no_accumulate
        run_arguments["sourcefile"] = args.sourcefile
        run_arguments["ncpus"] = args.ncpus
        run_arguments["threads"] = args.threads
        run_arguments["output_extension"] = args.output_extension
        run_arguments["mpi"] = args.mpi
        run_arguments["hpe-cray-ex"] = args.hpe_cray_ex
        run_arguments["hpe-cray-cs"] = args.hpe_cray_cs
        run_arguments["cray-xt"] = args.cray_xt
        run_arguments["cray-xc"] = args.cray_xc
        run_arguments["stderr"] = args.stderr
        run_arguments["skip_log"] = args.skip_log
        run_arguments["no_skip"] = args.no_skip
        run_arguments["queue"] = args.queue
        run_arguments["pre"] = args.pre
        run_arguments["post"] = args.post
        run_arguments["bwrap"] = args.bwrap

        return run_arguments

    # Used by $GMS_DIR/tests/checkgms.py
    run_arguments["mode"] = args.mode
    run_arguments["arch"] = args.arch
    run_arguments["system"] = args.system
    run_arguments["naccelerators"] = args.naccelerators
    run_arguments["array"] = args.array
    run_arguments["json_create"] = args.json_create
    run_arguments["exit_on_fail"] = args.exit_on_fail
    run_arguments["group"] = args.group
    run_arguments["verbose_parsing"] = args.verbose_parsing
    run_arguments["verbose_validation"] = args.verbose_validation
    run_arguments["skip_json_create"] = args.skip_json_create
    run_arguments["test_path"] = args.test_path
    if len(args.fail_rename) > 0:
        run_arguments["fail_rename"] = True
        run_arguments["fail_extension"] = args.fail_rename
    else:
        run_arguments["fail_rename"] = False
    run_arguments["pass_delete"] = args.pass_delete
    run_arguments["no_dump"] = args.no_dump
    run_arguments["print_resource_usage"] = args.print_resource_usage

    # Raise exception if certain flags are missing
    if run_arguments["mode"] == "regression" or run_arguments["mode"] == "both":
        if run_arguments["arch"] is None:
            get_available_regression_baselines()
            parser.error(
                "-a, --arch is required when running in 'regression' mode")
        if run_arguments["system"] is None:
            get_available_regression_baselines()
            parser.error(
                "-s, --system is required when running in 'regression' mode")

    return run_arguments

# Prints the available regression baselines based on the prescences of
# *.info in $GMS_DIR/info-regression-baseline


def get_available_regression_baselines():
    print("Available regression baselines:\n")
    for (directory_path, directory_name, directory_files) in os.walk(
            "info-regression-baseline"):
        directory_files.sort()
        for file_name in directory_files:
            if ".info" in file_name:
                # file_name=<arch>_<system>_<ncpus>-<naccelerators>.info
                baseline = file_name.split("_")
                # baseline=[<arch>,<system>,<ncpus>-<naccelerators>.info]
                baseline[-1] = baseline[-1].split(".")[0]
                # baseline=[<arch>,<system>,<ncpus>-<naccelerators>]
                baseline.append(baseline[-1].split("-")[-1])
                # baseline=[<arch>,<system>,<ncpus>-<naccelerators>,<naccelerators>]
                baseline[-2] = baseline[-2].split("-")[0]
                # baseline=[<arch>,<system>,<ncpus>,<naccelerators>]
                print("\t--arch={} --system={} --ncpus={} --naccelerators={}".format(
                    baseline[0], baseline[1], baseline[2], baseline[3]))
    print()

# Returns an array containing file paths to log files as elements


def get_log_file_paths(filepath_string_match="",folder_string_match="", file_string_match="",
                       folder_string_skip="", file_string_skip="", file_extension=".log", test_path=None):
    logFiles = []
    filepath_string_match_array = filepath_string_match.split(",")
    folder_string_match_array = folder_string_match.split(",")
    file_string_match_array = file_string_match.split(",")
    folder_string_skip_array = folder_string_skip.split(",")
    file_string_skip_array = file_string_skip.split(",")
    if test_path is None:
        test_path = "."
    for (directory_path, directory_name,
         directory_files) in os.walk(test_path):
        directory_files.sort()
        for file_name in directory_files:
            # Skip if it is not a log file
            if file_extension not in file_name[-4:]:
                continue
            # Skip if filepath does not match when provided
            if len(filepath_string_match_array[0]) > 0:
                skip = False
                file_path = os.path.join(directory_path,file_name)
                for match_string in filepath_string_match_array:
                    if match_string[0] == ".":
                        match_string=match_string[1:]
                    if match_string not in file_path:
                        skip = True
                if skip:
                    continue
            else:
                # Skip folders containing folder_string_skip
                if len(folder_string_skip_array[0]) > 0:
                    skip = False
                    for skip_string in folder_string_skip_array:
                        if skip_string in directory_path:
                            skip = True
                    if skip:
                        continue
                # Skip files containing folder_string_skip
                if len(file_string_skip_array[0]) > 0:
                    skip = False
                    for skip_string in file_string_skip_array:
                        if skip_string in file_name:
                            skip = True
                    if skip:
                        continue
                # Skip folders not containing folder_string_match
                if len(folder_string_match_array[0]) > 0:
                    skip = False
                    for match_string in folder_string_match_array:
                        if match_string not in directory_path:
                            skip = True
                    if skip:
                        continue
                # Skip files not containing file_string_match
                if len(file_string_match_array[0]) > 0:
                    skip = False
                    for match_string in file_string_match_array:
                        if match_string not in file_name:
                            skip = True
                    if skip:
                        continue
                # Skip files not containing ".log"
                if file_extension not in file_name[-4:]:
                    continue
            logFiles.append(directory_path + "/" + file_name)

    return sorted(logFiles)

# Returns an array containing file paths to input files as elements


def get_input_file_paths(filepath_string_match="",folder_string_match="", file_string_match="",
                         folder_string_skip="", file_string_skip="", script_path="."):
    inputFiles = []
    filepath_string_match_array = filepath_string_match.split(",")
    folder_string_match_array = folder_string_match.split(",")
    file_string_match_array = file_string_match.split(",")
    folder_string_skip_array = folder_string_skip.split(",")
    file_string_skip_array = file_string_skip.split(",")
    for (directory_path, directory_name,
         directory_files) in os.walk(script_path):
        directory_files.sort()
        for file_name in directory_files:
       	    # Skip if it is not an input file
       	    if ".inp" not in file_name[-4:]:
       	       	continue
            # Skip if filepath does not match when provided
            if len(filepath_string_match_array[0]) > 0:
                skip = True
                file_path = os.path.join(directory_path,file_name)
                for match_string in filepath_string_match_array:
                    if match_string[0] == ".":
                        match_string=match_string[1:]
                    if match_string in file_path:
                        skip = False
                        break
                if skip:
                    continue
            else:
                # Skip folders containing folder_string_skip
                if len(folder_string_skip_array[0]) > 0:
                    skip = False
                    for skip_string in folder_string_skip_array:
                        if skip_string in directory_path.replace(script_path, ""):
                            skip = True
                           #print(skip_string, directory_path, "skip")
                    if skip:
                        continue
                # Skip files containing file_string_skip
                if len(file_string_skip_array[0]) > 0:
                    skip = False
                    for skip_string in file_string_skip_array:
                        if skip_string in file_name:
                            skip = True
                           #print(skip_string, file_name, "skip")
                        if skip_string in os.path.join(directory_path,file_name):
                            skip = True
                    if skip:
                        continue
                # Skip folders not containing folder_string_match
                if len(folder_string_match_array[0]) > 0:
                    skip = False
                    for match_string in folder_string_match_array:
                        if match_string not in directory_path.replace(
                                script_path, ""):
                            skip = True
                    if skip:
                        continue
                # Skip files not containing file_string_match
                if len(file_string_match_array[0]) > 0:
                    skip = False
                    for match_string in file_string_match_array:
                        if match_string not in file_name:
                            skip = True
                    if skip:
                        continue
                # Skip files not containing ".inp"
                if '.inp' not in file_name:
                    continue
            inputFiles.append(os.path.join(directory_path,file_name))

    return sorted(inputFiles)

# Creates a JSON for a JSON for a validation value


def create_JSON(value_name=None, value_header="Validation Values", value_value=None, value_type=None, value_precision=10,
                value_tolerance=None, create_mode=None, value_arch=None, value_system=None, value_ncpus=None, value_naccelerators=None):
    value = {}
    value["header"] = value_header
    value["name"] = value_name
    value["value"] = value_value
    value["type"] = value_type
    value["precision"] = value_precision
    if (value_type == "float") or (value_type == "float_array"):
        if value_tolerance is None:
            value["tolerance"] = 5.0 * 10**(-(value_precision - 2))
        else:
            value["tolerance"] = value_tolerance
    elif value_type == "integer":
        value["tolerance"] = 0
        value["precision"] = 0

    if create_mode == "regression":
        value["architecture"] = value_arch
        value["system"] = value_system
        value["ncpus"] = int(value_ncpus)
        value["naccelerators"] = int(value_naccelerators)
    return value

# Prints the JSON of a parsed log file


def print_parsed_JSON(parsed_json=None, run_arguments=None):
    if parsed_json is None:
        print(lr_box("parsed_json is undefined!"))
        return

    if run_arguments is None:
        print(lr_box("run_arguments is undefined!"))
        return

    current_header = ""
    for validation in parsed_json["validation"]:
        if validation["header"] != current_header:
            current_header = validation["header"]
            if run_arguments["group"]:
                print(lb_box(current_header))

        if hasattr(validation["value"], '__len__') and (
                not isinstance(validation["value"], str)):
            if ("array" in validation["type"]):
                if run_arguments["array"]:
                    print(r_box(validation["name"]))
                    for index, parsed_element_value in zip(
                            range(len(validation["value"])), validation["value"]):
                        print(
                            r_box(index), n_(
                                parsed_element_value, validation["precision"]))
                else:
                    print(r_box("LENGTH OF", validation["name"]), n_(
                        len(validation["value"]), 0))
                    print(
                        r_box(
                            "FIRST ELEMENT OF", validation["name"]), n_(
                            validation["value"][0], validation["precision"]))
                    print(r_box("LAST ELEMENT OF", validation["name"]), n_(
                        validation["value"][-1], validation["precision"]))
        else:
            print(
                r_box(
                    validation["name"]), n_(
                    validation["value"], validation["precision"]))

# Prints the JSON of a validated log file


def print_validated_JSON(validated_json=None, run_arguments=None):
    if validated_json is None:
        print(lr_box("validated_json is undefined!"))
        return

    if run_arguments is None:
        print(lr_box("run_arguments is undefined!"))
        return

    # For situations where validation entries did not match up - we skip
    # printing
    if validated_json["result"] == "skip":
        return

    current_header = ""
    for validation in validated_json["validation"]:
        if validation["header"] != current_header:
            current_header = validation["header"]
            if run_arguments["group"]:
                print(lb_box(current_header))

        if hasattr(validation["value"], '__len__') and (
                not isinstance(validation["value"], str)):
            if run_arguments["array"]:
                print(r_box(validation["name"]))
                for index, validated_element_value, validated_element_error in zip(
                        range(len(validation["value"])), validation["value"], validation["error"]):
                    print(
                        r_box(index), n_(
                            validated_element_value, validation["precision"]), n_(
                            validated_element_error, validation["precision"]))
            if validation["result"] == "pass":
                print(
                    r_box(
                        validation["name"]), rs_(""), n_(
                        validation["max_error"], validation["precision"]), pass_box())
            elif validation["result"] == "skip":
                print(
                    r_box(
                        validation["name"]), rs_(""), n_(
                        validation["max_error"], validation["precision"]), skip_box())
            else:
                print(
                    r_box(
                        validation["name"]), rs_(""), n_(
                        validation["max_error"], validation["precision"]), fail_box())

        else:
            if validation["result"] == "pass":
                print(
                    r_box(
                        validation["name"]), n_(
                        validation["value"], validation["precision"]), n_(
                        validation["error"], validation["precision"]), pass_box())
            elif validation["result"] == "skip":
                print(
                    r_box(
                        validation["name"]), n_(
                        validation["value"], validation["precision"]), n_(
                        validation["error"], validation["precision"]), skip_box())
            else:
                print(
                    r_box(
                        validation["name"]), n_(
                        validation["value"], validation["precision"]), n_(
                        validation["error"], validation["precision"]), fail_box())

    if validated_json["result"] == "pass":
        print(seperator(), pass_box())
    else:
        print(seperator(), fail_box())

# Get parser groups


def get_parse_groups(run_arguments=None, mode=None):
    if mode is None:
        print(lr_box("mode is undefined!"))
        return

    parse_groups = []

    # Read parser.inp file and write contents to parse_file_json
    parse_file_json = {}
    parse_input = os.path.join(os.path.realpath(sys.path[0]), "./parse.inp")

    if os.path.isfile(os.path.join(
            os.path.realpath(sys.path[0]), parse_input)):
        with open(os.path.join(os.path.realpath(sys.path[0]), parse_input), 'r', encoding="utf-8", errors='ignore') as parser_file:
            parse_file_json = json.load(parser_file)
    else:
        print(lr_box("parse.input file was not found!"))

    # Generate parser_group_file_paths array
    parser_group_file_paths = []
    parser_group_files = parse_file_json["parse_group_files"]

    try:
        for parser_group_file in parser_group_files:
            if os.path.isfile(os.path.join(
                    os.path.realpath(sys.path[0]), parser_group_file)):
                # Only files in the _parsegroups_ folder that DO NOT contain 'regression'
                # in the filename are used for validation testing
                if mode == "validation":
                    if "regression" not in parser_group_file:
                        parser_group_file_paths.append(
                            os.path.join(
                                os.path.realpath(
                                    sys.path[0]),
                                parser_group_file))

                # Only files in the _parsegroups_ folder that contain 'regression'
                # in the filename are used for regression testing
                if mode == "regression":
                    if "regression" in parser_group_file:
                        parser_group_file_paths.append(parser_group_file)
            else:
                if run_arguments["debug"]:
                    print(
                        lr_box("parse group file was not found!"),
                        os.path.join(
                            os.path.realpath(
                                sys.path[0]),
                            parser_group_file))
    except BaseException:
        print(lr_box("invalid JSON in parse.inp!"))
        print(lr_box("parser_group_files key was not found!"))

    # Load parse group and append to parse_groups array
    # Note that the order of the parse group in the parse.inp is preserved
    # Order matters!
    for parser_group_file in parser_group_file_paths:
        if run_arguments["debug"]:
            print(lb_box("Debug : ", parser_group_file))
        with open(os.path.join(os.path.realpath(sys.path[0]), parser_group_file), 'r', encoding="utf-8", errors='ignore') as parser_file:
            parse_file_json = json.load(parser_file)
            parse_groups.append(parse_file_json)

    return parse_groups

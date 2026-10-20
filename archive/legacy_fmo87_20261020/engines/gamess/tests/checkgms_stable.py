import json
import os

import mmap

from checkgms_utils import *
from checkgms_parsers import parser
from shutil import copyfile

# Driver for testing


def checkgms(filenum=None, log_file_path=None, log_file_count=0,
             run_arguments={}, file_extension=".log", json_protect=None):

    if log_file_path is None:
        print(lr_box("log_file_path is undefined!"))
        return

    validated_json = {}

    log_filename = log_file_path.split("/")[-1]
    base_filename = log_filename.replace(file_extension, "")

    with open(log_file_path, 'r', encoding="utf-8", errors='ignore') as opened_file:
        if run_arguments["verbose_parsing"] or run_arguments["debug"]:
            print(l_box("Parsing"), file_progress(
                filenum, log_file_count), l_(log_file_path))

        # Validation Testing
        if run_arguments["mode"] == "validation" or run_arguments["mode"] == "both":
            validation_file_path = log_file_path.replace(
                file_extension, ".json")
            validated_json = get_results(filenum=filenum, log_file_path=log_file_path, log_file_count=log_file_count, opened_file=opened_file, run_arguments=run_arguments,
                                         file_extension=file_extension, json_protect=json_protect, mode="validation", log_filename=log_filename, validation_file_path=validation_file_path)

        # Performance Regression Testing
        if run_arguments["mode"] == "regression" or run_arguments["mode"] == "both":
            regression_baseline_info = "{}_{}_{}-{}.info".format(
                run_arguments["arch"].replace(
                    " ",
                    "-"),
                run_arguments["system"].replace(
                    " ",
                    "-"),
                run_arguments["ncpus"],
                run_arguments["naccelerators"])
            regression_baseline_ldd = "{}_{}_{}-{}.ldd".format(
                run_arguments["arch"].replace(
                    " ",
                    "-"),
                run_arguments["system"].replace(
                    " ",
                    "-"),
                run_arguments["ncpus"],
                run_arguments["naccelerators"])
            regression_baseline_info_absolute_path = "info-regression-baseline/{}".format(
                regression_baseline_info)
            regression_baseline_ldd_absolute_path = "info-regression-baseline/{}".format(
                regression_baseline_ldd)

            if not os.path.isfile(regression_baseline_info_absolute_path):
                try:
                    copyfile(
                        "../install.info",
                        regression_baseline_info_absolute_path)
                except BaseException:
                    raise IOError(
                        "File does not exist: {}".format("../install.info"))

            if not os.path.isfile(regression_baseline_ldd_absolute_path):
                try:
                    os.system(
                        "ldd ../gamess.00.x > {}".format(regression_baseline_ldd_absolute_path))
                    if not os.path.isfile(
                            regression_baseline_ldd_absolute_path):
                        raise IOError(
                            "Shared object dependencies capture (ldd) not found!")
                except BaseException:
                    raise IOError(
                        "File does not exist for shared object dependencies capture (ldd): {}".format("../gamess.00.x"))

            validation_file_path = log_file_path.replace(
                file_extension, ".perf")
            validated_json = get_results(filenum=filenum, log_file_path=log_file_path, log_file_count=log_file_count, opened_file=opened_file, run_arguments=run_arguments,
                                         file_extension=file_extension, json_protect=json_protect, mode="regression", log_filename=log_filename, validation_file_path=validation_file_path)

    return validated_json


# Return the results of either validation or performance regression testing

def get_results(filenum=None, log_file_path=None, log_file_count=0, opened_file=None, run_arguments={
}, file_extension=None, json_protect=None, mode=None, log_filename=None, validation_file_path=None):
    if filenum is None:
        print(lr_box("filenum is undefined!"))
        return

    if log_file_path is None:
        print(lr_box("log_file_path is undefined!"))
        return

    if opened_file is None:
        print(lr_box("opened_file is undefined!"))
        return

    if file_extension is None:
        print(lr_box("file_extension is undefined!"))
        return

    if mode is None:
        print(lr_box("mode is undefined!"))
        return

    if log_filename is None:
        print(lr_box("log_filename is undefined!"))
        return

    if validation_file_path is None:
        print(lr_box("log_filename is undefined!"))
        return

    # Populate the parser_groups array
    parse_groups = []
    parse_groups = get_parse_groups(run_arguments=run_arguments, mode=mode)

    parsed_json = {}
    validation_json = {}
    validated_json = {}

    if not run_arguments["dryrun"]:
        # Create a memory map of the log file for faster searching if a given
        # regex_string exists within the document
        try:
            parse_memory_map = mmap.mmap(
                opened_file.fileno(), 0, access=mmap.ACCESS_READ)
            # Parse the log file and return a JSON representing the parsed log file
            values_json = parse(file_handle=opened_file, parse_memory_map=parse_memory_map,
                                run_arguments=run_arguments, parse_groups=parse_groups, mode=mode)
        except ValueError:
             validated_json["result"] = "empty"
             return validated_json["result"]
        except:
             validated_json["result"] = "unknown"
             return validated_json["result"]
    else:
        values_json = {}

    parsed_json["name"] = log_filename.replace(file_extension, ".log")
    parsed_json["validation"] = values_json

    if run_arguments["verbose_parsing"]:
        # Formatted printing of the parsed JSON
        print_parsed_JSON(parsed_json=parsed_json, run_arguments=run_arguments)

    if mode == "regression":
        regression_baseline_exist = False

    # Check if validation file exists. Yes, validate. No, save the parsed JSON
    # as a validation file.
    if os.path.isfile(
            validation_file_path) and not run_arguments["json_create"]:
        if run_arguments["verbose_validation"] or run_arguments["debug"]:
            print(l_box("Validating"), file_progress(
                filenum, log_file_count), l_(log_file_path))

        # Early return if its a dry run
        if run_arguments["dryrun"]:
            validated_json["result"] = "skip"
            return validated_json["result"]

        if not run_arguments["dryrun"]:
            with open(validation_file_path, 'r', encoding="utf-8", errors='ignore') as validation_file:
                validation_json = json.load(validation_file)

                if mode == "regression":
                    # Check if the regression baseline exist in the validation
                    # file
                    for i in range(len(validation_json["validation"])):
                        # If validation value name matches, proceed
                        if validation_json["validation"][i]["architecture"] == parsed_json["validation"][-1]["architecture"]:
                            if validation_json["validation"][i]["system"] == parsed_json["validation"][-1]["system"]:
                                if validation_json["validation"][i]["ncpus"] == parsed_json["validation"][-1]["ncpus"]:
                                    if validation_json["validation"][i]["naccelerators"] == parsed_json["validation"][-1]["naccelerators"]:
                                        regression_baseline_exist = True

                    # If the regression baseline does not exist for the parsed architecture + system + npcu combination
                    # then add it to the validation file
                    if not regression_baseline_exist:
                        validated_json["result"] = "skip"
                        print(l_box("Appending validation file"), file_progress(
                            filenum, log_file_count), l_(validation_file_path))
                        if not run_arguments["dryrun"]:
                            with open(validation_file_path, 'w', encoding="utf-8", errors='ignore') as validation_file:
                                # Append the regression baseline from the
                                # parsed log file to the existing validation
                                # JSON and save
                                validation_json["validation"].append(
                                    parsed_json["validation"][-1])
                                validation_file.write(
                                    json.dumps(validation_json, indent=2))
                        return validated_json["result"]

            # Validate parsed JSON and return a modified JSON containing
            # calculated errors and validation results
            # ("pass/fail/skip/validation")
            validated_json = validate(
                validation_json=validation_json, parsed_json=parsed_json, mode=mode, run_arguments=run_arguments)
            if run_arguments["verbose_validation"] and validated_json["result"] != "validation":
                # Formatted printing of the validated JSON
                print_validated_JSON(
                    validated_json=validated_json, run_arguments=run_arguments)

        # Simplified formatted printing for validation results
        if validated_json["result"] == "pass":
            print(l_box("Validation result"), file_progress(
                filenum, log_file_count), l_large(log_file_path), pass_box())
        else:
            print(l_box("Validation result"), file_progress(
                filenum, log_file_count), l_large(log_file_path), fail_box())
        return validated_json["result"]

    # Check if validation file exists AND if it is protected
    elif os.path.isfile(validation_file_path) and json_protect is not None and json_protect in validation_file_path:
        validated_json["result"] = "skip"
        print(l_box("Skipping validation file create"), file_progress(
            filenum, log_file_count), l_(validation_file_path))
        return validated_json["result"]
    # Check if we are skipping validation file creation
    elif run_arguments["skip_json_create"]:
        validated_json["result"] = "skip"
        print(l_box("Skipping validation file create"), file_progress(
            filenum, log_file_count), l_(validation_file_path))
        return validated_json["result"]
    else:
        validated_json["result"] = "skip"
        print(l_box("Creating validation file"), file_progress(
            filenum, log_file_count), l_(validation_file_path))
        if not run_arguments["dryrun"]:
            with open(validation_file_path, 'w', encoding="utf-8", errors='ignore') as validation_file:
                validation_file.write(json.dumps(parsed_json, indent=2))
        return validated_json["result"]


# Compare parsed JSON file (obtained from the current log file) with the
# validation JSON file (obtained from the initial/original/baseline log
# file)

def validate(validation_json=None, parsed_json=None, mode=None, run_arguments=None):
    if validation_json is None:
        print(lr_box("validation_json is undefined!"))
        return

    if parsed_json is None:
        print(lr_box("parsed_json is undefined!"))
        return

    if mode is None:
        print(lr_box("mode is undefined!"))
        return

    # If filename matches, proceed
    if validation_json["name"] == parsed_json["name"]:

        if mode == "regression":
            # Find regression baseline entry in validation JSON that matches the parsed JSON and update
            # validation_json["validation"] so that it ONLY contains this entry so that the next condition
            # that checks the lengths of validation_json["validation"] and parsed_json["validation"]
            # matches and the validation can proceed.
            #
            # We need to do this so that we can value multiple entries in validation_json["validation"] for
            # various architecture + system + ncpus combination.
            for i in range(len(validation_json["validation"])):
                # If validation value name matches, proceed
                if validation_json["validation"][i]["architecture"] == parsed_json["validation"][-1]["architecture"]:
                    if validation_json["validation"][i]["system"] == parsed_json["validation"][-1]["system"]:
                        if validation_json["validation"][i]["ncpus"] == parsed_json["validation"][-1]["ncpus"]:
                            if validation_json["validation"][i]["naccelerators"] == parsed_json["validation"][-1]["naccelerators"]:
                                validation_json["validation"] = [
                                    validation_json["validation"][i]]
                                break

        # If the number of validation values are the same, proceed
        if len(validation_json["validation"]) == len(
                parsed_json["validation"]):
            parsed_json["result"] = "pass"

            # Loop through validation values
            for i in range(len(parsed_json["validation"])):
                # If validation value name matches, proceed
                if validation_json["validation"][i]["name"] == parsed_json["validation"][i]["name"]:

                    # Using absolute values during error calculation. Allowing
                    # sign changes.
                    if hasattr(validation_json["validation"][i]["value"], '__len__') and (
                            not isinstance(validation_json["validation"][i]["value"], str)):
                        error = []
                        for validation_element_value, parsed_element_value in zip(
                                validation_json["validation"][i]["value"], parsed_json["validation"][i]["value"]):
                            element_error = abs(
                                float(validation_element_value)) - abs(float(parsed_element_value))
                            error.append(element_error)
                    else:
                        if validation_json["validation"][i]["type"] == "float" or validation_json["validation"][i]["type"] == "skip_float":
                            if mode == "validation":
                                # error = | baseline value | - | observed value
                                # |
                                error = abs(float(validation_json["validation"][i]["value"])) - abs(
                                    float(parsed_json["validation"][i]["value"]))

                            if mode == "regression":
                                # error = observed value - baseline value
                                error = float(parsed_json["validation"][i]["value"]) - float(
                                    validation_json["validation"][i]["value"])

                        elif validation_json["validation"][i]["type"] == "integer" or validation_json["validation"][i]["type"] == "skip_integer":
                            # error = | baseline value | - | observed value |
                            error = abs(int(validation_json["validation"][i]["value"])) - abs(
                                int(parsed_json["validation"][i]["value"]))

                    parsed_json["validation"][i]["error"] = error

                    # Find largest absolute if we are dealing with an array of
                    # values and set that as our value for error
                    if hasattr(validation_json["validation"][i]["value"], '__len__') and (
                            not isinstance(validation_json["validation"][i]["value"], str)):
                        min_error = min(error)
                        max_error = max(error)
                        if abs(min_error) > abs(max_error):
                            error = min_error
                        else:
                            error = max_error
                        parsed_json["validation"][i]["max_error"] = error

                    if mode == "validation":
                        if "skip" in validation_json["validation"][i]["type"]:
                            parsed_json["validation"][i]["result"] = "skip"
                        elif abs(
                                error) <= validation_json["validation"][i]["tolerance"]:
                            parsed_json["validation"][i]["result"] = "pass"
                        else:
                            parsed_json["validation"][i]["result"] = "fail"
                            parsed_json["result"] = "fail"

                    if mode == "regression":
                        # if error is less than tolerance % then we pass else,
                        # we fail
                        if error <= float(validation_json["validation"][i]["value"]) * (
                                validation_json["validation"][i]["tolerance"] / 100.00):
                            parsed_json["validation"][i]["result"] = "pass"
                        else:
                            parsed_json["validation"][i]["result"] = "fail"
                            parsed_json["result"] = "fail"
                else:
                    print(lr_box("Validation entries do not match!"),
                          validation_json["validation"][i]["name"], 'vs.', parsed_json["validation"][i]["name"])
                    parsed_json["validation"][i]["result"] = "fail"
                    parsed_json["result"] = "validation"

        else:
            if (run_arguments["supress"] == False):
                print(lr_box("Number of validations does not match!"), len(
                    validation_json["validation"]), 'vs.', len(parsed_json["validation"]))
            validation_headers = []
            parsed_headers = []
            if run_arguments["verbose_validation"]:
                for i in validation_json["validation"]:
                    validation_headers.append((i['header'], i['name']))
                for i in parsed_json["validation"]:
                    parsed_headers.append((i['header'], i['name']))
                for i in list(set(parsed_headers) - set(validation_headers)):
                    print(lr_box("Found in .log, but not in .json!"), i[1])
                for i in list(set(validation_headers) - set(parsed_headers)):
                    print(lr_box("Found in .json, but not in .log!"), i[1])
            parsed_json["result"] = "validation"

    else:
        print(lr_box("Validating incorrect files!"),
              validation_json["name"], 'vs.', parsed_json["name"])
        parsed_json["result"] = "validation"

    return parsed_json


# Parses the log file and returns a JSON of containing parsed values
# defined by the parse groups in parse.inp

def parse(file_handle=None, parse_memory_map=None,
          run_arguments=None, parse_groups=None, mode=None):
    if mode is None:
        print(lr_box("mode is undefined!"))
        return

    parsed_values = []

    for parse_group_json in parse_groups:
        # Populate parse header or resort to parser group filename
        try:
            parse_header = parse_group_json["parse_header"]
        except BaseException:
            parse_header = "Parse Group"

        if run_arguments["debug"]:
            print(lb_box("Debug : ", parse_header))

        # Loop over parse group
        for parse in parse_group_json["parse_group"]:

            # Process parse JSON
            try:
                regex_string = parse["regex_string"]
            except BaseException:
                print(lr_box("regex_string was not found!"))
                continue

            try:
                parse_name = parse["parse_name"]
            except BaseException:
                parse_name = "Value Name"

            try:
                parse_case_sensitive = parse["parse_case_sensitive"]
            except BaseException:
                parse_case_sensitive = True

            try:
                parse_type = parse["parse_type"]
            except BaseException:
                parse_type = "float"

            parse_integer = False
            if parse_type == "integer":
                parse_integer = True

            try:
                parse_precision = parse["parse_precision"]
            except BaseException:
                parse_precision = 10

            try:
                parse_tolerance = parse["parse_tolerance"]
            except BaseException:
                parse_tolerance = None

            try:
                parse_value_index = parse["parse_value_index"]
            except BaseException:
                parse_value_index = 0

            try:
                parse_line_offset = parse["parse_line_offset"]
            except BaseException:
                parse_line_offset = 0

            try:
                parse_instance = parse["parse_instance"]
            except BaseException:
                parse_instance = -1

            try:
                parse_required = parse["parse_required"]
            except BaseException:
                parse_required = False

            try:
                parse_recipe = parse["parse_recipe"]
            except BaseException:
                parse_recipe = None

            if parse_recipe is not None:
                if parse_recipe == "lastinstance_firstindex":
                    parse_instance = -1
                    parse_value_index = 0
                elif parse_recipe == "lastinstance_lastindex":
                    parse_instance = -1
                    parse_value_index = -1
                elif parse_recipe == "firstinstance_firstindex":
                    parse_instance = 1
                    parse_value_index = 0
                elif parse_recipe == "firstinstance_lastindex":
                    parse_instance = 1
                    parse_value_index = -1
                else:
                    print(lr_box("parse_recipe does not currenty exist!"))
                    continue

            # Set parse_multiline to True if we are parsing an array
            parse_multiline = False
            if "array" in parse_type:
                parse_multiline = True

            parse_memory_map_temp = parse_memory_map

            parse_multiple = False
            try:
                parse_multiple_indexes = parse["parse_multiple_indexes"]
                parse_multiple_label = parse["parse_multiple_label"]
                parse_multiple = True
            except BaseException:
                parse_multiple_instances = None

            try:
                parse_multiple_instances = parse["parse_multiple_instances"]
                parse_multiple_label = parse["parse_multiple_label"]
                parse_multiple = True
            except BaseException:
                parse_multiple_instances = None

            try:
                parse_multiple_lines = parse["parse_multiple_lines"]
                parse_multiple_label = parse["parse_multiple_label"]
                parse_multiple = True
            except BaseException:
                parse_multiple_lines = None

            if parse_multiple:

                # Parse multiple instances
                if parse_multiple_instances is not None:
                    # Last N instances
                    if parse_multiple_instances < 0:
                        # Loop over parse instances
                        for instance in range(parse_multiple_instances, 0):
                            parse_name_instance = parse_name + " : " + \
                                parse_multiple_label + " " + str((instance))
                            parse_instance = instance
                            file_handle.seek(0)
                            parsed_value = parser(
                                regex_string, file_handle,
                                parse_memory_map=parse_memory_map_temp,
                                parse_value_index=parse_value_index,
                                parse_line_offset=parse_line_offset,
                                parse_instance=parse_instance,
                                parse_integer=parse_integer,
                                parse_case_sensitive=parse_case_sensitive,
                                parse_multiline=parse_multiline,
                                parse_mode=mode)

                            if parsed_value is not None:

                                # Calculate precision
                                if parse_precision < 0:
                                    rhs = parsed_value.split(".")[1]
                                    parse_precision = len(rhs)

                                # Calculate tolerance
                                if parse_tolerance is not None:
                                    if parse_tolerance < 0:
                                        lhs = parsed_value.replace(
                                            "-", "").split(".")[0]
                                        rhs = parsed_value.replace(
                                            "-", "").split(".")[1]
                                        if len(lhs) > 4:
                                            parse_tolerance = 5.0 * \
                                                10**(len(lhs) - 12)
                                        else:
                                            parse_tolerance = 5.0 * \
                                                10**(-(len(rhs) - 2))

                                parsed_values.append(
                                    create_JSON(
                                        value_name=parse_name_instance,
                                        value_value=parsed_value,
                                        value_header=parse_header,
                                        value_type=parse_type,
                                        value_precision=parse_precision,
                                        value_tolerance=parse_tolerance,
                                        create_mode=mode,
                                        value_arch=run_arguments["arch"],
                                        value_system=run_arguments["system"],
                                        value_ncpus=run_arguments["ncpus"],
                                        value_naccelerators=run_arguments["naccelerators"]))
                            else:
                                # We don't break because we are looping in
                                # reverse
                                pass

                        if parsed_value is None:
                            if parse_required:
                                # Break out of the for loop over parse group
                                break

                    # First N instances
                    if parse_multiple_instances > 0:
                        for instance in range(1, parse_multiple_instances + 1):
                            parse_name_instance = parse_name + " : " + \
                                parse_multiple_label + " " + str(instance)
                            parse_instance = instance
                            file_handle.seek(0)
                            parsed_value = parser(
                                regex_string, file_handle,
                                parse_memory_map=parse_memory_map_temp,
                                parse_value_index=parse_value_index,
                                parse_line_offset=parse_line_offset,
                                parse_instance=parse_instance,
                                parse_integer=parse_integer,
                                parse_case_sensitive=parse_case_sensitive,
                                parse_mode=mode)

                            if parsed_value is not None:

                                # Calculate precision
                                if parse_precision < 0:
                                    rhs = parsed_value.split(".")[1]
                                    parse_precision = len(rhs)

                                # Calculate tolerance
                                if parse_tolerance is not None:
                                    if parse_tolerance < 0:
                                        lhs = parsed_value.replace(
                                            "-", "").split(".")[0]
                                        rhs = parsed_value.replace(
                                            "-", "").split(".")[1]
                                        if len(lhs) > 4:
                                            parse_tolerance = 5.0 * \
                                                10**(len(lhs) - 12)
                                        else:
                                            parse_tolerance = 5.0 * \
                                                10**(-(len(rhs) - 2))

                                parsed_values.append(
                                    create_JSON(
                                        value_name=parse_name_instance,
                                        value_value=parsed_value,
                                        value_header=parse_header,
                                        value_type=parse_type,
                                        value_precision=parse_precision,
                                        value_tolerance=parse_tolerance,
                                        create_mode=mode,
                                        value_arch=run_arguments["arch"],
                                        value_system=run_arguments["system"],
                                        value_ncpus=run_arguments["ncpus"],
                                        value_naccelerators=run_arguments["naccelerators"]))
                            else:
                                # Break out of the for loop over parse
                                # instances
                                break

                        if parsed_value is None:
                            if parse_required:
                                # Break out of the for loop over parse group
                                break

                elif parse_multiple_lines is not None:
                    for line in range(0, parse_multiple_lines):
                        parse_name_line = parse_name + " : " + \
                            parse_multiple_label + " " + str(line + 1)
                        parse_line_offset_ = parse_line_offset + line
                        file_handle.seek(0)
                        parsed_value = parser(
                            regex_string, file_handle,
                            parse_memory_map=parse_memory_map_temp,
                            parse_value_index=parse_value_index,
                            parse_line_offset=parse_line_offset_,
                            parse_instance=parse_instance,
                            parse_integer=parse_integer,
                            parse_case_sensitive=parse_case_sensitive,
                            parse_multiline=parse_multiline,
                            parse_mode=mode)

                        if parsed_value is not None:

                            # Calculate precision
                            if parse_precision < 0:
                                rhs = parsed_value.split(".")[1]
                                parse_precision = len(rhs)

                            # Calculate tolerance
                            if parse_tolerance is not None:
                                if parse_tolerance < 0:
                                    lhs = parsed_value.replace(
                                        "-", "").split(".")[0]
                                    rhs = parsed_value.replace(
                                        "-", "").split(".")[1]
                                    if len(lhs) > 4:
                                        parse_tolerance = 5.0 * \
                                            10**(len(lhs) - 12)
                                    else:
                                        parse_tolerance = 5.0 * \
                                            10**(-(len(rhs) - 2))

                            parsed_values.append(
                                create_JSON(
                                    value_name=parse_name_line,
                                    value_value=parsed_value,
                                    value_header=parse_header,
                                    value_type=parse_type,
                                    value_precision=parse_precision,
                                    value_tolerance=parse_tolerance,
                                    create_mode=mode,
                                    value_arch=run_arguments["arch"],
                                    value_system=run_arguments["system"],
                                    value_ncpus=run_arguments["ncpus"],
                                    value_naccelerators=run_arguments["naccelerators"]))
                        else:
                            # Break out of the for loop over parse instances
                            break

                    if parsed_value is None:
                        if parse_required:
                            # Break out of the for loop over parse group
                            break

                elif parse_multiple_indexes is not None:
                    for parse_value_index in range(
                            0, parse_multiple_indexes + 1):
                        parse_name_line = parse_name + " " + \
                            parse_multiple_label + " " + \
                            str(parse_value_index + 1)
                        file_handle.seek(0)
                        parsed_value = parser(
                            regex_string, file_handle,
                            parse_memory_map=parse_memory_map_temp,
                            parse_value_index=parse_value_index,
                            parse_line_offset=parse_line_offset,
                            parse_instance=parse_instance,
                            parse_integer=parse_integer,
                            parse_case_sensitive=parse_case_sensitive,
                            parse_multiline=parse_multiline,
                            parse_mode=mode)

                        if parsed_value is not None:

                            # Calculate precision
                            if parse_precision < 0:
                                rhs = parsed_value.split(".")[1]
                                parse_precision = len(rhs)

                            # Calculate tolerance
                            if parse_tolerance is not None:
                                if parse_tolerance < 0:
                                    lhs = parsed_value.replace(
                                        "-", "").split(".")[0]
                                    rhs = parsed_value.replace(
                                        "-", "").split(".")[1]
                                    if len(lhs) > 4:
                                        parse_tolerance = 5.0 * \
                                            10**(len(lhs) - 12)
                                    else:
                                        parse_tolerance = 5.0 * \
                                            10**(-(len(rhs) - 2))

                            parsed_values.append(
                                create_JSON(
                                    value_name=parse_name_line,
                                    value_value=parsed_value,
                                    value_header=parse_header,
                                    value_type=parse_type,
                                    value_precision=parse_precision,
                                    value_tolerance=parse_tolerance,
                                    create_mode=mode,
                                    value_arch=run_arguments["arch"],
                                    value_system=run_arguments["system"],
                                    value_ncpus=run_arguments["ncpus"],
                                    value_naccelerators=run_arguments["naccelerators"]))
                        else:
                            # Break out of the for loop over parse instances
                            break

                    if parsed_value is None:
                        if parse_required:
                            # Break out of the for loop over parse group
                            break

            else:
                file_handle.seek(0)
                parsed_value = parser(
                    regex_string, file_handle,
                    parse_memory_map=parse_memory_map_temp,
                    parse_value_index=parse_value_index,
                    parse_line_offset=parse_line_offset,
                    parse_instance=parse_instance,
                    parse_integer=parse_integer,
                    parse_case_sensitive=parse_case_sensitive,
                    parse_multiline=parse_multiline,
                    parse_mode=mode)

                if parsed_value is not None:

                    # Calculate precision
                    if parse_precision < 0:
                        rhs = parsed_value.split(".")[1]
                        parse_precision = len(rhs)

                    # Calculate tolerance
                    if parse_tolerance is not None:

                        if parse_tolerance < 0:

                            # If we have an array
                            if hasattr(parsed_value, '__len__') and (
                                    not isinstance(parsed_value, str)):
                                lhs = parsed_value[0].replace(
                                    "-", "").split(".")[0]
                                rhs = parsed_value[0].replace(
                                    "-", "").split(".")[1]
                            # Else we have a scalar
                            else:
                                lhs = parsed_value.replace(
                                    "-", "").split(".")[0]
                                rhs = parsed_value.replace(
                                    "-", "").split(".")[1]

                            if len(lhs) > 4:
                                parse_tolerance = 5.0 * 10**(len(lhs) - 12)
                            else:
                                parse_tolerance = 5.0 * 10**(-(len(rhs) - 2))

                    parsed_values.append(
                        create_JSON(
                            value_name=parse_name,
                            value_value=parsed_value,
                            value_header=parse_header,
                            value_type=parse_type,
                            value_precision=parse_precision,
                            value_tolerance=parse_tolerance,
                            create_mode=mode,
                            value_arch=run_arguments["arch"],
                            value_system=run_arguments["system"],
                            value_ncpus=run_arguments["ncpus"],
                            value_naccelerators=run_arguments["naccelerators"]))
                else:
                    if parse_required:
                        # Break out of the for loop over parse group
                        break

    return parsed_values

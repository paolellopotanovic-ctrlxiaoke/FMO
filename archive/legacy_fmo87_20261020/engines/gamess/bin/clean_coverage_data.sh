#!/bin/bash
find object -maxdepth 1 -type f -writable \( -name '*.gcda' -o -name '*.gcov' \) -delete

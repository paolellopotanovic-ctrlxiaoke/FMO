#!/usr/bin/python3

import re

def srcline(src, line):
  global srcmap
  return "{0}#L{1}".format(srcmap[src],line)

data = open("warnings.gcc").read().split("\n")
srcmaptxt = open("map.txt").read().split("\n")
new_lines = open("new-lines.txt").read().split("\n")

res = []
out = []
srcmap = {}
resmap = {}

for line in srcmaptxt:
  if len(line) == 0:
    break
  line = line.split(" ")
  srcmap[line[0]] = line[1]

# 0 - search file
# 1 - search warning-error
state = 0

for line in data:
  if "No such file or directory" in line:
    break
  if "file not found" in line:
    break
  if state == 0:
    res.append(line.split(":"))
    state = 1
  elif state == 1:
    if len(re.findall("^Warning", line)) == 0 and len(re.findall("^Error", line)) == 0:
      continue
    res[-1].append(line)
    tres = []
    state = 0

res = [ [srcline(x[0].split("/")[-1], x[1]), x[-1]] for x in res ]

for line in res:
  if "Conversion from HOLLERITH to REAL(8)" in line[1]:
    continue
  if "Legacy Extension: Hollerith constant at" in line[1]:
    continue
  if "Equality comparison for REAL(8) at (1)" in line[1]:
    continue
  if "Inequality comparison for REAL(8) at (1)" in line[1]:
    continue
  resmap[line[0]] = line[1]

out.append("Total warnings:    {0:6}".format(len(res)))
out.append("Selected warnings: {0:6}".format(len(resmap)))
for new in new_lines:
  try:
    val = resmap[new]
    out.append("{0} {1}".format(new, val))
  except:
    pass
out.insert(2, "New warnings:      {0:6}".format(len(out)-2))

open("all-warnings.txt", "w").write("\n".join([ " ".join(x) for x in res ]))
open("new-warnings.txt", "w").write("\n".join(out))

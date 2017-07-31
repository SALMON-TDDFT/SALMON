#!/usr/bin/env python
#
#   Copyright 2017 ARTED developers
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
import os

os.chdir(os.path.dirname(__file__))

with open("template.txt") as fh_temp:
    temp = fh_temp.read()

for arch in os.listdir("arch"):
    with open(os.path.join("arch", arch)) as fh_head:
        head = fh_head.read()
    
    with open("Makefile.%s" % arch, "w") as fh_make:
        fh_make.write("# %s\n\n" % arch)
        fh_make.write(head)
        fh_make.write(temp)
        print("# Generated: %s" % arch)

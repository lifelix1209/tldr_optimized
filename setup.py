#!/usr/bin/env python

"""
Legacy setup entrypoint.

Project metadata and dependencies are declared in pyproject.toml.
This file keeps optional runtime tool checks for users that still run:

    python setup.py install
"""

import os
import sys
import subprocess

from setuptools import setup


def check_python():
    return sys.version_info.major == 3 and sys.version_info.minor >= 9


def check_minimap2():
    p = subprocess.Popen(["minimap2", "--version"], stdout=subprocess.PIPE)
    for line in p.stdout:
        major, minor = line.decode().strip().split(".")
        if int(major) >= 2:
            return True
    return False


def check_samtools():
    p = subprocess.Popen(["samtools"], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()
        if line.startswith("Version:"):
            major, minor = line.strip().split()[1].split(".")[:2]
            minor = minor.split("-")[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_mafft():
    p = subprocess.Popen(["mafft", "--help"], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode().strip()
        if line.startswith("MAFFT"):
            return True
    return False


def check_exonerate():
    p = subprocess.Popen(["exonerate"], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode()
        if line.startswith("exonerate from exonerate"):
            major, minor = line.strip().split()[-1].split(".")[:2]
            minor = minor.split("-")[0]
            if int(major) >= 2 and int(minor) >= 2:
                return True
    return False


if __name__ == "__main__":
    if not check_python():
        sys.exit("Dependency problem: python >= 3.9 is required")

    # Optional runtime tool check:
    # TLDR_STRICT_INSTALL_CHECK=1 pip install .
    strict_install_check = os.environ.get("TLDR_STRICT_INSTALL_CHECK", "0") == "1"
    if strict_install_check:
        if not check_minimap2():
            sys.exit("Dependency problem: minimap2 >= 2.0 not found")
        if not check_samtools():
            sys.exit("Dependency problem: samtools >= 1.2 not found")
        if not check_mafft():
            sys.exit("Dependency problem: mafft not found")
        if not check_exonerate():
            sys.exit("Dependency problem: exonerate not found")

setup()

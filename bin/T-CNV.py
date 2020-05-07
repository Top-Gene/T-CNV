#!/usr/bin/env python3
# * coding: utf-8 *

import os,sys

if sys.version_info < (3,0):
	import platform
	raise Exception("T-CNV requires python3 version 3.0+ (version %s detected)" % (platform.python_version()))

scriptDir = os.path.abspath(os.path.dirname(__file__))
workflowDir = os.path.abspath(os.path.join(scriptDir,"../lib"))
sys.path.append(workflowDir)

from configureOptions import ConfigureWorkflowOptions
from runWorkflow import DetectByCNVect

def main():
	workflowOptions = ConfigureWorkflowOptions().get_run_options()

	CNVect = DetectByCNVect(workflowOptions)
	CNVect.runWorkflow()

if __name__ == "__main__" :
	main()

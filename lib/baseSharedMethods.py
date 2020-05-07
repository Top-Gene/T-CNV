import math
from collections import defaultdict

def computelog2(num):
	try:
		return math.log(num)/math.log(2)
	except ValueError as e:
		return False

def createDict():
	return defaultdict(createDict)

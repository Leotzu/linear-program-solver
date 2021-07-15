# Magic incantation:
#__import__('code').interact(local={k: v for ns in (globals(), locals()) for k, v in ns.items()})
import numpy as np
import sys
from scipy.optimize import linprog

def parse_input():
	# get input file name
	filename = sys.argv[-1]
	# open input file and read lines
	f = open(filename)
	lines = f.readlines()
	
	for iter, line in enumerate(lines):
		#print(line)
		if iter is 0:
			obj = np.array([float(x.strip()) for x in line.split()])
			# obj *= -1 # this is only necessary for running the scipy's linprog
			num_vars = len(obj)
			num_constraints = len(lines)-1
			# declare lhs and rhs matrices with proper dimensions: 
			lhs = np.zeros(shape=(num_constraints, num_vars))
			rhs = np.zeros(shape=(num_constraints, 1))
			# initialize objective function variables: x1=1, x2=2, ..., xn=n
			obj_vars = np.array([x+1 for x in range(num_vars)], dtype=int)
			# initialize basis variables (lhs_vars): w1=-1, w2=-2, ... ,wn=-n
			lhs_vars = np.zeros(shape=(num_constraints, num_vars+1), dtype=int)	
			for w in range(num_constraints):
				lhs_vars[w] = np.array([x for x in range(num_vars+1)])
				lhs_vars[w,0] = -(w+1)
			continue
		# for every line in the txt file, fill in that row for both lhs and rhs (this is the efficient numpy way)
		lhs[iter-1] = np.array([float(x.strip()) for x in line.split()[:num_vars]])		
		rhs[iter-1] = np.array([float(x.strip()) for x in line.split()[num_vars:]])	
		
	'''
	# print out the parsed data for bug inspection	
	print("==================================")
	print("objective function:\n", obj)
	print("==================================")
	print("Left hand side:\n", lhs)
	print("==================================")
	print("Constants (rhs):\n", rhs)
	print("==================================")
	'''

	f.close()
	return obj, obj_vars, lhs, lhs_vars, rhs

# Determine if the dictionary is at an optimal solution
def is_optimal(obj):
	# If all vals in the objective function are <= 0, return true
	for i in obj:
		if i > 0:
			return False
	return True


# Checks if there are 2 zeros in the dictionary.
# This doesn't necessarily mean there will be a degenerate pivot, 
# but it's enough to trigger the lexicographic method to avoid cycling.
def is_degenerate(rhs):
	# if there are 2 zeros in rhs, return true
	zero_count = 0
	for coeff in rhs:
		if coeff == 0:
			zero_count += 1
		if zero_count == 2:
			return True
	return False

def lexicographic_method(obj, lhs, rhs):
	# TODO: how to do this symbolically?? or do I actually add these tiny values in??
	return rhs

# Checks for and handles degeneracy, unboundedness, infeasibility, and auxiliary method necessity.
# If none of the above: Find largest coefficient and pivot it via the Simplex Method. Repeat until optimal. 
# Return optimal dictionary form.
def largest_coeff_pivot(obj, obj_vars, obj_val, lhs, lhs_vars, rhs):
	# repeat pivots until optimal, unbounded, or infeasible. Return optimal or exit if the latters. 
	while not is_optimal(obj):
		
		# TODO: update rhs with lex method and continue pivoting
		# perhaps set a flag to note that the epsilon offsets will need to be removed before the final answer is printed. 
		if is_degenerate(rhs):
			#rhs = lexicographic_method(obj, lhs, rhs)
			# until worked out, just exit the program
			print('DEGENERATE DICTIONARY')
			exit(0)

		for b in rhs:
			if b < 0:
				print("NEED AUXILIARY TO FIND FEASIBLE POINT")
				exit(0)
		
		# TODO: LARGEST_COEFF_PIVOT: find largest coeff, pivot it, return 
		# choose the largest positive obj coeff. if all basis function coeffs in that col are > 0, then print_unbounded()
	
		# find entering variable:
		max_coeff = 0
		# initialized to -1 so that if no entering position is found, an error can be called
		entering_pos = -1
		for pos, elem in enumerate(obj):
			if elem > max_coeff:
				max_coeff = elem
				entering_pos = pos

		# find leaving variable:
		min_ratio = np.inf # numpy infinity
		leaving_pos = -1
		# a flag variable to check if LP is unbounded
		is_unbounded = False
		for pos in range(np.shape(lhs)[0]):
			# if this if true, then there exists a constraint s.t. the entering variable cannot go to infinity
			if lhs[pos, entering_pos] > 0:
				is_unbounded = False
				if rhs[pos] / lhs[pos,entering_pos] < min_ratio and rhs[pos] > 0:
					min_ratio = rhs[pos] / lhs[pos, entering_pos]
					leaving_pos = pos
		
		# if the LP is unbounded, print and exit entire program
		if is_unbounded:
			print('unbounded')
			exit(0)
		'''
		print(obj)
		print(lhs)
		print(rhs)
		print('entering_pos: ', entering_pos)
		print('leaving_pos: ', leaving_pos)
		print('current obj_val: ', obj_val)
		print('***************')
		'''
		# perform pivot (2 parts: 1. pivot variable placement; 2. update coefficients):
		# 1. pivot just the variable location arrays:
		temp_var = 0
		for pos in range(len(obj_vars)):
			if pos == entering_pos:
				temp_var = obj_vars[pos]
				obj_vars[pos] = lhs_vars[leaving_pos,0]
		for row_num, row in enumerate(lhs_vars):           #range(np.shape(lhs_vars)[0]):
			for col_num, val in enumerate(row):
				if col_num == 0 and row_num == leaving_pos:
					lhs_vars[row_num, col_num] = temp_var
				if col_num == entering_pos+1:
					lhs_vars[row_num, col_num] = obj_vars[entering_pos]
		
		# 2. update coefficients:
		
		# objective value
		obj_val += obj[entering_pos] * rhs[leaving_pos] / lhs[leaving_pos, entering_pos]	
	
		# this is a special coefficient that is used a lot, so it gets a special name
		pivot_coeff = lhs[leaving_pos, entering_pos]

		# objection function coefficients:
		prev_entering_coeff = obj[entering_pos]
		obj[entering_pos] /= -pivot_coeff
		for pos in range(len(obj)):
			if pos == entering_pos:
				#obj[pos] /= -pivot_coeff # DELETEME
				continue
			obj[pos] -= prev_entering_coeff * lhs[leaving_pos, pos] / pivot_coeff

		# lhs and rhs coefficients:
		for row_num, row in enumerate(lhs):
			for col_num, val in enumerate(row):
				if row_num == leaving_pos:
					if col_num == entering_pos:
						lhs[row_num, col_num] = 1 / pivot_coeff
						continue
					lhs[row_num, col_num] /= pivot_coeff
					continue
		rhs[leaving_pos] /= pivot_coeff
		# second pass through the matrix allows us to use the previously calculated coefficients
		for row_num, row in enumerate(lhs):
			# to ensure we don't change this value before we're done referencing it
			row_coeff = lhs[row_num, entering_pos]
			if row_num == leaving_pos:
				continue
			rhs[row_num] -= row_coeff * rhs[leaving_pos]
			for col_num, val in enumerate(row):
				if col_num == entering_pos:
					lhs[row_num, col_num] /= -pivot_coeff
					continue
				lhs[row_num, col_num] -= row_coeff * lhs[leaving_pos, col_num]
				rhs[row_num]
	'''
	print('final objective coefficients:\n', obj)
	print('final lhs:\n', lhs)
	print('final rhs:\n', rhs)
	print('final objective vars:\n', obj_vars)
	print('final basis vars:\n', lhs_vars)
	'''

	return obj, obj_vars, obj_val, lhs, lhs_vars, rhs

def use_linprog(obj_fun, lhs, rhs):
	# use linprog() to calculate solution
	result = linprog(c=obj_fun, A_ub=lhs, b_ub=rhs, method="simplex")
	
	# print out results
	if result.success:
		# Note: result.fun will be negative, since we solved the minimization problem, not max.
		optimal_soln = -result.fun
		#__import__('code').interact(local={k: v for ns in (globals(), locals()) for k, v in ns.items()})
		print("optimal")
		# %0.7G prints 7 sig figs, without including trailing zeros
		print('%0.7G' % optimal_soln)
		for count, val in enumerate(result.x):
			if count is not 0:
				print(' ', end='')
			print('%0.7G' % val, end='')
			#print(('%0.7f' % val).rstrip('0').rstrip('.'), end=' ')
		print()
	elif result.status == 2:
		print("infeasible")
	elif result.status == 3:
		print("unbounded")
	else:
		print("UNKNOWN ERROR")

def main():

	if len(sys.argv) != 2:
		print("Error: this program requires exactly 1 argument -> an input file containing the linear program to be solved.")
		print("-----> Example of correct program call: $ python LP_solver.py input_file.txt")
		return 0
	
	# parses input file into numpy arrays
	obj, obj_vars, lhs, lhs_vars, rhs = parse_input()	
	# initialize objective value to 0
	obj_val = 0
	
	'''
	# use scipy's linear program solver (solves 29/32 test cases) 	
	use_linprog(obj, lhs, rhs)
	'''
	
	obj, obj_vars, obj_val, lhs, lhs_vars, rhs = largest_coeff_pivot(obj, obj_vars, obj_val, lhs, lhs_vars, rhs) 
	
	# ********** THIS IS WHERE ALL OF largest_coeff_pivot() WAS *****************

	# PRINTING OPTIMAL SOLUTION TODO: Put this into a seperate function
	# Get optimal coefficients:
	optimal_coeffs = np.zeros(len(obj))
	for var_pos, marker in enumerate(lhs_vars[:,0]):
		# negative marker values correspond to w1, w2, ..., so skip them
		if marker < 0:
			continue
		for rhs_pos, rhs_val in enumerate(rhs):
			if var_pos == rhs_pos:
				optimal_coeffs[marker-1] = rhs_val
	# Print optimal output:
	print('optimal')
	# %0.7G prints 7 sig figs, without including trailing zeros
	print('%0.7G' % obj_val)
	for count, x in enumerate(optimal_coeffs):
		if count != 0:
			print(' ', end='')
		print('%0.7G' % x, end='')
	print()

	return 0


if __name__ == "__main__":
	main()

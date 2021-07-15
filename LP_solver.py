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
			obj_fun = np.array([float(x.strip()) for x in line.split()])
			# obj_fun *= -1 # this is necessary for running the scipy's linprog
			num_vars = len(obj_fun)
			num_constraints = len(lines)
			# declare lhs and rhs matrices with proper dimensions  
			lhs = np.zeros(shape=(num_constraints-1, num_vars))
			rhs = np.zeros(shape=(num_constraints-1, 1))
			continue

		# for every line in the txt file, fill in that row for both lhs and rhs (Apparently this is the efficient numpy way)
		lhs[iter-1] = np.array([float(x.strip()) for x in line.split()[:num_vars]])		
		rhs[iter-1] = np.array([float(x.strip()) for x in line.split()[num_vars:]])		
		
	#'''
	# print out the parsed data for bug inspection	
	print("==================================")
	print("objective function:\n", obj_fun)
	print("==================================")
	print("Left hand side:\n", lhs)
	print("==================================")
	print("Constants (rhs):\n", rhs)
	print("==================================")
	#'''
	f.close()
	return obj_fun, lhs, rhs

# Determine if the dictionary is at an optimal solution
def is_optimal(obj_fun):
	# If all vals in obj_fun are negative, return true
	for i in obj_fun:
		if i > 0:
			return False
	return True

def print_infeasible():
	# TODO: ???
	return False

def print_unbounded():
	# TODO: ???
	return False

def is_degenerate():
	# TODO: if there are 2 zeros in rhs, return true
	return False

def lexicographic_method(obj_fun, lhs, rhs):
	# TODO: how to do this symbolically??
	return rhs

def largest_coeff_pivot(obj_fun, lhs, rhs):
	# TODO: find largest coeff, pivot it, return 
	# choose the largest positive obj_fun coeff. if all basis function coeffs in that col are > 0, then print_unbounded()
	
	# find entering variable:
	max_var = 0
	# initialized to -1 so that if no entering var is found, an error can be called
	entering_var = -1
	for pos, elem in enumerate(obj_fun):
		if elem > max_entering:
			max_var = elem
			entering_var = pos

	# find leaving variable:
	min_ratio = np.inf # numpy infinity
	leaving_var = -1
	for pos, elem in enumerate(lhs):
		if elem > 0:
			if rhs[pos]/lhs[pos] < min_ratio:
				min_ratio = rhs[pos]/lhs[pos]
				leaving_var = pos
	
	# now perform pivot:
	# TODO:
	
	return obj_fun, lhs, rhs

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
	obj_fun, lhs, rhs = parse_input()
	
	'''
	# use scipy's linear program solver (solves 29/32 test cases) 	
	use_linprog(obj_fun, lhs, rhs)
	'''
	return 0
	
	while not is_optimal():
		''' # these are both detected while running largest_coeff_pivot()
		if is_infeasible():
			print('infeasible')
			return 0
		if is_unbounded():
			print('unbounded')
			return 0
		'''
		if is_degenerate():
			rhs = lexicographic_method(obj_fun, lhs, rhs)
		obj_fun, lhs, rhs = largest_coeff_pivot(obj_fun, lhs, rhs)

	return 0


if __name__ == "__main__":
	main()

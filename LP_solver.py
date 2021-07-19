# Magic debugging incantation:
#__import__('code').interact(local={k: v for ns in (globals(), locals()) for k, v in ns.items()})

import numpy as np
import sys
# used to verify results in the beginning
from scipy.optimize import linprog

def parse_input():
	
	#filename = sys.argv[-1]
	filename = "/dev/stdin"
	# open input file and read lines
	f = open(filename)
	lines = f.readlines()	

	for iter, line in enumerate(lines):
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
def is_optimal(obj, rhs):
	# If all vals in the objective function are <= 0 and 0,0,...,0 is feasible, return true
	for i in obj:
		if i > 0:
			return False
	for b in rhs:
		if b < 0:
			return False
	return True


# Checks if there are 2 zeros in the dictionary.
# This doesn't necessarily mean there will be a degenerate pivot, 
# but it's enough to trigger the lexicographic method to avoid possible future cycling.
def is_degenerate(rhs):
	# if there are 2 zeros in rhs, return true
	zero_count = 0
	for coeff in rhs:
		if coeff == 0:
			zero_count += 1
		if zero_count == 2:
			return True
	return False


def pivot_dictionary(obj, obj_vars, obj_val, lhs, lhs_vars, rhs, entering_pos, leaving_pos):
	# perform pivot (2 parts: 1. pivot variable placement; 2. update coefficients):
	
	# 1) pivot just the variable location arrays:
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
	
	# 2) update coefficients:
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
	return obj, obj_vars, obj_val, lhs, lhs_vars, rhs


# There are 2 parts: 1) determining if LP is feasible; and 2) finding feasible dictionary form and return this LP.
def auxiliary_method(obj, obj_vars, lhs, lhs_vars, rhs):
	print('AUXILIARY METHOD CALLED')
	# 1) check if the LP is feasible (does optimal Omega = 0?) 
	
	# create new aux_LP (a copy of the original which will be modified)
	aux_obj = np.copy(obj)
	aux_obj_vars = np.copy(obj_vars)
	aux_lhs = np.copy(lhs)
	aux_lhs_vars = np.copy(lhs_vars)
	aux_rhs = np.copy(rhs)
	# Find the row that contains the most negative aux_rhs value (this will be pivoted)
	most_neg_val = 0
	for row, b in enumerate(aux_rhs):
		if b < most_neg_val:
			most_neg_row = row
			most_neg_val = b
	# set all previous objective function coeffs to 0 in aux_LP
	for i in range(len(aux_obj)):
		aux_obj[i] = 0

	# Add a new variable/coeff (omega) to the end of each line in aux_obj, aux_obj_vars, aux_lhs, and aux_lhs_vars:
	aux_obj = np.append(aux_obj, -1)
	# add col of -1's to aux_lhs
	col_of_ones = np.ones((np.shape(aux_lhs)[0], 1))
	aux_lhs = np.hstack((aux_lhs, -col_of_ones)) 
	# omega=0 in the array of variable positions (recall: w_n=-n and x_n=1)
	aux_obj_vars = np.append(aux_obj_vars, 0)
	# add col of 0's to aux_lhs_vars (this represents the position of omega in lhs)
	col_of_zeros = np.zeros((np.shape(aux_lhs_vars)[0], 1))
	aux_lhs_vars = np.hstack((aux_lhs_vars, col_of_zeros))
	
	# now pivot, with entering=omega and leaving=(the variable in the most_neg_row)
	aux_entering_pos = len(aux_obj)-1
	aux_leaving_pos = most_neg_row
	aux_obj_val = 0
	# one-time pivot:
	aux_obj, aux_obj_vars, aux_obj_val, aux_lhs, aux_lhs_vars, aux_rhs = pivot_dictionary(aux_obj, aux_obj_vars, aux_obj_val, aux_lhs, aux_lhs_vars, aux_rhs, aux_entering_pos, aux_leaving_pos)
	# this aux_LP must now feasible, so it can be sent into largest_coeff_pivot():

	# Send aux_LP into largest_coeff_pivot() solve it and find the aux_obj_val
	#TODO: think about possible errors if this aux_LP is unbounded or something....
	'''
	print('LP BEING SENT INTO LARGEST_COEFF_PIVOT FROM WITHIN AUX_METHOD:')
	print('aux_obj:\n', aux_obj)
	print('aux_obj_val: ', aux_obj_val)
	print('aux_lhs:\n', aux_lhs)
	print('aux_rhs:\n', aux_rhs)
	print('aux_obj_vars:\n', aux_obj_vars)
	print('aux_lhs_vars:\n', aux_lhs_vars)
	print('==============================')
	'''
	
	aux_obj, aux_obj_vars, aux_obj_val, aux_lhs, aux_lhs_vars, aux_rhs = largest_coeff_pivot(aux_obj, aux_obj_vars, aux_obj_val, aux_lhs, aux_lhs_vars, aux_rhs)
	'''
	print('LP RETURNED TO AUX_METHOD FROM LARGEST_COEFF_PIVOT:')
	print('aux_obj:\n', aux_obj)
	print('aux_obj_val: ', aux_obj_val)
	print('aux_lhs:\n', aux_lhs)
	print('aux_rhs:\n', aux_rhs)
	print('aux_obj_vars:\n', aux_obj_vars)
	print('aux_lhs_vars:\n', aux_lhs_vars)
	print('==============================')
	'''


	# Check if original LP is feasible. Print infeasible and terminate program if not.
	if aux_obj_val != 0:
		print('AUXILIARY INFEASIBLE')
		print('infeasible')
		exit(0)

	# 2) Modify original LP so that it's feasible, then return it. 	
	# at this point, the original LP should not have been modified in any way, only the aux_LP

	# remove all the Omega variables (find the col to delete):
	col_to_delete = -1
	for row in aux_lhs_vars:
		for pos, val in enumerate(row):
			if val == 0:
				col_to_delete = pos
				break
		break
	# if omega is a basis variable, remove that row and keep all columns:
	if col_to_delete == -1:
		for row, var in enumerate(aux_lhs_vars[:,0]):
			if int(var) == 0:
				# delete this row instead of deleting a column
				lhs_vars = np.delete(aux_lhs_vars, row, 0)
				lhs = np.delete(aux_lhs, row, 0)
				break
	else:
		lhs_vars = np.delete(aux_lhs_vars, col_to_delete, 1)
		# there's one less col in lhs, so the Omega coeffs col are shifted over by one 
		lhs = np.delete(aux_lhs, col_to_delete-1, 1)
	
	# Set obj_vars to have the same form as the updated lhs_vars
	old_obj_vars = np.copy(obj_vars)
	obj_vars = lhs_vars[0,1:]
	
	
	# Update obj coeffs. Take original obj function, and sub in the new x1, x2, ... eqns.
	# Note: we know that the original obj_vars form was just 1,2,3,...,n, since aux_method() can only be called right at the beginning.
	old_obj_coeffs = np.copy(obj)
	# reset all the obj coeffs and sum them up while looping through lhs to collect the coeffs
	obj = np.zeros(len(old_obj_coeffs))
	# initialize obj_val in order to add to it when merging aux_LP and original_LP
	obj_val = 0
	nonbasic_x_vals = np.array([int(x) for x in range(len(obj))])
	# print('BEFORE**********nonbasic_x_vals:\n', nonbasic_x_vals)
	# loop through the first col of lhs_vars to find rows that looks like x_n = b_n + ...
	for row_num, var_num in enumerate(lhs_vars[:,0]):
		# if we're looking at an x1, x2, ..., or xn row
		if var_num > 0:
			# update obj_val and cross it off the list of x_n values that were in the basis (crossing off by setting to -1)
			nonbasic_x_vals[int(var_num-1)] = -1
			obj_val += aux_rhs[row_num] * old_obj_coeffs[int(var_num-1)]
			# update obj coeffs:
			for col_num in range(len(obj)):
				# ************ TODO warninig: it's possible that lhs[row,col] should be negative... I'm not sure ****************
				obj[col_num] -= old_obj_coeffs[int(var_num-1)] * lhs[row_num, col_num]
	
	# print('AFTER**********nonbasic_x_vals:\n', nonbasic_x_vals)
	
	# for all x_n values that weren't in the basis, just sub them back in as variables. 
	# BUT! Make sure you're adding the old_coeff to the old obj_vars variable, and not adding the old_coeff to the new obj_vars position.
	for col_num in range(len(obj)):
		if col_num in nonbasic_x_vals:
			# now find the old_obj_var associated with this col_num
			for new_pos, var in enumerate(obj_vars):
				if col_num+1 == int(var):
					#print('OBJ BEFORE:')
					#print(obj)
					# TODO: This could be += or -=..... Not sure. ******************
					obj[new_pos] += old_obj_coeffs[col_num]
					#print('OBJ AFTER:')
					#print(obj)

					

	# update rhs:
	rhs = np.copy(aux_rhs)

	# return original LP that's now organised to be feasible
	print('AUXILIARY COMPLETE')
	return obj, obj_vars, obj_val, lhs, lhs_vars, rhs


# Checks for and handles degeneracy, unboundedness, infeasibility, and auxiliary method necessity.
# If none of the above: Find largest coefficient and pivot it via the Simplex Method. Repeat until optimal. 
# Return optimal dictionary form.
def largest_coeff_pivot(obj, obj_vars, obj_val, lhs, lhs_vars, rhs):
	
	# set a cycle flag to false to check for cycling later	
	cycle_mode = False
	# repeat pivots until optimal, unbounded, or infeasible. Return optimal or exit if the latters. 
	while not is_optimal(obj, rhs):
		
		# Check if any rhs vals are negative (i.e. infeasible dict). If so, call auxiliary_method()
		# Note: if no val in rhs is negative to start, it never will become negative through future Simplex pivots
		for b in rhs:
			if b < 0:
				obj, obj_vars, obj_val, lhs, lhs_vars, rhs = auxiliary_method(obj, obj_vars, lhs, lhs_vars, rhs)
				'''
				print("LP COMING OUT OF AUXILIARY:")
				print('obj:\n', obj)
				print('obj_val: ', obj_val)
				print('lhs:\n', lhs)
				print('rhs:\n', rhs)
				print('obj_vars:\n', obj_vars)
				print('lhs_vars:\n', lhs_vars)
				'''
				break
		# need to check if the rearranged LP coming from the Aux_method is now optimal
		if is_optimal(obj, rhs):	
			break

		# Check if cycling is possible (if there exist 2 zeros vals in rhs):
		# if not in cycle_mode already, check if we should be
		if not cycle_mode:
			zero_count = 0
			for coeff in rhs:
				if coeff == 0:
					zero_count += 1
				if zero_count == 2:
					#print('ENTERING CYCLING_MODE')
					cycle_mode = True
					# set e_0=0, e_1=1, ... e_n=n for all epsilon offsets
					rhs_epsilon = np.array([int(e) for e in range(len(rhs))])	
					break

		# find entering variable:
		max_coeff = 0
		# initialized to -1 so that if no entering position is found, an error can be called
		entering_pos = -1
		for pos, elem in enumerate(obj):
			if elem > max_coeff:
				max_coeff = elem
				entering_pos = pos

		# find leaving variable:
		if cycle_mode:
			# Use Lexicographic method to choose leaving variable to avoid cycling:
			# 1) lex_method searches for a viable lhs coeff (i.e. > 0), then selects the one with the SMALLEST epsilon value (e_2 < e_1)
			# 2) during the pivot, e in rhs_epsilon must be moved accordingly. However, we only care about the LARGEST e value, so if during a pivot a smaller e value is added to a cell in rhs_epsilon that already contains an e value that's equal or greater to the new e, we ignore the new e. If the new e is greater than the current e, then we replace the old e with the new e.
			# 2.1) TODO: But what if one possible leaving var has e_1 + e_1 and another only have e_1? Shouldn't I keep track of the magnitude of the largest epsilon in each rhs_epsilon row??? (Come back to this later once most cases are working)
		
			# Find leaving_pos for lex method
			min_ratio = np.inf # numpy infinity
			largest_epsilon = -1 # used to find smallest lexicographic epsilon value (which marked as the largest epsilon label)
			leaving_pos = -1
			# a flag variable to check if LP is unbounded
			is_unbounded = True # TODO: Should the unbounded check be used while using the Lex method?? <--- POSSIBLE CAUSE OF ERROR
			for pos in range(np.shape(lhs)[0]):
				# if this if true, then there exists a constraint s.t. the entering variable cannot go to infinity
				if lhs[pos, entering_pos] > 0:
					is_unbounded = False
					# if there are any 0 values in rhs, then we find the one with the largest associated e value and chose that for the leaving pos.
					if rhs[pos] == 0 and rhs_epsilon[pos] > largest_epsilon:
						min_ratio = 0 # now no rhs values greater than 0 can be the leaving_pos
						largest_epsilon = rhs_epsilon[pos]
						leaving_pos = pos
						continue
					# if no 0 values in rhs, then do the following:
					if rhs[pos] / lhs[pos,entering_pos] < min_ratio and rhs[pos] > 0:
						min_ratio = rhs[pos] / lhs[pos, entering_pos]
						leaving_pos = pos

		# if not in cycle_mode, find the leaving variable like this:
		else:
			min_ratio = np.inf # numpy infinity
			leaving_pos = -1
			# a flag variable to check if LP is unbounded
			is_unbounded = True
			for pos in range(np.shape(lhs)[0]):
				# if this if true, then there exists a constraint s.t. the entering variable cannot go to infinity
				if lhs[pos, entering_pos] > 0:
					is_unbounded = False
					if rhs[pos] / lhs[pos,entering_pos] < min_ratio and rhs[pos] >= 0: # ******************************* Previously this was just >0
						min_ratio = rhs[pos] / lhs[pos, entering_pos]
						leaving_pos = pos
		
		# if the LP is unbounded, print and exit entire program
		if is_unbounded:
			print('unbounded')
			exit(0)
		'''
		print("LP BEFORE PIVOTING:")
		print(obj)
		print(lhs)
		print(rhs)
		print('entering_pos: ', entering_pos)
		print('leaving_pos: ', leaving_pos)
		print('current obj_val: ', obj_val)
		print('***************')
		'''
		# pivot rhs_epsilon values
		if cycle_mode:
			'''
			print('*************')
			print('leaving_pos: ', leaving_pos)
			print('entering_pos: ', entering_pos)
			print('obj_val: ', obj_val)
			print('obj:\n', obj)
			print('lhs:\n', lhs)
			print('rhs:\n', rhs)
			print('rhs_epsilon pre-pivot:\n', rhs_epsilon)
			'''
			# TODO: appropriately pivot rhs_epsilon, then the normal pivot can be called as usual.
			for pos in range(len(rhs_epsilon)):
				if pos == leaving_pos:
					# epsilon value remains the same for the leaving_pos row
					continue
				if lhs[leaving_pos, entering_pos] > 0 and rhs_epsilon[leaving_pos] < rhs_epsilon[pos]:
					rhs_epsilon[pos] = rhs_epsilon[leaving_pos]
			#print('rhs_epsilon post-pivot:\n', rhs_epsilon)
			
		# pivot with found entering and leaving vars:
		obj, obj_vars, obj_val, lhs, lhs_vars, rhs = pivot_dictionary(obj, obj_vars, obj_val, lhs, lhs_vars, rhs, entering_pos, leaving_pos) 
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
	
	'''
	if len(sys.argv) != 2:
		print("Error: this program requires exactly 1 argument -> an input file containing the linear program to be solved.")
		print("-----> Example of correct program call: $ python LP_solver.py input_file.txt")
		return 0
	'''	

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
				optimal_coeffs[int(marker-1)] = rhs_val
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

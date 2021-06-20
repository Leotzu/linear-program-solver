import numpy as np
import sys
from scipy.optimize import linprog


#def get_input():
	#nada

def main():
	if len(sys.argv) != 2:
		print("Error: this program requires exactly 1 argument -> an input file containing the linear program to be solved.")
		print("-----> Example of correct program call: $ python LP_solver.py input_file.txt")
		return 0
	
	# get input file name
	filename = sys.argv[-1]
	# open input file and read lines
	f = open(filename)
	lines = f.readlines()
	
	for iter, line in enumerate(lines):
		#print(line)
		if iter is 0:
			obj_fun = np.array([float(x.strip()) for x in line.split()])
			obj_fun *= -1	
			num_vars = len(obj_fun)
			num_constraints = len(lines)
			# declare lhs and rhs matrices with proper dimensions  
			lhs = np.zeros(shape=(num_constraints-1, num_vars))
			rhs = np.zeros(shape=(num_constraints-1, 1))
			continue

		#__import__('code').interact(local={k: v for ns in (globals(), locals()) for k, v in ns.items()})
		#lhs = np.append(lhs, [np.array([float(x.strip()) for x in line.split()])], axis=0)
		
		# for every line in the txt file, fill in that row for both lhs and rhs (Apparently this is the efficient numpy way)
		lhs[iter-1] = np.array([float(x.strip()) for x in line.split()[:num_vars]])		
		rhs[iter-1] = np.array([float(x.strip()) for x in line.split()[num_vars:]])		
		
		#__import__('code').interact(local={k: v for ns in (globals(), locals()) for k, v in ns.items()})
	'''
	# print out the parsed data for bug inspection	
	print("==================================")
	print("objective function:\n", obj_fun)
	print("==================================")
	print("Left hand side:\n", lhs)
	print("==================================")
	print("Constants:\n", rhs)
	print("==================================")
	'''
	#__import__('code').interact(local={k: v for ns in (globals(), locals()) for k, v in ns.items()})
	
	# use linprog() to calculate solution
	result = linprog(c=obj_fun, A_ub=lhs, b_ub=rhs, method="simplex")
	
	# print out results
	if result.success:
		# Note: result.fun will be negative, since we solved the minimization problem, not max.
		optimal_soln = -result.fun
		#__import__('code').interact(local={k: v for ns in (globals(), locals()) for k, v in ns.items()})
		print("optimal")	
		# this strips unecessary zeros and prints up to 7 decimals if they exist
		print(('%0.7f' % optimal_soln).rstrip('0').rstrip('.'))
		for i in result.x:
			print(('%0.7f' % i).rstrip('0').rstrip('.'), end=' ')
		print()
	elif result.status == 2:
		print("infeasible")
	elif result.status == 3:
		print("unbounded")
	else:
		print("UNKNOWN ERROR")
	
	#print(result)
	
	f.close()
	return 0

if __name__ == "__main__":
    main()

import numpy as np
from scipy.optimize import linprog

#def get_input():
	#nada

def main():
	print("I'm a file!")

	f = open('input_01.txt')
	lines = f.readlines()
	
	for iter, line in enumerate(lines):
		#print(line)
		if iter is 0:
			obj_fun = np.array([float(x.strip()) for x in line.split()])
			
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
	
		
	print("==================================")
	print("objective function:\n", obj_fun)
	print("==================================")
	print("Left hand side:\n", lhs)
	print("==================================")
	print("Constants:\n", rhs)
	print("==================================")
	
	f.close()
	return 0

if __name__ == "__main__":
    main()

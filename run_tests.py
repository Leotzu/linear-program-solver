import os
import sys
import subprocess
import filecmp

def main():
	tests_count = 0
	pass_count = 0
	TEST_DIR = "test_LPs/input/"
	#begin = int(sys.argv[1])
	for TEST_FILE in os.listdir(TEST_DIR):
		tests_count += 1
		#os.system('python3 LP_solver.py "' + TEST_DIR + '"445k21_Lecture01_bakery.txt')
		#os.system('python3 LP_solver.py ' + TEST_DIR + '445k21_Lecture01_bakery.txt')
		#print('python3 LP_solver.py ' + TEST_DIR + '445k21_Lecture01_bakery.txt')
		os.system('python3 LP_solver.py < ' + TEST_DIR + TEST_FILE + " > terminal_output.txt")
		
		#command = 'python3 LP_solver.py ' + TEST_DIR + TEST_FILE
		#result = subprocess.run('python3 LP_solver.py ' + TEST_DIR + TEST_FILE, capture_output=True)
		#result = subprocess.run(command)
		#if (os.system('python3 LP_solver.py ' + TEST_DIR + TEST_FILE) == 
		#print(TEST_FILE + ': ')
		expected_out = 'test_LPs/output/' + TEST_FILE
		f_expected = open(expected_out, 'r')
		f_actual = open('terminal_output.txt', 'r')
		passed_flag = 1 # checks if "pass" needs to be printed after for loop
		line_count = 0 # used to detect "pass*" where dif optimal soln is found
		for line_exp, line_act in zip(f_expected, f_actual):
			if line_exp != line_act:
				# check which line didn't match. If both are optimal but with dif soln, then "pass*"
				if line_count >= 2:
					pass_count += 1
					print("pass* --> " + TEST_FILE)
				else:
					print("FAIL  --> " + TEST_FILE)
				passed_flag = 0
				break
			line_count += 1
		if passed_flag:
			pass_count += 1
			print("pass  --> " + TEST_FILE)
	print('==================\n' + str(pass_count) + '/' + str(tests_count) + ' tests passed\n==================')


if __name__ == "__main__":
	main()

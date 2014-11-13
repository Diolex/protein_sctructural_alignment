import numpy as np
 
def horn(a,b,w):

	w = w.reshape(-1)
	r = np.array([[1,0,0],[0,1,0],[0,0,1]])
	
	
	sum_of_wts = np.sum(w)
	weights=1.0*w/sum_of_wts
	sqrtwts = np.sqrt(weights)
	
	lc = np.dot(weights,a).reshape(-1,1)
	rc = np.dot(weights,b).reshape(-1,1)
	
	left = mat_minus_vector(a,lc)
	left = mat_times_vector(left,sqrtwts)
	
	
	right = mat_minus_vector(b,rc);
	right = mat_times_vector(right,sqrtwts)
	
	
	M = np.dot(left.T,  right)
	Sxx = M[0][0]
	Syx = M[1][0]
	Szx = M[2][0]
	Sxy = M[0][1]
	Syy = M[1][1]
	Szy = M[2][1]
	Sxz = M[0][2]
	Syz = M[1][2]
	Szz = M[2][2]
	
	N = np.array([[(Sxx+Syy+Szz),(Syz-Szy),(Szx-Sxz),(Sxy-Syx)],[(Syz-Szy),(Sxx-Syy-Szz),(Sxy+Syx),(Szx+Sxz)],[(Szx-Sxz),(Sxy+Syx),(-Sxx+Syy-Szz),(Syz+Szy)],[(Sxy-Syx),(Szx+Sxz),(Syz+Szy),(-Sxx-Syy+Szz)]])
	
	l,V = np.linalg.eig(N)
	D = np.diag(l)
	l_real = l.real
	emax = l_real.argmax()
	
	q = V[:,emax]
	q = q.real
	
	ii = abs(q).argmax()
	sgn = np.sign(q[ii])
	
	
	quat = q.T
	nrm = np.linalg.norm(quat)
	
	quat = quat/nrm;
	
	q0 = quat[0]
	qx = quat[1]
	qy = quat[2]
	qz = quat[3]
	v = quat[1:4]
	
	Z = np.array([[q0, -qz,qy], [qz,q0,-qx],[-qy,qx,q0]])
	
	r = np.outer(v,v.T) + np.dot(Z,Z)
	sss = 1
	temp = np.dot(r , lc)
	t = rc - temp
	return r,t

def mat_minus_vector(M, v):
	temp = np.ones(M.shape[0])
	temp = np.outer(temp, v)
	result = M - temp
	return result

def mat_plus_vector(M, v):
	temp = np.ones(M.shape[0])
	temp = np.outer(temp, v)
	result = M + temp
	return result	
	
def mat_times_vector(M, v):
	temp = np.ones(M.shape[1])
	temp = np.outer(v, temp)
	result = np.ones(M.shape)
	for i in range(M.shape[0]):
		for j in range(M.shape[1]):
			result[i,j] = M[i,j] * temp[i,j]
	return result

def main():

	print("Hello World!\n");	
	a = np.array ([[0.7060,0.0318,0.2769],[0.0462,0.0971,0.8235],[0.6948,0.3171,0.9502], [0.0344,0.4387,0.3816], [0.7655,0.7952, 0.1869]])
	b = np.array ([[0.4898,0.4456,0.6463],[0.7094,0.7547,0.2760],[0.6797,0.6551,0.1626],[0.1190,0.4984, 0.9597], [0.3404,0.5853, 0.2238]])
	w = np.array([0.751267059305653,0.255095115459269,0.505957051665142,0.699076722656686,0.890903252535799])
	r,t = horn(a,b,w)
	print (mat_plus_vector(np.dot(a,r.T), t) - b)
	#print (r)
	#print (t)
	
if __name__ == "__main__":
	main()

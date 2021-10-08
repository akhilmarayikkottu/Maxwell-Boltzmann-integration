import numpy as np
import math as m


xrange = range
global Temp,mass, kB, titer, dcx, dcy, dcz, f_opt, error,tempDim,T_fin,error_fin
tempDim=120
data=np.genfromtxt('temp.txt',skip_header=1)
Temp=data[:,2]
T_fin=np.zeros(tempDim)
error_fin=np.zeros(tempDim)
# ~ print Temp

def Mx(T,c):
    m = 28*1.66E-27
    k =1.38064852E-23
    A = (m/(2*3.141*k*T))
    B = 2.718281828459045**(-m*c*c/(2*k*T))
    return(A**0.5*B)
    

mass = 28*1.66E-27
kB = 1.38064852E-23
titer = 30
f_opt = 0.0
f_opt2 = 0.0

#data1=np.genfromtxt('Sampled_data_x.txt',skip_header=0)
#Cx = data1[:]
#data2=np.genfromtxt('Sampled_data_y.txt',skip_header=0)
#Cy = data2[:]
#data3=np.genfromtxt('Sampled_data_z.txt',skip_header=0)
#Cz = data3[:]
Cx = np.linspace(-10000,10000,50000)
Cy = np.linspace(-10000,10000,50000)
Cz = np.linspace(-10000,10000,50000)

dcx = (np.max(Cx)-np.min(Cx))/len(Cx)
dcy = (np.max(Cy)-np.min(Cy))/len(Cy)
dcz = (np.max(Cz)-np.min(Cz))/len(Cz)
T = np.zeros((tempDim, (titer+1)))
error = np.zeros((tempDim, (titer+1)))
for i in range(tempDim):
	T[i][0]=Temp[i]
error[0][0] = 0.0
h = 1e-5
num_sample = len(Cx)
 
def Sum(cx,cy,cz,T):
	Term1X = 0.5*mass*cx*cx
	Term1Y = 0.5*mass*cy*cy
	Term1Z = 0.5*mass*cz*cz
	Term2X = Mx(T,cx)
	Term2Y = Mx(T,cy)
	Term2Z = Mx(T,cz)

	Sum1 = np.dot(Term1X,Term2X)*dcx
	Sum2 = np.dot(Term1Y,Term2Y)*dcy
	Sum3 = np.dot(Term1Z,Term2Z)*dcz
	Sum4 = 2*(Sum1+Sum2+Sum3)/(3*kB)
	return(Sum4)
	
for k in xrange(tempDim):
	for i in xrange(titer):
		print (k,i,T[k][i])
		for j in xrange(num_sample):
			f_opt += Sum(Cx[j],Cy[j],Cz[j],T[k][i])
			f_opt2 += Sum(Cx[j],Cy[j],Cz[j],T[k][i]+h)
		der = (f_opt2-f_opt)/h
		f_opt = T[k][i]-f_opt
		T[k][i+1] = T[k][i]-f_opt/der
		error[k][i+1]=T[k][i+1]-T[k][i]
		print (k,i,T[k][i+1])
		f_opt2 = 0.0

# bisection method
# ~ def func(x):
	# ~ f = 0.0
	# ~ for j in xrange(50000):
		# ~ f += Sum(Cx[j],Cy[j],Cz[j],x)
	# ~ f = x-f
	# ~ return(f)
	
# ~ def bisection(a,b):
  
    # ~ if (func(a) * func(b) >= 0):
        # ~ print("You have not assumed right a and b\n")
        # ~ return
   
    # ~ c = a
    # ~ while ((b-a) >= 0.01):
  
        # ~ # Find middle point
        # ~ c = (a+b)/2
   
        # ~ # Check if middle point is root
        # ~ if (func(c) == 0.0):
            # ~ break
   
        # ~ # Decide the side to repeat the steps
        # ~ if (func(c)*func(a) < 0):
            # ~ b = c
        # ~ else:
            # ~ a = c
              
    # ~ print("The value of root is : ","%.4f"%c)	

# ~ a = 1
# ~ b = 1000
# ~ bisection(a,b)

for i in range(tempDim):
	T_fin[i]=T[i][titer]
	error_fin[i]=error[i][titer]

res = open("Iter_T_err.txt","w+")
res.write("bin number"+"Iteration"+"\t"+"Temp"+"\t"+"Fraction Error"+"\n")
for k in range(tempDim):
	for i in range(titer):
		res.write(str(k)+"\t"+str(i)+"\t"+str(T[k][i])+"\t"+str(error[k][i])+"\n")
plt.figure()
plt.subplot(211)
plt.subplot(211).set_title("Temperature")	
plt.plot(T_fin,'bo')
plt.subplot(212)
plt.subplot(212).set_title("Error")	

plt.plot(error_fin,'r--')

plt.show()
plt.figure()
plt.plot().set_title("Error")	
plt.plot(error_fin,'r--')

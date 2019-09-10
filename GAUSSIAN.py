import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import time
import sys
sys.setrecursionlimit(50000)
#THE SOLUTION OF THE MATRIX BY ELIMINATION
#1)Definition of 3 matrices:matrix of gamma coefs, solution matrix and initial state matrix
#2)Solution of spacial temperature distribution for specific temporal point
#3)Insertion of the solutions into the 'solution' matrix and proceeding to the next temporal cordinate

class Matrix:
    def __init__(self,dim_s,dim_t,g,T1,T2):
        self.g = g
        self.dim = dim_s #defining the dimension of tridiagonal matrix
        self.dim2 = dim_t #defining time marching dimensions
        self.init_state = np.zeros(dim_s) #defining the known column for 0 state, initial state matrix
        self.init_state[0] = T1
        self.init_state[dim_s-1] = T2
        self.solution = np.zeros((dim_t,dim_s))#creating an empty matrix for solutions
        self.solution[0] = self.init_state
        self.matrix = np.zeros((dim_s,dim_s+1)) #creating a gamma tridiagonal matrix
        for i in range(dim_s):
            if i == 0 or i == dim_s-1:
                self.matrix[i,i] = 1
            else:
                for j in range(i-1,i+2):
                    if j == i-1 or j == i+1:
                        self.matrix[i,j] = -g
                    elif j==i:
                        self.matrix[i,j] = 2*g+1
        self.matrix[:,dim_s] = self.init_state

    def definition(self):
        for i in range(self.dim):
            print('Please define the ',i+1,' row\n')
            for j in range(self.dim+1):
                self.matrix[i,j] = float(input())
        self.echo()

    def implicit_Thomas(self):
        n = self.dim
        b = np.zeros(self.dim)#main diagonal
        f = np.zeros(self.dim)#lower
        c = np.zeros(self.dim)#upper
        #the process of filling up the matrices
        for i in range(n):
            if i == 0:
                b[i] = self.matrix[i,i]
                c[i] = self.matrix[i,i+1]
                f[i] = 0
            elif i == (self.dim-1):
                b[i] = self.matrix[i,i]
                f[i] = self.matrix[i,i-1]
                c[i] = 0
            else:
                b[i] = self.matrix[i,i]
                c[i] = self.matrix[i,i+1]
                f[i] = self.matrix[i,i-1]
        i = 0
        for k in range(1,self.dim2):
            d = self.solution[k-1]#solutions
            A = np.zeros(self.dim-1)#Matrix holding new up diagonal terms
            B = np.zeros(self.dim)
            for i in range(n):
                if i==0:
                    A[i] = c[i]/b[i]
                    B[i] = d[i]/b[i]
                elif i==n-1:
                    B[i] = (d[i]-f[i]*B[i-1])/(b[i]-f[i]*A[i-1])
                else:
                    A[i] = c[i]/(b[i]-f[i]*A[i-1])
                    B[i] = (d[i]-f[i]*B[i-1])/(b[i]-f[i]*A[i-1])
            self.solution[k,n-1] = B[n-1]
            for i in range(1,n):
                self.solution[k,n-1-i] = B[n-1-i]-A[n-i-1]*self.solution[k,n-i]



    def implicit(self):
        for k in range(self.dim2):#for time marching
            for i in range(self.dim-1): #the first loop is for marching along the main row
                for j in range(i+1,self.dim): #the second loop is for perfoming ellimination
                    self.matrix[j] = -(self.matrix[j,i]/self.matrix[i,i])*self.matrix[i] + self.matrix[j]
            self.matrix[self.dim-1] = self.matrix[self.dim-1]/self.matrix[self.dim-1,self.dim-1]
            #backwards substitution
            for i in range(1,self.dim+1):
                self.matrix[self.dim-i] = self.matrix[self.dim-i]/self.matrix[self.dim-i,self.dim-i]#row initialization
                for j in range(1,i):
                    self.matrix[self.dim-i,self.dim] = self.matrix[self.dim-i,self.dim] - self.matrix[self.dim-j,self.dim]*self.matrix[self.dim-i,self.dim-j]
                    self.matrix[self.dim-i,self.dim-j] = 0
#THE 'SOLUTION' MATRIX PART
            self.init_state = self.matrix[:,self.dim]#redefining the initial state matrix,current initial state is always
            #the 'self.matrix' solution
            self.solution[k] = self.init_state#filling up the 'solution' matrix from the solved 'self.matrix'
            self.matrix = np.zeros((self.dim,self.dim+1)) #recreating a gamma tridiagonal matrix
    #is there a possibility to avoid the previous state, in order to keep the last column unaltered and not reinitialize it twice
            for i in range(0,self.dim):
                if  i== 0 or  i== self.dim-1:
                    self.matrix[i,i] = 1
                else:
                    for j in range(i-1,i+2):
                        if j == i-1 or j == i+1:
                            self.matrix[i,j] = -self.g
                        elif j==i:
                            self.matrix[i,j] = 2*self.g+1
            self.matrix[:,dim_s] = self.init_state

#Solving PDE using implicit scheme
    def explicit(self):
        for i in range(1,self.dim2):
            for j in range(0,self.dim):
                if j == 0:
                    self.solution[i,j]=self.solution[0,0]
                elif j == (self.dim-1):
                    self.solution[i,j]=self.solution[0,j]
                else:
                    self.solution[i,j] = self.g*self.solution[i-1,j-1]+(1-2*self.g)*self.solution[i-1,j]+self.g*self.solution[i-1,j+1]


    def echo(self,object):
        print(object)

#choice of solution method

def sol_method():
    print('Please choose the solution scheme:\'implicit\' or \'explicit\'')
    m = input()
    if m == 'implicit':
        c.implicit_Thomas()
    elif m == 'explicit':
        c.explicit()

#Defining a function for the display of the solution
def display(t,dt,L,dx,T1,T2,a):
    print('I have allowed for two options for you to visualize the solution obtained:\n\'Statical with intervals\' and \'Dynamic with animation\'')
    print('Please if you prefer option one type \'statical\' if option two type \'dynamic\'')
    dis = input()
    if dis=='statical':
        #plotting initial condition and steady state solution
        #!!!BEGINNING OF COMMENTING OUT!!!
        #defining a figure
        fig, axs = plt.subplots(nrows=3, ncols=2,
                                sharex = True,sharey = True,
                                figsize=(6, 6))
        x = np.linspace(0,L,dim_s)
        #creating suplots
        axs[0,0].plot(x,c.solution[0],color = 'blue',marker = 'o')
        axs[0,0].set_title('T = 0 sec')

        axs[0,1].plot(x,c.solution[int(0.5/dt)],color = 'green',marker = 'o')
        axs[0,1].set_title('T = 0.5 sec')

        axs[1,0].plot(x,c.solution[int(1/dt)],color = 'b',marker = 'o')
        axs[1,0].set_title('T = 1 sec')

        axs[1,1].plot(x,c.solution[int(5/dt)],color = 'g',marker = 'o')
        axs[1,1].set_title('T = 5 sec')

        axs[2,0].plot(x,c.solution[int(7.5/dt)],color = 'b',marker = 'o')
        axs[2,0].set_title('T = 7.5 sec')

        axs[2,1].plot(x,c.solution[dim_t-1],color = 'g',marker = 'o')
        axs[2,1].set_title('T = 10 sec')
        #setting labels
        fig.text(0.5, 0.04, 'Distance', ha='center')
        fig.text(0.04, 0.5, 'Temperature', va='center', rotation='vertical')

        fig.tight_layout()
        plt.subplots_adjust(left=0.15,bottom=0.12)#really nice command which
                                                  #defines the plot proporsion on the figure

        plt.show()
        #THE END
    elif dis=='dynamic':
        T = 0
        if T2>T1:
            T=T2
        elif T2<T1:
            T=T1
        #A small scripto to animate the system's temperature evolution with time
        fig = plt.figure()
        ax = plt.axes(xlim=(-0.02*L,1.05*L), ylim=(-0.05*T,T+0.05*T))
        line, = ax.plot([], [], lw=2,color = 'green')

        def init():
            line.set_data([],[])
            return line,

        def animate(i):
            x = np.linspace(0,L,dim_s)
            y = c.solution[i]
            line.set_data(x,y)
            return line,
        #IMPORTANT!!!
        #If animation takes to long please adjust the frames and interval values
        #If it doesn't fully converge to the linear dependance please
        #increase the t variable
        anim = animation.FuncAnimation(fig, animate,init_func=init,
                                       frames=int(0.03*t/dt),
                                       interval=20,blit=True)
        plt.title('Parabolic PDE heat diffusion solution')
        plt.xlabel('Distance')
        plt.ylabel('Temperature')
        fig.tight_layout()
        plt.show()

#Implicit method for the Parbolic PDE
#t-total time,L-length,dt-time step,dx-spacial step,T1&T2 - coundary conditions
#FUNCTION TO INITIALIZE THE PROGRAM. ALL INITIAL CONDITIONS
def initialization():
    print('Please define all the variables and boundary conditions neccessary for the PDE solution\n')
    print('Time t = ')
    t = float(input())
    print('Time step dt = ')
    dt = float(input())
    print('Length L = ')
    L = float(input())
    print('Spacial step dx = ')
    dx = float(input())
    print('Boundary condition at x(0) T1 = ')
    T1 = float(input())
    print('Boundary condition at x(L) T2 = ')
    T2 = float(input())
    print('Alpha coefficient = ')
    a = float(input())
    gamma = double_check(t,dt,L,dx,T1,T2,a)
    return t,dt,L,dx,T1,T2,a,gamma

def double_check(t,dt,L,dx,T1,T2,a):
    #checking for stability criterion
    gamma = a*dt/(dx*dx)
    print('Please check if values are correct')
    print('t = ',t)
    print('dt = ',dt)
    print('L = ',L)
    print('dx = ',dx)
    print('T1 = ',T1)
    print('T2 = ',T2,'\n')
    print('alpha = ',a,'\n')
    print('For stability gamma <= :\n','gamma = ',gamma)
    print('If everything is correct type \'yes\' if not \'no\':')
    #answer variable
    check = input()
    if check =="yes":
        print('++')
        return gamma
    elif check == 'no':
        initialization()

#calling for initialization of the further variables
t,dt,L,dx,T1,T2,a,gamma = initialization()
dim_s = int(L/dx)
dim_t = int(t/dt)

c = Matrix(dim_s,dim_t,gamma,T1,T2)
start_time = time.time()
sol_method()
#PLEASE LEAVE UNCOMMENTED THE SCHEME YOU WOULD LIKE TO USE
#NOTE!!!c.implicit uses full GAUSSIAN elimination, which increases the computational time exponantialy
#c.implicit_Thomas()
#c.explicit()
#c.implicit()
#c.echo(c.solution)
print("COMPILING TIME\n--- %s seconds ---" % (time.time() - start_time))

display(t,dt,L,dx,T1,T2,a)

#1)perform error analysis
#2)brush up the graphs (DONE)
#3)write the report(DONE)
#4)solve the tridiagonal matrix using Thomas algorithm (DONE)

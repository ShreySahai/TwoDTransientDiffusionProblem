import 	numpy 					as 		np 
import 	matplotlib.animation 	as 		ani
import 	matplotlib.pyplot 		as 		plt
from 	matplotlib.pyplot 		import 	imshow, show, colorbar
class TwoDHeatTransfer:														#simulates transient 2d heat conduction using ADI method
	def __init__(self,parameters):
		"""parameters ={
			"Length"		:
			"T0"			:
			"Tin"			:
			"alpha"			:
			"Nelements"		:
			"timestep"		:
			"spacestep"		:
			"Niteration"	:
		}"""
		self.L 		=	parameters["Length"]
		self.T0		=	parameters["T0"]
		self.Tin	=	parameters["Tin"]
		self.alpha	=	parameters["alpha"]
		self.N 		=	parameters["Nelements"]
		self.dt 	=	parameters["timestep"]
		self.dx 	=	round(parameters["Length"]/parameters["Nelements"],1)
		self.nit 	= 	parameters["Niteration"]
		self.r 		=	self.dt/(self.dx * self.dx)							#varaiables handy while the TDM method 
		self.b 		=	2*(1/self.r+1)										#varaiables handy while the TDM method 
		self.f 		=	2*(1/self.r-1)										#varaiables handy while the TDM method 
	def formMatrix(self,d,N):												#forms matrix which arises as result of Alternating Directional Implicit method 
		A=np.zeros((N,N))													# [ b -2					]
		b=self.b															# [-1  b  -1				]
		for k in range(0,N):												# [   -1   b  -1			]
			A[k][k]=b														# [		  -1   b  -1		]
			if(k==0):														# [		      -1   b  -1	]
				A[k][k+1]=-2												# [   			  -1   b  -1]
				continue													# [   				   -1  b]
			A[k][k-1]=-1													
			if(k!=N-1):														
				A[k][k+1]=-1												
		return self.TDMA(A,d)
	def updatevalue(self,A,B,N):											#just copies values of one array onto other
		for i in range(0,N):
			for j in range(0,N):
				A[i][j]=B[i][j]
	def dirichelet(self,thetha,N):											#assigns dirichelet boundary condition of thetha = 1 at boundaries
		for x in range(0,N+1):
			thetha[N][x]=1
			thetha[x][N]=1													
	def TDMA(self,A,b):														#solve Ax = b using TDM algorithm
	    p=[]
	    q=[]
	    p.append(-1*A[0][1]/A[0][0])
	    q.append(b[0]/A[0][0])
	    x=np.zeros(len(b))
	    for i in range(1,len(b)):
	        if i!= len(b)-1:
	            p.append(-1*A[i][i+1]/(A[i][i]+A[i][i-1]*p[i-1]))
	        q.append((b[i]-A[i][i-1]*q[i-1])/(A[i][i]+A[i][i-1]*p[i-1]))
	    x[len(b)-1]=q[len(b)-1]
	    for i in range(len(b)-2,-1,-1):
	        x[i]=p[i]*x[i+1]+q[i]
	    return x													
	def Halfstepi(self,thetha):												#Implicit iteration over x-coordinate
		N=self.N
		f=self.f
		j=0
		d=[]
		thethastar=np.full((N+1,N+1),0,dtype='f')
		for i in range(0,N-1):
			d.append(2*thetha[i][j+1]+f*thetha[i][j])
		d.append(2*thetha[N-1][j+1]+f*thetha[N-1][j]+thetha[N][j])
		u=self.formMatrix(d,N)
		for i in range(0,N):
			thethastar[i][j]=u[i]
		for j in range(1,N):
			d=[]
			for i in range(0,N-1):
				d.append(thetha[i][j-1]+f*thetha[i][j]+thetha[i][j+1])
			d.append(thetha[N-1][j-1]+f*thetha[N-1][j]+thetha[N][j+1]+thetha[N][j])
			u=self.formMatrix(d,N)
			for i in range(0,N):
				thethastar[i][j]=u[i]
		return thethastar
	def Halfstepj(self,thetha):												#Implicit iteration over y-coordinate
		N=self.N
		f=self.f
		i=0
		d=[]
		thethastar=np.full((N+1,N+1),0,dtype='f')
		for j in range(0,N-1):
			d.append(2*thetha[i+1][j]+f*thetha[i][j])
		d.append(2*thetha[i+1][N-1]+f*thetha[i][N-1]+thetha[i][N])
		u=self.formMatrix(d,N)
		for j in range(0,N):
			thethastar[i][j]=u[j]
		for i in range(1,N):
			d=[]
			for j in range(0,N-1):
				d.append(thetha[i-1][j]+f*thetha[i][j]+thetha[i+1][j])
			d.append(thetha[i-1][N-1]+f*thetha[i][N-1]+thetha[i+1][N]+thetha[i][N])
			u=self.formMatrix(d,N)
			for j in range(0,N):
				thethastar[i][j]=u[j]
		return thethastar
	def extendsquare(self,A,N):												#Extended the matrix in all direction(Symmetry)
		box=np.full((2*N-1,2*N-1),0,dtype='f')
		for i in range(0,N):
			for j in range(0,N):
				box[N-1-i][N-1-j]=A[i][j]
				box[N-1-i][N-1+j]=A[i][j]
				box[N-1+i][N-1-j]=A[i][j]
				box[N-1+i][N-1+j]=A[i][j]
		return box
	def solve(self):														#Simulates the diffusion using ADI method
		N=self.N
		fig = plt.figure()													
		ax = fig.add_subplot(111)
		thetha=np.full((N+1,N+1),0,dtype='f')								#Dimensionless temperature
		ims=[]																#stores thetha at different timesteps later used to visualise
		im=imshow(self.extendsquare(thetha,N+1),animated=True)
		ims.append([im])
		self.dirichelet(thetha,N)
		for t in range(self.nit):											#time marching
			thethastar=self.Halfstepi(thetha)								#updates thetha after  half time step
			self.updatevalue(thetha,thethastar,N)
			thethastar=self.Halfstepj(thetha)
			self.updatevalue(thetha,thethastar,N)							#updates thetha after second half time step					
			box=self.extendsquare(thetha,N+1)
			temp=int(thetha[0][0]*(self.T0-self.Tin)+self.Tin)				#next three lines are for the showing timestep and temperature  
			time = plt.text(0.1, 0.1, "Timestep: "+str(t), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
			temp = ax.text(0.98,0.95, "Core Temperature: "+str(temp)+"K", horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
			im=imshow(box,vmin=0,vmax=1,animated=True)
			if t % 2==0:													#only visualising at every other time step
				ims.append([im,time,temp])
		anim =ani.ArtistAnimation(fig, ims, blit=True,interval=600,repeat_delay=10)
		show()

import numpy as np
import math

"""
Esta libreria es la libreria de Quaterniones generada por y para el Instituto
de Automatica INAUT:
    link
"""
class Quaternion:
    """
    Objeto cuaternion creado ad hoc para este caso
    Se toma:
    q=w+xi+yj+zj=[w,a]
    q=cos(fi/2)+sin(fi/2)*[x,y,z]
    
    """
    def __init__(self,q=np.array([1,0,0,0]),unitary=True):
        self.q=q
        if unitary==True:
            self.is_unit()
        
        self.unitary=unitary
        self.w=q[0]
        self.n=q[1:]
        self.m=np.array([[q[0],-q[1],-q[2],-q[3]],
                         [q[1],q[0],-q[3],q[2]],
                         [q[2],q[3],q[0],-q[1]],
                         [q[3],-q[2],q[1],q[0]]])



    def is_unit(self):
        norm=self.norm()
        if norm!=1.0:
            #print("Warning, is not unitary quaternion. For avoid this warning set unitary=False. Norm is: ",norm)
            if norm==0:
                print("Error, quaternion is null")
                return

            self.q=self.q*(1/norm)
            return 
         


    def conj(self):
        cq=-self.q
        cq[0]=-cq[0]
        return Quaternion(cq,unitary=self.unitary)

    def inner_product(self,p ):
        # Según taller 1 J. Gimenez
        q=self.q
        p=p.q
        ip=0
        for i in range(4):
            ip+=p[i]*q[i]

        return ip

    def norm(self)->float:
        return np.sqrt(self.inner_product(self))

    def cross(self,p):
        """
        p=Quaternion()
        p.cross(q)=pq
        """
        return Quaternion(np.matmul(self.m,p.q),unitary=p.unitary)
    
    def rot(self,p):
        """
        p is a quaternion / -conj(p)=p
        p rot = q * p *q^-1
        """
        if self.unitary==False:
            print("Error: unitary quaternion is required")
            return None
        
        p_=np.array([0,p[0],p[1],p[2]])
        p=Quaternion(p_,unitary=False)
        qp=self.cross(p)
        qpq=qp.cross(self.conj())
        return qpq.q[1:]

    def inverse(self):
        q=self.conj()
        norm=q.norm()
        return Quaternion(q.q* 1/norm**2,unitary=self.unitary)

    def euler_from_quaternion(self, R):
            """
            Convert a quaternion into euler angles (roll, pitch, yaw)
            roll is rotation around x in radians (counterclockwise)
            pitch is rotation around y in radians (counterclockwise)
            yaw is rotation around z in radians (counterclockwise)
            Se toma q=[w,x,y,z]
            """
            #x=R[0]
            #y=R[1]
            #z=R[2]
            #w=R[3]
            
            w=R[0]
            x=R[1]
            y=R[2]
            z=R[3]
            
            t0 = +2.0 * (w * x + y * z)
            t1 = +1.0 - 2.0 * (x * x + y * y)
            roll_x = math.atan2(t0, t1)
         
            t2 = +2.0 * (w * y - z * x)
            t2 = +1.0 if t2 > +1.0 else t2
            t2 = -1.0 if t2 < -1.0 else t2
            pitch_y = math.asin(t2)
         
            t3 = +2.0 * (w * z + x * y)
            t4 = +1.0 - 2.0 * (y * y + z * z)
            yaw_z = math.atan2(t3, t4)
         
            return roll_x, pitch_y, yaw_z # in radians
    
    def euler_to_quaternion(self,r):
        (roll, pitch, yaw) = (r[0], r[1], r[2])
        qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
        qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
        qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
        qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
        return [qw, qx, qy, qz]

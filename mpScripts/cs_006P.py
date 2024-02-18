# -*- coding: utf-8 -*-

# https://github.com/whucyb/computational_physics_N2014301020067/blob/master/Exercise_08/Exercise_08.py


import matplotlib.pyplot as plt
import numpy as np
class driven_nonlinear_pendulum(object):
    def __init__(self,F_D):
        self.omega = [0]
        self.theta = [0.2]
        self.dt = (2 * np.pi / (2.0 / 3.0)) / 600
        self.time = [0]
        self.omega_in_poincare_section = []
        self.theta_in_poincare_section = []
        self.F_D = F_D
    def calculate(self):
        l = 9.8
        g = 9.8
        q = 0.5
        Omega_D = 2.0 / 3.0
        for i in range(600000):
            temp_omega = self.omega[-1] - ((g / l) * np.sin(self.theta[-1]) + q * self.omega[-1] - self.F_D * np.sin(Omega_D * self.time[-1])) * self.dt
            temp_theta = self.theta[-1] + temp_omega * self.dt
            while temp_theta > np.pi:
                temp_theta -= 2 * np.pi
            while temp_theta < -np.pi:
                temp_theta += 2 * np.pi
            self.omega.append(temp_omega)
            self.theta.append(temp_theta)
            self.time.append(self.dt * i)
        for i in range(500,1000):
            n = i * 600
            self.omega_in_poincare_section.append(self.omega[n])
            self.theta_in_poincare_section.append(self.theta[n])
    def show_results(self):
        plt.subplot(121)
        p1,=plt.plot(self.time,self.theta)
        plt.title(r'$\theta$ versus time')
        plt.legend([p1],[r'$F_D$ = %.3f'%self.F_D],loc='best')
        plt.xlabel('time (s)')
        plt.ylabel(r'$\theta$ (radians)')
        plt.xlim(0,200)
        plt.ylim(-4,4)
        plt.subplot(122)
        p2,=plt.plot(self.theta_in_poincare_section,self.omega_in_poincare_section,'.')
        plt.title(u'Poincaré section' + '\n' + r'$\omega$ versus $\theta$')
        plt.legend([p2],[r'$F_D$ = %.3f'%self.F_D],loc='best')
        plt.xlabel(r'$\theta$ (radians)')
        plt.ylabel(r'$\omega$ (radians/s)')
        plt.xlim(-4,4)
        plt.ylim(-3,3)
        plt.show()
class bifurcation_diagram(object):
    def __init__(self):
        self.omega = []
        self.theta = []
        self.dt = (2 * np.pi / (2.0 / 3.0)) / 600
        self.time = []
        self.theta_in_bifurcation_diagram = []
        self.F_D = np.arange(1.35,1.5,0.001)
        self.theta_versus_F_D = []
    def calculate(self):
        l = 9.8
        g = 9.8
        q = 0.5
        Omega_D = 2.0 / 3.0
        for f_d in self.F_D:
            self.omega.append([0])
            self.theta.append([0.2])
            self.time.append([0])
            for i in range(600000):
                temp_omega = self.omega[-1][-1] - ((g / l) * np.sin(self.theta[-1][-1]) + q * self.omega[-1][-1] - f_d * np.sin(Omega_D * self.time[-1][-1])) * self.dt
                temp_theta = self.theta[-1][-1] + temp_omega * self.dt
                while temp_theta > np.pi:
                    temp_theta -= 2 * np.pi
                while temp_theta < -np.pi:
                    temp_theta += 2 * np.pi
                self.omega[-1].append(temp_omega)
                self.theta[-1].append(temp_theta)
                self.time[-1].append(self.dt * i)
            for i in range(500,1000):
                n = i * 600
                self.theta_in_bifurcation_diagram.append(self.theta[-1][n])
                self.theta_versus_F_D.append(f_d)
                print(self.theta_versus_F_D.append(f_d))
    def show_results(self):
        plt.plot(self.theta_versus_F_D,self.theta_in_bifurcation_diagram,'.')
        plt.title('Bifurcation diagram' + '\n' + r'$\theta$ versus $F_D$')
        plt.xlabel(r'$F_D$')
        plt.ylabel(r'$\theta$ (radians)')
        plt.xlim(1.35,1.5)
        plt.ylim(1,3)
        plt.show()



#%%
x = np.array([1,2,3,4,5,6])

z = []
for i in range(6):
    z.append(x[-1][i])
    
print(z)
    
    
    
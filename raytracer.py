# -*- coding: utf-8 -*-
"""
T Molnar 25/01/2018

"""

"""
3-D Optical Ray Tracer
Model investigates the propagation of Ray objects through various optical element objects
"""

"""Importing necessary modules
"""

import numpy as np
import matplotlib.pyplot as plt

from config import *


class Ray:
    def __init__(self, p = np.array([0,0,0]), k = np.array([0,0,1])):
        """Initialising the Ray class with parameters:
            
            p: Starting Position (3D Vector)
            k: Initial Direction (3D Vector)
            """
            
        if len(p) == 3:
            self.__position = [p]
            
        else:
            raise Exception('Position parameter size')
            
        if len(k) == 3:
            self.__direction = [k/np.linalg.norm(k)]
        
        else:
            raise Exception('Direction parameter size')
            
    def p(self):
        "Returns the most recent position of the Ray"
        return self.__position[-1] 
    
    
    def k(self):
        "Returns the most recent direction of the Ray"
        return self.__direction[-1]
    
    def append(self,p,k):
        "Adding new position and direction to the Ray"
        self.__position.append(p)
        self.__direction.append(k)
    
    def vertices(self):
        "Returns all positions of the Ray"
        return self.__position
    
    def directions(self):
        "Returns all directions of the Ray"
        return self.__direction
    
    def Plot(self, focal, radius):
        """ Plots: - Propagation of the Ray in the x-z plane (y=0)
                   - Input Plane of the Ray in the x-y plane (z=0)
                   - Output Plane of the Ray in the x-y plane (z=0)
            
            with the following parameters:
            
            focal: Focal Point of optical system
            radius: Radius of outermost rays in beam 
        """
        self.__radius = radius
        self.__focal = focal
        self.__x, self.__y, self.__z = [], [], []
        for i in self.vertices():
            self.__x.append(i[0])
            self.__y.append(i[1])
            self.__z.append(i[2])
        
        """Input Plane plot"""
        plt.figure(1)
        plt.plot(self.__x[0], self.__y[0], 'ro', markersize = 4)
        plt.title('Input Plane at z = 0 mm (Input Plane)')
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        plt.grid(linestyle='--' )
        plt.gca().set_aspect('equal', adjustable='box')
        
        """Output Plane plot"""
        plt.figure(2)
        plt.plot(self.__x[-1], self.__y[-1], 'ro', markersize = 2.5)
        plt.title('Output Plane at z = %i mm (Paraxial Focal Plane)' % (self.__focal))
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        plt.xlim(-0.06, 0.06)
        plt.ylim(-0.06, 0.06)
        plt.grid(linestyle='--' )
        plt.gca().set_aspect('equal', adjustable='box')
        
        """Ray Propagation plot"""
        plt.figure(3)
        plt.plot(self.__z, self.__x, 'r-')
        plt.title('Ray Propagation for y = 0 mm')
        plt.xlabel('z (mm)')
        plt.ylabel('x (mm)')
        plt.xlim(0, self.__focal)
        plt.ylim(-self.__radius, self.__radius)
        plt.grid(linestyle='--' )
        
class OpticalElement:
    
  def propagate_ray(self, Ray):
    "Propagate a Ray through the optical element"
    raise NotImplementedError()
        
class SphericalRefraction(OpticalElement):
    
    def __init__(self, z0, n1, n2, ra, C):
        """Initialising refractive surface with five parameters:
            
            z0: Intercept of surface with z-axis
            n1: Refractive index of material 1
            n2: Refractive index of material 2
            ra: Aperture radius
            C: Curvature of refractive object
            """
            
        self.__z0 = z0
        self.__C = C
        self.__n1, self.__n2 = n1,n2
        self.__ra = ra
        
        if self.__C != 0:
            self.__R = 1.0/self.__C
            self.__centre = np.array([0,0, self.__z0 + self.__R])
        
    def intercept(self, Ray):
        """Returns intercept of a Ray with refractive surface 
        with a parameter from the Ray class
        
        """
        
        if self.__C == 0:
            l = (self.__z0-Ray.p()[2])/Ray.k()[2]
            s = Ray.p() + l*Ray.k()
            
        else:
            r = Ray.p() - self.__centre
            discriminant = (np.dot(r, Ray.k()))**2-((np.linalg.norm(r))**2 - (self.__R)**2) 

            if discriminant > 0:
            
                l1 = (-np.dot(r, Ray.k())) + np.sqrt(discriminant)
                l2 = (-np.dot(r, Ray.k())) - np.sqrt(discriminant)
                l = [l1,l2]
            
                if self.__C > 0:
                    s = Ray.p() + min(l)*Ray.k()
            
                elif self.__C < 0:
                    s = Ray.p() + max(l)*Ray.k()
            
            else:
                return None
            
        if np.abs(s[0]) < self.__ra or np.abs(s[1]) < self.__ra:
            return s
            
        else:
            return None
            
    def Snell(self, k1, norm, n1, n2):
        """Snell's law implemented for the Ray, with the following parameters:
            
            k1: Incident direction of Ray
            norm: Surface normal of refractive surface (note: points in opposite direction to Ray propagation)
            n1: Refractive index of material 1
            n2: Refractive index of material 2
            
            """
            
        theta1 = np.arccos(-np.dot(k1,norm))
        r = n1/n2
        
        if np.sin(theta1) > n2/n1:
            return None 
        
        else:
            theta2 = np.arcsin(r*np.sin(theta1))
            k2 = r*k1 + (r*np.cos(theta1) - np.cos(theta2))*norm
            return k2
        
    def norm(self, Ray):
        """Returns the normal to the refractive surface
        with a parameter from the Ray class
        
        """
        s = self.intercept(Ray)
        
        if self.__C == 0:
            norm = np.array([0,0,-1])
            
        else:
            sign = self.__C/np.abs(self.__C)
            norm = sign*(s - self.__centre)/np.linalg.norm(s - self.__centre)
        
        return norm
        
        
    def propagate_ray(self, Ray):
        """Propagates a Ray through the optical element
        with a parameter from the Ray class
        
        """
        s = self.intercept(Ray)
        norm = self.norm(Ray)
        k2 = self.Snell(Ray.k(), norm, self.__n1, self.__n2)
        
        if s is None: 
            return None
        
        elif k2 is None:            
            return None

        else:
            Ray.append(s,k2)
            
class OutputPlane(SphericalRefraction):
    
    def __init__(self, z1):
        """Initializng Output Plane class as a SphericalRefraction object
        
        """
        self.__z1 = z1
        self.__n1, self.__n2 = 1, 1 #Arbitrary Constants as ray does not propagate through Output Plane
        self.__ra = 10000 # Near Infinite Aperture
        self.__C = 0
        
        
    def intercept(self, Ray):
        """ Returns intercept with the planar output plane
        with a parameter from the Ray class 
        
        """
        l = (self.__z1-Ray.p()[2])/Ray.k()[2]
        s = Ray.p() + l*Ray.k()
        return s
        
    def propagate_ray(self, Ray):
        s = self.intercept(Ray)
        k2 = Ray.k()
        Ray.append(s, k2)
       
        
def  ParaxialFocus(elements):
    """Returns the paraxial focus for a spherical refraction optical element 
       with the following parameter:
        
        elements: A list of optical elements through which the rays are propagated
    
    """
    test = Ray([0.0001, 0.0, 0.0], [0.0, 0.0, 1.0])
    for elem in elements:   
        s = elem.propagate_ray(test)
    s , k = test.p(), test.k()
    focal = s[2] + k[2]*(s[0]/(abs(k[0])))
    return focal

        
def Plot_Sum(rays, elements, radius):
    """Plots the positions of the rays at the Input and Output Plane and the Trajectories of the rays
       Also returns the RMS deviation of the ray positions from the optical axis
       with the following parameters:
        
        rays: A list of rays which are propagated through optical elements
        elements: A list of optical elements through which the rays are propagated
        radius: Radius of outermost rays in beam
        
        """
    focal = ParaxialFocus(elements)
    Output = OutputPlane(focal)
    elements.append(Output)
    Sum = 0.0
    for i in rays:
        for elem in elements:
            elem.propagate_ray(i)
        s = i.p()
        d_squared = s[0]**2 + s[1]**2
        Sum += d_squared
        i.Plot(focal, radius)
    #plt.figure(1).savefig('InputPlane2.pdf') used to save figures
    #plt.figure(2).savefig('OutputPlane2.pdf')
    #plt.figure(3).savefig('Propagation2.pdf')
    RMS = np.sqrt(Sum/len(rays))
    return RMS
    
def Bundle(radius, n, f):
    """Creates a bundle of rays in the form of a list
    with the following parameters:
        
        radius: Radius of outermost rays in beam
        n: number of different radii in bundle
        f: Factor by which the number of rays increase from neighbour radii
        
        """
    bundle = [Ray([0,0,0], [0,0,1])]
    for i in range(n+1):
        for j in range(f*i):
            angle = 2*j*(np.pi/(f*i))
            r = i*radius/n
            ray = Ray([r*np.cos(angle), r*np.sin(angle), 0], [0,0,1])
            bundle.append(ray)
    return bundle

def Run(radius=5, elements, n_radii = 8, factor = 10):
    """Runs the Program for the following parameters:
        
        radius: Radius of outermost rays in beam
        elements: A list of optical elements through which the rays are propagated 
        n_radii: number of different radii in bundle (equals 8, by default)
        factor: Factor by which the number of rays increase from neighbour radii (equals 10, by default)
        
    """
    rays = Bundle(radius, n_radii, factor)
    Plot_Sum(rays, elements, radius)
    
if __name__=="__main__":

    config = configuration()
    
    """Defining elements in a plano-convex lens
    """
    #Orientation 1     
    curv_plano1 = SphericalRefraction(105.0, 1.5168, 1.0, 1000, -0.02)
    flat_plano1 = SphericalRefraction(100.0, 1.0, 1.5168, 1000, 0.0)
    planoconvex1 = [flat_plano1, curv_plano1]
    #Orientation 2
    curv_plano2 = SphericalRefraction(100.0, 1.0, 1.5168, 1000, 0.02)
    flat_plano2 = SphericalRefraction(105.0, 1.5168, 1.0, 1000, 0.0)
    planoconvex2 = [curv_plano2, flat_plano2]

    """Running program
    """
    Run(radius=config['radius'], elements=planoconvex1,
        n_radii=config['n_radii'], factor=config['factor'])

        

    
    
        

from math import sin,cos,pi

import Angpow

class FuncType0(Angpow.ClassFunc1D_get_value):
  def get_value(self,a_x):    
    return a_x * (a_x-1.)*(a_x-1.)

import math
  
class FuncType1(Angpow.ClassFunc1D_get_value):
  m_ell = 0
  m_R = 0
  def __init__(self,a_ell,a_R):
    Angpow.ClassFunc1D_get_value.__init__(self)
    self.m_ell = a_ell
    self.m_R = a_R
  def get_value(self,a_x):
    return math.cos(a_x*self.m_R - self.m_ell*math.pi*0.5 - math.pi*0.25)

############
def pow(x,a):
  return x**a

##############
def truth_test0(l,R1,R2,a,b):

  M_PI = pi
  res = 0
  
  if R1 != R2:
     res = ((-3*(-2 + pow(a,2.)*pow(R1 - R2,2.))*cos(a*(R1 - R2)))/pow(R1 - R2,4.) - 
               cos(a*(R1 - R2))/pow(R1 - R2,2.) + (4*a*cos(a*(R1 - R2)))/pow(R1 - R2,2.) + 
         (3*(-2 + pow(b,2.)*pow(R1 - R2,2.))*cos(b*(R1 - R2)))/pow(R1 - R2,4.) + 
         cos(b*(R1 - R2))/pow(R1 - R2,2.) - (4*b*cos(b*(R1 - R2)))/pow(R1 - R2,2.) - 
         (a*(-6 + pow(a,2.)*pow(R1 - R2,2.))*sin(a*(R1 - R2)))/pow(R1 - R2,3.) + 
         (2*(-2 + pow(a,2.)*pow(R1 - R2,2.))*sin(a*(R1 - R2)))/pow(R1 - R2,3.) - 
         (a*sin(a*(R1 - R2)))/(R1 - R2) + 
         (b*(-6 + pow(b,2.)*pow(R1 - R2,2.))*sin(b*(R1 - R2)))/pow(R1 - R2,3.) - 
         (2*(-2 + pow(b,2.)*pow(R1 - R2,2.))*sin(b*(R1 - R2)))/pow(R1 - R2,3.) + 
         (b*sin(b*(R1 - R2)))/(R1 - R2) + 
         (a*(R1 + R2)*cos(l*M_PI - a*(R1 + R2)) + sin(l*M_PI - a*(R1 + R2)))/pow(R1 + R2,2.) - 
         (2*((-2 + pow(a,2.)*pow(R1 + R2,2.))*cos(l*M_PI - a*(R1 + R2)) + 
              2*a*(R1 + R2)*sin(l*M_PI - a*(R1 + R2))))/pow(R1 + R2,3.) + 
         (a*(R1 + R2)*(-6 + pow(a,2.)*pow(R1 + R2,2.))*cos(l*M_PI - a*(R1 + R2)) + 
            3*(-2 + pow(a,2.)*pow(R1 + R2,2.))*sin(l*M_PI - a*(R1 + R2)))/pow(R1 + R2,4.) - 
         (b*(R1 + R2)*cos(l*M_PI - b*(R1 + R2)) + sin(l*M_PI - b*(R1 + R2)))/pow(R1 + R2,2.) + 
         (2*((-2 + pow(b,2.)*pow(R1 + R2,2.))*cos(l*M_PI - b*(R1 + R2)) + 
              2*b*(R1 + R2)*sin(l*M_PI - b*(R1 + R2))))/pow(R1 + R2,3.) - 
         (b*(R1 + R2)*(-6 + pow(b,2.)*pow(R1 + R2,2.))*cos(l*M_PI - b*(R1 + R2)) + 
          3*(-2 + pow(b,2.)*pow(R1 + R2,2.))*sin(l*M_PI - b*(R1 + R2)))/pow(R1 + R2,4.))/2.;

  else:
      res = (-2*pow(a,2.)*(6 - 8*a + 3*pow(a,2.))*pow(R1,4.) + 
       2*pow(b,2.)*(6 - 8*b + 3*pow(b,2.))*pow(R1,4.) + 
       6*R1*(2 - 4*pow(a,2.)*pow(R1,2.) + 2*pow(a,3.)*pow(R1,2.) + 
          a*(-3 + 2*pow(R1,2.)))*cos(l*M_PI - 2*a*R1) - 
       6*R1*(2 - 4*pow(b,2.)*pow(R1,2.) + 2*pow(b,3.)*pow(R1,2.) + 
          b*(-3 + 2*pow(R1,2.)))*cos(l*M_PI - 2*b*R1) + 
       3*(-3 + (2 - 8*a + 6*pow(a,2.))*pow(R1,2.))*sin(l*M_PI - 2*a*R1) - 
             3*(-3 + (2 - 8*b + 6*pow(b,2.))*pow(R1,2.))*sin(l*M_PI - 2*b*R1))/(48.*pow(R1,4.)); 
    
  return res
#########

def test0():
  ell = 20
  R1 = 2000.
  R2 = 2200.

  chebyshev_order_1 = 8
  chebyshev_order_2 =  chebyshev_order_1
  n_sub_intervals = 5

  print("ell=",ell," , Nintervales=",n_sub_intervals)
  
  # k-integral bounds
  kMin = 0.
  kMax = 1.0 #Mpc^(-1)
  klp = Angpow.std_vector_double(n_sub_intervals+1)
  dK = kMax-kMin
  for i in range(0,n_sub_intervals+1):
    klp[i] = kMin + dK * i/n_sub_intervals

  f1 = FuncType1(ell,R1)
  f2 = FuncType1(ell,R2)
  f0 = FuncType0()
  
  iOrd0 = 2
  iOrd1 = chebyshev_order_1
  iOrd2 = chebyshev_order_2
  
  farr = Angpow.std_vector_CheFunc()
  che_fun_0 = Angpow.CheFunc(f1, iOrd1)
  farr.push_back(che_fun_0)

  che_fun_1 = Angpow.CheFunc(f2, iOrd2)
  farr.push_back(che_fun_1)

  che_fun_2 = Angpow.CheFunc(f0, iOrd0)
  farr.push_back(che_fun_2)
  
  # Initialisation of the Clenshow-Curtis quadrature
  cheAlgo = Angpow.CheAlgo(farr)

  # Integration
  integral = 0.

  for p in range(1,n_sub_intervals+1):
    # get the bounds
    lowBound = klp[p-1]
    uppBound = klp[p]
    
    if lowBound > uppBound:
      print('KIntegrator::Compute uppBound < lowBound Fatal')
      return
    
    # Loop on each function to compute their  Foward Chebyshev coefficents
    for i in range(0,farr.size()):
      farr[i].ChebyshevTransform(lowBound, uppBound)

    # Compute the sampling of all the functions in the final space dimension
    cheAlgo.InverseChebyshevTransform()

    # Compute the integral thanks to CC quadrature and the function sampling 
    integral += (uppBound - lowBound) * cheAlgo.ComputeIntegralUnscaled()
    
  print("Approx. Integ = ",integral)

  trueInt = truth_test0(ell,R1,R2,kMin,kMax)
  diff = trueInt - integral
  print(f"True Integ ={trueInt}, diff = {diff}")
  
  




if __name__ == '__main__':

  test0()

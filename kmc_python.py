#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 21:12:30 2018

@author: rd
"""

#import sys
import numpy as np
import warnings 
import scipy.constants as cnst
import math
import random
import csv
from itertools import starmap

class organic_kmc(object):
  
  global lamda_i,lamda_s,h_omega,T,n,numatom
  lamda_i=0.245;lamda_s=0.269;numatom=57
  T=300;h_omega=0.2;n=1001
  
  
  
  def __init__(self,po=1,time=0,lifetime=2.0E-13,dist_xyz=None,box_abc=0):
    if (po is None):
      raise ValueError('You must provide initial position of exciton.')
    if type(po)==int:
      self.po=po
      self.time=time
      self.lifetime=lifetime
      self.x = dist_xyz
      self.y = dist_xyz
      self.z = dist_xyz
      self.boxa=box_abc
      self.boxb=box_abc
      self.boxc=box_abc
    else:
      raise ValueError('Please type one integer number between 1 to 1000')
    return


#We are going to read energy file and restore them as a list.

  def read_energy(self):
    get_col = []
    get_energy = []
    #get_col = lambda col: (line.split('\t')[col-1] for line in open('results.dat')) 
    with open('/home/r/Desktop/kmc_py/saamph/saamph_E') as f:
      get_col=[line for line in f]
    for i in range(len(get_col)):
      get_en = float(get_col[i])
      get_energy.append(get_en)
    get_energy.insert(0,0)
    size_mol = np.size(get_energy)
    
    return size_mol,get_energy
  

#We are going to read top of the prj file
  
  def read_parameter(self):
    get_head = []
    with open('/home/r/Desktop/kmc_py/saamph/saamph.prj') as para:
      lines = [line.rstrip() for line in para]
    para.close()
    
    for i in range(8):
     if lines[i].startswith('N')==True:
       get_splt = lines[i]
     else:
    #       _size=np.size(lines[i].split())
       get_splt = [float(elmnt) for elmnt in lines[i].split()]
     get_head.append(get_splt) 
    box_abc=get_head[7][0:3]
        
    return get_head,box_abc,lines
  
#  @property
#  def box_abc(self):
#    self.boxa=self.read_parameter()[1][0]
#    self.boxb=self.read_parameter()[1][1]
#    self.boxc=self.read_parameter()[1][2]
#    return self.boxa,self.boxb,self.boxc
  
  def tot_dimer_coupl(self):
    with open('/home/r/Desktop/kmc_py/saamph/coupl_saamph.dat') as para:
      lines_cou = [line.rstrip() for line in para]
    para.close()
    dim_amph1=[]; dim_amph2=[]; amph_coupl=[]
    for i in range(np.size(lines_cou)):
      aa, bb, cc = lines_cou[i].split()
      dim_amph1.append(int(aa))
      dim_amph2.append(int(bb))
      amph_coupl.append(float(cc))
    
    tot_coupl=np.zeros([n,n])
    for i in range(len(dim_amph1)):
      tot_coupl[dim_amph1[i]][dim_amph2[i]]=amph_coupl[i]
      tot_coupl[dim_amph2[i]][dim_amph1[i]]=tot_coupl[dim_amph1[i]][dim_amph2[i]]
    return tot_coupl,dim_amph1,dim_amph2
  
  def dimer_dict(self):
    dimer1=self.tot_dimer_coupl()[1];dimer2=self.tot_dimer_coupl()[2]
    dimer=[];dimer_lst=[];
    dimer1.extend([0])
    for i in range(len(dimer1)-1):
      if dimer1[i]==dimer1[i+1]:
        dimer.append(dimer2[i])
      else:
        dimer.append(dimer2[i])
        dimer_lst.append(dimer)
        dimer=[]
    dimer_lst.insert(0,0) 
    return dimer_lst
  
  def read_position(self):
    with open('/home/r/Desktop/kmc_py/saamph/SA321_amorphous_center.cen') as para:
      lines_cou = [line.rstrip() for line in para]
    para.close()
    x_mol=[];y_mol=[];z_mol=[]
    for j in range (len(lines_cou)):
       x_,y_,z_=lines_cou[j].split()
       x_mol.append(float(x_));y_mol.append(float(y_));z_mol.append(float(z_))
    x_mol.insert(0,0);y_mol.insert(0,0);z_mol.insert(0,0)
    return x_mol,y_mol,z_mol
  
     
  def cal_dis(self,*argv):
    argv1=int(argv[0]);argv2=int(argv[1])
#    print('---------',argv1,argv2)
    if np.size(argv)<2:
      warnings.warn('The arguments are at least two or more!!!')
    else:
      dis_x=self.read_position()[0][argv2] - self.read_position()[0][argv1]
      dis_y=self.read_position()[1][argv2] - self.read_position()[1][argv1]
      dis_z=self.read_position()[2][argv2] - self.read_position()[2][argv1]
      self.x,self.y,self.z = self.periodic_xyz(dis_x,dis_y,dis_z)
      dis_r=np.sqrt(self.x**2 + self.y**2 + self.z**2)
    return dis_r
  
  def periodic_xyz(self,a,b,c):
    boxa=self.read_parameter()[1][0]
    boxb=self.read_parameter()[1][1]
    boxc=self.read_parameter()[1][2]
    
    if abs(a) > 0.5*boxa:
      tmp=boxa-abs(a)
      if a > 0:
        a=-tmp
      else:
        a=tmp
    
    if abs(b) > 0.5*boxb:
      tmp=boxb-abs(b)
      if b > 0:
        b=-tmp
      else:
        b=tmp
        
    if abs(c) > 0.5*boxc:
      tmp=boxc-abs(c)
      if c > 0:
        c=-tmp
      else:
        c=tmp
  
    return a,b,c
  
  def cal_marcus(self):
    global k_marcus
    f=open('/home/r/Desktop/kmc_py/saamph/saamph_rates.dat', 'w')
    k_marcus=np.zeros((n,n));nozero_rate=[];_dimer1=[];_dimer2=[]
#    k_B=cnst.physical_constants['Boltzmann constant in eV/K'][0]
    k_B=8.617e-05
    h_plnk=cnst.physical_constants['Planck constant in eV s'][0]
    h_plnk=4.137e-15
    E=self.read_energy()[1];t=self.tot_dimer_coupl()[0]
    cnst1=4*cnst.pi**2/(h_plnk*np.sqrt(4*cnst.pi*lamda_s*k_B*T))
    s=lamda_i/h_omega
    for i in range(n):
      for j in range(n):
        
        for kk in range(50):
          k_marcus[i,j]=k_marcus[i,j]+math.exp(-s)*(s**kk/math.factorial(kk))* \
             math.exp(-(E[j]-E[i]+lamda_s+kk*h_omega)**2/(4*lamda_s*k_B*T))
        k_marcus[i,j]=cnst1*(t[i,j])**2*k_marcus[i,j]
        if k_marcus[i,j] !=0:
          nozero_rate.append(k_marcus[i,j])
          _dimer1.append(i)
          _dimer2.append(j)
          row=[i,j,'%00.5E'% k_marcus[i,j]]
          writer = csv.writer(f, delimiter="\t")
          writer.writerow(row)
    return k_marcus,nozero_rate

  
  
  def moving_exiton(self):

    time=0;duf_d=0;times=[];energies=[];krates=[]
    while(time < self.lifetime):
      list_jump=self.dimer_dict()[int(self.po)]
      dtlist=[]
      for i in range(np.size(list_jump)):
        rnd=random.uniform(0,1)
        *argv,=(self.po,list_jump[i])
#        print('-------',argv[0],argv[1])
        k_rate=k_marcus[argv[0]][argv[1]]
        
        dt=-math.log(rnd)/k_rate
        krates.append([argv[0],argv[1],k_rate,dt,-math.log(rnd)])
        dtlist.append(dt)
      mindt=min(dtlist)
#      print(self.po)
#    times will not be same or better to check     
      ind_po=dtlist.index(mindt)
      current_po=list_jump[ind_po]

      dist_=self.cal_dis(self.po,current_po)
      self.po=current_po
      curr_tm_po =(mindt,self.po) 
      time=time+mindt
      duf_d=duf_d+dist_**2
      en_=self.read_energy()[1][self.po]
      energies.append(en_)
      times.append(time)
#      with open('ee', 'a') as f:
#        writer = csv.writer(f, delimiter='\n')
#        writer.writerow(krates)

#    return    
    return en_,time, energies,times

def main():
  organic=organic_kmc()
  
  return organic

if __name__=='__main__':
  dd=main()
  

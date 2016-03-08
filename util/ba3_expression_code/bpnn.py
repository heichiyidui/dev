#!/usr/bin/env python3
import sys
while '/share/python/lib64/python' in sys.path :
    sys.path.remove('/share/python/lib64/python') 
# the cluster got this wrong dir
from numpy import *
import time

################################################################################
#                                                                              #
# class Bpnn, a simple neural network with one hidden layer.                   #
#                                                                              #
# Use the tanh activity function for the hidden and output layer.              #
# Often only one output unit only, for binary problems.                        #
#                                                                              #
# Supervised-training with the back-propagation algorithm.                     #
#                                                                              #
# Partly based on code by Neil Schemenauer <nas@arctrix.com>                   #
################################################################################

class Bpnn:
    
    def __init__(self, ni=1, nh=1, no=1, eta = 0.001):
        # number of input, hidden, and output nodes
        self.ni = ni + 1 # +1 for the bias node
        self.nh = nh + 1 # +1 for the bias node 
        self.no = no     # often only one output unit!
        
        # activations of nodes
        self.act_i = ones(self.ni)
        self.act_h = ones(self.nh)
        self.act_o = ones(self.no)
        
        # error vectors 
        self.hid_error = zeros(self.nh-1)
        self.out_error = zeros(self.no  )
        
        # connection weights
        self.wei_i     = random.uniform(-1,1,(self.ni, self.nh-1))
        self.wei_o     = random.uniform(-1,1,(self.nh, self.no))
        
        # last changes in weights for momentum   
        self.chg_i     = zeros((self.ni, self.nh-1))
        self.chg_o     = zeros((self.nh, self.no)  )
        
        # temp change matrices
        self.tmp_chg_i = zeros((self.ni, self.nh-1))
        self.tmp_chg_o = zeros((self.nh, self.no  ))
        
        # the learning rate
        self.eta = eta 
        # the momentum factor can be set in the back_propagation function
    
    def save_net(self,file_name):
        ofile=open(file_name,'w')
        ofile.write(str(self.ni-1)+' '+str(self.nh-1)+' '+str(self.no)+' ')
        ofile.write(str(self.eta)+'\n') 
        
        for i in range(self.ni):
            for j in range(self.nh-1):
                ofile.write(str(self.wei_i[i][j])+' ')
        ofile.write('\n')
        
        for i in range(self.nh):
            for j in range(self.no):
                ofile.write(str(self.wei_o[i][j])+' ')
        ofile.write('\n')
        ofile.close()
    
    def read_net(self,file_name):
        ifile=open(file_name,'r')
        cols=ifile.readline().split()
        (ni,nh,no,eta)=(int(cols[0]),int(cols[1]),int(cols[2]),float(cols[3]))
        self.__init__(ni,nh,no,eta) 
        
        cols=ifile.readline().split()
        index=0;
        for i in range(self.ni):
            for j in range(self.nh-1):
                self.wei_i[i][j]=float(cols[index])
                index+=1
        
        cols=ifile.readline().split()
        index=0
        for i in range(self.nh):
            for j in range(self.no):
                self.wei_o[i][j]=float(cols[index])
                index+=1
        
        ifile.close()
    
    # the sigmoid function is y = tanh(x)
    def propagate(self,inputs):
        # the first units of the input and hidden layer are 1.0 after __init__
        self.act_i[1:] = inputs 
        # hidden activations
        self.act_h[1:] = tanh(dot(self.act_i,self.wei_i))
        # output activations
        self.act_o     = tanh(dot(self.act_h,self.wei_o))
        return self.act_o

    # squared error 
    def se_error(self,target):
        return ((self.act_o-target)**2).sum()
    
    # the derivative of the sigmoid function is 1 - y**2
    def back_propagate(self,target,alpha=0.9):
        # alpha is the momentum factor of training
        
        # calculate error terms for output
        self.out_error = (1.0 - self.act_o ** 2) * (target - self.act_o)
        
        # calculate error terms for hidden
        # the first unit, bias is allways 1.0 , so the error will allways be 0 
        self.hid_error = ((1.0-self.act_h**2) * dot(self.wei_o,self.out_error))\
                         [1:]
        
        # update output weights
        self.tmp_chg_o = self.out_error * \
                            reshape(self.act_h,(self.act_h.shape[0],1)) 
        self.tmp_chg_o *= self.eta 
        self.tmp_chg_o += alpha * self.chg_o
        self.wei_o += self.tmp_chg_o
        self.chg_o =  self.tmp_chg_o
        
        # update input weights
        self.tmp_chg_i = self.hid_error * \
                          reshape(self.act_i,(self.act_i.shape[0],1))
        self.tmp_chg_i *= self.eta
        self.tmp_chg_i += alpha * self.chg_i
        self.wei_i += self.tmp_chg_i
        self.chg_i =  self.tmp_chg_i

################################################################################
# end of class Bpnn
# test function    

if __name__ == '__main__':
    # the usual XOR test set
    xor_inputs  = (( 1, 1),(-1, 1),( 1,-1),(-1,-1))
    xor_targets = (-1, 1, 1, -1)
    
    random.seed()
    net=Bpnn(2,2,1,0.01)
    
    begin_time=time.clock()
    for i in range(1000):
        for j in range(4):
            net.propagate     (xor_inputs[j])
            net.back_propagate(xor_targets[j])
            print(i,net.se_error(xor_targets[j]))
    
    end_time=time.clock()
    print('Total training time: ',end_time-begin_time)
    
    avg_se_err=0;
    for i in range(4):
        print(net.propagate(xor_inputs[i]))
        avg_se_err+=net.se_error(xor_targets[i])
    avg_se_err/=4
    print('se_error: ',avg_se_err)
    
    # test save and load network file 
    net.save_net('t.net')
    net2=Bpnn()
    net2.read_net('t.net')
    print('\nload network from the t.net file')
    avg_se_err=0;
    for i in range(4):
        net2.propagate(xor_inputs[i])
        avg_se_err+=net2.se_error(xor_targets[i])
        print(xor_inputs[i],net2.propagate(xor_inputs[i]),\
              net2.se_error(xor_targets[i]))
    avg_se_err/=4
    print('se_error: ',avg_se_err)

#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk

################################################################################
#                               options                                        #
################################################################################

CXB_VERSION = '0.01'
CXB_WEBSITE = 'http://code.google.com/p/cxb/'

from options import *

################################################################################
#                           The Tab class                                      #
################################################################################

class CxbTab(ttk.Frame):
    
    def __init__(self, parent=None, file=None):
        self.fileName='';
        ttk.Frame.__init__(self, parent,border=0,relief=FLAT)
        self.pack(expand=YES, fill=BOTH)                 
        self.make_widgets()
    
    def make_widgets(self):
        
        statusBar = ttk.Frame(self)
        lineLabel = ttk.Label(statusBar,text="line1")
        lineLabel.pack(side=RIGHT)
        
        middlePanel = ttk.Frame(self)
        ysbar = ttk.Scrollbar(middlePanel)
        
        self.text = Text(middlePanel, 
            font=EDITOR['FONT'],wrap=EDITOR['WRAP'],
            height=EDITOR['HEIGHT'],width=EDITOR['WIDTH'],
            background=EDITOR['BACKGROUND_COLOR'],
            foreground=EDITOR['FOREGROUND_COLOR'],
            insertwidth=EDITOR['CURSOR_WIDTH'],
            insertbackgroun=EDITOR['CURSOR_COLOR'],
            relief=FLAT,bd=0,padx=3,pady=3,setgrid=1)
            
        ysbar.pack(side=RIGHT,expand=NO,fill=Y)
        self.text.pack(side=LEFT,expand=YES,fill=BOTH)
        
        bottomPanel=ttk.Frame(self)
        xsbar = ttk.Scrollbar(bottomPanel,orient=HORIZONTAL)
        grip  = ttk.Sizegrip(bottomPanel)
        
        grip.pack(side=RIGHT,expand=NO,fill=NONE)
        xsbar.pack(side=LEFT,expand=YES,fill=X)
        
        xsbar.config(command=self.text.xview)
        ysbar.config(command=self.text.yview)
        
        self.text.config(xscrollcommand=xsbar.set,yscrollcommand=ysbar.set)           
        
        statusBar.pack(side=TOP,expand=NO,fill=X)
        bottomPanel.pack(side=BOTTOM,expand=NO,fill=X)
        middlePanel.pack(side=LEFT,expand=YES,fill=BOTH)
        self.pack()
        
################################################################################
#                           the main script                                    #
################################################################################

if __name__ == '__main__':
    root = Tk()
    root.title(string='程序本v'+ CXB_VERSION)
    
    tabs = ttk.Notebook(root)
    
    tab1 = CxbTab(tabs,None); # first tab,
    tabs.add(tab1, text=tab1.fileName)
    
    tabs.pack(expand=YES,fill=BOTH)
    
    root.mainloop()
    
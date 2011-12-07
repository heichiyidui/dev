#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk

# how the hell can I get the c version of this thing?
################################################################################
#                           global options                                     #
################################################################################

CXB_VERSION = '0.01'
CXB_WEBSITE = 'http://code.google.com/p/cxb/'

from options import *

################################################################################
#                           The Tab class                                      #
################################################################################

class CxbTab(ttk.Frame):
    
    def __init__(self, parent=None, file=None):
        self.ifileName='';
        ttk.Frame.__init__(self, parent,border=0,relief=SUNKEN)
        self.pack(expand=YES, fill=BOTH)                 
        self.makewidgets()
        self.settext(file)
    
    def makewidgets(self):
        
        xsbar = ttk.Scrollbar(self,orient=HORIZONTAL)
        ysbar = ttk.Scrollbar(self)
#        statusbar = ttk.Label(self, text='untitled',
#                             font=DEFAULT_FONT,width=DEFAULT_WIDTH)
        
        text = Text(self, 
            font=EDITOR['FONT'],wrap=EDITOR['WRAP'],
            height=EDITOR['HEIGHT'],width=EDITOR['WIDTH'],
            background=EDITOR['BACKGROUND_COLOR'],
            foreground=EDITOR['FOREGROUND_COLOR'],
            insertwidth=EDITOR['CURSOR_WIDTH'],
            insertbackgroun=EDITOR['CURSOR_COLOR'],
            
            relief=FLAT,bd=0,
            padx=3,pady=3,setgrid=True)
        
        xsbar.config(command=text.xview)
        ysbar.config(command=text.yview)
        
        text.config(xscrollcommand=xsbar.set)           
        text.config(yscrollcommand=ysbar.set)           
        
#        statusbar.pack(side=TOP,fill=X)
        xsbar.pack(side=BOTTOM,fill=X)
        ysbar.pack(side=RIGHT, fill=Y)                 
        
        text.pack(side=LEFT, expand=YES, fill=BOTH)  
        self.text = text
    
    def settext(self, file=None):
        text=''
        if file: 
            text = open(file, 'r').read()
        self.text.delete('1.0', END)                 
        self.text.insert('1.0', text)               
        self.text.mark_set(INSERT, '1.0')          
        self.text.focus()                           
    
    def gettext(self):                             
        return self.text.get('1.0', END+'-1c') 

################################################################################
#                           the main script                                    #
################################################################################

if __name__ == '__main__':
    root = Tk()
    root.title(string='cxb version'+ CXB_VERSION)
    
    tabs = ttk.Notebook(root)
    
    tab1 = CxbTab(tabs,None); # first tab,
    tabs.add(tab1, text=tab1.ifileName)
    
    tabs.pack(expand=YES,fill=BOTH)
    
    root.mainloop()

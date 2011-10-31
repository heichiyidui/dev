#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk

################################################################################
#                           global options                                     #
################################################################################

CXB_VERSION - '0.01'
CXV_WEBSITE = 'http://code.google.com/p/cxb/'

from default import *

################################################################################
#                           The Tab class                                      #
################################################################################

class CxbTab(ttk.Frame):
    
    def __init__(self, parent=None, text='', file=None):
        ttk.Frame.__init__(self, parent)
        self.pack(expand=YES, fill=BOTH)                 
        self.makewidgets()
        self.settext(text, file)

    def makewidgets(self):
        
        xsbar = ttk.Scrollbar(self,orient=HORIZONTAL)
        ysbar = ttk.Scrollbar(self)
        statusbar = ttk.Label(self, text='untitled',
                             font=DEFAULT_FONT,width=DEFAULT_WIDTH)
        
        text = Text(self, wrap='none',
            font=DEFAULT_FONT, height=DEFAULT_HEIGHT,
            background=DEFAULT_BACKGROUND_COLOR,
            foreground=DEFAULT_FOREGROUND_COLOR,
            width=DEFAULT_WIDTH, relief=SUNKEN,
            insertwidth=DEFAULT_CURSOR_WIDTH,
            insertbackgroun=DEFAULT_CURSOR_COLOR,
            padx=4,pady=4,setgrid=True)
        
        xsbar.config(command=text.xview)
        ysbar.config(command=text.yview)
        
        text.config(xscrollcommand=xsbar.set)           
        text.config(yscrollcommand=ysbar.set)           
        
        statusbar.pack(side=TOP,fill=X)
        xsbar.pack(side=BOTTOM,fill=X)
        ysbar.pack(side=RIGHT, fill=Y)                 
        
        text.pack(side=LEFT, expand=YES, fill=BOTH)  
        self.text = text

    def settext(self, text='', file=None):
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
    root.title(string='程序本 v0.01')
    
    tabs = ttk.Notebook(root)
    tab1 = CxbTab(tabs); # first tab,
    tabs.add(tab1, text='One')
    
    tabs.pack(expand=YES,fill=BOTH)
    
    root.mainloop()

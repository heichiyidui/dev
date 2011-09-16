#!/usr/bin/env python3 
#...................................................................................

from tkinter import *
from tkinter import ttk

class cxb_tab(ttk.Frame):

    def __init__(self, parent=None, text='', file=None):
        ttk.Frame.__init__(self, parent)
        self.pack(expand=YES, fill=BOTH)                 
        self.makewidgets()
        self.settext(text, file)

    def makewidgets(self):
        
        xsbar = ttk.Scrollbar(self,orient=HORIZONTAL)
        ysbar = ttk.Scrollbar(self)
        label = ttk.Label(self, text = 'untitled')
        
        text = Text(self,height=42, width=82, relief=SUNKEN)
        
        xsbar.config(command=text.xview)
        ysbar.config(command=text.yview)
        
        text.config(xscrollcommand=xsbar.set)           
        text.config(yscrollcommand=ysbar.set)           
        
        label.pack(side=BOTTOM,fill=X)
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

root = Tk()
root.title('程序本')

tabs=ttk.Notebook(root)

tab1=cxb_tab(tabs)
tab2=cxb_tab(tabs)

tabs.add(tab1,text='tab1')
tabs.add(tab2,text='tab2')

tabs.pack()
root.mainloop()

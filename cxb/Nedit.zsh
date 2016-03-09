################################################################################
#                                                                              #
#                           NEdit, or Gedit?                                   #
#                                                                              #
################################################################################

################################################################################
# 1. The open file problem                                                     #
################################################################################

# Use gedit to open files it always opens an untitled document as well.
# sudo vi /usr/share/applications/gedit.desktop
# replace Exec=gedit %U
# with Exec=gedit $1 < /dev/null
# OK, now with version 3.14.0 it's fixed.

################################################################################
# 2. the GtkSourceView style                                                   #
################################################################################

mkdir -p ~/.local/share/gtksourceview-3.0/styles/
cp green.xml ~/.local/share/gtksourceview-3.0/styles/
# go preference and choose the Green color scheme

mkdir -p ~/.local/share/gtksourceview-3.0/language-specs/
cp c.lang python.lang zsh.lang ~/.local/share/gtksourceview-3.0/language-specs/

# at last, 
# got the regular expression lookahead right there for funtion coloring

# for cursor width and color 
mkdir -p  ~/.config/gtk-3.0/
cp gtk.css ~/.config/gtk-3.0/

################################################################################
# 3. the plugins                                                               #
################################################################################

# Install via YaST (I'm a opensuse user) the gedit-plugins package
# for the code comment plugin
# select the code comment plugin in Edit/Preferences/Plugins
# the default Ctrl+M and Ctrl+Shift+M shortcuts are difficult to use
# in /usr/lib64/gedit/plugins/codecomment.py
# change them to Ctrl+E and Ctrl+Shift+E

# use the smart space, snippet and word completion plugin as well

################################################################################
# 4. preference                                                                #
################################################################################

# View: select all 
# Editor: Tab 4 spaces
#    insert spaces
#    enable auto indention
#    no local backup
# Colour: Green
#
################################################################################
#                                                                              #
#                      Gedit is too slow in KDE? Kate?                         #
#                                                                              #
################################################################################

################################################################################
# 1. language syntax                                                           #
################################################################################

sudo cp cpp.xml    /usr/share/kde4/apps/katepart/syntax/ 
sudo cp python.xml /usr/share/kde4/apps/katepart/syntax/ 

################################################################################
#  2 kate setting                                                              #
################################################################################

mkdir -p ~/.kde4/share/apps/ktexteditor_snippets/data/
cp 0_C++\ Snippets-0.6.xml 0_Python\ Snippets-0.1.xml \
    ~/.kde4/share/apps/ktexteditor_snippets/data/

# switch off startup help

# setting -> show toolbar off

# config kate
#    application
#        general 
#            warn about foreign process modification -> on
#        plugin
#            build plugin -> on
#            kate snippets -> on
#            symbol viewer -> on
#        Kate Snippets
#            C++ snippets -> on
#            Python Snippets -> on
#        Symbol Viewer
#            Display functions parameters -> on
#    Editor
#        appearance
#            Show indentation lines -> on
#            Highlight range between selected brackets -> on
#            Show line numbers -> on
#        Font And Colour
#            set font to "Liberation Mono 11" 
#            import scheme kate_Normal.kateschema
#        Editing
#            using 4 spaces for indentation
#            backspace in leading blank space unindent on
#        Open/Save
#            local backup off

################################################################################
#                                                                              #
#                            sublime text editor?                              #
#                                                                              #
################################################################################
#
# $70 is a bit too much?
# It seems to be a damn good one. 
# 
# sublime text editor build 3047 Release Date: 27 June 2013
#
# download and unzip it into ~/bin/ 
ln -s ~/bin/sublime_text_3/sublime_text ~/bin/

# put this line into ~/.alias 
alias se='sublime_text '

#############
mkdir -p ~/.config/sublime-text-3/Packages/User/
cp Python.sublime-settings Green.tmTheme Preferences.sublime-settings \
   ~/.config/sublime-text-3/Packages/User/

################################################################################

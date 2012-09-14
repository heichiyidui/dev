
################################################################################
#                                                                              #
#                           NEdit, or Gedit?                                   #
#                                                                              #
################################################################################

################################################################################
# 1. the GtkSourceView style                                                   #
################################################################################

mkdir -p ~/.local/share/gtksourceview-3.0/styles/
cp green.xml ~/.local/share/gtksourceview-3.0/styles/

mkdir -p ~/.local/share/gtksourceview-3.0/language-specs/
cp c.lang python.lang zsh.lang ~/.local/share/gtksourceview-3.0/language-specs/

# at last, 
# got the regular expression lookahead right there for funtion coloring

################################################################################
# 2. the plugins                                                               #
################################################################################

# 2.1 install via YaST (I'm a opensuse user) the gedit-plugins package
# for the code comment plugin
# select the code comment plugin in Edit/Preferences/Plugins
# the default Ctrl+M and Ctrl+Shift+M shortcuts are difficult to use
# in /usr/lib64/gedit/plugins/codecomment.py
# change them to Ctrl+E and Ctrl+Shift+E

# Use gedit to open files it always opens an untitled document as well.
sudo vi /usr/share/applications/gedit.desktop
# replace Exec=gedit %U
# with Exec=gedit $1 < /dev/null

preference: 
Text Wrapping off

################################################################################
#                                                                              #
#                       Gedit is too slow in KDE? Kate?                        #
#                                                                              #
################################################################################

1. switch off startup help

2. setting -> show toolbar off

3. config kate
    application
        general 
            warn about foreign process modification -> on
        plugin
            build plugin -> on
            kate snippets -> on
            Tab bar -> on
            symbol viewer -> on
        Kate Snippets
            cp 0_C++ Snippets-0.6.xml and 0_Python Snippets-0.1.xml to \
            ~/.kde4/share/apps/ktexteditor_snippets/data/        
            C++ snippets -> on
            Python Snippets -> on
        Symbol Viewer
            Display functions parameters -> on
    Editor
        appearance
            Show indentation lines -> on
            Highlight range between selected brackets -> on
            Show line numbers -> on
        Font And Colour
            set font to "Efont Fixed Regular 12"

 in file /usr/share/kde4/apps/katepart/syntax/bash.xml
   add zsh to extensions in line 11
 in file /usr/share/kde4/apps/katepart/syntax/cpp.xml
   add size_t and string to types
   add <RegExpr attribute="Function" String="[_a-zA-Z][_a-zA-Z0-9]*\s*[(]" />
   to <contexts> <context attribute="Normal Text">
   add <itemData name="Function"    defStyleNum="dfFuntion" spellChecking="false"/>
   to <itemDatas>
 in file /usr/share/kde4/apps/katepart/syntax/python.xml
   in <contexts> <context name="Normal">, 
   BEFORE the normal regExpr line
   add <RegExpr attribute="Function" String="[_a-zA-Z][_a-zA-Z0-9]*\s*[(]" context="#stay"/>
   in <itemDatas>, after the normal text line,
   add <itemData name="Function" defStyleNum="dfFuntion" spellChecking="false"/>

	import scheme kate_Normal.kateschema

        Editing
            indentation 4 spaces
            backspace in leading to unindent
        Open/Save
            Don't autosave


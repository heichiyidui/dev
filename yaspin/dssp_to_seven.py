
MAP_8_TO_3={'X':'X',
            'E':'E','B':'E',
            'H':'H','G':'H',
            'C':'C','S':'C','T':'C','I':'C'}

MAP_7_to_i={'X':-1, 'C':0,'D':1,'E':2,'F':3,'G':4,'H':5,'I':6}

def dssp_to_seven(dssp):
    three_state_string = ''
    for dssp_char in dssp:
        three_state_string+=MAP_8_TO_3[dssp_char]
    
    seven_state_string = ''
    seven_state_string += three_state_string[0]
    for i in range(len(three_state_string)-2):
        ss_char = three_state_string[i+1]
        if three_state_string[i:i+3] in ('CHH','EHH'):
            ss_char='G'
        if three_state_string[i:i+3] in ('HHC','HHE'):
            ss_char='I'
        if three_state_string[i:i+3] in ('HEE','CEE'):
            ss_char='D'
        if three_state_string[i:i+3] in ('EEH','EEC'):
            ss_char='F'
        seven_state_string += ss_char
    seven_state_string += three_state_string[-1]
    
    seven_state_def = []
    for ss_char in seven_state_string:
        seven_state_def.append(MAP_7_to_i[ss_char])
    
    return seven_state_def
    


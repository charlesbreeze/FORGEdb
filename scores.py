index = 0

eqtl = 0

abc = 0

motif = 0

cato = 0

motifcato = 0

dnase = 0

hmark = 0

fdbscore = 0

import sys
sys.argv[1]

import fileinput
#with open('test.txt2') as a_file:
#for line in a_file:
for line in fileinput.input(sys.argv[1]):
            index = index + 1
            #print line,  # The comma to suppress the extra new line char
            if "rs" in line:
                index = 0
                print line,
            if index == 4:
                eqtl = float(line.rstrip())
            if index == 5:
                eqtl = float(line.rstrip()) + float(eqtl)
                if eqtl > 0.1:
                    eqtl = 2    
                #print eqtl
            if index == 6:
                abc = float(line.rstrip())
                if abc > 0.1:
                    abc = 2
                #print abc
            if index == 7:
                motif = float(line.rstrip())
                #print motif
            if index == 8:
                cato = float(line.rstrip())
                #print cato
                motifcato = float(motif) + float(cato)
                #print motifcato
            if index == 9:
                dnase = float(line.rstrip())
            if index == 10:
                hmark = float(line.rstrip())
                if hmark > 0.1:
                    hmark = 2
                #print hmark
            if index == 11:
                dnase = float(line.rstrip()) + float(dnase)
            if index == 12:
                dnase = float(line.rstrip()) + float(dnase)
            if index == 13:
                dnase = float(line.rstrip()) + float(dnase)
                if dnase > 0.1:
                    dnase = 2
                #print dnase
                fdbscore = float(eqtl) + float(abc) + float(motifcato) + float(hmark) + float(dnase)
                print fdbscore
                
                
                
                
                
                            

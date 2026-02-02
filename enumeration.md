Here I would like to keep an explanation of what I am doing in enumeration.ipynb.

I have reran BruteSwitching but now with only half the graphs. This is about twice as fast.

GM4: time for 7 vertices 2.2667720317840576
GM6: time for 7 vertices 0.4853804111480713
WQH6: time for 7 vertices 3.5667808055877686
AH6: time for 7 vertices 10.407105445861816

GM4: time for 8 vertices 66.30962491035461
GM6: time for 8 vertices 26.474975109100342
GM8: time for 8 vertices 1.3516762256622314
WQH6: time for 8 vertices 214.730149269104
WQH8: time for 8 vertices 27.68884539604187
AH^: time for 8 vertices 629.846931219101

Potential idea is to first generate all the graphs using nauty and then reading them from graph6 format, but I don't think that can be much faster.


With the new way generating way, generating 114048 graphs only takes 5 sec:
        switching sets          | cospectralmates | total    | time
WQH8+2: 2768256 = 72^2*267*2          5757             11514    14m40   
GM8+2:  114048 = 72^2*11*2            3049             6098     1m15  --> now with complements earlier it takes about 1m45, 40 sec of which are reducelisttoiso

GM6+3   170368 = 22^3*4*4            1244             2488     2m0 --> 1m11 zonder QTAQ, 58 zonder canonical label
        (5200 distinct graphs)
WQH6+3: 553696 = 22^3*13*4           3167            6212      4m39
AH8+0:  136= 136                     14               20       7m34  
AH8+1:  38080 = 32*1190               713             1226?      19s (if already have all switchingsets)
AH8+2:  2437120 = 32^2*1190*2         25302          50604?      26m20 (if already have all switchingsets)
AH6+2:  7680 = 16^2*15*2               216            414       4s
AH6+3:  245760 = 16^3*15*4             6215          11846      2m39

On cluster we got the following
See nds files. Reminder that if it says cluster behind it then there are still isomorphic ones in that list.


        switching sets          | cospectralmates | compiledlist   | time
GM6+4:  7,496,192 = 22^4*8*11          132,000 x           5,465,356        less than half an hour split in 2 groups of 100 - GM advantage: I generated 59636 per thing: 240000 (all enumerationN.g6 are from this one)
WQH8+2:   2*(2,768,256 = 72^2*267*2)?      11,514 x           47440           less than an hour split in 20
GM8+3:    32,845,824 = 72^3*22*4      589,084               9,202,394      less than two hours split in 72
GM4+6:  163,577,856 = 8^6*4*156         1,977,190       55,177,120      up to 1.5 hours per job when split in 800 jobs.
WQH6+4: 41,229,056 = 22^4*16*11          407,770           7,839,458
AH6+4:  10,813,440 = 16^4*15*11 
irred   1,441,792 = 16^4*2*11
Fano+2: 2048 = 4*16^2*2                                                         1 second
Fano+3: 65536 = 4*16^3*4                                                         few seconds


Again:
AH6+3 (done, got slightly smaller 11846 -> 11784)
+irred AH6 (done)
AH8+2


Compare 1/a(Gamma,Q) for WQH6:
computer | my counts
72              72 - empty
12              12 - 3K_2
12              12 - C_6
8               8 - K_4+K_2              
4               4 - C_4+2K_1 
2               2 - weird one 
8               8 - H
72              72 - K3,3
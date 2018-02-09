FUNCTION mwa_get_pointing_number,obs,string=string

delays=*obs.delays
delays=Fix(delays)
CASE 1 OF
    Min(delays EQ [0,5,10,15,1,6,11,16,2,7,12,17,3,8,13,18]):pointing_number=-5
    Min(delays EQ [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]):pointing_number=-4
    Min(delays EQ [0,3,6,9,0,3,6,9,0,3,6,9,0,3,6,9]):pointing_number=-3
    Min(delays EQ [0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]):pointing_number=-2
    Min(delays EQ [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]):pointing_number=-1
    Min(delays EQ [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]):pointing_number=0
    Min(delays EQ [15,10,5,0,16,11,6,1,17,12,7,2,18,13,8,3]):pointing_number=5
    Min(delays EQ [12,8,4,0,13,9,5,1,14,10,6,2,15,11,7,3]):pointing_number=4
    Min(delays EQ Reverse([0,3,6,9,0,3,6,9,0,3,6,9,0,3,6,9])):pointing_number=3
    Min(delays EQ Reverse([0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6])):pointing_number=2
    Min(delays EQ Reverse([0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3])):pointing_number=1
    ELSE: message, "Pointing number not within 5 pointings around zenith"
ENDCASE

IF Keyword_Set(string) THEN pointing_number=Strn(pointing_number) 
RETURN,pointing_number
END

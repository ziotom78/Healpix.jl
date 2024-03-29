resol = Healpix.Resolution(4)
highresol = Healpix.Resolution(2^29)

@test Healpix.pix2xyfRing(resol, 1) == (3, 3, 0)
@test Healpix.pix2xyfRing(resol, 2) == (3, 3, 1)
@test Healpix.pix2xyfRing(resol, 3) == (3, 3, 2)
@test Healpix.pix2xyfRing(resol, 4) == (3, 3, 3)
@test Healpix.pix2xyfRing(resol, 5) == (2, 3, 0)
@test Healpix.pix2xyfRing(resol, 6) == (3, 2, 0)
@test Healpix.pix2xyfRing(resol, 7) == (2, 3, 1)
@test Healpix.pix2xyfRing(resol, 8) == (3, 2, 1)
@test Healpix.pix2xyfRing(resol, 9) == (2, 3, 2)
@test Healpix.pix2xyfRing(resol, 10) == (3, 2, 2)
@test Healpix.pix2xyfRing(resol, 11) == (2, 3, 3)
@test Healpix.pix2xyfRing(resol, 12) == (3, 2, 3)
@test Healpix.pix2xyfRing(resol, 13) == (1, 3, 0)
@test Healpix.pix2xyfRing(resol, 14) == (2, 2, 0)
@test Healpix.pix2xyfRing(resol, 15) == (3, 1, 0)
@test Healpix.pix2xyfRing(resol, 16) == (1, 3, 1)
@test Healpix.pix2xyfRing(resol, 17) == (2, 2, 1)
@test Healpix.pix2xyfRing(resol, 18) == (3, 1, 1)
@test Healpix.pix2xyfRing(resol, 19) == (1, 3, 2)
@test Healpix.pix2xyfRing(resol, 20) == (2, 2, 2)
@test Healpix.pix2xyfRing(resol, 21) == (3, 1, 2)
@test Healpix.pix2xyfRing(resol, 22) == (1, 3, 3)
@test Healpix.pix2xyfRing(resol, 23) == (2, 2, 3)
@test Healpix.pix2xyfRing(resol, 24) == (3, 1, 3)
@test Healpix.pix2xyfRing(resol, 25) == (0, 3, 0)
@test Healpix.pix2xyfRing(resol, 26) == (1, 2, 0)
@test Healpix.pix2xyfRing(resol, 27) == (2, 1, 0)
@test Healpix.pix2xyfRing(resol, 28) == (3, 0, 0)
@test Healpix.pix2xyfRing(resol, 29) == (0, 3, 1)
@test Healpix.pix2xyfRing(resol, 30) == (1, 2, 1)
@test Healpix.pix2xyfRing(resol, 31) == (2, 1, 1)
@test Healpix.pix2xyfRing(resol, 32) == (3, 0, 1)
@test Healpix.pix2xyfRing(resol, 33) == (0, 3, 2)
@test Healpix.pix2xyfRing(resol, 34) == (1, 2, 2)
@test Healpix.pix2xyfRing(resol, 35) == (2, 1, 2)
@test Healpix.pix2xyfRing(resol, 36) == (3, 0, 2)
@test Healpix.pix2xyfRing(resol, 37) == (0, 3, 3)
@test Healpix.pix2xyfRing(resol, 38) == (1, 2, 3)
@test Healpix.pix2xyfRing(resol, 39) == (2, 1, 3)
@test Healpix.pix2xyfRing(resol, 40) == (3, 0, 3)
@test Healpix.pix2xyfRing(resol, 41) == (3, 3, 4)
@test Healpix.pix2xyfRing(resol, 42) == (0, 2, 0)
@test Healpix.pix2xyfRing(resol, 43) == (1, 1, 0)
@test Healpix.pix2xyfRing(resol, 44) == (2, 0, 0)
@test Healpix.pix2xyfRing(resol, 45) == (3, 3, 5)
@test Healpix.pix2xyfRing(resol, 46) == (0, 2, 1)
@test Healpix.pix2xyfRing(resol, 47) == (1, 1, 1)
@test Healpix.pix2xyfRing(resol, 48) == (2, 0, 1)
@test Healpix.pix2xyfRing(resol, 49) == (3, 3, 6)
@test Healpix.pix2xyfRing(resol, 50) == (0, 2, 2)
@test Healpix.pix2xyfRing(resol, 51) == (1, 1, 2)
@test Healpix.pix2xyfRing(resol, 52) == (2, 0, 2)
@test Healpix.pix2xyfRing(resol, 53) == (3, 3, 7)
@test Healpix.pix2xyfRing(resol, 54) == (0, 2, 3)
@test Healpix.pix2xyfRing(resol, 55) == (1, 1, 3)
@test Healpix.pix2xyfRing(resol, 56) == (2, 0, 3)
@test Healpix.pix2xyfRing(resol, 57) == (3, 2, 4)
@test Healpix.pix2xyfRing(resol, 58) == (0, 1, 0)
@test Healpix.pix2xyfRing(resol, 59) == (1, 0, 0)
@test Healpix.pix2xyfRing(resol, 60) == (2, 3, 5)
@test Healpix.pix2xyfRing(resol, 61) == (3, 2, 5)
@test Healpix.pix2xyfRing(resol, 62) == (0, 1, 1)
@test Healpix.pix2xyfRing(resol, 63) == (1, 0, 1)
@test Healpix.pix2xyfRing(resol, 64) == (2, 3, 6)
@test Healpix.pix2xyfRing(resol, 65) == (3, 2, 6)
@test Healpix.pix2xyfRing(resol, 66) == (0, 1, 2)
@test Healpix.pix2xyfRing(resol, 67) == (1, 0, 2)
@test Healpix.pix2xyfRing(resol, 68) == (2, 3, 7)
@test Healpix.pix2xyfRing(resol, 69) == (3, 2, 7)
@test Healpix.pix2xyfRing(resol, 70) == (0, 1, 3)
@test Healpix.pix2xyfRing(resol, 71) == (1, 0, 3)
@test Healpix.pix2xyfRing(resol, 72) == (2, 3, 4)
@test Healpix.pix2xyfRing(resol, 73) == (2, 2, 4)
@test Healpix.pix2xyfRing(resol, 74) == (3, 1, 4)
@test Healpix.pix2xyfRing(resol, 75) == (0, 0, 0)
@test Healpix.pix2xyfRing(resol, 76) == (1, 3, 5)
@test Healpix.pix2xyfRing(resol, 77) == (2, 2, 5)
@test Healpix.pix2xyfRing(resol, 78) == (3, 1, 5)
@test Healpix.pix2xyfRing(resol, 79) == (0, 0, 1)
@test Healpix.pix2xyfRing(resol, 80) == (1, 3, 6)
@test Healpix.pix2xyfRing(resol, 81) == (2, 2, 6)
@test Healpix.pix2xyfRing(resol, 82) == (3, 1, 6)
@test Healpix.pix2xyfRing(resol, 83) == (0, 0, 2)
@test Healpix.pix2xyfRing(resol, 84) == (1, 3, 7)
@test Healpix.pix2xyfRing(resol, 85) == (2, 2, 7)
@test Healpix.pix2xyfRing(resol, 86) == (3, 1, 7)
@test Healpix.pix2xyfRing(resol, 87) == (0, 0, 3)
@test Healpix.pix2xyfRing(resol, 88) == (1, 3, 4)
@test Healpix.pix2xyfRing(resol, 89) == (2, 1, 4)
@test Healpix.pix2xyfRing(resol, 90) == (3, 0, 4)
@test Healpix.pix2xyfRing(resol, 91) == (0, 3, 5)
@test Healpix.pix2xyfRing(resol, 92) == (1, 2, 5)
@test Healpix.pix2xyfRing(resol, 93) == (2, 1, 5)
@test Healpix.pix2xyfRing(resol, 94) == (3, 0, 5)
@test Healpix.pix2xyfRing(resol, 95) == (0, 3, 6)
@test Healpix.pix2xyfRing(resol, 96) == (1, 2, 6)
@test Healpix.pix2xyfRing(resol, 97) == (2, 1, 6)
@test Healpix.pix2xyfRing(resol, 98) == (3, 0, 6)
@test Healpix.pix2xyfRing(resol, 99) == (0, 3, 7)
@test Healpix.pix2xyfRing(resol, 100) == (1, 2, 7)
@test Healpix.pix2xyfRing(resol, 101) == (2, 1, 7)
@test Healpix.pix2xyfRing(resol, 102) == (3, 0, 7)
@test Healpix.pix2xyfRing(resol, 103) == (0, 3, 4)
@test Healpix.pix2xyfRing(resol, 104) == (1, 2, 4)
@test Healpix.pix2xyfRing(resol, 105) == (1, 1, 4)
@test Healpix.pix2xyfRing(resol, 106) == (2, 0, 4)
@test Healpix.pix2xyfRing(resol, 107) == (3, 3, 8)
@test Healpix.pix2xyfRing(resol, 108) == (0, 2, 5)
@test Healpix.pix2xyfRing(resol, 109) == (1, 1, 5)
@test Healpix.pix2xyfRing(resol, 110) == (2, 0, 5)
@test Healpix.pix2xyfRing(resol, 111) == (3, 3, 9)
@test Healpix.pix2xyfRing(resol, 112) == (0, 2, 6)
@test Healpix.pix2xyfRing(resol, 113) == (1, 1, 6)
@test Healpix.pix2xyfRing(resol, 114) == (2, 0, 6)
@test Healpix.pix2xyfRing(resol, 115) == (3, 3, 10)
@test Healpix.pix2xyfRing(resol, 116) == (0, 2, 7)
@test Healpix.pix2xyfRing(resol, 117) == (1, 1, 7)
@test Healpix.pix2xyfRing(resol, 118) == (2, 0, 7)
@test Healpix.pix2xyfRing(resol, 119) == (3, 3, 11)
@test Healpix.pix2xyfRing(resol, 120) == (0, 2, 4)
@test Healpix.pix2xyfRing(resol, 121) == (1, 0, 4)
@test Healpix.pix2xyfRing(resol, 122) == (2, 3, 8)
@test Healpix.pix2xyfRing(resol, 123) == (3, 2, 8)
@test Healpix.pix2xyfRing(resol, 124) == (0, 1, 5)
@test Healpix.pix2xyfRing(resol, 125) == (1, 0, 5)
@test Healpix.pix2xyfRing(resol, 126) == (2, 3, 9)
@test Healpix.pix2xyfRing(resol, 127) == (3, 2, 9)
@test Healpix.pix2xyfRing(resol, 128) == (0, 1, 6)
@test Healpix.pix2xyfRing(resol, 129) == (1, 0, 6)
@test Healpix.pix2xyfRing(resol, 130) == (2, 3, 10)
@test Healpix.pix2xyfRing(resol, 131) == (3, 2, 10)
@test Healpix.pix2xyfRing(resol, 132) == (0, 1, 7)
@test Healpix.pix2xyfRing(resol, 133) == (1, 0, 7)
@test Healpix.pix2xyfRing(resol, 134) == (2, 3, 11)
@test Healpix.pix2xyfRing(resol, 135) == (3, 2, 11)
@test Healpix.pix2xyfRing(resol, 136) == (0, 1, 4)
@test Healpix.pix2xyfRing(resol, 137) == (0, 0, 4)
@test Healpix.pix2xyfRing(resol, 138) == (1, 3, 8)
@test Healpix.pix2xyfRing(resol, 139) == (2, 2, 8)
@test Healpix.pix2xyfRing(resol, 140) == (3, 1, 8)
@test Healpix.pix2xyfRing(resol, 141) == (0, 0, 5)
@test Healpix.pix2xyfRing(resol, 142) == (1, 3, 9)
@test Healpix.pix2xyfRing(resol, 143) == (2, 2, 9)
@test Healpix.pix2xyfRing(resol, 144) == (3, 1, 9)
@test Healpix.pix2xyfRing(resol, 145) == (0, 0, 6)
@test Healpix.pix2xyfRing(resol, 146) == (1, 3, 10)
@test Healpix.pix2xyfRing(resol, 147) == (2, 2, 10)
@test Healpix.pix2xyfRing(resol, 148) == (3, 1, 10)
@test Healpix.pix2xyfRing(resol, 149) == (0, 0, 7)
@test Healpix.pix2xyfRing(resol, 150) == (1, 3, 11)
@test Healpix.pix2xyfRing(resol, 151) == (2, 2, 11)
@test Healpix.pix2xyfRing(resol, 152) == (3, 1, 11)
@test Healpix.pix2xyfRing(resol, 153) == (0, 3, 8)
@test Healpix.pix2xyfRing(resol, 154) == (1, 2, 8)
@test Healpix.pix2xyfRing(resol, 155) == (2, 1, 8)
@test Healpix.pix2xyfRing(resol, 156) == (3, 0, 8)
@test Healpix.pix2xyfRing(resol, 157) == (0, 3, 9)
@test Healpix.pix2xyfRing(resol, 158) == (1, 2, 9)
@test Healpix.pix2xyfRing(resol, 159) == (2, 1, 9)
@test Healpix.pix2xyfRing(resol, 160) == (3, 0, 9)
@test Healpix.pix2xyfRing(resol, 161) == (0, 3, 10)
@test Healpix.pix2xyfRing(resol, 162) == (1, 2, 10)
@test Healpix.pix2xyfRing(resol, 163) == (2, 1, 10)
@test Healpix.pix2xyfRing(resol, 164) == (3, 0, 10)
@test Healpix.pix2xyfRing(resol, 165) == (0, 3, 11)
@test Healpix.pix2xyfRing(resol, 166) == (1, 2, 11)
@test Healpix.pix2xyfRing(resol, 167) == (2, 1, 11)
@test Healpix.pix2xyfRing(resol, 168) == (3, 0, 11)
@test Healpix.pix2xyfRing(resol, 169) == (0, 2, 8)
@test Healpix.pix2xyfRing(resol, 170) == (1, 1, 8)
@test Healpix.pix2xyfRing(resol, 171) == (2, 0, 8)
@test Healpix.pix2xyfRing(resol, 172) == (0, 2, 9)
@test Healpix.pix2xyfRing(resol, 173) == (1, 1, 9)
@test Healpix.pix2xyfRing(resol, 174) == (2, 0, 9)
@test Healpix.pix2xyfRing(resol, 175) == (0, 2, 10)
@test Healpix.pix2xyfRing(resol, 176) == (1, 1, 10)
@test Healpix.pix2xyfRing(resol, 177) == (2, 0, 10)
@test Healpix.pix2xyfRing(resol, 178) == (0, 2, 11)
@test Healpix.pix2xyfRing(resol, 179) == (1, 1, 11)
@test Healpix.pix2xyfRing(resol, 180) == (2, 0, 11)
@test Healpix.pix2xyfRing(resol, 181) == (0, 1, 8)
@test Healpix.pix2xyfRing(resol, 182) == (1, 0, 8)
@test Healpix.pix2xyfRing(resol, 183) == (0, 1, 9)
@test Healpix.pix2xyfRing(resol, 184) == (1, 0, 9)
@test Healpix.pix2xyfRing(resol, 185) == (0, 1, 10)
@test Healpix.pix2xyfRing(resol, 186) == (1, 0, 10)
@test Healpix.pix2xyfRing(resol, 187) == (0, 1, 11)
@test Healpix.pix2xyfRing(resol, 188) == (1, 0, 11)
@test Healpix.pix2xyfRing(resol, 189) == (0, 0, 8)
@test Healpix.pix2xyfRing(resol, 190) == (0, 0, 9)
@test Healpix.pix2xyfRing(resol, 191) == (0, 0, 10)
@test Healpix.pix2xyfRing(resol, 192) == (0, 0, 11)

@test Healpix.pix2xyfRing(highresol, 1) == (536870911, 536870911, 0)
@test Healpix.pix2xyfRing(highresol, 2) == (536870911, 536870911, 1)
@test Healpix.pix2xyfRing(highresol, 3) == (536870911, 536870911, 2)
@test Healpix.pix2xyfRing(highresol, 4) == (536870911, 536870911, 3)
@test Healpix.pix2xyfRing(highresol, 5) == (536870910, 536870911, 0)
@test Healpix.pix2xyfRing(highresol, 6) == (536870911, 536870910, 0)
@test Healpix.pix2xyfRing(highresol, 7) == (536870910, 536870911, 1)
@test Healpix.pix2xyfRing(highresol, 8) == (536870911, 536870910, 1)
@test Healpix.pix2xyfRing(highresol, 9) == (536870910, 536870911, 2)
@test Healpix.pix2xyfRing(highresol, 10) == (536870911, 536870910, 2)
@test Healpix.pix2xyfRing(highresol, 864691128455135233) == (469762048, 469762047, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135234) == (469762049, 469762046, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135235) == (469762050, 469762045, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135236) == (469762051, 469762044, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135237) == (469762052, 469762043, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135238) == (469762053, 469762042, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135239) == (469762054, 469762041, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135240) == (469762055, 469762040, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135241) == (469762056, 469762039, 6)
@test Healpix.pix2xyfRing(highresol, 864691128455135242) == (469762057, 469762038, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270465) == (268435456, 268435455, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270466) == (268435457, 268435454, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270467) == (268435458, 268435453, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270468) == (268435459, 268435452, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270469) == (268435460, 268435451, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270470) == (268435461, 268435450, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270471) == (268435462, 268435449, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270472) == (268435463, 268435448, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270473) == (268435464, 268435447, 6)
@test Healpix.pix2xyfRing(highresol, 1729382256910270474) == (268435465, 268435446, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405697) == (67108864, 67108863, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405698) == (67108865, 67108862, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405699) == (67108866, 67108861, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405700) == (67108867, 67108860, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405701) == (67108868, 67108859, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405702) == (67108869, 67108858, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405703) == (67108870, 67108857, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405704) == (67108871, 67108856, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405705) == (67108872, 67108855, 6)
@test Healpix.pix2xyfRing(highresol, 2594073385365405706) == (67108873, 67108854, 6)
@test Healpix.pix2xyfRing(highresol, 3458764513820540918) == (1, 0, 8)
@test Healpix.pix2xyfRing(highresol, 3458764513820540919) == (0, 1, 9)
@test Healpix.pix2xyfRing(highresol, 3458764513820540920) == (1, 0, 9)
@test Healpix.pix2xyfRing(highresol, 3458764513820540921) == (0, 1, 10)
@test Healpix.pix2xyfRing(highresol, 3458764513820540922) == (1, 0, 10)
@test Healpix.pix2xyfRing(highresol, 3458764513820540923) == (0, 1, 11)
@test Healpix.pix2xyfRing(highresol, 3458764513820540924) == (1, 0, 11)
@test Healpix.pix2xyfRing(highresol, 3458764513820540925) == (0, 0, 8)
@test Healpix.pix2xyfRing(highresol, 3458764513820540926) == (0, 0, 9)
@test Healpix.pix2xyfRing(highresol, 3458764513820540927) == (0, 0, 10)

@test Healpix.xyf2pixRing(resol, 3, 3, 0) == 1
@test Healpix.xyf2pixRing(resol, 3, 3, 1) == 2
@test Healpix.xyf2pixRing(resol, 3, 3, 2) == 3
@test Healpix.xyf2pixRing(resol, 3, 3, 3) == 4
@test Healpix.xyf2pixRing(resol, 2, 3, 0) == 5
@test Healpix.xyf2pixRing(resol, 3, 2, 0) == 6
@test Healpix.xyf2pixRing(resol, 2, 3, 1) == 7
@test Healpix.xyf2pixRing(resol, 3, 2, 1) == 8
@test Healpix.xyf2pixRing(resol, 2, 3, 2) == 9
@test Healpix.xyf2pixRing(resol, 3, 2, 2) == 10
@test Healpix.xyf2pixRing(resol, 2, 3, 3) == 11
@test Healpix.xyf2pixRing(resol, 3, 2, 3) == 12
@test Healpix.xyf2pixRing(resol, 1, 3, 0) == 13
@test Healpix.xyf2pixRing(resol, 2, 2, 0) == 14
@test Healpix.xyf2pixRing(resol, 3, 1, 0) == 15
@test Healpix.xyf2pixRing(resol, 1, 3, 1) == 16
@test Healpix.xyf2pixRing(resol, 2, 2, 1) == 17
@test Healpix.xyf2pixRing(resol, 3, 1, 1) == 18
@test Healpix.xyf2pixRing(resol, 1, 3, 2) == 19
@test Healpix.xyf2pixRing(resol, 2, 2, 2) == 20
@test Healpix.xyf2pixRing(resol, 3, 1, 2) == 21
@test Healpix.xyf2pixRing(resol, 1, 3, 3) == 22
@test Healpix.xyf2pixRing(resol, 2, 2, 3) == 23
@test Healpix.xyf2pixRing(resol, 3, 1, 3) == 24
@test Healpix.xyf2pixRing(resol, 0, 3, 0) == 25
@test Healpix.xyf2pixRing(resol, 1, 2, 0) == 26
@test Healpix.xyf2pixRing(resol, 2, 1, 0) == 27
@test Healpix.xyf2pixRing(resol, 3, 0, 0) == 28
@test Healpix.xyf2pixRing(resol, 0, 3, 1) == 29
@test Healpix.xyf2pixRing(resol, 1, 2, 1) == 30
@test Healpix.xyf2pixRing(resol, 2, 1, 1) == 31
@test Healpix.xyf2pixRing(resol, 3, 0, 1) == 32
@test Healpix.xyf2pixRing(resol, 0, 3, 2) == 33
@test Healpix.xyf2pixRing(resol, 1, 2, 2) == 34
@test Healpix.xyf2pixRing(resol, 2, 1, 2) == 35
@test Healpix.xyf2pixRing(resol, 3, 0, 2) == 36
@test Healpix.xyf2pixRing(resol, 0, 3, 3) == 37
@test Healpix.xyf2pixRing(resol, 1, 2, 3) == 38
@test Healpix.xyf2pixRing(resol, 2, 1, 3) == 39
@test Healpix.xyf2pixRing(resol, 3, 0, 3) == 40
@test Healpix.xyf2pixRing(resol, 3, 3, 4) == 41
@test Healpix.xyf2pixRing(resol, 0, 2, 0) == 42
@test Healpix.xyf2pixRing(resol, 1, 1, 0) == 43
@test Healpix.xyf2pixRing(resol, 2, 0, 0) == 44
@test Healpix.xyf2pixRing(resol, 3, 3, 5) == 45
@test Healpix.xyf2pixRing(resol, 0, 2, 1) == 46
@test Healpix.xyf2pixRing(resol, 1, 1, 1) == 47
@test Healpix.xyf2pixRing(resol, 2, 0, 1) == 48
@test Healpix.xyf2pixRing(resol, 3, 3, 6) == 49
@test Healpix.xyf2pixRing(resol, 0, 2, 2) == 50
@test Healpix.xyf2pixRing(resol, 1, 1, 2) == 51
@test Healpix.xyf2pixRing(resol, 2, 0, 2) == 52
@test Healpix.xyf2pixRing(resol, 3, 3, 7) == 53
@test Healpix.xyf2pixRing(resol, 0, 2, 3) == 54
@test Healpix.xyf2pixRing(resol, 1, 1, 3) == 55
@test Healpix.xyf2pixRing(resol, 2, 0, 3) == 56
@test Healpix.xyf2pixRing(resol, 3, 2, 4) == 57
@test Healpix.xyf2pixRing(resol, 0, 1, 0) == 58
@test Healpix.xyf2pixRing(resol, 1, 0, 0) == 59
@test Healpix.xyf2pixRing(resol, 2, 3, 5) == 60
@test Healpix.xyf2pixRing(resol, 3, 2, 5) == 61
@test Healpix.xyf2pixRing(resol, 0, 1, 1) == 62
@test Healpix.xyf2pixRing(resol, 1, 0, 1) == 63
@test Healpix.xyf2pixRing(resol, 2, 3, 6) == 64
@test Healpix.xyf2pixRing(resol, 3, 2, 6) == 65
@test Healpix.xyf2pixRing(resol, 0, 1, 2) == 66
@test Healpix.xyf2pixRing(resol, 1, 0, 2) == 67
@test Healpix.xyf2pixRing(resol, 2, 3, 7) == 68
@test Healpix.xyf2pixRing(resol, 3, 2, 7) == 69
@test Healpix.xyf2pixRing(resol, 0, 1, 3) == 70
@test Healpix.xyf2pixRing(resol, 1, 0, 3) == 71
@test Healpix.xyf2pixRing(resol, 2, 3, 4) == 72
@test Healpix.xyf2pixRing(resol, 2, 2, 4) == 73
@test Healpix.xyf2pixRing(resol, 3, 1, 4) == 74
@test Healpix.xyf2pixRing(resol, 0, 0, 0) == 75
@test Healpix.xyf2pixRing(resol, 1, 3, 5) == 76
@test Healpix.xyf2pixRing(resol, 2, 2, 5) == 77
@test Healpix.xyf2pixRing(resol, 3, 1, 5) == 78
@test Healpix.xyf2pixRing(resol, 0, 0, 1) == 79
@test Healpix.xyf2pixRing(resol, 1, 3, 6) == 80
@test Healpix.xyf2pixRing(resol, 2, 2, 6) == 81
@test Healpix.xyf2pixRing(resol, 3, 1, 6) == 82
@test Healpix.xyf2pixRing(resol, 0, 0, 2) == 83
@test Healpix.xyf2pixRing(resol, 1, 3, 7) == 84
@test Healpix.xyf2pixRing(resol, 2, 2, 7) == 85
@test Healpix.xyf2pixRing(resol, 3, 1, 7) == 86
@test Healpix.xyf2pixRing(resol, 0, 0, 3) == 87
@test Healpix.xyf2pixRing(resol, 1, 3, 4) == 88
@test Healpix.xyf2pixRing(resol, 2, 1, 4) == 89
@test Healpix.xyf2pixRing(resol, 3, 0, 4) == 90
@test Healpix.xyf2pixRing(resol, 0, 3, 5) == 91
@test Healpix.xyf2pixRing(resol, 1, 2, 5) == 92
@test Healpix.xyf2pixRing(resol, 2, 1, 5) == 93
@test Healpix.xyf2pixRing(resol, 3, 0, 5) == 94
@test Healpix.xyf2pixRing(resol, 0, 3, 6) == 95
@test Healpix.xyf2pixRing(resol, 1, 2, 6) == 96
@test Healpix.xyf2pixRing(resol, 2, 1, 6) == 97
@test Healpix.xyf2pixRing(resol, 3, 0, 6) == 98
@test Healpix.xyf2pixRing(resol, 0, 3, 7) == 99
@test Healpix.xyf2pixRing(resol, 1, 2, 7) == 100
@test Healpix.xyf2pixRing(resol, 2, 1, 7) == 101
@test Healpix.xyf2pixRing(resol, 3, 0, 7) == 102
@test Healpix.xyf2pixRing(resol, 0, 3, 4) == 103
@test Healpix.xyf2pixRing(resol, 1, 2, 4) == 104
@test Healpix.xyf2pixRing(resol, 1, 1, 4) == 105
@test Healpix.xyf2pixRing(resol, 2, 0, 4) == 106
@test Healpix.xyf2pixRing(resol, 3, 3, 8) == 107
@test Healpix.xyf2pixRing(resol, 0, 2, 5) == 108
@test Healpix.xyf2pixRing(resol, 1, 1, 5) == 109
@test Healpix.xyf2pixRing(resol, 2, 0, 5) == 110
@test Healpix.xyf2pixRing(resol, 3, 3, 9) == 111
@test Healpix.xyf2pixRing(resol, 0, 2, 6) == 112
@test Healpix.xyf2pixRing(resol, 1, 1, 6) == 113
@test Healpix.xyf2pixRing(resol, 2, 0, 6) == 114
@test Healpix.xyf2pixRing(resol, 3, 3, 10) == 115
@test Healpix.xyf2pixRing(resol, 0, 2, 7) == 116
@test Healpix.xyf2pixRing(resol, 1, 1, 7) == 117
@test Healpix.xyf2pixRing(resol, 2, 0, 7) == 118
@test Healpix.xyf2pixRing(resol, 3, 3, 11) == 119
@test Healpix.xyf2pixRing(resol, 0, 2, 4) == 120
@test Healpix.xyf2pixRing(resol, 1, 0, 4) == 121
@test Healpix.xyf2pixRing(resol, 2, 3, 8) == 122
@test Healpix.xyf2pixRing(resol, 3, 2, 8) == 123
@test Healpix.xyf2pixRing(resol, 0, 1, 5) == 124
@test Healpix.xyf2pixRing(resol, 1, 0, 5) == 125
@test Healpix.xyf2pixRing(resol, 2, 3, 9) == 126
@test Healpix.xyf2pixRing(resol, 3, 2, 9) == 127
@test Healpix.xyf2pixRing(resol, 0, 1, 6) == 128
@test Healpix.xyf2pixRing(resol, 1, 0, 6) == 129
@test Healpix.xyf2pixRing(resol, 2, 3, 10) == 130
@test Healpix.xyf2pixRing(resol, 3, 2, 10) == 131
@test Healpix.xyf2pixRing(resol, 0, 1, 7) == 132
@test Healpix.xyf2pixRing(resol, 1, 0, 7) == 133
@test Healpix.xyf2pixRing(resol, 2, 3, 11) == 134
@test Healpix.xyf2pixRing(resol, 3, 2, 11) == 135
@test Healpix.xyf2pixRing(resol, 0, 1, 4) == 136
@test Healpix.xyf2pixRing(resol, 0, 0, 4) == 137
@test Healpix.xyf2pixRing(resol, 1, 3, 8) == 138
@test Healpix.xyf2pixRing(resol, 2, 2, 8) == 139
@test Healpix.xyf2pixRing(resol, 3, 1, 8) == 140
@test Healpix.xyf2pixRing(resol, 0, 0, 5) == 141
@test Healpix.xyf2pixRing(resol, 1, 3, 9) == 142
@test Healpix.xyf2pixRing(resol, 2, 2, 9) == 143
@test Healpix.xyf2pixRing(resol, 3, 1, 9) == 144
@test Healpix.xyf2pixRing(resol, 0, 0, 6) == 145
@test Healpix.xyf2pixRing(resol, 1, 3, 10) == 146
@test Healpix.xyf2pixRing(resol, 2, 2, 10) == 147
@test Healpix.xyf2pixRing(resol, 3, 1, 10) == 148
@test Healpix.xyf2pixRing(resol, 0, 0, 7) == 149
@test Healpix.xyf2pixRing(resol, 1, 3, 11) == 150
@test Healpix.xyf2pixRing(resol, 2, 2, 11) == 151
@test Healpix.xyf2pixRing(resol, 3, 1, 11) == 152
@test Healpix.xyf2pixRing(resol, 0, 3, 8) == 153
@test Healpix.xyf2pixRing(resol, 1, 2, 8) == 154
@test Healpix.xyf2pixRing(resol, 2, 1, 8) == 155
@test Healpix.xyf2pixRing(resol, 3, 0, 8) == 156
@test Healpix.xyf2pixRing(resol, 0, 3, 9) == 157
@test Healpix.xyf2pixRing(resol, 1, 2, 9) == 158
@test Healpix.xyf2pixRing(resol, 2, 1, 9) == 159
@test Healpix.xyf2pixRing(resol, 3, 0, 9) == 160
@test Healpix.xyf2pixRing(resol, 0, 3, 10) == 161
@test Healpix.xyf2pixRing(resol, 1, 2, 10) == 162
@test Healpix.xyf2pixRing(resol, 2, 1, 10) == 163
@test Healpix.xyf2pixRing(resol, 3, 0, 10) == 164
@test Healpix.xyf2pixRing(resol, 0, 3, 11) == 165
@test Healpix.xyf2pixRing(resol, 1, 2, 11) == 166
@test Healpix.xyf2pixRing(resol, 2, 1, 11) == 167
@test Healpix.xyf2pixRing(resol, 3, 0, 11) == 168
@test Healpix.xyf2pixRing(resol, 0, 2, 8) == 169
@test Healpix.xyf2pixRing(resol, 1, 1, 8) == 170
@test Healpix.xyf2pixRing(resol, 2, 0, 8) == 171
@test Healpix.xyf2pixRing(resol, 0, 2, 9) == 172
@test Healpix.xyf2pixRing(resol, 1, 1, 9) == 173
@test Healpix.xyf2pixRing(resol, 2, 0, 9) == 174
@test Healpix.xyf2pixRing(resol, 0, 2, 10) == 175
@test Healpix.xyf2pixRing(resol, 1, 1, 10) == 176
@test Healpix.xyf2pixRing(resol, 2, 0, 10) == 177
@test Healpix.xyf2pixRing(resol, 0, 2, 11) == 178
@test Healpix.xyf2pixRing(resol, 1, 1, 11) == 179
@test Healpix.xyf2pixRing(resol, 2, 0, 11) == 180
@test Healpix.xyf2pixRing(resol, 0, 1, 8) == 181
@test Healpix.xyf2pixRing(resol, 1, 0, 8) == 182
@test Healpix.xyf2pixRing(resol, 0, 1, 9) == 183
@test Healpix.xyf2pixRing(resol, 1, 0, 9) == 184
@test Healpix.xyf2pixRing(resol, 0, 1, 10) == 185
@test Healpix.xyf2pixRing(resol, 1, 0, 10) == 186
@test Healpix.xyf2pixRing(resol, 0, 1, 11) == 187
@test Healpix.xyf2pixRing(resol, 1, 0, 11) == 188
@test Healpix.xyf2pixRing(resol, 0, 0, 8) == 189
@test Healpix.xyf2pixRing(resol, 0, 0, 9) == 190
@test Healpix.xyf2pixRing(resol, 0, 0, 10) == 191
@test Healpix.xyf2pixRing(resol, 0, 0, 11) == 192

@test Healpix.xyf2pixRing(highresol, 536870911, 536870911, 0) == 1
@test Healpix.xyf2pixRing(highresol, 536870911, 536870911, 1) == 2
@test Healpix.xyf2pixRing(highresol, 536870911, 536870911, 2) == 3
@test Healpix.xyf2pixRing(highresol, 536870911, 536870911, 3) == 4
@test Healpix.xyf2pixRing(highresol, 536870910, 536870911, 0) == 5
@test Healpix.xyf2pixRing(highresol, 536870911, 536870910, 0) == 6
@test Healpix.xyf2pixRing(highresol, 536870910, 536870911, 1) == 7
@test Healpix.xyf2pixRing(highresol, 536870911, 536870910, 1) == 8
@test Healpix.xyf2pixRing(highresol, 536870910, 536870911, 2) == 9
@test Healpix.xyf2pixRing(highresol, 536870911, 536870910, 2) == 10
@test Healpix.xyf2pixRing(highresol, 469762048, 469762047, 6) == 864691128455135233
@test Healpix.xyf2pixRing(highresol, 469762049, 469762046, 6) == 864691128455135234
@test Healpix.xyf2pixRing(highresol, 469762050, 469762045, 6) == 864691128455135235
@test Healpix.xyf2pixRing(highresol, 469762051, 469762044, 6) == 864691128455135236
@test Healpix.xyf2pixRing(highresol, 469762052, 469762043, 6) == 864691128455135237
@test Healpix.xyf2pixRing(highresol, 469762053, 469762042, 6) == 864691128455135238
@test Healpix.xyf2pixRing(highresol, 469762054, 469762041, 6) == 864691128455135239
@test Healpix.xyf2pixRing(highresol, 469762055, 469762040, 6) == 864691128455135240
@test Healpix.xyf2pixRing(highresol, 469762056, 469762039, 6) == 864691128455135241
@test Healpix.xyf2pixRing(highresol, 469762057, 469762038, 6) == 864691128455135242
@test Healpix.xyf2pixRing(highresol, 268435456, 268435455, 6) == 1729382256910270465
@test Healpix.xyf2pixRing(highresol, 268435457, 268435454, 6) == 1729382256910270466
@test Healpix.xyf2pixRing(highresol, 268435458, 268435453, 6) == 1729382256910270467
@test Healpix.xyf2pixRing(highresol, 268435459, 268435452, 6) == 1729382256910270468
@test Healpix.xyf2pixRing(highresol, 268435460, 268435451, 6) == 1729382256910270469
@test Healpix.xyf2pixRing(highresol, 268435461, 268435450, 6) == 1729382256910270470
@test Healpix.xyf2pixRing(highresol, 268435462, 268435449, 6) == 1729382256910270471
@test Healpix.xyf2pixRing(highresol, 268435463, 268435448, 6) == 1729382256910270472
@test Healpix.xyf2pixRing(highresol, 268435464, 268435447, 6) == 1729382256910270473
@test Healpix.xyf2pixRing(highresol, 268435465, 268435446, 6) == 1729382256910270474
@test Healpix.xyf2pixRing(highresol, 67108864, 67108863, 6) == 2594073385365405697
@test Healpix.xyf2pixRing(highresol, 67108865, 67108862, 6) == 2594073385365405698
@test Healpix.xyf2pixRing(highresol, 67108866, 67108861, 6) == 2594073385365405699
@test Healpix.xyf2pixRing(highresol, 67108867, 67108860, 6) == 2594073385365405700
@test Healpix.xyf2pixRing(highresol, 67108868, 67108859, 6) == 2594073385365405701
@test Healpix.xyf2pixRing(highresol, 67108869, 67108858, 6) == 2594073385365405702
@test Healpix.xyf2pixRing(highresol, 67108870, 67108857, 6) == 2594073385365405703
@test Healpix.xyf2pixRing(highresol, 67108871, 67108856, 6) == 2594073385365405704
@test Healpix.xyf2pixRing(highresol, 67108872, 67108855, 6) == 2594073385365405705
@test Healpix.xyf2pixRing(highresol, 67108873, 67108854, 6) == 2594073385365405706
@test Healpix.xyf2pixRing(highresol, 1, 0, 8) == 3458764513820540918
@test Healpix.xyf2pixRing(highresol, 0, 1, 9) == 3458764513820540919
@test Healpix.xyf2pixRing(highresol, 1, 0, 9) == 3458764513820540920
@test Healpix.xyf2pixRing(highresol, 0, 1, 10) == 3458764513820540921
@test Healpix.xyf2pixRing(highresol, 1, 0, 10) == 3458764513820540922
@test Healpix.xyf2pixRing(highresol, 0, 1, 11) == 3458764513820540923
@test Healpix.xyf2pixRing(highresol, 1, 0, 11) == 3458764513820540924
@test Healpix.xyf2pixRing(highresol, 0, 0, 8) == 3458764513820540925
@test Healpix.xyf2pixRing(highresol, 0, 0, 9) == 3458764513820540926
@test Healpix.xyf2pixRing(highresol, 0, 0, 10) == 3458764513820540927

@test Healpix.pix2xyfNest(resol, 1) == (0, 0, 0)
@test Healpix.pix2xyfNest(resol, 2) == (1, 0, 0)
@test Healpix.pix2xyfNest(resol, 3) == (0, 1, 0)
@test Healpix.pix2xyfNest(resol, 4) == (1, 1, 0)
@test Healpix.pix2xyfNest(resol, 5) == (2, 0, 0)
@test Healpix.pix2xyfNest(resol, 6) == (3, 0, 0)
@test Healpix.pix2xyfNest(resol, 7) == (2, 1, 0)
@test Healpix.pix2xyfNest(resol, 8) == (3, 1, 0)
@test Healpix.pix2xyfNest(resol, 9) == (0, 2, 0)
@test Healpix.pix2xyfNest(resol, 10) == (1, 2, 0)
@test Healpix.pix2xyfNest(resol, 11) == (0, 3, 0)
@test Healpix.pix2xyfNest(resol, 12) == (1, 3, 0)
@test Healpix.pix2xyfNest(resol, 13) == (2, 2, 0)
@test Healpix.pix2xyfNest(resol, 14) == (3, 2, 0)
@test Healpix.pix2xyfNest(resol, 15) == (2, 3, 0)
@test Healpix.pix2xyfNest(resol, 16) == (3, 3, 0)
@test Healpix.pix2xyfNest(resol, 17) == (0, 0, 1)
@test Healpix.pix2xyfNest(resol, 18) == (1, 0, 1)
@test Healpix.pix2xyfNest(resol, 19) == (0, 1, 1)
@test Healpix.pix2xyfNest(resol, 20) == (1, 1, 1)
@test Healpix.pix2xyfNest(resol, 21) == (2, 0, 1)
@test Healpix.pix2xyfNest(resol, 22) == (3, 0, 1)
@test Healpix.pix2xyfNest(resol, 23) == (2, 1, 1)
@test Healpix.pix2xyfNest(resol, 24) == (3, 1, 1)
@test Healpix.pix2xyfNest(resol, 25) == (0, 2, 1)
@test Healpix.pix2xyfNest(resol, 26) == (1, 2, 1)
@test Healpix.pix2xyfNest(resol, 27) == (0, 3, 1)
@test Healpix.pix2xyfNest(resol, 28) == (1, 3, 1)
@test Healpix.pix2xyfNest(resol, 29) == (2, 2, 1)
@test Healpix.pix2xyfNest(resol, 30) == (3, 2, 1)
@test Healpix.pix2xyfNest(resol, 31) == (2, 3, 1)
@test Healpix.pix2xyfNest(resol, 32) == (3, 3, 1)
@test Healpix.pix2xyfNest(resol, 33) == (0, 0, 2)
@test Healpix.pix2xyfNest(resol, 34) == (1, 0, 2)
@test Healpix.pix2xyfNest(resol, 35) == (0, 1, 2)
@test Healpix.pix2xyfNest(resol, 36) == (1, 1, 2)
@test Healpix.pix2xyfNest(resol, 37) == (2, 0, 2)
@test Healpix.pix2xyfNest(resol, 38) == (3, 0, 2)
@test Healpix.pix2xyfNest(resol, 39) == (2, 1, 2)
@test Healpix.pix2xyfNest(resol, 40) == (3, 1, 2)
@test Healpix.pix2xyfNest(resol, 41) == (0, 2, 2)
@test Healpix.pix2xyfNest(resol, 42) == (1, 2, 2)
@test Healpix.pix2xyfNest(resol, 43) == (0, 3, 2)
@test Healpix.pix2xyfNest(resol, 44) == (1, 3, 2)
@test Healpix.pix2xyfNest(resol, 45) == (2, 2, 2)
@test Healpix.pix2xyfNest(resol, 46) == (3, 2, 2)
@test Healpix.pix2xyfNest(resol, 47) == (2, 3, 2)
@test Healpix.pix2xyfNest(resol, 48) == (3, 3, 2)
@test Healpix.pix2xyfNest(resol, 49) == (0, 0, 3)
@test Healpix.pix2xyfNest(resol, 50) == (1, 0, 3)
@test Healpix.pix2xyfNest(resol, 51) == (0, 1, 3)
@test Healpix.pix2xyfNest(resol, 52) == (1, 1, 3)
@test Healpix.pix2xyfNest(resol, 53) == (2, 0, 3)
@test Healpix.pix2xyfNest(resol, 54) == (3, 0, 3)
@test Healpix.pix2xyfNest(resol, 55) == (2, 1, 3)
@test Healpix.pix2xyfNest(resol, 56) == (3, 1, 3)
@test Healpix.pix2xyfNest(resol, 57) == (0, 2, 3)
@test Healpix.pix2xyfNest(resol, 58) == (1, 2, 3)
@test Healpix.pix2xyfNest(resol, 59) == (0, 3, 3)
@test Healpix.pix2xyfNest(resol, 60) == (1, 3, 3)
@test Healpix.pix2xyfNest(resol, 61) == (2, 2, 3)
@test Healpix.pix2xyfNest(resol, 62) == (3, 2, 3)
@test Healpix.pix2xyfNest(resol, 63) == (2, 3, 3)
@test Healpix.pix2xyfNest(resol, 64) == (3, 3, 3)
@test Healpix.pix2xyfNest(resol, 65) == (0, 0, 4)
@test Healpix.pix2xyfNest(resol, 66) == (1, 0, 4)
@test Healpix.pix2xyfNest(resol, 67) == (0, 1, 4)
@test Healpix.pix2xyfNest(resol, 68) == (1, 1, 4)
@test Healpix.pix2xyfNest(resol, 69) == (2, 0, 4)
@test Healpix.pix2xyfNest(resol, 70) == (3, 0, 4)
@test Healpix.pix2xyfNest(resol, 71) == (2, 1, 4)
@test Healpix.pix2xyfNest(resol, 72) == (3, 1, 4)
@test Healpix.pix2xyfNest(resol, 73) == (0, 2, 4)
@test Healpix.pix2xyfNest(resol, 74) == (1, 2, 4)
@test Healpix.pix2xyfNest(resol, 75) == (0, 3, 4)
@test Healpix.pix2xyfNest(resol, 76) == (1, 3, 4)
@test Healpix.pix2xyfNest(resol, 77) == (2, 2, 4)
@test Healpix.pix2xyfNest(resol, 78) == (3, 2, 4)
@test Healpix.pix2xyfNest(resol, 79) == (2, 3, 4)
@test Healpix.pix2xyfNest(resol, 80) == (3, 3, 4)
@test Healpix.pix2xyfNest(resol, 81) == (0, 0, 5)
@test Healpix.pix2xyfNest(resol, 82) == (1, 0, 5)
@test Healpix.pix2xyfNest(resol, 83) == (0, 1, 5)
@test Healpix.pix2xyfNest(resol, 84) == (1, 1, 5)
@test Healpix.pix2xyfNest(resol, 85) == (2, 0, 5)
@test Healpix.pix2xyfNest(resol, 86) == (3, 0, 5)
@test Healpix.pix2xyfNest(resol, 87) == (2, 1, 5)
@test Healpix.pix2xyfNest(resol, 88) == (3, 1, 5)
@test Healpix.pix2xyfNest(resol, 89) == (0, 2, 5)
@test Healpix.pix2xyfNest(resol, 90) == (1, 2, 5)
@test Healpix.pix2xyfNest(resol, 91) == (0, 3, 5)
@test Healpix.pix2xyfNest(resol, 92) == (1, 3, 5)
@test Healpix.pix2xyfNest(resol, 93) == (2, 2, 5)
@test Healpix.pix2xyfNest(resol, 94) == (3, 2, 5)
@test Healpix.pix2xyfNest(resol, 95) == (2, 3, 5)
@test Healpix.pix2xyfNest(resol, 96) == (3, 3, 5)
@test Healpix.pix2xyfNest(resol, 97) == (0, 0, 6)
@test Healpix.pix2xyfNest(resol, 98) == (1, 0, 6)
@test Healpix.pix2xyfNest(resol, 99) == (0, 1, 6)
@test Healpix.pix2xyfNest(resol, 100) == (1, 1, 6)
@test Healpix.pix2xyfNest(resol, 101) == (2, 0, 6)
@test Healpix.pix2xyfNest(resol, 102) == (3, 0, 6)
@test Healpix.pix2xyfNest(resol, 103) == (2, 1, 6)
@test Healpix.pix2xyfNest(resol, 104) == (3, 1, 6)
@test Healpix.pix2xyfNest(resol, 105) == (0, 2, 6)
@test Healpix.pix2xyfNest(resol, 106) == (1, 2, 6)
@test Healpix.pix2xyfNest(resol, 107) == (0, 3, 6)
@test Healpix.pix2xyfNest(resol, 108) == (1, 3, 6)
@test Healpix.pix2xyfNest(resol, 109) == (2, 2, 6)
@test Healpix.pix2xyfNest(resol, 110) == (3, 2, 6)
@test Healpix.pix2xyfNest(resol, 111) == (2, 3, 6)
@test Healpix.pix2xyfNest(resol, 112) == (3, 3, 6)
@test Healpix.pix2xyfNest(resol, 113) == (0, 0, 7)
@test Healpix.pix2xyfNest(resol, 114) == (1, 0, 7)
@test Healpix.pix2xyfNest(resol, 115) == (0, 1, 7)
@test Healpix.pix2xyfNest(resol, 116) == (1, 1, 7)
@test Healpix.pix2xyfNest(resol, 117) == (2, 0, 7)
@test Healpix.pix2xyfNest(resol, 118) == (3, 0, 7)
@test Healpix.pix2xyfNest(resol, 119) == (2, 1, 7)
@test Healpix.pix2xyfNest(resol, 120) == (3, 1, 7)
@test Healpix.pix2xyfNest(resol, 121) == (0, 2, 7)
@test Healpix.pix2xyfNest(resol, 122) == (1, 2, 7)
@test Healpix.pix2xyfNest(resol, 123) == (0, 3, 7)
@test Healpix.pix2xyfNest(resol, 124) == (1, 3, 7)
@test Healpix.pix2xyfNest(resol, 125) == (2, 2, 7)
@test Healpix.pix2xyfNest(resol, 126) == (3, 2, 7)
@test Healpix.pix2xyfNest(resol, 127) == (2, 3, 7)
@test Healpix.pix2xyfNest(resol, 128) == (3, 3, 7)
@test Healpix.pix2xyfNest(resol, 129) == (0, 0, 8)
@test Healpix.pix2xyfNest(resol, 130) == (1, 0, 8)
@test Healpix.pix2xyfNest(resol, 131) == (0, 1, 8)
@test Healpix.pix2xyfNest(resol, 132) == (1, 1, 8)
@test Healpix.pix2xyfNest(resol, 133) == (2, 0, 8)
@test Healpix.pix2xyfNest(resol, 134) == (3, 0, 8)
@test Healpix.pix2xyfNest(resol, 135) == (2, 1, 8)
@test Healpix.pix2xyfNest(resol, 136) == (3, 1, 8)
@test Healpix.pix2xyfNest(resol, 137) == (0, 2, 8)
@test Healpix.pix2xyfNest(resol, 138) == (1, 2, 8)
@test Healpix.pix2xyfNest(resol, 139) == (0, 3, 8)
@test Healpix.pix2xyfNest(resol, 140) == (1, 3, 8)
@test Healpix.pix2xyfNest(resol, 141) == (2, 2, 8)
@test Healpix.pix2xyfNest(resol, 142) == (3, 2, 8)
@test Healpix.pix2xyfNest(resol, 143) == (2, 3, 8)
@test Healpix.pix2xyfNest(resol, 144) == (3, 3, 8)
@test Healpix.pix2xyfNest(resol, 145) == (0, 0, 9)
@test Healpix.pix2xyfNest(resol, 146) == (1, 0, 9)
@test Healpix.pix2xyfNest(resol, 147) == (0, 1, 9)
@test Healpix.pix2xyfNest(resol, 148) == (1, 1, 9)
@test Healpix.pix2xyfNest(resol, 149) == (2, 0, 9)
@test Healpix.pix2xyfNest(resol, 150) == (3, 0, 9)
@test Healpix.pix2xyfNest(resol, 151) == (2, 1, 9)
@test Healpix.pix2xyfNest(resol, 152) == (3, 1, 9)
@test Healpix.pix2xyfNest(resol, 153) == (0, 2, 9)
@test Healpix.pix2xyfNest(resol, 154) == (1, 2, 9)
@test Healpix.pix2xyfNest(resol, 155) == (0, 3, 9)
@test Healpix.pix2xyfNest(resol, 156) == (1, 3, 9)
@test Healpix.pix2xyfNest(resol, 157) == (2, 2, 9)
@test Healpix.pix2xyfNest(resol, 158) == (3, 2, 9)
@test Healpix.pix2xyfNest(resol, 159) == (2, 3, 9)
@test Healpix.pix2xyfNest(resol, 160) == (3, 3, 9)
@test Healpix.pix2xyfNest(resol, 161) == (0, 0, 10)
@test Healpix.pix2xyfNest(resol, 162) == (1, 0, 10)
@test Healpix.pix2xyfNest(resol, 163) == (0, 1, 10)
@test Healpix.pix2xyfNest(resol, 164) == (1, 1, 10)
@test Healpix.pix2xyfNest(resol, 165) == (2, 0, 10)
@test Healpix.pix2xyfNest(resol, 166) == (3, 0, 10)
@test Healpix.pix2xyfNest(resol, 167) == (2, 1, 10)
@test Healpix.pix2xyfNest(resol, 168) == (3, 1, 10)
@test Healpix.pix2xyfNest(resol, 169) == (0, 2, 10)
@test Healpix.pix2xyfNest(resol, 170) == (1, 2, 10)
@test Healpix.pix2xyfNest(resol, 171) == (0, 3, 10)
@test Healpix.pix2xyfNest(resol, 172) == (1, 3, 10)
@test Healpix.pix2xyfNest(resol, 173) == (2, 2, 10)
@test Healpix.pix2xyfNest(resol, 174) == (3, 2, 10)
@test Healpix.pix2xyfNest(resol, 175) == (2, 3, 10)
@test Healpix.pix2xyfNest(resol, 176) == (3, 3, 10)
@test Healpix.pix2xyfNest(resol, 177) == (0, 0, 11)
@test Healpix.pix2xyfNest(resol, 178) == (1, 0, 11)
@test Healpix.pix2xyfNest(resol, 179) == (0, 1, 11)
@test Healpix.pix2xyfNest(resol, 180) == (1, 1, 11)
@test Healpix.pix2xyfNest(resol, 181) == (2, 0, 11)
@test Healpix.pix2xyfNest(resol, 182) == (3, 0, 11)
@test Healpix.pix2xyfNest(resol, 183) == (2, 1, 11)
@test Healpix.pix2xyfNest(resol, 184) == (3, 1, 11)
@test Healpix.pix2xyfNest(resol, 185) == (0, 2, 11)
@test Healpix.pix2xyfNest(resol, 186) == (1, 2, 11)
@test Healpix.pix2xyfNest(resol, 187) == (0, 3, 11)
@test Healpix.pix2xyfNest(resol, 188) == (1, 3, 11)
@test Healpix.pix2xyfNest(resol, 189) == (2, 2, 11)
@test Healpix.pix2xyfNest(resol, 190) == (3, 2, 11)
@test Healpix.pix2xyfNest(resol, 191) == (2, 3, 11)
@test Healpix.pix2xyfNest(resol, 192) == (3, 3, 11)

@test Healpix.xyf2pixNest(resol, 0, 0, 0) == 1
@test Healpix.xyf2pixNest(resol, 1, 0, 0) == 2
@test Healpix.xyf2pixNest(resol, 0, 1, 0) == 3
@test Healpix.xyf2pixNest(resol, 1, 1, 0) == 4
@test Healpix.xyf2pixNest(resol, 2, 0, 0) == 5
@test Healpix.xyf2pixNest(resol, 3, 0, 0) == 6
@test Healpix.xyf2pixNest(resol, 2, 1, 0) == 7
@test Healpix.xyf2pixNest(resol, 3, 1, 0) == 8
@test Healpix.xyf2pixNest(resol, 0, 2, 0) == 9
@test Healpix.xyf2pixNest(resol, 1, 2, 0) == 10
@test Healpix.xyf2pixNest(resol, 0, 3, 0) == 11
@test Healpix.xyf2pixNest(resol, 1, 3, 0) == 12
@test Healpix.xyf2pixNest(resol, 2, 2, 0) == 13
@test Healpix.xyf2pixNest(resol, 3, 2, 0) == 14
@test Healpix.xyf2pixNest(resol, 2, 3, 0) == 15
@test Healpix.xyf2pixNest(resol, 3, 3, 0) == 16
@test Healpix.xyf2pixNest(resol, 0, 0, 1) == 17
@test Healpix.xyf2pixNest(resol, 1, 0, 1) == 18
@test Healpix.xyf2pixNest(resol, 0, 1, 1) == 19
@test Healpix.xyf2pixNest(resol, 1, 1, 1) == 20
@test Healpix.xyf2pixNest(resol, 2, 0, 1) == 21
@test Healpix.xyf2pixNest(resol, 3, 0, 1) == 22
@test Healpix.xyf2pixNest(resol, 2, 1, 1) == 23
@test Healpix.xyf2pixNest(resol, 3, 1, 1) == 24
@test Healpix.xyf2pixNest(resol, 0, 2, 1) == 25
@test Healpix.xyf2pixNest(resol, 1, 2, 1) == 26
@test Healpix.xyf2pixNest(resol, 0, 3, 1) == 27
@test Healpix.xyf2pixNest(resol, 1, 3, 1) == 28
@test Healpix.xyf2pixNest(resol, 2, 2, 1) == 29
@test Healpix.xyf2pixNest(resol, 3, 2, 1) == 30
@test Healpix.xyf2pixNest(resol, 2, 3, 1) == 31
@test Healpix.xyf2pixNest(resol, 3, 3, 1) == 32
@test Healpix.xyf2pixNest(resol, 0, 0, 2) == 33
@test Healpix.xyf2pixNest(resol, 1, 0, 2) == 34
@test Healpix.xyf2pixNest(resol, 0, 1, 2) == 35
@test Healpix.xyf2pixNest(resol, 1, 1, 2) == 36
@test Healpix.xyf2pixNest(resol, 2, 0, 2) == 37
@test Healpix.xyf2pixNest(resol, 3, 0, 2) == 38
@test Healpix.xyf2pixNest(resol, 2, 1, 2) == 39
@test Healpix.xyf2pixNest(resol, 3, 1, 2) == 40
@test Healpix.xyf2pixNest(resol, 0, 2, 2) == 41
@test Healpix.xyf2pixNest(resol, 1, 2, 2) == 42
@test Healpix.xyf2pixNest(resol, 0, 3, 2) == 43
@test Healpix.xyf2pixNest(resol, 1, 3, 2) == 44
@test Healpix.xyf2pixNest(resol, 2, 2, 2) == 45
@test Healpix.xyf2pixNest(resol, 3, 2, 2) == 46
@test Healpix.xyf2pixNest(resol, 2, 3, 2) == 47
@test Healpix.xyf2pixNest(resol, 3, 3, 2) == 48
@test Healpix.xyf2pixNest(resol, 0, 0, 3) == 49
@test Healpix.xyf2pixNest(resol, 1, 0, 3) == 50
@test Healpix.xyf2pixNest(resol, 0, 1, 3) == 51
@test Healpix.xyf2pixNest(resol, 1, 1, 3) == 52
@test Healpix.xyf2pixNest(resol, 2, 0, 3) == 53
@test Healpix.xyf2pixNest(resol, 3, 0, 3) == 54
@test Healpix.xyf2pixNest(resol, 2, 1, 3) == 55
@test Healpix.xyf2pixNest(resol, 3, 1, 3) == 56
@test Healpix.xyf2pixNest(resol, 0, 2, 3) == 57
@test Healpix.xyf2pixNest(resol, 1, 2, 3) == 58
@test Healpix.xyf2pixNest(resol, 0, 3, 3) == 59
@test Healpix.xyf2pixNest(resol, 1, 3, 3) == 60
@test Healpix.xyf2pixNest(resol, 2, 2, 3) == 61
@test Healpix.xyf2pixNest(resol, 3, 2, 3) == 62
@test Healpix.xyf2pixNest(resol, 2, 3, 3) == 63
@test Healpix.xyf2pixNest(resol, 3, 3, 3) == 64
@test Healpix.xyf2pixNest(resol, 0, 0, 4) == 65
@test Healpix.xyf2pixNest(resol, 1, 0, 4) == 66
@test Healpix.xyf2pixNest(resol, 0, 1, 4) == 67
@test Healpix.xyf2pixNest(resol, 1, 1, 4) == 68
@test Healpix.xyf2pixNest(resol, 2, 0, 4) == 69
@test Healpix.xyf2pixNest(resol, 3, 0, 4) == 70
@test Healpix.xyf2pixNest(resol, 2, 1, 4) == 71
@test Healpix.xyf2pixNest(resol, 3, 1, 4) == 72
@test Healpix.xyf2pixNest(resol, 0, 2, 4) == 73
@test Healpix.xyf2pixNest(resol, 1, 2, 4) == 74
@test Healpix.xyf2pixNest(resol, 0, 3, 4) == 75
@test Healpix.xyf2pixNest(resol, 1, 3, 4) == 76
@test Healpix.xyf2pixNest(resol, 2, 2, 4) == 77
@test Healpix.xyf2pixNest(resol, 3, 2, 4) == 78
@test Healpix.xyf2pixNest(resol, 2, 3, 4) == 79
@test Healpix.xyf2pixNest(resol, 3, 3, 4) == 80
@test Healpix.xyf2pixNest(resol, 0, 0, 5) == 81
@test Healpix.xyf2pixNest(resol, 1, 0, 5) == 82
@test Healpix.xyf2pixNest(resol, 0, 1, 5) == 83
@test Healpix.xyf2pixNest(resol, 1, 1, 5) == 84
@test Healpix.xyf2pixNest(resol, 2, 0, 5) == 85
@test Healpix.xyf2pixNest(resol, 3, 0, 5) == 86
@test Healpix.xyf2pixNest(resol, 2, 1, 5) == 87
@test Healpix.xyf2pixNest(resol, 3, 1, 5) == 88
@test Healpix.xyf2pixNest(resol, 0, 2, 5) == 89
@test Healpix.xyf2pixNest(resol, 1, 2, 5) == 90
@test Healpix.xyf2pixNest(resol, 0, 3, 5) == 91
@test Healpix.xyf2pixNest(resol, 1, 3, 5) == 92
@test Healpix.xyf2pixNest(resol, 2, 2, 5) == 93
@test Healpix.xyf2pixNest(resol, 3, 2, 5) == 94
@test Healpix.xyf2pixNest(resol, 2, 3, 5) == 95
@test Healpix.xyf2pixNest(resol, 3, 3, 5) == 96
@test Healpix.xyf2pixNest(resol, 0, 0, 6) == 97
@test Healpix.xyf2pixNest(resol, 1, 0, 6) == 98
@test Healpix.xyf2pixNest(resol, 0, 1, 6) == 99
@test Healpix.xyf2pixNest(resol, 1, 1, 6) == 100
@test Healpix.xyf2pixNest(resol, 2, 0, 6) == 101
@test Healpix.xyf2pixNest(resol, 3, 0, 6) == 102
@test Healpix.xyf2pixNest(resol, 2, 1, 6) == 103
@test Healpix.xyf2pixNest(resol, 3, 1, 6) == 104
@test Healpix.xyf2pixNest(resol, 0, 2, 6) == 105
@test Healpix.xyf2pixNest(resol, 1, 2, 6) == 106
@test Healpix.xyf2pixNest(resol, 0, 3, 6) == 107
@test Healpix.xyf2pixNest(resol, 1, 3, 6) == 108
@test Healpix.xyf2pixNest(resol, 2, 2, 6) == 109
@test Healpix.xyf2pixNest(resol, 3, 2, 6) == 110
@test Healpix.xyf2pixNest(resol, 2, 3, 6) == 111
@test Healpix.xyf2pixNest(resol, 3, 3, 6) == 112
@test Healpix.xyf2pixNest(resol, 0, 0, 7) == 113
@test Healpix.xyf2pixNest(resol, 1, 0, 7) == 114
@test Healpix.xyf2pixNest(resol, 0, 1, 7) == 115
@test Healpix.xyf2pixNest(resol, 1, 1, 7) == 116
@test Healpix.xyf2pixNest(resol, 2, 0, 7) == 117
@test Healpix.xyf2pixNest(resol, 3, 0, 7) == 118
@test Healpix.xyf2pixNest(resol, 2, 1, 7) == 119
@test Healpix.xyf2pixNest(resol, 3, 1, 7) == 120
@test Healpix.xyf2pixNest(resol, 0, 2, 7) == 121
@test Healpix.xyf2pixNest(resol, 1, 2, 7) == 122
@test Healpix.xyf2pixNest(resol, 0, 3, 7) == 123
@test Healpix.xyf2pixNest(resol, 1, 3, 7) == 124
@test Healpix.xyf2pixNest(resol, 2, 2, 7) == 125
@test Healpix.xyf2pixNest(resol, 3, 2, 7) == 126
@test Healpix.xyf2pixNest(resol, 2, 3, 7) == 127
@test Healpix.xyf2pixNest(resol, 3, 3, 7) == 128
@test Healpix.xyf2pixNest(resol, 0, 0, 8) == 129
@test Healpix.xyf2pixNest(resol, 1, 0, 8) == 130
@test Healpix.xyf2pixNest(resol, 0, 1, 8) == 131
@test Healpix.xyf2pixNest(resol, 1, 1, 8) == 132
@test Healpix.xyf2pixNest(resol, 2, 0, 8) == 133
@test Healpix.xyf2pixNest(resol, 3, 0, 8) == 134
@test Healpix.xyf2pixNest(resol, 2, 1, 8) == 135
@test Healpix.xyf2pixNest(resol, 3, 1, 8) == 136
@test Healpix.xyf2pixNest(resol, 0, 2, 8) == 137
@test Healpix.xyf2pixNest(resol, 1, 2, 8) == 138
@test Healpix.xyf2pixNest(resol, 0, 3, 8) == 139
@test Healpix.xyf2pixNest(resol, 1, 3, 8) == 140
@test Healpix.xyf2pixNest(resol, 2, 2, 8) == 141
@test Healpix.xyf2pixNest(resol, 3, 2, 8) == 142
@test Healpix.xyf2pixNest(resol, 2, 3, 8) == 143
@test Healpix.xyf2pixNest(resol, 3, 3, 8) == 144
@test Healpix.xyf2pixNest(resol, 0, 0, 9) == 145
@test Healpix.xyf2pixNest(resol, 1, 0, 9) == 146
@test Healpix.xyf2pixNest(resol, 0, 1, 9) == 147
@test Healpix.xyf2pixNest(resol, 1, 1, 9) == 148
@test Healpix.xyf2pixNest(resol, 2, 0, 9) == 149
@test Healpix.xyf2pixNest(resol, 3, 0, 9) == 150
@test Healpix.xyf2pixNest(resol, 2, 1, 9) == 151
@test Healpix.xyf2pixNest(resol, 3, 1, 9) == 152
@test Healpix.xyf2pixNest(resol, 0, 2, 9) == 153
@test Healpix.xyf2pixNest(resol, 1, 2, 9) == 154
@test Healpix.xyf2pixNest(resol, 0, 3, 9) == 155
@test Healpix.xyf2pixNest(resol, 1, 3, 9) == 156
@test Healpix.xyf2pixNest(resol, 2, 2, 9) == 157
@test Healpix.xyf2pixNest(resol, 3, 2, 9) == 158
@test Healpix.xyf2pixNest(resol, 2, 3, 9) == 159
@test Healpix.xyf2pixNest(resol, 3, 3, 9) == 160
@test Healpix.xyf2pixNest(resol, 0, 0, 10) == 161
@test Healpix.xyf2pixNest(resol, 1, 0, 10) == 162
@test Healpix.xyf2pixNest(resol, 0, 1, 10) == 163
@test Healpix.xyf2pixNest(resol, 1, 1, 10) == 164
@test Healpix.xyf2pixNest(resol, 2, 0, 10) == 165
@test Healpix.xyf2pixNest(resol, 3, 0, 10) == 166
@test Healpix.xyf2pixNest(resol, 2, 1, 10) == 167
@test Healpix.xyf2pixNest(resol, 3, 1, 10) == 168
@test Healpix.xyf2pixNest(resol, 0, 2, 10) == 169
@test Healpix.xyf2pixNest(resol, 1, 2, 10) == 170
@test Healpix.xyf2pixNest(resol, 0, 3, 10) == 171
@test Healpix.xyf2pixNest(resol, 1, 3, 10) == 172
@test Healpix.xyf2pixNest(resol, 2, 2, 10) == 173
@test Healpix.xyf2pixNest(resol, 3, 2, 10) == 174
@test Healpix.xyf2pixNest(resol, 2, 3, 10) == 175
@test Healpix.xyf2pixNest(resol, 3, 3, 10) == 176
@test Healpix.xyf2pixNest(resol, 0, 0, 11) == 177
@test Healpix.xyf2pixNest(resol, 1, 0, 11) == 178
@test Healpix.xyf2pixNest(resol, 0, 1, 11) == 179
@test Healpix.xyf2pixNest(resol, 1, 1, 11) == 180
@test Healpix.xyf2pixNest(resol, 2, 0, 11) == 181
@test Healpix.xyf2pixNest(resol, 3, 0, 11) == 182
@test Healpix.xyf2pixNest(resol, 2, 1, 11) == 183
@test Healpix.xyf2pixNest(resol, 3, 1, 11) == 184
@test Healpix.xyf2pixNest(resol, 0, 2, 11) == 185
@test Healpix.xyf2pixNest(resol, 1, 2, 11) == 186
@test Healpix.xyf2pixNest(resol, 0, 3, 11) == 187
@test Healpix.xyf2pixNest(resol, 1, 3, 11) == 188
@test Healpix.xyf2pixNest(resol, 2, 2, 11) == 189
@test Healpix.xyf2pixNest(resol, 3, 2, 11) == 190
@test Healpix.xyf2pixNest(resol, 2, 3, 11) == 191
@test Healpix.xyf2pixNest(resol, 3, 3, 11) == 192

@test Healpix.xyf2pixNest(highresol, 0, 0, 0) == 1
@test Healpix.xyf2pixNest(highresol, 1, 0, 0) == 2
@test Healpix.xyf2pixNest(highresol, 0, 1, 0) == 3
@test Healpix.xyf2pixNest(highresol, 1, 1, 0) == 4
@test Healpix.xyf2pixNest(highresol, 2, 0, 0) == 5
@test Healpix.xyf2pixNest(highresol, 3, 0, 0) == 6
@test Healpix.xyf2pixNest(highresol, 2, 1, 0) == 7
@test Healpix.xyf2pixNest(highresol, 3, 1, 0) == 8
@test Healpix.xyf2pixNest(highresol, 0, 2, 0) == 9
@test Healpix.xyf2pixNest(highresol, 1, 2, 0) == 10
@test Healpix.xyf2pixNest(highresol, 0, 0, 3) == 864691128455135233
@test Healpix.xyf2pixNest(highresol, 1, 0, 3) == 864691128455135234
@test Healpix.xyf2pixNest(highresol, 0, 1, 3) == 864691128455135235
@test Healpix.xyf2pixNest(highresol, 1, 1, 3) == 864691128455135236
@test Healpix.xyf2pixNest(highresol, 2, 0, 3) == 864691128455135237
@test Healpix.xyf2pixNest(highresol, 3, 0, 3) == 864691128455135238
@test Healpix.xyf2pixNest(highresol, 2, 1, 3) == 864691128455135239
@test Healpix.xyf2pixNest(highresol, 3, 1, 3) == 864691128455135240
@test Healpix.xyf2pixNest(highresol, 0, 2, 3) == 864691128455135241
@test Healpix.xyf2pixNest(highresol, 1, 2, 3) == 864691128455135242
@test Healpix.xyf2pixNest(highresol, 0, 0, 6) == 1729382256910270465
@test Healpix.xyf2pixNest(highresol, 1, 0, 6) == 1729382256910270466
@test Healpix.xyf2pixNest(highresol, 0, 1, 6) == 1729382256910270467
@test Healpix.xyf2pixNest(highresol, 1, 1, 6) == 1729382256910270468
@test Healpix.xyf2pixNest(highresol, 2, 0, 6) == 1729382256910270469
@test Healpix.xyf2pixNest(highresol, 3, 0, 6) == 1729382256910270470
@test Healpix.xyf2pixNest(highresol, 2, 1, 6) == 1729382256910270471
@test Healpix.xyf2pixNest(highresol, 3, 1, 6) == 1729382256910270472
@test Healpix.xyf2pixNest(highresol, 0, 2, 6) == 1729382256910270473
@test Healpix.xyf2pixNest(highresol, 1, 2, 6) == 1729382256910270474
@test Healpix.xyf2pixNest(highresol, 0, 0, 9) == 2594073385365405697
@test Healpix.xyf2pixNest(highresol, 1, 0, 9) == 2594073385365405698
@test Healpix.xyf2pixNest(highresol, 0, 1, 9) == 2594073385365405699
@test Healpix.xyf2pixNest(highresol, 1, 1, 9) == 2594073385365405700
@test Healpix.xyf2pixNest(highresol, 2, 0, 9) == 2594073385365405701
@test Healpix.xyf2pixNest(highresol, 3, 0, 9) == 2594073385365405702
@test Healpix.xyf2pixNest(highresol, 2, 1, 9) == 2594073385365405703
@test Healpix.xyf2pixNest(highresol, 3, 1, 9) == 2594073385365405704
@test Healpix.xyf2pixNest(highresol, 0, 2, 9) == 2594073385365405705
@test Healpix.xyf2pixNest(highresol, 1, 2, 9) == 2594073385365405706
@test Healpix.xyf2pixNest(highresol, 536870911, 536870908, 11) == 3458764513820540918
@test Healpix.xyf2pixNest(highresol, 536870910, 536870909, 11) == 3458764513820540919
@test Healpix.xyf2pixNest(highresol, 536870911, 536870909, 11) == 3458764513820540920
@test Healpix.xyf2pixNest(highresol, 536870908, 536870910, 11) == 3458764513820540921
@test Healpix.xyf2pixNest(highresol, 536870909, 536870910, 11) == 3458764513820540922
@test Healpix.xyf2pixNest(highresol, 536870908, 536870911, 11) == 3458764513820540923
@test Healpix.xyf2pixNest(highresol, 536870909, 536870911, 11) == 3458764513820540924
@test Healpix.xyf2pixNest(highresol, 536870910, 536870910, 11) == 3458764513820540925
@test Healpix.xyf2pixNest(highresol, 536870911, 536870910, 11) == 3458764513820540926
@test Healpix.xyf2pixNest(highresol, 536870910, 536870911, 11) == 3458764513820540927

@test Healpix.pix2xyfNest(highresol, 1) == (0, 0, 0)
@test Healpix.pix2xyfNest(highresol, 2) == (1, 0, 0)
@test Healpix.pix2xyfNest(highresol, 3) == (0, 1, 0)
@test Healpix.pix2xyfNest(highresol, 4) == (1, 1, 0)
@test Healpix.pix2xyfNest(highresol, 5) == (2, 0, 0)
@test Healpix.pix2xyfNest(highresol, 6) == (3, 0, 0)
@test Healpix.pix2xyfNest(highresol, 7) == (2, 1, 0)
@test Healpix.pix2xyfNest(highresol, 8) == (3, 1, 0)
@test Healpix.pix2xyfNest(highresol, 9) == (0, 2, 0)
@test Healpix.pix2xyfNest(highresol, 10) == (1, 2, 0)
@test Healpix.pix2xyfNest(highresol, 864691128455135233) == (0, 0, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135234) == (1, 0, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135235) == (0, 1, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135236) == (1, 1, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135237) == (2, 0, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135238) == (3, 0, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135239) == (2, 1, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135240) == (3, 1, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135241) == (0, 2, 3)
@test Healpix.pix2xyfNest(highresol, 864691128455135242) == (1, 2, 3)
@test Healpix.pix2xyfNest(highresol, 1729382256910270465) == (0, 0, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270466) == (1, 0, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270467) == (0, 1, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270468) == (1, 1, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270469) == (2, 0, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270470) == (3, 0, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270471) == (2, 1, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270472) == (3, 1, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270473) == (0, 2, 6)
@test Healpix.pix2xyfNest(highresol, 1729382256910270474) == (1, 2, 6)
@test Healpix.pix2xyfNest(highresol, 2594073385365405697) == (0, 0, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405698) == (1, 0, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405699) == (0, 1, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405700) == (1, 1, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405701) == (2, 0, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405702) == (3, 0, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405703) == (2, 1, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405704) == (3, 1, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405705) == (0, 2, 9)
@test Healpix.pix2xyfNest(highresol, 2594073385365405706) == (1, 2, 9)
@test Healpix.pix2xyfNest(highresol, 3458764513820540918) == (536870911, 536870908, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540919) == (536870910, 536870909, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540920) == (536870911, 536870909, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540921) == (536870908, 536870910, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540922) == (536870909, 536870910, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540923) == (536870908, 536870911, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540924) == (536870909, 536870911, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540925) == (536870910, 536870910, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540926) == (536870911, 536870910, 11)
@test Healpix.pix2xyfNest(highresol, 3458764513820540927) == (536870910, 536870911, 11)

@test Healpix.ring2nest(Healpix.Resolution(32), 1643) == 436
@test Healpix.nest2ring(Healpix.Resolution(64), 1791) == 4933


################################################################################

(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 0)
@test z ≈ -4.3333333333e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 1)
@test z ≈ -4.3333333333e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 2)
@test z ≈ -4.3333333333e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 3)
@test z ≈ -4.3333333333e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 0)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 0)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 1)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 1)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 2)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 2)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 3)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 3)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 0)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 0)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 0)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 1)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 1)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 1)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 2)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 2)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 2)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 3)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 3)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 3)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 0)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 0)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 0)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 0)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 1)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 1)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 1)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 1)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 2)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 2)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 2)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 2)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 3)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 3)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 3)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 3)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 4)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 0)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 0)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 0)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 5)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 1)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 1)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 1)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 6)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 2)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 2)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 2)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 7)
@test z ≈ -2.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 3)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 3)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 3)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 4)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 0)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 0)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 5)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 5)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 1)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 1)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 6)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 6)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 2)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 2)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 7)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 7)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 3)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 3)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 4)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 4)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 4)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 0)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 7.8539816340e-01; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 5)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 5)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 5)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 1)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 2.3561944902e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 6)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 6)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 6)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 2)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 3.9269908170e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 7)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 7)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 7)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 3)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 5.4977871438e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 4)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 4)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 4)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 5)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 5)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 5)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 5)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 6)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 6)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 6)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 6)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 7)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 7)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 7)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 7)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 4)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 4)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 4)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 4)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 8)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 5)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 5)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 5)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 9)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 6)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 6)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 6)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 10)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 7)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 7)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 7)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 3, 11)
@test z ≈ -3.3333333333e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 4)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 4)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 7.8539816340e-01; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 8)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 8)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 5)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 7.8539816340e-01; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 5)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 2.3561944902e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 9)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 9)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 6)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 2.3561944902e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 6)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 3.9269908170e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 10)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 10)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 7)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 3.9269908170e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 7)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 5.4977871438e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 3, 11)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 2, 11)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 4)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 5.4977871438e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 4)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 8)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 8)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 8)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 5)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 9)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 9)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 9)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 6)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 10)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 10)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 10)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 7)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 3, 11)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 2, 11)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 1, 11)
@test z ≈ 1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 8)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 8)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 8)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 8)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 9)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 9)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 9)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 9)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 10)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 10)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 10)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 10)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 3, 11)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 2, 11)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 1, 11)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(3, 0, 11)
@test z ≈ 6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 8)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 5.4977871438e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 8)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 7.8539816340e-01; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 8)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 2.3561944902e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 9)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 7.8539816340e-01; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 9)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 2.3561944902e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 9)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 3.9269908170e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 10)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 2.3561944902e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 10)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 3.9269908170e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 10)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 5.4977871438e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 2, 11)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 3.9269908170e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 1, 11)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 5.4977871438e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(2, 0, 11)
@test z ≈ 0.0000000000e+00 ; @test phi ≈ 7.8539816340e-01; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 8)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 8)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 9)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 1.5707963268e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 9)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 10)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 3.1415926536e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 10)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 1, 11)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 4.7123889804e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(1, 0, 11)
@test z ≈ -6.6666666667e-01 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 8)
@test z ≈ -1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 9)
@test z ≈ -1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 10)
@test z ≈ -1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
(z, phi, sintheta, have_sintheta) = Healpix.xyf2loc(0, 0, 11)
@test z ≈ -1.0000000000e+00 ; @test phi ≈ 0.0000000000e+00; have_sintheta && (@test sintheta ≈ 0.0000000000e+00)
